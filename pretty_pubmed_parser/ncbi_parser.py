import logging, sys
from typing import List, Dict, Optional 
from datetime import datetime
from dateutil.relativedelta import relativedelta

from pydantic import BaseModel, HttpUrl
from Bio import Entrez

from configs import ENTREZ_EMAIL, BATCH_SIZE, TIME_BATCH_MONTHS, START_YEAR

class ReadyAuthor(BaseModel):

    forename: str
    lastname: str

class ReadyPaper(BaseModel):

    title: str
    pmid: int
    doi: str
    link: HttpUrl

    abstract: str
    journal: str

    authors: Optional[List[ReadyAuthor]]

def setup_ncbi() -> logging.Logger:
    """Setups NCBI Entrez module and returns logger."""
    Entrez.email = ENTREZ_EMAIL

    pubmed_logger = logging.getLogger('pubmed')
    pubmed_logger.addHandler(logging.StreamHandler(sys.stdout))
    pubmed_logger.setLevel(logging.INFO)

    return pubmed_logger



def search(query: str, results_num: int = BATCH_SIZE, sort: str = 'pub_date') -> dict:
    handle = Entrez.esearch(
        db='pubmed', 
        sort=sort, 
        retmax=str(results_num),
        retmode='xml', 
        term=query
    )
    results = dict(Entrez.read(handle))

    return results

def search_by_journal(journal_name: str, results_num: int = BATCH_SIZE, sort: str = 'pub_date') -> dict:
    handle = Entrez.esearch(
        db='pubmed', 
        sort=sort, 
        retmax=str(results_num),
        retmode='xml', 
        term=f'"{journal_name}"[Journal]'
    )
    results = dict(Entrez.read(handle))

    return results

def search_num_pages(query: str) -> int:
    handle = Entrez.esearch(
        db='pubmed', 
        retmax=str(10),
        retmode='xml', 
        term=query
    )
    results = dict(Entrez.read(handle))

    return  (int(results['Count']) // BATCH_SIZE) + 1

def search_by_journal_num_pages(journal_name: str) -> int:
    handle = Entrez.esearch(
        db='pubmed', 
        retmax=str(10),
        retmode='xml', 
        term=f'"{journal_name}"[Journal]'
    )
    results = dict(Entrez.read(handle))

    return  (int(results['Count']) // BATCH_SIZE) + 1



# TODO DRY, probably decorators or OOP?

# def get_all_pmids(query: str, sort='pub_date') -> dict:

#     today = datetime.today()
#     date = datetime(year=START_YEAR, month=1, day=1)
#     pmids = []

#     while date <= today:

#         new_date = date + relativedelta(months=TIME_BATCH_MONTHS)

#         min_month = date.strftime('%Y/%m')
#         max_month = new_date.strftime('%Y/%m')

#         batch_results = search(query, results_num=10000)

#         pmids += list(batch_results)        
        

#     handle = Entrez.esearch(
#         db='pubmed', 
#         sort=sort,
#         retmode='uilist', 
#         term=query
#     )
#     results = dict(Entrez.read(handle))

#     return results


def fetch_details(id_list):
    ids = ','.join(id_list)
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    details = Entrez.read(handle)
    return details

def clean_text(s: str) -> str:
    
    s = ' '.join(s.replace('/', ' / ').split())
    s = s.replace('!', '\!')
    
    return s


def get_doi(article_id_list) -> str:
    for string_element in article_id_list:
        if string_element.__dict__['attributes']['IdType'] == 'doi':
            return string_element.__str__().lower()


def get_abstract(paper) -> str:
    
    string_list = []

    try:
        abstract = paper['Abstract']['AbstractText']
    except KeyError:
        return 'No abstract available for this paper :('
    
    for part in abstract:
        if part.__dict__['attributes']:
            string_list.append('*' + part.__dict__['attributes']['Label'] + '*')
        string_list.append(part.__str__())
        string_list.append('')

    cleaned_list = [clean_text(s) for s in string_list]
    
    return '\n'.join(cleaned_list)

def parse_paper(paper) -> ReadyPaper:

    title = clean_text(paper['MedlineCitation']['Article']['ArticleTitle'])
    pmid = int(paper['MedlineCitation']['PMID'])
    doi = get_doi(paper['PubmedData']['ArticleIdList'])
    link = f'https://doi.org/{doi}'

    abstract = get_abstract(paper['MedlineCitation']['Article'])
    journal = paper['MedlineCitation']['Article']['Journal']['Title']

    authors = []
    if 'AuthorList' in paper['MedlineCitation']['Article'].keys():
        for author in paper['MedlineCitation']['Article']['AuthorList']:
            if author.attributes['ValidYN'] == 'Y':
                authors.append(
                    ReadyAuthor(
                        forename=author['ForeName'],
                        lastname=author['LastName']
                    )
                )
                

    ready_paper = ReadyPaper(
        title=title,
        pmid=pmid,
        doi=doi,
        link=link,
        abstract=abstract,
        journal=journal,
        authors = None if not authors else authors
    )

    if authors:
        ready_paper.authors = authors=authors

    return ready_paper

def get_ready_paper_by_id(id: int):
    papers = fetch_details([id])['PubmedArticle']
    return parse_paper(papers[0])


def get_ready_papers(
        pubmed_logger: logging.Logger,
        query: str,
        max_results_num: int = BATCH_SIZE,
        ids: Optional[List[int]] = None,
        esearch_args: Optional[Dict[str: str]] = None,
) -> Dict[int: ReadyPaper]:
    
    if esearch_args is None:
        esearch_args = {}
    
    if not ids:
        results = search(query, max_results_num)
        ids = results['IdList']

    output = {}
    papers = fetch_details(ids)['PubmedArticle']      
    for paper in papers:
        ready_paper = None
        try:
            ready_paper = parse_paper(paper)
        except Exception as e:
            pubmed_logger.error(f'Error while parsing paper matching query "{query}"')
            pubmed_logger.error(e)
            continue
        if ready_paper:
            output[ready_paper.pmid] = ready_paper
            


    


    

def get_papers_by_query(pubmed_logger: logging.Logger, query: str, results_num: int = BATCH_SIZE) -> List[ReadyPaper]:

    results = search(query, results_num)
    ids = results['IdList']
    papers = fetch_details(ids)['PubmedArticle']

    if len(papers) < results_num:
        results = search(query, results_num * 2)
        ids = results['IdList']
        papers = fetch_details(ids)['PubmedArticle']

    outuput = []

    for paper in papers:
        try:
            outuput.append(parse_paper(paper))
        except Exception as e:
            pubmed_logger.error(f'Error while parsing paper about {query}!')
            pubmed_logger.error(e)
            continue

    if len(outuput) < results_num:
        results = search(query, results_num * 2)
        ids = results['IdList']
        papers = fetch_details(ids)['PubmedArticle']

    for paper in papers:
        try:
            outuput.append(parse_paper(paper))
        except Exception as e:
            pubmed_logger.error(f'Error while parsing paper about {query}!')
            pubmed_logger.error(e)
            continue

    if len(outuput) < results_num:
        pubmed_logger.warning(f'Not enough papers about {query}, found only {len(outuput)}!')

    if len(outuput) > results_num:
        outuput = outuput[:results_num]

    return outuput

def get_papers_by_journal(pubmed_logger: logging.Logger, journal_name: str, results_num = BATCH_SIZE):

    results = search_by_journal(journal_name, results_num)
    ids = results['IdList']
    papers = fetch_details(ids)['PubmedArticle']

    if len(papers) < results_num:
        results = search_by_journal(journal_name, results_num * 2)
        ids = results['IdList']
        papers = fetch_details(ids)['PubmedArticle']

    outuput = []

    for paper in papers:
        try:
            outuput.append(parse_paper(paper))
        except Exception as e:
            pubmed_logger.error(f'Error while parsing papers from journal "{journal_name}"!')
            pubmed_logger.error(e)
            continue

    if len(outuput) < results_num:
        results = search_by_journal(journal_name, results_num * 2)
        ids = results['IdList']
        papers = fetch_details(ids)['PubmedArticle']

    for paper in papers:
        try:
            outuput.append(parse_paper(paper))
        except Exception as e:
            pubmed_logger.error(f'Error while parsing papers from journal "{journal_name}"!')
            pubmed_logger.error(e)
            continue

    if len(outuput) < results_num:
        pubmed_logger.warning(f'Not enough papers from journal "{journal_name}", found only {len(outuput)}!')

    if len(outuput) > results_num:
        outuput = outuput[:results_num]

    return outuput

    # TODO DRY

    # TODO
    # DATETIME
    # PAGINATION:
    # - create function to get the number of pages
    # - add page_num param to the get_papers_ funcs
    # - divide into batches, matching by ids
    # - query these papers (separate id list creation and further operations into separate funcs)
    # - if there are some errors, some pages will be just less than 100, it's ok