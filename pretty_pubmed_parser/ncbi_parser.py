import logging, sys
from typing import List, Dict, Optional 
from datetime import datetime
from dateutil.relativedelta import relativedelta
from functools import lru_cache
from frozendict import frozendict

from pydantic import BaseModel, HttpUrl
from fastapi import HTTPException
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
    date: str

    authors: Optional[List[ReadyAuthor]]

def setup_ncbi() -> logging.Logger:
    """Setups NCBI Entrez module and returns logger."""
    Entrez.email = ENTREZ_EMAIL

    pubmed_logger = logging.getLogger('pubmed')
    pubmed_logger.addHandler(logging.StreamHandler(sys.stdout))
    pubmed_logger.setLevel(logging.INFO)

    return pubmed_logger

def generate_journal_query(journal_name: str) -> str:
    return f'"{journal_name}"[Journal]'

def generate_handle_args(
        query: str, max_results_num: str, esearch_args: Optional[Dict[str, str]] = None,
        min_date: Optional[str] = None, max_date: Optional[str] = None
):
    
    if min_date is None:
        date = datetime(year=START_YEAR, month=1, day=1)
        min_date = date.strftime('%Y/%m/%d')

    if max_date is None:
        today = datetime.today()
        max_date = today.strftime('%Y/%m/%d')

    handle_args = {
        'db': 'pubmed',
        'retmode': 'xml',
        'sort': 'pub_date',
        'term': query,
        'retmax': str(max_results_num),
        'mindate': min_date,
        'maxdate': max_date
    }

    if esearch_args:
        handle_args.update(esearch_args)

    return handle_args


def extract_num_pages(handle_args: dict) -> int:

    handle_args['retmax'] = str(5)
    handle = Entrez.esearch(**handle_args)
    results = dict(Entrez.read(handle))

    print(results)

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


def extract_doi(article_id_list) -> str:
    for string_element in article_id_list:
        if string_element.__dict__['attributes']['IdType'] == 'doi':
            return string_element.__str__().lower()


def extract_abstract(paper) -> str:
    
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

def extract_date(paper):

    date_dict = dict(paper['MedlineCitation']['Article']['ArticleDate'][0])
    date_str= f'{date_dict["Year"]}-{date_dict["Month"]}-{date_dict["Day"]}'

    return date_str

def parse_paper(paper) -> ReadyPaper:

    title = clean_text(paper['MedlineCitation']['Article']['ArticleTitle'])
    pmid = int(paper['MedlineCitation']['PMID'])
    doi = extract_doi(paper['PubmedData']['ArticleIdList'])
    link = f'https://doi.org/{doi}'

    abstract = extract_abstract(paper['MedlineCitation']['Article'])
    journal = paper['MedlineCitation']['Article']['Journal']['Title']
    date = extract_date(paper)

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
        date=date,
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
        esearch_args: Optional[Dict[str, str]] = None,
) -> Dict[int, ReadyPaper]:
    
    if esearch_args and ids:
        raise ValueError("You cannot pass id list and esearch arguments at the same time!")
    
    if esearch_args is None:
        esearch_args = {}
    
    if not ids:

        handle_args = generate_handle_args(query, max_results_num, esearch_args=esearch_args)

        results = dict(
            Entrez.read(Entrez.esearch(**handle_args))
        )

        ids = results['IdList']

    output = {}
    papers = fetch_details(ids)['PubmedArticle']      
    for paper in papers:
        ready_paper = None
        try:
            ready_paper = parse_paper(paper)
        except Exception as e:
            pubmed_logger.error(f'Error while parsing paper matching query "{query}. We will try to parse it again."')
            pubmed_logger.error(e)
            continue
        if ready_paper:
            output[ready_paper.pmid] = ready_paper

    # try again for not found ids
    for not_parsed_id in set(ids) - set(output.keys()):
        try:
            ready_paper = get_ready_paper_by_id(not_parsed_id)
        except Exception as e:
            pubmed_logger.error(f'Error while parsing paper with id "{not_parsed_id}" matching query "{query}". No further attempts will be performed.')
            pubmed_logger.error(e)
            continue

    if len(output) < len(ids):
        pubmed_logger.warning(f'Not enough papers for query "{query}", found only {len(output)}!')


    return output

@lru_cache(maxsize=32)
def get_batch_by_num_pages(handle_args: frozendict, num_pages: int) -> list:

    if num_pages * BATCH_SIZE + 1 < 10000:
        results = dict(
            Entrez.read(Entrez.esearch(**handle_args))
        )

        id_list = list(results['IdList'])

    else:
        raise HTTPException(status_code=400, detail='Too much results found! Please, specify your query (you can narrow the dates diapason).')
    
    return id_list

    # TODO
    # DATETIME - DONE!
    # PAGINATION:
    # - create function to get the number of pages - DONE!

    # - add page_num param to the get_papers_ funcs
    # - divide into batches, matching by ids
    # - query these papers (separate id list creation and further operations into separate funcs)
    # - if there are some errors, some pages will be just less than 100, it's ok

    # NEWER IMPLEMENTATION:
    # - Create ids param to get_ready_papers - DONE
    # - Generate ids and save into Python list?

    # HOW TO ? WE