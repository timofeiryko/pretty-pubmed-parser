import logging, sys

from fastapi import FastAPI
from pydantic import BaseModel, HttpUrl
from Bio import Entrez

from configs import ENTREZ_EMAIL

def setup_ncbi():
    """Setups NCBI Entrez module and logging."""
    Entrez.email = ENTREZ_EMAIL

    pubmed_logger = logging.getLogger('pubmed')
    pubmed_logger.addHandler(logging.StreamHandler(sys.stdout))
    pubmed_logger.setLevel(logging.INFO)

app = FastAPI()

class ReadyPaper(BaseModel):

    title: str
    pmid: int
    link: HttpUrl

    abstract: str
    text: str
    journal: str

    authors: list[str]

def search(query: str, results_num: int = 10, sort='pub+date') -> dict:
    handle = Entrez.esearch(
        db='pubmed', 
        sort='pub+date', 
        retmax=str(results_num),
        retmode='xml', 
        term=query
    )
    results = Entrez.read(handle)

    return results

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
    ready_paper = ReadyPaper(
        title=clean_text(paper['MedlineCitation']['Article']['ArticleTitle']),
        doi=get_doi(paper['PubmedData']['ArticleIdList']),
        pmid=int(paper['MedlineCitation']['PMID']),
        abstract=get_abstract(paper['MedlineCitation']['Article'])
    )
    return ready_paper

def get_ready_papers(query: str, results_num: int = 10):

    results = search(query, results_num)
    ids = results['IdList']
    papers = fetch_details(ids)['PubmedArticle']

    if len(papers) < results_num:
        results = search(query, results_num * 2)
        ids = results['IdList']
        papers = fetch_details(ids)['PubmedArticle']

    if len(papers) > results_num:
        papers = papers[:results_num]

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

    return outuput

@app.get("/")
async def root():
    return {"message": "Pretty PubMed Parser API"}