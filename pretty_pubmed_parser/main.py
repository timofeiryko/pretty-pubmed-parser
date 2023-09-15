from typing import Dict

import uvicorn
from fastapi import FastAPI

from ncbi_parser import ReadyPaper, setup_ncbi, generate_journal_query, generate_handle_args, get_ready_papers, extract_num_pages
from configs import BATCH_SIZE

app = FastAPI()

@app.get("/")
async def root():
    return {"message": "Pretty PubMed Parser API"}

@app.get("/search/{query}")
async def search(query: str, results_num: int = BATCH_SIZE, response_model=Dict[int, ReadyPaper]):

    pubmed_logger = setup_ncbi()
        
    result = get_ready_papers(pubmed_logger, query, results_num)

    return result

@app.get("/search/pages_num/{query}")
async def search_num_pages(query: str, response_model=int):

    setup_ncbi()

    handle_args = generate_handle_args(query, max_results_num=BATCH_SIZE)
    result = extract_num_pages(handle_args)

    return result


@app.get("/journal/{journal_name}")
async def search_journal(journal_name: str, results_num: int = BATCH_SIZE, response_model=Dict[int, ReadyPaper]):

    pubmed_logger = setup_ncbi()
    journal_query = generate_journal_query(journal_name)
    result = get_ready_papers(pubmed_logger, journal_query, results_num)

    return result

@app.get("/journal/pages_num/{journal_name}")
async def journal_num_pages(journal_name: str, response_model=int):

    setup_ncbi()

    journal_query = generate_journal_query(journal_name)
    handle_args = generate_handle_args(journal_query, max_results_num=BATCH_SIZE)
    result = extract_num_pages(handle_args)

    return result



if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)