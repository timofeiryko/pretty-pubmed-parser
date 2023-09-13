from typing import List

import uvicorn
from fastapi import FastAPI

from ncbi_parser import ReadyPaper, setup_ncbi, get_papers_by_query, get_papers_by_journal
from configs import BATCH_SIZE

app = FastAPI()

@app.get("/")
async def root():
    return {"message": "Pretty PubMed Parser API"}

@app.get("/search/{query}")
async def search(query: str, results_num: int = BATCH_SIZE, response_model=List[ReadyPaper]):
    pubmed_logger = setup_ncbi()
    results = get_papers_by_query(pubmed_logger, query, results_num)
    return results

@app.get("/journal/{journal_name}")
async def search_journal(journal_name: str, results_num: int = BATCH_SIZE, response_model=List[ReadyPaper]):
    pubmed_logger = setup_ncbi()
    results = get_papers_by_journal(pubmed_logger, journal_name, results_num)
    return results

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)