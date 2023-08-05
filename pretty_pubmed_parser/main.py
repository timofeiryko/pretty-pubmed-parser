from typing import List

from fastapi import FastAPI

from ncbi_parser import ReadyPaper, setup_ncbi, get_ready_papers

app = FastAPI()

@app.get("/")
async def root():
    return {"message": "Pretty PubMed Parser API"}

@app.get("/search/{query}")
async def search(query: str, results_num: int = 10, response_model=List[ReadyPaper]):
    pubmed_logger = setup_ncbi()
    results = get_ready_papers(pubmed_logger, query)
    return results

if __name__ == '__main__':

    pubmed_logger = setup_ncbi()

    results = get_ready_papers(pubmed_logger, 'cancer', 10)

    for result in results:
        print(result)
        print()