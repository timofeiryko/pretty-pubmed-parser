# Using

## Getting started

### Download and install

Download this repository:

```bash
git clone git@github.com:timofeiryko/pretty-pubmed-parser.git
```

Then locate to the project directory:

```bash
cd pretty-pubmed-parser
```

You can install dependencies using poetry:

```bash
poetry install
```

Alternatively, you can use pip:

```bash
pip install -r requirements.txt
```

### Run

You can run this API using UVicorn:

```bash
uvicorn main:app --reload
```

To try it out, you can use the built-in FastAPI docs at `http://localhost:8000/docs`.

You can search by some keyword or by journal, you can also liimit the number of retrived results or use pagination. Filtering by publication date is also possible. Example query:

```bash
curl -X 'GET' \
  'http://localhost:8000/journal/Nature?results_num=100&p=2&min_date=2020%2F01%2F01&max_date=2022%2F01%2F01' \
  -H 'accept: application/json'
```

### Customization

You can change the parameters of the API like `ENTREZ_EMAIL` in the `configs.py` file.
