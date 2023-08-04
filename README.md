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

### Set up

You can change the parameters of the API like `ENTREZ_EMAIL` in the `configs.py` file.