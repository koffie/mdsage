import json
from pathlib import Path
from sage.all import cached_function

SRC_DIR = Path(__file__).parent.resolve()
DATA_DIR = SRC_DIR.joinpath("data_files")


@cached_function
def load_json_data(filename):
    filename = DATA_DIR.joinpath(filename)
    with open(filename, "r") as file:
        data = file.read()
    return json.loads(data)
