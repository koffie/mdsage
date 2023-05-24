import json
from pathlib import Path
from sage.all import cached_function

SRC_DIR = Path(__file__).parent.resolve()
DATA_DIR = SRC_DIR.joinpath("data_files")


@cached_function
def load_json_data(filename):
    """
    loads json data from the data_files directory

    TESTS::

        sage: from mdsage.utilities import load_json_data
        sage: data = load_json_data("small_class_number_fundamental_discriminant.json")
        sage: data["1"]
        [-3, -4, -7, -8, -11, -19, -43, -67, -163]

    """
    filename = DATA_DIR.joinpath(filename)
    with open(filename, "r") as file:
        data = file.read()
    return json.loads(data)
