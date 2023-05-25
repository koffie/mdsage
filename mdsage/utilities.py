import json
import logging
from pathlib import Path
from sage.all import cached_function, cputime

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


class CPUTimeLogger:
    """
    A class for logging how long different steps of an algorithm take.

    EXAMPLES::

        sage: from mdsage.utilities import CPUTimeLogger
        sage: import logging, time
        sage: logger = CPUTimeLogger("my algorithm")
        sage: # do some work
        sage: i = 0
        sage: for j in srange(10**5):
        ....:     i += j
        sage: logger.log("loop completed", level = logging.WARNING) # random
        0.250791 my algorithm loop completed




    Note that only the cputime is logged, so time.sleep is logged as not taking any time.

        sage: time.sleep(0.1)
        sage: logger.log("sleep completed", level = logging.WARNING) # random
        0.0016680000000000583 my algorithm sleep completed


    If the logging level passed to logger.log is lower than the current loging
    level then nothing is logged:

        sage: logger.log("second sleep completed", level = logging.INFO)
    """

    def __init__(self, name: str, logger=None):
        if logger is None:
            logger = logging.getLogger()

        self.t = cputime()
        self.logger = logger
        self.name = name

    def log(self, message: str, level=logging.INFO):
        self.logger.log(level, f"{cputime(self.t):.6f} {self.name} {message}")
        self.t = cputime()
