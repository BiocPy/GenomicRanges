import logging

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def access_gtf_ucsc(genome: str, type: str = "refGene") -> str:
    """Generate a path to a genome gtf file from UCSC e.g. (http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/genes/)

    Args:
        genome (str): genome shortcode; e.g. hg19, hg38, mm10 etc
        type (str): One of refGene, ensGene, knownGene or ncbiRefSeq

    Raises:
        Exception: TypeError, when `type` does not match with a valid input

    Returns:
        str: returns the path to the file
    """
    base_path = (
        f"http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/genes/"
    )

    if type not in ["refGene", "ensGene", "knownGene", "ncbiRefSeq"]:
        logging.error(f"Provided value for type: {type}")
        logging.error(
            "type must be one of refGene, ensGene, knownGene or ncbiRefSeq"
        )
        raise Exception("TypeError")

    full_path = f"{base_path}/{genome}.{type}.gtf.gz"

    return full_path
