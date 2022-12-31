from ..GenomicRanges import GenomicRanges
from .gtf import parse_gtf
from .pdf import fromPandas

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def access_gtf_ucsc(genome: str, type: str = "refGene") -> str:
    """Generate a path to a genome gtf file from UCSC e.g. (http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/genes/)

    Args:
        genome (str): genome shortcode; e.g. hg19, hg38, mm10 etc
        type (str): One of refGene, ensGene, knownGene or ncbiRefSeq

    Raises:
        Exception: ValueError, when `type` does not match with a valid input

    Returns:
        str: returns the path to the file
    """
    base_path = f"http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/genes/"

    if type not in ["refGene", "ensGene", "knownGene", "ncbiRefSeq"]:
        raise ValueError(
            f"type must be one of refGene, ensGene, knownGene or ncbiRefSeq, provided {type}"
        )

    full_path = f"{base_path}/{genome}.{type}.gtf.gz"

    return full_path


def readUCSC(genome: str, type: str = "refGene") -> GenomicRanges:
    """Load a genome annotation from UCSC as `GenomicRanges`.

    Args:
        genome (str): genome shortcode; e.g. hg19, hg38, mm10 etc
        type (str): One of refGene, ensGene, knownGene or ncbiRefSeq

    Returns:
        GenomicRanges:  a new `GenomicRanges` representing the gene model
    """
    path = access_gtf_ucsc(genome, type=type)
    compressed = True
    data = parse_gtf(path, compressed=compressed)

    return fromPandas(data)
