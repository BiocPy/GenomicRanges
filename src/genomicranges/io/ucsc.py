from typing import Literal

from .gtf import parse_gtf

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def access_gtf_ucsc(
    genome: str,
    type: Literal["refGene", "ensGene", "knownGene", "ncbiRefSeq"] = "refGene",
) -> str:
    """Generate a path to a genome gtf file from UCSC,
    e.g. for `hg19 genome <http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/genes/>`_.

    Args:
        genome:
            Genome shortcode; e.g. hg19, hg38, mm10 etc.

        type:
            Defaults to "refGene".

    Raises:
        Exception, ValueError:
            When ``type`` does not match with a valid input.

    Returns:
        The URI to the file.
    """
    base_path = f"http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/genes/"

    if type not in ["refGene", "ensGene", "knownGene", "ncbiRefSeq"]:
        raise ValueError(
            f"type must be one of refGene, ensGene, knownGene or ncbiRefSeq, provided {type}"
        )

    full_path = f"{base_path}/{genome}.{type}.gtf.gz"

    return full_path


def read_ucsc(
    genome: str,
    type: Literal["refGene", "ensGene", "knownGene", "ncbiRefSeq"] = "refGene",
) -> "GenomicRanges":
    """Load a genome annotation from UCSC as :py:class:`~genomicranges.GenomicRanges.GenomicRanges`.

    Args:
        genome:
            Genome shortcode; e.g. hg19, hg38, mm10 etc.

        type:
            Defaults to "refGene".

    Returns:
        The gene model from UCSC.
    """
    path = access_gtf_ucsc(genome, type=type)
    compressed = True
    data = parse_gtf(path, compressed=compressed)

    from ..GenomicRanges import GenomicRanges

    return GenomicRanges.from_pandas(data)
