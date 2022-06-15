import pandas as pd
import ncls

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def split_pandas_df(data: pd.DataFrame) -> tuple:
    # add index and group by
    data["_index"] = range(0, data.shape[0])
    groups = data.groupby("seqnames")

    # generate NCLS indexes for each seqname
    indexes = {}
    for group, rows in groups:
        indexes[group] = ncls.NCLS(
            rows.starts.astype(int).values,
            rows.ends.astype(int).values,
            rows._index.astype(int).values,
        )

    ranges = pd.DataFrame(
        {"seqnames": data.seqnames, "starts": data.starts, "ends": data.ends}
    )

    metadata = data.drop(["seqnames", "starts", "ends", "_index"], axis=1)
    if metadata.empty:
        metadata = None

    return (indexes, data["_index"].tolist(), ranges, metadata)
