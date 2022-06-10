import pandas as pd
import ncls

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def split_pandas_df(data: pd.DataFrame) -> tuple:
    # add index and group by
    data["_index"] = range(0, data.shape[0])
    groups = data.groupby("seqname")

    # generate NCLS indexes for each seqname
    indexes = {}
    for group, rows in groups:
        indexes[group] = ncls.NCLS(
            rows.start.values, rows.end.values, rows._index.values
        )

    ranges = pd.DataFrame(
        {"seqnames": data.seqname, "starts": data.start, "ends": data.end}
    )

    metadata = data.drop(["seqname", "start", "end", "_index"], axis=1)
    if metadata.empty:
        metadata = None

    return (indexes, data["_index"].tolist(), ranges, metadata)
