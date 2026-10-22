#!/usr/bin/env python3
"""Extract all Mycobacterium tuberculosis runs from the AWS SRA metadata (anonymous S3).

Usage: fetch_sra_metadata.py <out.parquet>
Source: s3://sra-pub-metadata-us-east-1/sra/metadata/ (parquet, refreshed by NCBI).
Needs: s3fs, pyarrow.
"""
import os, sys
import pyarrow as pa, pyarrow.compute as pc, pyarrow.parquet as pq, s3fs

kw = {"proxies": {"https": os.environ["HTTPS_PROXY"]}} if "HTTPS_PROXY" in os.environ else {}
fs = s3fs.S3FileSystem(anon=True, client_kwargs={"region_name": "us-east-1"}, config_kwargs=kw)
cols = ["acc", "sample_acc", "biosample", "sample_name", "organism", "bioproject", "sra_study", "mbases",
        "avgspotlen", "librarylayout", "platform", "instrument", "assay_type", "geo_loc_name_country_calc",
        "sample_name_sam", "center_name", "releasedate"]
out = []
for f in fs.ls("sra-pub-metadata-us-east-1/sra/metadata/"):
    pf = pq.ParquetFile(fs.open(f, block_size=8 * 2**20))
    for rg in range(pf.num_row_groups):
        t = pf.read_row_group(rg, columns=cols)
        out.append(t.filter(pc.fill_null(pc.match_substring(t["organism"], "Mycobacterium tuberculosis"), False)))
t = pa.concat_tables(out)
t = t.set_column(t.schema.get_field_index("sample_name_sam"), "sample_name_sam", pc.binary_join(t["sample_name_sam"], ";"))
pq.write_table(t, sys.argv[1])
