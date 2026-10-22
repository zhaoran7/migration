"""Export public atlas coordinates without loading its expression matrix."""
import csv
import gzip
from pathlib import Path
import sys
import anndata as ad

file, group, output = sys.argv[1:]
x = ad.read_h5ad(file, backed="r")
key = next(k for k in ("X_umap", "X_UMAP", "UMAP") if k in x.obsm)
xy = x.obsm[key]
labels = x.obs[group].astype(str)
Path(output).parent.mkdir(parents=True, exist_ok=True)
with gzip.open(output, "wt", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["cell", "umap1", "umap2", "label"])
    writer.writerows(zip(x.obs_names, xy[:, 0], xy[:, 1], labels))
x.file.close()
