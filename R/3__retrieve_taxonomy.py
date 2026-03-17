import pandas as pd
from ete3 import NCBITaxa
ncbi = NCBITaxa()

names = pd.read_csv("data//processed//common_nodes.csv", index_col=0)
names.columns = ["common_nodes"]
names = names["common_nodes"]

sp_id = []
for n in names:
    try:
        idx = ncbi.get_name_translator([n])[n][0]
    except:
        idx = -1
    sp_id.append(idx)

results = pd.DataFrame({"sp_name": names, "sp_id": sp_id})
results.drop(results.index[results["sp_id"] == -1], inplace=True)

sp_id = results["sp_id"]
gen_id = [ncbi.get_lineage(idx)[-2] for idx in sp_id]
fam_id = [ncbi.get_lineage(idx)[-3] for idx in sp_id]
ord_id = [ncbi.get_lineage(idx)[-4] for idx in sp_id]
cls_id = [ncbi.get_lineage(idx)[-5] for idx in sp_id]
gen_name = [ncbi.get_taxid_translator([idx])[idx] for idx in gen_id]
fam_name = [ncbi.get_taxid_translator([idx])[idx] for idx in fam_id]
ord_name = [ncbi.get_taxid_translator([idx])[idx] for idx in ord_id]
cls_name = [ncbi.get_taxid_translator([idx])[idx] for idx in cls_id]

results = pd.DataFrame({
    "sp_name": results["sp_name"], "sp_id": sp_id,
    "gen_id": gen_id, "fam_id": fam_id, 
    "ord_id": ord_id, "cl_id": cls_id,
    "gen_name": gen_name, "fam_name": fam_name, 
    "ord_name": ord_name, "cl_name": cls_name
})
results.to_csv("sp_name_idx.csv")
