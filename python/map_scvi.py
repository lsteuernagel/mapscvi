## run anndata script

# Import relevant modules
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
import torch
import scvi
from os import listdir
from os.path import isfile, join
print(scvi.__version__)

# get parameters
query_path = sys.argv[1]
model_path = sys.argv[2]
latent_dim_path = sys.argv[3]
max_epochs = sys.argv[4]

# load query anndata
#print("load query anndata")
print("using anndata stored at: "+query_path)
adata_query = sc.read_h5ad(query_path)
#ensure there are no bytestrings which can crash scvi
str_df = adata_query.obs
str_df = str_df.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)
adata_query.obs = str_df

#print("online update")
model_path = str(model_path)
print("using model stored at: "+model_path)
# online update:
# https://docs.scvi-tools.org/en/stable/user_guide/notebooks/scarches_scvi_tools.html
vae_q = scvi.model.SCVI.load_query_data(
    adata_query,
    model_path, # get model directly from disk!
    inplace_subset_query_vars=False
)
# train --> how many epochs ?
print("Training with disabled error bar.")
vae_q.train(max_epochs=int(max_epochs), plan_kwargs=dict(weight_decay=0.0),progress_bar_refresh_rate=0) # use same epochs and weight_decay = 0
# results
adata_query.obsm["X_scVI"] = vae_q.get_latent_representation() # get laten dim

print("save")
#save in latent_dim_path
output = pd.DataFrame(adata_query.obsm["X_scVI"])
output = output.set_index(adata_query.obs_names)
output2 = output.set_axis(["scVI_" + str(s) for s in output.axes[1].to_list()], axis=1, inplace=False)
output2.to_csv(latent_dim_path, sep='\t',index=True)

print("finalized")
