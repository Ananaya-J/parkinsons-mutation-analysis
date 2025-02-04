import pandas as pd
data = pd.read_csv("annotated.hg19_multianno.txt", sep="\t")
pd_genes = ["SNCA", "LRRK2", "PARK2", "PINK1", "DJ1"]
pd_mutations = data[data["Gene.refGene"].isin(pd_genes)]
pd_mutations.to_csv("pd_mutations.csv", index=False)
