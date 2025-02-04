import pandas as pd
import requests

# Load the pd_mutations.csv file
input_file = "pd_mutations.csv"
output_file = "pd_mutations_annotated.csv"

# Read the input file
df = pd.read_csv(input_file)

# Ensure the input file has the required columns
required_columns = ["Chr", "Start", "Ref", "Alt"]
if not all(col in df.columns for col in required_columns):
    raise ValueError(f"Input file must contain the following columns: {required_columns}")

# Function to query ClinVar
def query_clinvar(chrom, pos, ref, alt):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&term={chrom}[Chromosome]+AND+{pos}[Position]+AND+{ref}[Ref]+AND+{alt}[Alt]&retmode=json"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    return None

# Function to query dbSNP
def query_dbsnp(chrom, pos, ref, alt):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&term={chrom}[Chromosome]+AND+{pos}[Position]+AND+{ref}[Ref]+AND+{alt}[Alt]&retmode=json"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    return None

# Function to query gnomAD
def query_gnomad(chrom, pos, ref, alt):
    url = f"https://gnomad.broadinstitute.org/api/?query=variant(chrom:{chrom},pos:{pos},ref:{ref},alt:{alt})"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    return None

# Annotate each variant
clinvar_data = []
dbsnp_data = []
gnomad_data = []

for index, row in df.iterrows():
    chrom = row["Chr"]
    pos = row["Start"]
    ref = row["Ref"]
    alt = row["Alt"]

    # Query ClinVar
    clinvar_result = query_clinvar(chrom, pos, ref, alt)
    clinvar_data.append(clinvar_result)

    # Query dbSNP
    dbsnp_result = query_dbsnp(chrom, pos, ref, alt)
    dbsnp_data.append(dbsnp_result)

    # Query gnomAD
    gnomad_result = query_gnomad(chrom, pos, ref, alt)
    gnomad_data.append(gnomad_result)

# Add results to the DataFrame
df["ClinVar_Additional"] = clinvar_data
df["dbSNP_Additional"] = dbsnp_data
df["gnomAD_Additional"] = gnomad_data

# Save the annotated DataFrame to a new CSV file
df.to_csv(output_file, index=False)

print(f"Annotated data saved to {output_file}")
