from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import time
STOP_AT_RECORD = 5000
class NCBIRetriever:
    def __init__(self, email, api_key):
        Entrez.email = email #zmienilam na przypisanie od razu
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'
    def search_taxid(self, taxid):
        try:
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            self.organism_name = records[0]["ScientificName"]
            print(f"Organism: {self.organism_name} (TaxID: {taxid})")

            search_term = f"txid{taxid}[Organism]"
            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            results = Entrez.read(handle)

            self.count = int(results["Count"])

            if self.count == 0:
                print(f"No records found for {self.organism_name}")
                return None
            print(f"Found {self.count} records")

            self.webenv = results["WebEnv"]
            self.query_key = results["QueryKey"]
            print(f"Found {self.count} records.")
            return self.count
        except Exception as e:
            print(f"Search error: {e}")
            return 0

    def fetch_filtered_records(self, min_len, max_len):
        records = []
        fetched = 0
        batch_size =500

        for start in range(0, min(self.count, STOP_AT_RECORD), batch_size):
            print(f"Fetching in progress :{start}")

            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=self.webenv,
                query_key=self.query_key
            )

            for record in SeqIO.parse(handle, "gb"):
                seq_len = len(record.seq)
                if min_len <= seq_len <= max_len:
                    records.append({
                        "accession": record.id,
                        "length": seq_len,
                        "description": record.description
                    })
            handle.close()
            if len(records) >= 500:
                break
            time.sleep(0.1) # z kluczem 10 zapytan
        return records

def generate_csv(data, filename):
    df=pd.DataFrame(data)
    df.to_csv(filename, index=False)

def generate_plot(data, filename):
    df= pd.DataFrame(data)
    df_sorted = df.sort_values(by="length", ascending=False)
    plt.figure(figsize=(10, 6))
    plt.plot(df_sorted["accession"], df_sorted["length"], marker='o')
    plt.xticks(rotation=90, fontsize=6)
    plt.xlabel("Accession")
    plt.ylabel("Sequence Length")
    plt.title("Sequence Lengths (sorted)")
    plt.tight_layout()
    plt.savefig(filename)

def main():
    email = input("Enter your email address for NCBI: ")
    api_key = input("Enter your NCBI API key: ")
    taxid = input("Enter taxonomic ID (taxid) of the organism: ")
    min_len = int(input("Enter minimum sequence length: "))
    max_len = int(input("Enter maximum sequence length: "))

    retriever = NCBIRetriever(email, api_key)
    if retriever.search_taxid(taxid) == 0:
        print("No records found")
        return

    records = retriever.fetch_filtered_records(min_len, max_len)
    if not records:
        print("No sequences matched given restrictions")
        return

    csv_file = f"{taxid}_filtered.csv"
    plot_file = f"{taxid}_plot.png"
    generate_csv(records, csv_file)
    generate_plot(records, plot_file)

if __name__ == "__main__":
    main()