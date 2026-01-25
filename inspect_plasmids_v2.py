from Bio import Entrez, SeqIO

Entrez.email = "test@example.com"

def inspect_record(accession):
    print(f"Fetching {accession}...")
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        print(f"Length: {len(record.seq)}")
        found = False
        for feature in record.features:
            if feature.type == "CDS":
                qualifiers = feature.qualifiers
                gene = qualifiers.get("gene", [""])[0]
                product = qualifiers.get("product", [""])[0]

                if "fosA" in gene or "fosA" in product or "fosfomycin" in product:
                     print(f"Found putative fos gene:")
                     print(f"  Gene: {gene}")
                     print(f"  Product: {product}")
                     print(f"  Location: {feature.location}")
                     found = True

        if not found:
            print("No explicit 'fosA' feature found. Listing all CDS products...")
            for feature in record.features:
                 if feature.type == "CDS":
                    print(f"  {feature.qualifiers.get('product', [''])[0]} ({feature.location})")
    except Exception as e:
        print(f"Error fetching {accession}: {e}")

inspect_record("KP324830.1")
inspect_record("KY270852.1")
