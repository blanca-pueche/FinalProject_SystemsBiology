import requests
from Bio import Entrez

# Configurar correo para NCBI Entrez
Entrez.email = "tu_email@example.com"

# Función para buscar genes en NCBI
def get_ncbi_genes(mesh_term):
    search_handle = Entrez.esearch(db="gene", term=mesh_term, retmax=10)
    record = Entrez.read(search_handle)
    search_handle.close()
    gene_ids = record["IdList"]
    return gene_ids

# Función para buscar genes en DisGeNET
def get_disgenet_genes(mesh_term):
    url = f"https://www.disgenet.org/api/gda/disease/{mesh_term}?format=json"
    headers = {"Authorization": "Bearer TU_API_KEY"}  # Reemplazar con tu clave de API
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        data = response.json()
        return [entry["gene_symbol"] for entry in data]  # Extrae los genes
    return []

# Función para buscar genes en OMIM
def get_omim_genes(mesh_term):
    url = f"https://api.omim.org/api/entry/search?search={mesh_term}&include=geneMap&apiKey=TU_API_KEY"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        genes = []
        for entry in data.get("omim", {}).get("entryList", []):
            if "geneMap" in entry["entry"]:
                genes.append(entry["entry"]["geneMap"]["geneSymbols"])
        return genes
    return []

# Ejecutar búsqueda
mesh_term = "Diabetes Mellitus"  # Reemplazar con el término MeSH deseado
ncbi_genes = get_ncbi_genes(mesh_term)
disgenet_genes = get_disgenet_genes(mesh_term)
omim_genes = get_omim_genes(mesh_term)

# Mostrar resultados
print("NCBI Genes:", ncbi_genes)
print("DisGeNET Genes:", disgenet_genes)
print("OMIM Genes:", omim_genes)
