import json
import numpy as np
import os
import sys

def convert_sim_to_bin(json_path, out_folder="data"):
    """
    Convierte la salida del simulador externo (JSON) al formato binario.
    """
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    if not os.path.exists(json_path):
        print(f"Error: No se encuentra el archivo {json_path}")
        return

    print(f"Cargando datos desde: {json_path}")
    
    with open(json_path, 'r') as f:
        try:
            data = json.load(f)
        except json.JSONDecodeError:
            print("Error: El archivo no es un JSON válido.")
            return

    # Extraer metadatos
    n_snps = data.get('num_snps')
    n_inds = data.get('num_inds')
    
    print(f"Dimensiones: {n_snps} SNPs x {n_inds} Individuos")

    genotypes = np.array(data['genotype'], dtype=np.int8)
    
    if genotypes.shape != (n_snps, n_inds):
        print(f"Aviso: Ajustando forma de genotipos de {genotypes.shape} a ({n_snps}, {n_inds})")
        genotypes = genotypes.reshape((n_snps, n_inds))

    phenotypes = np.array(data['phenotype'], dtype=np.float64)

    
    geno_path = os.path.join(out_folder, "genotypes_simulador.bin")
    pheno_path = os.path.join(out_folder, "phenotypes_simulador.bin")
    
    genotypes.tofile(geno_path)
    phenotypes.tofile(pheno_path)

    print("\nCONVERSIÓN")
    print(f"Archivos generados en '{out_folder}':")
    print(f"   - genotypes.bin  ({os.path.getsize(geno_path)} bytes)")
    print(f"   - phenotypes.bin ({os.path.getsize(pheno_path)} bytes)")
    print("-" * 40)

if __name__ == "__main__":
    archivo_json = sys.argv[1] if len(sys.argv) > 1 else "sim/0_1_ASW.json"
    convert_sim_to_bin(archivo_json)