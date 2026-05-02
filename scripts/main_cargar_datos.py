import subprocess
import os
import sys

def run_command(cmd_list):
    print(f">> Ejecutando: {' '.join(cmd_list)}")
    result = subprocess.run(cmd_list, capture_output=True, text=True)
    
    if result.returncode == 0:
        print(result.stdout)
        return True
    else:
        print("ERROR:")
        print(result.stderr)
        return False

def main():
   
    # IMPORTANTE: n_individuals debe ser divisible entre los procesos MPI que uses para generar
    N_INDIV = 500  
    N_SNPS = 1500

    # Crear carpeta data si no existe
    if not os.path.exists("data"):
        os.makedirs("data")

    
    for f in ["data/genotypes.csv", "data/phenotypes.csv", "data/genotype.bin", "data/phenotype.bin"]:
        if os.path.exists(f):
            os.remove(f)

    print(f"=== EXPERIMENTO: {N_INDIV} Indiv, {N_SNPS} SNPs ===")

    gen_cmd = [
        "mpirun", "-np", "2", 
        sys.executable, "src/generate_data_mpi.py", 
        "--n_individuals", str(N_INDIV), 
        "--n_snps", str(N_SNPS)
    ]
    if not run_command(gen_cmd): return
    conv_cmd = [sys.executable, "scripts/convertir_archivos.py"]
    if not run_command(conv_cmd): return

    print("=== Datos listos ===")

if __name__ == "__main__":
    main()