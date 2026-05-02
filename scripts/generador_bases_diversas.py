import subprocess
import os
import sys

def run_command(cmd_list):
    
    print(f">> Ejecutando: {' '.join(cmd_list)}")
    result = subprocess.run(cmd_list, capture_output=True, text=True)
    if result.returncode == 0:
        return True
    else:
        print(f"ERROR en comando: {result.stderr}")
        return False

def preparar_escenario(n_indiv, n_snps, carpeta_final):
    
    print(f"\n--- Generar base de datos: {n_snps} SNPs ---")
    gen_cmd = [
        "mpirun", "-np", "15", 
        sys.executable, "src/generate_data_mpi.py", 
        "--n_individuals", str(n_indiv), 
        "--n_snps", str(n_snps)
    ]
    if not run_command(gen_cmd): return False

    conv_cmd = [sys.executable, "scripts/convertir_archivos.py"]
    if not run_command(conv_cmd): return False

    os.makedirs(carpeta_final, exist_ok=True)
    os.rename("data/genotypes.bin", f"{carpeta_final}/genotypes_{n_snps}.bin")
    os.rename("data/phenotypes.bin", f"{carpeta_final}/phenotypes_{n_snps}.bin")
    
    if os.path.exists("data/genotypes.csv"): os.remove("data/genotypes.csv")
    if os.path.exists("data/phenotypes.csv"): os.remove("data/phenotypes.csv")

    print(f"Base de datos {n_snps} listo en {carpeta_final}")
    return True

def main():
    N_INDIV = 7500
    
    tamanos_snps = [3750, 5303, 7500, 10606, 15000, 21213, 30000]
    
    CARPETA_DATASET = "data/dataset_experimentacion"
    
    if not os.path.exists("data"): os.makedirs("data")

    print(f"---- GENERACIÓN DE BASES DE DATOS ----")
    
    for s in tamanos_snps:
        exito = preparar_escenario(N_INDIV, s, CARPETA_DATASET)
        if not exito:
            print(f"Terminando prematuramente por error en tamaño {s}")
            break

    print(f"\n---- PROCESO FINALIZADO: Todos los binarios están en {CARPETA_DATASET} ----")

if __name__ == "__main__":
    main()