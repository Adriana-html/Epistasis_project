import subprocess
import os
import sys

def run_command(cmd_list):
    """Ejecuta un comando (lista de strings) y maneja la salida."""
    print(f">> Ejecutando: {' '.join(cmd_list)}")
    # Usamos shell=False por seguridad al pasar listas
    result = subprocess.run(cmd_list, capture_output=True, text=True)
    
    if result.returncode == 0:
        print(result.stdout)
        return True
    else:
        print("ERROR:")
        print(result.stderr)
        return False

def main():
    # 1. Definir los parámetros que quieres para tu prueba
    # IMPORTANTE: n_individuals debe ser divisible entre los procesos MPI que uses para generar
    N_INDIV = 500  
    N_SNPS = 1500

    # 2. Crear carpeta data si no existe
    if not os.path.exists("data"):
        os.makedirs("data")

    # 3. Limpiar archivos previos (tu generador da error si ya existen)
    for f in ["data/genotypes.csv", "data/phenotypes.csv", "data/genotype.bin", "data/phenotype.bin"]:
        if os.path.exists(f):
            os.remove(f)

    print(f"=== PREPARANDO EXPERIMENTO: {N_INDIV} Indiv, {N_SNPS} SNPs ===")

    # EJECUTAR GENERADOR (MPI)
    # Aquí pasamos los argumentos --n_individuals y --n_snps
    gen_cmd = [
        "mpirun", "-np", "2", 
        sys.executable, "src/generate_data_mpi.py", 
        "--n_individuals", str(N_INDIV), 
        "--n_snps", str(N_SNPS)
    ]
    if not run_command(gen_cmd): return
    # CONVERTIR A BINARIO
    conv_cmd = [sys.executable, "scripts/convertir_archivos.py"]
    if not run_command(conv_cmd): return

    print("=== Datos listos ===")

if __name__ == "__main__":
    main()