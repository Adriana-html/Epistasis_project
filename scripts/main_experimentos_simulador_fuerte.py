import subprocess
import re
import csv
import os


PACIENTES = 7500
SNPS_FIJOS = 30000 
CORES_A_PROBAR = [1, 2, 4, 8, 16, 24, 32, 64]
REPETICIONES_PARALELO = 3 
REPETICIONES_SERIAL = 1
CARPETA_DATA = os.path.abspath("data")

def ejecutar_y_validar(comando, nombre_prueba, num_iter):
    print(f"\n>>> {nombre_prueba} - Repeticion {num_iter+1}:")
    try:
        process = subprocess.run(comando, shell=True, capture_output=True, text=True)
        if process.returncode != 0:
            print(f"ERROR: {process.stderr}")
            return None, None

        salida = process.stdout.strip()
        print(salida)

        t_match = re.search(r"Tiempo:\s+([\d.]+)", salida)
        tiempo = float(t_match.group(1)) if t_match else None

        s_match = re.search(r"(?:Max Score|Score):\s+([\d.]+)", salida)
        score = float(s_match.group(1)) if s_match else None
        
        return tiempo, score
    except Exception as e:
        print(f"    Error inesperado: {e}")
        return None, None

def main():
    os.makedirs("results", exist_ok=True)
    resultados_strong = []
    tiempo_base = 0.0
    score_referencia = None

    print(f"DATASET TERCEROS (Fijo {SNPS_FIJOS} SNPs) ===")

    
    os.system("rm -f data/genotypes.bin data/phenotypes.bin")
    os.system(f"ln -sf {CARPETA_DATA}/genotypes_simulador_{SNPS_FIJOS}.bin data/genotypes.bin")
    os.system(f"ln -sf {CARPETA_DATA}/phenotypes_simulador_{SNPS_FIJOS}.bin data/phenotypes.bin")

    for cores in CORES_A_PROBAR:
        tiempos_corrida = []
        actual_reps = REPETICIONES_SERIAL if cores == 1 else REPETICIONES_PARALELO
        nombre_prueba = f"Strong_{cores}Cores"

        if cores == 1:
            comando = f"./serial_v3 {PACIENTES} {SNPS_FIJOS}"
        else:
            comando = f"mpirun -np {cores} --bind-to core ./mpi_v3 {PACIENTES} {SNPS_FIJOS}"

        for i in range(actual_reps):
            t, s = ejecutar_y_validar(comando, nombre_prueba, i)
            if t is not None:
                tiempos_corrida.append(t)
                
                # Validación de consistencia del resultado
                if cores == 1 and i == 0:
                    score_referencia = s
                elif score_referencia is not None and s is not None:
                    if abs(s - score_referencia) > 1e-6:
                        print(f"AVISO: Score {s} difiere de referencia {score_referencia}")

        if tiempos_corrida:
            t_promedio = sum(tiempos_corrida) / len(tiempos_corrida)
            if cores == 1:
                tiempo_base = t_promedio
            
            # Speedup = T(1) / T(n)
            speedup = tiempo_base / t_promedio
            # Eficiencia = Speedup / n
            eficiencia = speedup / cores

            resultados_strong.append({
                "Cores": cores,
                "Tiempo_Promedio": round(t_promedio, 6),
                "Speedup": round(speedup, 2),
                "Eficiencia": round(eficiencia, 2)
            })

    # Guardar en CSV
    csv_path = "results/analisis_strong_scaling_terceros.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Cores", "Tiempo_Promedio", "Speedup", "Eficiencia"])
        writer.writeheader()
        writer.writerows(resultados_strong)
    
    print(f"\nEXPERIMENTACIÓN FINALIZADA: {csv_path}")

if __name__ == "__main__":
    main()