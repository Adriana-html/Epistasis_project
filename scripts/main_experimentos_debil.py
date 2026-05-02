import subprocess
import re
import csv
import os
import math

PACIENTES = 7500
CONFIGURACIONES = {
    1: 3750,
    2: 5303,
    4: 7500,
    8: 10606,
    16: 15000,
    32: 21213,
    64: 30000
}
REPETICIONES_PARALELO = 3 
REPETICIONES_SERIAL = 1
# ruta absoluta para evitar que el link se rompa
CARPETA_DATA = os.path.abspath("data/dataset_experimentacion")

def ejecutar_y_validar(comando, nombre_prueba, num_iter):
    print(f"\n>>> {nombre_prueba} - Repeticion {num_iter+1}:")
    try:
        process = subprocess.run(comando, shell=True, capture_output=True, text=True)
        if process.returncode != 0:
            print(f"    ERROR: {process.stderr}")
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
    # Asegurar que las carpetas existan
    os.makedirs("results", exist_ok=True)
    os.makedirs("data", exist_ok=True) 

    resultados_weak = []
    tiempo_base = 0.0
    score_referencia = None

    print(f"Base: 3,750 SNPs | Máx: 30,000 SNPs | Pacientes: {PACIENTES}")

    cores_ordenados = sorted(CONFIGURACIONES.keys())

    for cores in cores_ordenados:
        snps = CONFIGURACIONES[cores]
        tiempos_corrida = []
        actual_reps = REPETICIONES_SERIAL if cores == 1 else REPETICIONES_PARALELO
        
        nombre_prueba = f"Weak_{cores}Cores_{snps}SNPs"
        
       
        # Borrar enlaces anteriores
        os.system("rm -f data/genotypes.bin data/phenotypes.bin")
        
        # rutas absolutas (Ruta origen -> Nombre destino)
        os.system(f"ln -s {CARPETA_DATA}/genotypes_{snps}.bin data/genotypes.bin")
        os.system(f"ln -s {CARPETA_DATA}/phenotypes_{snps}.bin data/phenotypes.bin")

        if cores == 1:
            comando = f"./serial_epistasis {PACIENTES} {snps}"
        else:
            comando = f"mpirun -np {cores} --bind-to core ./mpi_epistasis {PACIENTES} {snps}"

        for i in range(actual_reps):
            t, s = ejecutar_y_validar(comando, nombre_prueba, i)
            if t is not None:
                tiempos_corrida.append(t)
                if cores == 1 and i == 0:
                    score_referencia = s
                elif score_referencia is not None and s is not None:
                    if abs(s - score_referencia) > 1e-5:
                        print(f"ALERTA: Score {s} difiere de referencia {score_referencia}")
                    else:
                        print(f"Score validado.")

        if tiempos_corrida:
            t_promedio = sum(tiempos_corrida) / len(tiempos_corrida)
            if cores == 1: tiempo_base = t_promedio
            eficiencia = (tiempo_base / t_promedio) * 100
            resultados_weak.append({
                "Configuracion": nombre_prueba,
                "Cores": cores,
                "SNPs": snps,
                "Tiempo_Promedio": round(t_promedio, 6),
                "Eficiencia_%": round(eficiencia, 2)
            })

    csv_path = "results/analisis_weak_scaling_int8.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Configuracion", "Cores", "SNPs", "Tiempo_Promedio", "Eficiencia_%"])
        writer.writeheader()
        writer.writerows(resultados_weak)

    print(f"\n\n----- EXPERIMENTACIÓN FINALIZADA: {csv_path} ----")

if __name__ == "__main__":
    main()