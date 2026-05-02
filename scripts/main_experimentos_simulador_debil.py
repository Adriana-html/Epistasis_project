import subprocess
import re
import csv
import os


PACIENTES = 7500
CONFIGURACIONES = {
    1: 3750, 2: 5303, 4: 7500, 8: 10606, 16: 15000, 32: 21213, 64: 30000
}
REPETICIONES_PARALELO = 3 
REPETICIONES_SERIAL = 1
# Ajusta esta ruta a donde están los archivos
CARPETA_DATA = os.path.abspath("data") 

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
    os.makedirs("results", exist_ok=True)
    resultados_weak = []
    tiempo_base = 0.0
    score_referencia = None

    cores_ordenados = sorted(CONFIGURACIONES.keys())

    for cores in cores_ordenados:
        snps = CONFIGURACIONES[cores]
        tiempos_corrida = []
        actual_reps = REPETICIONES_SERIAL if cores == 1 else REPETICIONES_PARALELO
        nombre_prueba = f"Weak_{cores}Cores_{snps}SNPs"
        
        # GESTIÓN DE ENLACES
        os.system("rm -f data/genotypes.bin data/phenotypes.bin")
        # Asegúrate de que el nombre coincida con tu imagen: genotypes_simulador_{snps}.bin
        os.system(f"ln -sf {CARPETA_DATA}/genotypes_simulador_{snps}.bin data/genotypes.bin")
        os.system(f"ln -sf {CARPETA_DATA}/phenotypes_simulador_{snps}.bin data/phenotypes.bin")

        # CAMBIO DE EJECUTABLE (Debe ser el que soporta double) ---
        if cores == 1:
            comando = f"./serial_v3 {PACIENTES} {snps}"
        else:
            # mpi_v2 (la versión para el simulador de terceros)
            comando = f"mpirun -np {cores} --bind-to core ./mpi_v3 {PACIENTES} {snps}"

        for i in range(actual_reps):
            t, s = ejecutar_y_validar(comando, nombre_prueba, i)
            if t is not None:
                tiempos_corrida.append(t)
                if cores == 1 and i == 0:
                    score_referencia = s
                elif score_referencia is not None and s is not None:
                    if abs(s - score_referencia) > 1e-6:
                        print(f"ALERTA: Score {s} difiere de referencia {score_referencia}")

        if tiempos_corrida:
            t_promedio = sum(tiempos_corrida) / len(tiempos_corrida)
            if cores == 1: tiempo_base = t_promedio
            eficiencia = (tiempo_base / t_promedio) * 100
            resultados_weak.append({
                "Cores": cores, "SNPs": snps,
                "Tiempo_Promedio": round(t_promedio, 6),
                "Eficiencia_%": round(eficiencia, 2)
            })

    # Guardar resultados específicos
    csv_path = "results/analisis_weak_scaling_terceros.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Cores", "SNPs", "Tiempo_Promedio", "Eficiencia_%"])
        writer.writeheader()
        writer.writerows(resultados_weak)

if __name__ == "__main__":
    main()