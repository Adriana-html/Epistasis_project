import subprocess
import re
import csv
import os

# --- Configuración ---
PACIENTES = 500
SNPS = 1500
CORES_A_PROBAR = [1, 2, 4] 
REPETICIONES = 1 

def ejecutar_y_validar(comando, nombre_prueba, num_iter):
    """Ejecuta el comando y extrae tiempo y score."""
    print(f"\n>>> {nombre_prueba} - Repeticion {num_iter+1}:")
    
    try:
        process = subprocess.run(comando, shell=True, capture_output=True, text=True)
        if process.returncode != 0:
            print(f"    ERROR: {process.stderr}")
            return None, None

        salida = process.stdout.strip()
        print(salida)

        # Extraer Tiempo
        t_match = re.search(r"Tiempo:\s+([\d.]+)", salida)
        tiempo = float(t_match.group(1)) if t_match else None

        # Extraer Max Score (Pearson)
        s_match = re.search(r"(?:Max Score|Score):\s+([\d.]+)", salida)
        score = float(s_match.group(1)) if s_match else None
        
        return tiempo, score
    except Exception as e:
        print(f"    Error inesperado: {e}")
        return None, None

def main():
    os.makedirs("results", exist_ok=True)
    resultados_finales = []
    tiempo_serial_base = 0.0
    score_serial_referencia = None

    print(f"=== INICIANDO EXPERIMENTACION (SNPs: {SNPS}, Pacientes: {PACIENTES}) ===")

    for cores in CORES_A_PROBAR:
        tiempos_corrida = []
        
        if cores == 1:
            nombre_prueba = "SERIAL"
            comando = f"./serial_epistasis {PACIENTES} {SNPS}"
        else:
            nombre_prueba = f"MPI ({cores} cores)"
            comando = f"mpirun -np {cores} ./mpi_epistasis {PACIENTES} {SNPS}"

        for i in range(REPETICIONES):
            t, s = ejecutar_y_validar(comando, nombre_prueba, i)
            
            if t is not None:
                tiempos_corrida.append(t)
                
                # Validación de Score contra el Serial
                if cores == 1 and i == 0:
                    score_serial_referencia = s
                elif score_serial_referencia is not None and s is not None:
                    # Usamos math.isclose o una tolerancia pequeña por precisión de punto flotante
                    if abs(s - score_serial_referencia) > 1e-5:
                        print(f"ALERTA: El score ({s}) difiere del Serial ({score_serial_referencia})")
                    else:
                        print(f"Score validado correctamente.")

        if tiempos_corrida:
            t_promedio = sum(tiempos_corrida) / len(tiempos_corrida)
            if cores == 1: tiempo_serial_base = t_promedio
            
            speedup = tiempo_serial_base / t_promedio
            eficiencia = speedup / cores

            resultados_finales.append({
                "Configuracion": nombre_prueba,
                "Cores": cores,
                "Tiempo_Promedio": round(t_promedio, 6),
                "Speedup": round(speedup, 2),
                "Eficiencia": round(eficiencia, 2)
            })

    # Guardar en CSV
    csv_path = "results/analisis_rendimiento.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Configuracion", "Cores", "Tiempo_Promedio", "Speedup", "Eficiencia"])
        writer.writeheader()
        writer.writerows(resultados_finales)

    print(f"\n\n=== EXPERIMENTACION FINALIZADA: {csv_path} ===")

if __name__ == "__main__":
    main()