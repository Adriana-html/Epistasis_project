import numpy as np
import pandas as pd

def generar_binario_limpio(csv_input, bin_output):
    try:
        # index_col=0 le dice a pandas que la columna 'Indiv_X' es el índice (no datos)
        # Esto  deja fuera la primera columna y la primera fila
        df = pd.read_csv(csv_input, index_col=0)
        
        # Convertimos a matriz numérica (solo los 0s, 1s y 2s)
        matriz = df.values.astype(np.int32)
        
        # Guardamos en binario
        matriz.tofile(bin_output)
        
        print(f"Éxito: {csv_input} -> {bin_output}")
        print(f"Dimensiones extraídas: {matriz.shape}")
    except Exception as e:
        print(f"Error procesando {csv_input}: {e}")

# Ejecución
generar_binario_limpio('data/genotypes.csv', 'data/genotypes.bin')
generar_binario_limpio('data/phenotypes.csv', 'data/phenotypes.bin')