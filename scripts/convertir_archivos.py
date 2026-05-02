import numpy as np
import pandas as pd
import os

def generar_binario(csv_input, bin_output):
    try:
        
        df = pd.read_csv(csv_input, index_col=0)      
        matriz = df.values.astype(np.int8)
        matriz.tofile(bin_output)
        peso_real = os.path.getsize(bin_output)
        
        print(f"--- Procesamiento de {csv_input} ---")
        print(f"Dimensiones: {matriz.shape}")
        print(f"Archivo generado: {bin_output} ({peso_real / (1024**2):.2f} MB)")
            
    except Exception as e:
        print(f"Error procesando {csv_input}: {e}")

if __name__ == "__main__":
    
    generar_binario('data/genotypes.csv', 'data/genotypes.bin')
    generar_binario('data/phenotypes.csv', 'data/phenotypes.bin')