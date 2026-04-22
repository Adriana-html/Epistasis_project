#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <algorithm>

using namespace std;

// Actualizado para usar punteros int* para que coincida con el main
double calculate_pearson(const int* snp1, const int* snp2, const int* phenotype, int M) {
    double sum_X = 0, sum_Y = 0, sum_X2 = 0, sum_Y2 = 0, sum_XY = 0;

    for (int i = 0; i < M; ++i) {
        double X = (double)snp1[i] * snp2[i]; 
        double Y = (double)phenotype[i];
        
        sum_X += X;
        sum_Y += Y;
        sum_X2 += X * X;
        sum_Y2 += Y * Y;
        sum_XY += X * Y;
    }

    double numerator = ((double)M * sum_XY) - (sum_X * sum_Y);
    double denominator_sq = ((double)M * sum_X2 - (sum_X * sum_X)) * ((double)M * sum_Y2 - (sum_Y * sum_Y));
    
    if (denominator_sq <= 0.0) return 0.0;
    return numerator / sqrt(denominator_sq);
}

int main(int argc, const char **argv) { 
    if (argc < 3) {
        cout << "Uso: " << argv[0] << " <num_pacientes> <num_snps>" << endl;
        return 1;
    }
    const int num_patients = atoi(argv[1]);
    const int num_snps = atoi(argv[2]);

    // 1. Carga de Fenotipos
    vector<int> phenotypes(num_patients);
    // Nota: Asegúrate de que el nombre coincida con tu archivo (phenotype.bin vs phenotypes.bin)
    FILE* f_pheno = fopen("data/phenotypes.bin", "rb"); 
    if (!f_pheno) {
        cerr << "Error: No se pudo abrir phenotype.bin" << endl;
        return 1;
    }
    if (fread(phenotypes.data(), sizeof(int), num_patients, f_pheno) != (size_t)num_patients) {
            cerr << "Error leyendo phenotypes.bin" << endl;
    }
    fclose(f_pheno);

    // 2. Carga de Genotipos
    vector<int> snps_flat((long long)num_snps * num_patients);
    FILE* f_geno = fopen("data/genotypes.bin", "rb");
    if (!f_geno) {
        cerr << "Error: No se pudo abrir genotypes.bin" << endl;
        return 1;
    }
    size_t total_elements = (size_t)num_snps * num_patients;
    if (fread(snps_flat.data(), sizeof(int), total_elements, f_geno) != total_elements) {
        cerr << "Error leyendo genotypes.bin" << endl;
    }
    fclose(f_geno);

    cout << "Datos cargados correctamente. Iniciando calculo serial..." << endl;
    
    double global_max = -1.0;
    int best_i = -1, best_j = -1;

    // 3. Medicion de tiempo y calculo
    auto start_time = chrono::high_resolution_clock::now();
    
    for (int i = 0; i < num_snps; ++i) {
        for (int j = i + 1; j < num_snps; ++j) {
            // Pasamos la direccion de memoria del inicio de cada SNP
            double r = calculate_pearson(&snps_flat[i * num_patients], 
                                         &snps_flat[j * num_patients], 
                                         phenotypes.data(), num_patients);
            
            if (abs(r) > global_max) {
                global_max = abs(r);
                best_i = i;
                best_j = j;
            }
        }
    }
    
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end_time - start_time;

    // 4. Escritura de resultados
    ofstream out_file("results/serial_result.txt");
    out_file << "--- Resultado Serial Epistasis ---" << endl;
    out_file << "SNPs: " << num_snps << ", Pacientes: " << num_patients << endl;
    out_file << "Mejor Par: (rs" << best_i << ", rs" << best_j << ")" << endl;
    out_file << "Max Score (Pearson): " << fixed << setprecision(6) << global_max << endl;
    out_file << "Tiempo de ejecucion: " << diff.count() << " segundos." << endl;
    out_file.close();

    cout << "Calculo finalizado. Resultado guardado en 'results/serial_result.txt'" << endl;
    cout << "Tiempo: " << diff.count() << " s" << endl;

    return 0;
}