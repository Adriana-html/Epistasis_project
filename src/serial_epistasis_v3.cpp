//misma implementacion serial pero con el cambio de tipo de datos a int8_t para genotipos y double para fenotipos, ademas de algunos ajustes menores en la lectura de archivos y escritura de resultados. Este código es la versión usada para el simulador de terceros
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <cstdint>

using namespace std;


double calculate_pearson(const int8_t* snp1, const int8_t* snp2, const double* phenotype, int M) {
    double sum_X = 0, sum_Y = 0, sum_X2 = 0, sum_Y2 = 0, sum_XY = 0;

    for (int i = 0; i < M; ++i) {
        
        double X = (double)(snp1[i] * snp2[i]); 
        double Y = phenotype[i]; 
        
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

    
    vector<double> phenotypes(num_patients);
    
    
    FILE* f_pheno = fopen("data/phenotypes.bin", "rb"); 
    if (!f_pheno) {
        cerr << "Error: No se pudo abrir data/phenotypes.bin" << endl;
        return 1;
    }
    
    
    if (fread(phenotypes.data(), sizeof(double), num_patients, f_pheno) != (size_t)num_patients) {
        cerr << "Error leyendo phenotypes.bin" << endl;
    }
    fclose(f_pheno);

    
    vector<int8_t> snps_flat((long long)num_snps * num_patients);
    FILE* f_geno = fopen("data/genotypes.bin", "rb");
    if (!f_geno) {
        cerr << "Error: No se pudo abrir data/genotypes.bin" << endl;
        return 1;
    }
    size_t total_elements = (size_t)num_snps * num_patients;
    if (fread(snps_flat.data(), sizeof(int8_t), total_elements, f_geno) != total_elements) {
        cerr << "Error leyendo genotypes.bin" << endl;
    }
    fclose(f_geno);

    cout << "Datos cargados correctamente (Phenotypes: double, Genotypes: int8). Iniciando..." << endl;
    
    double global_max = -1.0;
    int best_i = -1, best_j = -1;

    auto start_time = chrono::high_resolution_clock::now();
    
    for (int i = 0; i < num_snps; ++i) {
        for (int j = i + 1; j < num_snps; ++j) {
            double r = calculate_pearson(&snps_flat[(long long)i * num_patients], 
                                         &snps_flat[(long long)j * num_patients], 
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

    // Escritura de resultados
    ofstream out_file("results/serial_result_v3.txt");
    out_file << "--- Resultado Serial Epistasis V3 (Terceros) ---" << endl;
    out_file << "SNPs: " << num_snps << ", Pacientes: " << num_patients << endl;
    out_file << "Max Score (Pearson): " << fixed << setprecision(6) << global_max << endl;
    out_file << "Tiempo: " << diff.count() << " s" << endl;
    out_file.close();

    cout << "Tiempo: " << diff.count() << " s" << endl;

    return 0;
}