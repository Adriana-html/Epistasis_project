//codigo usado para el simulador de terceros.
#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <limits>
#include <fstream>

using namespace std;


double calculate_pearson(const int8_t* a, const int8_t* b, const double* y, int M) {
    double sx = 0, sy = 0, sx2 = 0, sy2 = 0, sxy = 0;
    for (int i = 0; i < M; i++) {
        
        double X = (int)a[i] * (int)b[i];
        double Y = y[i]; // Fenotipo cuantitativo
        sx += X; sy += Y;
        sx2 += X * X; sy2 += Y * Y;
        sxy += X * Y;
    }

    double num = (double)M * sxy - sx * sy;
    double den = ((double)M * sx2 - sx * sx) * ((double)M * sy2 - sy * sy);
    if (den <= 0) return 0.0;
    return abs(num / sqrt(den)); // Retornamos valor absoluto de la correlación
}

// k → (i,j) 
inline void k_to_ij(long long k, int N, int &i, int &j) {
    long long T = (long long)N * (N - 1) / 2;
    i = N - 2 - floor((sqrt(8 * (T - 1 - k) + 1) - 1) / 2);
    j = k - ((long long)N * i - (long long)i * (i + 1) / 2) + i + 1;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 3) {
        if (rank == 0)
            cout << "Uso: mpirun -np N ./mpi_epistasis <num_patients> <num_snps>\n";
        MPI_Finalize();
        return 0;
    }

    int M = atoi(argv[1]); // num_patients
    int N = atoi(argv[2]); // num_snps

    // Total combinaciones (Carga de trabajo)
    long long T = (long long)N * (N - 1) / 2;
    long long start = (rank * T) / size;
    long long end = ((rank + 1) * T) / size;

    if (rank == 0) {
        cout << "Total tareas (T): " << T << endl;
        cout << "Procesando con " << size << " CPUs..." << endl;
    }

    // Identificar SNPs necesarios para este proceso
    unordered_set<int> needed;
    for (long long k = start; k < end; k++) {
        int i, j;
        k_to_ij(k, N, i, j);
        needed.insert(i);
        needed.insert(j);
    }
    vector<int> snp_ids(needed.begin(), needed.end());
    sort(snp_ids.begin(), snp_ids.end());
    int S = snp_ids.size();

    // Buffer para genotipos
    vector<int8_t> snp_data(S * M);
    MPI_File fh;
    int err = MPI_File_open(MPI_COMM_WORLD, "data/genotypes.bin",
                           MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    
    if(err != MPI_SUCCESS) {
        if(rank == 0) cout << "Error al abrir genotypes_epigen.bin" << endl;
        MPI_Finalize();
        return 1;
    }

   
    int start_block = 0;
    while (start_block < S) {
        int end_block = start_block;
        while (end_block + 1 < S && snp_ids[end_block + 1] == snp_ids[end_block] + 1) {
            end_block++;
        }
        int block_size = end_block - start_block + 1;

        
        MPI_Offset offset = (MPI_Offset)snp_ids[start_block] * M * sizeof(int8_t);

        MPI_File_read_at_all(
            fh,
            offset,
            &snp_data[start_block * M],
            block_size * M,
            MPI_INT8_T, 
            MPI_STATUS_IGNORE
        );
        start_block = end_block + 1;
    }
    MPI_File_close(&fh);

    // Lectura de Fenotipo
    vector<double> phenotype(M);
    if (rank == 0) {
        MPI_File fh2;
        MPI_File_open(MPI_COMM_SELF, "data/phenotypes.bin",
                      MPI_MODE_RDONLY, MPI_INFO_NULL, &fh2);
        
        MPI_File_read(fh2, phenotype.data(), M, MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_close(&fh2);
    }
    
    MPI_Bcast(phenotype.data(), M, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    
    vector<int> index_map(N, -1);
    for (int i = 0; i < S; i++) {
        index_map[snp_ids[i]] = i;
    }

    double local_max = -1.0;
    int best_i_local = -1, best_j_local = -1;

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    for (long long k = start; k < end; k++) {
        int i, j;
        k_to_ij(k, N, i, j);

        int li = index_map[i];
        int lj = index_map[j];

        double score = calculate_pearson(
            &snp_data[li * M],
            &snp_data[lj * M],
            phenotype.data(),
            M
        );

        if (score > local_max) {
            local_max = score;
            best_i_local = i;
            best_j_local = j;
        }
    }

    double t1 = MPI_Wtime();

    
    string filename = "results/simulador/cpu_results_" + to_string(rank) + ".txt";
    ofstream outfile(filename);
    outfile << "Rank: " << rank << endl;
    outfile << "Mejor SNP i local: " << best_i_local << endl;
    outfile << "Mejor SNP j local: " << best_j_local << endl;
    outfile << "Score local: " << local_max << endl;
    outfile << "Tiempo local: " << (t1 - t0) << " s" << endl;
    outfile.close();

    
    struct {
        double val;
        int rank;
    } local_out, global_out;

    local_out.val = local_max;
    local_out.rank = rank;

    MPI_Reduce(&local_out, &global_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "================================\n";
        cout << "Tiempo: " << (t1 - t0) << " s\n";
        cout << "Max score global: " << global_out.val << "\n";
        cout << "Encontrado por CPU rank: " << global_out.rank << "\n";
        cout << "================================\n";
    }

    MPI_Finalize();
    return 0;
}