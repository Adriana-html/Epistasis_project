#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <fstream>

using namespace std;

double calculate_pearson(const int* a, const int* b, const int* y, int M){
    double sx =0, sy=0, sx2=0, sy2=0, sxy=0;
    for (int i=0;i<M;i++){
        double X = a[i]*b[i];
        double Y = y[i];
        sx+=X; sy+=Y;
        sx2+=X*X; sy2+=Y*Y;
        sxy+=X*Y;
    }

    double num = M*sxy - sx*sy;
    double den = (M*sx2 - sx*sx)*(M*sy2 - sy*sy);
    if (den<=0) return 0.0;
    return num / sqrt(den);
}
// k → (i,j)
inline void k_to_ij(long long k, int N, int &i, int &j){
    long long T = (long long)N*(N-1)/2;
    i = N - 2 - floor((sqrt(8*(T-1-k)+1)-1)/2);
    j = k - ( (long long)N*i - (long long)i*(i+1)/2 ) + i + 1;
}

int main(int argc, char** argv){

    MPI_Init(&argc,&argv);
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    if(argc<3){
        if(rank==0)
            cout<<"Usage: ./app <num_patients> <num_snps>\n";
        MPI_Finalize();
        return 0;
    }

    int M = atoi(argv[1]);
    int N = atoi(argv[2]);

    //Total combinaciones 
    long long T = (long long)N*(N-1)/2;
    long long start = (rank*T)/size;
    long long end   = ((rank+1)*T)/size;
    if (rank == 0) {
        cout << "Total tareas (T): " << T << endl;
        cout << "Tareas por proceso aprox: " << T / size << endl;
    }
    //identiticar SNPs necesarios
    unordered_set<int> needed;
    for(long long k=start;k<end;k++){
        int i,j;
        k_to_ij(k,N,i,j);
        needed.insert(i);
        needed.insert(j);
    }
    vector<int> snp_ids(needed.begin(), needed.end());
    sort(snp_ids.begin(), snp_ids.end());
    int S = snp_ids.size();

    vector<int> snp_data(S * M);
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD,"genotypes.bin",
                  MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

    // leer por bloques contiguos
    int start_block = 0;

    while(start_block < S){
        int end_block = start_block;

        // agrupar consecutivos
        while(end_block+1 < S &&
              snp_ids[end_block+1] == snp_ids[end_block]+1){
            end_block++;
        }

        int block_size = end_block - start_block + 1;

        MPI_Offset offset =
            (MPI_Offset)snp_ids[start_block] * M * sizeof(int);

        MPI_File_read_at_all(
            fh,
            offset,
            &snp_data[start_block*M],
            block_size*M,
            MPI_INT,
            MPI_STATUS_IGNORE
        );

        start_block = end_block + 1;
    }

    MPI_File_close(&fh);


    //Fenotipo
    vector<int> phenotype(M);

    if(rank==0){
        MPI_File fh2;
        MPI_File_open(MPI_COMM_SELF,"phenotypes.bin",
                      MPI_MODE_RDONLY,MPI_INFO_NULL,&fh2);

        MPI_File_read(fh2,phenotype.data(),M,MPI_INT,MPI_STATUS_IGNORE);
        MPI_File_close(&fh2);
    }

    MPI_Bcast(phenotype.data(),M,MPI_INT,0,MPI_COMM_WORLD);

    //índice rápido (id -> posicion)
    vector<int> index_map(N,-1);
    for(int i=0;i<S;i++){
        index_map[snp_ids[i]] = i;
    }
    double local_max = -numeric_limits<double>::infinity();
    int best_i_local = -1, best_j_local = -1;

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();
    

    if (rank == 0) {
        cout << "Check SNP_DATA[0]: " << snp_data[0] << endl;
        cout << "Check PHENOTYPE[0]: " << phenotype[0] << endl;
    }
    for(long long k=start;k<end;k++){
        int i,j;
        k_to_ij(k,N,i,j);

        int li = index_map[i];
        int lj = index_map[j];

        double score = calculate_pearson(
            &snp_data[li*M],
            &snp_data[lj*M],
            phenotype.data(),
            M
        );

        if(score > local_max){
            local_max = score;
            best_i_local = i;
            best_j_local = j;
        }
    }

    // --- ESCRITURA DE ARCHIVOS INDIVIDUALES POR CPU ---
    string filename = "cpu_results_" + to_string(rank) + ".txt";
    ofstream outfile(filename);
    outfile << "Rank: " << rank << endl;
    outfile << "Rango de tareas: " << start << " a " << end << endl;
    outfile << "Mejor SNP i local: " << best_i_local << endl;
    outfile << "Mejor SNP j local: " << best_j_local << endl;
    outfile << "Score local: " << local_max << endl;
    outfile.close();

    // --- REDUCCIÓN PARA EL GANADOR GLOBAL ---
    struct {
        double val;
        int rank;
    } local_out, global_out;

    local_out.val = local_max;
    local_out.rank = rank;
    

    MPI_Reduce(&local_out, &global_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
    

    double t1 = MPI_Wtime();

    if(rank==0){
        cout << "================================\n";
        cout << "Tiempo: " << (t1 - t0) << " s\n";
        cout << "Max score global: " << global_out.val << "\n";
        cout << "Encontrado por CPU rank: " << global_out.rank << "\n";
        cout << "================================\n";
        cout << "Revisa los archivos cpu_results_X.txt para ver el detalle de cada nucleo.\n";
    }

    MPI_Finalize();
    return 0;
}