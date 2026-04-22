#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <chrono>

using namespace std;

// Calculates the Pearson correlation for the interaction of two SNPs
double calculate_pearson(const vector<unsigned int>& snp1, const vector<unsigned int>& snp2, const vector<unsigned int>& phenotype)
{
 // The size of the phenotype, number of samples
 int M = phenotype.size();
 
 double sum_X = 0, sum_Y = 0, sum_X2 = 0, sum_Y2 = 0, sum_XY = 0;

 // Loop over all the samples
 for (unsigned int i = 0; i < M; ++i)
  {
  // The epistatic interaction is the product of the two SNP values
  double X = snp1[i] * snp2[i]; //AND logic for the interaction
  double Y = phenotype[i];
  
  sum_X += X;
  sum_Y += Y;
  sum_X2 += X * X;
  sum_Y2 += Y * Y;
  sum_XY += X * Y;
 }

 // Pearson correlation
 double numerator = (M * sum_XY) - (sum_X * sum_Y);
 double denominator_sq = ((M * sum_X2) - (sum_X * sum_X)) * ((M * sum_Y2) - (sum_Y * sum_Y));
 
 // Handle edge cases like zero variance to avoid NaN results
 if (denominator_sq <= 0.0) return 0.0; //caso divicion por cero, no hay correlacion
 
 return numerator / sqrt(denominator_sq);
 
}

int main(int argc, const char **argv)
{
 const unsigned int num_patients = 500;
 const unsigned int num_snps = 1500;

 // Temporal variables to store the reading values from files
 string line, val, id;
 
 // 1. Read Phenotypes
 vector<unsigned int> phenotypes(num_patients); //vector de fenotipos, cada fenotipo es un valor binario (0 o 1) que indica la presencia o ausencia de la enfermedad en el paciente correspondiente
 ifstream pheno_file("phenotypes.csv"); //validamos que el archivo exista y se puede abrir correctamente.

 // Check for error at opening the file
 if (!pheno_file.is_open())
  {
   cerr << "Error opening phenotypes.csv" << endl;
   return 1;
  }

 // Skip the first line
 getline(pheno_file, line); 

 // Load data into memory
 for (unsigned int i = 0; i < num_patients; ++i)
  {
   // Get a line and parse it
   getline(pheno_file, line);
   stringstream ss(line);

   getline(ss, id, ','); // Ignore the identifier o ID del paciente.
   getline(ss, val); // Get the actual binary phenotype, valor que indica si el paciente tiene la enfermedad (1) o no (0).
   
   phenotypes.push_back(stoi(val)); //convertimos el valor de string a entero y lo almacenamos en el vector de fenotipos.
   //sugerencia
   //phenotypes[i] = stoi(val); //alternativamente, podríamos asignar directamente al índice correspondiente del vector, ya que sabemos que el número de pacientes es fijo y coincide con el tamaño del vector.
  }
 pheno_file.close();
 
 // 2. Read Genotypes
 // We store as snps[snp_index][patient_index] for easier access later
 //snps[i][j] representa el genotipo del paciente j para el SNP i, donde el valor es 0, 1 o 2 dependiendo de la cantidad de alelos de riesgo presentes en ese SNP para ese paciente.
 vector<vector<unsigned int> > snps(num_snps, vector<unsigned int>(num_patients));
 ifstream geno_file("genotypes.csv");
 
 // Check for errors at opening the files
 if (!geno_file.is_open()) {
  cerr << "Error opening genotypes.csv" << endl;
  return 1;
 }
 
 // Skip the first line
 getline(geno_file, line);
 
 // Load data into memory
 for (unsigned int i = 0; i < num_patients; ++i)
  {
   // Read a line
   getline(geno_file, line);
   stringstream ss(line);
   
   getline(ss, id, ','); // Ignore the identifier
   
   for (unsigned int s = 0; s < num_snps; ++s)
    {
     // Parse the reading
     getline(ss, val, ',');
     snps[s][i] = stoi(val); // s es el índice del SNP, i es el índice del paciente.
    }
   
  }
 geno_file.close();
 
 // 3. Compute Epistatic Interactions (All unique pairs)
 cout << "--- Epistasis Interaction Scores ---" << endl;
 cout << fixed << setprecision(5);

 long long total_pairs = 0;
 unsigned int print_limit = 15; // Limit output to avoid flooding the console
 
 // N*(N-1)/2 comparisons.
 //pares unicos, es decir, sin repetir (rs1, rs2) y (rs2, rs1) y sin comparar un SNP consigo mismo (rs1, rs1).
 cout <<"Iniciando cálculo..." << endl;
 auto start_time = chrono::high_resolution_clock::now(); // Iniciamos el temporizador para medir el tiempo de ejecución del cálculo de interacciones.
 for (unsigned int i = 0; i < num_snps; ++i)
  {
  for (unsigned int j = i + 1; j < num_snps; ++j)
   {
    //para cada par se calcula la correlacion de Pearson entre la interaccion de los dos SNPs
    double r = calculate_pearson(snps[i], snps[j], phenotypes);
    total_pairs++;
    
    if (total_pairs <= print_limit)
     {
      cout << "Pair (rs" << i << ", rs" << j << "): r = " << r << endl;
     } //evitamos saturar la consola
    else if (total_pairs == print_limit + 1)
     {
      cout << "... continuing computations ..." << endl;
     }
   }
  
  }
 auto end_time = chrono::high_resolution_clock::now(); // Finalizamos el temporizador después de completar el cálculo de interacciones.
 std::chrono::duration<double> diff = end_time - start_time;
 cout << "Time taken for epistasis computation: " << diff.count() << " seconds." << endl;

 cout << "\nSuccessfully computed " << total_pairs << " interaction pairs." << endl;
 
 return 0;

}
