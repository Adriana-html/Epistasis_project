from mpi4py import MPI
import numpy as np
import pandas as pd
import argparse
import sys
import os

def write_mpi_io_csv(comm, rank, filename, data_bytes, header_df):
    """
    Writes variable-length byte chunks to a shared file using MPI-IO.
    """
    # 1. Rank 0 generates the header and broadcasts its length to all processors
    header_bytes = b""
    if rank == 0: #rank 0 es el proceso maestro
        header_str = header_df.to_csv()
        header_bytes = header_str.encode('utf-8')
    header_len = comm.bcast(len(header_bytes), root=0)

    # 2. Calculate offsets using Exscan
    local_size = len(data_bytes)
    offset_arr = np.array([0], dtype=np.int64)
    local_size_arr = np.array([local_size], dtype=np.int64)

    # Exscan computes the sum of local_sizes BEFORE the current rank
    comm.Exscan(local_size_arr, offset_arr, op=MPI.SUM)
    if rank == 0:
        offset_arr[0] = 0  # Rank 0 always starts data exactly after the header

    actual_offset = header_len + offset_arr[0]

    # 3. Open and Write via MPI-IO
    amode = MPI.MODE_WRONLY | MPI.MODE_CREATE
    fh = MPI.File.Open(comm, filename, amode)

    if rank == 0:
        fh.Write_at(0, header_bytes) # Write the header at the very beginning

    # Collective write: all processors write their chunks at their computed offsets
    fh.Write_at_all(actual_offset, data_bytes)
    fh.Close()

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # 1. Parse arguments
    parser = argparse.ArgumentParser(description="Parallel GWAS Simulator")
    parser.add_argument("--n_individuals", type=int, default=1000, help="Total number of individuals")
    parser.add_argument("--n_snps", type=int, default=5000, help="Total number of SNPs")
    args = parser.parse_args()

    n_individuals = args.n_individuals # Total number of individuals to simulate, M value
    n_snps = args.n_snps # Total number of SNPs to simulate, N value

    # 2. Check workload divisibility and logic requirements
    if n_individuals % size != 0: #verify that the total number of individuals can be perfectly divided by the number of processors, otherwise stop immediately
        if rank == 0:
            print(f"Error: Total individuals ({n_individuals}) must be perfectly divisible by the number of MPI processors ({size}).")
        sys.exit(1)

    if n_snps <= 1200: #verify that the total number of SNPs is greater than 1200, otherwise the hidden interaction logic will not work, and stop immediately
        if rank == 0:
            print(f"Error: n_snps ({n_snps}) must be greater than 1200 for the hidden interaction logic to work.")
        sys.exit(1)

    # Make sure the files does not exist, otherwise stop immediatelly
    if os.path.exists("genotypes.csv"):
        if rank == 0:
            print(f"Error: the file genotypes.csv already exist.")
        sys.exit(1)

    # Make sure the files does not exist, otherwise stop immediatelly
    if os.path.exists("phenotypes.csv"):
        if rank == 0:
            print(f"Error: the file phenotypes.csv already exist.")
        sys.exit(1)

    comm.Barrier() # Ensure files doesnt exist to continue
        
    # 3. Partition the workload
    chunk_size = n_individuals // size
    start_idx = rank * chunk_size
    end_idx = start_idx + chunk_size

    # Set a unique random seed per processor so they don't generate identical data!
    np.random.seed(42 + rank) 

    if rank == 0:
        print(f"Starting simulation for {n_individuals} individuals and {n_snps} SNPs across {size} processors...")

    # 4. Generate Local Genotypes
    local_genotypes = np.random.choice([0, 1, 2], size=(chunk_size, n_snps), p=[0.64, 0.32, 0.04])
    snp_ids = [f"rs{i}" for i in range(n_snps)]
    local_ind_ids = [f"Indiv_{i}" for i in range(start_idx, end_idx)]

    # 5. Generate Local Phenotypes
    local_phenotypes = np.zeros(chunk_size, dtype=int)
    for i in range(chunk_size):
        if local_genotypes[i, 500] == 2 and local_genotypes[i, 1200] >= 1:
            local_phenotypes[i] = np.random.choice([0, 1], p=[0.1, 0.9])
        else:
            local_phenotypes[i] = np.random.choice([0, 1], p=[0.95, 0.05])

    # 6. Convert local arrays to Pandas DataFrames, then to CSV byte strings
    df_gen = pd.DataFrame(local_genotypes, index=local_ind_ids, columns=snp_ids)
    df_pheno = pd.DataFrame(local_phenotypes, index=local_ind_ids, columns=["Trait_Status"])

    # We skip headers here because Rank 0 will write the master header separately
    gen_bytes = df_gen.to_csv(header=False).encode('utf-8')
    pheno_bytes = df_pheno.to_csv(header=False).encode('utf-8')

    # Create empty dataframes just to generate the CSV headers efficiently on Rank 0
    header_gen = pd.DataFrame(columns=snp_ids)
    header_pheno = pd.DataFrame(columns=["Trait_Status"])

    # 7. Write to shared files via MPI-IO
    write_mpi_io_csv(comm, rank, "data/genotypes.csv", gen_bytes, header_gen)
    write_mpi_io_csv(comm, rank, "data/phenotypes.csv", pheno_bytes, header_pheno)

    if rank == 0:
        print("Done! Both files successfully generated and written in parallel via MPI-IO.")

if __name__ == "__main__":
    main()
