#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <mpi.h>

#include "minimizerstandard.h"

#include "wavestandard.h"
#include "wavesimple.h"
#include "montecarlostandard.h"
#include "matrix.h"
#include "hamiltonianstandard.h"
#include "hamiltoniansimple.h"

MinimizerStandard::MinimizerStandard(int rank, int nProcesses) :
    Minimizer(rank, nProcesses),
    my_rank(rank),
    numprocs(nProcesses),
    dimension(3),
    charge(2),
    step_length(1.0),
    number_particles(2)
{
}

void MinimizerStandard::run()
{
    WaveSimple *wave = new WaveSimple(number_particles, dimension);
    HamiltonianSimple *hamiltonian = new HamiltonianSimple(number_particles, dimension, charge);
    string outfilename;
    int total_number_cycles, i;
    double *cumulative_e, *cumulative_e2;
    double *total_cumulative_e, *total_cumulative_e2;
    double *all_energies;
    double  time_start, time_end, total_time;
    int number_cycles = 1000000; //  default number of cycles
    int max_variations = 11;    //  default number of variations
    double alpha, variance, energy, error;

    time_start = MPI_Wtime();

    if (my_rank == 0) {
        outfilename = "output.dat";
        ofile.open(outfilename.c_str());
    }

    // Setting output file name for this rank:
    ostringstream ost;
    ost << "blocks_rank" << my_rank << ".dat";
    // Open file for writing:
    blockofile.open(ost.str().c_str(), ios::out | ios::binary);

    total_cumulative_e = new double[max_variations+1];
    total_cumulative_e2 = new double[max_variations+1];
    cumulative_e = new double[max_variations+1];
    cumulative_e2 = new double[max_variations+1];

    //  initialize the arrays  by zeroing them
    for( i=1; i <= max_variations; i++){
        cumulative_e[i] = cumulative_e2[i]  = total_cumulative_e[i] = total_cumulative_e2[i]  = 0.0;
    }

    // broadcast the total number of  variations
    MPI_Bcast (&max_variations, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&number_cycles, 1, MPI_INT, 0, MPI_COMM_WORLD);

    total_number_cycles = number_cycles*numprocs;

    // array to store all energies for last variation of alpha
    all_energies = new double[number_cycles+1];

    //  Do the mc sampling  and accumulate data with MPI_Reduce
    MonteCarloStandard *monteCarlo = new MonteCarloStandard(wave, hamiltonian, number_particles, dimension, charge, my_rank, step_length);
    monteCarlo->sample(max_variations, number_cycles, cumulative_e, cumulative_e2,
                all_energies);
    //  Collect data in total averages
    for( i=1; i <= max_variations; i++){
        MPI_Reduce(&cumulative_e[i], &total_cumulative_e[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&cumulative_e2[i], &total_cumulative_e2[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    time_end = MPI_Wtime();
    total_time = time_end-time_start;
    // Print out results
    if ( my_rank == 0) {
        cout << "Time = " <<  total_time  << " on number of processors: "  << numprocs  << endl;
        alpha = 0.5*charge;
        for( i=1; i <= max_variations; i++){
            energy = total_cumulative_e[i]/total_number_cycles;
            variance = total_cumulative_e2[i]/total_number_cycles-energy*energy;
            error=sqrt(variance/(total_number_cycles-1));
            ofile << setiosflags(ios::showpoint | ios::uppercase);
            ofile << setw(15) << setprecision(8) << alpha;
            ofile << setw(15) << setprecision(8) << energy;
            ofile << setw(15) << setprecision(8) << variance;
            ofile << setw(15) << setprecision(8) << error << endl;
            alpha += 0.1;
        }
        ofile.close();  // close output file
    }
    blockofile.write((char*)(all_energies+1),
                     number_cycles*sizeof(double));
    blockofile.close();
    delete [] total_cumulative_e; delete [] total_cumulative_e2;
    delete [] cumulative_e; delete [] cumulative_e2; delete [] all_energies;
}
