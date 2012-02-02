// Variational Monte Carlo for atoms with up to two electrons

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <mpi.h>

#include "wavestandard.h"
#include "montecarlostandard.h"
#include "matrix.h"
using namespace  std;
// output file as global variable
ofstream ofile, blockofile;

//  Here we define global variables  used in various functions
//  These can be changed by reading from file the different parameters
int dimension = 3; // three-dimensional system
int charge = 2;  //  we fix the charge to be that of the helium atom
int my_rank, numprocs;  //  these are the parameters used by MPI  to define which node and how many
double step_length = 1.0;  //  we fix the brute force jump to 1 Bohr radius
int number_particles  = 2;  //  we fix also the number of electrons to be 2


// Begin of main program

int main(int argc, char* argv[])
{
    WaveStandard *wave = new WaveStandard(number_particles, dimension);
    char *outfilename;
    char *blockoutfilename;
    int total_number_cycles, i;
    double *cumulative_e, *cumulative_e2;
    double *total_cumulative_e, *total_cumulative_e2;
    double *all_energies;
    double  time_start, time_end, total_time;
    int number_cycles = 1000000; //  default number of cycles
    int max_variations = 10;    //  default number of variations
    double alpha, variance, energy, error;

    //  MPI initializations
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    time_start = MPI_Wtime();

    if (my_rank == 0 && argc <= 1) {
        cout << "Bad Usage: " << argv[0] <<
                " read also output file on same line" << endl;
        //    exit(1);
    }
    if (my_rank == 0 && argc > 1) {
        outfilename=argv[1];
        ofile.open(outfilename);
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
    MonteCarloStandard *monteCarlo = new MonteCarloStandard(wave, number_particles, dimension, charge, my_rank, step_length);
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
            alpha += 0.1;
            energy = total_cumulative_e[i]/total_number_cycles;
            variance = total_cumulative_e2[i]/total_number_cycles-energy*energy;
            error=sqrt(variance/(total_number_cycles-1));
            ofile << setiosflags(ios::showpoint | ios::uppercase);
            ofile << setw(15) << setprecision(8) << alpha;
            ofile << setw(15) << setprecision(8) << energy;
            ofile << setw(15) << setprecision(8) << variance;
            ofile << setw(15) << setprecision(8) << error << endl;
        }
        ofile.close();  // close output file
    }
    blockofile.write((char*)(all_energies+1),
                     number_cycles*sizeof(double));
    blockofile.close();
    delete [] total_cumulative_e; delete [] total_cumulative_e2;
    delete [] cumulative_e; delete [] cumulative_e2; delete [] all_energies;
    // End MPI
    MPI_Finalize ();

    return 0;
}  //  end of main function

