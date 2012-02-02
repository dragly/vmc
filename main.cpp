// Variational Monte Carlo for atoms with up to two electrons

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "mpi.h"
using namespace  std;
// output file as global variable
ofstream ofile, blockofile;
// the step length and its squared inverse for the second derivative
#define h 0.001
#define h2 1000000

// declaration of functions

// The Mc sampling for the variational Monte Carlo
void  mc_sampling(int, int, double *, double *, double *);

// The variational wave function
double  wave_function(double **, double);

// The local energy
double  local_energy(double **, double, double);

//  allocate space for a matrix
void  **matrix(int, int, int);

//  free space for  a matrix
void free_matrix(void **);

// ran2 for uniform deviates, initialize with negative seed.
double ran2(long *);


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

    if (my_rank == 0 && argc <= 2) {
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
    mc_sampling(max_variations, number_cycles, cumulative_e, cumulative_e2,
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


// Monte Carlo sampling with the Metropolis algorithm

void mc_sampling(int max_variations, int number_cycles, double *cumulative_e, double *cumulative_e2, double *all_energies)
{
    int cycles, variate, accept, dim, i, j, k;
    long idum;
    double wfnew, wfold, alpha, energy, energy2, delta_e;
    double **r_old, **r_new;
    alpha = 0.5*charge;

    // every node has its own seed for the random numbers
    idum = -1-my_rank;
    // allocate matrices which contain the position of the particles
    r_old = (double **) matrix( number_particles, dimension, sizeof(double));
    r_new = (double **) matrix( number_particles, dimension, sizeof(double));
    for (i = 0; i < number_particles; i++) {
        for ( j=0; j < dimension; j++) {
            r_old[i][j] = r_new[i][j] = 0;
        }
    }
    // loop over variational parameters
    for (variate=1; variate <= max_variations; variate++){
        // initialisations of variational parameters and energies
        alpha += 0.1;
        energy = energy2 = 0; accept =0; delta_e=0;
        //  initial trial position, note calling with alpha
        for (i = 0; i < number_particles; i++) {
            for ( j=0; j < dimension; j++) {
                r_old[i][j] = step_length*(ran2(&idum)-0.5);
            }
        }
        wfold = wave_function(r_old, alpha);
        // loop over monte carlo cycles
        for (cycles = 1; cycles <= number_cycles; cycles++){
            // new position
            for (i = 0; i < number_particles; i++) {
                for ( j=0; j < dimension; j++) {
                    r_new[i][j] = r_old[i][j]+step_length*(ran2(&idum)-0.5);
                }
                //  for the other particles we need to set the position to the old position since
                //  we move only one particle at the time
                for (k = 0; k < number_particles; k++) {
                    if ( k != i) {
                        for ( j=0; j < dimension; j++) {
                            r_new[k][j] = r_old[k][j];
                        }
                    }
                }
                wfnew = wave_function(r_new, alpha);
                // The Metropolis test is performed by moving one particle at the time
                if(ran2(&idum) <= wfnew*wfnew/wfold/wfold ) {
                    for ( j=0; j < dimension; j++) {
                        r_old[i][j]=r_new[i][j];
                    }
                    wfold = wfnew;
                }
            }  //  end of loop over particles
            // compute local energy
            delta_e = local_energy(r_old, alpha, wfold);
            // save all energies on last variate
            if(variate==max_variations){
                all_energies[cycles] = delta_e;
            }
            // update energies
            energy += delta_e;
            energy2 += delta_e*delta_e;
        }   // end of loop over MC trials
        // update the energy average and its squared
        cumulative_e[variate] = energy;
        cumulative_e2[variate] = energy2;
    }    // end of loop over variational  steps
    free_matrix((void **) r_old); // free memory
    free_matrix((void **) r_new); // free memory
}   // end mc_sampling function


// Function to compute the squared wave function, simplest form

double  wave_function(double **r, double alpha)
{
    int i, j, k;
    double wf, argument, r_single_particle, r_12;

    argument = wf = 0;
    for (i = 0; i < number_particles; i++) {
        r_single_particle = 0;
        for (j = 0; j < dimension; j++) {
            r_single_particle  += r[i][j]*r[i][j];
        }
        argument += sqrt(r_single_particle);
    }
    wf = exp(-argument*alpha) ;
    return wf;
}


// Function to calculate the local energy with num derivative

double  local_energy(double **r, double alpha, double wfold)
{
    int i, j , k;
    double e_local, wfminus, wfplus, e_kinetic, e_potential, r_12,
            r_single_particle;
    double **r_plus, **r_minus;

    // allocate matrices which contain the position of the particles
    // the function matrix is defined in the progam library
    r_plus = (double **) matrix( number_particles, dimension, sizeof(double));
    r_minus = (double **) matrix( number_particles, dimension, sizeof(double));
    for (i = 0; i < number_particles; i++) {
        for ( j=0; j < dimension; j++) {
            r_plus[i][j] = r_minus[i][j] = r[i][j];
        }
    }
    // compute the kinetic energy
    e_kinetic = 0;
    for (i = 0; i < number_particles; i++) {
        for (j = 0; j < dimension; j++) {
            r_plus[i][j] = r[i][j]+h;
            r_minus[i][j] = r[i][j]-h;
            wfminus = wave_function(r_minus, alpha);
            wfplus  = wave_function(r_plus, alpha);
            e_kinetic -= (wfminus+wfplus-2*wfold);
            r_plus[i][j] = r[i][j];
            r_minus[i][j] = r[i][j];
        }
    }
    // include electron mass and hbar squared and divide by wave function
    e_kinetic = 0.5*h2*e_kinetic/wfold;
    // compute the potential energy
    e_potential = 0;
    // contribution from electron-proton potential
    for (i = 0; i < number_particles; i++) {
        r_single_particle = 0;
        for (j = 0; j < dimension; j++) {
            r_single_particle += r[i][j]*r[i][j];
        }
        e_potential -= charge/sqrt(r_single_particle);
    }
    // contribution from electron-electron potential
    for (i = 0; i < number_particles-1; i++) {
        for (j = i+1; j < number_particles; j++) {
            r_12 = 0;
            for (k = 0; k < dimension; k++) {
                r_12 += (r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
            }
            e_potential += 1/sqrt(r_12);
        }
    }
    free_matrix((void **) r_plus); // free memory
    free_matrix((void **) r_minus);
    e_local = e_potential+e_kinetic;
    return e_local;
}


/*
 * The function
 *      void  **matrix()
 * reserves dynamic memory for a two-dimensional matrix
 * using the C++ command new . No initialization of the elements.
 * Input data:
 *  int row      - number of  rows
 *  int col      - number of columns
 *  int num_bytes- number of bytes for each
 *                 element
 * Returns a void  **pointer to the reserved memory location.
 */

void **matrix(int row, int col, int num_bytes)
{
    int      i, num;
    char     **pointer, *ptr;

    pointer = new(nothrow) char* [row];
    if(!pointer) {
        cout << "Exception handling: Memory allocation failed";
        cout << " for "<< row << "row addresses !" << endl;
        return NULL;
    }
    i = (row * col * num_bytes)/sizeof(char);
    pointer[0] = new(nothrow) char [i];
    if(!pointer[0]) {
        cout << "Exception handling: Memory allocation failed";
        cout << " for address to " << i << " characters !" << endl;
        return NULL;
    }
    ptr = pointer[0];
    num = col * num_bytes;
    for(i = 0; i < row; i++, ptr += num )   {
        pointer[i] = ptr;
    }

    return  (void **)pointer;

} // end: function void **matrix()

/*
   * The function
   *      void free_matrix()
   * releases the memory reserved by the function matrix()
   *for the two-dimensional matrix[][]
   * Input data:
   *  void far **matr - pointer to the matrix
   */

void free_matrix(void **matr)
{

    delete [] (char *) matr[0];
    delete [] matr;

}  // End:  function free_matrix()


/*
** The function
**         ran2()
** is a long periode (> 2 x 10^18) random number generator of
** L'Ecuyer and Bays-Durham shuffle and added safeguards.
** Call with idum a negative integer to initialize; thereafter,
** do not alter idum between sucessive deviates in a
** sequence. RNMX should approximate the largest floating point value
** that is less than 1.
** The function returns a uniform deviate between 0.0 and 1.0
** (exclusive of end-point values).
*/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
    int            j;
    long           k;
    static long    idum2 = 123456789;
    static long    iy=0;
    static long    iv[NTAB];
    double         temp;

    if(*idum <= 0) {
        if(-(*idum) < 1) *idum = 1;
        else             *idum = -(*idum);
        idum2 = (*idum);
        for(j = NTAB + 7; j >= 0; j--) {
            k     = (*idum)/IQ1;
            *idum = IA1*(*idum - k*IQ1) - k*IR1;
            if(*idum < 0) *idum +=  IM1;
            if(j < NTAB)  iv[j]  = *idum;
        }
        iy=iv[0];
    }
    k     = (*idum)/IQ1;
    *idum = IA1*(*idum - k*IQ1) - k*IR1;
    if(*idum < 0) *idum += IM1;
    k     = idum2/IQ2;
    idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
    if(idum2 < 0) idum2 += IM2;
    j     = iy/NDIV;
    iy    = iv[j] - idum2;
    iv[j] = *idum;
    if(iy < 1) iy += IMM1;
    if((temp = AM*iy) > RNMX) return RNMX;
    else return temp;g
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

// End: function ran2()

