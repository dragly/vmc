#include "onerun.h"

#include "../inih/ini.h"

#include "../wavefunction/wavefunction.h"
#include "../hamiltonian/hamiltonian.h"
#include "../montecarlo/montecarlo.h"
#include "../config.h"
#include "../blocker.h"

// disable annoying unused parameter warnings from the MPI library which we don't have any control over
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>

OneRun::OneRun(Config *config_) :
    myRank(config_->myRank()),
    m_nProcesses(config_->m_nProcesses()),
    wave(config_->wave()),
    monteCarlo(config_->monteCarlo()),
    hamiltonian(config_->hamiltonian()),
    config(config_),
    onlyBlocking(false)
{
    monteCarlo->setStoreEnergies(true);
    blocker = new Blocker(config);
}

OneRun::~OneRun() {
    delete blocker;
}

void OneRun::loadConfiguration(INIParser *settings) {
    blocker->loadConfiguration(settings);
    nSamples = settings->GetDouble("OneRun","nCycles", 1000);
    alpha = settings->GetDouble("OneRun","alpha", 0);
    beta = settings->GetDouble("OneRun","beta", 0);
    onlyBlocking = settings->GetBoolean("OneRun","onlyBlocking", onlyBlocking);
    scratchDir = settings->GetString("Blocking", "scratchDir", "/scratch/blocking");
}

void OneRun::run() {
    if(!onlyBlocking) {
        double timeStart;
        double timeEnd;
        double totalTime;

        //    WaveIdeal *wave = new WaveIdeal(nParticles, config->nDimensions());
        //    HamiltonianIdeal *hamiltonian = new HamiltonianIdeal(nParticles, config->nDimensions(), charge);


        timeStart = MPI_Wtime();
        // broadcast the total number of  variations
        MPI_Bcast (&nSamples, 1, MPI_INT, 0, MPI_COMM_WORLD);

        double totalCumulativeEnergy;
        double totalCumulativeEnergySquared;
        double cumulativeEnergy;
        double cumulativeEnergySquared;
        double nTotalSamples = nSamples*m_nProcesses;

        double parameters[2];

        parameters[0] = alpha;
        parameters[1] = beta;

        wave->setParameters(parameters);
        monteCarlo->setOutputEnergies(true);
        monteCarlo->sample(nSamples);

        cumulativeEnergy = monteCarlo->energy();
        cumulativeEnergySquared = monteCarlo->energySquared();
        allEnergies = monteCarlo->allEnergies();

        double totalCumulativeEnergyDummy;
        MPI_Reduce(&cumulativeEnergy, &totalCumulativeEnergyDummy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        totalCumulativeEnergy = totalCumulativeEnergyDummy / m_nProcesses;

        double totalCumulativeEnergySquaredDummy;
        MPI_Reduce(&cumulativeEnergySquared, &totalCumulativeEnergySquaredDummy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        totalCumulativeEnergySquared = totalCumulativeEnergySquaredDummy / m_nProcesses;

        timeEnd = MPI_Wtime();
        totalTime = timeEnd-timeStart;

        // Print out results
        if ( myRank == 0) {
            std::cout << "Writing data to file" << std::endl;
            // output file as global variable
            ofstream dataFile;

            double energy = totalCumulativeEnergy ;
            double variance = totalCumulativeEnergySquared - energy*energy;
            double error = sqrt(variance / (nTotalSamples-1));
            dataFile.open("onerun.dat");
            dataFile << setprecision(20) << config->omega() << " ";
            dataFile << setprecision(20) << alpha << " ";
            dataFile << setprecision(20) << beta << " ";
            dataFile << setprecision(20) << energy << " ";
            dataFile << setprecision(20) << variance << " ";
            dataFile << setprecision(20) << error << std::endl;

            dataFile.close();


            ofstream energyPartsFile;
            energyPartsFile.open("energyparts.dat");
            energyPartsFile << hamiltonian->totalsString() << std::endl;
            energyPartsFile.close();

        }
        writeBlockData();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if ( myRank == 0) {
        std::cout << "Starting blocking " << std::endl;
        blocker->runBlocking();
    }
}

void OneRun::writeBlockData() {
    std::cout << "Writing blocking data" << std::endl;
    ofstream blockofile;
    // Setting output file name for this myRank:
    ostringstream fileName;
    ostringstream path;
    path << scratchDir << "/" << config->nParticles() << "p-omega" << config->omega();

    if(myRank == 0) {
        struct stat st;
        if(stat(scratchDir.c_str(), &st) != 0) {
            if(mkdir(scratchDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0) {
                std::cerr << "Error creating directory " << scratchDir << std::endl;
                exit(948);
            }
        }
        if(stat(path.str().c_str(), &st) != 0) {
            if(mkdir(path.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0) {
                std::cerr << "Error creating directory " << path.str() << std::endl;
                exit(948);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    fileName << path.str() << "/blocks_rank" << myRank << ".dat";
    // Open file for writing:
    blockofile.open(fileName.str().c_str(), ios::out | ios::binary);
    if(blockofile.is_open()) {
        blockofile.write((char*)(allEnergies+1), nSamples*sizeof(double));
        blockofile.close();
    } else {
        std::cerr << "Error opening file " << fileName.str() << std::endl;
        exit(948);
    }
}
