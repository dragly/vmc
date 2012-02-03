// Variational Monte Carlo for atoms with up to two electrons

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <mpi.h>
#include <QApplication>

#include "wavestandard.h"
#include "wavesimple.h"
#include "montecarlostandard.h"
#include "matrix.h"
#include "hamiltonianstandard.h"
#include "hamiltoniansimple.h"
#include "minimizerstandard.h"
#include "mainwindow.h"
using namespace  std;


// Begin of main program

int main(int argc, char* argv[])
{
    //  MPI initializations
    MPI_Init (0, 0);
    int rank, nProcesses;
    MPI_Comm_size (MPI_COMM_WORLD, &nProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    int execVal;
    if(rank == 0) {
        QApplication a(argc, argv);
        MainWindow *mainWindow = new MainWindow(rank, nProcesses);
        mainWindow->showNormal();
        execVal = a.exec();
    } else {
        execVal = 0;
    }

    // End MPI
    MPI_Finalize ();

    return execVal;
}  //  end of main function

