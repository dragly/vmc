/*!
  * \mainpage FYS4411 Diffusion and Variational Monte Carlo
  *
  * These pages contain the documentation for the code used in the FYS4411 project for variational and diffusion Monte Carlo calculations on quantum dots.
  */

// Begin of main program

// TODO: Convert every vector to armadillo code ...

#include <iostream>
#include <stdlib.h>
#include "mainapplication.h"

#include <armadillo>

using namespace std;

using namespace arma;

int main(int argc, char* argv[])
{
    MainApplication *app = new MainApplication(&argc, &argv);
    app->runConfiguration();
    app->finalize();

    // A little output to help notice when the program is finished executing
    int ret = system("kdialog --title 'vmc' --passivepopup 'Execution finished' 3");
    (void)ret;
    return 0;
}  //  end of main function

