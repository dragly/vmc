// Begin of main program

// TODO: Convert every vector to armadillo code ...

#include <iostream>
#include <stdlib.h>
#include "mainapplication.h"

using namespace std;

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

