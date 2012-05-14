// Begin of main program

// TODO: Convert every vector to armadillo code ...

#include <iostream>

#include "mainapplication.h"

using namespace std;

int main(int argc, char* argv[])
{
    cout << "Starting application" << endl;
    MainApplication *app = new MainApplication(&argc, &argv);
    cout << "Running config" << endl;
    app->runConfiguration();
    cout << "Finalizing" << endl;
    app->finalize();
    return 0;
}  //  end of main function

