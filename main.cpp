// Begin of main program

#include <iostream>

#include "mainapplication.h"

using namespace std;

int main(int argc, char* argv[])
{
    cout << "Running main application..." << endl;
    MainApplication *app = new MainApplication(&argc, &argv);

    app->loadConfiguration();
    app->runConfiguration();
    app->finalize();
    return 0;
}  //  end of main function

