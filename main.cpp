// Begin of main program

// TODO: Convert every vector to armadillo code ...

#include <iostream>

#include "mainapplication.h"

using namespace std;

int main(int argc, char* argv[])
{
    MainApplication *app = new MainApplication(&argc, &argv);
    app->runConfiguration();
    app->finalize();
    return 0;
}  //  end of main function

