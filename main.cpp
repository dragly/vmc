// Begin of main program

#include <iostream>

#include "mainapplication.h"

using namespace std;

int main(int argc, char* argv[])
{
    MainApplication *app = new MainApplication(&argc, &argv);
    cout << "Running main application..." << endl;
    return app->runApplication();;
}  //  end of main function

