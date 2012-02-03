#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "minimizer.h"

#include <QMainWindow>

class MinimizerQObject;

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(int rank, int nProcesses, QWidget *parent = 0);
    ~MainWindow();

private slots:

private:
    Ui::MainWindow *ui;

    MinimizerQObject *minimizerQObject;
    QThread *minimizerThread;
    int rank;
    int nProcesses;
};

#endif // MAINWINDOW_H
