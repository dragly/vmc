#include <QFuture>
#include <QtConcurrentRun>

#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "minimizerstandard.h"
#include "minimizerqobject.h"

MainWindow::MainWindow(int rank, int nProcesses, QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    rank(rank),
    nProcesses(nProcesses)
{
    ui->setupUi(this);

    minimizerQObject = new MinimizerQObject(rank, nProcesses);
    minimizerThread = new QThread(this);
    minimizerQObject->moveToThread(minimizerThread);
    minimizerThread->start();

    connect(ui->runButton, SIGNAL(clicked()), minimizerQObject, SLOT(runMinimizer()));
}

MainWindow::~MainWindow()
{
    delete ui;
}
