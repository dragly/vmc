#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "minimizerstandard.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    minimizer = new MinimizerStandard(rank, nProcesses);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_runButton_clicked()
{

}
