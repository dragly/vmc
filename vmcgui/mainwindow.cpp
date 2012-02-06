#include <QSettings>

#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    settings = new QSettings("../config.ini", QSettings::IniFormat);
    ui->nCyclesTxt->setText(settings->value("MinimizerStandard/nCycles", 1000).toString());
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_nCyclesTxt_textChanged(const QString &arg1)
{
    settings->setValue("MinimizerStandard/nCycles", ui->nCyclesTxt->text());
}
