#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "minimizer.h"

#include <QMainWindow>

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_runButton_clicked();

private:
    Ui::MainWindow *ui;

    Minimizer *minimizer;
};

#endif // MAINWINDOW_H
