#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QSettings>

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
    void on_nCyclesTxt_textChanged(const QString &arg1);

private:
    Ui::MainWindow *ui;

    QSettings *settings;
};

#endif // MAINWINDOW_H
