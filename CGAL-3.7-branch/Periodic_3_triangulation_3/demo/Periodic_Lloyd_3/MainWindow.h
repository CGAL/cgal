#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include "ui_MainWindow.h"
#include <CGAL/Qt/DemosMainWindow.h>
#include <QFileDialog>
#include <QInputDialog>
#include <QSlider>
#include <QTimer>

#include <QtGui>
#include <QProcess>

class QWidget;

class MainWindow : public CGAL::Qt::DemosMainWindow, private Ui::MainWindow
{
  Q_OBJECT

public:
  MainWindow(QWidget* = 0);
  ~MainWindow() { process->close(); delete(process); }

  void connectActions();

  Scene scene;
  QTimer * qtimer;

private:
  QProcess* process;

public slots:
  void newPoints(int i);
  void newPointSet();
  void loadPoints();
  void savePoints();
  void speedChanged(int i);
  void togglePause(bool p);
  void toggle8Copies(bool on);
  void toggle2D(bool on);
  void lloydStep();
  void help();

  signals:
  void sceneChanged();
  void speedChanged();
};




#endif
