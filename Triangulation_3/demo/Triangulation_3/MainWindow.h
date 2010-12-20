#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include "ui_MainWindow.h"
#include <CGAL/Qt/DemosMainWindow.h>
#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>

#include <QtGui>

#include "Scene.h"

class QWidget;

class MainWindow : public CGAL::Qt::DemosMainWindow, private Ui::MainWindow
{
  Q_OBJECT

public:
  MainWindow(QWidget* parent = 0);
  ~MainWindow() {}

public slots:
  // file menu
  void on_actionLoad_Points_triggered();
  void on_actionSave_Points_triggered();

  // edit menu
  void on_actionGenerate_Points_triggered();
  void stopAnimation();

  // mode menu
  void setMode(QAction *a);

  // show menu
  void on_actionClear_Scene_triggered();

  // about menu
  void popupAboutCGAL();

  signals:
  void sceneChanged();

protected:
  void closeEvent(QCloseEvent *event);

private:
  void connectActions();

private:
  Scene m_scene;
};

#endif
