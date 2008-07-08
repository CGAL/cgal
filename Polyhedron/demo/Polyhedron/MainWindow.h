#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "ui_MainWindow.h"
#include <CGAL/Qt/DemosMainWindow.h>

class QDragEnterEvent;
class QDropEvent;
class Scene;

class MainWindow : 
  public CGAL::Qt::DemosMainWindow,
  public Ui::MainWindow
{
  Q_OBJECT
public:
  MainWindow(QWidget* parent = 0);

public slots:
  void updateViewerBBox();
  void open(QString filename);

protected slots:
  void on_actionLoadPolyhedron_triggered();
  void on_actionErasePolyhedron_triggered();
  void on_actionDuplicatePolyhedron_triggered();

protected:
  void dragEnterEvent(QDragEnterEvent *event);
  void dropEvent(QDropEvent *event);

private:
  Scene* scene;
};

#endif // ifndef MAINWINDOW_H
