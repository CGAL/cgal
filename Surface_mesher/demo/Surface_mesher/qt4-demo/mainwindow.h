#ifndef _MAINWINDOW_H
#define _MAINWINDOW_H

#include "ui_mainwindow.h"

#include <QMainWindow>

class QDragEnterEvent;
class QDropEvent;
class Surface;
class QGLViewer;

class MainWindow : public  QMainWindow, private Ui::MainWindow 
{
  Q_OBJECT
public:
  MainWindow();
  void dragEnterEvent(QDragEnterEvent *);
  void dropEvent(QDropEvent *event);
  void surface_open(const QString& filename);
  
private slots:
  void on_action_Open_triggered();
  void on_action_Quit_triggered();
  void on_action_Options_triggered();
  
signals:
  void new_sharp_edges_angle_bounds(double, double);

private:
  Surface* surface;
  QGLViewer* viewer_ptr;
  double sharp_edges_angle_lower_bound;
  double sharp_edges_angle_upper_bound;
};


#endif // _MAINWINDOW_H
