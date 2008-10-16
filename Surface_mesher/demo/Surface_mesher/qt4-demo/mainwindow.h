#ifndef _MAINWINDOW_H
#define _MAINWINDOW_H

#include <QMainWindow>
#include "ui_mainwindow.h"
#include <CGAL/Qt/DemosMainWindow.h>

class QDragEnterEvent;
class QDropEvent;
class Surface;
class QGLViewer;
class QDoubleSpinBox;
class QCloseEvent;

class MainWindow : public CGAL::Qt::DemosMainWindow, public Ui::MainWindow
{
  Q_OBJECT
public:
  MainWindow(MainWindow* other_window = 0);
  void dragEnterEvent(QDragEnterEvent *);
  void dropEvent(QDropEvent *event);

public slots:
  void show_only(QString);
  void surface_open(const QString& filename);

private slots:
  void on_action_Open_triggered();
  void on_action_Quit_triggered();
  void on_action_Clone_triggered();
  
private:
  void closeEvent(QCloseEvent *event);
  Surface* surface;
  double sharp_edges_angle_lower_bound;
  double sharp_edges_angle_upper_bound;
  QDoubleSpinBox* spinbox_isovalue;
  QDoubleSpinBox* spinbox_radius_bound;
  QDoubleSpinBox* spinbox_distance_bound;
};


#endif // _MAINWINDOW_H
