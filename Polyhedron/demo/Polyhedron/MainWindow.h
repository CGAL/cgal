#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtOpenGL/qgl.h>
#include "ui_MainWindow.h"
#include <CGAL/Qt/DemosMainWindow.h>

class QDragEnterEvent;
class QDropEvent;
class Scene;
struct Polyhedron;

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
  void on_treeView_itemSelectionChanged();

	void on_actionConvexHull_triggered();
  void on_actionLoadPolyhedron_triggered();
  void on_actionErasePolyhedron_triggered();
  void on_actionDuplicatePolyhedron_triggered();

  // subdivision methods are defined in MainWindow_subdivision_methods.cpp
  void on_actionCatmullClark_triggered();
  void on_actionSqrt3_triggered();

protected:
  void dragEnterEvent(QDragEnterEvent *event);
  void dropEvent(QDropEvent *event);

  inline bool onePolygonIsSelected() const;
  inline int getSelectedPolygonIndex() const;

  Polyhedron* getSelectedPolygon();
private:
  Scene* scene;
};

bool MainWindow::onePolygonIsSelected() const
{
  return treeView->selectionModel()->selectedRows().size() == 1;
}

int MainWindow::getSelectedPolygonIndex() const
{
  QModelIndexList selectedRows = treeView->selectionModel()->selectedRows();
  if(selectedRows.empty())
    return -1;
  else
    return selectedRows.first().row();
}

#endif // ifndef MAINWINDOW_H
