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
  void selectionChanged();
  void openRecentFile();
  void setCurrentFile(const QString &fileName);
  void updateRecentFileActions();

  // settings
  void quit();
  void readSettings();
  void writeSettings();

  void on_actionLoadPolyhedron_triggered();
  void on_actionErasePolyhedron_triggered();
  void on_actionDuplicatePolyhedron_triggered();

  // defined in MainWindow_simplify.cpp
  void on_actionSimplify_triggered();

  // defined in MainWindow_convex_hull.cpp
  void on_actionConvexHull_triggered();

  // defined in MainWindow_kernel.cpp
  void on_actionKernel_triggered();

  // subdivision methods are defined in MainWindow_subdivision_methods.cpp
  void on_actionLoop_triggered();
  void on_actionCatmullClark_triggered();
  void on_actionSqrt3_triggered();

  // PCA
  void on_actionFitPlane_triggered();
  void on_actionFitLine_triggered();

protected:
  void dragEnterEvent(QDragEnterEvent *event);
  void dropEvent(QDropEvent *event);
  void closeEvent(QCloseEvent *event);

  void selectPolyhedron(int i);
  bool onePolygonIsSelected() const;
  int getSelectedPolygonIndex() const;

  Polyhedron* getSelectedPolygon();
private:
  QString strippedName(const QString &fullFileName);
  QAction* recentFilesSeparator;

  Scene* scene;

  enum { MaxRecentFiles = 10 };
  QAction *recentFileActs[MaxRecentFiles];
};

#endif // ifndef MAINWINDOW_H
