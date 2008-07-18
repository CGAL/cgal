#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtOpenGL/qgl.h>
#include <CGAL/Qt/DemosMainWindow.h>

class QDragEnterEvent;
class QDropEvent;
class Scene;
class Viewer;
class QTreeView;
namespace Ui {
  class MainWindow;
}

#include "Polyhedron_type_fwd.h"

class MainWindow : 
  public CGAL::Qt::DemosMainWindow
{
  Q_OBJECT
public:
  MainWindow(QWidget* parent = 0);
  ~MainWindow();

public slots:
  void updateViewerBBox();
  void open(QString filename);

  void selectPolyhedron(int i);

  void setAddKeyFrameKeyboardModifiers(Qt::KeyboardModifiers);

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
  bool on_actionErasePolyhedron_triggered();
  void on_actionEraseAll_triggered();
  void on_actionDuplicatePolyhedron_triggered();

  void on_actionActivatePolyhedron_triggered();
  void on_actionSetPolyhedronA_triggered();
  void on_actionSetPolyhedronB_triggered();

  // save
  // TODO: save all, save current (do we store the current file name?)
  void on_actionSaveAs_triggered(); // save selected polyhedron as...

  // merge (TODO)

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

  // Boolean operations defined in MainWindow_boolean_operations.cpp
  void on_actionUnion_triggered();
  void on_actionIntersection_triggered();
  void on_actionDifference_triggered();

  // curvature estimation, in MainWindow_curvature_estimation.cpp
  void on_actionEstimateCurvature_triggered();

  // PCA, in MainWindow_pca.cpp
  void on_actionFitPlane_triggered();
  void on_actionFitLine_triggered();

  // self intersection, in MainWindow_self_intersection.cpp
  void on_actionSelfIntersection_triggered();

protected:
  enum  Boolean_operation { BOOLEAN_UNION,
                            BOOLEAN_INTERSECTION,
                            BOOLEAN_DIFFERENCE };
  // define in MainWindow_boolean_operations.cpp
  void boolean_operation(const Boolean_operation operation);

  void dragEnterEvent(QDragEnterEvent *event);
  void dropEvent(QDropEvent *event);
  void closeEvent(QCloseEvent *event);

  bool onePolygonIsSelected() const;
  int getSelectedPolygonIndex() const;

  Polyhedron* getSelectedPolygon();
private:
  QString strippedName(const QString &fullFileName);
  QAction* recentFilesSeparator;

  Scene* scene;
  Viewer* viewer;
  QTreeView* treeView;
  Ui::MainWindow* ui;

  enum { MaxRecentFiles = 10 };
  QAction *recentFileActs[MaxRecentFiles];
};

#endif // ifndef MAINWINDOW_H
