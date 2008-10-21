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

  // settings
  void quit();
  void readSettings();
  void writeSettings();

	// load, erase, duplicate
  void on_actionEraseAll_triggered();
  void on_actionLoadPolyhedron_triggered();
  bool on_actionErasePolyhedron_triggered();
  void on_actionDuplicatePolyhedron_triggered();

	// activate
  void on_actionActivatePolyhedron_triggered();
  void on_actionSetPolyhedronA_triggered();
  void on_actionSetPolyhedronB_triggered();

	// inside out
  void on_actionInsideOut_triggered();

  // defined in MainWindow_parameterization.cpp
  void on_actionMVC_triggered();
  void on_actionDCP_triggered();

	// save selected polyhedron as...
	void on_actionSaveAs_triggered(); 

  // defined in MainWindow_simplify.cpp
  void on_actionSimplify_triggered();

  // defined in MainWindow_convex_hull.cpp
  void on_actionConvexHull_triggered();

  // defined in MainWindow_kernel.cpp
  void on_actionKernel_triggered();

  // subdivision methods, in MainWindow_subdivision_methods.cpp
  void on_actionLoop_triggered();
  void on_actionCatmullClark_triggered();
  void on_actionSqrt3_triggered();

  // Boolean operations, in MainWindow_boolean_operations.cpp
  void on_actionUnion_triggered();
  void on_actionDifference_triggered();
  void on_actionIntersection_triggered();

  // curvature estimation, in MainWindow_curvature_estimation.cpp
  void on_actionEstimateCurvature_triggered();

  // PCA, in MainWindow_pca.cpp
  void on_actionFitLine_triggered();
  void on_actionFitPlane_triggered();

  // self intersection, in MainWindow_self_intersection.cpp
  void on_actionSelfIntersection_triggered();

  // remeshing, in MainWindow_remeshing.cpp
  void on_actionRemeshing_triggered();


protected:
  enum  Boolean_operation { BOOLEAN_UNION,
                            BOOLEAN_INTERSECTION,
                            BOOLEAN_DIFFERENCE };
  // defined in MainWindow_boolean_operations.cpp
  void boolean_operation(const Boolean_operation operation);

#ifdef CGAL_TAUCS_ENABLED
  enum  Parameterization_method  { PARAM_MVC,
                                   PARAM_DCP};
  // defined in MainWindow_parameterization.cpp
  void parameterize(const Parameterization_method method);
#endif

  void dragEnterEvent(QDragEnterEvent *event);
  void dropEvent(QDropEvent *event);
  void closeEvent(QCloseEvent *event);

  bool onePolygonIsSelected() const;
  int getSelectedPolygonIndex() const;

  Polyhedron* getSelectedPolygon();
private:
  QString strippedName(const QString &fullFileName);

  Scene* scene;
  Viewer* viewer;
  QTreeView* treeView;
  Ui::MainWindow* ui;
};

#endif // ifndef MAINWINDOW_H
