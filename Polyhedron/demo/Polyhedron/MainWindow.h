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

  void on_actionConvexHull_triggered();
  void on_actionLoadPolyhedron_triggered();
  void on_actionErasePolyhedron_triggered();
  void on_actionDuplicatePolyhedron_triggered();
	void on_actionSimplify_triggered();

  // subdivision methods are defined in MainWindow_subdivision_methods.cpp
  void on_actionCatmullClark_triggered();
  void on_actionSqrt3_triggered();

protected:
  void dragEnterEvent(QDragEnterEvent *event);
  void dropEvent(QDropEvent *event);

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
