#ifndef CGAL_QT_DEMOS_MAIN_WINDOW_H
#define CGAL_QT_DEMOS_MAIN_WINDOW_H

#include <QMainWindow>

// forward declaration
class QLabel;
class QGraphicsView;
class QAction;
class QMenu;

namespace CGAL {
namespace Qt {

// forward declaration
class GraphicsViewNavigation;

class DemosMainWindow : public QMainWindow {

  Q_OBJECT

private:
  QMenu* getMenu(QString objectName, QString title);
  void popupAboutBox(QString title, QString html_resource_name);
  QMenu* getHelpMenu();

protected:
  DemosMainWindow (QWidget * parent = 0, ::Qt::WindowFlags flags = 0 );
  void setupStatusBar();
  void addNavigation(QGraphicsView*);
  void setupOptionsMenu(QMenu* menu  = NULL);
  void addAboutCGAL(QMenu* menu  = NULL);
  void addAboutDemo(QString htmlResourceName, QMenu* menu  = NULL);

protected slots:
  void setUseAntialiasing(bool checked);
  void setUseOpenGL(bool checked);
  void popupAboutCGAL();
  void popupAboutDemo();
protected:
  QGraphicsView* view;
  GraphicsViewNavigation* navigation;
  QLabel* xycoord ;

  QAction *actionUse_OpenGL;
  QAction *actionUse_Antialiasing;
  QAction *actionAbout;
  QAction *actionAboutCGAL;

  QString aboutHtmlResource;
}; // end class DemosMainWindow
 
} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_DEMOS_MAIN_WINDOW_H
