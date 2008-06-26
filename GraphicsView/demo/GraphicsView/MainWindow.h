#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QtGui>
#include <QMainWindow>
#include <QString>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include "ui_MainWindow.h"

namespace CGAL {
  namespace Qt {
    class GraphicsViewNavigation;
    template <class Delaunay> class TriangulationGraphicsItem;
    template <class Delaunay> class ConstrainedTriangulationGraphicsItem;
    template <class Delaunay> class TriangulationMovingPoint;
    template <class Delaunay> class TriangulationCircumcircle;
    template <class K> class GraphicsViewPolylineInput;
  } // namespace Qt
} // namespace CGAL

class QLabel;
class QWidget;
class QGLWidget;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Triangulation_vertex_base_2<K>  Vertex_base;
typedef CGAL::Constrained_triangulation_face_base_2<K> Face_base;
typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base>  TDS;
typedef CGAL::Exact_predicates_tag              Itag;


typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> Delaunay;

class MainWindow : public QMainWindow, public Ui::MainWindow
{
  Q_OBJECT
  
private:  
  Delaunay dt; 
  QGraphicsScene scene;  

  CGAL::Qt::GraphicsViewNavigation* navigation;

  CGAL::Qt::ConstrainedTriangulationGraphicsItem<Delaunay> * dgi;

  QLabel* xycoord ;

  CGAL::Qt::TriangulationMovingPoint<Delaunay> * mp;
  CGAL::Qt::GraphicsViewPolylineInput<K> * pi;
  CGAL::Qt::TriangulationCircumcircle<Delaunay> *tcc;
public:
  MainWindow();

private:
  void setupStatusBar();

  template <typename Iterator> 
  void insert_polyline(Iterator b, Iterator e)
  {
    Point_2 p, q;
    Delaunay::Vertex_handle vh, wh;
    Iterator it = b;
    vh = dt.insert(*it);
    p = *it;
    ++it;
    for(; it != e; ++it){
      q = *it;
      if(p != q){
        wh = dt.insert(*it);
        dt.insert_constraint(vh,wh);
        vh = wh;
        p = q;
      } else {
        std::cout << "duplicate point: " << p << std::endl; 
      }
    }
    emit(changed());
  }

public slots:

  void processInput(CGAL::Object o);

  void on_actionUse_OpenGL_toggled(bool checked);

  void on_actionMovingPoint_toggled(bool checked);

  void on_actionShowDelaunay_toggled(bool checked);

  void on_actionInsertPolyline_toggled(bool checked);
  
  void on_actionCircumcenter_toggled(bool checked);

  void on_actionClear_triggered();

  void on_actionRecenter_triggered();

  void on_actionLoadConstraints_triggered();

  void loadConstraints(QString);

  void on_actionSaveConstraints_triggered();

  void saveConstraints(QString);

  void on_actionInsertRandomPoints_triggered();

  void on_actionAbout_triggered();

  void on_actionAboutCGAL_triggered();

  signals:

  void changed();
};

#endif // MAIN_WINDOW_H

