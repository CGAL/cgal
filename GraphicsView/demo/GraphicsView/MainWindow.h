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
  class QtNavigation;
  template <class Delaunay> class QtTriangulationGraphicsItem;
  template <class Delaunay> class QtVoronoiGraphicsItem;
  template <class Delaunay> class QtConstrainedTriangulationGraphicsItem;
  template <class Delaunay> class QTriangulationMovingPoint_2;
  template <class Delaunay> class QTriangulationCircumcenter_2;
  template <class K> class QtPolylineInput;
} // end namespace CGAL

class QLabel;
class QWidget;
class QGLWidget;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Triangulation_vertex_base_2<K>  Vertex_base;
typedef CGAL::Constrained_triangulation_face_base_2<K> Face_base;
typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base>  TDS;
typedef CGAL::Exact_predicates_tag              Itag;


#define DELAUNAY_VORONOI

#ifdef DELAUNAY_VORONOI
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;
#else 
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> Delaunay;
#endif

class MainWindow : public QMainWindow, public Ui::MainWindow
{
  Q_OBJECT
  
private:  
  Delaunay dt; 
  QGraphicsScene scene;  

  CGAL::QtNavigation* navigation;

#ifdef DELAUNAY_VORONOI 
  CGAL::QtTriangulationGraphicsItem<Delaunay> * dgi; 
  CGAL::QtVoronoiGraphicsItem<Delaunay> * vgi;
#else
  CGAL::QtConstrainedTriangulationGraphicsItem<Delaunay> * dgi;
#endif

  QLabel* xycoord ;

  CGAL::QTriangulationMovingPoint_2<Delaunay> * mp;
  CGAL::QtPolylineInput<K> * pi;
  CGAL::QTriangulationCircumcenter_2<Delaunay> *tcc;
public:
  MainWindow();

private:
  void connectActions();
  void setupStatusBar();

#ifndef DELAUNAY_VORONOI 
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

#endif // ifndef DELAUNAY_VORONOI
  
public slots:

  void process(CGAL::Object o);

  void on_actionUse_OpenGL_toggled(bool checked);

  void on_actionMovingPoint_toggled(bool checked);

  void on_actionShowDelaunay_toggled(bool checked);

  void on_actionShowVoronoi_toggled(bool checked);

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

