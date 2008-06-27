#ifndef CONSTRAINED_DELAUNAY_TRIANGULATION_MAIN_WINDOW_H
#define CONSTRAINED_DELAUNAY_TRIANGULATION_MAIN_WINDOW_H

#include <QtGui>
#include <QString>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include "ui_Constrained_Delaunay_triangulation_MainWindow.h"
#include "DemosMainWindow.h"

// forward declarations
namespace CGAL {
  namespace Qt {
    template <class Delaunay> class ConstrainedTriangulationGraphicsItem;
    template <class Delaunay> class TriangulationMovingPoint;
    template <class Delaunay> class TriangulationCircumcircle;
    template <class K> class GraphicsViewPolylineInput;
  } // namespace Qt
} // namespace CGAL

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Triangulation_vertex_base_2<K>  Vertex_base;
typedef CGAL::Constrained_triangulation_face_base_2<K> Face_base;
typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base>  TDS;
typedef CGAL::Exact_predicates_tag              Itag;


typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> Delaunay;

class Constrained_Delaunay_triangulation_MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Constrained_Delaunay_triangulation_MainWindow
{
  Q_OBJECT
  
private:  
  Delaunay dt; 
  QGraphicsScene scene;  

  CGAL::Qt::ConstrainedTriangulationGraphicsItem<Delaunay> * dgi;

  CGAL::Qt::TriangulationMovingPoint<Delaunay> * mp;
  CGAL::Qt::GraphicsViewPolylineInput<K> * pi;
  CGAL::Qt::TriangulationCircumcircle<Delaunay> *tcc;
public:
  Constrained_Delaunay_triangulation_MainWindow();

private:
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

signals:
  void changed();
};

#endif // CONSTRAINED_DELAUNAY_TRIANGULATION_MAIN_WINDOW_H

