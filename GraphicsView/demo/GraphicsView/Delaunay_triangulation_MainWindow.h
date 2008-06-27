#ifndef DELAUNAY_TRIANGULATION_MAIN_WINDOW_H
#define DELAUNAY_TRIANGULATION_MAIN_WINDOW_H

#include <QtGui>
#include <QString>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include "ui_Delaunay_triangulation_MainWindow.h"
#include "DemosMainWindow.h"

// forward declarations
namespace CGAL {
  namespace Qt {
    template <class Delaunay> class TriangulationGraphicsItem;
    template <class Delaunay> class VoronoiGraphicsItem;
    template <class Delaunay> class TriangulationMovingPoint;
    template <class Delaunay> class TriangulationCircumcircle;
    template <class K> class GraphicsViewPolylineInput;
  } // namespace Qt
} // namespace CGAL

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;

typedef CGAL::Delaunay_triangulation_2<K> Delaunay;

class Delaunay_triangulation_MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Delaunay_triangulation_MainWindow
{
  Q_OBJECT
  
private:  
  Delaunay dt; 
  QGraphicsScene scene;  

  CGAL::Qt::TriangulationGraphicsItem<Delaunay> * dgi;
  CGAL::Qt::VoronoiGraphicsItem<Delaunay> * vgi;

  CGAL::Qt::TriangulationMovingPoint<Delaunay> * mp;
  CGAL::Qt::GraphicsViewPolylineInput<K> * pi;
  CGAL::Qt::TriangulationCircumcircle<Delaunay> *tcc;
public:
  Delaunay_triangulation_MainWindow();

public slots:

  void processInput(CGAL::Object o);

  void on_actionMovingPoint_toggled(bool checked);

  void on_actionShowDelaunay_toggled(bool checked);

  void on_actionShowVoronoi_toggled(bool checked);

  void on_actionInsertPoint_toggled(bool checked);
  
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

#endif // DELAUNAY_TRIANGULATION_MAIN_WINDOW_H

