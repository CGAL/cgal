#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QtGui>
#include <QMainWindow>
#include <QLabel>
#include <QString>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include "QTriangulationMovingPoint_2.h"
#include "QPolylineInput_2.h"
#include "QTriangulationCircumcenter_2.h"
#include "QTriangulationGraphicsItem_2.h"
#include "QConstrainedTriangulationGraphicsItem_2.h"
#include "QVoronoiGraphicsItem_2.h"
#include "QNavigation.h"
#include "QNavigation2.h"

#include "ui_MainWindow.h"


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

  CGAL::QNavigation* navigation;
  CGAL::QNavigation2* navigation2;

#ifdef DELAUNAY_VORONOI 
  CGAL::QTriangulationGraphicsItem_2<Delaunay> * dgi; 
  CGAL::QVoronoiGraphicsItem_2<Delaunay> * vgi;
#else
  CGAL::QConstrainedTriangulationGraphicsItem_2<Delaunay> * dgi;
#endif

  QLabel* xycoord ;

  CGAL::QTriangulationMovingPoint_2<Delaunay> * mp;
  CGAL::QPolylineInput_2<K> * pi;
  CGAL::QTriangulationCircumcenter_2<Delaunay> *tcc;
public:
  MainWindow();

private:
  void connectActions();
  void setupStatusBar();

public slots:

  void process(CGAL::Object o);

  void on_actionMovingPoint_toggled(bool checked);

  void on_actionShowDelaunay_toggled(bool checked);

  void on_actionShowVoronoi_toggled(bool checked);

  void on_actionInsertPolyline_toggled(bool checked);
  
  void on_actionCircumcenter_toggled(bool checked);

  void on_actionClear_triggered();

  void on_actionLoadConstraints_triggered();

  void loadConstraints(QString);

  void on_actionSaveConstraints_triggered();

  void saveConstraints(QString);

  void on_actionInsertRandomPoints_triggered();

  void on_actionAbout_triggered();

  void on_actionAboutCGAL_triggered();

  void updateMouseCoordinates(QString s);
  
  signals:

  void changed();
};

#endif // MAIN_WINDOW_H

