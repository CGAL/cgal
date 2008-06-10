#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QtGui>
#include <QMainWindow>
#include <QLabel>
#include <QString>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include "QConstrainedTriangulation_2.h"
#include "QTriangulation_2.h"
#include "TriangulationMovingPoint_2.h"
#include "PolylineInput_2.h"
#include "TriangulationCircumcenter_2.h"
#include "TriangulationGraphicsItem_2.h"
#include "ConstrainedTriangulationGraphicsItem_2.h"
#include "VoronoiGraphicsItem_2.h"
#include "Navigation.h"
#include "Navigation2.h"

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

  CGAL::Navigation* navigation;
  CGAL::Navigation2* navigation2;

#ifdef DELAUNAY_VORONOI 
  CGAL::TriangulationGraphicsItem_2<Delaunay> * dgi; 
  CGAL::QTriangulation_2<Delaunay>  * sdt;
  CGAL::VoronoiGraphicsItem_2<Delaunay> * vgi;
#else
  CGAL::ConstrainedTriangulationGraphicsItem_2<Delaunay> * dgi;
  CGAL::QConstrainedTriangulation_2<Delaunay>  * sdt;
#endif

  QLabel* xycoord ;

  CGAL::TriangulationMovingPoint_2<Delaunay> * mp;
  CGAL::PolylineInput_2<K> * pi;
  CGAL::TriangulationCircumcenter_2<Delaunay> *tcc;
public:
  MainWindow();

private:
  void connectActions();
  void setupStatusBar();

public slots:

  void movingPoint(bool checked);

  void showDelaunay(bool checked);

  void on_actionShowVoronoi_toggled(bool checked);

  void insertPolyline(bool checked);
  
  void circumcenter(bool checked);

  void clear();

  void loadConstraints();

  void loadConstraints(QString);

  void saveConstraints();

  void saveConstraints(QString);

  void insertRandomPoints();

  void about();

  void aboutCGAL();

  void updateMouseCoordinates(QString s);
  
};

#endif // MAIN_WINDOW_H

