#ifndef __VIEWER_H__
#define __VIEWER_H__

#include <qwidget.h>
#include <qmainwindow.h>
#include <qmenubar.h>
#include <qscrollview.h>
#include <qtoolbar.h>
#include <qpushbutton.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qstatusbar.h>
#include <list>
//#include <iostream>
//#include <math.h>
#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
// // #include <CGAL/Triangulation_euclidean_traits_2.h>
// #include <CGAL/Constrained_triangulation_2.h>
// #include <CGAL/Constrained_Delaunay_triangulation_2.h>
// #include <CGAL/Constrained_triangulation_face_base_2.h> 
// #include <CGAL/triangulation_assertions.h>
// #include <CGAL/Triangulation_short_names_2.h>
// #include <CGAL/Arithmetic_filter.h>
//#include <CGAL/double.h>
//#include <LEDA/real.h>
// #include "Mesh.h"


class TrViewer;

class TrFrame : public QMainWindow {
  Q_OBJECT  // Qt object
    
 public:
  TrFrame();
  void setStatus(int x, int y);
 protected:
  QScrollView *scr;
  TrViewer *trv;
  QLineEdit *editW, *editH;
  QLabel *lblStatus;
 private slots:
   void clearTriangulation();
   void openTriangulation();
   void saveTriangulation();
   void onChangeSizes();
   void mesh();
 public:
   friend class TrViewer;
};

// struct Point {
//   int x, y;
//   Point(int a, int b):x(a), y(b){}
// };

// struct Line {
//   int x1,y1,x2,y2;
//   Line(int a, int b, int c, int d):x1(a), y1(b), x2(c), y2(d) {}
// };

#ifdef CGAL_MESH_H

#include <CGAL/intersections.h>

typedef CGAL::Simple_cartesian<double>  K1;
typedef CGAL::Filtered_kernel<K1>       Kernel;
struct K : public Kernel {}; 

typedef K::Point_2 Point;
typedef K::Line_2 Line;

typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds,
  CGAL::Exact_predicates_tag> Tr;
//typedef CGAL::Constrained_triangulation_2<K, Tds> Tr;

typedef CGAL::Mesh_2<Tr> Msh;

class TrViewer : public QWidget {
  Q_OBJECT
 public:
  TrViewer(TrFrame *f, QWidget *parent);
 protected:
  list<Point> points;
  list<Line> lines;
  Msh *mesh;
  //  bool auto_mesh;
  TrFrame *frame;
  void paintEvent(QPaintEvent *pe);
  void mouseMoveEvent(QMouseEvent *me);
  bool dragging;
  int startx, starty, endx, endy, oldx, oldy;
  void mousePressEvent(QMouseEvent *me);
  void mouseReleaseEvent(QMouseEvent *me);
  void keyPressEvent(QKeyEvent *ke);
 public slots:
  void onMesh();
  //  void toggleAutoMesh(bool);
 public:
  friend class TrFrame;
};

#endif

#endif
