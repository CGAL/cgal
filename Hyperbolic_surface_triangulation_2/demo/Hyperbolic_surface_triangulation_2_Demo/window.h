#ifndef CGAL_HYPERBOLIC_SURFACE_TRIANGULATION_2_DEMO_WINDOW
#define CGAL_HYPERBOLIC_SURFACE_TRIANGULATION_2_DEMO_WINDOW

// Qt headers
#include <CGAL/Qt/utility.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/DemosMainWindow.h>

// UI generated header
#include "ui_drawing_window_description.h"

// CGAL headers
#include <CGAL/Gmpq.h>
#include <CGAL/Complex_without_sqrt.h>
#include <CGAL/Hyperbolic_surface_triangulation_2.h>
#include <CGAL/Hyperbolic_surfaces_traits_2.h>

typedef CGAL::Complex_without_sqrt<CGAL::Gmpq>                                Complex;
typedef CGAL::Hyperbolic_surfaces_traits_2<CGAL::Gmpq>                        Traits;
typedef CGAL::Hyperbolic_surface_triangulation_2<Traits>                      Triangulation;

class DemoWindowItem :
    public CGAL::Qt::GraphicsItem
{
    Q_OBJECT // Qt macro for Qt objects
    // (Q_OBJECT does not support templates)
private:
  typedef CGAL::Bbox_2                                Bbox_2; // "Bounding box" : just a box type used for drawing

  // Edges to draw
  std::vector<std::pair<Complex,Complex>>             _edges;

  // Pens for drawing
  QPen                                                _poincare_disk_pen;
  QPen                                                _edges_pen;

  // radius of the poincaré disk
  const int _poincare_disk_radius_in_pixels = 600;
  // Approximation treshold : used to decide when to simplify a computation (ex : draw a line instead of an arc if an hyperbolic segment is very small)
  const double computation_treshold = 0.001;
  const double computation_treshold_squared = computation_treshold*computation_treshold;

public:
  // Initializer
  DemoWindowItem();

  // Qt methods
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  QRectF boundingRect() const;
  void modelChanged();

  // Drawing method
  void draw_triangulation(Triangulation& triangulation);

private:
  // Sub-methods for drawing edges and vertices
  void draw_point(QPainter* painter, Complex position);

  void draw_edge(QPainter* painter, Complex source, Complex target);
  void draw_line(QPainter* painter, double point_1_x, double point_1_y, double point_2_x, double point_2_y);
  void draw_arc(QPainter* painter, double point_1_x, double point_1_y, double point_2_x, double point_2_y, double center_x, double center_y);

  double deg_angle(double x, double y);
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DemoWindow :
    public CGAL::Qt::DemosMainWindow, public Ui::MainWindow
{
    Q_OBJECT // Qt macro for Qt objects
    // (Q_OBJECT does not support templates)
private:
  QGraphicsScene                  _scene;
  DemoWindowItem*                 _item;

public:
  DemoWindow();
  DemoWindowItem& item();

  // Events handling
  void keyPressEvent(QKeyEvent* event);
};

#endif // CGAL_HYPERBOLIC_SURFACE_TRIANGULATION_2_DEMO_WINDOW
