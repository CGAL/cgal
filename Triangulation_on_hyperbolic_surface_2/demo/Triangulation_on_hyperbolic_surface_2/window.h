// Copyright (c) 2024
// INRIA Nancy (France), and Université Gustave Eiffel Marne-la-Vallee (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Vincent Despré, Loïc Dubois, Monique Teillaud

#ifndef CGAL_TRIANGULATION_ON_HYPERBOLIC_SURFACE_DEMO_WINDOW
#define CGAL_TRIANGULATION_ON_HYPERBOLIC_SURFACE_DEMO_WINDOW

// Qt headers
#include <CGAL/Qt/utility.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/DemosMainWindow.h>

#include <QSvgGenerator>

// UI generated header
#include "ui_drawing_window_description.h"

#include <CGAL/Exact_rational.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>
#include <CGAL/Hyperbolic_surface_traits_2.h>
#include <CGAL/Triangulation_on_hyperbolic_surface_2.h>
#include <CGAL/Delaunay_triangulation_on_hyperbolic_surface_2.h>
#include <CGAL/Hyperbolic_Dirichlet_domain_2.h>

typedef CGAL::Exact_rational                                                      NumberType;
// typedef CGAL::Simple_cartesian<NumberType>                                        Kernel;
typedef CGAL::Circular_kernel_2<CGAL::Simple_cartesian<NumberType>,CGAL::Algebraic_kernel_for_circles_2_2<NumberType>> Kernel;
typedef CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<Kernel>                  ParentTraits;
typedef CGAL::Hyperbolic_surface_traits_2<ParentTraits>                           Traits;
typedef Traits::Hyperbolic_point_2                                                Point;
typedef typename Traits::Hyperbolic_Voronoi_point_2                               Voronoi_point;
typedef CGAL::Hyperbolic_fundamental_domain_2<Traits>                             Domain;
typedef CGAL::Triangulation_on_hyperbolic_surface_2<Traits>                       Triangulation;
typedef CGAL::Delaunay_triangulation_on_hyperbolic_surface_2<Traits>              Delaunay_triangulation;
typedef typename Delaunay_triangulation::Anchor                                   Anchor;

class DemoWindowItem
  : public CGAL::Qt::GraphicsItem
{
  Q_OBJECT // Qt macro for Qt objects
  // (Q_OBJECT does not support templates)

private:
  typedef CGAL::Bbox_2                                Bbox_2; // "Bounding box": just a box type used for drawing

  // Edges to draw
  std::vector<std::pair<Point,Point> > edges_;
  std::vector<std::pair<Point,Point> > dirichlet_edges_;

  // Pens for drawing
  QPen poincare_disk_pen_;
  QPen edges_pen_;
  QPen dirichlet_pen_;

  // radius of the poincaré disk
  const int poincare_disk_radius_in_pixels_ = 600;

  // Approximation threshold: used to decide when to simplify a computation (ex: draw a line
  // instead of an arc if an hyperbolic segment is very small)
  const double computation_threshold = 0.001;
  const double computation_threshold_squared = computation_threshold*computation_threshold;

public:
  // Initializer
  DemoWindowItem();

  // Qt methods
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  QRectF boundingRect() const;
  void modelChanged();

  // Drawing methods
  void draw_triangulation(Triangulation& triangulation);
  void draw_dirichlet(Domain& domain);
  void draw_triangles(std::vector<Anchor> anchors);
  void draw_triangulation(Delaunay_triangulation& triangulation, Anchor& anchor);

private:
  // Sub-methods for drawing edges and vertices
  void draw_point(QPainter* painter, Point position);

  void draw_edge(QPainter* painter, Point source, Point target);
  void draw_line(QPainter* painter, double point_1_x, double point_1_y, double point_2_x, double point_2_y);
  void draw_arc(QPainter* painter, double point_1_x, double point_1_y, double point_2_x, double point_2_y, double center_x, double center_y);

  double deg_angle(double x, double y);
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DemoWindow
  : public CGAL::Qt::DemosMainWindow,
    public Ui::MainWindow
{
  Q_OBJECT // Qt macro for Qt objects
  // (Q_OBJECT does not support templates)

private:
  QGraphicsScene scene_;
  DemoWindowItem* item_;
  QString         path_;

public:
  DemoWindow();
  DemoWindowItem& item();

  // Events handling
  void keyPressEvent(QKeyEvent* event);

public Q_SLOTS:
  void on_actionSave_as_SVG_triggered();
  void on_actionSave_as_PNG_triggered();
};

#endif // CGAL_TRIANGULATION_ON_HYPERBOLIC_SURFACE_DEMO_WINDOW
