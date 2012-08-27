// Copyright (c) 2003-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Rineau

#ifndef MESH_2_SHOW_CLUSTERS_H
#define MESH_2_SHOW_CLUSTERS_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>

#include <CGAL/Delaunay_triangulation_2.h>

#include <qpixmap.h>

class Show_clusters_aux : public CGAL::Qt_widget_layer
{
  Q_OBJECT
private:
  virtual void reinit_clusters() {}

public slots:
  void reinitClusters()
  {
    reinit_clusters();
  }
public:
  Show_clusters_aux(QObject* parent, const char* name)
    : Qt_widget_layer(parent, name)
  {
  }
};

template <class Mesher>
class Show_clusters : public Show_clusters_aux
{
public:
  typedef typename Mesher::Triangulation Tr;
  typedef typename Tr::Point		Point;
  typedef typename Tr::Geom_traits         Geom_traits;
  typedef CGAL::Delaunay_triangulation_2<Geom_traits> DT;
  typedef typename DT::Finite_vertices_iterator Vertices_iterator;
  typedef typename DT::Vertex_handle DT_vertex_handle;
  typedef typename Tr::Segment		Segment;
  typedef typename Tr::Face_handle		Face_handle;
  typedef typename Tr::Vertex_handle	Vertex_handle;
  typedef typename Tr::Geom_traits::FT	FT;
  typedef typename Mesher::Clusters Clusters;
  typedef typename Clusters::Cluster Cluster;
  typedef typename Clusters::Cluster_vertices_iterator CVIt;
  typedef typename Clusters::Vertices_in_cluster_iterator ViCIt;
  typedef typename Clusters::const_iterator Clusters_const_iterator;
  typedef std::list<Point> List_of_points;
  typedef typename List_of_points::const_iterator Point_iterator;

  Show_clusters(Mesher* m,
		CGAL::Color color_ = CGAL::GREEN,
		int pointsize = 3,
		CGAL::PointStyle pointstyle = CGAL::DISC,
		CGAL::Color lc = CGAL::RED,
                CGAL::Color reduced_line_color_ = CGAL::BLUE,
		int linewidth = 2,
                QObject* parent = 0, const char* name = 0)
    : Show_clusters_aux(parent, name),
      mesher(m), dt(), color(color_),
      size(pointsize), style(pointstyle), line_color(lc),
      reduced_line_color(reduced_line_color_),
      width(linewidth)
    {
      reinit_clusters();
    }

  void activating()
  {
    reinit_clusters();
  }

  void change_mesher(Mesher* m)
  {
    mesher = m;
    reinit_clusters();
  }

  void reinit_clusters()
  {
    if(!is_active()) return;

    dt.clear();

    if( mesher != 0 )
      for(CVIt it = mesher->clusters().clusters_vertices_begin();
          it != mesher->clusters().clusters_vertices_end();
          ++it)
        dt.push_back( (*it)->point() );
  }

  void draw()
  {
    widget->lock();

    QColor oldColor = widget->color();
    int oldPointSize = widget->pointSize();
    CGAL::PointStyle oldStyle = widget->pointStyle();

    *widget << color << CGAL::PointStyle(style)
    	    << CGAL::PointSize(size);

    for(Vertices_iterator it = dt.finite_vertices_begin();
        it != dt.finite_vertices_end();
        ++it)
      *widget << it->point();

    widget->setPointStyle(oldStyle);
    widget->setPointSize(oldPointSize);
    widget->setColor(oldColor);

    widget->unlock();
    oldPixmap = widget->get_pixmap();
  }

  void mouseMoveEvent(QMouseEvent *e)
  {
    if( mesher == 0 ) return;

    FT x, y;
    widget->x_real(e->x(), x);
    widget->y_real(e->y(), y);
    Point p(x, y);

    DT_vertex_handle v = dt.nearest_vertex(p);

    if(v == NULL) return;
    if(v == oldVertex) return;

    oldVertex = v;

    QColor oldColor = widget->color();
    int oldWidth = widget->lineWidth();

    widget->lock();

    widget->get_painter().drawPixmap(0, 0, oldPixmap);

    *widget << CGAL::LineWidth(width);

    typename Tr::Locate_type lt;
    int i;
    Face_handle fh = mesher->triangulation().locate(v->point(), lt, i);
    CGAL_assertion( lt == Tr::VERTEX );

    Vertex_handle v2 = fh->vertex(i);

    int n = mesher->clusters().number_of_clusters_at_vertex(v2);

    for(int j = 0; j < n; ++j)
      {
	std::pair<ViCIt,ViCIt> seq =
          mesher->clusters().vertices_in_cluster_sequence(v2, j);
        Cluster c;
        Clusters_const_iterator dummy_c_it;
        mesher->clusters().get_cluster(v2, *(seq.first), c, dummy_c_it);
        if( c.is_reduced() )
          *widget << reduced_line_color;
        else
          *widget << line_color;
	for(ViCIt it = seq.first;
	    it != seq.second;
	    ++it)
	  *widget << Segment(v2->point(), (*it)->point());
      }

    widget->setLineWidth(oldWidth);
    widget->setColor(oldColor);

    widget->unlock();
  }

  void leaveEvent(QEvent *)
  {
    widget->get_painter().drawPixmap(0, 0, oldPixmap);
    widget->update();
  }

private:
  Mesher* mesher;
  DT dt;
  DT_vertex_handle oldVertex;
  QPixmap oldPixmap;
  bool  should_restore_pixmap;
  CGAL::Color color;
  int size;
  CGAL::PointStyle style;
  CGAL::Color line_color;
  CGAL::Color reduced_line_color;
  int width;
};

#endif
