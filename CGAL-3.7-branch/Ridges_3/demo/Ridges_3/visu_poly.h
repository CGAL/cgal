#ifndef CGAL_VISU_POLY_H_
#define CGAL_VISU_POLY_H_

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>
//marc
#include "enriched_polyhedron.h"

typedef  CGAL::Cartesian<double>               Kernel;
typedef  Kernel::Point_3                       Point;
//marc
/* typedef  CGAL::Polyhedron_3<Kernel>            Mesh; */
typedef Enriched_polyhedron<Kernel,Enriched_items> Mesh;

typedef Kernel::FT	FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef Mesh::Vertex_handle Vertex_handle;
typedef Mesh::Facet_handle Facet_handle;
typedef Mesh::Vertex Vertex;
typedef Mesh::Facet Facet;
typedef Mesh::Face_handle Face_handle;
typedef Mesh::Halfedge_handle Halfedge_handle;
typedef Mesh::Facet_iterator Facet_iterator;
typedef Mesh::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
typedef Mesh::Halfedge_iterator Halfedge_iterator;
typedef Mesh::Point_iterator Point_iterator;
typedef Mesh::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
typedef Mesh::Vertex_iterator Vertex_iterator;
typedef Mesh::Vertex_const_iterator   Vertex_const_iterator;
typedef Mesh::Edge_iterator Edge_iterator ;
typedef Mesh::Halfedge_around_facet_const_circulator Halfedge_around_facet_const_circulator;

typedef CGAL::Inverse_index < Vertex_const_iterator > Vertex_index;

#endif
