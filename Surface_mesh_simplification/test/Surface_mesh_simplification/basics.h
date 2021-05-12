#ifndef BASICS_H
#define BASICS_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <exception>
#include <map>

#include <boost/format.hpp>

#define CGAL_CHECK_EXPENSIVE

//#define CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE   4
//#define CGAL_SURFACE_SIMPLIFICATION_ENABLE_LT_TRACE 4

void Surface_simplification_external_trace(std::string s)
{
  static std::ofstream out("log.txt");
  out << s << std::endl;
}

#include <CGAL/Real_timer.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>

#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <CGAL/assertions_behaviour.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT FT;

typedef Kernel::Vector_3     Vector;
typedef Kernel::Point_3      Point;

typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3> Surface;

typedef Surface::Vertex                                  Vertex;
typedef Surface::Vertex_iterator                         Vertex_iterator;
typedef Surface::Vertex_handle                           Vertex_handle;
typedef Surface::Vertex_const_handle                     Vertex_const_handle;
typedef Surface::Halfedge_handle                         Halfedge_handle;
typedef Surface::Halfedge_const_handle                   Halfedge_const_handle;
typedef Surface::Edge_iterator                           Edge_iterator;
typedef Surface::Facet_iterator                          Facet_iterator;
typedef Surface::Facet_const_iterator                    Facet_const_iterator;
typedef Surface::Facet_const_handle                      Facet_const_handle;
typedef Surface::Halfedge_around_vertex_const_circulator HV_circulator;
typedef Surface::Halfedge_around_facet_circulator        HF_circulator;
typedef Surface::size_type                               size_type;

using namespace std;
using namespace boost;
using namespace CGAL;



#endif // BASICS_H
