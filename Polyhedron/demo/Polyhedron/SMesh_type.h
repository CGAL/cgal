#ifndef SMESH_TYPE_H
#define SMESH_TYPE_H

#include "properties.h"

#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/properties.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel EPICK;
typedef EPICK::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> SMesh;
typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<SMesh>::halfedge_descriptor halfedge_descriptor;


namespace boost {

template <typename P>
struct property_map<CGAL::Surface_mesh<P>, CGAL::vertex_selection_t>
{

  typedef typename boost::graph_traits<CGAL::Surface_mesh<P> >::vertex_descriptor vertex_descriptor;

  typedef typename CGAL::Surface_mesh<P>::template Property_map<vertex_descriptor, int> type;
  typedef type const_type;
};


template <typename P>
struct property_map<CGAL::Surface_mesh<P>, CGAL::face_selection_t>
{

  typedef typename boost::graph_traits<CGAL::Surface_mesh<P> >::face_descriptor face_descriptor;

  typedef typename CGAL::Surface_mesh<P>::template Property_map<face_descriptor, int> type;
  typedef type const_type;
};

} // namespace boost

#endif // SMESH_TYPE_H
