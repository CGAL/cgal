#ifndef SMESH_TYPE_H
#define SMESH_TYPE_H

#include "properties.h"

#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


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


namespace CGAL {

template <typename P, typename Property_tag>
struct Get_pmap_of_surface_mesh_ {
  typedef typename boost::property_map<Surface_mesh<P>, Property_tag >::type type;
};

#define CGAL_PROPERTY_SURFACE_MESH_RETURN_TYPE(Tag) \
  typename boost::lazy_disable_if<                      \
     boost::is_const<P>,                                \
     Get_pmap_of_surface_mesh_<P, Tag >                  \
   >::type

template <typename P>
CGAL_PROPERTY_SURFACE_MESH_RETURN_TYPE(CGAL::face_selection_t)
 inline get(CGAL::face_selection_t, Surface_mesh<P> & smesh)
{
 typedef typename boost::graph_traits<Surface_mesh<P> >::face_descriptor face_descriptor;
  return smesh. template add_property_map<face_descriptor,int>("f:selection").first;
}



template <typename P>
CGAL_PROPERTY_SURFACE_MESH_RETURN_TYPE(CGAL::vertex_selection_t)
inline get(CGAL::vertex_selection_t, Surface_mesh<P> & smesh)
{
  typedef typename boost::graph_traits<Surface_mesh<P> >::vertex_descriptor vertex_descriptor;
  return smesh. template add_property_map<vertex_descriptor,int>("v:selection").first;
}
#undef CGAL_PROPERTY_SURFACE_MESH_RETURN_TYPE
} // namespace CGAL

#endif // SMESH_TYPE_H
