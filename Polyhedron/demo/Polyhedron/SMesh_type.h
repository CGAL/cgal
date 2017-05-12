#ifndef SMESH_TYPE_H
#define SMESH_TYPE_H

#include <boost/scoped_ptr.hpp>
#include <boost/array.hpp>
#include <set>
#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>
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



namespace CGAL{
SMesh::Property_map<boost::graph_traits<SMesh>::halfedge_descriptor,bool>
inline get(CGAL::halfedge_is_feature_t, SMesh& smesh)
{
  typedef boost::graph_traits<SMesh >::halfedge_descriptor halfedge_descriptor;
  return smesh.add_property_map<halfedge_descriptor,bool>("h:is_feature").first;
}




template < typename I>
SMesh::Property_map< boost::graph_traits<SMesh >::face_descriptor,I>
inline get(CGAL::face_patch_id_t<I>, SMesh& smesh)
{
  typedef boost::graph_traits<SMesh >::face_descriptor face_descriptor;
  return smesh.add_property_map<face_descriptor,I>("f:patch_id").first;
}


 CGAL::static_property_map<boost::graph_traits<SMesh> ::face_descriptor,std::pair<int,int> >
   inline get(CGAL::face_patch_id_t<void>, SMesh & smesh)
 {
   typedef CGAL::static_property_map<boost::graph_traits<SMesh>::face_descriptor,
       std::pair<int,int> > Pmap;

   return Pmap(std::make_pair(0,1));
 }

 SMesh::Property_map< boost::graph_traits<SMesh >::face_descriptor,int>
 inline get(CGAL::face_selection_t, SMesh& smesh)
{
  typedef  boost::graph_traits<SMesh >::face_descriptor face_descriptor;
  return smesh.add_property_map<face_descriptor,int>("f:selection").first;
}



 SMesh::Property_map< boost::graph_traits<SMesh >::vertex_descriptor,int>
inline get(CGAL::vertex_selection_t, SMesh& smesh)
{
  typedef  boost::graph_traits<SMesh >::vertex_descriptor vertex_descriptor;
  return smesh.add_property_map<vertex_descriptor,int>("v:selection").first;
}

 SMesh::Property_map< boost::graph_traits<SMesh >::vertex_descriptor,int>
 inline get(CGAL::vertex_num_feature_edges_t, SMesh& smesh)
{
  typedef  boost::graph_traits<SMesh >::vertex_descriptor vertex_descriptor;
  return smesh.add_property_map<vertex_descriptor,int>("v:nfe").first;
}


 template < typename Patch_ID>
 SMesh::Property_map< boost::graph_traits<SMesh >::vertex_descriptor,std::set<int> >
 inline get(CGAL::vertex_incident_patches_t<Patch_ID>, SMesh& smesh)
{
  typedef  boost::graph_traits<SMesh >::vertex_descriptor vertex_descriptor;
  return smesh.add_property_map<vertex_descriptor,std::set<int> >("v:ip").first;
}


 SMesh::Property_map< boost::graph_traits<SMesh >::vertex_descriptor,std::size_t>
 inline get(CGAL::vertex_time_stamp_t, SMesh & smesh)
 {
   typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;
   return smesh.add_property_map<vertex_descriptor,std::size_t>("v:time_stamp").first;
 }



 SMesh::Property_map< boost::graph_traits<SMesh>::halfedge_descriptor,std::size_t>
 inline get(CGAL::halfedge_time_stamp_t, SMesh& smesh)
 {
   typedef boost::graph_traits<SMesh>::halfedge_descriptor halfedge_descriptor;
   return smesh.add_property_map<halfedge_descriptor,std::size_t>("h:time_stamp").first;
 }



 SMesh::Property_map< boost::graph_traits<SMesh>::face_descriptor,std::size_t>
 inline get(CGAL::face_time_stamp_t, SMesh& smesh)
 {
   typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
   return smesh.add_property_map<face_descriptor,std::size_t>("v:time_stamp").first;
 }

}
namespace boost
{
template<>
struct property_map<SMesh, CGAL::halfedge_is_feature_t>
{
  typedef boost::graph_traits<SMesh>::halfedge_descriptor halfedge_descriptor;

  typedef SMesh::Property_map<halfedge_descriptor, bool> type;
  typedef type const_type;
};

template<typename I>
struct property_map<SMesh, CGAL::face_patch_id_t<I> >
{
  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;

  typedef SMesh::Property_map<face_descriptor, I> type;
  typedef type const_type;
};

template <>
 struct property_map<SMesh, CGAL::face_patch_id_t<void> >
 {

   typedef typename boost::graph_traits<SMesh>::face_descriptor face_descriptor;

   typedef CGAL::static_property_map< boost::graph_traits<SMesh>::face_descriptor,std::pair<int,int> > type;
   typedef type const_type;
 };


template<>
struct property_map<SMesh, CGAL::vertex_selection_t>
{

  typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;

  typedef SMesh::Property_map<vertex_descriptor, int> type;
  typedef type const_type;
};

template<>
struct property_map<SMesh, CGAL::face_selection_t>
{

  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;

  typedef SMesh::Property_map<face_descriptor, int> type;
  typedef type const_type;
};

template<>
struct property_map<SMesh, CGAL::vertex_num_feature_edges_t>
{

  typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;

  typedef SMesh::Property_map<vertex_descriptor, int> type;
  typedef type const_type;
};

template<typename Patch_ID>
struct property_map<SMesh, CGAL::vertex_incident_patches_t<Patch_ID>>
{

  typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;

  typedef SMesh::Property_map<vertex_descriptor, std::set<int>> type;
  typedef type const_type;
};



template<>
struct property_map<SMesh, CGAL::vertex_time_stamp_t>
{

  typedef  boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;

  typedef  SMesh::Property_map<vertex_descriptor, std::size_t> type;
  typedef type const_type;
};

template<>
struct property_map<SMesh, CGAL::halfedge_time_stamp_t>
{

  typedef  boost::graph_traits<SMesh >::halfedge_descriptor halfedge_descriptor;

  typedef  SMesh::Property_map<halfedge_descriptor, std::size_t> type;
  typedef type const_type;
};


template<>
struct property_map<SMesh, CGAL::face_time_stamp_t>
{

  typedef boost::graph_traits<SMesh >::face_descriptor face_descriptor;

  typedef SMesh::Property_map<face_descriptor, std::size_t> type;
  typedef type const_type;
};
} //boost


#endif // SMESH_TYPE_H
