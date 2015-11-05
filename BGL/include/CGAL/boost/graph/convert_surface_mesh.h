

#ifndef CGAL_BOOST_GRAPH_CONVERT_SURFACE_MESH_H
#define CGAL_BOOST_GRAPH_CONVERT_SURFACE_MESH_H

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/iterator.h>


namespace CGAL {

  template <typename SourceMesh, typename TargetMesh, typename V2V, typename H2H>
  void convert_surface_mesh(const SourceMesh& sm, TargetMesh& tm, V2V& v2v, H2H& h2h)
{
  typedef boost::graph_traits<SourceMesh>::vertex_descriptor sm_vertex_descriptor;
  typedef boost::graph_traits<TargetMesh>::vertex_descriptor tm_vertex_descriptor;

  typedef boost::graph_traits<SourceMesh>::face_descriptor sm_face_descriptor;
  typedef boost::graph_traits<TargetMesh>::face_descriptor tm_face_descriptor;

  typedef boost::graph_traits<SourceMesh>::halfedge_descriptor sm_halfedge_descriptor;
  typedef boost::graph_traits<TargetMesh>::halfedge_descriptor tm_halfedge_descriptor;

  typedef typename boost::property_map<SourceMesh, vertex_point_t>::const_type sm_PMap;
  typedef typename boost::property_map<TargetMesh, vertex_point_t>::type tm_PMap;

  sm_PMap sm_pmap = get(vertex_point, sm);
  tm_PMap tm_pmap = get(vertex_point, tm);


  BOOST_FOREACH(sm_vertex_descriptor svd, vertices(sm)){
    tm_vertex_descriptor tvd = add_vertex(tm);
    v2v.insert(std::make_pair(svd, tvd));
    put(tm_pmap, tvd, get(sm_pmap, svd));
  }

  boost::unordered_map<sm_face_descriptor, tm_face_descriptor> f2f;
  BOOST_FOREACH(sm_face_descriptor sfd, faces(sm)){
    std::vector<tm_vertex_descriptor> tv;
    BOOST_FOREACH(sm_vertex_descriptor svd, vertices_around_face(halfedge(sfd,sm),sm)){
      tv.push_back(v2v.at(svd));
    }
    f2f[sfd] = Euler::add_face(tv,tm);
  }
  
  BOOST_FOREACH(sm_face_descriptor sfd, faces(sm)){
    sm_halfedge_descriptor shd = halfedge(sfd,sm), done(shd);
    tm_halfedge_descriptor thd = halfedge(f2f[sfd],tm);
    tm_vertex_descriptor tvd = v2v.at(target(shd,sm));
    while(target(thd,tm) != tvd){
      thd = next(thd,tm);
    }
    do {
      h2h.insert(std::make_pair(shd, thd));
      shd = next(shd,sm);
      thd = next(thd,tm);
    }while(shd != done);
  }
  
}

} // namespace CGAL

#endif //  CGAL_BOOST_GRAPH_CONVERT_SURFACE_MESH_H
