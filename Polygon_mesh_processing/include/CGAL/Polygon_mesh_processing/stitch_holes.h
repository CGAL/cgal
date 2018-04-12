#ifndef CGAL_STITCH_HOLES_H
#define CGAL_STITCH_HOLES _H


#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>


namespace CGAL{

namespace Polygon_mesh_processing{


template <typename PolygonMesh, typename OutputIterator>
void extract_connected_components(PolygonMesh& mesh,
                                  OutputIterator out)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  std::set<halfedge_descriptor> border_halfedges;

  BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh))
  {
    if(is_border(h, mesh))
      border_halfedges.insert(h);
  }

  std::set<halfedge_descriptor> connected_component;
  BOOST_FOREACH(halfedge_descriptor h, border_halfedges)
  {
    if(connected_component.insert(h).second)
    {
      halfedge_descriptor start = h;
      do{
        h = next(h, mesh);
        connected_component.insert(h);
      } while(h != start);

      *out++=connected_component;
    }
  }
}


template <typename PolygonMesh, typename ConnectedComponents>
std::size_t count_identical_points(PolygonMesh& mesh,
                                   std::vector<ConnectedComponents> cc_list)
{
  // cc is a std::vector<std::set<halfedge_descriptor> >

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  typedef typename boost::property_map<PolygonMesh, boost::vertex_point_t>::const_type Vpm;
  Vpm vpm = get(boost::vertex_point, mesh);

  for(std::set<halfedge_descriptor> i_cc : cc_list)
  {
    // for each cc
    BOOST_FOREACH(halfedge_descriptor h, i_cc)
    {
      vertex_descriptor vs = source(h, mesh);
      vertex_descriptor vt = target(h, mesh);

      // find identicals
    }
  }
  return 0; // how many found
}

/// \ingroup PMP_repairing_grp
/// merges two vertices into one
///
/// @tparam TriangleMesh a model of `FaceListGraph`
///
/// @param mesh the input triangle mesh
/// @param v_keep the vertex to be kept
/// @param v_rm the vertex to be removed
template <typename PolygonMesh>
void merge_identical_points(PolygonMesh& mesh,
                            typename boost::graph_traits<PolygonMesh>::vertex_descriptor v_keep,
                            typename boost::graph_traits<PolygonMesh>::vertex_descriptor v_rm)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor h = halfedge(v_rm, mesh);
  halfedge_descriptor start = h;

  do{
    set_target(h, v_keep, mesh);
    h = opposite(next(h, mesh), mesh);
  } while( h != start );

  remove_vertex(v_rm, mesh);
}




}
}

#endif //CGAL_STITCH_HOLES_H
