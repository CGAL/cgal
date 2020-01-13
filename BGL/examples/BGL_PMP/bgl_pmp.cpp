
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <pmp/SurfaceMesh.h>

#include <CGAL/boost/graph/graph_traits_SurfaceMesh.h>
#include <CGAL/boost/graph/properties_SurfaceMesh.h>

#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <CGAL/boost/graph/helpers.h>
#include <iostream>
#include <fstream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;

typedef boost::graph_traits<pmp::SurfaceMesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<pmp::SurfaceMesh>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<pmp::SurfaceMesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<pmp::SurfaceMesh>::face_descriptor face_descriptor;

int main()
{
  typedef boost::property_map<pmp::SurfaceMesh, boost::vertex_index_t>::type VIM;
  typedef boost::property_map<pmp::SurfaceMesh, boost::edge_index_t>::type EIM;
  typedef boost::property_map<pmp::SurfaceMesh, CGAL::vertex_point_t>::type VPM;
  pmp::SurfaceMesh sm;
  pmp::Vertex v0, v1, v2, v3;
  v0 = sm.add_vertex(pmp::Point(0,0,0));
  v1 = sm.add_vertex(pmp::Point(1,0,0));
  v2 = sm.add_vertex(pmp::Point(0,2,0));

  sm.add_triangle(v0,v1,v2);

  std::list<edge_descriptor> mst;

  boost::kruskal_minimum_spanning_tree(sm, 
                                       std::back_inserter(mst));

  edge_descriptor e = *(edges(sm).begin());
  std::cout << e << std::endl;
  halfedge_descriptor h = halfedge(e,sm);
  std::cout << h << std::endl;
  edge_descriptor e2 = edge(h,sm);
  std::cout << e2 << std::endl;
  halfedge_descriptor h2 = halfedge(e2,sm);
  std::cout << h2 << std::endl;

  h = opposite(h,sm);
  e = edge(h,sm);
  std::cout << e << std::endl;
  h = halfedge(e,sm);
  std::cout << h << std::endl;
  
  
  VIM vim = get(boost::vertex_index, sm);
  VPM vpm = get(CGAL::vertex_point, sm);
  for(auto v : vertices(sm)){
    std::cout << degree(v,sm) << " " << get(vim,v) << " " << get(vpm,v) << std::endl;
  }

  boost::property_traits<VPM>::value_type vt;

  
  EIM eim = get(boost::edge_index,sm);
  
  for(auto e : edges(sm)){
    std::cout << e  << " " <<  get(eim,e) << std::endl;
  }

  typedef boost::property_map<pmp::SurfaceMesh, CGAL::dynamic_vertex_property_t<int>>::type V_index_map;
  V_index_map dvim;
  dvim = get(CGAL::dynamic_vertex_property_t<int>(), sm);

#if 1
  // Incrementally fill the holes
  unsigned int nb_holes = 0;
  for(halfedge_descriptor h : halfedges(sm))
  {
    if(CGAL::is_border(h,sm))
    {
      std::vector<face_descriptor>  patch_facets;
      std::vector<vertex_descriptor> patch_vertices;
      bool success = std::get<0>(
        CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
                  sm,
                  h,
                  std::back_inserter(patch_facets),
                  std::back_inserter(patch_vertices),
     CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, sm)).
                  geom_traits(Kernel())) );

      CGAL_assertion(CGAL::is_valid_polygon_mesh(sm));

      std::cout << "* FILL HOLE NUMBER " << ++nb_holes << std::endl;
      std::cout << "  Number of facets in constructed patch: " << patch_facets.size() << std::endl;
      std::cout << "  Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
      std::cout << "  Is fairing successful: " << success << std::endl;
    }
  }
#endif
  
  return 0;
}
