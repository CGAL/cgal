// #define USE_SURFACE_MESH

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#ifdef USE_SURFACE_MESH
#include <CGAL/Surface_mesh.h>
#else
#include <CGAL/Polyhedron_3.h>
#endif
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>

#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polygon_2.h>
#include <boost/foreach.hpp>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Epec;

template <typename K>
int test_slicer()
{
#ifdef USE_SURFACE_MESH
  typedef CGAL::Surface_mesh<K::Point_3> Mesh;
#else
  typedef CGAL::Polyhedron_3<K> Mesh;
#endif
  
  typedef CGAL::AABB_halfedge_graph_segment_primitive<Mesh> HGSP;
  typedef CGAL::AABB_traits<K, HGSP>    AABB_traits;
  typedef CGAL::AABB_tree<AABB_traits>  AABB_tree;
  typedef std::vector<typename K::Point_3> Polyline_type;
  typedef std::list< Polyline_type > Polylines;
  
  //API test
  {
    std::ifstream input("data/U.off");
    Mesh m;
    
    if (!input || !(input >> m)){
      std::cerr << "Error: can not read file.\n";
      return 1;
    }
    
    AABB_tree tree(edges(m).first, edges(m).second, m);
    
    CGAL::Polygon_mesh_slicer<Mesh, K> slicer(m, tree);
    Polylines polylines;
    slicer(typename K::Plane_3(1,1,0,0), std::back_inserter(polylines));
    assert(polylines.size()==1);
  }
  
  std::ifstream input("data_slicer/open_cube_meshed.off");
  Mesh m;
  
  if (!input || !(input >> m)){
    std::cerr << "Error: can not read file.\n";
    return 1;
  }
  
  CGAL::Polygon_mesh_slicer<Mesh, K> slicer(m);
  
  Polylines polylines;
  
  // test isolated vertex
  typename K::Plane_3 plane(0,1,0,0);
  slicer(plane, std::back_inserter(polylines));
  assert(polylines.size()==2); // two polylines
  assert( (polylines.front().size()==1) != (polylines.back().size()==1)); //only one isolated vertex
  
  CGAL::Polygon_2<K> polygon;
  BOOST_FOREACH(const typename K::Point_3& p, polylines.back())
  {
    polygon.push_back(plane.to_2d(p));
  }
  
  assert(polygon.is_counterclockwise_oriented());
  
  //test two nested polylines, one open and one closed
  polylines.clear();
  polygon.clear();
  plane = typename K::Plane_3(0,1,0,0.5); 
  slicer(plane, std::back_inserter(polylines));
  assert(polylines.size()==2);// two polylines
  assert( (polylines.front().front()==polylines.front().back()) !=
      (polylines.back().front()==polylines.back().back()) ); //one open and one closed polyline
  
  
  BOOST_FOREACH(const Polyline_type& polyline, polylines)
  {
    BOOST_FOREACH(const typename K::Point_3& p, polyline)
    {
      polygon.push_back(plane.to_2d(p));
    }
    if(polyline.front() == polyline.back())
      polygon.erase(polygon.vertices_end()-1);
    assert(polygon.is_counterclockwise_oriented());
    polygon.clear();
  }
  polylines.clear();
  
  // test only coplanar edges
  plane = typename K::Plane_3(0,0,1,1);
  slicer(plane, std::back_inserter(polylines));
  assert(polylines.size()==1); // one polyline
  assert(polylines.front().front()==polylines.front().back()); // that is closed
  
  
  BOOST_FOREACH(const Polyline_type& polyline, polylines)
  {
    BOOST_FOREACH(const typename K::Point_3& p, polyline)
    {
      polygon.push_back(plane.to_2d(p));
    }
    if(polyline.front() == polyline.back())
      polygon.erase(polygon.vertices_end()-1);
    assert(polygon.is_counterclockwise_oriented());
    polygon.clear();
  }
  polylines.clear();
  
  //test only coplanar border edges
  plane = typename K::Plane_3(0,0,1,-1);
  slicer(plane, std::back_inserter(polylines));
  assert(polylines.size()==1); // one polyline
  assert(polylines.front().front()!=polylines.front().back()); // that is closed
  
  BOOST_FOREACH(const Polyline_type& polyline, polylines)
  {
    BOOST_FOREACH(const typename K::Point_3& p, polyline)
    {
      polygon.push_back(plane.to_2d(p));
    }
    assert(polygon.is_counterclockwise_oriented());
    polygon.clear();
  }
  polylines.clear();
  //test no intersection
  slicer(typename K::Plane_3(0,0,1,333), std::back_inserter(polylines));
  assert(polylines.empty());
  
  return 0;
}


int main()
{
  assert(test_slicer<Epic>() == 0);
  assert(test_slicer<Epec>() == 0);
  
  return 0;
}
