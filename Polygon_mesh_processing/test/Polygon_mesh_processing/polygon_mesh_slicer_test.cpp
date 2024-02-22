// #define USE_SURFACE_MESH

#include <CGAL/Polygon_mesh_slicer.h>

#ifdef USE_SURFACE_MESH
#include <CGAL/Surface_mesh.h>
#else
#include <CGAL/Polyhedron_3.h>
#endif

#include <CGAL/AABB_halfedge_graph_segment_primitive.h>

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <fstream>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Epec;

template<class K, class Polyline_type>
bool is_ccw(int xi, int yi,
    const Polyline_type& polyline)
{
  CGAL::Polygon_2<K> polygon;
  if(polyline.front() == polyline.back())
  {
    for(const typename K::Point_3& p : polyline)
    {
      polygon.push_back(typename K::Point_2(p[xi], p[yi]));
    }
    polygon.erase(polygon.vertices_end()-1);
    return polygon.is_counterclockwise_oriented();
  }
  return true;
}
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
  typedef std::vector< Polyline_type > Polylines;

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
  slicer(typename K::Plane_3(0,1,0,0), std::back_inserter(polylines));
  assert(polylines.size()==2); // two polylines
  assert( (polylines.front().size()==1) != (polylines.back().size()==1)); //only one isolated vertex


  //test two nested polylines, one open and one closed
  polylines.clear();
  slicer(typename K::Plane_3(0,1,0,0.5), std::back_inserter(polylines));
  assert(polylines.size()==2);// two polylines
  assert( (polylines.front().front()==polylines.front().back()) !=
          (polylines.back().front()==polylines.back().back()) ); //one open and one closed polyline

  // test only coplanar edges
  polylines.clear();
  slicer(typename K::Plane_3(0,0,1,1), std::back_inserter(polylines));
  assert(polylines.size()==1); // one polyline
  assert(polylines.front().front()==polylines.front().back()); // that is closed


  //test only coplanar border edges
  polylines.clear();
  slicer(typename K::Plane_3(0,0,1,-1), std::back_inserter(polylines));
  assert(polylines.size()==1); // one polyline
  assert(polylines.front().front()!=polylines.front().back()); // that is closed

  //test no intersection
  polylines.clear();
  slicer(typename K::Plane_3(0,0,1,333), std::back_inserter(polylines));
  assert(polylines.empty());

// Now test the orientation of polylines
  polylines.clear();
  slicer(typename K::Plane_3(0,1,0,0.5), std::back_inserter(polylines));
  assert(polylines.size()==2); // two polylines
  int closed_id = polylines.front().front()==polylines.front().back() ? 0 : 1;

  assert( is_ccw<K>(0, 2 , polylines[closed_id]) );

  polylines.clear();
  slicer(typename K::Plane_3(0,-1,0,-0.5), std::back_inserter(polylines));
  assert(polylines.size()==2); // two polylines
  closed_id = polylines.front().front()==polylines.front().back() ? 0 : 1;
  assert( !is_ccw<K>(0, 2, polylines[closed_id]) );

  polylines.clear();
  slicer(typename K::Plane_3(0,0,1,1), std::back_inserter(polylines));
  assert(polylines.size()==1); // one polyline
  assert( is_ccw<K>(0, 1 , polylines[0]) );

  polylines.clear();
  slicer(typename K::Plane_3(0,0,-1,-1), std::back_inserter(polylines));
  assert(polylines.size()==1); // one polyline
  assert( !is_ccw<K>(0, 1 , polylines[0]) );

  // reverse face orientation (no need to rebuild the tree)
  CGAL::Polygon_mesh_processing::reverse_face_orientations(m);
  polylines.clear();
  slicer(typename K::Plane_3(0,1,0,0.5), std::back_inserter(polylines));
  assert(polylines.size()==2); // two polylines
  closed_id = polylines.front().front()==polylines.front().back() ? 0 : 1;
  assert( !is_ccw<K>(0, 2 , polylines[closed_id]) );

  polylines.clear();
  slicer(typename K::Plane_3(0,-1,0,-0.5), std::back_inserter(polylines));
  assert(polylines.size()==2); // two polylines
  closed_id = polylines.front().front()==polylines.front().back() ? 0 : 1;
  assert( is_ccw<K>(0, 2, polylines[closed_id]) );

  polylines.clear();
  slicer(typename K::Plane_3(0,0,1,1), std::back_inserter(polylines));
  assert(polylines.size()==1); // one polyline
  assert( !is_ccw<K>(0, 1 , polylines[0]) );

  polylines.clear();
  slicer(typename K::Plane_3(0,0,-1,-1), std::back_inserter(polylines));
  assert(polylines.size()==1); // one polyline
  assert( is_ccw<K>(0, 1 , polylines[0]) );



  return 0;
}


int main()
{
  assert(test_slicer<Epic>() == 0);
  assert(test_slicer<Epec>() == 0);

  return 0;
}
