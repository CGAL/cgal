//#define USE_SURFACE_MESH

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#ifdef USE_SURFACE_MESH
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#else
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#endif
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>

#include <CGAL/Polygon_mesh_slicer_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Timer.h>

#include <functional>
#include <boost/foreach.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef K::Iso_cuboid_3 Iso_cuboid_3;

#ifdef USE_SURFACE_MESH
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
#else
typedef CGAL::Polyhedron_3<K> Mesh;
#endif

typedef CGAL::AABB_halfedge_graph_segment_primitive<Mesh> HGSP;
typedef CGAL::AABB_traits<K, HGSP>    AABB_traits;
typedef CGAL::AABB_tree<AABB_traits>  AABB_tree;
typedef std::vector<K::Point_3> Polyline;
typedef std::list< Polyline > Polylines;

typedef CGAL::Timer Timer;

template <typename PolygonMesh>
 class Point_projector
  : public std::unary_function<typename boost::graph_traits<PolygonMesh>::vertex_descriptor,
                               typename boost::property_traits<typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type>::value_type> {
  typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type Ppmap;
public:
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor value_type;
  typedef typename Ppmap::value_type result_type;
private:
  Ppmap* ppmap;

public:
  Point_projector()
    : ppmap(NULL)
  {}

  Point_projector(PolygonMesh& pm)
    : ppmap(&get(CGAL::vertex_point, pm))
  {}
  
  result_type operator()(const value_type& v) const 
  {
    assert(ppmap != NULL);
    return (*ppmap)[v];
  }
};

template <typename PolygonMesh>
class Point_iterator :
public boost::transform_iterator<Point_projector<PolygonMesh>, typename boost::graph_traits<PolygonMesh>::vertex_iterator>
{

  typedef typename boost::transform_iterator<Point_projector<PolygonMesh>, typename boost::graph_traits<PolygonMesh>::vertex_iterator> Base;

public:
  Point_iterator()
  {}

  Point_iterator(typename boost::graph_traits<PolygonMesh>::vertex_iterator it, const PolygonMesh& pm)
    : Base(it,Point_projector<PolygonMesh>(const_cast<PolygonMesh&>(pm)))
  {}
};


template <typename PolygonMesh>
CGAL::Iterator_range<Point_iterator<PolygonMesh> > points(const PolygonMesh& m)
{ 
  return CGAL::make_range(Point_iterator<PolygonMesh>(vertices(m).begin(),m), Point_iterator<PolygonMesh>(vertices(m).end(),m));
}


int main(int argc, char* argv[])
{
  std::ifstream input(argv[1]);
  Mesh m;
  int N = 10;
  if (!input || !(input >> m)){
    std::cerr << "Error: can not read file.\n";
    return 1;
  } 

  Timer t;
  t.start();
  std::cerr << "bbox"<< std::endl;
  Iso_cuboid_3 ic = CGAL::bounding_box(points(m).begin(), points(m).end());
  Point_3 p = midpoint(ic.min(), ic.max());
  double zmin = ic.min().z();
  double zmax = ic.max().z();
  double delta = (zmax - zmin)/N;
 
  std::cerr << "slicer"<< std::endl;
  CGAL::Polygon_mesh_slicer_3<Mesh, K> slicer(m);


  int polycount = 0;
  int vertex_count = 0;
  for(int i=0; i < N; i++){
    Polylines polylines;
    slicer(K::Plane_3(Point_3(0,0,zmin+delta*i), Vector_3(0,0,1)), std::back_inserter(polylines));
    t.stop();
    polycount += polylines.size();
    BOOST_FOREACH(Polyline pl, polylines){
      vertex_count += pl.size();
      std::cout << pl.size();
      BOOST_FOREACH(Point_3 p, pl){
        std::cout << " " << p;
      }
      std::cout << std::endl;
    }
    t.start();
    
  }
  t.stop();
  std::cerr << N << " layers in a model with " << num_faces(m) << " triangles"<< std::endl; 
  std::cerr << polycount << " polylines with in total " << vertex_count << " vertices computed in "<< t.time() << " sec." << std::endl;
  
  return 0;
}
