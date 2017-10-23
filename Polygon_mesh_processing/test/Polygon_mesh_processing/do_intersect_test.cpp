
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Epic;
typedef CGAL::Exact_predicates_exact_constructions_kernel       Epec;

template<class Point>
int load_polyline(std::ifstream& input,
                   std::vector<Point>& points)
{
  std::size_t n;
  input >> n;
  points.reserve(n);
  while(n--){
    Point p;
    input >> p;
    points.push_back(p);
    if(!input.good())
    {
      std::cerr << "Error: cannot read file: " << std::endl;
      return 1;
    }
  }
  return 0;
}

template <typename K>
int
test_faces_intersections(const char* filename1,
                         const char* filename2,
                         const bool expected)
{
  typedef CGAL::Surface_mesh<typename K::Point_3>                Mesh;
  typedef typename boost::graph_traits<Mesh>::face_descriptor    face_descriptor;

  std::ifstream input1(filename1);
  std::ifstream input2(filename2);
  Mesh m1,m2;

  if ( !input1 || !(input1 >> m1) ) {
    std::cerr << "Error: cannot read file: " << filename1 << std::endl;
    return 1;
  }
  if ( !input2 || !(input2 >> m2) ) {
    std::cerr << "Error: cannot read file: " << filename2 << std::endl;
    return 1;
  }

  std::cout << "Reading files: " << filename1 <<", "<< filename2 << std::endl;

  CGAL::Timer timer;
  timer.start();

  std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
  CGAL::internal::compute_face_face_intersection(
        m1, m2,
        std::back_inserter(intersected_tris),
        CGAL::Polygon_mesh_processing::parameters::all_default(),
        CGAL::Polygon_mesh_processing::parameters::all_default());

  bool intersecting_1 = !intersected_tris.empty();

  std::cout << "intersections test took " << timer.time() << " sec." << std::endl;
  std::cout << intersected_tris.size() << " pairs of triangles are intersecting." << std::endl;

  timer.reset();
  bool intersecting_2 = CGAL::Polygon_mesh_processing::do_intersect(m1, m2);

  std::cout << "does_intersect test took " << timer.time() << " sec." << std::endl;
  std::cout << (intersecting_2 ? "There are intersections." :
                                 "There are no intersections.") << std::endl;

  assert(intersecting_1 == intersecting_2);
  assert(intersecting_1 == expected);

  std::cout << filename1 << "and " <<filename2  << " passed the tests." << std::endl << std::endl;

  return 0;
}

template <typename K>
int
test_faces_polyline_intersections(const char* filename1,
                                  const char* filename2,
                                  const bool expected)
{
  typedef typename K::Point_3                                    Point;
  typedef typename CGAL::Surface_mesh<Point>                     Mesh;
  typedef typename boost::graph_traits<Mesh>::face_descriptor    face_descriptor;

  std::ifstream input1(filename1);
  std::ifstream input2(filename2);
  Mesh m;

  if ( !input1 || !(input1 >> m) ) {
    std::cerr << "Error: cannot read file: " << filename1 << std::endl;
    return 1;
  }

  if ( !input2 ) {
    std::cerr << "Error: cannot read file: " << filename2 << std::endl;
    return 1;
  }
  std::vector<Point> points;
  if(load_polyline(input2, points) >0)
    return 1;

  std::cout << "Reading files: " << filename1 <<", "<< filename2 << std::endl;

  CGAL::Timer timer;
  timer.start();

  std::vector<std::pair<std::size_t, std::size_t> > intersected_tris;
  CGAL::internal::compute_face_polyline_intersection(
        m, points,
        std::back_inserter(intersected_tris),
        CGAL::Polygon_mesh_processing::parameters::all_default());

  bool intersecting_1 = !intersected_tris.empty();

  std::cout << "intersections test took " << timer.time() << " sec." << std::endl;
  std::cout << intersected_tris.size() << " intersections." << std::endl;

  timer.reset();
  bool intersecting_2 = CGAL::Polygon_mesh_processing::do_intersect(m,points);

  std::cout << "does_intersect test took " << timer.time() << " sec." << std::endl;
  std::cout << (intersecting_2 ? "There are intersections." :
                                 "There are no intersections.") << std::endl;

  assert(intersecting_1 == intersecting_2);
  assert(intersecting_1 == expected);

  std::cout << filename1 << "and " <<filename2  << " passed the tests." << std::endl << std::endl;

  return 0;
}

template <typename K>
int
test_polylines_intersections(const char* filename1,
                                  const char* filename2,
                                  const bool expected)
{
  typedef typename K::Point_3                                    Point;
  typedef typename CGAL::Surface_mesh<Point>                     Mesh;
  typedef typename boost::graph_traits<Mesh>::face_descriptor    face_descriptor;

  std::ifstream input1(filename1);
  std::ifstream input2(filename2);

  if ( !input1 ) {
    std::cerr << "Error: cannot read file: " << filename1 << std::endl;
    return 1;
  }

  if ( !input2 ) {
    std::cerr << "Error: cannot read file: " << filename2 << std::endl;
    return 1;
  }

  std::vector<Point> points1;
  if(load_polyline(input1, points1)>0)
    return 1;
  std::vector<Point> points2;
  if(load_polyline(input2, points2)>0)
    return 1;


  std::cout << "Reading files: " << filename1 <<", "<< filename2 << std::endl;

  CGAL::Timer timer;
  timer.start();

  std::vector<std::pair<std::size_t, std::size_t> > intersected_polys;
  CGAL::internal::compute_polyline_polyline_intersection(
        points1, points2,
        std::back_inserter(intersected_polys),
        K());

  bool intersecting_1 = !intersected_polys.empty();

  std::cout << "intersections test took " << timer.time() << " sec." << std::endl;
  std::cout << intersected_polys.size() << " intersections." << std::endl;

  timer.reset();

  bool intersecting_2 = CGAL::Polygon_mesh_processing::do_intersect(points1, points2);

  std::cout << "does_intersect test took " << timer.time() << " sec." << std::endl;
  std::cout << (intersecting_2 ? "There are intersections." :
                                 "There are no intersections.") << std::endl;

  assert(intersecting_1 == intersecting_2);
  assert(intersecting_1 == expected);

  std::cout << filename1 << "and " <<filename2  << " passed the tests." << std::endl << std::endl;

  return 0;
}

int main()
{


  bool expected = true;
  const char* filename1 =  "data/tetra1.off";
  const char* filename2 =  "data/tetra2.off";
  const char* filename3 =  "data/triangle.polylines.txt";
  const char* filename4 =  "data/planar.polylines.txt";


  std::cout << "First test (Epic):" << std::endl;
  int r = test_faces_intersections<Epic>(filename1, filename2,  expected);

  std::cout << "Second test (Epec):" << std::endl;
  r += test_faces_intersections<Epec>(filename1, filename2, expected);
  std::cout << "Third test (Polyline):" << std::endl;
  r += test_faces_polyline_intersections<Epic>(filename1, filename3, expected);
  std::cout << "Fourth test (Polylines):" << std::endl;
  r += test_polylines_intersections<Epic>(filename3, filename4, expected);

  return r;
}
