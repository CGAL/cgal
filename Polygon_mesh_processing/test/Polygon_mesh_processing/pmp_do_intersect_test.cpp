#include <CGAL/Polygon_mesh_processing/intersection.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Timer.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

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

template<class Point>
int load_polylines(std::ifstream& input,
                   std::vector<std::vector<Point> >& points)
{
  std::size_t n;
  while(input >> n) {
    std::vector<Point> new_polyline;
    points.push_back(new_polyline);
    std::vector<Point>&polyline = points.back();
    polyline.reserve(n);
    while(n--){
      Point p;
      input >> p;
      polyline.push_back(p);
      if(!input.good()) return 1;
    }
    std::string line_remainder;
    std::getline(input, line_remainder);

    if(input.bad() || input.fail()) return 1;
    }
  return 0;
}

template <typename K>
int
test_faces_intersections(const std::string filename1,
                         const std::string filename2,
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
  CGAL::Polygon_mesh_processing::internal::compute_face_face_intersection(
        m1, m2,
        std::back_inserter(intersected_tris),
        CGAL::parameters::default_values(),
        CGAL::parameters::default_values());

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
test_faces_polyline_intersections(const std::string filename1,
                                  const std::string filename2,
                                  const bool expected)
{
  typedef typename K::Point_3                                    Point;
  typedef typename CGAL::Surface_mesh<Point>                     Mesh;

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
  CGAL::Polygon_mesh_processing::internal::compute_face_polyline_intersection(
        m, points,
        std::back_inserter(intersected_tris),
        CGAL::parameters::default_values());

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
test_faces_polylines_intersections(const std::string filename1,
                                  const std::string filename2,
                                  const bool expected)
{
  typedef typename K::Point_3                                    Point;
  typedef typename CGAL::Surface_mesh<Point>                     Mesh;

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
  std::vector<std::vector<Point> > points;
  if(load_polylines(input2, points) >0)
    return 1;

  std::cout << "Reading files: " << filename1 <<", "<< filename2 << std::endl;

  CGAL::Timer timer;
  timer.start();

  std::vector<std::pair<std::size_t,
      std::pair<std::size_t, std::size_t> > > intersected_tris;

  CGAL::Polygon_mesh_processing::internal::compute_face_polylines_intersection(
        faces(m),
        points,
        m,
        std::back_inserter(intersected_tris),
        CGAL::parameters::default_values());

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
test_polylines_polylines_intersections(const std::string filename1,
                                  const std::string filename2,
                                  const bool expected)
{
  typedef typename K::Point_3                                    Point;

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
  std::vector<std::vector<Point> > polys1;
  if(load_polylines(input1, polys1) >0)
    return 1;
  std::vector<std::vector<Point> > polys2;
  if(load_polylines(input2, polys2) >0)
    return 1;

  std::cout << "Reading files: " << filename1 <<", "<< filename2 << std::endl;

  CGAL::Timer timer;
  timer.start();

  std::vector<
        std::pair<
          std::pair<std::size_t, std::size_t>,
          std::pair<std::size_t, std::size_t>
        >
      > intersected_segs;

  CGAL::Polygon_mesh_processing::internal::compute_polylines_polylines_intersection(
        polys1,
        polys2,
        std::back_inserter(intersected_segs),
        K());

  bool intersecting_1 = !intersected_segs.empty();

  std::cout << "intersections test took " << timer.time() << " sec." << std::endl;
  std::cout << intersected_segs.size() << " intersections." << std::endl;

  timer.reset();
  bool intersecting_2 = CGAL::Polygon_mesh_processing::do_intersect(polys1, polys2);

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
test_polylines_intersections(const std::string filename1,
                                  const std::string filename2,
                                  const bool expected)
{
  typedef typename K::Point_3                                    Point;

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
  CGAL::Polygon_mesh_processing::internal::compute_polyline_polyline_intersection(
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

template <typename K>
int test_inter_in_range(const std::vector<std::string>& filenames, std::size_t expected, bool volume)
{
  typedef typename K::Point_3                                    Point;
  typedef typename CGAL::Surface_mesh<Point>                     Mesh;

  // reading input meshes
  std::vector<Mesh> meshes(filenames.size());
  for (std::size_t i=0; i< filenames.size(); ++i)
  {
    std::ifstream input(filenames[i]);
    input >> meshes[i];
  }
  std::vector<std::pair<std::size_t, std::size_t> > output;
  CGAL::Polygon_mesh_processing::intersecting_meshes(meshes, std::back_inserter(output),
    CGAL::parameters::do_overlap_test_of_bounded_sides(volume));
  std::cout << output.size() <<" pairs."<<std::endl;
  if(output.size() != expected)
    return 1;
  return 0;
}

int main()
{

  bool expected = true;
  const std::string filename1 =  "data/tetra1.off";
  const std::string filename2 =  CGAL::data_file_path("meshes/reference_tetrahedron.off");
  const std::string filename3 =  "data/triangle.polylines.txt";
  const std::string filename4 =  "data/planar.polylines.txt";
  const std::string filename5 =  "data/tetra3_inter.polylines.txt";
  const std::string filename6 =  "data/polylines_inter.polylines.txt";
  const std::string filename7 =  "data/tetra2.off";
  const std::string filename8 =  "data/tetra4.off";
  const std::string filename9 =  "data/small_spheres.off";
  const std::string filename10 = "data/hollow_sphere.off";
  const std::string filename11 = CGAL::data_file_path("meshes/sphere.off");


  std::cout << "First test (Epic):" << std::endl;
  int r = test_faces_intersections<Epic>(filename1, filename2,  expected);

  std::cout << "Second test (Epec):" << std::endl;
  r += test_faces_intersections<Epec>(filename1, filename2, expected);
  std::cout << "Third test (Polyline):" << std::endl;
  r += test_faces_polyline_intersections<Epic>(filename1, filename3, expected);
  std::cout << "Fourth test (Polylines):" << std::endl;
  r += test_polylines_intersections<Epic>(filename3, filename4, expected);
  std::cout << "Fifth test (Polyline Range and Faces):" << std::endl;
  r += test_faces_polylines_intersections<Epic>(filename2, filename5, expected);
  std::cout << "Sixth test (Polyline Ranges):" << std::endl;
  r += test_polylines_polylines_intersections<Epic>(filename5, filename6, expected);
  std::cout << "Seventh test (number of intersecting meshes (surface) ):" << std::endl;
  std::vector<std::string> names;
  names.push_back(filename1);
  names.push_back(filename2);
  names.push_back(filename7);
  names.push_back(filename8);
  r += test_inter_in_range<Epic>(names, 4, false);

  names.clear();
  names.push_back(filename9);
  names.push_back(filename10);
  r += test_inter_in_range<Epic>(names, 0, true);
  names.clear();
  names.push_back(filename9);
  names.push_back(filename11);
  r += test_inter_in_range<Epic>(names, 1, true);
  return r;
}
