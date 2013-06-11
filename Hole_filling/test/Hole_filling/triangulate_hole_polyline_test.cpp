#include <CGAL/Hole_filling.h>
#include <CGAL/Simple_cartesian.h>

#include <cassert>
#include <vector>
#include <boost/tuple/tuple.hpp>

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Simple_cartesian<double>   Kernel;
typedef Kernel::Point_3                  Point_3;
typedef CGAL::Polyhedron_3<Kernel>       Polyhedron;

// to debug visually, construct polyhedron from patch
template<class HDS>
class Polyhedron_builder : public CGAL::Modifier_base<HDS> {
public:
  Polyhedron_builder(std::vector<boost::tuple<int, int, int> >* triangles, 
    std::vector<Point_3>* polyline) 
    : triangles(triangles), polyline(polyline) 
  { }

  void operator()(HDS& hds) {
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
    B.begin_surface(polyline->size() -1, triangles->size());

    for(std::vector<Point_3>::iterator it = polyline->begin();
      it != --polyline->end(); ++it) {
        B.add_vertex(*it);
    }

    for(std::vector<boost::tuple<int, int, int> >::iterator it = triangles->begin();
      it != triangles->end(); ++it) {
        B.begin_facet();
        B.add_vertex_to_facet(it->get<0>());
        B.add_vertex_to_facet(it->get<1>());
        B.add_vertex_to_facet(it->get<2>());
        B.end_facet();
    }

    B.end_surface();
  }

private:
  std::vector<boost::tuple<int, int, int> >* triangles;
  std::vector<Point_3>* polyline;
};


// it reads .polylines.txt (there should be one polyline with last point repeated)
bool read_polyline_one_line(const char* file_name, std::vector<Point_3>& points) {
  std::ifstream stream(file_name);
  if(!stream) { return false; }

  int count;
  if(!(stream >> count)) { return false; }
  while(count-- > 0) {
    Point_3 p;
    if(!(stream >> p)) { return false; }
    points.push_back(p);
  }
  return true;
}

// last point should be be repeated
bool read_polyline_with_extra_points(const char* file_name, 
  std::vector<Point_3>& points,
  std::vector<Point_3>& extras)
{
  std::ifstream stream(file_name);
  if(!stream) { return false; }

  for(int i =0; i < 2; ++i) {
    int count;
    if(!(stream >> count)) { return false; }
    while(count-- > 0) {
      Point_3 p;
      if(!(stream >> p)) { return false; }
      i == 0 ? points.push_back(p) : extras.push_back(p);
    }
  }
  return true;
}

bool check_triangles(std::vector<Point_3>& points, std::vector<boost::tuple<int, int, int> >& tris) {
  if(points.size() - 3 != tris.size()) {
    std::cerr << "Error: there should be n-2 triangles in generated patch." << std::endl;
    return false;
  }

  const int max_index = points.size()-1;
  for(std::vector<boost::tuple<int, int, int> >::iterator it = tris.begin(); it != tris.end(); ++it) {
    if(it->get<0>() == it->get<1>() ||
      it->get<0>() == it->get<2>() ||
      it->get<1>() == it->get<2>() ) 
    {
      std::cerr << "Error: indices of triangles should be all different." << std::endl;
      return false; 
    }  

    if(it->get<0>() >= max_index ||
      it->get<1>() >= max_index ||
      it->get<2>() >= max_index ) 
    {
      std::cerr << "Error: max possible index check failed." << std::endl;
      return false; 
    } 
  }

  return true;
}

bool check_constructed_polyhedron(const char* file_name,
  std::vector<boost::tuple<int, int, int> >* triangles, 
  std::vector<Point_3>* polyline) 
{
  Polyhedron poly;
  Polyhedron_builder<Polyhedron::HalfedgeDS> patch_builder(triangles, polyline);
  poly.delegate(patch_builder);

  if(!poly.is_valid()) {
    std::cerr << "Error: constructed patch does not constitute a valid polyhedron." << std::endl;
    return false;
  }
  std::string out_file_name;
  out_file_name.append(file_name).append(".off");
  std::ofstream out(out_file_name);
  out << poly; out.close();
  return true;
}

bool test_1(const char* file_name) {
  std::cerr << "test_1 with '" << file_name << "' file..." << std::endl;
  std::vector<Point_3> points; // this will contain n and +1 repeated point
  if(!read_polyline_one_line(file_name, points)) {
    std::cerr << "Error: can not read file." << std::endl;
    return false;
  }

  std::vector<boost::tuple<int, int, int> > tris;
  CGAL::triangulate_hole_polyline(points.begin(), --points.end(), std::back_inserter(tris));

  if(!check_triangles(points, tris)) { return false; }
  if(!check_constructed_polyhedron(file_name, &tris, &points)) { return false; }

  std::cerr << "Done!" << std::endl;
  return true;
}

bool test_2(const char* file_name) {
  std::cerr << "test_2 with '" << file_name << "' file..." << std::endl;
  std::vector<Point_3> points; // this will contain n and +1 repeated point
  std::vector<Point_3> extras;
  if(!read_polyline_with_extra_points(file_name, points, extras)) {
    std::cerr << "Error: can not read file." << std::endl;
    return false;
  }

  std::vector<boost::tuple<int, int, int> > tris;
  CGAL::triangulate_hole_polyline(points.begin(), points.end(),
    extras.begin(), extras.end(), std::back_inserter(tris));

  if(!check_triangles(points, tris)) { return false; }
  if(!check_constructed_polyhedron(file_name, &tris, &points)) { return false; }
  std::cerr << "Done!" << std::endl;
  return true;
}

int main() {
  std::vector<std::string> input_files_1;
  input_files_1.push_back("data/triangle.polylines.txt");
  input_files_1.push_back("data/quad.polylines.txt");
  input_files_1.push_back("data/U.polylines.txt");

  for(std::vector<std::string>::iterator it = input_files_1.begin(); it != input_files_1.end(); ++it) {
    assert(test_1(it->c_str()));
  }

  std::vector<std::string> input_files_2;
  input_files_2.push_back("data/hole1.txt");
  input_files_2.push_back("data/hole2.txt");
  input_files_2.push_back("data/hole3.txt");
  input_files_2.push_back("data/hole4.txt");

  for(std::vector<std::string>::iterator it = input_files_2.begin(); it != input_files_2.end(); ++it) {
    assert(test_2(it->c_str()));
  }
  std::cerr << "All Done!" << std::endl;
}