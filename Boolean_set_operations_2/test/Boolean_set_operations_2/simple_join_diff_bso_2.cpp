#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_set_2.h>
#include <list>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                                   Point_2;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2;
typedef CGAL::Polygon_set_2<Kernel>                       Polygon_set_2;

void print_polygon(std::ostream& out, const Polygon_2& polygon)
{
  out << polygon.size()+1 << " ";
  for (Polygon_2::Vertex_const_iterator vit=polygon.vertices_begin();
                                        vit!=polygon.vertices_end();++vit)
    out << *vit <<" 0 ";
  out << *polygon.vertices_begin() << " 0\n";
}

void print_polygons(std::ostream& out, const Polygon_set_2& polygon_set)
{
  std::vector<Polygon_with_holes_2> polygons_wh(polygon_set.number_of_polygons_with_holes());
  polygon_set.polygons_with_holes(&polygons_wh[0]);
  for(Polygon_with_holes_2& polygon_wh : polygons_wh)
  {
    print_polygon(out, polygon_wh.outer_boundary());
    for(Polygon_with_holes_2::Hole_const_iterator it=polygon_wh.holes_begin();
                                                  it!=polygon_wh.holes_end();++it)
    {
      print_polygon(out, *it);
    }
  }
}


void join(const char* filename)
{
  std::cout << "Running join on " << filename << "\n";
  std::ifstream input_file(filename);
  if (!input_file.is_open()) {
    std::cerr << "Failed to open the " << filename <<std::endl;
    exit(-1);
  }

  Polygon_set_2                   ps;
  std::list<Polygon_with_holes_2> pwhs;
  size_t                          num_pwhs;

  input_file >> num_pwhs;
  for ( unsigned int i=0; i<num_pwhs; i++ ) {
      Polygon_with_holes_2 pwh;
      input_file >> pwh;
      pwhs.push_back( pwh );
  }
  input_file.close();
  ps.join(pwhs.begin(), pwhs.end());

  //std::ofstream out("out.cgal");
  //print_polygons(out, ps);
}

void difference(const char* filename)
{
  std::cout << "Running difference on " << filename << "\n";
  std::ifstream input_file(filename);
  if (!input_file.is_open()) {
    std::cerr << "Failed to open the " << filename <<std::endl;
    exit(-1);
  }

  Polygon_set_2                   ps;
  std::list<Polygon_with_holes_2> pwhs;
  size_t                          num_pwhs;

  input_file >> num_pwhs;

  for ( unsigned int i=0; i<num_pwhs; i++ ) {
      Polygon_with_holes_2 pwh;
      input_file >> pwh;
      pwhs.push_back( pwh );
  }
  input_file.close();

  std::list<Polygon_with_holes_2>::const_iterator it = pwhs.begin();
  ps.join(*it);
  ++it;
  while ( it != pwhs.end() ) {
    ps.difference( (*it) );
    ++it;
  }
  std::ofstream out("out.cgal");
  print_polygons(out, ps);
}

int main()
{
  join("data/union_01.txt");
  join("data/union_02.txt");
  join("data/union_03.txt");
  difference("data/diff_01.txt");
  return 0;
}
