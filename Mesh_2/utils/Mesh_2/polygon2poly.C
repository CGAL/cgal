#include <CGAL/basic.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/Conforming_Delaunay_triangulation_2.h>
#include <CGAL/Read_write.h>

typedef CGAL::Simple_cartesian<double> K1;
typedef CGAL::Filtered_kernel<K1> K2;
struct K : public K2 {};

typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds,
  CGAL::Exact_predicates_tag> Tr;
typedef CGAL::Conforming_Delaunay_triangulation_2<Tr> Conform;
typedef Conform::Vertex_handle Vertex_handle;

typedef K::Point_2 Point;

Conform conform;

int main(int, char**)
{
  unsigned int number_of_points;
  std::cin >> number_of_points;
  CGAL::skip_until_EOL(std::cin);
  CGAL::skip_comment_OFF(std::cin);
  
  // read vertices
  std::vector<Vertex_handle> vertices(number_of_points);
  for(unsigned int i = 0; i < number_of_points; ++i)
    {
      Point p;
      std::cin >> p;
      CGAL::skip_until_EOL(std::cin); CGAL::skip_comment_OFF(std::cin);
      vertices[i] = conform.insert(p);
    }
  
  // les segments sont formes de 2 points consecutifs
  // read segments
  unsigned int k;
  for(k = 0; k < number_of_points-1; ++k)
    {
      conform.insert_constraint(vertices[k], vertices[k+1]);
    }

  // le dernier segment, qui relie le dernier point au premier
  conform.insert_constraint(vertices[k], vertices[0]);

  
  write_poly(conform, std::cout);
  return 0;
};
