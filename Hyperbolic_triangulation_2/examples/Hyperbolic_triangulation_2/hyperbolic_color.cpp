#include <fstream>

// CGAL headers
#include <CGAL/IO/io.h>

#include <CGAL/Exact_circular_kernel_2.h>

#include <CGAL/Hyperbolic_triangulation_face_base_with_info_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>

#include <CGAL/IO/Color.h>

typedef CGAL::Exact_circular_kernel_2 K;
typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2< K > Gt;

typedef Gt::Point_2 Point_2;

typedef CGAL::Hyperbolic_triangulation_face_base_with_info_2<CGAL::Color, Gt> Fb;
typedef CGAL::Triangulation_data_structure_2 <
                        CGAL::Triangulation_vertex_base_2<Gt>, Fb > TDS;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Gt, TDS> Dt;

int main()
{
  std::vector<Point_2> pts;
  Point_2 p;

  std::ifstream ifs("input-file");
  while(ifs >> p) {
    pts.push_back(p);
  }

  Dt dt;
  
  dt.insert(pts.begin(),pts.end());
  Dt::Vertex_handle vo = dt.insert(Point_2(0,0));
  
  int origin_faces = 0;
  Dt::Hyperbolic_faces_iterator fit;
  for (fit = dt.hyperbolic_faces_begin(); fit != dt.hyperbolic_faces_end(); ++fit)
    if (fit->has_vertex(vo))
      {
          fit->info() = CGAL::RED;
	  origin_faces++;
      }

  int red_faces = 0;
  for (fit = dt.hyperbolic_faces_begin(); fit != dt.hyperbolic_faces_end(); ++fit)
    if (fit->info() == CGAL::RED)
      red_faces++;

  assert(red_faces == origin_faces);

  std::cout << "number of points " << std::distance(pts.begin(),pts.end())+1 << std::endl;
  std::cout << "Number of (finite) vertices: " << dt.number_of_vertices() << std::endl;
  std::cout << "number of (finite) Euclidean faces: " << dt.number_of_faces() << std::endl;
  std::cout << "number of hyperbolic faces: " << dt.number_of_hyperbolic_faces() << std::endl;
  std::cout << "number of faces having the origin as vertex: " << origin_faces << std::endl;

  return 0;
}
