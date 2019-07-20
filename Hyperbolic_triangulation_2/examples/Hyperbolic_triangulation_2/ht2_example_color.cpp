#include <CGAL/Hyperbolic_triangulation_face_base_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/IO/Color.h>

#include <iostream>
#include <vector>

typedef CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<>           Gt;
typedef CGAL::Hyperbolic_triangulation_face_base_2<Gt>                  Hyperbolic_face_base;
typedef CGAL::Triangulation_face_base_with_info_2<CGAL::Color, Gt,
                                                  Hyperbolic_face_base> Face_base_with_info;
typedef CGAL::Triangulation_vertex_base_2<Gt>                           Vertex_base;
typedef CGAL::Triangulation_data_structure_2<Vertex_base,
                                             Face_base_with_info>       TDS;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Gt, TDS>              Dt;
typedef Dt::Point                                                       Point_2;
typedef CGAL::Creator_uniform_2<Gt::FT, Point_2>                        Creator;

int main(int argc, char** argv)
{
  int N;
  if(argc < 2) {
    std::cout << "usage: " << argv[0] << " [number_of_points]" << std::endl;
    std::cout << "Defaulting to 100k points..." << std::endl;
    N = 100000;
  } else {
    N = atoi(argv[1]);
  }

  std::cout << "Number of points: " << N << std::endl;

  CGAL::Random_points_in_disc_2<Point_2, Creator> in_disc;
  std::vector<Point_2> pts;
  for(int i=0; i<N; ++i)
    pts.push_back(*(in_disc++));

  Dt dt;
  dt.insert(pts.begin(), pts.end());
  Dt::Vertex_handle vo = dt.insert(Point_2(0,0));

  int origin_faces = 0;
  Dt::Hyperbolic_faces_iterator fit;
  for(fit = dt.hyperbolic_faces_begin(); fit != dt.hyperbolic_faces_end(); ++fit)
  {
    if(fit->has_vertex(vo))
    {
      fit->info() = CGAL::red();
      origin_faces++;
    }
  }

  int red_faces = 0;
  for(fit = dt.hyperbolic_faces_begin(); fit != dt.hyperbolic_faces_end(); ++fit)
  {
    if(fit->info() == CGAL::red())
    {
      red_faces++;
    }
  }

  assert(red_faces == origin_faces);

  std::cout << "Number of points                        " << N << std::endl;
  std::cout << "Number of vertices:                     " << dt.number_of_vertices() << std::endl;
  std::cout << "Number of hyperbolic faces:             " << dt.number_of_hyperbolic_faces() << std::endl;
  std::cout << "Number of faces incident to the origin: " << origin_faces << std::endl;

  return EXIT_SUCCESS;
}
