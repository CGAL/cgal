#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Simple_cartesian.h>
#include <CGAL/MP_Float.h>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Constrained_regular_triangulation_3.h>
#include <CGAL/Constrained_triangulation_vertex_base_3.h>
#include <CGAL/Constrained_triangulation_cell_base_3.h>

#include <CGAL/Triangulation_2_traits_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
//#include <CGAL/Conforming_Delaunay_triangulation_2.h>
#include <CGAL/Conforming_Delaunay_triangulation_2_traits_3.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>

//#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>

//typedef CGAL::Quotient<CGAL::MP_Float> FT;
struct K : public CGAL::Exact_predicates_inexact_constructions_kernel {};
//struct K : public CGAL::Simple_cartesian<FT> {};

typedef CGAL::Conforming_Delaunay_triangulation_2_traits_3<K> Traits_2;
typedef CGAL::Triangulation_vertex_base_2<Traits_2> Vb_2;
typedef CGAL::Delaunay_mesh_face_base_2<Traits_2> Fb_2;
typedef CGAL::Triangulation_data_structure_2<Vb_2, Fb_2> Tds_2;
typedef CGAL::Constrained_Delaunay_triangulation_2<Traits_2, Tds_2> CDT_2;
typedef CGAL::Constrained_triangulation_plus_2<CDT_2> CDTP_2;
//typedef CGAL::Conforming_Delaunay_triangulation_2<CDTP_2> Conf_DT_2;

typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Traits_3;
typedef CGAL::Constrained_triangulation_cell_base_3<Traits_3> Cb_3;
typedef CGAL::Constrained_triangulation_vertex_base_3<Traits_3> Vb_3;
typedef CGAL::Triangulation_data_structure_3<Vb_3, Cb_3> Tds_3;

typedef Traits_3::RT                                          Weight;
typedef Traits_3::Bare_point                                  Point_3;
typedef Traits_3::Weighted_point                              Weighted_point_3;

typedef CGAL::Regular_triangulation_3<Traits_3, Tds_3>          Rt_3;
typedef CGAL::Constrained_regular_triangulation_3<Rt_3, CDT_2>  Crt;

typedef Crt::Vertex_iterator                                 Vertex_iterator;
typedef Crt::Vertex_handle                                   Vertex_handle;
typedef Crt::Cell_handle                                   Cell_handle;


int main(int argc, char** argv)
{
  Crt T;

  // insertion of points on a 3D grid
//   std::vector<Vertex_handle> V;

//   for (int z=0 ; z<5 ; z++)
//     for (int y=0 ; y<5 ; y++)
//       for (int x=0 ; x<5 ; x++) {
// 	  Point_3 p(x, y, z);
//           Weight w = (x+y-z*y*x)*2.0; // let's say this is the weight.
// 	  Weighted_point_3 wp(p, w);
// 	  V.push_back(T.insert(wp));
//       }

//   //  assert( T.is_valid() );
//   assert( T.dimension() == 3 );

//   std::cout << "Number of vertices : " << T.number_of_vertices() << std::endl;

//   T.clear();

  std::cerr.precision(20);

  if ( argc == 2 )
    {
      typedef CGAL::PLC_loader<Crt, CDT_2> Loader;

      std::ifstream f(argv[1]);

      T.clear();

      Loader(T).load_triangulation(f);
    }

//   T.insert(Point_3(1,1,0));
//   T.insert(Point_3(-1,1,0));
//   T.insert(Point_3(-1,-1,0));
//   T.insert(Point_3(1,-1,0));
//   Vertex_handle va = T.insert(Point_3(0,0,-0.1));
//   Vertex_handle vb = T.insert(Point_3(0,0,0.1));
  
  //  int i1, i2;
  //  Cell_handle c;

  //  CGAL_assertion( T.is_edge(va, vb, c, i1, i2) );
  T.off_file_output(std::cout);
}
