#include "cgal_types.h"

#include <CGAL/Conforming_Delaunay_triangulation_2_traits_3.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>

//#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>

typedef CGAL::Conforming_Delaunay_triangulation_2_traits_3<K> Traits_2_3;
typedef CGAL::Special_mesh_traits_2<Rt_3, Traits_2_3> Traits_2;
typedef CGAL::Triangulation_vertex_base_2<Traits_2> Vb_2;
typedef CGAL::Delaunay_mesh_face_base_2<Traits_2> Fb_2;
typedef CGAL::Triangulation_data_structure_2<Vb_2, Fb_2> Tds_2;
typedef CGAL::Constrained_Delaunay_triangulation_2<Traits_2, Tds_2> CDT_2;
typedef CGAL::Constrained_triangulation_plus_2<CDT_2> CDTP_2;
//typedef CGAL::Conforming_Delaunay_triangulation_2<CDTP_2> Conf_DT_2;

typedef Traits_3::RT                                          Weight;
typedef Traits_3::Bare_point                                  Point_3;
typedef Traits_3::Weighted_point                              Weighted_point_3;

typedef CGAL::Constrained_regular_triangulation_3<Rt_3, CDT_2>  Crt;

typedef Crt::Vertex_iterator                                 Vertex_iterator;
typedef Crt::Vertex_handle                                   Vertex_handle;
typedef Crt::Cell_handle                                   Cell_handle;


int main(int argc, char** argv)
{
  Crt T;

  typedef CGAL::Mesh_3::Refine_edges<Crt> Edges_level;
  Edges_level edges_level(T);

  typedef CGAL::Mesh_3::Refine_facets<Crt> Facets_level;
  Facets_level facets_level(T, &edges_level);

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

      std::string arg = argv[1];

      std::ifstream f(arg.data());

      T.clear();

      if( arg.find(".off") != 0 )
	CGAL::Off_loader<Crt, CDT_2>(T).load_triangulation(f);
      else
	CGAL::PLC_loader<Crt, CDT_2>(T).load_triangulation(f);
    }

  typedef CGAL::Mesh_3::facets::Refine_facets_visitor<Crt, Facets_level>
    Facets_level_visitor;
  Facets_level_visitor facets_visitor(&facets_level);

  edges_level.scan_triangulation();
  edges_level.refine(CGAL::Null_mesh_visitor());

  facets_level.scan_triangulation();
  facets_level.refine(CGAL::Null_mesh_visitor());//facets_visitor);

  //  T.fill_edges_to_be_conformed();
  //  T.conform_edges();
  //  T.fill_facets_to_be_conformed();
  //  T.conform_facets();

  T.off_file_output(std::cout);
}
