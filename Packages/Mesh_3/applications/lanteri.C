#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// vertex and cell bases
#include <CGAL/Complex_2_in_triangulation_vertex_base_3.h>
#include <CGAL/Complex_2_in_triangulation_cell_base_3.h>
#include <CGAL/Mesh_3/Complex_2_in_triangulation_cell_base_3.h>

// c2t3
#include <CGAL/Complex_2_in_triangulation_3.h>

// delaunay
#include <CGAL/Delaunay_triangulation_3.h>

// traits class for reading meshes
#include <CGAL/Point_with_surface_index_geom_traits.h>

// radius_ratio
#include <CGAL/Mesh_3/Slivers_exuder.h>

// traits class
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Point_with_surface_index_geom_traits<K> My_traits;

// vertex and cell types
typedef CGAL::Triangulation_vertex_base_3<My_traits> Vb1;
typedef CGAL::Complex_2_in_triangulation_vertex_base_3<My_traits, Vb1> Vb;
typedef CGAL::Triangulation_cell_base_3<My_traits> Cb1;
typedef CGAL::Complex_2_in_triangulation_cell_base_3<My_traits, Cb1> Cb2;
typedef CGAL::Mesh_3::Complex_2_in_triangulation_cell_base_3<My_traits, Cb2> Cb3;

template <class GT, class Cb = CGAL::Triangulation_ds_cell_base_3 <> >
class Cell_with_volume_index : public Cb
{
private:
  int volume;

public:
  typedef typename GT::Point_3 Point;
  typedef typename Cb::Triangulation_data_structure Tds;
  typedef typename Tds::Vertex_handle Vertex_handle;
  typedef typename Tds::Cell_handle Cell_handle;

  template < class TDS3 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS3>::Other  Cb3;
    typedef Cell_with_volume_index<GT, Cb3> Other;
  };
    
  // Constructors
  Cell_with_volume_index() : Cb(), volume(-1)
  {
  }
  Cell_with_volume_index (Vertex_handle v0,
                          Vertex_handle v1,
                          Vertex_handle v2,
                          Vertex_handle v3)
    :  Cb (v0, v1, v2, v3), volume(-1)
  {
  }
  Cell_with_volume_index (Vertex_handle v0,
                          Vertex_handle v1,
                          Vertex_handle v2,
                          Vertex_handle v3,
                          Cell_handle n0,
                          Cell_handle n1,
                          Cell_handle n2,
                          Cell_handle n3)
    : Cb (v0, v1, v2, v3, n0, n1, n2, n3), volume(-1) 
  {
  }

  // volume index
  int volume_index() const
  {
    return volume;
  }
      
  void set_volume_index(const int i)
  {
    volume = i;
  }
}; // end template class Cell_with_volume_index

typedef Cell_with_volume_index<My_traits, Cb3> Cb;

// triangulation
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;

typedef CGAL::Delaunay_triangulation_3<My_traits, Tds> Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2T3;

// ios
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/File_medit.h>

#include "utils.h"
#include "distribution.h"

bool process_cells(const std::vector<Qualities>& volume_cells_quality);
bool process_volume_edges(const std::vector<Qualities>& volume_edges_lenght);
bool process_surface_edges(const std::vector<Qualities>& surface_edges_lenght);

#include "lanteri_utils.h"

typedef Tr::Vertex_handle Vertex_handle;

int main(int , char**)
{
  Tr tr;
  C2T3 c2t3(tr);

  double r1, r2, r3, r4, r5;
  std::cout << "Input r1, r2, r3, r4, r5:" << std::endl;
  std::cin >> r1 >> r2 >> r3 >> r4 >> r5;

  std::string filename;
  std::cout << "Input filename:" << std::endl;
  std::cin >> filename;
  std::ifstream ifs(filename.c_str());
  if( !ifs || !cin)
  {
    return EXIT_FAILURE;
  }

  std::cout << "  Reading " << filename << std::endl;
  if( ! CGAL::input_pslg_from_medit(ifs,
                                    c2t3,
                                    true,         // debug
                                    &std::cout) ) // debug to cout
    return EXIT_FAILURE;
  
  display_faces_counts(tr, "    ", &std::cout);

  std::cout << "\n  Statistics:\n";

  std::cout << "(vertices)\n";
  display_vertices_by_surface_indices_statistics(tr, "    ", &std::cout);

  std:: cout << "(facets)\n";
  display_facets_by_surface_indices_statistics(c2t3, "    ", &std::cout);

  // sets volume indices
  for(Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
      cit != tr.finite_cells_end();
      ++cit)
    if(cit->is_in_domain())
    {
      const double sq_r = 
        squared_distance(K::Point_3(0, 0, 0),
                         static_cast<K::Point_3>(tr.dual(cit)));
    
      if( sq_r < r1*r1 )
        cit->set_volume_index(1);
      else if( sq_r < r2*r2 )
        cit->set_volume_index(2);
      else if( sq_r < r3*r3 )
        cit->set_volume_index(3);
      else if( sq_r < r4*r4 )
        cit->set_volume_index(4);
      else if( sq_r < r5*r5 )
        cit->set_volume_index(5);
    }

  std::cout << "(cells)\n";
  display_cells_by_volume_indices_statistics(tr, "    ", &std::cout);

  std::cout << "\n"
            << "  Edges lengths:\n";
  if(!scan_edges_and_process(tr, "    ", &std::cout))
    return EXIT_FAILURE;
    
  if(!scan_cells_and_process(tr))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
