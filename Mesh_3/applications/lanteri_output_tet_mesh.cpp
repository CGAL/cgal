#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// regular
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>

// IO.h must be included before vertex and cell bases.
#include <CGAL/Mesh_3/IO.h>

// vertex and cell bases
#include <CGAL/Surface_mesh_vertex_base_3.h>
#include <CGAL/Surface_mesh_cell_base_3.h>
#include <CGAL/Volume_mesher_cell_base_3.h>

// c2t3
#include <CGAL/Complex_2_in_triangulation_3.h>

// traits class for reading meshes
#include <CGAL/Weighted_point_with_surface_index_geom_traits.h>

// radius_ratio
#include <CGAL/Mesh_3/Slivers_exuder.h>

// Point_traits
#include <CGAL/Point_traits.h>

#include <sstream>

// traits class
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_filtered_traits_3<K> Regular_traits;
typedef CGAL::Weighted_point_with_surface_index_geom_traits<Regular_traits> My_traits;

// vertex and cell types
typedef CGAL::Surface_mesh_vertex_base_3<My_traits> Vb;
typedef CGAL::Regular_triangulation_cell_base_3<My_traits> Cb1;
typedef CGAL::Surface_mesh_cell_base_3<My_traits, Cb1> Cb2;
typedef CGAL::Volume_mesher_cell_base_3<My_traits, Cb2> Cb3;

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
  Cell_with_volume_index() : Cb(), volume(0)
  {
  }
  Cell_with_volume_index (Vertex_handle v0,
                          Vertex_handle v1,
                          Vertex_handle v2,
                          Vertex_handle v3)
    :  Cb (v0, v1, v2, v3), volume(0)
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
    : Cb (v0, v1, v2, v3, n0, n1, n2, n3), volume(0) 
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

#ifdef CGAL_MESH_3_IO_H
  static
  std::string io_signature()
  {
    return CGAL::Get_io_signature<Cb>()();
  }
#endif
}; // end template class Cell_with_volume_index

template <class Gt, class Cb>
inline
std::istream&
operator>>(std::istream &is, Cell_with_volume_index<Gt, Cb>& c)
{
  return is >> static_cast<Cb&>(c);
}

template <class Gt, class Cb>
inline
std::ostream&
operator<<(std::ostream &os, const Cell_with_volume_index<Gt, Cb>& c)
{
  return os <<  static_cast<const Cb&>(c);
}

typedef Cell_with_volume_index<My_traits, Cb3> Cb;

// triangulation
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Regular_triangulation_3<My_traits, Tds> Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2T3;

// ios
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>


int main(int , char**)
{
  Tr tr;
  C2T3 c2t3(tr);

  double r1, r2, r3, r4, r5;
  std::vector<double> size_bounds(5);
  std::vector<double> radii(5);
  
  std::cout << "Input r1, r2, r3, r4, r5:" << std::endl;
  std::cin >> r1 >> r2 >> r3 >> r4 >> r5;
  std::cout << "Input the corresponding 5 size bounds:" << std::endl;
  std::cin >> size_bounds[0]
           >> size_bounds[1]
           >> size_bounds[2]
           >> size_bounds[3]
           >> size_bounds[4];
  if(!std::cin)
    return EXIT_FAILURE;

  std::string file_prefix;
  std::cout << "Input filename prefix:" << std::endl;
  std::cin >> file_prefix;

  // input file
  std::ifstream ifs((file_prefix + ".cgal").c_str());

  // ouput tets mesh file
  std::ofstream ofs_maillage((file_prefix + ".maillage").c_str());

  // output five surface meshes files
  std::vector<std::ofstream*> ofs_surfaces(6); 
  // ofs_surfaces[0] will not be used.

  for(int i = 1; i <= 5; ++i)
  {
    std::stringstream str_stream;
    str_stream << file_prefix << ".surface" << i;

    ofs_surfaces[i] = new std::ofstream(str_stream.str().c_str());
  }

  std::cout << "  Reading file " << (file_prefix + ".cgal") << std::endl;

  if( CGAL::Mesh_3::input_mesh(ifs, c2t3,
                               true,
                               &std::cerr) )
  {
    std::cout << "  Writing file " << (file_prefix + ".maillage") << std::endl
              << "  and files " << (file_prefix + ".surface*") << std::endl;

    typedef C2T3::Triangulation Tr;
    typedef Tr::Finite_cells_iterator Finite_cells_iterator;
    typedef Tr::Finite_facets_iterator Finite_facets_iterator;
    typedef Tr::Finite_vertices_iterator Finite_vertices_iterator;
    typedef Tr::Vertex_handle Vertex_handle;
    typedef Tr::Point Point;
    typedef CGAL::Point_traits<Point> P_traits;
    typedef P_traits::Bare_point Bare_point;

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
          cit->set_volume_index(1); // brain
        else if( sq_r < r2*r2 )
          cit->set_volume_index(2); // LCR
        else if( sq_r < r3*r3 )
          cit->set_volume_index(3); // head skull
        else if( sq_r < r4*r4 )
          cit->set_volume_index(4); // skin
        else if( sq_r < r5*r5 )
          cit->set_volume_index(-1);// air
      }

    const Tr& tr = c2t3.triangulation();

    // Headers

    ofs_maillage << tr.number_of_vertices() << " "
                 << CGAL::number_of_cells_in_domain(tr) << " "
                 << CGAL::number_of_facets_on_surface_with_index(c2t3, 5)
                 << std::endl;

    for(int i = 1; i <= 5; ++i)
    {
      *ofs_surfaces[i] << tr.number_of_vertices() << " "
                       << CGAL::number_of_facets_on_surface_with_index(c2t3, i)
                       << std::endl;
    }

    // precision
    ofs_maillage << std::setprecision(20);
    for(int i = 1; i <= 5; ++i)
    {        
      *ofs_surfaces[i] << std::setprecision(20);
    }
 
    // Vertices

    std::map<Vertex_handle, int> V; // vertices are counter from 1
    int inum = 1;
    for( Finite_vertices_iterator vit = tr.finite_vertices_begin();
         vit != tr.finite_vertices_end();
         ++vit)
    {
      V[vit] = inum++;
      Point p = static_cast<Point>(vit->point());
      ofs_maillage << p.x() << " " << p.y() << " " << p.z()
                   << std::endl;
      for(int i = 1; i <= 5; ++i)
      {
        *ofs_surfaces[i] << p.x() << " " << p.y() << " " << p.z()
                         << std::endl;
      }
    }

    // Tetrahedra

    for( Finite_cells_iterator cit = tr.finite_cells_begin(); 
       cit != tr.finite_cells_end(); ++cit)
      if( cit->is_in_domain() )
      {
        for (int i=0; i<4; i++)
          ofs_maillage << V[cit->vertex(i)] << " ";
        
        ofs_maillage << cit->volume_index() 
                     << std::endl;
      }
  
    // Facets
    for( Finite_facets_iterator fit = tr.finite_facets_begin(); 
         fit != tr.finite_facets_end(); ++fit)
    {
      int surface_index = 0;
      if (c2t3.face_status(fit->first,fit->second)
          != C2T3::NOT_IN_COMPLEX)
      {
        for (int i=0; i<4; i++)
          if (i != (*fit).second)
          {
            const Vertex_handle& vh = (*fit).first->vertex(i);
            surface_index = vh->point().surface_index();
            if(surface_index == 5)
            {
              ofs_maillage << V[vh] << " ";
            }
            *ofs_surfaces[surface_index] << V[vh] << " ";
          }
        if(surface_index == 5)
          ofs_maillage << "4"
                       << std::endl;
        *ofs_surfaces[surface_index] << surface_index
                                     << std::endl;
      }
    }

    for(int i = 1; i <= 5; ++i)
      if( ofs_surfaces[i]->bad() )
        return EXIT_FAILURE;
    
    if( ofs_maillage.good() )
    {
      for(int i = 1; i <= 5; ++i)
        delete ofs_surfaces[i];
      return EXIT_SUCCESS;
    }
    else 
      return EXIT_FAILURE;
  }
  else
    return EXIT_FAILURE;
}
