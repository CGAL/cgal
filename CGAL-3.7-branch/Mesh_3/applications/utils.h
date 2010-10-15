#include <iostream>
#include <string>
#include <vector>
#include <CGAL/Complex_2_in_triangulation_3.h>

template <class Tr>
void display_faces_counts(const Tr& tr,
                          std::string prefix = "",
                          // prefix to each line output
                          std::ostream* out_stream = &std::cout
                          // output stream
                          )
{
  *out_stream << prefix << "Vertices: " << tr.number_of_vertices() << std::endl
              << prefix << "Facets on surface: " 
              << CGAL::Surface_mesher::number_of_facets_on_surface(tr) 
              << std::endl
              << prefix << "Cells: "
              << tr.number_of_cells() << std::endl
              << prefix << "Cells in volume: " 
              << CGAL::number_of_cells_in_domain(tr) << std::endl;
}

template <class C2T3>
void
display_facets_by_surface_indices_statistics(C2T3& c2t3,
                                             std::string prefix = "",
                                             // prefix to each line output
                                             std::ostream* out_stream =
                                               &std::cout
                                             // output stream
                                             )
{
  typedef typename C2T3::Triangulation Tr;

  const Tr& tr = c2t3.triangulation();

  std::vector<int> number_of_facets_by_surface_index;
  int number_of_hybrid_facets = 0;

  for(typename Tr::Finite_facets_iterator fit = tr.finite_facets_begin();
      fit != tr.finite_facets_end();
      ++fit)
  {
    const typename Tr::Cell_handle& cell = fit->first;
    const int index = fit->second;

    if(c2t3.face_status(cell, index) != C2T3::NOT_IN_COMPLEX)
    {
      const typename Tr::Vertex_handle& va = cell->vertex((index+1)&3);
      const typename Tr::Vertex_handle& vb = cell->vertex((index+1)&3);
      const typename Tr::Vertex_handle& vc = cell->vertex((index+1)&3);

      const unsigned int va_index = va->point().surface_index();
      const unsigned int vb_index = vb->point().surface_index();
      const unsigned int vc_index = vc->point().surface_index();

      if(va_index != vb_index || va_index != vc_index )
        ++number_of_hybrid_facets;
      else
      {
        if(number_of_facets_by_surface_index.size() <= va_index)
          number_of_facets_by_surface_index.resize(va_index+1);
        ++number_of_facets_by_surface_index[va_index];
      }
    }
  }
    
  *out_stream << prefix << "Hybrid facets: "
              << number_of_hybrid_facets << std::endl;

  std::vector<int>::size_type vector_size = 
    number_of_facets_by_surface_index.size();
  for(unsigned int i = 0; i < vector_size; ++i)
  {
    *out_stream << prefix << "Facets on surface #" << i << ": "
                << number_of_facets_by_surface_index[i] << std::endl;
  }
}

template <class Tr>
void
display_vertices_by_surface_indices_statistics(const Tr& tr,
                                               std::string prefix = "",
                                               // prefix to each line output
                                               std::ostream* out_stream =
                                                 &std::cout
                                               // output stream
                                               )
{
  std::vector<int> number_of_vertices_by_surface_index;

  for(typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
      vit != tr.finite_vertices_end();
      ++vit)
  {
    const unsigned int index = vit->point().surface_index();
    
    if(number_of_vertices_by_surface_index.size() <= index)
      number_of_vertices_by_surface_index.resize(index+1);
    ++number_of_vertices_by_surface_index[index];
  }
    
  std::vector<int>::size_type vector_size = 
    number_of_vertices_by_surface_index.size();
  for(unsigned int i = 0; i < vector_size; ++i)
  {
    *out_stream << prefix << "Vertices with index #" << i << ": "
                << number_of_vertices_by_surface_index[i] << std::endl;
  }
}

template <class Tr>
void
display_cells_by_volume_indices_statistics(const Tr& tr,
                                           std::string prefix = "",
                                           // prefix to each line output
                                           std::ostream* out_stream =
                                              &std::cout
                                           // output stream
                                           )
{
  std::vector<int> number_of_cells_by_volume_index;
  int number_of_cells_with_default_volume_index = 0;

  for(typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
      cit != tr.finite_cells_end();
      ++cit)
    if(cit->is_in_domain())
    {
      const int index = cit->volume_index();

      if(index < 0)
        ++number_of_cells_with_default_volume_index;
      else
      {
        const unsigned int positive_index = index;

        if(number_of_cells_by_volume_index.size() <= positive_index)
          number_of_cells_by_volume_index.resize(index+1);
        ++number_of_cells_by_volume_index[index];
      }
    }

  *out_stream << prefix << "Cells with default index: "
              << number_of_cells_with_default_volume_index << std::endl;
    
  std::vector<int>::size_type vector_size = 
    number_of_cells_by_volume_index.size();
  for(unsigned int i = 0; i < vector_size; ++i)
  {
    *out_stream << prefix << "Cells with index #" << i << ": "
                << number_of_cells_by_volume_index[i] << std::endl;
  }
}
