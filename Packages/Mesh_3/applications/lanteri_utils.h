template <class Tr>
bool
scan_edges_and_process(const Tr& tr,
                       std::string prefix = "",
                       // prefix to each line output
                       std::ostream* out_stream =
                         &std::cout
                       // output stream
                       )
{
  std::vector<Qualities> surface_edges_length;
  std::vector<Qualities> volume_edges_length;

  std::vector<double> surface_edges_length_max;   
  std::vector<double> volume_edges_length_max;

  for(typename Tr::Finite_edges_iterator fit = tr.finite_edges_begin();
      fit!=tr.finite_edges_end();
      ++fit)
  {
    const typename Tr::Vertex_handle& va = fit->first->vertex(fit->second);
    const typename Tr::Vertex_handle& vb = fit->first->vertex(fit->third);

    const double length = 
      CGAL::sqrt(CGAL::to_double(squared_distance(va->point(),
                                                  vb->point())));

    const unsigned int& index_a = va->point().surface_index();
    const unsigned int& index_b = vb->point().surface_index();

    if( index_a != 0 && index_a == index_b ) // surface edge
    {
      if( surface_edges_length.size() <= index_a )
      {
        surface_edges_length.resize(index_a+1);
        surface_edges_length_max.resize(index_a+1);
      }
      if( length > surface_edges_length_max[index_a] )
        surface_edges_length_max[index_a] = length;
        
      surface_edges_length[index_a].push_back(length);
    }
    else                                     // volume edge
    {
      const int index = fit->first->volume_index();
      
      if(index >=0)
      {
        const unsigned int positive_index = index;
        if( volume_edges_length.size() <= positive_index )
        {
          volume_edges_length.resize(positive_index+1);
          volume_edges_length_max.resize(positive_index+1);
        }
        if( length > volume_edges_length_max[positive_index] )
          volume_edges_length_max[positive_index] = length;
        volume_edges_length[positive_index].push_back(length);
      }
    }
  }

  const typename Qualities::size_type surface_vector_size = 
    surface_edges_length_max.size();
  for(unsigned int i = 0; i < surface_vector_size; ++i)
  {
    *out_stream << prefix
                << "Maximum length for edges on surface #" << i << ": "
                << surface_edges_length_max[i] << std::endl;
  }

  const typename Qualities::size_type volume_vector_size = 
    volume_edges_length_max.size();
  for(unsigned int i = 0; i < volume_vector_size; ++i)
  {
    *out_stream << prefix
                << "Maximum length for edges in volume #" << i << ": "
                << volume_edges_length_max[i] << std::endl;
  }

  return process_surface_edges(surface_edges_length) &&
    process_volume_edges(volume_edges_length);
}

template <class Tr>
bool
scan_cells_and_process(const Tr& tr)
{
  std::vector<Qualities> volume_cells_quality;

  for(typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
      cit != tr.finite_cells_end();
      ++cit)
    if(cit->is_in_domain())
    {
      const double quality = 
        CGAL::to_double(radius_ratio(tr.tetrahedron(cit)));
      // radius ratio is in common namespace, in Slivers_exuder.h
      int index = cit->volume_index();
      if(index < 0)
        index = 0;

      const unsigned int positive_index = index;

      if( volume_cells_quality.size() <= positive_index )
        volume_cells_quality.resize(positive_index+1);
      volume_cells_quality[positive_index].push_back(quality);
    }

  return process_cells(volume_cells_quality);
}
