#include <string>
#include <iomanip>

// analyse edges
template <class Tr>
bool
scan_edges_and_process(const Tr& tr,
                       std::vector<double> length_bounds,
                       std::string filename_prefix,
                       std::string prefix = "",
                       // prefix to each line output
                       std::ostream* out_stream =
                         &std::cout
                       // output stream
                       )
{
  // reminder: Qualities is std::vector<double> (voir "distribution.h")
  std::vector<Qualities> surface_edges_length;
  std::vector<Qualities> volume_edges_length;

  std::vector<double> surface_edges_length_max;   
  std::vector<double> volume_edges_length_max;

  std::vector<double> surface_edges_length_min;   
  std::vector<double> volume_edges_length_min;
  
  std::vector<double> surface_edges_length_sum;   
  std::vector<double> volume_edges_length_sum;

  std::vector<int> surface_edges_length_num;
  std::vector<int> volume_edges_length_num;

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
      // resize vectors
      if( surface_edges_length.size() <= index_a )
      {
        
        surface_edges_length.resize(index_a+1);
        
        surface_edges_length_max.resize(index_a+1);
                  
        surface_edges_length_min.resize(index_a+1, std::numeric_limits<double>::infinity());
        
        surface_edges_length_sum.resize(index_a+1);
        surface_edges_length_num.resize(index_a+1);
        
      }
      
      // update statistics
      if( length > surface_edges_length_max[index_a] )
        surface_edges_length_max[index_a] = length;
      
      if( length < surface_edges_length_min[index_a] )
        surface_edges_length_min[index_a] = length;
        
      surface_edges_length_sum[index_a] += length;
      
      surface_edges_length_num[index_a]++;
        
        
      surface_edges_length[index_a].push_back(length);
    }
    else                                     // volume edge
    {
      const int index = fit->first->volume_index();
      
      if(index >=0)
      {
        const unsigned int positive_index = index;

        // resize vectors
        if( volume_edges_length.size() <= positive_index )
        {
          volume_edges_length.resize(positive_index+1);
          volume_edges_length_max.resize(positive_index+1, 0);
          volume_edges_length_min.resize(positive_index+1, std::numeric_limits<double>::infinity());
          volume_edges_length_sum.resize(positive_index+1, 0);
          volume_edges_length_num.resize(positive_index+1, 0);
        }

        // update statistics
        if( length > volume_edges_length_max[positive_index] )
          volume_edges_length_max[positive_index] = length;
          
        if( length < volume_edges_length_min[positive_index] )
          volume_edges_length_min[positive_index] = length;
          
        volume_edges_length_sum[positive_index] += length;

        volume_edges_length_num[positive_index]++;
     
        
        volume_edges_length[positive_index].push_back(length);
      }
    }
  }

  *out_stream << std::setprecision(3)
              << std::setw(46) << "max"
              << std::setw(13) << "min"
              << std::setw(13) << "avg"
              << std::endl;

  const typename Qualities::size_type surface_vector_size = 
    surface_edges_length_max.size();
  for(unsigned int i = 0; i < surface_vector_size; ++i)
  {
    *out_stream << prefix
                << "length for edges on surface #" << i << ": "
                << std::setw(10) << surface_edges_length_max[i] << ",  "
                << std::setw(10) << surface_edges_length_min[i] << ",  " 
                << std::setw(10) << surface_edges_length_sum[i]/surface_edges_length_num[i]
                << std::endl;
  }

  *out_stream << std::endl;

  const typename Qualities::size_type volume_vector_size = 
    volume_edges_length_max.size();
  for(unsigned int i = 0; i < volume_vector_size; ++i)
  {
    *out_stream << prefix
                << "length for edges in volume #" << i << ":  "
                << std::setw(10) << volume_edges_length_max[i]  << ",  "
                << std::setw(10) << volume_edges_length_min[i]  << ",  "
                << std::setw(10) << volume_edges_length_sum[i]/volume_edges_length_num[i]                
                << std::endl;
  }

  return process_surface_edges(surface_edges_length,
                               length_bounds,
                               filename_prefix) &&
         process_volume_edges(volume_edges_length,
                               length_bounds,
                              filename_prefix);
}

// analyse cells
template <class Tr>
bool
scan_cells_and_process(const Tr& tr, std::string filename_prefix)
{
  std::vector<Qualities> volume_cells_quality;
  
  // global data
  double v_max = 0, v_min = std::numeric_limits<double>::infinity(), v_sum = 0;
  int v_num = 0;

  // surfacebased data
  std::vector<double> v_max_l, v_min_l, v_sum_l;
  std::vector<int> v_num_l;


  for(typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
      cit != tr.finite_cells_end();
      ++cit)
    if(cit->is_in_domain())
    {
  
      // analyse cells' quality
      const double quality = 
        CGAL::to_double(Pierre::radius_ratio(tr.tetrahedron(cit)));
      // radius ratio is in common namespace, in Slivers_exuder.h
      int index = cit->volume_index();
      if(index < 0)
        index = 0;

      const unsigned int positive_index = index;

      if( volume_cells_quality.size() <= positive_index )
        volume_cells_quality.resize(positive_index+1);
      volume_cells_quality[positive_index].push_back(quality);

      // analyse cells' volume globaly
      double v = tr.tetrahedron(cit).volume();
      if (v_max < v)
        v_max = v;
      
      if (v_min > v)
        v_min = v;
      
      v_sum += v;

      v_num++;
      
      
      // analyse cells' volume localy

      const int idx = cit->volume_index();
      
      if(idx >=0)
      {
        const unsigned int pidx = idx;

        // resize vectors
        if( v_max_l.size() <= pidx )
        {
          v_max_l.resize(pidx+1, 0);
          v_min_l.resize(pidx+1, std::numeric_limits<double>::infinity());
          v_sum_l.resize(pidx+1, 0);
          v_num_l.resize(pidx+1, 0);
        }


        // update statistics
        if( v > v_max_l[pidx] )
          v_max_l[pidx] = v_max;
          
        if( v < v_min_l[pidx] )
          v_min_l[pidx] = v_min;
          
        v_sum_l[pidx] += v;

        v_num_l[pidx]++;
      }

    }
  
    // global volume output
    std::cout << std::setprecision(3)
              << "\n\n  Tetrahedras volume:\n" 
              << "\n    min: "<< v_min 
              << "\n    avg: "<< v_sum / v_num  
              << "\n    max: "<< v_max << "\n"
              << std::endl;
  
    // local volume output
    for(unsigned int i = 0; i < v_max_l.size(); i++)
    {
      std::cout << std::setprecision(3)
                << "    volume of cell in volume #" << i << ":    "
                << std::setw(10) << v_max_l[i]  << ",  "
                << std::setw(10) << v_min_l[i]  << ",  "
                << std::setw(10) << v_sum_l[i]/v_num_l[i]                
                << std::endl;
    }

  return process_cells(volume_cells_quality, filename_prefix);
}
