#include <string>
#include <sstream>
#include <iomanip> // std::setprecision, std::setw
#include <algorithm> // std::maximum, std::minimum
#include <numeric> // std::accumulate

#include <boost/tuple/tuple.hpp>

template <typename Iterator> // better be RandomAccessIterator, because
                          // std::distance() is used
boost::tuple<
  typename std::iterator_traits<Iterator>::value_type,
  typename std::iterator_traits<Iterator>::value_type,
  typename std::iterator_traits<Iterator>::value_type,
  typename std::iterator_traits<Iterator>::difference_type
>
compute_max_min_sum_size(Iterator begin, Iterator end)
{
  typedef typename std::iterator_traits<Iterator>::value_type T;
  typedef typename std::iterator_traits<Iterator>::difference_type size_type;

  T maximum = std::numeric_limits<T>::infinity();
  T minimum = - maximum;

  const Iterator pos_min = std::min_element(begin, end);

  if( pos_min != end ) // non-empty range
    minimum = *pos_min;

  const Iterator pos_max = std::max_element(begin, end);

  if( pos_max != end ) // non-empty range
    maximum = *pos_max;

  const T sum = std::accumulate(begin, end, 0.);

  const size_type size = std::distance(begin, end);

	return boost::make_tuple(maximum, minimum,  sum, size);
}

template <typename Iterator> // better be RandomAccessIterator, because
                          // std::distance() is used
std::string output_max_min_average(Iterator begin, Iterator end)
{
  typedef typename std::iterator_traits<Iterator>::value_type T;
  typedef typename std::iterator_traits<Iterator>::difference_type size_type;

  T minimum;
  T maximum;
  T sum;
  size_type size;
  boost::tie(maximum, minimum, sum, size) = compute_max_min_sum_size(begin, end);

  return format_max_min_sum_size(maximum, minimum, sum, size);
}

template <typename T, typename size_type>
std::string format_max_min_sum_size(const T maximum,
                                    const T minimum,
                                    const T sum,
                                    const size_type size)
{
  std::stringstream output_stream;

  output_stream << std::setw(10) << maximum << ",  "
                << std::setw(10) << minimum << ",  ";

  if( size == 0 )
    output_stream << std::setw(10) << "nan";
  else
    output_stream << std::setw(10) << sum / size;

  return output_stream.str();
}

// analyse edges
template <class Tr>
bool
scan_edges_and_process(const Tr& tr,
                       std::vector<double> length_bounds,
                       std::string filename_prefix,
                       std::string prefix = "",
                       // prefix to each line output
                       std::ostream* out_stream = &std::cout
                       // output stream
                       )
{
  // reminder: Qualities is std::vector<double> (voir "distribution.h")
  std::vector<Qualities> surface_edges_length;
  std::vector<Qualities> volume_edges_length;

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
      }
      
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
        }

        volume_edges_length[positive_index].push_back(length);
      }
    }
  }

  *out_stream << std::setprecision(3)
              << prefix
              << std::setw(42) << "max"
              << std::setw(13) << "min"
              << std::setw(13) << "avg"
              << std::endl;

  const typename Qualities::size_type surface_vector_size = 
    surface_edges_length.size();
  for(unsigned int i = 0; i < surface_vector_size; ++i)
  {
    *out_stream << prefix
                << "length for edges on surface #" << i << ": "
                << output_max_min_average(surface_edges_length[i].begin(),
                                          surface_edges_length[i].end())
                << std::endl;
  }

  *out_stream << std::endl;

  const typename Qualities::size_type volume_vector_size = 
    volume_edges_length.size();
  for(unsigned int i = 0; i < volume_vector_size; ++i)
  {
    *out_stream << prefix
                << "length for edges in volume #" << i << ":  "
                << output_max_min_average(volume_edges_length[i].begin(),
                                          volume_edges_length[i].end())
                << std::endl;
  }

  return process_surface_edges(surface_edges_length,
                               length_bounds,
                               filename_prefix,
                               out_stream) &&
         process_volume_edges(volume_edges_length,
                              length_bounds,
                              filename_prefix,
                              out_stream);
}

// analyse cells
template <class Tr>
bool
scan_cells_and_process(const Tr& tr,
                       std::string filename_prefix,
                       std::string prefix = "",
                       // prefix to each line output
                       std::ostream* out_stream = &std::cout
                       // output stream
)
{
  std::vector<Qualities> cells_quality;
  std::vector<Qualities> cells_volume;
  
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

      if( cells_quality.size() <= positive_index )
      {
        cells_quality.resize(positive_index+1);
        cells_volume.resize(positive_index+1);
      }

      cells_quality[positive_index].push_back(quality);
      cells_volume[positive_index].push_back(tr.tetrahedron(cit).volume());
    }

  const typename Qualities::size_type vectors_size = 
    cells_quality.size();

  std::vector<double> maximum(vectors_size);
  std::vector<double> minimum(vectors_size);
  std::vector<double> sum(vectors_size);
  std::vector<unsigned int> size(vectors_size);

  for(unsigned int i = 0; i < vectors_size; ++i)
  {
    boost::tie(maximum[i], minimum[i], sum[i], size[i]) = 
      compute_max_min_sum_size(cells_volume[i].begin(),
                               cells_volume[i].end());
  }

  // global volume output
  *out_stream << std::setprecision(3)
              << prefix << "min tetrahedron volume: "
                    << *(std::min_element(++minimum.begin(), minimum.end())) << "\n"
              << prefix << "avg tetrahedron volume: "
    // ++minimum.begin(), because we do not want to take into account cells
    // of index 0, which should be an empty set.

              << std::accumulate(sum.begin(), sum.end(), 0.) / 
                       std::accumulate(++size.begin(), size.end(), 0)
               // the division may be "not a number"

              << "\n"
              << prefix << "max tetrahedron volume: " 
              << *(std::max_element(++maximum.begin(), maximum.end())) << "\n"
              << std::endl;
  
  *out_stream << std::setprecision(3)
              << prefix
              << std::setw(42) << "max"
              << std::setw(13) << "min"
              << std::setw(13) << "avg"
              << std::endl;

  // local volume output
  for(unsigned int i = 0; i < vectors_size; i++)
  {
    *out_stream << std::setprecision(3)
                << prefix << "volume of cells in volume #" << i << ":   "
                << format_max_min_sum_size(maximum[i], minimum[i], sum[i], size[i])
                << std::endl;
  }


  *out_stream << std::setprecision(3)
              << prefix
              << std::setw(42) << "max"
              << std::setw(13) << "min"
              << std::setw(13) << "avg"
              << std::endl;

  for(unsigned int i = 0; i < vectors_size; i++)
  {
    *out_stream << std::setprecision(3)
                << prefix << "quality of cells in volume #" << i << ":  "
                << output_max_min_average(cells_quality[i].begin(),
                                          cells_quality[i].end())
                << std::endl;
  }

  return process_cells(cells_quality, filename_prefix);
}
