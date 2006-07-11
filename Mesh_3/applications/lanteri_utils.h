From Marc.Scherfenberg@sophia.inria.fr Tue Jul 11 17:57:16 2006
X-Spam-Checker-Version: SpamAssassin 3.1.0 (2005-09-13) on reale.ens.fr
X-Spam-Level: 
X-Spam-Status: No, score=0.0 required=5.0 tests=none autolearn=disabled 
	version=3.1.0
Received: from nef2.ens.fr (nef2 [129.199.96.40]) by clipper.ens.fr (8.13.1/jb-1.1)
	id k6BFxRBv023558 for <rineau@clipper.ens.fr>; Tue, 11 Jul 2006 17:59:27 +0200 (MEST)
Return-Path: <Marc.Scherfenberg@sophia.inria.fr>
Received: from sophia.inria.fr (sophia.inria.fr [138.96.64.20])
          by nef2.ens.fr (8.13.6/1.01.28121999) with ESMTP id k6BFxRud064892
          for <rineau@clipper.ens.fr>; Tue, 11 Jul 2006 17:59:27 +0200 (CEST)
X-Envelope-To: <rineau@clipper.ens.fr>
Received: from localhost (localhost [127.0.0.1])
	by sophia.inria.fr (8.13.7/8.13.4) with ESMTP id k6BFxRCd007910
	for <rineau@clipper.ens.fr>; Tue, 11 Jul 2006 17:59:27 +0200
Received: from [138.96.211.57] (historis.inria.fr [138.96.211.57])
	by sophia.inria.fr (8.13.7/8.13.4) with ESMTP id k6BFvGJ5007139
	for <rineau@clipper.ens.fr>; Tue, 11 Jul 2006 17:57:16 +0200
Message-ID: <44B3CA5C.9070307@sophia.inria.fr>
Date: Tue, 11 Jul 2006 17:57:16 +0200
From: Marc Scherfenberg <Marc.Scherfenberg@sophia.inria.fr>
User-Agent: Thunderbird 1.5.0.4 (X11/20060614)
MIME-Version: 1.0
To: Laurent Rineau <rineau@clipper.ens.fr>
Subject: angle
Content-Type: multipart/mixed;
  boundary="------------080605060100030707060101"
X-Greylist: Recipient e-mail whitelisted, not delayed by milter-greylist-1.5.10 (nef2.ens.fr [129.199.96.32]); Tue, 11 Jul 2006 17:59:27 +0200 (CEST)
X-Greylist: Sender IP whitelisted, not delayed by milter-greylist-2.0.2 (sophia.inria.fr [138.96.64.20]); Tue, 11 Jul 2006 17:57:16 +0200 (MEST)
X-Virus-Scanned: by amavisd-milter (http://amavis.org/)
X-Virus-Scanned: by amavisd-new at sophia.inria.fr
X-UID: 29230
X-Length: 14242
Status: R
X-Status: NT
X-KMail-EncryptionState:  
X-KMail-SignatureState:  
X-KMail-MDN-Sent:  

This is a multi-part message in MIME format.
--------------080605060100030707060101
Content-Type: text/plain; charset=ISO-8859-1; format=flowed
Content-Transfer-Encoding: 7bit

thanks for reviewing

--------------080605060100030707060101
Content-Type: text/x-chdr;
 name="lanteri_utils_angle.h"
Content-Transfer-Encoding: 7bit
Content-Disposition: inline;
 filename="lanteri_utils_angle.h"

#include <string>
#include <sstream>
#include <iomanip> // std::setprecision, std::setw
#include <algorithm> // std::maximum, std::minimum
#include <numeric> // std::accumulate
#include <cmath> // arcsin

#include <boost/tuple/tuple.hpp>


template<typename Triangulation>
class Compute_min_angle
{  
  public:
    
    typedef typename Triangulation::Cell_handle Cell_handle;
    typedef typename Triangulation::Tetrahedron Tetrahedron;  
    typedef typename Triangulation::Point Point;
  
    // constructor
    Compute_min_angle(Triangulation _tr) : tr(_tr) {}  
    
    // computes the minimum angle between all 6 faces of a tetrahedra
    double 
    operator()(const Cell_handle cell) const
    {
      const Tetrahedron tet = tr.tetrahedron(cell);
      
      double min_quotient = compute_quotient(cell, tet, 0, 1, 2, 3);
      min_quotient = std::min(min_quotient, compute_quotient(cell, tet, 0, 2, 1, 3));
      min_quotient = std::min(min_quotient, compute_quotient(cell, tet, 0, 3, 1, 2));
      min_quotient = std::min(min_quotient, compute_quotient(cell, tet, 1, 2, 0, 3));  
      min_quotient = std::min(min_quotient, compute_quotient(cell, tet, 1, 3, 0, 2));  
      min_quotient = std::min(min_quotient, compute_quotient(cell, tet, 2, 3, 0, 1));  
      
      const double V = CGAL::to_double(tet.volume());
      
      return asin( 1.5 * V * min_quotient) * 180 / 3.14159265; 
    }      
  
  
  private:  
    
    Triangulation tr;    
    
    double compute_quotient(const Cell_handle cell, const Tetrahedron tet, const int i, const int j, const int k, const int l) const
    {
      
      const Point pi = tet.vertex(i);
      const Point pj = tet.vertex(j);
            
      const double le = CGAL::to_double(CGAL::sqrt(CGAL::squared_distance(pi, pj)));
        
      const double Ak = CGAL::to_double(CGAL::sqrt(tr.triangle(cell, k).squared_area()));   
      
      const double Al = CGAL::to_double(CGAL::sqrt(tr.triangle(cell, l).squared_area()));    
      
      return le / Ak / Al;
    }   
};



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



void output_legend(std::ostream* out_stream = &std::cout, std::string prefix = "")
{
  *out_stream << std::endl;

  *out_stream << std::setprecision(3)
      << prefix
      << std::setw(42) << "min"
      << std::setw(13) << "avg"
      << std::setw(13) << "max"
      << std::endl;    
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

  output_stream << std::setw(10) << minimum << ",  ";

  if( size == 0 )
    output_stream << std::setw(10) << "nan" << ",  ";
  else
    output_stream << std::setw(10) << sum / size << ",  ";

  output_stream << std::setw(10)  << maximum; 
  
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
  
   
  // local edge-lengths output (surface)  
  output_legend(out_stream, prefix);

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

  // local edge-lengths output (volume)
  output_legend(out_stream, prefix); 
  
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

  *out_stream << std::endl;
  
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
  std::vector<Qualities> cells_min_angle;
  
  const Compute_min_angle<Tr> compute_min_angle(tr); 
  
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
        cells_min_angle.resize(positive_index+1);    
      }

      cells_quality[positive_index].push_back(quality);
      cells_volume[positive_index].push_back(tr.tetrahedron(cit).volume());
      cells_min_angle[positive_index].push_back(compute_min_angle(cit));   
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

              << std::accumulate(sum.begin(), sum.end(), 0.) / 
                       std::accumulate(++size.begin(), size.end(), 0)
               // the division may be "not a number"

              << "\n"
              << prefix << "max tetrahedron volume: " 
              << *(std::max_element(++maximum.begin(), maximum.end())) << "\n"
              << std::endl;
  
  // local volume output
  output_legend(out_stream, prefix);  
  
  for(unsigned int i = 0; i < vectors_size; i++)
  {
    *out_stream << std::setprecision(3)
                << prefix << "volume of cells in volume #" << i << ":   "
                << format_max_min_sum_size(maximum[i], minimum[i], sum[i], size[i])
                << std::endl;
  }


  // local quality output    
  output_legend(out_stream, prefix); 
  
  for(unsigned int i = 0; i < vectors_size; i++)
  {
    *out_stream << std::setprecision(3)
                << prefix << "quality of cells in volume #" << i << ":  "
                << output_max_min_average(cells_quality[i].begin(),
                                          cells_quality[i].end())
                << std::endl;
  }

  
  // local angle output 
  output_legend(out_stream, prefix);
  
  for(unsigned int i = 0; i < vectors_size; i++)
  {
    *out_stream << std::setprecision(3)
        << prefix << "angles of cells in volume #" << i << ":   "
        << output_max_min_average(cells_min_angle[i].begin(),
                                  cells_min_angle[i].end())
        << std::endl;
  } 
    
  *out_stream << std::endl; 
  
  return process_cells(cells_quality, filename_prefix);
}
   



--------------080605060100030707060101--

