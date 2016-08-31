#ifndef CGAL_DATA_CLASSIFICATION_NEIGHBORHOOD_H
#define CGAL_DATA_CLASSIFICATION_NEIGHBORHOOD_H

#include <vector>

#include <boost/iterator/counting_iterator.hpp>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Default_diagonalize_traits.h>
#include <CGAL/centroid.h>

#include <CGAL/Data_classification/Image.h>

namespace CGAL {

namespace Data_classification {

  /*!
    \ingroup PkgDataClassification

    \brief 

    \tparam Kernel The geometric kernel used.

  */
template <typename Kernel, typename RandomAccessIterator, typename PointPMap>
class Neighborhood
{
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point;
  
  class My_point_property_map{
    RandomAccessIterator begin;
    PointPMap point_pmap;
    
  public:
    typedef Point value_type;
    typedef const value_type& reference;
    typedef std::size_t key_type;
    typedef boost::lvalue_property_map_tag category;
    My_point_property_map () { }
    My_point_property_map (RandomAccessIterator begin, PointPMap point_pmap)
      : begin (begin), point_pmap (point_pmap) { }
    reference operator[] (key_type k) const { return get(point_pmap, begin[k]); }
    friend inline reference get (const My_point_property_map& ppmap, key_type i) 
    { return ppmap[i]; }
  };

  typedef Search_traits_3<Kernel> SearchTraits_3;
  typedef Search_traits_adapter <std::size_t, My_point_property_map, SearchTraits_3> Search_traits;
  typedef Sliding_midpoint<Search_traits> Splitter;
  typedef Distance_adapter<std::size_t, My_point_property_map, Euclidean_distance<SearchTraits_3> > Distance;
  typedef Kd_tree<Search_traits, Splitter, Tag_true> Tree;
  typedef Fuzzy_sphere<Search_traits> Sphere;
  typedef Orthogonal_k_neighbor_search<Search_traits, Distance, Splitter, Tree> Knn;


  Tree* m_tree;
  Distance m_distance;
  RandomAccessIterator m_begin;
  RandomAccessIterator m_end;
  PointPMap m_point_pmap;

  std::vector<std::vector<std::size_t> > m_precomputed_neighbors;
  
public:

  Neighborhood () : m_tree (NULL) { }
  
  Neighborhood (const RandomAccessIterator& begin,
                const RandomAccessIterator& end,
                PointPMap point_pmap)
    : m_tree (NULL), m_begin (begin), m_end (end),
      m_point_pmap (point_pmap)
  {
    std::size_t size = end - begin;
    
    My_point_property_map pmap (begin, point_pmap);
    m_tree = new Tree (boost::counting_iterator<std::size_t> (0),
                       boost::counting_iterator<std::size_t> (size),
                       Splitter(),
                       Search_traits (pmap));
    m_distance = Distance (pmap);
  }

  ~Neighborhood ()
  {
    if (m_tree != NULL)
      delete m_tree;
  }

  template <typename OutputIterator>
  void range_neighbors (std::size_t index, const FT radius_neighbors, OutputIterator output) const
  {
    CGAL_assertion (m_tree != NULL);
    Sphere fs (get(m_point_pmap, m_begin[index]), radius_neighbors, 0, m_tree->traits());
    m_tree->search (output, fs);
  }

  template <typename OutputIterator>
  void k_neighbors (std::size_t index, const std::size_t k, OutputIterator output) const
  {
    CGAL_assertion (m_tree != NULL);
    Knn search (*m_tree, get(m_point_pmap, m_begin[index]), k, 0, true, m_distance);
    for (typename Knn::iterator it = search.begin(); it != search.end(); ++ it)
      *(output ++) = it->first;
  }
};
  

}
  
}


#endif // CGAL_DATA_CLASSIFICATION_NEIGHBORHOOD_H
