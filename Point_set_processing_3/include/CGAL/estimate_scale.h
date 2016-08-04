#ifndef CGAL_ESTIMATE_SCALE_H
#define CGAL_ESTIMATE_SCALE_H

#include <CGAL/Search_traits_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/assertions.h>
#include <CGAL/hierarchy_simplify_point_set.h>

#include <iterator>
#include <list>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>  
#endif // CGAL_LINKED_WITH_TBB

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL
namespace internal {

template <class Kernel>
class Quick_multiscale_approximate_knn_distance
{
  typedef typename Kernel::FT FT;
  typedef Search_traits_3<Kernel> Tree_traits;
  typedef Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Iterator;

  std::size_t m_cluster_size;
  std::vector<Tree*> m_trees;

public:

  template <typename InputIterator, typename PointPMap>
  Quick_multiscale_approximate_knn_distance (InputIterator first,
                                             InputIterator beyond,
                                             PointPMap point_pmap,
                                             std::size_t cluster_size = 10)
    : m_cluster_size (cluster_size)
  {
    m_trees.push_back (new Tree ());
    std::size_t nb_pts = 0;
    for (InputIterator it = first; it != beyond; ++ it)
      {
        m_trees[0]->insert (get(point_pmap, *it));
        ++ nb_pts;
      }
    
    std::size_t nb_trees = 0;
    while (nb_pts > m_cluster_size)
      {
        nb_trees ++;
        nb_pts /= m_cluster_size;
      }

    m_trees.reserve (nb_trees);

    InputIterator first_unused = beyond;

    for (std::size_t i = 1; i < nb_trees; ++ i)
      {
        m_trees.push_back (new Tree());
        first_unused
          = CGAL::hierarchy_simplify_point_set (first, first_unused, point_pmap, cluster_size, 1./3.);

        for (InputIterator it = first; it != first_unused; ++ it)
          m_trees.back()->insert (get(point_pmap, *it));
      }
  }

  ~Quick_multiscale_approximate_knn_distance()
  {
    for (std::size_t i = 0; i < m_trees.size(); ++ i)
      delete m_trees[i];
  }

  template <typename InputIterator, typename PointPMap>
  std::size_t compute_scale (InputIterator query, PointPMap point_pmap)
  {
    std::size_t out = 0;

    std::size_t weight = 1;
    FT dist_min = std::numeric_limits<FT>::max();
    FT sum_sq_distances = 0.;
    std::size_t nb = 0;
    
    for (std::size_t t = 0; t < m_trees.size(); ++ t)
      {
        Neighbor_search search (*(m_trees[t]), get(point_pmap, *query), m_cluster_size);
        Iterator it = search.begin();
        
        if (t != 0) // Skip first point except on first scale
          ++ it;

        for (; it != search.end(); ++ it)
          {
            sum_sq_distances += weight * it->second;
            nb += weight;
            
            if (nb < 6) // do not consider values under 6
              continue;

            FT dist = std::sqrt (sum_sq_distances / nb)
              / std::pow (nb, 0.375); // nb^(5/12)

            if (dist < dist_min)
              {
                dist_min = dist;
                out = nb;
              }
          }
        weight *= m_cluster_size;
      }
    return out;
  }

};

} /* namespace internal */
/// \endcond



// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

template <typename SamplesInputIterator,
          typename SamplesPointPMap,
          typename QueriesInputIterator,
          typename QueriesPointPMap,
          typename OutputIterator,
          typename Kernel
>
void
estimate_local_k_neighbor_scales(
  SamplesInputIterator first,
  SamplesInputIterator beyond,
  SamplesPointPMap samples_pmap,
  QueriesInputIterator first_query,
  QueriesInputIterator beyond_query,
  QueriesPointPMap queries_pmap,
  OutputIterator output,
  const Kernel& /*kernel*/) ///< geometric traits.
{
  // Build multi-scale KD-tree
  internal::Quick_multiscale_approximate_knn_distance<Kernel> kdtree (first, beyond, samples_pmap);

  // Compute local scales everywhere
  for (QueriesInputIterator it = first_query; it != beyond_query; ++ it)
    *(output ++) = kdtree.compute_scale (it, queries_pmap);

}

template <typename InputIterator,
          typename PointPMap,
          typename Kernel
>
std::size_t
estimate_global_k_neighbor_scale(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map: value_type of InputIterator -> Point_3
  const Kernel& kernel) ///< geometric traits.
{
  std::vector<std::size_t> scales;
  estimate_local_k_neighbor_scales (first, beyond, point_pmap,
                                    first, beyond, point_pmap,
                                    std::back_inserter (scales),
                                    kernel);
  std::sort (scales.begin(), scales.end());
  return scales[scales.size() / 2];
}

template <typename InputIterator,
          typename PointPMap,
          typename Kernel
>
typename Kernel::FT
estimate_global_range_scale(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map: value_type of InputIterator -> Point_3
  const Kernel& /*kernel*/) ///< geometric traits.
{

  return 0.;
}



} //namespace CGAL

#endif // CGAL_ESTIMATE_SCALE_3_H
