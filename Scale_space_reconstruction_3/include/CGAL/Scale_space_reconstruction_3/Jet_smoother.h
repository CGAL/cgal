#ifndef CGAL_SCALE_SPACE_RECONSTRUCTION_3_JET_SMOOTHER_H
#define CGAL_SCALE_SPACE_RECONSTRUCTION_3_JET_SMOOTHER_H

#include <CGAL/jet_smooth_point_set.h>

#ifdef CGAL_LINKED_WITH_TBB
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#endif // CGAL_LINKED_WITH_TBB

namespace CGAL
{

namespace Scale_space_reconstruction_3
{
  
template <typename Geom_traits,
#ifdef CGAL_LINKED_WITH_TBB
          typename ConcurrencyTag = CGAL::Parallel_tag>
#else
          typename ConcurrencyTag = CGAL::Sequential_tag>
#endif
class Jet_smoother
{
public:
  typedef typename Geom_traits::FT FT; ///< defines the point type.
  typedef typename Geom_traits::Point_3 Point; ///< defines the point typ.e
private:

  unsigned int m_k;
  unsigned int m_degree_fitting;
  unsigned int m_degree_monge;

public:

  Jet_smoother (unsigned int k = 12,
                unsigned int degree_fitting = 2,
                unsigned int degree_monge = 2)
    : m_k (k), m_degree_fitting (degree_fitting), m_degree_monge (degree_monge)
  { }

  template <typename InputIterator>
  void operator() (InputIterator begin, InputIterator end)
  {
    CGAL::jet_smooth_point_set<ConcurrencyTag> (begin, end, m_k, m_degree_fitting, m_degree_monge);
  }
  
};


} // namespace Scale_space_reconstruction_3

} // namespace CGAL

#endif // CGAL_SCALE_SPACE_RECONSTRUCTION_3_JET_SMOOTHER_H
