// TODO: Copyright info

#ifndef CGAL_POINTMATCHER_REGISTER_POINT_SETS_H
#define CGAL_POINTMATCHER_REGISTER_POINT_SETS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/pointmatcher/compute_registration_transformation.h>

#include <boost/type_traits/is_same.hpp>

#include <Eigen/Dense>

namespace CGAL {

namespace pointmatcher {

// TODO: Document
template <class PointRange1, class PointRange2,
          class NamedParameters1, class NamedParameters2>
void // TODO: Can we return some sort of score?
register_point_sets (const PointRange1& point_set_1, PointRange2& point_set_2,
                     const NamedParameters1& np1, const NamedParameters2& np2)
{
  using boost::choose_param;
  using boost::get_param;

  namespace PSP = CGAL::Point_set_processing_3;
  typedef typename PSP::GetK<PointRange1, NamedParameters1>::Kernel Kernel;

  // compute registration transformation
  typename Kernel::Aff_transformation_3 res =
    compute_registration_transformation(point_set_1, point_set_2, np1, np2);

  // property map type of point_set_2 
  typedef typename PSP::GetPointMap<PointRange2, NamedParameters2>::type PointMap2;
  PointMap2 point_map2 = choose_param(get_param(np2, internal_np::point_map), PointMap2());

  // update CGAL points
  for (typename PointRange2::iterator it=point_set_2.begin(),
                                      end=point_set_2.end(); it!=end; ++it)
  {
    put(point_map2, *it, get(point_map2, *it).transform(res));
  }
}

// convenience overloads
template <class PointRange1, class PointRange2,
          class NamedParameters1>
void // TODO: Can we return some sort of score?
register_point_sets(const PointRange1& point_set_1, PointRange2& point_set_2,
                    const NamedParameters1& np1)
{
  namespace params = CGAL::Point_set_processing_3::parameters;
  return register_point_sets(point_set_1, point_set_2, np1, params::all_default(point_set_1));
}

template <class PointRange1, class PointRange2>
void // TODO: Can we return some sort of score?
register_point_sets(const PointRange1& point_set_1, PointRange2& point_set_2)
{
  namespace params = CGAL::Point_set_processing_3::parameters;
  return register_point_sets(point_set_1, point_set_2,
                             params::all_default(point_set_1),
                             params::all_default(point_set_2));
}

} } // end of namespace CGAL::pointmatcher

#endif // CGAL_POINTMATCHER_REGISTER_POINT_SETS_H