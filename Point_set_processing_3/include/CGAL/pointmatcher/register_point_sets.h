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

namespace internal {

template <class Kernel,
          class PointRange1,
          class PointRange2,
          class PointMap1,
          class PointMap2,
          class VectorMap1,
          class VectorMap2>
void // TODO: Can we return some sort of score?
register_point_sets(const PointRange1& range1, PointRange2& range2,
                    PointMap1 point_map1,   PointMap2 point_map2,
                    VectorMap1 vector_map1, VectorMap2 vector_map2,
                    ICP<typename Kernel::FT> icp)
{
  typename Kernel::Aff_transformation_3 res =
    compute_registration_transformation<Kernel>(range1, range2,
                                                point_map1, point_map2,
                                                vector_map1, vector_map2,
                                                icp);

  // update CGAL points
  for (typename PointRange2::iterator it=range2.begin(),
                                      end=range2.end(); it!=end; ++it)
  {
    put(point_map2, *it, get(point_map2, *it).transform(res));
  }
}

} // end of namespace internal

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

  // property map types
  typedef typename PSP::GetPointMap<PointRange1, NamedParameters1>::type PointMap1;
  typedef typename PSP::GetPointMap<PointRange2, NamedParameters2>::type PointMap2;
  CGAL_static_assertion_msg((boost::is_same< typename boost::property_traits<PointMap1>::value_type,
                                             typename boost::property_traits<PointMap2>::value_type> ::value),
                            "The point type of input ranges must be the same");

  typedef typename PSP::GetNormalMap<PointRange1, NamedParameters1>::type NormalMap1;
  typedef typename PSP::GetNormalMap<PointRange2, NamedParameters2>::type NormalMap2;
  CGAL_static_assertion_msg((boost::is_same< typename boost::property_traits<NormalMap1>::value_type,
                                             typename boost::property_traits<NormalMap2>::value_type> ::value),
                            "The vector type of input ranges must be the same");

  typedef typename PSP::GetK<PointRange1, NamedParameters1>::Kernel Kernel;
  typedef typename Kernel::FT Scalar;

  PointMap1 point_map1 = choose_param(get_param(np1, internal_np::point_map), PointMap1());
  NormalMap1 normal_map1 = choose_param(get_param(np1, internal_np::normal_map), NormalMap1());
  PointMap2 point_map2 = choose_param(get_param(np2, internal_np::point_map), PointMap2());
  NormalMap2 normal_map2 = choose_param(get_param(np2, internal_np::normal_map), NormalMap2());

  internal::register_point_sets<Kernel>(point_set_1, point_set_2,
                                        point_map1, point_map2,
                                        normal_map1, normal_map2,
                                        internal::construct_icp<Scalar>(np1, np2));
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