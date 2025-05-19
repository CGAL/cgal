namespace CGAL {

/*!
 * \ingroup PkgConvexHull3Traits
 *
 * The class `Extreme_points_traits_adapter_3` serves as a traits class for the function
 * `extreme_points_3()`. It permits the use of this function for accessing vertices, indices,
 * or anything that can be used as `key_type` for `PointPropertyMap`.
 *
 * \tparam PointPropertyMap a model of `ReadablePropertyMap` with `CGAL::Point_3`
 * as value_type.
 * \tparam Base_traits a model of `ConvexHullTraits_3` and `IsStronglyConvexTraits_3`.
 * If the kernel `R` of the points from `PointPropertyMap`
 * is a kernel with exact predicates but inexact constructions
 * (in practice we check `R::Has_filtered_predicates_tag` is `Tag_true` and `R::FT` is a floating point type),
 * then the default traits class used is `Convex_hull_traits_3<R>`, and `R` otherwise.
 *
 * \cgalModels{ConvexHullTraits_3,IsStronglyConvexTraits_3}
 */
template<class PointPropertyMap, class Base_traits=Default>
class Extreme_points_traits_adapter_3
{
public:
  /*!
   * Constructor for the adapter.
   *
   * It uses the functors from `Base_traits` with a call to `get()` for each argument.
   *
   * \param pmap the propertymap used for retrieving the data.
   * \param base the base ConvexHullTraits_3 used for `extreme_points_3()`.
   */
  Extreme_points_traits_adapter_3(const PointPropertyMap& pmap, Base_traits base=Base_traits());
  /*!

   */
  typedef boost::property_traits<PointPropertyMap>::key_type Point_3;

  /*!
   * \return the `CGAL::Point_3` associated with `k`
   */
  boost::property_traits<PointPropertyMap>::reference
  get_point(const typename boost::property_traits<PointPropertyMap>::key_type& k) const;

};

/*!
 * \ingroup PkgConvexHull3Functions
 * Returns `Extreme_points_traits_adapter_3<PointPropertyMap, Base_traits>(pmap, traits)`.
 */
template<class PointPropertyMap,class Base_traits>
Extreme_points_traits_adapter_3<PointPropertyMap, Base_traits>
make_extreme_points_traits_adapter(const PointPropertyMap& pmap, Base_traits traits);
}//end CGAL
