namespace CGAL {

/*!
 * \ingroup PkgConvexHull3Traits
 * 
 * The class `Extreme_points_traits_adapter_3` serves as a traits class for the function 
 * `extreme_points_3()`. It permits the use of this function for accessing vertices, indices,
 * or anything that can be used as `key_type` for `PointPropertyMap`.
 * 
 * \tparam Base_traits a model of `ConvexHullTraits_3` and `IsStronglyConvexTraits_3`
 * \tparam PointPropertyMap a model of `ReadablePropertyMap` with `CGAL::Point_3` 
 * as value_type.
 * 
 * \cgalModels `ConvexHullTraits_3`
 * \cgalModels `IsStronglyConvexTraits_3` 
 */
template<class Base_traits,class PointPropertyMap>
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
}//end CGAL
