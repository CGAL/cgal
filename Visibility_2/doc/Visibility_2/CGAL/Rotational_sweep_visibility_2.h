namespace CGAL {
/*!
\ingroup PkgVisibility_2Classes

\brief This class is a model of the concept `Visibility_2` can answer visibility queries within a polygon that may have holes.


\details The algorithm does not require preprocessing. It relies on the algorithm of T. Asano \cite ta-aeafvpprh-85 based on angular plane sweep, with a time complexity of \f$O (n \log n)\f$ in the number of vertices.


\tparam Arrangement_2_ is the type used to represent the input environment.
It must be an instance of `CGAL::Arrangement_2`, where its `CGAL::Arrangement_2::Traits_2` must be an instance of 
`CGAL::Arr_segment_traits_2`, or of `CGAL::Arr_non_caching_segment_traits_2`. 

\tparam RegularizationCategory indicates whether the output should be regularized. It can be
specified by one of the following: `#Tag_true` or `#Tag_false`, where `#Tag_false` is the default value.



\cgalModels `Visibility_2` 

\sa `CGAL::Simple_polygon_visibility_2`
\sa `CGAL::Triangular_expansion_visibility_2`

*/
template <typename Arrangement_2_, typename RegularizationCategory = Tag_true>
class Rotational_sweep_visibility_2 {
public:

/// \name Types 
/// @{
   
 /*!
  The type of the input arrangement.
  */
  typedef Arrangement_2  Arrangement_2;

/// @}

/// \name Tags 
/// @{
  /*!
  \copydoc ::Visibility_2::Regularization_category
  */
  typedef RegularizationCategory Regularization_category;
  
  /*! 
    See \ref Visibility_2::Supports_general_polygon_category.
  */
  typedef Tag_true Supports_general_polygon_category; 

  /*! 
    See \ref Visibility_2::Supports_simple_polygon_category.
  */
  typedef Tag_true Supports_simple_polygon_category; 
/// @}


/// \name Functions 
/// @{


/*!
Attaches the given arrangement to the visibility object. 

This operation takes \f$O(1)\f$ as the class does no pre-processing. 

In case the object is already attached to another arrangement, 
the visibility object gets detached before being attached to `arr`.
*/
  void attach(const Arrangement_2& arr);
  

/// @}

}; /* end Visibility_2 */
}
