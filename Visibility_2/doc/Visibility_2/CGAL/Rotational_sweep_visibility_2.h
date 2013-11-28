namespace CGAL {
/*!
\ingroup PkgVisibility_2Classes

\brief This class is a model of the concept `Visibility_2` can answer visibility queries within a polygon that may have holes.


\details The algorithm does not require preprocessing. It relies on the algorithm of T. Asano \cite ta-aeafvpprh-85 based on angular plane sweep, with a time complexity of \f$O (n \log n)\f$ in the number of vertices.

\tparam Arrangement_2 is the type of input polygonal environment and output visibility polygon.
Arrangement_2::Traits_2 must be an instance of Arr_segment_traits_2.

\tparam RegularizationTag indicates whether the output should be regularized. It can be
specified by one of the following: ::Tag_true or ::Tag_false, where ::Tag_false is the default value.


\cgalModels `Visibility_2` 

\sa `CGAL::Simple_polygon_visibility_2<Arrangement_2, RegularizationTag>`
\sa `CGAL::Triangular_expansion_visibility_2<Arrangement_2, RegularizationTag>`

*/
template <typename Arrangement_2, typename RegularizationTag = Tag_false>
class Rotational_sweep_visibility_2 {
public:

/// \name Types 
/// @{
   
 /*!
  The type of the input arrangement.
  */
  typedef Arrangement_2  Input_arrangement_2;

   /*!
    The type of the output arrangement.
    */
  typedef Arrangement_2 Output_arrangement_2;

/// @}

/// \name Tags 
/// @{
  /*! 
    Tag identifying whether the regularized visibility area is computed. 
  */
  typedef RegularizationTag Regularization_tag;
  
  /*! 
    The class supports general polygons (i.e.\ with holes). 
  */
  typedef ::Tag_true Supports_general_polygon_tag; 

  /*! 
    The class supports general simple polygons. 
  */
  typedef ::Tag_true Supports_simple_polygon_tag; 
/// @}


/// \name Functions 
/// @{


/*!
Attaches the given arrangement to the visibility object. 

This operation takes \f$O(1)\f$ as the class does no pre-processing. 

In case the object is already attached to another arrangement, 
the visibility object gets detached before being attached to `arr`.
*/
  void attach(const Input_arrangement_2& arr);
  

/// @}

}; /* end Visibility_2 */
}
