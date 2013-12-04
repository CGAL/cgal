namespace CGAL {
/*!
\ingroup PkgVisibility_2Classes

\brief This class is a model of the concept `Visibility_2` can answer visibility queries within
a simple polygon with no holes.

\details This class implements the algorithm of B.Joe and R.B.Simpson \cite bjrb-clvpa-87. 
The algorithm is a modification and extension of the  linear time algorithm of Lee \cite dtl-voasp-83.
It computes the visibility region from a viewpoint that is in the interior or on the boundary of the polygon. 

While scanning the boundary the algorithm uses a stack to manipulate the vertices, and ultimately 
yields the visibility region. For each scanned edge, at most 2 points are pushed onto the stack. 
Overall, at most 2\f$ n \f$ points are pushed or popped. Thus, the time and space complexities of the
algorithm are \f$ O(n) \f$ even in case of degeneracies such as needles, where \f$ n \f$ 
is the number of the vertices of the polygon.

\tparam Arrangement_2_ is the type used to represent the input environment.
It must be an instance of CGAL::Arrangement_2, where its CGAL::Arrangement_2::Traits_2 must be an instance of 
CGAL::Arr_segment_traits_2. 

\tparam RegularizationTag indicates whether the output should be regularized. It can be
specified by one of the following: ::Tag_true or ::Tag_false, where ::Tag_false is the default value.


\cgalModels `Visibility_2` 

\sa `CGAL::Rotational_sweep_visibility_2`
\sa `CGAL::Triangular_expansion_visibility_2`
*/
template <typename Arrangement_2_, typename RegularizationTag = Tag_true>
class Simple_polygon_visibility_2 {
public:

/// \name Types 
/// @{

 /*!
   The arrangement type is used for input.
 */
  typedef Arrangement_2 Arrangement_2;
 
/// @}

/// \name Tags 
/// @{
  /*! 
    Tag identifying whether the regularized visibility area is computed. 
  */
  typedef RegularizationTag Regularization_tag;
  
  /*! 
    The class does not support general polygons (i.e.\ with holes). 
  */
  typedef ::Tag_false Supports_general_polygon_tag; 

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
  void attach(const Arrangement_2& arr);



/// @}

}; /* end Visibility_2 */
}
