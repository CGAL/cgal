namespace CGAL {
/*!
\ingroup PkgVisibility_2Classes

This class is a model of the concept `Visibility_2` offers visibility queries within 
a simple polygon with no holes. Uses the linear algorithm of B.Joe and R.B.Simpson [1] to obtain 
the visibility polygon.

\sa `Visibility_2` 

*/
template <typename ArrExtensionTraits_2, typename Regularization_tag>
class Simple_visibility_2 {
public:

/// \name Types 
/// @{

  typedef ArrExtensionTraits_2 Arr_extension_traits_2;
 /*!
   The Arrangement type is used for input.
 */
  typedef Arr_extension_traits_2::Input_Arrangement_2 Input_Arrangement_2;

 /*!
  *The Arrangement type is used for output.
  */
  typedef Arr_extension_traits_2::Output_Arrangement_2 Output_Arrangement_2;

 /*! 
   The Point_2 type which is used for queries. 
 */ 
  typedef Input_Arrangement_2::Point_2 Point_2;



  
/// @}

/// \name Constructors 
/// @{

/*! 
Constructs a `Visibility_2` object from a given `Input_Arrangement_2`
*/ 
Simple_visibility_2(const Input_Arrangement_2& arr); 

/// @}


/// \name functions 
/// @{

/*!

Computes the visibility region for the given query point q. 
\pre face is a face of this->arr() with no holes
\pre p is in the interior of face
\pre out_arr is the output arrangement 

*/ 
  void visibility_region(const Point_2& q, const Face& face, Output_Arrangement_2& out_arr); 

/*! 
Computes for the given query point q the visibility region that is on the side of the given halfedge.   
\pre half_edge is a half edge of  this->arr()
\pre p is on halfedge  
\pre out_arr is the output arrangement 

*/ 
  void visibility_region(const Point_2& q, const Halfedge& halfedge, Output_Arrangement_2& out_arr); 

/// @}


}; /* end Visibility_2 */

/*!
[1] B. Joe, R. B. Simpson, "Corrections to Lee's visibility polygon algorithm", BIT Numerical Mathematics 
Volume 27, Issue 4 , pp 458-473
*/
}
