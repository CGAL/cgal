namespace CGAL {
/*!
\ingroup PkgVisibility_2Classes

This class is a model of the concept `Visibility_2` offers visibility queries within 
a simple polygon with no holes. Uses a linear algorithm to obtain visibility.

\sa `Visibility_2` 

*/
template <typename Arrangement_2, typename Traits>
class Simple_visibility_2 {
public:

/// \name Types 
/// @{


 /*! 
   The Point_2 type which is used for queries. 
 */ 
  typedef Arrangement_2::Point_2 Point_2; 

  /*! 
    Tag identifying whether `Visibility_2` computes the regularized visibility area. 
  */
  typedef bool Regularization_tag; 

  /*!
    The Arrangement type which is used for output.
  */
  typedef notknown Output_Arrangement_2;
  
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


}
