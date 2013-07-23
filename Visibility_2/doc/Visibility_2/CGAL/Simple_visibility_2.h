namespace CGAL {
/*!
\ingroup PkgVisibility_2Classes

This class is a model of the concept `Visibility_2` offers visibility queries within 
a simple polygon with no holes. Uses the algorithm of B.Joe and R.B.Simpson \cite bjrb-clvpa-87 to 
obtain the visibility polygon, based on a scan of the boundary of the polygon and the notion 
of angular displacement as a control variable. The algorithm is a modification and extension 
of the  linear time algorithm of Lee \cite TODO. It computes the visibility polygon from a 
viewpoint that is in the interior or on the boundary of the polygon. 

The algorithm uses a stack to manipulate the vertices, and ultimately yields the visibility
polygon. For each scanned edge, at most 2 points are pushed on the stack. Overall, it
will have at most 2n points pushed and popped, thus the time and space complexities of the 
algorithm are O(n) even in case of degeneracies such as needles, where n is the number of 
the vertices of the polygon.

The class also supports a regularization tag which allows the user to obtain an output
with all the isolated vertices and edges that have the same face on both sides removed.

\cgalModels `Visibility_2` 

\sa 'CGAL::Visibility_2'
\sa `CGAL::Naive_visibility_2<Arrangement_2, Regularization_tag>`
\sa `CGAL::Preprocessed_visibility_2<Arrangement_2, Regularization_tag>`

*/
template <typename Arrangement_2, typename Regularization_tag>
class Simple_visibility_2 {
public:

/// \name Types 
/// @{

  typedef Arrangement_2 Arr_traits_2;
 /*!
   The Arrangement type is used for input.
 */
  typedef Arr_traits_2::Input_Arrangement_2 Input_Arrangement_2;

 /*!
  *The Arrangement type is used for output.
  */
  typedef Arr_traits_2::Output_Arrangement_2 Output_Arrangement_2;

 /*! 
   The Point_2 type which is used for queries. 
 */ 
  typedef Input_Arrangement_2::Point_2 Point_2;
  
  
/// @}

/// \name Constructors 
/// @{

/*!
Default constructor creates an empty `Visibility_2` object
*/
Simple_visibility_2();

/*! 
Constructs a `Visibility_2` object from a given `Input_Arrangement_2`
*/ 
Simple_visibility_2(const Input_Arrangement_2& arr); 

/// @}


/// \name functions 
/// @{

/*!
Returns whether an arrangement is attached to the visibility object
*/
  bool is_attached ();

/*!
Attaches the given arrangement to the visibility object
*/
  void attach (const Input_Arrangement_2 &arr);

/*!
Detaches the arrangement from the visibility object it is currently attached to
*/
  void detach();

/*!
Access to the attached arrangement
*/
  Input_Arrangement_2 arr();

/*!
Computes the visibility region for the given query point q. 
\pre face is a face of this->arr() with no holes
\pre p is in the interior of face
\pre out_arr is the output arrangement 
\return a face handle to the bounded face of the output
*/ 
  Face_handle visibility_region(const Point_2& q, const Face& face, Output_Arrangement_2& out_arr); 

/*! 
Computes for the given query point q the visibility region that is on the side of the given halfedge.   
\pre half_edge is a half edge of this->arr()
\pre p is located on the halfedge  
\pre out_arr is the output arrangement 
\return a face handle to the bounded face of the output
*/ 
  Face_handle visibility_region(const Point_2& q, const Halfedge& halfedge, Output_Arrangement_2& out_arr); 

/// @}


}; /* end Visibility_2 */
}
