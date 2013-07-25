namespace CGAL {
/*!
\ingroup PkgVisibility_2Classes

This class is a model of the concept `Visibility_2` offers visibility queries within 
a simple polygon with no holes. This class implements the algorithm of B.Joe and R.B.Simpson \cite bjrb-clvpa-87 to 
obtain the visibility region, based on a scan of the boundary of the polygon and the notion 
of angular displacement as a control variable. The algorithm is a modification and extension 
of the  linear time algorithm of Lee \cite dtl-voasp-83. It computes the visibility region from a 
viewpoint that is in the interior or on the boundary of the polygon. 

The algorithm uses a stack to manipulate the vertices, and ultimately yields the visibility
region. For each scanned edge, at most 2 points are pushed on the stack. Overall, it
will have at most 2n points pushed and popped, thus the time and space complexities of the 
algorithm are O(n) even in case of degeneracies such as needles, where n is the number of 
the vertices of the polygon.

The class offers to either compute the visibility region or the visibility polygon, which can be chosen 
at compile time via the second template argument Regularization_tag. The default for the Regularization_tag
is CGAL::Tag_false, which means the visibility region will be computed. Setting the template argument
to CGAL::Tag_true will produce the output as a visibility polygon.

\cgalModels `Visibility_2` 

\sa `CGAL::Naive_visibility_2<Arrangement_2, Regularization_tag>`
\sa `CGAL::Preprocessed_visibility_2<Arrangement_2, Regularization_tag>`

*/
template <typename Arrangement_2, typename Regularization_tag>
class Simple_visibility_2 {
public:

/// \name Types 
/// @{

 /*!
   The Arrangement type is used for input.
 */
  typedef Arrangement_2 Input_Arrangement_2;

 /*!
  *The Arrangement type is used for output.
  */
  typedef Arrangement_2 Output_Arrangement_2;

 /*! 
   The Point_2 type which is used for queries. 
 */ 
  typedef Input_Arrangement_2::Point_2 Point_2;

  /*!
   Face_handle type of input Arrangement.
   */
  typedef Input_Arrangement_2::Face_handle Face_handle;

  /*!
   Halfedge_handle type of input Arrangement.
   */
  typedef Input_Arrangement_2::Halfedge_handle Halfedge_handle;
   
/// @}

/// \name Constructors 
/// @{

/*!
Default constructor creates an empty 'Simple_visibility_2' object, that is not
attached to any  arrangement yet.
*/
Simple_visibility_2();

/*! 
Constructs a `Simple_visibility_2` object from a given `Input_Arrangement_2` and attaches it to `arr`.
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
  const Input_Arrangement_2& arr();

/*! 
Computes the visibility region for the given query point `q` in the
face `f` of the arrangement that is attached to the visibility object. 
The visibility region of `q` will be stored in `out_arr`.
\param out_arr is the output arrangement 
\pre `f` is a face of  this->arr()
\pre q is in the interior or on the foundary of the given face `f`
\return the face handle to the face in `out_arr` that represents the visibility region
*/ 
  Face_handle visibility_region(const Point_2& q, const Face& face, Output_Arrangement_2& out_arr); 

/*! 
Computes for the given query point `q` the visibility region that is on the side of `halfedge`. 
The visibility region of `q` will be stored in `out_arr`.
\param out_arr is the output arrangement  
\pre half_edge is a half edge of  this->arr()
\pre q is on halfedge
\return the face handle to the face in `out_arr` that represents the visibility region
*/ 
  Face_handle visibility_region(const Point_2& q, const Halfedge& halfedge, Output_Arrangement_2& out_arr); 

/// @}

}; /* end Visibility_2 */
}
