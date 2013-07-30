namespace CGAL {
/*!
\ingroup PkgVisibility_2Classes

\brief This class is a model of the concept `Visibility_2` offering visibility queries within
a simple polygon with no holes.

\details This class implements the algorithm of B.Joe and R.B.Simpson \cite bjrb-clvpa-87 to
obtain the visibility region, based on a scan of the boundary of the polygon and the notion 
of angular displacement as a control variable. The algorithm is a modification and extension 
of the  linear time algorithm of Lee \cite dtl-voasp-83. It computes the visibility region from a 
viewpoint that is in the interior or on the boundary of the polygon. 

The algorithm uses a stack to manipulate the vertices, and ultimately yields the visibility
region. For each scanned edge, at most 2 points are pushed on the stack. Overall, it
will have at most 2n points pushed and popped, thus the time and space complexities of the 
algorithm are \f$ O(n) \f$ even in case of degeneracies such as needles, where n is the number of 
the vertices of the polygon.

The class offers the option to either compute the visibility region or the visibility polygon, which can be chosen 
at compile time via the second template argument Regularization_tag. The default for the Regularization_tag
is ::Tag_false, which means the visibility region will be computed. Setting the template argument
to ::Tag_true will produce the output as a visibility polygon.

\tparam Arrangement_2 is the type of input polygonal environment and output visibility polygon.

\tparam Regularization_tag indicates whether the output should be regularized. It can be
specified by one of the following: ::Tag_true or ::Tag_false.

\cgalModels `Visibility_2` 

\sa `CGAL::Rotational_sweep_visibility_2<Arrangement_2, Regularization_tag>`
\sa `CGAL::Preprocessed_rotational_sweep_visibility_2<Arrangement_2, Regularization_tag>`
\sa `CGAL::Triangular_expansion_visibility_2<Arrangement_2, Regularization_tag>`
*/
template <typename Arrangement_2, typename Regularization_tag>
class Simple_visibility_2 {
public:

/// \name Types 
/// @{

 /*!
   The arrangement type is used for input.
 */
  typedef Arrangement_2 Input_Arrangement_2;

 /*!
  *The arrangement type is used for output.
  */
  typedef Arrangement_2 Output_Arrangement_2;

 /*! 
   The Point_2 type which is used for queries. 
 */ 
  typedef Input_Arrangement_2::Point_2 Point_2;

  /*!
   Face_handle type of input arrangement.
   */
  typedef Input_Arrangement_2::Face_handle Face_handle;

  /*!
   Halfedge_handle type of input arrangement.
   */
  typedef Input_Arrangement_2::Halfedge_handle Halfedge_handle;
   
/// @}



/// \name Tags 
/// @{
  /*! 
    Tag identifying whether the regularized visibility area is computed. 
  */
  typedef Regularization_tag Regularization_tag;
  
  /*! 
    Tag identifying that the class does not support general polygons (i.e. with holes). 
  */
  typedef Tag_false Supports_general_polygon_tag; 

  /*! 
    Tag identifying that the class supports general simple polygons. 
  */
  typedef ::Tag_true Supports_simple_polygon_tag; 
/// @}


/// \name Constructors 
/// @{

/*!
Default constructor creates an empty 'Simple_visibility_2' object, that is not
attached to any arrangement yet.
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
\param q is the query point from which the visibility region is computed
\param f is the face of the arrangement in which the visibility region is computed
\pre `f` is a face of  `this->arr()`
\pre `q` is in the interior or on the boundary of the given face `f`
\return the face handle to the face in `out_arr` that represents the visibility region
*/ 
  Face_handle visibility_region(const Point_2& q, const Face& f, Output_Arrangement_2& out_arr); 

/*! 
Computes the visibility region for the given query point `q` that is on the side of `halfedge`.
The visibility region of `q` will be stored in `out_arr`.
\param q is the query point from which the visibility region is computed
\param halfedge the halfedge on which `q` is located
\param out_arr is the output arrangement  
\pre `half_edge` is a half edge of  `this->arr()`
\pre `q` is on halfedge
\return the face handle to the face in `out_arr` that represents the visibility region
*/ 
  Face_handle visibility_region(const Point_2& q, const Halfedge& halfedge, Output_Arrangement_2& out_arr); 

/// @}

}; /* end Visibility_2 */
}
