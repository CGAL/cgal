namespace CGAL {
/*!
\ingroup PkgVisibility_2Classes

\brief This class is a model of the concept `Visibility_2` can answer visibility queries within
a simple polygon with no holes.

\details This class implements the algorithm of B.Joe and R.B.Simpson \cite bjrb-clvpa-87 to
obtain the visibility region, based on a scan of the boundary of the polygon and the notion 
of angular displacement as a control variable. The algorithm is a modification and extension 
of the  linear time algorithm of Lee \cite dtl-voasp-83. It computes the visibility region from a 
viewpoint that is in the interior or on the boundary of the polygon. 

The algorithm uses a stack to manipulate the vertices, and ultimately yields the visibility
region. For each scanned edge, at most 2 points are pushed onto the stack. Overall, it
will have at most 2\f$ n \f$ points pushed and popped, thus the time and space complexities of the
algorithm are \f$ O(n) \f$ even in case of degeneracies such as needles, where n is the number of 
the vertices of the polygon.

The class offers the option to either compute the visibility region or the visibility polygon, which can be chosen 
at compile time via the second template argument RegularizationTag. The default for the RegularizationTag
is ::Tag_false, which implies that the visibility region will be computed. Setting the template argument
to ::Tag_true will produce the output as a visibility polygon.

\tparam Arrangement_2 is the type of input polygonal environment and output visibility polygon.

\tparam RegularizationTag indicates whether the output should be regularized. It can be
specified by one of the following: ::Tag_true or ::Tag_false, where ::Tag_false is the default value.

\cgalModels `Visibility_2` 

\sa `CGAL::Rotational_sweep_visibility_2<Arrangement_2, RegularizationTag>`
\sa `CGAL::Preprocessed_rotational_sweep_visibility_2<Arrangement_2, RegularizationTag>`
\sa `CGAL::Triangular_expansion_visibility_2<Arrangement_2, RegularizationTag>`
*/
template <typename Arrangement_2, typename RegularizationTag = Tag_false>
class Simple_polygon_visibility_2 {
public:

/// \name Types 
/// @{

 /*!
   The arrangement type is used for input.
 */
  typedef Arrangement_2 Input_arrangement_2;

 /*!
  *The arrangement type is used for output.
  */
  typedef Arrangement_2 Output_arrangement_2;

 /*! 
   The 2D point type used for the queries. 
 */ 
  typedef Input_arrangement_2::Point_2 Point_2;

  /*!
   Face_const_handle type of input arrangement.
   */
  typedef Input_arrangement_2::Face_const_handle Face_const_handle;

  /*!
   Halfedge_const_handle type of input arrangement.
   */
  typedef Input_arrangement_2::Halfedge_const_handle Halfedge_const_handle;
   
/// @}



/// \name Tags 
/// @{
  /*! 
    Tag identifying whether the regularized visibility area is computed. 
  */
  typedef RegularizationTag Regularization_tag;
  
  /*! 
    Tag identifying that the class does not support general polygons (i.e.\ with holes). 
  */
  typedef ::Tag_false Supports_general_polygon_tag; 

  /*! 
    Tag identifying that the class supports general simple polygons. 
  */
  typedef ::Tag_true Supports_simple_polygon_tag; 
/// @}


/// \name Constructors 
/// @{

/*!
Default constructor creates an empty 'Simple_polygon_visibility_2' object that is not
attached to any arrangement yet.
*/
Simple_polygon_visibility_2();

/*! 
Constructs a `Simple_polygon_visibility_2` object that is attached to `arr`.
*/ 
Simple_polygon_visibility_2(const Input_arrangement_2& arr); 

/// @}

/// \name functions 
/// @{

/*!
Returns whether an arrangement is attached to the visibility object
*/
  bool is_attached ();

/*!
Attaches the given arrangement to the visibility object.
In case the object is already attached to another arrangement, 
the visibility object gets detached before being attached to 'arr'.
*/
  void attach (const Input_arrangement_2& arr);

/*!
Detaches the arrangement from the visibility object it is currently attached to
*/
  void detach();

/*!
Access to the attached arrangement
*/
  const Input_arrangement_2& arr();

/*! 
Computes the visibility region for the given query point `q` in the
face `f` of the arrangement that is attached to the visibility object. 
The visibility region of `q` will be stored in `out_arr`.
\param out_arr is the output arrangement 
\param q is the query point
\param f is the face of the arrangement in which the visibility region is computed
\pre `f` is a face of  `this->arr()` and represents a valid polygon. 
\pre `q` is in the interior of the given face `f`
\return a handle to the face in `out_arr` that represents the visibility region
*/ 
  Face_handle visibility_region(const Point_2& q, const Face_const_handle f, Output_arrangement_2& out_arr);


/*!
Computes the visibility region for the given query point `q` that is on `e`.If `q` is an interior point of `e`, the computed visibility region is restricted to the halfplane indicated by `e`. If `q` is an endpoint of `e`, the visibility region is restricted by `e` and its next.
The visibility region of `q` will be stored in `out_arr`.
\param q is the query point
\param e the halfedge on which `q` is located
\param out_arr is the output arrangement
\pre `e` is a halfedge of  `this->arr()`
\pre `q` is on `e`
\pre `q` equals to `e->target()->point()` if `q` is an endpoint of `e`
\return a handle to the face in `out_arr` that represents the visibility region
*/
  Face_handle visibility_region(const Point_2& q, const Halfedge_const_handle e, Output_arrangement_2& out_arr);


/// @}

}; /* end Visibility_2 */
}
