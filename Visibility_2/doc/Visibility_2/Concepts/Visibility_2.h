
/*!
\ingroup PkgVisibility_2Concepts
\cgalConcept

A model of the concept `Visibility_2` offers visibility queries within 
an Arrangement. 

\cgalHasModel `CGAL::Simple_visibility_2<Arrangement_2, Regularization_tag>`
\cgalHasModel `CGAL::Rotational_sweep_visibility_2<Arrangement_2, Regularization_tag>`
\cgalHasModel `CGAL::Preprocessed_visibility_2<Arrangement_2, Regularization_tag>`

\sa `CGAL::Simple_visibility_2<Arrangement_2, Regularization_tag>`
\sa `CGAL::Rotational_sweep_visibility_2<Arrangement_2, Regularization_tag>`
\sa `CGAL::Preprocessed_visibility_2<Arrangement_2, Regularization_tag>`

*/
class Visibility_2 {
public:

/// \name Types 
/// @{
  
 /*! 
   The supported arrangement type of input. 
 */ 
  typedef Hidden_type Input_Arrangement_2; 

  /*!
    The supported arrangement type of output.
   */
  typedef Hidden_type Output_Arrangement_2;

 /*! 
   The supported Point_2 type of the input type, used for queries. 
 */ 
  typedef Input_Arrangement_2::Point_2 Point_2;

  /*!
   * The supported Face handle type of the input Arrangement.
   */
  typedef Input_Arrangement_2::Face_handle Face_handle;

  /*!
   * The supported Halfedge handle type of the input Arrangement.
   */
  typedef Input_Arrangement_2::Halfedge_handle Halfedge_handle;

/// @}

/// \name Tags 
/// @{
  /*! 
    Tag identifying whether the regularized visibility area is computed (either `Tag_true` or `Tag_false`). 
  */
  typedef Hidden_type Regularization_tag;
  
  /*! 
    Tag identifying whether general polygons, that is, with holes are supported (either `Tag_true` or `Tag_false`). 
  */
  typedef Hidden_type Supports_general_polygon_tag; 

  /*! 
    Tag identifying whether general simple polygons are supported (either `Tag_true` or `Tag_false`). 
  */
  typedef Hidden_type Supports_simple_polygon_tag; 
/// @}

/// \name Constructors 
/// @{

/*! 
Default constructor creates an empty `Visibility_2` object, that is not
attached to any arrangement yet.
*/ 
Visibility_2(); 

/*! 
Constructs a `Visibility_2` object from a given `Arrangement_2` and attaches it to `arr`.
*/ 
Visibility_2(const Input_Arrangement_2 &arr);

/// @}


/// \name Functions 
/// @{

/*!
Returns whether an arrangement is attached to the visibility object
*/
  bool is_attached ();

/*!
Attaches the given arrangement `arr` to the visibility object.
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
The visibility region of `q` will be saved to `out_arr`.
\param q is the query point from which the visibility region is computed
\param f is the face of the arrangement in which the visibility region is computed
\param out_arr is the output arrangement 
\pre `f` is a face of  `this->arr()`
\pre `q` is in the interior or on the boundary of the given face `f`
\return the face handle to the face in `out_arr` that represents the visibility region
*/ 
  Face_handle visibility_region(const Point_2& q, const Face_handle& f, Output_Arrangement_2& out_arr);

/*! 
Computes for the given query point `q` the visibility region that is on the side of `halfedge`. 
The visibility region of `q` will be saved to `out_arr`.
\param q is the query point from which the visibility region is computed
\param halfedge the halfedge on which 'q' is located
\param out_arr is the output arrangement  
\pre `half_edge` is a half edge of  `this->arr()`
\pre `q` is on halfedge
\return the face handle to the face in `out_arr` that represents the visibility region
*/ 
  Face_handle visibility_region(const Point_2& q, const Halfedge_handle& halfedge, Output_Arrangement_2& out_arr);

/// @}

}; /* end Visibility_2 */

