
/*!
\ingroup PkgVisibility_2Concepts
\cgalConcept

A model of the `Visibility_2` concept can be attached to an `Arrangement_2` instance to 
answer visibility queries within the faces of this arrangement.

\cgalHasModel `CGAL::Simple_polygon_visibility_2<Arrangement_2, RegularizationTag>`
\cgalHasModel `CGAL::Rotational_sweep_visibility_2<Arrangement_2, RegularizationTag>`
\cgalHasModel `CGAL::Triangular_expansion_visibility_2<Arrangement_2, RegularizationTag>`

*/
class Visibility_2 {
public:

/// \name Types 
/// @{
  
 /*! 
   The arrangement type of the input.
 */ 
  typedef Hidden_type Input_arrangement_2; 

  /*!
    The arrangement type of the output.
   */
  typedef Hidden_type Output_arrangement_2;

 /*! 
   The 2D point type used for the queries. 
 */ 
  typedef Input_arrangement_2::Point_2 Point_2;

  /*!
   * The face handle type of the input arrangement.
   */
  typedef Input_arrangement_2::Face_const_handle Face_const_handle;

  /*!
   * The halfedge handle type of the input arrangement.
   */
  typedef Input_arrangement_2::Halfedge_const_handle Halfedge_const_handle;

/// @}

/// \name Tags 
/// @{
  /*! 
    Tag identifying whether the regularized visibility area is computed (either \ref CGAL::Tag_true or \ref CGAL::Tag_false). 
  */
  typedef Hidden_type Regularization_tag;
  
  /*! 
    Tag identifying whether general polygons (with holes) are supported (either \ref CGAL::Tag_true or \ref CGAL::Tag_false). 
  */
  typedef Hidden_type Supports_general_polygon_tag; 

  /*! 
    Tag identifying whether general simple polygons are supported (either \ref CGAL::Tag_true or \ref CGAL::Tag_false). 
  */
  typedef Hidden_type Supports_simple_polygon_tag; 
/// @}

/// \name Constructors 
/// @{

/*! 
Default constructor creates an empty `Visibility_2` object that is not
attached to any arrangement yet.
*/ 
Visibility_2(); 

/*! 
Constructs a `Visibility_2` object that is attached to `arr`.
*/ 
Visibility_2(const Input_arrangement_2& arr);

/// @}


/// \name Functions 
/// @{

/*!
Returns whether an arrangement is attached to the visibility object
*/
  bool is_attached() const;

/*!
Attaches the given arrangement `arr` to the visibility object.
In case the object is already attached to another arrangement, 
the visibility object gets detached before being attached to `arr`.
*/
  void attach(const Input_arrangement_2& arr);

  
/*!
Detaches the arrangement from the visibility object it is currently attached to
*/
  void detach();

/*!
Access to the attached arrangement
*/
  const Input_arrangement_2& arr() const;

/*! 
Computes the visibility region for the given query point `q` in the
face \f$ f \f$ of the arrangement that is attached to the visibility object.
The visibility region of `q` will be saved to `out_arr`, that is, 
all features but the unbounded face of `out_arr` represent the visibility region. 
\param q is the query point
\param f is the face of the arrangement in which the visibility region is computed
\param out_arr is the output arrangement 
\pre `f` is a face of  `this->arr()`
\pre `q` is in the interior of the given face `f`
\return the face handle to the face in `out_arr` that represents iterior of the visibility region
*/ 
  typename Output_arrangement_2::Face_handle compute_visibility(const Point_2& q, const Face_const_handle f, Output_arrangement_2& out_arr) const;

/*!
Computes the visibility region in `e->face()` for the given query 
point `q` which must be located on `e`. 
If `q` is an interior point of `e`, the computed visibility 
region is restricted to the side indicated by the halfedge `e`.
If `q` is an endpoint of `e`, the visibility region is restricted by `e` and `e->next()`.
The visibility region of `q` will be stored in `out_arr`, that is, 
all features but the unbounded face of `out_arr` represent the visibility region.
\param q is the query point
\param e the halfedge on which `q` is located
\param out_arr is the output arrangement
\pre `e` is a halfedge of `this->arr()`
\pre `q` is on `e`
\pre `q` equals to `e->target()->point()` if `q` is an endpoint of `e`
\return a handle to the face in `out_arr` that represents the iterior of the visibility region
*/
  typename Output_arrangement_2::Face_handle compute_visibility(const Point_2& q, const Halfedge_const_handle e, Output_arrangement_2& out_arr) const;

/// @}

}; /* end Visibility_2 */

