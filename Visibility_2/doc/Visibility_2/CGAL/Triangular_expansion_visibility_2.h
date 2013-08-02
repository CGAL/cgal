namespace CGAL {
/*!
\ingroup PkgVisibility_2Classes

\brief This class is a model of the concept `Visibility_2` offering visibility queries within a polygon that may have holes.

\details The algorithm obtains a constrained Delaunay triangulation from input arrangement then computes visibility by 
expanding the triangle that contains the query point. Preprocessing takes \f$ O(nh)\f$ time and \f$ O(n) \f$ space, 
where \f$ n \f$ and \f$ h \f$ are the numbers  of vertices and holes of input polygon respectively. The worst query 
time is \f$O(n^2) \f$, but it is expected to work significantly faster for cases with a small number of holes.

\tparam Arrangement_2 is the type of input polygonal environment and output visibility polygon.

\tparam RegularizationTag indicates whether the output should be regularized. It can be
specified by one of the following: ::Tag_true or ::Tag_false, which is the default value.

\cgalModels `Visibility_2` 

\sa `CGAL::Simple_polygon_visibility_2<Arrangement_2, RegularizationTag>`
\sa `CGAL::Rotational_sweep_visibility_2<Arrangement_2, RegularizationTag>`
\sa `CGAL::Preprocessed_rotational_sweep_visibility_2<Arrangement_2, RegularizationTag>`


*/
template <typename Arrangement_2, typename RegularizationTag>
class Triangular_expansion_visibility_2 {
public:

/// \name Types 
/// @{

 /*!
  The type of the input arrangement.
  */
   typedef Arrangement_2 Input_arrangement_2;

 /*!
  The type of the output arrangement.
  */
   typedef Arrangement_2 Output_arrangement_2;

 /*! 
   The Point_2 type, which is used for queries.
 */ 
  typedef Input_arrangement_2::Point_2 Point_2; 

  /*!
   Face_handle type of the input arrangement.
   */
  typedef Input_arrangement_2::Face_handle Face_handle;

  /*!
   Halfedge_handle type of the input arrangement.
   */
  typedef Input_arrangement_2::Halfedge_handle Halfedge_handle;

/// @}


/// \name Tags 
/// @{
  /*! 
    Tag identifying whether the regularized visibility area is computed. 
  */
  typedef RegularizationTag Regularization_tag;
  
  /*! 
    Tag identifying that the class supports general polygons (i.e.\ with holes). 
  */
  typedef ::Tag_true Supports_general_polygon_tag; 

  /*! 
    Tag identifying that the class supports general simple polygons. 
  */
  typedef ::Tag_true Supports_simple_polygon_tag; 
/// @}



/// \name Constructors 
/// @{

/*!
Default constructor creates an empty 'Triangular_expansion_visibility_2' object, that is not
attached to any arrangement yet.
*/
Triangular_expansion_visibility_2();

/*! 
Constructs a `Triangular_expansion_visibility_2` object from a given `Input_arrangement_2` instance and attaches it to `arr` and does preprocessing.
*/ 
Triangular_expansion_visibility_2(const Input_arrangement_2& arr);

/// @}

/// \name functions 
/// @{

/*!
Returns whether an arrangement is attached to the visibility object
*/
  bool is_attached();

/*!
Attaches the given arrangement to the visibility object and does preprocessing.
In case the object is already attached to another arrangement, 
the visibility object gets detached before being attached to 'arr'.
*/
  void attach(const Input_arrangement_2 &arr);

/*!
Detaches the arrangement from the visibility object it is currently attached to
*/
  void detach ();

/*!
Access to the attached arrangement
*/
  const Input_arrangement_2& arr();

/*! 
Computes the visibility region for the given query point `q` in the
face `f` of the arrangement that is attached to the visibility object. 
The visibility region of `q` will be stored in `out_arr`.
\param q is the query point from which the visibility region is computed
\param f is the face of the arrangement in which the visibility region is computed
\param out_arr is the output arrangement 
\pre `f` is a face of  `this->arr()`, defined as a regular polygon 
\pre `q` is in the interior or on the boundary of the given face `f`
\return a handle to the face in `out_arr` that represents the visibility region
*/ 
  Face_handle visibility_region(const Point_2& q, const Face_handle& f, Output_arrangement_2& out_arr);


/*!
Computes the visibility region for the given query point `q` that is on `e`.If `q` is an interior point of `e`, the computed visibility region is restricted to the halfplane indicated by `e`. If `q` is an endpoint of `e`, the visibility region is restricted by `e` and its next.
The visibility region of `q` will be stored in `out_arr`.
\param q is the query point from which the visibility region is computed
\param e the halfedge on which `q` is located
\param out_arr is the output arrangement
\pre `e` is a halfedge of  `this->arr()`
\pre `q` is on `e`
\pre `q` equals to `e->target()->point()` if `q` is an endpoint of `e`
\return a handle to the face in `out_arr` that represents the visibility region
*/
  Face_handle visibility_region(const Point_2& q, const Halfedge_handle& e, Output_arrangement_2& out_arr);

/// @}

}; /* end Visibility_2 */
}  /* namespace CGAL */
