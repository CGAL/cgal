namespace CGAL {
/*!
\ingroup PkgVisibility_2Classes

\brief This class is a model of the concept `Visibility_2` can answer visibility queries within a polygon that may have holes.

\details The algorithm obtains a constrained triangulation from input arrangement, then computes visibility by 
expanding the triangle that contains the query point. 
Preprocessing takes \f$ O(n)\f$ time and \f$ O(n) \f$ space, where \f$ n \f$ is the number of vertices of input polygon. 
The query time is \f$ O(nh)\f$, where \f$ h \f$ is the number of holes+1 of input polygon. Thus, for simple polygons 
the algorithm is even linear but it can also be  \f$ O(n^2)\f$ in the worst case as the number of holes can be linear in \f$ n \f$. 

\tparam Arrangement_2 is the type of input polygonal environment and output visibility polygon.

\tparam RegularizationTag indicates whether the output should be regularized. It can be
specified by one of the following: ::Tag_true or ::Tag_false, where ::Tag_false is the default value.

\cgalModels `Visibility_2` 

\sa `CGAL::Simple_polygon_visibility_2<Arrangement_2, RegularizationTag>`
\sa `CGAL::Rotational_sweep_visibility_2<Arrangement_2, RegularizationTag>`


*/
template <typename Arrangement_2, typename RegularizationTag = Tag_false>
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
   The 2D point type used for the queries. 
 */ 
  typedef Input_arrangement_2::Point_2 Point_2; 

  /*!
   Face_const_handle type of the input arrangement.
   */
  typedef Input_arrangement_2::Face_const_handle Face_const_handle;

  /*!
   Halfedge_const_handle type of the input arrangement.
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
    The class supports general polygons (i.e.\ with holes). 
  */
  typedef ::Tag_true Supports_general_polygon_tag; 

  /*! 
    The class supports general simple polygons. 
  */
  typedef ::Tag_true Supports_simple_polygon_tag; 
/// @}



/// \name Constructors 
/// @{

/*!
Default constructor creates an empty `Triangular_expansion_visibility_2` object that is not
attached to any arrangement yet.
*/
Triangular_expansion_visibility_2();

/*! 
Constructs a `Triangular_expansion_visibility_2` object that is attached to `arr`.
*/ 
Triangular_expansion_visibility_2(const Input_arrangement_2& arr);

/// @}

/// \name functions 
/// @{

/*!
Returns whether an arrangement is attached to the visibility object
*/
  bool is_attached() const;

/*!
Attaches the given arrangement to the visibility object and computes the restricted triangulation. 
This takes \f$ O(n) \f$ time, where \f$ n \f$ of vertices. Modifying the attached arrangement 
also changes the stored restricted triangulation which in the worst may again take \f$ O(n) \f$ time.

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
Computes the visibility region of `q` in the
face `f` of the arrangement that is attached to the visibility object. 
The visibility region of `q` will be stored in `out_arr`.
\param q is the query point
\param f is the face of the arrangement in which the visibility region is computed
\param out_arr is the output arrangement 
\pre `f` is a face of `arr()` and represents a valid polygon. 
\pre `q` is in the interior of the given face `f`
\return a handle to the face in `out_arr` that represents the visibility region
*/ 
  typename Output_arrangement_2::Face_handle compute_visibility(const Point_2& q, const Face_const_handle f, Output_arrangement_2& out_arr) const;


/*!
Computes the visibility region of `q` that is on `e`. If `q` is an interior point of `e`, the computed visibility region is restricted to the halfplane indicated by `e`. If `q` is an endpoint of `e`, the visibility region is restricted by `e` and its next.
The visibility region of `q` will be stored in `out_arr`.
\param q is the query point
\param e the halfedge on which `q` is located
\param out_arr is the output arrangement
\pre `e` is a halfedge of `arr()`
\pre `q` is on `e`
\pre `q` equals to `e->target()->point()` if `q` is an endpoint of `e`
\return a handle to the face in `out_arr` that represents the visibility region
*/
  typename Output_arrangement_2::Face_handle compute_visibility(const Point_2& q, const Halfedge_const_handle e, Output_arrangement_2& out_arr) const;

/// @}

}; /* end Visibility_2 */
}  /* namespace CGAL */
