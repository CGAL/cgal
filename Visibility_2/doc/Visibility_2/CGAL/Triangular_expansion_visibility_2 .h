namespace CGAL {
/*!
\ingroup PkgVisibility_2Classes

This class is a model of the concept `Visibility_2` offering visibility queries within 
an Arrangement. The algorithm converts the input into restricted Delaunay triangulation and computes visibility by expanding the triangle that contains the query point. Preprocessing takes \f$ O(n^2)\f$ time and \f$ O(n) \f$ space. The worst query time is \f$O(n^2) \f$ but it works significantly fast for the case with a small number of holes.

Arrangement_2 gives information about the input/output Arrangements and the extension between them.
The Regularization_tag indicates whether the output Arrangement result should be regularized. It can be
specified by one of the following: CGAL::Tag_true or CGAL::Tag_false.

\cgalModels `Visibility_2` 

\sa `CGAL::Simple_visibility_2<Arrangement_2, Regularization_tag>`
\sa `CGAL::Rotational_sweep_visibility_2<Arrangement_2, Regularization_tag>`
\sa `CGAL::Preprocessed_Rotational_sweep_visibility_2<Arrangement_2, Regularization_tag>`


*/
template <typename Arrangement_2, typename Regularization_tag>
class Triangular_expansion_visibility_2 {
public:

/// \name Types 
/// @{

 /*!
  The type of input arrangement.
  */
   typedef Arrangement_2 Input_Arrangement_2;

 /*!
  The type of output arrangement.
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
    Tag identifying that the class supports general polygons (i.e. with holes). 
  */
  typedef Tag_true Supports_general_polygon_tag; 

  /*! 
    Tag identifying that the class supports general simple polygons. 
  */
  typedef Tag_true Supports_simple_polygon_tag; 
/// @}



/// \name Constructors 
/// @{

/*!
Default constructor creates an empty 'Preprocessed_visibility_2' object, that is not
attached to any arrangement yet.
*/
Triangular_expansion_visibility_2();

/*! 
Constructs a `Preprocessed_visibility_2` object from a given `Input_Arrangement_2` and attaches it to `arr` and does preprocessing.
*/ 
Triangular_expansion_visibility_2(const Input_Arrangement_2& arr);

/// @}

/// \name functions 
/// @{

/*!
Returns whether an arrangement is attached to the visibility object
*/
  bool is_attached();

/*!
Attaches the given arrangement to the visibility object and does preprocessing.
*/
  void attach(const Input_Arrangement_2 &arr);

/*!
Detaches the arrangement from the visibility object it is currently attached to
*/
  void detach ();

/*!
Access to the attached arrangement
*/
  const Input_Arrangement_2& arr();

/*! 
Computes the visibility region for the given query point `q` in the
face `f` of the arrangement that is attached to the visibility object. 
The visibility region of `q` will be stored in `out_arr`.
\param q is the query point from which the visibility region is computed
\param f is the face of the arrangement in which the visibility region is computed
\param out_arr is the output arrangement 
\pre `f` is a face of  `this->arr()`
\pre `q` is in the interior or on the boundary of the given face `f`
\return the face handle to the face in `out_arr` that represents the visibility region
*/ 
  Face_handle visibility_region(const Point_2& q, const Face& f, Output_Arrangement_2& out_arr); 

/*! 
Computes for the given query point `q` the visibility region that is on the side of `halfedge`. 
The visibility region of `q` will be stored in `out_arr`.
\param q is the query point from which the visibility region is computed
\param halfedge the halfedge on which 'q' is located
\param out_arr is the output arrangement  
\pre `half_edge` is a half edge of  `this->arr()`
\pre `q` is on halfedge
\return the face handle to the face in `out_arr` that represents the visibility region
*/ 
  Face_handle visibility_region(const Point_2& q, const Halfedge& halfedge, Output_Arrangement_2& out_arr); 

/// @}

}; /* end Visibility_2 */
}  /* namespace CGAL */
