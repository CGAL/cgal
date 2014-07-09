namespace CGAL {
/*!
\ingroup PkgVisibility_2Classes

\brief This class is a model of the concept `Visibility_2` can answer visibility queries within a polygon that may have holes.

\details The algorithm obtains a constrained triangulation from input arrangement, then computes visibility by 
expanding the triangle that contains the query point. 
Preprocessing takes \f$ O(n)\f$ time and \f$ O(n) \f$ space, where \f$ n \f$ is the number of vertices of input polygon. 
The query time is \f$ O(nh)\f$, where \f$ h \f$ is the number of holes+1 of input polygon. Thus, for simple polygons 
the algorithm is even linear but it can also be  \f$ O(n^2)\f$ in the worst case as the number of holes can be linear in \f$ n \f$. 


\tparam Arrangement_2_ is the type used to represent the input environment.
It must be an instance of CGAL::Arrangement_2, where its CGAL::Arrangement_2::Traits_2 must be an instance of 
CGAL::Arr_segment_traits_2. 

\tparam RegularizationTag indicates whether the output should be regularized. It can be
specified by one of the following: ::Tag_true or ::Tag_false, where ::Tag_false is the default value.

\cgalModels `Visibility_2` 

\sa `CGAL::Simple_polygon_visibility_2`
\sa `CGAL::Rotational_sweep_visibility_2`


*/
template <typename Arrangement_2_, typename RegularizationTag = Tag_true>
class Triangular_expansion_visibility_2 {
public:

/// \name Types 
/// @{

 /*!
  The type of the input arrangement.
  */
   typedef Arrangement_2_ Arrangement_2;

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


/// \name Functions 
/// @{

/*!
Attaches the given arrangement to the visibility object and computes the restricted triangulation. 
This takes \f$ O(n) \f$ time, where \f$ n \f$ is the number of vertices. 

From this moment on the class observes changes in the arrangement. If the arrangement changes 
the a new restricted triangulation is computed right before a new query. It is also possible
to force a re-computation by re-attaching the current arrangement. 

In case the object is already attached to another arrangement, 
the visibility object gets detached before being attached to `arr`.
*/
void attach(const Arrangement_2& arr);


/// @}

}; /* end Visibility_2 */
}  /* namespace CGAL */
