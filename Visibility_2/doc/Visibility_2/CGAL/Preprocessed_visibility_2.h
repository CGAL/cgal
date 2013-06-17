namespace CGAL {
/*!
\ingroup PkgVisibility_2Classes

This class is a model of the concept `Visibility_2` offers visibility queries within 
an Arrangement. The algorithm it applies to obtain visibility is using preprocessing.

\sa `Visibility_2` 
\sa `CGAL::Naive_visibility_2<Arrangement_2, Traits>`

*/
template <typename Arrangement_2, typename Traits>
class Preprocessed_visibility_2 {
public:

/// \name Types 
/// @{


 /*! 
   The Point_2 type which is used for queries. 
 */ 
  typedef Arrangement_2::Point_2 Point_2; 

  /*! 
    Tag identifying whether `Visibility_2` computes regularized visibility area. 
  */
  typedef bool Regularization_tag; 

  /*!
    The Arrangement type which is used for output.
  */
  typedef notknown Ouput_Arrangement_2;
  
/// @}

/// \name Constructors 
/// @{

/*! 
Constructs a `Visibility_2` object from a given `Input_Arrangement_2`
*/ 
Naive_visibility_2(const Input_Arrangement_2& arr); 

/// @}


/// \name functions 
/// @{



/*!
Return whether the object is attached to an arrangement.
*/
  bool is_attached ();

/*!
Attaches visibility object to the given arrangement arr.
*/
  void attach ( Input_Arrangement_2 arr);

  
/*!
Detaches the object from the arrangement it is currently attached to.
*/
  void detach ();

/*!
Access to the attached Arrangement_2.
*/
  Input_Arrangement_2 arr();

/*! 
Computes the visibility region for the given query point q. 
\pre face is a face of  this->arr()
\pre p is in the interior of face 
\pre out_arr is the output arrangement

*/ 
  void visibility_region(const Point_2& q, const Face& face, Output_Arrangement_2& out_arr); 

/*! 
Computes for the given query point q the visibility region that is on the side of the given halfedge.   
\pre half_edge is a half edge of  this->arr()
\pre p is on halfedge
\pre out_arr is the output arrangement  

*/ 
  void visibility_region(const Point_2& q, const Halfedge& halfedge, Output_Arrangement_2& out_arr); 

/// @}


}; /* end Visibility_2 */


}  /* namespace CGAL */
