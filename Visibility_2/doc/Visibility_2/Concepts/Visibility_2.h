
/*!
\ingroup PkgVisibility_2Concepts
\cgalConcept

A model of the concept `Visibility_2` offers visibility queries within 
an Arrangement. 

\cgalHasModel `CGAL::Naive_visibility_2<Arrangement_2, Traits>`
\cgalHasModel `CGAL::Preprocessed_visibility_2<Arrangement_2, Traits>`

\sa `CGAL::Naive_visibility_2<Arrangement_2, Traits>` 
\sa `CGAL::Preprocessed_visibility_2<Arrangement_2, Traits>`

*/
class Visibility_2 {
public:

/// \name Types 
/// @{
  
 /*! 
   The supported Arrangement type. 
 */ 
  typedef Hidden_type Arrangement_2; 

 /*! 
   The supported Point_2 type which is used for queries . 
 */ 
  typedef Hidden_type Point_2; 

  /*! 
    Tag identifying whether `Visibility_2` computes regularized visbility area. 
  */
  typedef Hidden_type Regularization_tag; 

  /*!
    Traits showing what kind of arrangement extension the output will have
  */
  typedef Hidden_type Traits;
  
/// @}

/// \name Constructors 
/// @{

/*! 
Constructs a `Visibility_2` object from a given `Arrangement_2` and a given Regularization tag.
*/ 
Visibility_2(const Arrangement_2& arr, Regularization_tag Rt); 

/// @}


/// \name functions 
/// @{



/*!
Return whether the object is attachted to an arrangement.
*/
  bool is_attached ();

/*!
Attaches visibility object to the given arrangement arr.
*/
  void attach ( Arrangement_2 arr);

  
/*!
Detaches the object from the arrangement it is currently attached to.
*/
  void detach ();

/*!
Access to the attached Arrangement_2.
*/
  Arrangement_2 arr();

/*! 
Computes the visibility region for the given query point q. 
\pre face is a face of  this->arr()
\pre p is in the interior of face 

*/ 
  Arrangement_2 visibility_region(const Point_2& q, const Face& face); 

/*! 
Computes for the given query point q the visibility region that is on the side of the given halfedge.   
\pre half_edge is a half edge of  this->arr()
\pre p is on halfedge  

*/ 
  Arrangement_2 visibility_region(const Point_2& q, const Halfedge& halfedge); 

/// @}


}; /* end Visibility_2 */

