
/*!
\ingroup PkgVisibility_2Concepts
\cgalConcept

A model of the concept `Visibility_2` offers visibility queries within 
an Arrangement. 

\cgalHasModel `CGAL::Naive_visibility_2<ArrExtensionTraits_2, Regularization_tag>`
\cgalHasModel `CGAL::Preprocessed_visibility_2<ArrExtensionTraits_2, Regularization_tag>`

\sa `CGAL::Naive_visibility_2<ArrExtensionTraits_2, Regularization_tag>`
\sa `CGAL::Preprocessed_visibility_2<ArrExtensionTraits_2, Regularization_tag>`
\sa `CGAL::Simple_visibility_2<ArrExtensionTraits_2, Regularization_tag>`

*/
class Visibility_2 {
public:

/// \name Types 
/// @{
  
 /*! 
   The supported Arrangement type of input. 
 */ 
  typedef Hidden_type Input_Arrangement_2; 

  /*!
    The supported Arrangement type of output.
   */
  typedef Hidden_type Output_Arrangement_2;

 /*! 
   The supported Point_2 type which is used for queries . 
 */ 
  typedef Hidden_type Point_2; 

  /*!
   * The supported Face handle type of Arrangement.
   */
  typedef Hidden_type Face_handle;

  /*!
   * The supported Halfedge handle type of Arrangement.
   */
  typedef Hidden_type Halfedge_handle;

  /*! 
    Tag identifying whether `Visibility_2` computes regularized visbility area. 
  */
  typedef Hidden_type Regularization_tag; 



  
/// @}

/// \name Constructors 
/// @{

/*! 
Constructs a `Visibility_2` object from a given `Arrangement_2` and a given Regularization tag.
*/ 
Visibility_2(const Input_Arrangement_2& arr);

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

*/ 
  void visibility_region(const Point_2& q, const Face_handle& face, Output_Arrangement_2& out_arr); 

/*! 
Computes for the given query point q the visibility region that is on the side of the given halfedge.   
\pre half_edge is a half edge of  this->arr()
\pre p is on halfedge  

*/ 
  void visibility_region(const Point_2& q, const Halfedge_handle& halfedge, Output_Arrangement_2& out_arr); 

/// @}


}; /* end Visibility_2 */

