
/*!
\ingroup PkgVisibility_2Concepts
\cgalConcept

A model of the concept `Visibility_2` offers visibility queries within 
an Arrangement. 

\cgalHasModel `CGAL::Simple_visibility_2<Arrangement_2, Regularization_tag>`
\cgalHasModel `CGAL::Naive_visibility_2<Arrangement_2, Regularization_tag>`
\cgalHasModel `CGAL::Preprocessed_visibility_2<Arrangement_2, Regularization_tag>`

\sa `CGAL::Simple_visibility_2<Arrangement_2, Regularization_tag>`
\sa `CGAL::Naive_visibility_2<Arrangement_2, Regularization_tag>`
\sa `CGAL::Preprocessed_visibility_2<Arrangement_2, Regularization_tag>`

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
  /*!
   * The supported Face_handle type of output Arrangement.
   */
  typedef Output_Arrangement_2::Face_handle Out_face_handle;
	


  /*! 
    Tag identifying whether `Visibility_2` computes regularized visbility area. 
  */
  typedef Hidden_type Regularization_tag; 



  
/// @}

/// \name Constructors 
/// @{

/*! 
default constructor. 
*/ 
Visibility_2(); 

/*! 
Constructs a `Visibility_2` object from a given `Arrangement_2` and a given Regularization tag.
*/ 
Visibility_2(const Input_Arrangement_2 &arr);

/// @}


/// \name functions 
/// @{

/*!
Returns whether an arrangement is attached to the visibility object
*/
  bool is_attached ();

/*!
Attaches the given arrangement to the visibility object.
*/
  void attach (const Input_Arrangement_2 &arr);

  
/*!
Detaches the arrangement from the visibility object it is currently attached to
*/
  void detach();

/*!
Access to the attached arrangement
*/
  Input_Arrangement_2 arr();

/*! 
Computes the visibility region for the given query point q. Visibility polygon of q and its face will be saved to out_arr and out_face.
\pre face is a face of  this->arr()
\pre p is in the interior of face
\pre out_arr is the output arrangement 
\return a face handle to the bounded face of the output
*/ 
  Face_handle visibility_region(const Point_2& q, const Face_handle& face, Output_Arrangement_2& out_arr);

/*! 
Computes for the given query point q the visibility region that is on the side of the given halfedge. Visibility polygon of q and its face will be saved to out_arr and out_face.
\pre half_edge is a half edge of  this->arr()
\pre p is on halfedge
\pre out_arr is the output arrangement  
\return a face handle to the bounded face of the output
*/ 
  Face_handle visibility_region(const Point_2& q, const Halfedge_handle& halfedge, Output_Arrangement_2& out_arr);

/// @}

}; /* end Visibility_2 */

