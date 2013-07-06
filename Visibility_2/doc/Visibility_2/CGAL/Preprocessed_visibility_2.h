namespace CGAL {
/*!
\ingroup PkgVisibility_2Classes

This class is a model of the concept `Visibility_2` offers visibility queries within 
an Arrangement. The algorithm it applies to obtain visibility is using preprocessing.

\cgalModels `Visibility_2` 

\sa 'CGAL::Visibility_2'
\sa `CGAL::Simple_visibility_2<ArrExtensionTraits_2, Regularization_tag>`
\sa `CGAL::Naive_visibility_2<ArrExtensionTraits_2, Regularization_tag>`

*/
template <typename ArrExtensionTraits_2, typename Regularization_tag>
class Preprocessed_visibility_2 {
public:

/// \name Types 
/// @{

  /*!
   The type of Arrangement extension traits.
  */
   typedef ArrExtensionTraits_2 Arr_extension_traits_2;

 /*!
  The type of input Arrangement.
  */
   typedef Arr_extension_traits_2::Input_Arrangement_2  Input_Arrangement_2;

 /*!
  The type of output Arrangement.
  */
   typedef Arr_extension_traits_2::Output_Arrangement_2 Output_Arrangement_2;

 /*! 
   The Point_2 type which is used for queries. 
 */ 
  typedef Input_Arrangement_2::Point_2 Point_2; 

  /*!
   Face_handle type of input Arrangement.
   */
  typedef Input_Arrangement_2::Face_handle Face_handle;

  /*!
   Halfedge_handle type of input Arrangement.
   */
  typedef Input_Arrangement_2::Halfedge_handle Halfedge_handle;

  
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
