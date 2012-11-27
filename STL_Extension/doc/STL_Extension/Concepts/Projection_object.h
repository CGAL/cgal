


/*!
\ingroup PkgStlExtensionConcepts
\cgalConcept

\anchor sectionProjectionFunctionObjects 


The concept `Projection_object` is modeled after the STL 
concept `UnaryFunction`, but takes also care of (const) 
references. 

\cgalHasModel CGAL::Identity
\cgalHasModel CGAL::Dereference
\cgalHasModel CGAL::Get_address
\cgalHasModel CGAL::Cast_function_object
\cgalHasModel CGAL::Project_vertex
\cgalHasModel CGAL::Project_facet
\cgalHasModel CGAL::Project_point
\cgalHasModel CGAL::Project_normal
\cgalHasModel CGAL::Project_plane
\cgalHasModel CGAL::Project_next
\cgalHasModel CGAL::Project_prev
\cgalHasModel CGAL::Project_next_opposite
\cgalHasModel CGAL::Project_opposite_prev


*/

class Projection_object {
public:


/// \name Definition 
/// @{ 
/*! 
argument type. 
*/ 
typedef Hidden_type argument_type; 
/// @} 


/// \name Definition 
/// @{ 
/*! 
result type. 
*/ 
typedef Hidden_type result_type; 
/// @} 

/// \name Creation 
/// @{ 
/*! 
default constructor. 
*/ 
Projection_object(); 



/// @} 


/// \name Operations 
/// @{ 
/*! 

*/ 
result_type& operator()(argument_type &) const; 



/// @} 


/// \name Operations 
/// @{ 
/*! 

*/ 
const result_type& operator()(const argument_type &) const; 



/// @} 



}; /* end Projection_object */

