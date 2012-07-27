


/*!
\ingroup PkgStlExtensionConcepts
\cgalconcept

\anchor sectionProjectionFunctionObjects 


The concept `Projection_object` is modeled after the STL 
concept `UnaryFunction`, but takes also care of (const) 
references. 

\hasModel CGAL::Identity
\hasModel CGAL::Dereference
\hasModel CGAL::Get_address
\hasModel CGAL::Cast_function_object
\hasModel CGAL::Project_vertex
\hasModel CGAL::Project_facet
\hasModel CGAL::Project_point
\hasModel CGAL::Project_normal
\hasModel CGAL::Project_plane
\hasModel CGAL::Project_next
\hasModel CGAL::Project_prev
\hasModel CGAL::Project_next_opposite
\hasModel CGAL::Project_opposite_prev


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

