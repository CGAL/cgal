


/*!
\ingroup PkgSTLExtensionConcepts
\cgalConcept

\anchor sectionProjectionFunctionObjects


The concept `ProjectionObject` is modeled after the STL
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

class ProjectionObject {
public:


/// \name Definition
/// @{
/*!
argument type.
*/
typedef unspecified_type argument_type;
/// @}


/// \name Definition
/// @{
/*!
result type.
*/
typedef unspecified_type result_type;
/// @}

/// \name Creation
/// @{
/*!
default constructor.
*/
ProjectionObject();



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



}; /* end ProjectionObject */

