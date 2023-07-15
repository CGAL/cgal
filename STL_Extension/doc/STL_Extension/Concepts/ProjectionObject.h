


/*!
\ingroup PkgSTLExtensionConcepts
\cgalConcept

\anchor sectionProjectionFunctionObjects


The concept `ProjectionObject` is modeled after the STL
concept `UnaryFunction`, but takes also care of (const)
references.

\cgalHasModelsBegin
\cgalModels{CGAL::Identity}
\cgalModels{CGAL::Dereference}
\cgalModels{CGAL::Get_address}
\cgalModels{CGAL::Cast_function_object}
\cgalModels{CGAL::Project_vertex}
\cgalModels{CGAL::Project_facet}
\cgalModels{CGAL::Project_point}
\cgalModels{CGAL::Project_normal}
\cgalModels{CGAL::Project_plane}
\cgalModels{CGAL::Project_next}
\cgalModels{CGAL::Project_prev}
\cgalModels{CGAL::Project_next_opposite}
\cgalModels{CGAL::Project_opposite_prev}
\cgalHasModelsEnd


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

