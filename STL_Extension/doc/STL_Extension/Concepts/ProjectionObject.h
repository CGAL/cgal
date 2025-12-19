


/*!
\ingroup PkgSTLExtensionConcepts
\cgalConcept

\anchor sectionProjectionFunctionObjects


The concept `ProjectionObject` is modeled after the STL
concept `UnaryFunction`, but takes also care of (const)
references.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Identity}
\cgalHasModels{CGAL::Dereference}
\cgalHasModels{CGAL::Get_address}
\cgalHasModels{CGAL::Cast_function_object}
\cgalHasModels{CGAL::Project_vertex}
\cgalHasModels{CGAL::Project_facet}
\cgalHasModels{CGAL::Project_point}
\cgalHasModels{CGAL::Project_normal}
\cgalHasModels{CGAL::Project_plane}
\cgalHasModels{CGAL::Project_next}
\cgalHasModels{CGAL::Project_prev}
\cgalHasModels{CGAL::Project_next_opposite}
\cgalHasModels{CGAL::Project_opposite_prev}
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

