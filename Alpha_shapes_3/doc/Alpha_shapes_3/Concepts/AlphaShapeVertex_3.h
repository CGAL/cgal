
/*!
\ingroup PkgAlphaShapes3Concepts
\cgalconcept

This concept describe the requirements for the base vertex of an alpha shape. 

\refines `TriangulationVertexBase_3.` 


\hasModel `CGAL::Alpha_shape_vertex_base_3`

\sa `CGAL::Alpha_status` 

*/

class AlphaShapeVertex_3 {
public:

/// \name Types 
/// @{

/*! 
Must be the same as the point type provided by 
the geometric traits class of the triangulation. 
*/ 
typedef Hidden_type Point; 

/*! 
Must be `CGAL::Alpha_status<NT>` 
where `NT` is the number type used in the geometric traits class 
of the triangulation. 
*/ 
typedef Hidden_type Alpha_status; 

/// @} 

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
AlphaShapeVertex_3(); 

/*! 
constructor setting 
the point associated to. 
*/ 
AlphaShapeVertex_3(Point p); 

/*! 
constructor setting the point associated to and an incident cell. 
*/ 
AlphaShapeVertex_3(Point p, const Cell_handle& c); 

/// @} 

/// \name Modifiers 
/// @{

/*! 
Returns a pointer the alpha status of the vertex. 
*/ 
Alpha_status* get_alpha_status(); 

/// @}

}; /* end AlphaShapeVertex_3 */

