
/*!
\ingroup PkgAlphaShapes3Concepts
\cgalConcept

This concept describes the requirements for the base vertex of a alpha shape with a fixed value alpha. 

\cgalRefines `TriangulationVertexBase_3` 


\cgalHasModel `CGAL::Fixed_alpha_shape_vertex_base_3`
*/

class FixedAlphaShapeVertex_3 {
public:

/// \name Types 
/// @{

/*!
Must be the same as the point type provided by 
the geometric traits class of the triangulation. 
*/ 
typedef unspecified_type Point; 

/// @} 

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
FixedAlphaShapeVertex_3(); 

/*!
constructor setting the point. 
*/ 
FixedAlphaShapeVertex_3(Point p); 

/*!
constructor setting the point and an incident cell. 
*/ 
FixedAlphaShapeVertex_3(Point p, const Cell_handle& c); 

/// @} 

/// \name Access Functions 
/// @{

/*!
Returns a boolean indicating whether the point is on the convex hull of the point of the triangulation. 
*/ 
bool is_on_chull(); 

/*!
Returns the classification of the vertex. 
*/ 
Classification_type get_classification_type(); 

/// @} 

/// \name Modifiers 
/// @{

/*!
Sets the classification of the vertex. 
*/ 
void set_classification_type(Classification_type type); 

/*!
Sets whether the vertex is on the convex hull. 
*/ 
void is_on_chull(bool b); 

/// @}

}; /* end FixedAlphaShapeVertex_3 */

