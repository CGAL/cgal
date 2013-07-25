
/*!
\ingroup PkgAlphaShape2Concepts
\cgalConcept

\cgalRefines `TriangulationVertexBase_2` 

\cgalHasModel `CGAL::Alpha_shape_vertex_base_2`
*/
class AlphaShapeVertex_2 {
public:

/// \name Types 
/// @{

/*!
A coordinate type. 
The type must provide a copy constructor, assignment, comparison 
operators, negation, multiplication, division and allow the 
declaration and initialization with a small integer constant 
(cf. requirements for number types). An obvious choice would be 
coordinate type of the point class. 
*/ 
typedef unspecified_type FT; 

/// @} 

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
AlphaShapeVertex_2(); 

/*!
constructor setting 
the point. 
*/ 
AlphaShapeVertex_2(Point p); 

/*!
constructor setting the point associated to and an incident face. 
*/ 
AlphaShapeVertex_2(Point p, const Face_handle& ff); 

/// @} 

/// \name Access Functions 
/// @{

/*!
returns two alpha values \f$ \alpha_1 \leq\alpha_2\f$, such as for 
\f$ \alpha\f$ between \f$ \alpha_1\f$ and \f$ \alpha_2\f$, the vertex is attached but singular, and 
for \f$ \alpha\f$ upper \f$ \alpha_2\f$, the vertex is regular. 
*/ 
std::pair< FT, FT > get_range(); 

/// @} 

/// \name Modifiers 
/// @{

/*!
sets the alpha values \f$ \alpha_1 \leq\alpha_2\f$, such as for 
\f$ \alpha\f$ between \f$ \alpha_1\f$ and \f$ \alpha_2\f$, the vertex is attached but singular, and 
for \f$ \alpha\f$ upper \f$ \alpha_2\f$, the vertex is regular. 
*/ 
void set_range(std::pair< FT, FT > I); 

/// @}

}; /* end AlphaShapeVertex_2 */

