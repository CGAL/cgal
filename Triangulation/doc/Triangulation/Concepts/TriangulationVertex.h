
/*!
\ingroup PkgTriangulationsConcepts
\cgalConcept

The concept `TriangulationVertex` describes the requirements on the type used by the 
class `CGAL::Triangulation<TriangulationTraits, TriangulationDataStructure>`, and its derived classes, to 
represent a vertex. 

\cgalRefines `TriangulationDSVertex` 
We only list below the additional specific requirements of ::TriangulationVertex. 
Compared to ::TriangulationDSVertex, the main difference is the addition of 
an association of the vertex with a geometric point. 

\cgalHasModel `CGAL::Triangulation_vertex<TriangulationTraits, Data, TriangulationDSVertex> `

Input/Output 
-------------- 

These operators can be used directly and are called by the I/O 
operator of class `Triangulation`. 

\sa `CGAL::Triangulation_vertex<TriangulationTraits, Data, TriangulationDSVertex>` 
\sa `TriangulationFullCell` 
\sa `CGAL::Triangulation<TriangulationTraits, TriangulationDataStructure>` 

*/

class TriangulationVertex {
public:

/// \name Types 
/// @{

/*! 
The type of the point stored in the vertex. It must be 
the same as the point type `TriangulationTraits::Point` (or its refined 
concepts) when the `TriangulationVertex` is used in the class 
`Triangulation<TriangulationTraits, TriangulationDataStructure>` (or its derived classes). 
*/ 
typedef unspecified_type Point; 

/// @} 

/// \name Creation 
/// @{

/*! 
Constructs a vertex with incident full cell `c`. The vertex is embedded at point `p`. 
*/ 
TriangulationVertex(Full_cell_handle c, const Point & p); 

/*! 
Same as above, but without incident full cell. 
*/ 
TriangulationVertex(const Point & p); 

/*! 
Same as above, but with a default-constructed `Point`. 
*/ 
TriangulationVertex(); 

/// @} 

/// \name Operations 
/// @{

/*! 
The parameter `p` becomes the new geometrical position of the vertex. 
*/ 
void set_point(const Point & p); 

/*! 
Returns the vertex's position. 
*/ 
const Point & point() const; 

/*! 
Inputs the non-combinatorial information given by the vertex, i.e., 
the point and other possible information. 
*/ 
std::istream & operator>>(std::istream & is, TriangulationVertex & v); 

/*! 
Outputs the non-combinatorial information given by the vertex, i.e., 
the point and other possible information. 
*/ 
std::ostream & operator<<(std::ostream & os, const TriangulationVertex & v); 

/// @}

}; /* end TriangulationVertex */

