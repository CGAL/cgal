
/*!
\ingroup PkgTriangulation2Concepts
\cgalConcept

The concept `TriangulationVertexBase_2` describes the requirements for the 
vertex base class of a triangulation data structure 
to be plugged in a basic, Delaunay or constrained 
triangulations. 

The concept `TriangulationVertexBase_2` refines the concept 
`TriangulationDSVertexBase_2` 
adding geometric information: 
the vertex base of a triangulation stores a point. 

\cgalRefines `TriangulationDSVertexBase_2` 

\cgalHasModel `CGAL::Triangulation_vertex_base_2<Traits,Vb>`

\sa `TriangulationDataStructure_2` 
\sa `TriangulationDataStructure_2::Vertex` 
\sa `CGAL::Triangulation_vertex_base_2<Traits>` 

*/

class TriangulationVertexBase_2 {
public:

/// \name Types 
/// @{

/*!
Must be the same as the point type `TriangulationTraits_2::Point_2` 
defined by the geometric traits class of the triangulation. 
*/ 
typedef unspecified_type Point; 

/// @} 

/// \name Creation 
/// @{

/*!
constructs a vertex embedded in point `p`. 
*/ 
TriangulationVertexBase_2(Point p); 

/*!
constructs a vertex embedded in point `p` and pointing on face `f`. 
*/ 
TriangulationVertexBase_2(Point p, Face_handle f); 

/// @} 

/// \name Access Functions 
/// @{

/*!
returns the point. 
*/ 
Point point() const; 

/// @} 

/// \name Setting 
/// @{

/*!
sets the point. 
*/ 
void set_point(Point p); 

/*!
Inputs the non-combinatorial information given by the vertex: 
the point and other possible information. 
*/ 
istream& operator>> (istream& is, TriangulationVertexBase_2 & v); 

/*!
Outputs the non combinatorial operation given by the vertex: the 
point and other possible information. 
*/ 
ostream& operator<< (ostream& os, const TriangulationVertexBase_2 & v); 

/// @}

}; /* end TriangulationVertexBase_2 */

