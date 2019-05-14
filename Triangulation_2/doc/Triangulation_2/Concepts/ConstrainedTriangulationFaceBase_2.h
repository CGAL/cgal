
/*!
\ingroup PkgTriangulation2Concepts
\cgalConcept

In a constrained triangulation, 
the information about constrained edges is stored in the 
faces of the triangulation. 
The base face of a constrained triangulation 
has to be a model of the concept 
`ConstrainedTriangulationFaceBase_2` which refines the concept 
`TriangulationFaceBase_2` 
providing functionalities to deal with 
constraints. 

\cgalRefines `TriangulationFaceBase_2` 

\cgalHeading{Types}

Defines the same types as the `TriangulationFaceBase_2` concept 

\cgalHasModel `CGAL::Constrained_triangulation_face_base_2<Traits>` 

\sa `TriangulationFaceBase_2` 
\sa `CGAL::Constrained_triangulation_2<Traits,Tds>` 
\sa `CGAL::Constrained_triangulation_face_base_2<Traits>` 

*/

class ConstrainedTriangulationFaceBase_2 {
public:

/// \name Access Functions 
/// @{

/*!
returns true if the edge between the face and its neighbor 
`neighbor(i)` is constrained. 
\pre \f$ 0\leq i \leq2\f$. 
*/ 
bool is_constrained(int i); 

/// @} 

/// \name Modifiers 
/// @{

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
sets the edge between the face and its neighbor `neighbor(i)` 
as a constrained or unconstrained edge according to `b`.
\cgalAdvancedEnd
*/ 
void set_constraint(int i, bool b); 

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
sets the status (constrained or unconstrained) of the three edges of the face. 
\cgalAdvancedEnd
*/ 
void set_constraints(bool c0, bool c1, bool c2); 

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
changes the orientation of the face by exchanging `vertex(0)` 
with `vertex(1)` and `neighbor(0)` with `neighbor(1)` 
and the corresponding constrained status.
\cgalAdvancedEnd
*/ 
void reorient(); 

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
performs a counterclockwise permutation of the 
vertices, neighbors and constrained status of the face. 
\cgalAdvancedEnd
*/ 
void ccw_permute(); 

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
performs a clockwise permutation of the 
vertices and neighbors and constrained status of the face. 
\cgalAdvancedEnd
*/ 
void cw_permute(); 

/// @} 

/// \name Miscellaneous 
/// @{

/*!

\cgalAdvancedFunction
\cgalAdvancedBegin
tests the validity of the face as a face of a plain
triangulation and additionally checks if the edges of the face are
consistently marked as constrained or unconstrained edges in the face
and its neighbors.
\cgalAdvancedEnd
*/ 
bool is_valid(); 

/// @}

}; /* end ConstrainedTriangulationFaceBase_2 */

