
/*!
\ingroup PkgTriangulation2Concepts
\cgalconcept

In a constrained triangulation, 
the information about constrained edges is stored in the 
faces of the triangulation. 
The base face of a constrained triangulation 
has to be a model of the concept 
`ConstrainedTriangulationFaceBase_2` which refines the concept 
`TriangulationFaceBase_2` 
providing functionalities to deal with 
constraints. 

\refines ::TriangulationFaceBase_2 

Types 
-------------- 

Defines the same types as the `TriangulationFaceBase_2` concept 

\hasModel `CGAL::Constrained_triangulation_face_base_2<Traits>` 

\sa `TriangulationFaceBase_2` 
\sa `CGAL::Constrained_triangulation_2<Traits,Tds>` 
\sa `CGAL::Constrained_triangulation_face_base_2<Traits>` 

*/

class ConstrainedTriangulationFaceBase_2 {
public:

/// \name Access Functions 
/// @{

/*! 
returns true if the edge between `f` and its neighbor 
`f`.`neighbor(i)` is constrained. 
\pre \f$ 0\leq i \leq2\f$. 
*/ 
bool is_constrained(int i); 

/// @} 

/// \name Modifiers 
/// @{

/*! 
\advanced sets the edge between `f` and its neighbor `f`.`neighbor(i)` 
as a constrained or unconstrained edge according to `b`. 
*/ 
void set_constraint(int i, bool b); 

/*! 
\advanced sets the status (constrained or unconstrained) of the three edges of `f`. 
*/ 
void set_constraints(bool c0, bool c1, bool c2); 

/*! 
\advanced changes the orientation of `f` by exchanging `vertex(0)` 
with `vertex(1)` and `neighbor(0)` with `neighbor(1)` 
and the corresponding constrained status. 
*/ 
void reorient(); 

/*! 
\advanced performs a counterclockwise permutation of the 
vertices, neighbors and constrained status of `f`. 
*/ 
void ccw_permute(); 

/*! 
\advanced performs a clockwise permutation of the 
vertices and neighbors and constrained status of `f`. 
*/ 
void cw_permute(); 

/// @} 

/// \name Miscellaneous 
/// @{

/*! 

\advanced tests the validity of face `f` as a face of a plain
triangulation and additionally checks if the edges of `f` are
consistently marked as constrained or unconstrained edges in face
`f`and its neighbors.
*/ 
bool is_valid(); 

/// @}

}; /* end ConstrainedTriangulationFaceBase_2 */

