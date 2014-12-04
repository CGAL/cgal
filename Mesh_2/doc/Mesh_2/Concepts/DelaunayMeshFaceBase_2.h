


/*!
\ingroup PkgMesh2Concepts
\cgalConcept


The concept `DelaunayMeshFaceBase_2` refines the concept 
`TriangulationFaceBase_2`. It adds two functions giving access 
to a Boolean marker, that indicates if the face is in the 
meshing domain or not. 

\cgalRefines `ConstrainedTriangulationFaceBase_2` 

\cgalHasModel `CGAL::Delaunay_mesh_face_base_2<Traits, Fb>` 


*/

class DelaunayMeshFaceBase_2 {
public:

/// \name Access Functions
/// @{

/*!
returns true if this face is in the domain to be refined.
*/
bool is_in_domain() const;
/*!
sets if this face is in the domain.
*/
void set_in_domain(const bool b); 


/*!
return true if this face has its circumcenter hidden
by a constrained edge. It does not "see" it,
following the Constrained Delaunay triangulation visibility criterion.
*/
bool is_blind() const;
/*!
sets if this face is blind or not.
*/
void set_blind(const bool b);

/*!
If this face is blind, this function returns
the first constrained edge that prevents it from "seeing"
its circumcenter.
\pre is_blind() returns true
*/
Edge blinding_constraint() const;

/*!
\pre is_blind() returns true
\pre e is a constrained edge
*/
void set_blinding_constraint(const Egde& e);


/// @}

}; /* end DelaunayMeshFaceBase_2 */

