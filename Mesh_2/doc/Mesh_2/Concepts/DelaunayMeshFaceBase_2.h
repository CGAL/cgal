


/*!
\ingroup PkgMesh2Concepts
\cgalconcept


The concept `DelaunayMeshFaceBase_2` refines the concept 
`TriangulationFaceBase_2`. It adds two functions giving access 
to a Boolean marker, that indicates if the face is in the 
meshing domain or not. 

\refines ::ConstrainedTriangulationFaceBase_2 

\hasModel `CGAL::Delaunay_mesh_face_base_2<Traits, Fb>` 


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





/// @}

}; /* end DelaunayMeshFaceBase_2 */

