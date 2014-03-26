
/*!
\ingroup PkgSurfaceParameterizationConcepts
\cgalConcept

`ParameterizationPatchableMesh_3` inherits from concept `ParameterizationMesh_3`, thus is a concept of a 3D surface mesh. 

`ParameterizationPatchableMesh_3` adds the ability to support patches and virtual seams. <I>Patches</I> are a subset of a 3D mesh. <I>Virtual seams</I> are the ability to behave exactly as if the surface was cut following a certain path. 

This mainly means that: vertices can be tagged as inside or outside the patch to parameterize. the fields specific to parameterizations (index, u, v, `is_parameterized`) can be set per <I>corner</I> (aka half-edge). 

The main purpose of this feature is to allow the `Surface_mesh_parameterization` package to parameterize any 3D surface by decomposing it as a list of topological disks. 

Design Pattern 
-------------- 

`ParameterizationPatchableMesh_3` is an Adaptor \cgalCite{cgal:ghjv-dpero-95} : it changes the interface of a 3D mesh to match the interface expected by class `Parameterization_mesh_patch_3`. 

\cgalRefines `ParameterizationMesh_3` 


Creation 
-------------- 

Construction and destruction are undefined. 

\cgalHasModel  Adaptator for `Polyhedron_3` is provided: `CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron_3_>`

\sa `ParameterizationMesh_3`

*/

class ParameterizationPatchableMesh_3 {
public:

/// \name Operations 
/// @{

/*!

%Get vertex seaming flag. Default value is undefined. 

*/ 
int get_vertex_seaming(Vertex_const_handle vertex) const; 

/*!
Set vertex seaming flag. Default value is undefined. 
*/ 
void set_vertex_seaming(Vertex_handle vertex, int seaming); 

/*!

%Get oriented edge's seaming flag, i.e.\ position of the oriented edge w.r.t.\ to the UNIQUE main border. 

*/ 
int get_halfedge_seaming(Vertex_const_handle source, Vertex_const_handle target) const; 

/*!
Set oriented edge's seaming flag, i.e.\ position of the oriented edge w.r.t.\ to the UNIQUE main border. 
*/ 
void set_halfedge_seaming(Vertex_handle source, Vertex_handle target, int seaming); 

/*!

%Get the 2D position (= (u, v) pair) of corners at the <I>right</I> of the `prev_vertex -> vertex -> next_vertex` line. Default value is undefined. 

*/ 
Point_2 get_corners_uv(Vertex_const_handle vertex, Vertex_const_handle prev_vertex, Vertex_const_handle next_vertex) const; 

/*!
Set the 2D position (= (u, v) pair) of corners at the <I>right</I> of the `prev_vertex -> vertex -> next_vertex` line. Default value is undefined. 
*/ 
void set_corners_uv(Vertex_handle vertex, Vertex_const_handle prev_vertex, Vertex_const_handle next_vertex, const Point_2& uv); 

/*!
%Get <I>is parameterized</I> field of corners at the <I>right</I> of the `prev_vertex -> vertex -> next_vertex` line. Default value is undefined. 
*/ 
bool are_corners_parameterized(Vertex_const_handle vertex, Vertex_const_handle prev_vertex, Vertex_const_handle next_vertex) const; 

/*!
Set <I>is parameterized</I> field of corners at the <I>right</I> of the `prev_vertex -> vertex -> next_vertex` line. Default value is undefined. 
*/ 
void set_corners_parameterized(Vertex_handle vertex, Vertex_const_handle prev_vertex, Vertex_const_handle next_vertex, bool parameterized); 

/*!
%Get index of corners at the <I>right</I> of the `prev_vertex -> vertex -> next_vertex` line. Default value is undefined. 
*/ 
int get_corners_index(Vertex_const_handle vertex, Vertex_const_handle prev_vertex, Vertex_const_handle next_vertex) const; 

/*!
Set index of corners at the <I>right</I> of the `prev_vertex -> vertex -> next_vertex` line. Default value is undefined. 
*/ 
void set_corners_index(Vertex_handle vertex, Vertex_const_handle prev_vertex, Vertex_const_handle next_vertex, int index); 

/*!

%Get all purpose tag of corners at the <I>right</I> of the `prev_vertex -> vertex -> next_vertex` line. Default value is undefined. 

*/ 
int get_corners_tag(Vertex_const_handle vertex, Vertex_const_handle prev_vertex, Vertex_const_handle next_vertex) const; 

/*!
Set all purpose tag of corners at the <I>right</I> of the `prev_vertex -> vertex -> next_vertex` line. Default value is undefined. 
*/ 
void set_corners_tag(Vertex_handle vertex, Vertex_const_handle prev_vertex, Vertex_const_handle next_vertex, int tag); 

/// @}

}; /* end ParameterizationPatchableMesh_3 */

