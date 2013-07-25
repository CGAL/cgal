
/*!
\ingroup PkgArrangement2ConceptsTraits
\cgalConcept

A model for the `OverlayTraits` should be able to operate on records (namely, 
vertices, halfedges and faces) of two input <span class="textsc">Dcel</span> classes, named 
`Dcel_A` and `Dcel_B`, and construct the records of an output <span class="textsc">Dcel</span> class, referred to as `Dcel_R`. 

Models for the concept are used by the global `overlay()` function to 
maintain the auxiliary data stored with the <span class="textsc">Dcel</span> records of the resulting 
overlaid arrangement, based on the contents of the input records. 

\cgalHasModel `CGAL::Arr_default_overlay_traits<Arrangement>` 
\cgalHasModel `CGAL::Arr_face_overlay_traits<Arr1,Arr2,ResArr,OvlFaceData>` 

\sa `overlay` 

*/

class OverlayTraits {
public:

/// \name Types 
/// @{

/*!
a constant handle a vertex in `Dcel_A`. 
*/ 
typedef unspecified_type Vertex_handle_A; 

/*!
a constant handle to a halfedge in `Dcel_A`. 
*/ 
typedef unspecified_type Halfedge_handle_A; 

/*!
a constant handle to a face `Dcel_A`. 
*/ 
typedef unspecified_type Face_handle_A; 

/*!
a constant handle to a vertex in `Dcel_B`. 
*/ 
typedef unspecified_type Vertex_handle_B; 

/*!
a constant handle to a halfedge in `Dcel_B`. 
*/ 
typedef unspecified_type Halfedge_handle_B; 

/*!
a constant handle to a face in `Dcel_B`. 
*/ 
typedef unspecified_type Face_handle_B; 

/*!
a handle to a vertex in `Dcel_R`. 
*/ 
typedef unspecified_type Vertex_handle_R; 

/*!
a handle to a halfedge in `Dcel_R`. 
*/ 
typedef unspecified_type Halfedge_handle_R; 

/*!
a handle to a faces in `Dcel_R`. 
*/ 
typedef unspecified_type Face_handle_R; 

/// @} 

/// \name Vertex Creation
/// Whenever a vertex in the overlaid arrangement is created, one of
/// the following functions is called in order to attach the
/// appropriate auxiliary data to this vertex:
/// @{

/*!
constructs the vertex `v` induced by the coinciding vertices 
`v1` and `v2`. 
*/ 
void create_vertex (Vertex_handle_A v1, 
Vertex_handle_B v2, 
Vertex_handle_R v) const; 

/*!
constructs the vertex `v` induced by the vertex `v1` that lies on 
the halfedge `e2`. 
*/ 
void create_vertex (Vertex_handle_A v1, 
Halfedge_handle_B e2, 
Vertex_handle_R v) const; 

/*!
constructs the vertex `v` induced by the vertex `v1` that lies 
inside the face `f2`. 
*/ 
void create_vertex (Vertex_handle_A v1, 
Face_handle_B f2, 
Vertex_handle_R v) const; 

/*!
constructs the vertex `v` induced by the vertex `v2` that lies on 
the halfedge `e1`. 
*/ 
void create_vertex (Halfedge_handle_A e1, 
Vertex_handle_B v2, 
Vertex_handle_R v) const; 

/*!
constructs the vertex `v` induced by the vertex `v2` that lies 
inside the face `f1`. 
*/ 
void create_vertex (Face_handle_A f1, 
Vertex_handle_B v2, 
Vertex_handle_R v) const; 

/*!
constructs the vertex `v` induced by the intersection of the 
halfedges `e1` and `e2`. 
*/ 
void create_vertex (Halfedge_handle_A e1, 
Halfedge_handle_B e2, 
Vertex_handle_R v) const; 
/// @}

/// \name Edge Creation
/// Whenever an edge in the overlaid arrangement is created, one of
/// the following functions is called in order to attach the
/// appropriate auxiliary data to this vertex. Note that an edge is
/// created after both its end-vertices are created, (and the
/// corresponding `create_vertex()` methods were invoked). In all
/// cases, the edge is represented by a halfedge `e` directed in
/// lexicographic decreasing order (from right to left). The
/// `create_edge()` method should attach auxiliary data to the twin
/// halfedge (namely to `e->twin()`) as well:
/// @{

/*!
constructs the halfedge `e` induced by an overlap between the 
halfedges `e1` and `e2`. 
*/ 
void create_edge (Halfedge_handle_A e1, 
Halfedge_handle_B e2, 
Halfedge_handle_R e) const; 

/*!
constructs the halfedge `e` induced by the halfedge `e1` that lies 
inside the face `f2`. 
*/ 
void create_edge (Halfedge_handle_A e1, 
Face_handle_B f2, 
Halfedge_handle_R e) const; 

/*!
constructs the halfedge `e` induced by the halfedge `e2` that lies 
inside the face `f1`. 
*/ 
void create_edge (Face_handle_A f1, 
Halfedge_handle_B e2, 
Halfedge_handle_R e) const; 

/// \name Face Creation
/// The following function is invoked whenever a new face is
/// created. It is guaranteed that all halfedges along the face
/// boundary have already been created an have their auxiliary data
/// fields attached to them:
/// @{

/*!
constructs the face `f` induced by the an overlap between the 
faces `f1` and `f2`. 
*/ 
void create_face (Face_handle_A f1, 
Face_handle_B f2, 
Face_handle_R f) const; 

/// @}

}; /* end OverlayTraits */

