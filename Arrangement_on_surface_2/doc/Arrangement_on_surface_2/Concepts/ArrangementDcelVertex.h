
/*!
\ingroup PkgArrangement2ConceptsDCEL
\cgalConcept

A vertex record in a <span class="textsc">Dcel</span> data structure. A vertex is always associated 
with a point. However, the vertex record only stores a pointer to the 
associated point, and the actual `Point` object is stored elsewhere. 

A vertex usually has several halfedges incident to it, such that it is 
possible to access one of these halfedges directly and to traverse all 
incident halfedges around the vertex. However, the <span class="textsc">Dcel</span> may also contain 
isolated vertices that have no incident halfedges. In this case, the vertex 
stores an isolated vertex-information record, indicating the face that 
contains this vertex in its interior. 

\sa `ArrangementDcel` 
\sa `ArrangementDcelHalfedge` 
\sa `ArrangementDcelIsolatedVertex` 

*/

class ArrangementDcelVertex {
public:

/// \name Types 
/// @{

/*!
the corresponding <span class="textsc">Dcel</span> halfedge type. 
*/ 
typedef unspecified_type Halfedge; 

/*!
the corresponding <span class="textsc">Dcel</span> isolated 
vertex-information type. 
*/ 
typedef unspecified_type Isolated_vertex; 

/*!
the point type associated with the vertex. 
*/ 
typedef unspecified_type Point; 

/// @} 

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Arr_dcel_vertex(); 

/*!
assigns `v` with the contents of the `other` vertex. 
*/ 
void assign (const Self& other); 

/// @} 

/// \name Access Functions 
/// All functions below also have `const` counterparts, returning
/// non-mutable pointers or references:
/// @{

/*!
returns whether the vertex is isolated (has no incident halfedges). 
*/ 
bool is_isolated() const; 

/*!
returns an incident halfedge that has `v` as its target. 
\pre `v` is <I>not</I> an isolated vertex. 
*/ 
Halfedge* halfedge(); 

/*!
returns the isolated vertex-information record. 
\pre `v` is an isolated vertex. 
*/ 
Isolated_vertex* isolated_vertex(); 

/*!
returns whether the vertex is not associated with a valid point (i.e.\ it 
lies at infinity). 
*/ 
bool has_null_point () const; 

/*!
returns the associated point. 
\pre `v`() is associated with a valid point. 
*/ 
Point& point(); 

/*!
returns the placement of the \f$ x\f$-coordinate in the parameter space, 
that is, either the left boundary-side, the interior, or the right 
boundary-side. 
*/ 
Arr_parameter_space parameter_space_in_x () const; 

/*!
returns the placement of the \f$ y\f$-coordinate in the parameter space, 
that is, either the bottom boundary-side, the interior, or the top 
boundary-side. 
*/ 
Arr_parameter_space parameter_space_in_y () const; 

/// @} 

/// \name Modifiers 
/// @{

/*!
sets the incident halfedge, marking `v` as a regular vertex. 
*/ 
void set_halfedge (Halfedge* e); 

/*!
sets the isolated vertex-information record, marking `v` 
as an isolated vertex. 
*/ 
void set_isolated_vertex (Isolated_vertex* iv); 

/*!
sets the associated point. 
*/ 
void set_point (Point* p); 

/*!
sets `v` as a vertex on a boundary side. 
\pre Either `inf_x` or `inf_y` is not `ARR_INTERIOR`. 
*/ 
void set_boundary (Arr_parameter_space inf_x, Arr_parameter_space inf_y); 

/// @}

}; /* end ArrangementDcelVertex */

