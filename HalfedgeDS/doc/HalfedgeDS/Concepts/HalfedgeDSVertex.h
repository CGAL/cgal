/*!
\ingroup PkgHDSConcepts
\cgalconcept

The concept `HalfedgeDSVertex` defines the requirements for the local `Vertex` 
type in the `HalfedgeDS` concept. It is also required in 
the `Vertex_wrapper<Refs,Traits>` member class template of an 
items class, see the `HalfedgeDSItems` concept. 

A vertex optionally stores a reference to an incident halfedge that 
points to the vertex. A type tag indicates whether the related member 
functions are supported. 
Figure \ref figureHalfedgeDSOptionalMethods 
depicts the relationship between a halfedge and its incident 
halfedges, vertices, and faces. 

For the protection of the integrity of the data structure classes such as 
`CGAL::Polyhedron_3` are allowed to redefine the modifying member 
functions to be private. In order to make them accessible for the 
halfedge data structure they must be derived from a base class `Base` 
where the modifying member functions are still public. (The protection 
could be bypassed by an user, but not by accident.) 

\hasModel `CGAL::HalfedgeDS_vertex_base<Refs>` 
\hasModel `CGAL::HalfedgeDS_vertex_min_base<Refs>` 

\sa `HalfedgeDS<Traits,Items,Alloc>` 
\sa `HalfedgeDSItems` 
\sa `HalfedgeDSHalfedge` 
\sa `HalfedgeDSFace` 

*/

class HalfedgeDSVertex {
public:

/// \name Types 
/// @{

/*! 
instantiated `HalfedgeDS` ( \f$ \equiv\f$ `Refs`). 
*/ 
typedef Hidden_type HalfedgeDS; 

/*! 
base class that allows modifications. 
*/ 
typedef Hidden_type Base; 

/*! 
model of `HalfedgeDSHalfedge`. 
*/ 
typedef Hidden_type Halfedge; 

/*! 
model of `HalfedgeDSFace`. 
*/ 
typedef Hidden_type Face; 

/*! 
handle to vertex. 
*/ 
typedef Hidden_type Vertex_handle; 

/*! 
handle to halfedge. 
*/ 
typedef Hidden_type Halfedge_handle; 

/*! 
handle to face. 
*/ 
typedef Hidden_type Face_handle; 

/*! 

*/ 
typedef Hidden_type Vertex_const_handle; 

/*! 

*/ 
typedef Hidden_type Halfedge_const_handle; 

/*! 

*/ 
typedef Hidden_type Face_const_handle; 

/*! 
`CGAL::Tag_true` or 
`CGAL::Tag_false`. 
*/ 
typedef Hidden_type Supports_vertex_halfedge; 

/// @} 

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
Vertex(); 

/// @} 

/// \name Operations available if `Supports_vertex_halfedge` \f$ \equiv\f$ `CGAL::Tag_true` 
/// @{

/*! 

*/ 
Halfedge_handle halfedge(); 
/*! 

incident halfedge that points to `v`. 
*/ 
Halfedge_const_handle halfedge() const; 
/*! 

sets incident halfedge to \f$ h\f$. 
*/ 
void set_halfedge( Halfedge_handle h); 

/// @}

}; /* end HalfedgeDSVertex */
