/*!
\ingroup PkgConvexHull3Concepts
\cgalConcept

Requirements of the polyhedron type built by the 
function `CGAL::convex_hull_3()`. 

\cgalHasModel `CGAL::Polyhedron_3` 

*/
class ConvexHullPolyhedron_3 {
public:

/// \name Types 
/// @{

/*!
type of point stored in a vertex 
*/ 
typedef unspecified_type Point_3; 

/*!
a model of `ConvexHullPolyhedronVertex_3` 
*/ 
typedef unspecified_type Vertex; 

/*!
a model of `ConvexHullPolyhedronHalfedge_3` 
*/ 
typedef unspecified_type Halfedge; 

/*!
a model of `ConvexHullPolyhedronFacet_3` 
*/ 
typedef unspecified_type Facet; 

/*!
halfedge data structure 
*/ 
typedef unspecified_type Halfedge_data_structure; 

/*!
handle to halfedge 
*/ 
typedef unspecified_type Halfedge_handle; 

/*!
iterator for halfedge 
*/ 
typedef unspecified_type Halfedge_iterator; 

/*!
handle to facet 
*/ 
typedef unspecified_type Facet_handle; 

/*!
iterator for facet 
*/ 
typedef unspecified_type Facet_iterator; 

/// @} 

/// \name Creation 
/// Only a default constructor is required. 
/// @{

/*!

*/ 
ConvexHullPolyhedron_3(); 

/// @} 

/// \name Operations 
/// @{

/*!
iterator over all facets (excluding holes). 
*/ 
Facet_iterator facets_begin(); 

/*!
past-the-end iterator. 
*/ 
Facet_iterator facets_end(); 

/*!
iterator over all halfedges. 
*/ 
Halfedge_iterator P.halfedges_begin(); 

/*!
past-the-end iterator. 
*/ 
Halfedge_iterator P.halfedges_end(); 

/*!

adds a new tetrahedron to the polyhedral surface with its 
vertices initialized with `p1`, `p2`, `p3` and `p4`. 
Returns that halfedge 
of the tetrahedron which incident vertex is initialized with `p1`, the 
incident vertex of the next halfedge with `p2`, and the vertex 
thereafter with `p3`. The remaining fourth vertex is initialized with 
`p4`. 

*/ 
Halfedge_handle make_tetrahedron(Point_3 p1, Point_3 p2, Point_3 p3, Point_3 p4); 

/*!
removes the incident facet of `h` 
and changes all halfedges incident to the facet into border edges or removes 
them from the polyhedral surface if they were already border edges. 
*/ 
void erase_facet(Halfedge_handle h); 

/*!

creates a new facet within the hole incident to `h` and `g` by 
connecting the tip of `g` with the tip of `h` with two new halfedges 
and a new vertex and filling this separated part of the hole with 
a new facet, such that the new facet is incident to `g`. Returns the 
halfedge of the new edge that is incident to the new facet and 
the new vertex. 

*/ 
Halfedge_handle 
add_vertex_and_facet_to_border(Halfedge_handle h, Halfedge_handle g); 

/*!

creates a new facet within the hole incident to `h` and `g` by 
connecting the tip of `g` with the tip of `h` with a new halfedge and 
filling this separated part of the hole with a new facet, such that 
the new facet is incident to `g`. Returns the halfedge of the new 
edge that is incident to the new facet. 

*/ 
Halfedge_handle 
add_facet_to_border(Halfedge_handle h, Halfedge_handle g); 

/*!

fills a hole with a newly created facet. Makes all border 
halfedges of the hole denoted by h incident to the new facet. 
Returns `h`. 

*/ 
Halfedge_handle fill_hole(Halfedge_handle h); 

/*!

calls the `operator()` of the modifier `m`. See `Modifier_base` 
for a description of modifier design and its usage. 

*/ 
void delegate(Modifier_base<Halfedge_data_structure>& m); 

/// @}

}; /* end ConvexHullPolyhedron_3 */
