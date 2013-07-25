
/*!
\ingroup PkgSurfaceMesher3Concepts
\cgalConcept

The concept `SurfaceMeshTriangulation_3` describes 
the triangulation type used by the surface mesher 
`CGAL::make_surface_mesh()` to represent 
the three dimensional triangulation 
embedding the surface mesh. 
Thus, this concept describes the requirements 
for the triangulation type `SurfaceMeshC2T3::Triangulation` 
nested in the model of `SurfaceMeshComplex_2InTriangulation_3` 
plugged as the template parameter `SurfaceMeshC2T3` of 
`CGAL::make_surface_mesh()`. It also describes 
the requirements for the triangulation type 
plugged in the class 
`CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr>`. 

\cgalHasModel Any 3D Delaunay triangulation class of \cgal 

\sa `CGAL::Triangulation_3<TriangulationTraits_3,TriangulationDataStructure_3>` 
\sa `CGAL::Delaunay_triangulation_3<DelaunayTriangulationTraits_3,TriangulationDataStructure_3>` 
\sa `SurfaceMeshComplex_2InTriangulation_3` 
\sa `CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr>` 
\sa `CGAL::make_surface_mesh()` 

*/

class SurfaceMeshTriangulation_3 {
public:

/// \name Types 
/// <I>Vertices</I> and <I>cells</I> of the triangulation are
/// manipulated via handles, which support the two dereference
/// operators and `operator->`. The following iterators allow one to
/// visit all finite vertices, edges and facets of the triangulation.
/// @{

/*!
The point type. It must be DefaultConstructible, CopyConstructible and 
Assignable. 
*/ 
typedef unspecified_type Point; 

/*!
Handle to a data representing a <I>vertex</I>. `Vertex_handle` must be 
a model of `Handle` and its <I>value type</I> must be model of 
`TriangulationDataStructure_3::Vertex`. 
*/ 
typedef unspecified_type Vertex_handle; 

/*!
Handle to a data representing a <I>cell</I>. `Cell_handle` must be a 
model of `Handle` and its <I>value type</I> must be model of 
`TriangulationDataStructure_3::Cell`. 
*/ 
typedef unspecified_type Cell_handle; 

/*!
The edge type. 
*/ 
typedef CGAL::Triple<Cell_handle, int, int> Edge; 

/*!
The facet type. 
*/ 
typedef std::pair<Cell_handle, int> Facet; 

/*!
Iterator over finite vertices 
*/ 
typedef unspecified_type Finite_vertices_iterator; 

/*!
Iterator over finite edges 
*/ 
typedef unspecified_type Finite_edges_iterator; 

/*!
Iterator over finite facets 
*/ 
typedef unspecified_type Finite_facets_iterator; 

/*!
The geometric traits class. Must be a model of 
`DelaunayTriangulationTraits_3`. 
*/ 
typedef unspecified_type Geom_traits; 

/// @} 

/// \name Creation 
/// @{

/*!
%Default constructor. 
*/ 
SurfaceMeshTriangulation_3(); 

/*!
Copy constructor. All vertices and faces are duplicated. 
*/ 
SurfaceMeshTriangulation_3(SurfaceMeshTriangulation_3 tr); 

/// @} 

/// \name Assignment 
/// @{

/*!
The triangulation `tr` is duplicated, and modifying the copy after the 
duplication does not modify the original. The previous triangulation held 
by `t` is deleted. 
*/ 
SurfaceMeshTriangulation_3 & 
operator=(const SurfaceMeshTriangulation_3 & tr); 

/*!
Deletes all finite vertices and all cells of `t`. 
*/ 
void clear(); 

/// @} 

/// \name Access Functions 
/// @{

/*!
Returns the dimension of the affine hull. 
*/ 
int dimension() const; 

/*!
Returns a const reference to a model of 
`DelaunayTriangulationTraits_3`. 
*/ 
const DelaunayTriangulationTraits_3 & geom_traits() const; 

/// @} 

/// \name Voronoi diagram 
/// @{

/*!
Returns the dual of facet `f`.

In dimension 3: either a segment, if the two cells incident to `f` 
are finite, or a ray, if one of them is infinite; 

In dimension 2: a point. 
*/ 
Object dual(Facet f) const; 

/// @} 

/// \name Queries 
/// A point `p` is said to be in conflict with a cell `c` in dimension
/// 3 (resp.\ a facet `f` in dimension 2) iff `t.side_of_sphere(c, p)` 
/// (resp.\ `t.side_of_circle(f, p)`) returns
/// `ON_BOUNDED_SIDE`. The set of cells (resp.\ facets in dimension 2)
/// which are in conflict with `p` is connected, and it forms a
/// hole. 
/// @{

/*!
Computes the conflict hole induced by `p`. The starting cell 
(resp.\ facet) `c` must be in conflict. 

Then this function returns respectively in the output iterators: 

- `cit`: the cells (resp.\ facets) in conflict. 

- `bfit`: the facets (resp.\ edges) on the boundary, that is, the facets 
(resp.\ edges) `(t, i)` where the cell (resp. facet) `t` is in 
conflict, but `t->neighbor(i)` is not. 

- `ifit`: the facets (resp.\ edges) inside the hole, that is, delimiting 
two cells (resp.\ facets) in conflict. 

Returns the `Triple` composed of the resulting output iterators. 
*/ 
template <class OutputIteratorBoundaryFacets, 
class OutputIteratorCells, 
class OutputIteratorInternalFacets> 
Triple<OutputIteratorBoundaryFacets, 
OutputIteratorCells, 
OutputIteratorInternalFacets> 
find_conflicts(Point p, Cell_handle c, 
OutputIteratorBoundaryFacets bfit, 
OutputIteratorCells cit, 
OutputIteratorInternalFacets ifit); 

/// @}

/// \name
/// The following iterators allow the user to visit facets,
/// edges and vertices of the triangulation.
/// @{

/*!
Starts at an arbitrary finite vertex. 
*/ 
Finite_vertices_iterator finite_vertices_begin() const; 

/*!
Past-the-end iterator 
*/ 
Finite_vertices_iterator finite_vertices_end() const; 

/*!
Starts at an arbitrary finite edge.
*/ 
Finite_edges_iterator finite_edges_begin() const; 

/*!
Past-the-end iterator 
*/ 
Finite_edges_iterator finite_edges_end() const; 

/*!
Starts at an arbitrary finite facet. 
*/ 
Finite_facets_iterator finite_facets_begin() const; 

/*!
Past-the-end iterator 
*/ 
Finite_facets_iterator finite_facets_end() const; 

/*!
Copies the `Cell_handle`s of all cells incident to `v` to the output 
iterator `cells`. If `t.dimension() < 3`, then do nothing. 
Returns the resulting output iterator. 
\pre `v != Vertex_handle()`, `t.is_vertex(v)`. 
*/ 
template <class OutputIterator> 
OutputIterator 
incident_cells(Vertex_handle v, OutputIterator cells) const; 

/*!
Copies the `Cell_handle`s of all cells incident to `v` to the output 
iterator `cells`. If `t.dimension() < 3`, then do nothing. 
Returns the resulting output iterator. 
*/ 
template <class OutputIterator> 
OutputIterator 
incident_cells(Vertex_handle v, OutputIterator cells) const; 

/*!
Tests whether `p` is a vertex of `t` by locating `p` in 
the triangulation. If `p` is found, the associated vertex `v` 
is given. 
*/ 
bool is_vertex(const Point & p, Vertex_handle & v) const; 

/*!
Tests whether `(u,v)` is an edge of `t`. If the edge is found, 
it gives a cell `c` having this edge and the indices `i` 
and `j` of the vertices `u` and `v` in `c`, in this order. 
\pre `u` and `v` are vertices of `t`. 
*/ 
bool is_edge(Vertex_handle u, Vertex_handle v, 
Cell_handle & c, int & i, int & j) const; 

/*!
Returns `true`, iff vertex `v` is the infinite vertex. 
*/ 
bool is_infinite(const Vertex_handle v) const; 

/*!
Returns `true`, iff `c` is incident to the infinite vertex. 
\pre `t.dimension() == 3`. 
*/ 
bool is_infinite(const Cell_handle c) const; 

/*!
Returns the same facet seen from the other adjacent cell. 
*/ 
Facet mirror_facet(Facet f) const; 

/*!
Return the indexes of the `j`th vertex of the facet of a cell 
opposite to vertex `i`. 
*/ 
int vertex_triple_index(const int i, const int j); 

/// @} 

/// \name Point location 
/// @{

/*!

If the point `query` lies inside the convex hull of the points, the cell 
that contains the query in its interior is returned. If `query` lies on a 
facet, an edge or on a vertex, one of the cells having `query` on 
its boundary is returned. 

If the point `query` lies outside the convex hull of the points, 
an infinite cell with vertices \f$ \{ p, q, r, \infty\}\f$ is returned such that 
the tetrahedron \f$ ( p, q, r, query )\f$ is positively oriented 
(the rest of the triangulation lies on the other side of facet 
\f$ ( p, q, r )\f$). 

Note that locate works even in degenerate dimensions: in dimension 2 
(resp. 1, 0) the `Cell_handle` returned is the one that represents 
the facet (resp. edge, vertex) containing the query point. 

The optional argument `start` is used as a starting place for the search. 

*/ 
Cell_handle 
locate(const Point & query, Cell_handle start = Cell_handle()) const; 

/*!
If `query` lies inside the affine hull of the points, the `k`-face 
(finite or infinite) that contains `query` in its interior is 
returned, by means of the cell returned together with `lt`, which 
is set to the locate type of the query (`VERTEX, EDGE, FACET, CELL`, or `OUTSIDE_CONVEX_HULL` if the cell is infinite and `query` 
lies strictly in it) and two indices `li` and `lj` that 
specify the `k`-face of the cell containing `query`. 

If the `k`-face is a cell, `li` and `lj` have no 
meaning; if it is a facet (resp. vertex), `li` gives the index of 
the facet (resp. vertex) and `lj` has no meaning; if it is and 
edge, `li` and `lj` give the indices of its vertices. 

If the point `query` lies outside the affine hull of the points, 
which can happen in case of degenerate dimensions, `lt` is set to 
`OUTSIDE_AFFINE_HULL`, and the cell returned has no meaning. 
As a particular case, if there is no finite vertex yet in the 
triangulation, `lt` is set to `OUTSIDE_AFFINE_HULL` and 
`locate` returns the default constructed handle. 

The optional argument `start` is used as a starting place for the search. 

*/ 
Cell_handle 
locate(const Point & query, Locate_type & lt, 
int & li, int & lj, Cell_handle start = Cell_handle() ) const; 

/*!
Creates a new vertex by starring a hole. It takes an iterator range 
`[cell_begin, cell_end)` of `Cell_handle`s which specifies 
a hole: a set of connected cells (resp. facets in dimension 2) which is 
star-shaped wrt `p`. 
(`begin`, `i`) is a facet (resp. an edge) on the boundary of the hole, 
that is, `begin` belongs to the set of cells (resp. facets) previously 
described, and `begin->neighbor(i)` does not. Then this function deletes 
all the cells (resp. facets) describing the hole, creates a new vertex 
`v`, and for each facet (resp. edge) on the boundary of the hole, creates 
a new cell (resp. facet) with `v` as vertex. Then `v->set_point(p)` 
is called and `v` is returned. 

\pre `t.dimension() >= 2`, the set of cells (resp. facets in dimension 2) is connected, its boundary is connected, and `p` lies inside the hole, which is star-shaped wrt `p`. 
*/ 
template <class CellIt> 
Vertex_handle insert_in_hole(Point p, CellIt cell_begin, CellIt cell_end, 
Cell_handle begin, int i); 

/// @}

}; /* end SurfaceMeshTriangulation_3 */

