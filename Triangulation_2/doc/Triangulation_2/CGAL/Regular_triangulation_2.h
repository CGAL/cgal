
namespace CGAL {
/*!
\ingroup PkgTriangulation2TriangulationClasses

The class `Regular_triangulation_2`
is designed to maintain the
regular triangulation of a set of weighted points.

Let \f$ { PW} = \{(p_i, w_i), i = 1, \ldots , n \}\f$ be a set of
weighted points where each \f$ p_i\f$ is a point and each \f$ w_i\f$
is a scalar called the weight of point \f$ p_i\f$.
Alternatively, each weighted point \f$ (p_i, w_i)\f$ can be regarded
as a two dimensional sphere with center \f$ p_i\f$ and radius \f$ r_i=\sqrt{w_i}\f$.

The power diagram of the set \f$ { PW}\f$ is a planar partition
such that each cell corresponds to sphere \f$ (p_i, w_i)\f$ of \f$ { PW}\f$
and is the locus of points \f$ p\f$ whose power with respect to \f$ (p_i, w_i)\f$
is less than its power with respect to any other sphere \f$ (p_j, w_j)\f$
in \f$ { PW}\f$.
The dual of this diagram is a triangulation
whose domain covers the convex hull of the set
\f$ { P}= \{ p_i, i = 1, \ldots , n \}\f$ of center points
and whose vertices are a subset of \f$ { P}\f$.
Such a triangulation is called a regular triangulation.
The three points \f$ p_i, p_j\f$ and \f$ p_k\f$ of \f$ { P}\f$
form a triangle in the regular triangulation of \f$ { PW}\f$
iff there is a point \f$ p\f$ of the plane whose
powers with respect to \f$ (p_i, w_i)\f$, \f$ (p_j, w_j)\f$
and \f$ (p_k, w_k)\f$ are equal and less than the power of \f$ p\f$
with respect to any other sphere in \f$ { PW}\f$.

Let us defined the power product of two weighted points
\f$ (p_i, w_i)\f$ and \f$ (p_j, w_j)\f$ as:
\f[ \Pi(p_i, w_i,p_j, w_j) = p_ip_j ^2 - w_i - w_j . \f]
\f$ \Pi(p_i, w_i,p_j, 0)\f$ is simply the power of point \f$ p_j\f$
with respect to the sphere \f$ (p_i, w_i)\f$, and two weighted points
are said to be orthogonal if their power product is null.
The power circle of three weighted points
\f$ (p_i, w_i)\f$, \f$ (p_j, w_j)\f$
and \f$ (p_k, w_k)\f$ is defined as the unique circle
\f$ (\pi, \omega)\f$ orthogonal to
\f$ (p_i, w_i)\f$, \f$ (p_j, w_j)\f$
and \f$ (p_k, w_k)\f$.

The regular triangulation of the sets \f$ { PW}\f$
satisfies the following <I>regular property</I> (which just reduces to the
Delaunay property when all the weights are null):
a triangle \f$ p_ip_jp_k\f$ of the regular triangulation
of \f$ { PW}\f$ is such that the power product of any weighted point
\f$ (p_l, w_l)\f$ of \f$ { PW}\f$ with the power circle of
\f$ (p_i, w_i)\f$, \f$ (p_j, w_j)\f$ is \f$ (p_k, w_k)\f$ is positive or null.
We call power test of the weighted point \f$ (p_l, w_l)\f$ with respect
to the face \f$ p_ip_jp_k\f$, the predicates testing
the sign of
the power product of \f$ (p_l, w_l)\f$ with respect to
the power circle of
\f$ (p_i, w_i)\f$, \f$ (p_j, w_j)\f$ is \f$ (p_k, w_k)\f$. This power product
is given by the following
determinant
\f[
\left| \begin{array}{cccc}
1 & x_i & y_i & x_i ^2 + y_i ^2 - w_i \\
1 & x_j & y_j & x_j ^2 + y_j ^2 - w_j \\
1 & x_k & y_k & x_k ^2 + y_k ^2 - w_k \\
1 & x_l & y_l & x_l ^2 + y_l ^2 - w_l
\end{array}
\right|
\f]

A pair of neighboring faces \f$ p_ip_jp_k\f$
and \f$ p_ip_jp_l\f$ is said to be locally regular
(with respect to the weights in \f$ { PW}\f$)
if the power test of \f$ (p_l,w_l)\f$ with respect to
\f$ p_ip_jp_k\f$ is positive.
A classical result of computational geometry
establishes that a triangulation of the convex hull of \f$ { P}\f$
such that any pair of neighboring faces is regular with respect
to \f$ { PW}\f$, is a
regular triangulation of \f$ { PW}\f$.

Alternatively, the regular triangulation
of the weighted points set \f$ { PW}\f$
can be obtained as the projection
on the two dimensional plane of the convex hull of the set of three
dimensional points
\f$ { P'}= \{ (p_i,p_i ^2 - w_i ), i = 1, \ldots , n \}\f$.

The vertices of the regular triangulation
of a set of weighted points \f$ { PW}\f$ form only a subset
of the set of center points of \f$ { PW}\f$.
Therefore the insertion of a weighted point in a regular triangulation
does not necessarily imply the creation of a new vertex.
If the new inserted point does not appear as a vertex in the
regular triangulation, it is said to be hidden.

Hidden points
are stored in special vertices called hidden vertices.
A hidden point is considered as hidden by
the facet of the triangulation where its point component is located :
in fact, the hidden point
can appear as vertex of the triangulation
only if this facet is removed.
Each face of a regular triangulation stores
the list of hidden vertices whose points are located in the facet.
When a facet is removed,
points hidden by this facet are reinserted in the triangulation.

\tparam Traits is the geometric traits parameter and must be a model of the concept
`RegularTriangulationTraits_2`.
The concept `RegularTriangulationTraits_2` refines the
concept `TriangulationTraits_2` by adding the type
`Weighted_point_2` to describe weighted points
and the type `Power_side_of_oriented_power_circle_2` to perform
power tests on weighted points.

\tparam Tds  must be a model of `TriangulationDataStructure_2`.
The face base of a regular triangulation has to be a model of the concept
`RegularTriangulationFaceBase_2`. while
the vertex base class has to be a model
of `RegularTriangulationVertexBase_2`.
\cgal provides a default
instantiation for the `Tds` parameter by the class
`Triangulation_data_structure_2 < Regular_triangulation_vertex_base_2<Traits>, Regular_triangulation_face_base_2<Traits> >`.

\sa `CGAL::Triangulation_2<Traits,Tds>`
\sa `TriangulationDataStructure_2`
\sa `RegularTriangulationTraits_2`
\sa `RegularTriangulationFaceBase_2`
\sa `RegularTriangulationVertexBase_2`
\sa `CGAL::Regular_triangulation_face_base_2<Traits>`
\sa `CGAL::Regular_triangulation_vertex_base_2<Traits>`

*/
template< typename Traits, typename Tds >
class Regular_triangulation_2 : public Triangulation_2<Traits,Tds> {
public:

/// \name Types
/// @{

/*!

*/
typedef Traits::Distance Distance;

/*!

*/
typedef Traits::Line Line;

/*!

*/
typedef Traits::Ray Ray;

/*!

*/
typedef Traits::Point_2 Bare_point;

/*!

*/
typedef Traits::Weighted_point_2 Weighted_point;

/*!
An iterator that allows to enumerate the
vertices that are not hidden.
*/
typedef unspecified_type All_vertices_iterator;

/*!
An iterator that allows to enumerate the
finite vertices that are not hidden.
*/
typedef unspecified_type Finite_vertices_iterator;

/*!
An iterator that allows to enumerate the
hidden vertices.
*/
typedef unspecified_type Hidden_vertices_iterator;

/// @}

/// \name Creation
/// @{

/*!
Introduces an empty regular triangulation.
*/
Regular_triangulation_2(const Traits& gt = Traits());

/*!
Copy constructor.
*/
Regular_triangulation_2(const Regular_triangulation_2 &rt);

/*!
Equivalent to constructing an empty triangulation with the optional traits class argument and calling insert(first,last).
*/
template < class InputIterator >
Regular_triangulation_2<Traits,Tds>
Regular_triangulation_2 ( InputIterator first, InputIterator last, Traits gt = Traits());


/// @}

/// \name Insertion and Removal
/// @{

/*!
inserts weighted point `p` in the regular triangulation.
If the point `p` does not appear as a vertex of the triangulation,
the returned vertex is a hidden vertex.
If given the parameter `f` is used as an hint
for the place to start the location process of point `p`.
*/
Vertex_handle insert(const Weighted_point& p, Face_handle f=Face_handle());

/*!
insert a weighted point `p` whose bare-point is assumed to be
located in `lt,loc,li`.
See the description of member function
`Triangulation_2::locate()`.

*/
Vertex_handle insert(const Weighted_point &p,
Locate_type lt,
Face_handle loc, int li );

/*!
Equivalent to `insert(p)`.
*/
Vertex_handle push_back(const Point& p);

/*!
inserts the weighted points in the range `[first,last)`.
It returns the difference of the number of vertices between after and
before the insertions (it may be negative due to hidden points).
Note that this function is not guaranteed to insert the weighted points
following the order of `InputIterator`, as `spatial_sort()`
is used to improve efficiency.
\tparam InputIterator must be an input iterator with the value type \link Regular_triangulation_2::Weighted_point `Weighted_point` \endlink.
*/
template < class InputIterator >
std::ptrdiff_t
insert(InputIterator first, InputIterator last);

/*!
inserts the weighted points in the  range `[first,last)`.
It returns the difference of the number of vertices between after and
before the insertions (it may be negative due to hidden points).
Note that this function is not guaranteed to insert the weighted points
following the order of `WeightedPointWithInfoInputIterator`, as `spatial_sort`
is used to improve efficiency.
Given a pair `(p,i)`, the vertex `v` storing `p` also stores `i`, that is
`v.point() == p` and `v.info() == i`. If several pairs have the same point,
only one vertex is created, one of the objects of type `Vertex::Info` will be stored in the vertex.
\pre `Vertex` must be model of the concept `TriangulationVertexBaseWithInfo_2`.

\tparam WeightedPointWithInfoInputIterator must be an input iterator with value type
`std::pair<%Weighted_point,Vertex::Info>`.
*/
template < class WeightedPointWithInfoInputIterator >
std::ptrdiff_t
insert(WeightedPointWithInfoInputIterator first, WeightedPointWithInfoInputIterator last);

/*!
removes the vertex from the triangulation.
*/
void remove(Vertex_handle v);

/// @}

/// \name Queries
/// @{

/*!
outputs  the faces, boundary edges, and hidden vertices of the
conflict zone of point `p` to output iterators.

\tparam OutputItFaces is an output iterator with `Face_handle` as
value type.
\tparam OutputItBoundaryEdges is an output
iterator with `Edge` as value type.
\tparam OutputItHiddenVertices is an output iterator with
`Vertex_handle` as value type.

This member function outputs in
the container pointed to by `fit` the faces which are in
conflict with point `p`, i.e., the faces whose power circles
have negative power wrt. `p`. It outputs in the container
pointed to by `eit` the boundary of the zone in conflict
with `p`. It inserts the vertices that would be hidden by `p`
into the container pointed to by `vit`. The boundary edges of
the conflict zone are output in counter-clockwise order and each edge
is described through its incident face which is not in conflict with
`p`. The function returns in a `CGAL::Triple` the resulting output
iterators.
\pre `dimension()==2`.
*/
template <class OutputItFaces, class OutputItBoundaryEdges,
class OutputItHiddenVertices>
CGAL::Triple<OutputItFaces,OutputItBoundaryEdges,OutputItHiddenVertices>
get_conflicts_and_boundary_and_hidden_vertices(const Weighted_point
&p, OutputItFaces fit, OutputItBoundaryEdges eit,
OutputItHiddenVertices vit, Face_handle start) const;

/*!
outputs the faces and boundary edges  of the
conflict zone of point `p` to output iterators.

See `get_conflicts_and_boundary_and_hidden_vertices()` for details.

The function returns in a `std::pair` the resulting output iterators.
\pre `dimension()==2`.
*/
template <class OutputItFaces, class OutputItBoundaryEdges>
std::pair<OutputItFaces,OutputItBoundaryEdges>
get_conflicts_and_boundary(const Weighted_point
&p, OutputItFaces fit, OutputItBoundaryEdges eit, Face_handle start) const;

/*!
outputs the faces and hidden vertices  of the
conflict zone of point `p` to output iterators.

See `get_conflicts_and_boundary_and_hidden_vertices()` for details.
The function returns
in a `std::pair` the resulting output iterators.
\pre `dimension()==2`.
*/
template <class OutputItFaces,
class OutputItHiddenVertices>
std::pair<OutputItFaces,OutputItHiddenVertices>
get_conflicts_and_hidden_vertices(const Weighted_point
&p, OutputItFaces fit, OutputItHiddenVertices vit, Face_handle start) const;

/*!
outputs the boundary edges and hidden vertices  of the
conflict zone of point `p` to output iterators.

See `get_conflicts_and_boundary_and_hidden_vertices()` for details.
The function returns in a `std::pair` the resulting output
iterators.
*/
template <class OutputItBoundaryEdges, class OutputItHiddenVertices>
std::pair<OutputItBoundaryEdges,OutputItHiddenVertices>
get_boundary_of_conflicts_and_hidden_vertices(const Weighted_point
&p, OutputItBoundaryEdges eit,
OutputItHiddenVertices vit, Face_handle start) const;

/*!
outputs the faces of the
conflict zone of point `p` to output iterators.
 The function returns the resulting output iterator.
\pre `dimension()==2`.
*/
template <class OutputItFaces>
OutputItFaces
get_conflicts (const Point &p,
OutputItFaces fit,
Face_handle start) const;

/*!
outputs the boundary edges
of the conflict zone of `p` in counterclockwise order
where each edge is described through the incident face
which is not in conflict with `p`.
The function returns the resulting output iterator.
*/
template <class OutputItBoundaryEdges>
OutputItBoundaryEdges
get_boundary_of_conflicts(const Point &p,
OutputItBoundaryEdges eit,
Face_handle start) const;

/*!
outputs the hidden vertices of the conflict zone of `p`
into an output iterator.
 The function returns the resulting output iterator.
*/
template <class OutputItHiddenVertices>
OutputItHiddenVertices
get_hidden_vertices(const Point &p,
OutputItHiddenVertices vit,
Face_handle start) const;

/*!
Returns the vertex of the triangulation which is nearest to `p`
with respect to the power distance. This means that the power of the
query point `p` with respect to the weighted point in the nearest
vertex is smaller than the power of `p` with respect to the
weighted point in any other vertex. Ties are broken arbitrarily. The
default constructed handle is returned if the triangulation is empty.
*/
Vertex_handle nearest_power_vertex(Bare_point p) const;

/// @}

/// \name Access Functions
/// @{

/*!
returns the number of finite vertices that are not hidden.
*/
int number_of_vertices() const;

/*!
returns the number of hidden vertices.
*/
int number_of_hidden_vertices() const;

/*!
starts at an arbitrary hidden vertex.
*/
Hidden_vertices_iterator hidden_vertices_begin() const;

/*!
past the end iterator for the sequence of hidden vertices.
*/
Hidden_vertices_iterator hidden_vertices_end() const;

/*!
starts at an arbitrary unhidden finite vertex
*/
Finite_vertices_iterator finite_vertices_begin() const;

/*!
Past-the-end iterator
*/
Finite_vertices_iterator finite_vertices_end()
const;

/*!
starts at an arbitrary unhidden vertex.
*/
All_vertices_iterator all_vertices_end() const;

/*!
past the end iterator.
*/
All_vertices_iterator all_vertices_begin() const;

/// @}

/// \name Dual Power Diagram
/// The following member functions provide the elements of the dual
/// power diagram.
/// @{

/*!
returns the center of the circle orthogonal to the three weighted
points corresponding to the vertices of face `f`.
\pre `f` is not infinite.
*/
Point weighted_circumcenter(const Face_handle &f) const;

/*!
same as weighted_circumcenter.
*/
Point dual(const Face_handle &f) const;

/*!
If both incident faces are finite, returns a segment whose endpoints are the
duals of each incident face. If only one incident face is finite, returns a
ray whose endpoint is the dual of the finite incident face and supported by
the line which is the bisector of the edge's endpoints. If both incident faces
are infinite, returns the line which is the bisector of the edge's endpoints
otherwise.
*/
Object dual(const Edge &e) const;

/*!
Idem
*/
Object dual(const Edge_circulator& ec) const;

/*!
Idem
*/
Object dual(const Edge_iterator& ei) const;

/*!
output the dual power diagram to stream `ps`.
*/
template < class Stream>
Stream& draw_dual(Stream & ps);

/// @}

/// \name Predicates
/// @{

/*!
Returns the power test of `p` with respect to the
power circle associated with `f`.
*/
Oriented_side
power_test(Face_handle f,
const Weighted_point& p) const;

/// @}

/// \name Miscellaneous
/// @{

/*!
Tests the validity of the triangulation as a
`Triangulation_2` and additionally test the regularity of the
triangulation. This method is useful to debug regular triangulation
algorithms implemented by the user.
*/
bool is_valid(bool verbose = false, int level = 0) const;

/// @}

}; /* end Regular_triangulation_2 */
} /* end namespace CGAL */
