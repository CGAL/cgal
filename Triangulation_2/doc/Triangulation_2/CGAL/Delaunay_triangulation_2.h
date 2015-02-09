
namespace CGAL {

/*!
\ingroup PkgTriangulation2TriangulationClasses

The class `Delaunay_triangulation_2` is designed to represent 
the Delaunay triangulation of a set of points in a plane. 
A Delaunay triangulation of a set of points 
is a triangulation of the sets of points that fulfills 
the following <I>empty circle property</I> 
(also called <I>Delaunay property</I>): the circumscribing 
circle of any facet 
of the triangulation contains no point of the set in its interior. 
For a point set with no case of co-circularity 
of more than three points, 
the Delaunay triangulation is unique, it is the dual 
of the Voronoi diagram of the points. 

\tparam Tds must be a model of `TriangulationDataStructure_2`. 
\cgal provides a default instantiation for this parameter, 
which is the class 
`CGAL::Triangulation_data_structure_2 < CGAL::Triangulation_vertex_base_2<Traits>, CGAL::Triangulation_face_base_2<Traits> >`. 

\tparam Traits must be a model of `DelaunayTriangulationTraits_2`. 
The concept `DelaunayTriangulationTraits_2` refines the 
concept `TriangulationTraits_2`, providing 
a predicate type to check the empty circle property. 

Changing this predicate type 
allows the user to build Delaunay triangulations for different metrics 
such that \f$ L_1\f$ or \f$ L_{\infty}\f$ or any metric defined by a 
convex object. However, the user of an exotic metric 
must be careful that the constructed triangulation 
has to be a triangulation of the convex hull 
which means that convex hull edges have to be Delaunay edges. 
This is granted for any smooth convex metric (like \f$ L_2\f$) 
and can be ensured for other metrics (like \f$ L_{\infty}\f$) 
by the addition to the point set of well chosen sentinel points.
The concept of `DelaunayTriangulationTraits_2` is described ::DelaunayTriangulationTraits_2.

When dealing with a large triangulations, the user is advised to
encapsulate the Delaunay triangulation class into a triangulation
hierarchy, which means to use the class
`Triangulation_hierarchy_2<Tr>` with the template parameter
instantiated with `Delaunay_triangulation_2`. The triangulation
hierarchy will then offer the same functionalities but be much more
for efficient for locations and insertions.

\cgalHeading{Types}

All the types defined in `Triangulation_2<Traits,Tds>` are inherited. 

\cgalHeading{Implementation}

Insertion is implemented by inserting in the triangulation, then 
performing a sequence of Delaunay flips. The number of flips is \f$ O(d)\f$ 
if the new vertex is of degree \f$ d\f$ in the new triangulation. For 
points distributed uniformly at random, insertion takes time \f$ O(1)\f$ on 
average. 

Removal calls the removal in the triangulation and then re-triangulates 
the hole in such a way that the Delaunay criterion is satisfied. Removal of a 
vertex of degree \f$ d\f$ takes time \f$ O(d^2)\f$. 
The degree \f$ d\f$ is \f$ O(1)\f$ for a random 
vertex in the triangulation. 

After a point location step, the nearest neighbor 
is found in time \f$ O(n)\f$ in the 
worst case, but in time \f$ O(1)\f$ 
for vertices distributed uniformly at random and any query point. 

\sa `CGAL::Triangulation_2<Traits,Tds>`
\sa `TriangulationDataStructure_2`
\sa `DelaunayTriangulationTraits_2`
\sa `Triangulation_hierarchy_2<Tr>`

*/
template< typename Traits, typename Tds >
class Delaunay_triangulation_2 : public Triangulation_2<Traits,Tds> {
public:

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Delaunay_triangulation_2(const Traits& gt = 
Traits()); 

/*!
copy constructor. All the vertices and faces are duplicated. 
*/ 
Delaunay_triangulation_2( 
const Delaunay_triangulation_2<Traits,Tds> &tr); 

/*!
Equivalent to constructing an empty triangulation with the optional traits class argument and calling insert(first,last).
*/
template < class InputIterator >
Delaunay_triangulation_2<Traits,Tds> dt ( InputIterator first, InputIterator last, Traits gt = Traits());


/// @} 

/// \name Insertion and Removal 
/// The following insertion and removal functions overwrite the
/// functions inherited from the class `Triangulation_2<Traits,Tds>`
/// to maintain the Delaunay property. In the degenerate case when
/// there are co-circular points, the Delaunay triangulation is known
/// not to be uniquely defined. In this case, \cgal chooses a
/// particular Delaunay triangulation using a symbolic perturbation
/// scheme \cgalCite{cgal:dt-pvr3d-03}. Note that the other modifier
/// functions of `Triangulation_2<Traits,Tds>` are not
/// overwritten. Thus a call to \link Triangulation_2::insert_in_face()\endlink  
/// \link Triangulation_2::insert_in_face() `insert_in_face()`\endlink 
/// \link Triangulation_2::insert_in_edge() `insert_in_edge()`\endlink,
/// \link Triangulation_2::insert_outside_convex_hull() `insert_outside_convex_hull()`\endlink, 
/// \link Triangulation_2::insert_outside_affine_hull() `insert_outside_affine_hull()`\endlink or
///  \link Triangulation_2::flip() `flip()`\endlink on a valid Delaunay triangulation might lead to a
/// triangulation which is no longer a Delaunay triangulation.  
/// @{

/*!
inserts point `p`. 
If point `p` coincides with an already existing vertex, this 
vertex is returned and the triangulation is not updated. 
Optional parameter `f` is used to initialize the location of `p`. 

*/ 
Vertex_handle insert(const Point& p, Face_handle f=Face_handle()); 

/*!
inserts a point `p` at the location given by `(lt,loc,li)`. 
\sa `Triangulation_2::locate()`
*/ 
Vertex_handle insert(const Point& p, Locate_type& lt, 
Face_handle loc, int li ); 

/*!
equivalent to `insert(p)`. 
*/ 
Vertex_handle push_back(const Point& p); 

/*!
inserts the points in the range `[first,last)`.
Returns the number of inserted points. 
Note that this function is not guaranteed to insert the points 
following the order of `PointInputIterator`, as `spatial_sort()` 
is used to improve efficiency. 
\tparam PointInputIterator must be an input iterator with the value type `Point`. 
*/ 
template < class PointInputIterator > 
std::ptrdiff_t 
insert(PointInputIterator first, PointInputIterator last); 

/*!

inserts the points in the iterator range `[first,last)`. Returns the number of inserted points. 
Note that this function is not guaranteed to insert the points 
following the order of `PointWithInfoInputIterator`, as `spatial_sort()` 
is used to improve efficiency. 
Given a pair `(p,i)`, the vertex `v` storing `p` also stores `i`, that is 
`v.point() == p` and `v.info() == i`. If several pairs have the same point, 
only one vertex is created, and one of the objects of type `Vertex::Info` will be stored in the vertex. 
\pre `Vertex` must be model of the concept `TriangulationVertexBaseWithInfo_2`. 

\tparam PointWithInfoInputIterator must be an input iterator with the value type `std::pair<Point,Vertex::Info>`. 

*/ 
template < class PointWithInfoInputIterator > 
std::ptrdiff_t 
insert(PointWithInfoInputIterator first, PointWithInfoInputIterator last); 

/*!
removes the vertex from the triangulation. 
*/ 
void remove(Vertex_handle v); 

/// @} 

/// \name Displacement 
/// @{

/*!
if there is not already another vertex placed on `p`, 
the triangulation is modified such that the new position of vertex `v` 
is `p`, and `v` is returned. Otherwise, the triangulation is not 
modified and the vertex at point `p` is returned. 
\pre Vertex `v` must be finite. 
*/ 
Vertex_handle move_if_no_collision(Vertex_handle v, const Point & p); 

/*!
same as `move_if_no_collision()`, if there is no collision. Otherwise, `v` 
is deleted and the vertex placed on `p` is returned. 
\pre Vertex `v` must be finite. 
*/ 
Vertex_handle move(Vertex_handle v, const Point & p); 

/// @} 

/// \name Queries 
/// @{

/*!
returns any nearest vertex of `p`. The implemented function 
begins with a location step and 
`f` may be used to initialize the location. 
*/ 
Vertex_handle 
nearest_vertex(const Point& p, Face_handle f=Face_handle()) const; 

/*!
outputs the faces and boundary edges of the conflict zone of point `p` into
output iterators.

This function outputs in the container pointed to by `fit` 
the faces which are in conflict with point `p`, 
i. e., the faces whose circumcircle contains `p`. 
It outputs in the container pointed to by `eit` the 
the boundary of the zone in conflict with `p`. 
The boundary edges 
of the conflict zone are output in counter-clockwise order 
and each edge is described through its incident face 
which is not in conflict with `p`. 
The function returns in a `std::pair` the resulting output iterators. 


\tparam OutItFaces is an output iterator with `Face_handle` as value type. 
\tparam OutItBoundaryEdges is an output iterator with `Edge` as value type. 
\pre `dimension()==2`. 
*/ 
template <class OutputItFaces, class OutputItBoundaryEdges> 
std::pair<OutputItFaces,OutputItBoundaryEdges> 
get_conflicts_and_boundary(const Point &p, 
OutputItFaces fit, 
OutputItBoundaryEdges eit, 
Face_handle start) const; 

/*!
outputs the faces of the conflict zone of point `p` into an output iterator.
same as `get_conflicts_and_boundary()` except that only the faces in conflict with `p` 
are output. The function returns the resulting output iterator. 
\pre `dimension()==2`. 
*/ 
template <class OutputItFaces> 
OutputItFaces 
get_conflicts (const Point &p, 
OutputItFaces fit, 
Face_handle start) const; 

/*!
outputs the boundary edges of the conflict zone of point `p` into an output iterator.

This function outputs in the container pointed to by `eit`, 
the boundary of the zone in conflict with `p`. The boundary edges 
of the conflict zone are output in counterclockwise order 
and each edge is described through the incident face 
which is not in conflict with `p`. 
The function returns the resulting output iterator. 
\tparam OutputItBoundaryEdges is an output iterator with 
`Edge` as value type. 
*/ 
template <class OutputItBoundaryEdges> 
OutputItBoundaryEdges 
get_boundary_of_conflicts(const Point &p, 
OutputItBoundaryEdges eit, 
Face_handle start) const; 

/// @} 

/// \name Voronoi Diagram 
/// The following member functions provide the elements of the dual
/// Voronoi diagram.
/// @{

/*!
Returns the center of the circle circumscribed to face `f`. 
\pre `f` is not infinite. 
*/ 
Point dual(const Face_handle &f) const; 

/*!
returns a segment, a ray or a line supported by the bisector of the 
endpoints of `e`. 
If faces incident to `e` are both finite, a segment whose endpoints are the 
duals of each incident face is returned. If only one incident face is 
finite, a 
ray whose endpoint is the dual of the finite incident face is returned. 
Otherwise both incident faces 
are infinite and the bisector line is returned. 
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
output the dual Voronoi diagram to stream `ps`. 
*/ 
template < class Stream> 
Stream& draw_dual(Stream & ps); 

/// @} 

/// \name Predicates 
/// @{

/*!
Returns the side of `p` with respect to the circle circumscribing 
the triangle associated with `f`. 
*/ 
Oriented_side side_of_oriented_circle(Face_handle f, const Point& p) const; 

/// @} 

/// \name Miscellaneous 
/// The checking function `is_valid()` is also overwritten to
/// additionally test the empty circle property.
/// @{

/*!
tests the validity of the triangulation as a
`Triangulation_2` and additionally tests the Delaunay property. This
method is mainly useful for debugging Delaunay triangulation
algorithms designed by the user.
*/ 
bool is_valid(bool verbose = false, int level = 0) const; 

/// @}

}; /* end Delaunay_triangulation_2 */
} /* end namespace CGAL */
