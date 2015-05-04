
namespace CGAL {

/*!
\ingroup PkgTriangulation2TriangulationClasses

The class `Constrained_triangulation_plus_2<Tr>`
provides a constrained triangulation with an additional data
structure
that keeps track of the input constraints and of their refinement
in the triangulation.
The class `Constrained_triangulation_plus_2<Tr>`
inherits from its template parameter Tr, which has to be instantiated
by a constrained or constrained Delaunay triangulation.
The intersection tag of the base class determines whether 
intersecting input constraints are supported or not.
When intersections of input constraints are supported,
the base class constructs a triangulation of the arrangement
of the constraints,
introducing new vertices at each proper intersection
point.


The data structure maintains for each input constraint
the sequence of vertices on this constraint. These vertices are
either vertices of the input constraint or intersection points.

\todo The following description does not match the code
Two consecutive vertices of an input constraint form a *subconstraint*.
A subconstraint is a pair of vertex handles and corresponds to a constrained edge of the
triangulation, which is a pair of a face handle and an index. 

The triangulation also enables the retrieval of the set
of subconstraints of the triangulation (not ordered along constraints).
It further enables the retrieval of the set of input constraints that induce a subconstraint,.
As it is straightforward to obtain a subconstraint from a constrained edge `e`,
one can obtain the input constraints which induce `e`.

\tparam Tr must be either a CGAL::Constrained_triangulation_2 or a CGAL::Constrained_Delaunay_triangulation_2

\sa `CGAL::Constrained_triangulation_2<Traits,Tds>` 
\sa `CGAL::Constrained_Delaunay_triangulation_2<Traits,Tds>` 
\sa `ConstrainedTriangulationTraits_2` 
\sa `ConstrainedDelaunayTriangulationTraits_2` 

*/
template< typename Tr >
class Constrained_triangulation_plus_2 : public Tr {
public:

/// \name Types 
/// @{

/*! 
The triangulation base class. 
*/ 
typedef Tr Triangulation; 

/*! 
The intersection tag as defined in `Tr`.
*/ 
  typedef Tr::Intersection_tag Intersection_tag; 
  
/*!
The identifier of a polyline constraint. 
The class is model of `Assignable`, `CopyConstructible`, `DefaultConstructible`, `LessThanComparable` and  `EqualityComparable`.  

A default constructed `Constraint_id` is a singular value that can not be the ID of a constraint.
*/
  typedef unspecified_type Constraint_id;

/*! 
An iterator to visit 
all the input constraints. The order of visit is undefined. 
The value type of this iterator is `Constraint_id`. 
*/ 
typedef unspecified_type Constraint_iterator; 

/*!
A subconstraint is a pair of vertices that correspond to an `Edge`.
\todo The documentation does not match the code
 */
typedef std::pair<Vertex_handle, Vertex_handle> Subconstraint;

/*! 
An iterator 
to visit all the subconstraints of the triangulation. 
The order of visit is undefined. 
The value type of this iterator is `Subconstraint` 
corresponding to the vertices of the 
subconstraint. 
*/ 
typedef unspecified_type Subconstraint_iterator; 

/*! 
An iterator on the 
vertices of the chain of subconstraints representing a 
constraint. The value type of this iterator is `Vertex_handle`. 
*/ 
typedef unspecified_type Vertices_in_constraint_iterator; 

/*! 
A context enables the access to the vertices of a constraint that pass 
through a subconstraint.

*/ 
  class Context {
  public:
    /*!
      returns the constraint id.
     */
    Constraint_id id() const;

    /*!
      returns the first vertex of the enclosing constraint.
     */
    Vertices_in_constraint_iterator vertices_begin() const;

  /*!
      returns the past-the-end of the vertices of the enclosing constraint.
     */
    Vertices_in_constraint_iterator vertices_end() const;

  /*!
    returns the iterator `vici`  of the enclosing constraint
    for which `*vici` and `*%std::next(vici)`
    correspond to the two vertices of the subconstraint.
     */
    Vertices_in_constraint_iterator current() const;
  }; 

/*! 
An iterator on 
constraints enclosing a given subconstraint. The value type of this 
iterator 
is `Context`. 
*/ 
typedef unspecified_type Context_iterator; 

/// @} 

/// \name Creation 
/// @{

/*! 
Introduces an empty triangulation.
*/ 
Constrained_triangulation_plus_2(const Geom_traits& gt=Geom_traits()); 

/*! 
Copy constructor. 
*/ 
Constrained_triangulation_plus_2(const 
Constrained_triangulation_plus_2& ct); 


/*! 
Introduces a constrained triangulation 
from the constraints in the range `[first,last)`.
\tparam ConstraintIterator must be an `InputIterator` with the value type `std::pair<Point,Point>` or `Segment`. 
*/ 
template<class ConstraintIterator> 
Constrained_triangulation_plus_2( 
ConstraintIterator first, 
ConstraintIterator last,
const Geom_traits& gt= Geom_traits()); 

/// @} 

/// \name Assignment 
/// @{

/*! 
Assignment. All the vertices and faces are duplicated. 
The bidirectional mapping between constraints and subconstraints is also duplicated.
*/ 
Constrained_triangulation_plus_2 operator=(const 
Constrained_triangulation_plus_2& tr); 

/*! 
The triangulations `tr` and this triangulation are swapped. 
This operation should be preferred to the assignment or to 
the copy constructor if `tr` is deleted after that. 
*/ 
void swap(Constrained_triangulation_plus_2 tr); 

/// @} 

/// \name Insertion and Removal 
/// The class `Constrained_triangulation_plus_2` overwrites the
/// following insertion and removal member functions for points and
/// constraints.
/// @{

/*! 
Inserts point `p` as a vertex of the triangulation. 
*/ 
Vertex_handle insert(const Point& p, 
Face_handle start = Face_handle() ); 

/*! 
Inserts point `p` in the triangulation at the location given by `(lt,loc,i)`. 
\sa `Triangulation_2::locate()`
*/ 
Vertex_handle insert(const Point& p, 
Locate_type lt, 
Face_handle loc, int li ); 

/*! 
Equivalent to `insert(p)`. 
*/ 
Vertex_handle push_back(const Point& p); 

/*! 
Inserts the points in the range `[first,last)`.
Returns the number of inserted points. 
\tparam PointIterator must be an `InputIterator` with the value type `Point`. 
*/ 
template < class PointIterator > 
size_type 
insert(PointIterator first, PointIterator last); 

/*! 
Inserts the constraint segment `ab` in the triangulation.
If the two points are equal the point is inserted but no constraint,
and the default constructed `Constraint_id` is returned.
*/ 
Constraint_id insert_constraint(Point a, Point b);

/*! 
Inserts the constraint `c`. 
*/ 
  void push_back(const std::pair<Point,Point>& c); 

/*! 
Inserts a constraint whose endpoints are the vertices 
pointed by `va` and `vb` in the triangulation.
If the two vertex handles are equal no constraint is inserted,
and the default constructed `Constraint_id` is returned.
*/ 
Constraint_id insert_constraint(Vertex_handle va, Vertex_handle vb);

/*!
Inserts a polyline defined by the points in the range `[first,last)`
and returns the constraint id.
The polyline is considered as a closed curve if the first and last point are equal or if  `close == true`. This enables for example passing the vertex range of a `Polygon_2`.
When traversing the vertices of a closed polyline constraint with a  `Vertices_in_constraint_iterator` the first and last vertex are the same.
In case the range is empty `Constraint_id()`is returned.
In case all points are equal the point is inserted but no constraint,
and `Constraint_id()`is returned.
\tparam PointIterator must be an `InputIterator` with the value type `Point`. 
*/
template < class PointIterator>
Constraint_id insert_constraint(PointIterator first, PointIterator last, bool close=false);

/*!
inserts the constraints in the range `[first,last)`.
Note that this function is not guaranteed to insert the constraints
following the order of `ConstraintIterator`, as `spatial_sort()`
is used to improve efficiency.
More precisely, all endpoints are inserted prior to the segments and according to the order provided by the spatial sort.
Once endpoints have been inserted, the segments are inserted in the order of the input iterator,
using the vertex handles of its endpoints.

\return the number of inserted points.
\tparam ConstraintIterator must be an `InputIterator` with the value type `std::pair<Point,Point>` or `Segment`.
*/
template <class ConstraintIterator>
std::size_t insert_constraints(ConstraintIterator first, ConstraintIterator last);

/*!
Same as above except that each constraint is given as a pair of indices of the points
in the range [points_first, points_last). The indices must go from 0 to `std::distance(points_first, points_last)`
\tparam PointIterator is an `InputIterator` with the value type `Point`.
\tparam IndicesIterator is an `InputIterator` with `std::pair<Int, Int>` where `Int` is an integral type implicitly convertible to `std::size_t`
*/
template <class PointIterator, class IndicesIterator>
std::size_t insert_constraints(PointIterator points_first, PointIterator points_last,
                               IndicesIterator indices_first, IndicesIterator indices_last);


/*! 
Removes the constraint `cid`, without removing the points from the triangulation.
*/ 
void remove_constraint(Constraint_id cid); 

/// @} 

/// \name Access 
/// @{

/*! 
Returns a `Constraint_iterator` that points at the first 
constraint of the triangulation. 
*/ 
Constraint_iterator constraints_begin() const; 

/*! 
Returns the past-the-end iterator of the constraints of the triangulation. 
*/ 
Constraint_iterator constraints_end() const; 

/*! 
Returns a `Subconstraint_iterator` pointing at the first 
subconstraint of the triangulation. 
*/ 
Subconstraint_iterator subconstraints_begin() const; 

/*! 
Returns the past-the-end iterator of the subconstraints of the triangulation. 
*/ 
Subconstraint_iterator subconstraints_end() const; 

/*! 
Returns the number of constraints enclosing the subconstraint 
`(va,vb)`. 
\pre `va` and `vb` refer to the vertices of a constrained edge of the triangulation. 
*/ 
int number_of_enclosing_constraints(Vertex_handle va, 
Vertex_handle vb) const; 

/*! 
Returns the `Context` relative to one of the constraints 
enclosing the subconstraint `(va,vb)`. 
\pre `va` and `vb` refer to the vertices of a constrained edge of the triangulation. 
*/ 
Context context(Vertex_handle va, Vertex_handle vb) const; 

/*! 
Returns an iterator pointing at the first `Context` 
of the sequence of contexts
corresponding to the constraints enclosing the subconstraint `(va,vb)`. 
\pre `va` and `vb` refer to the vertices of a constrained edge of the triangulation. 
*/ 
Context_iterator contexts_begin(Vertex_handle va, 
Vertex_handle vb) const; 

/*!
Returns an iterator past the end `Context` 
of the sequence of contexts
corresponding to the constraints enclosing the subconstraint `(va,vb)`. 
\pre `va` and `vb` refer to the vertices of a constrained edge of the triangulation. 
*/ 
Context_iterator contexts_end(Vertex_handle va, 
Vertex_handle vb) const; 

/*!
Returns an iterator on the first vertex on the constraint `cid`. 
*/ 
Vertices_in_constraint_iterator 
vertices_in_constraint_begin(Constraint_id cid) const; 

/*!
Returns an iterator past the last vertex on the constraint `cid`. 
*/ 
Vertices_in_constraint_iterator 
vertices_in_constraint_end(Constraint_id cid) const; 

/*!
\deprecated Returns an iterator on the first vertex on the constraint `cid`.
\attention This function only works for constraints that were inserted
as a pair of points or pair of vertex handles.
*/
Vertices_in_constraint_iterator
vertices_in_constraint_begin(Vertex_handle va, Vertex_handle vb) const; 

/*!
\deprecated Returns an iterator past the last vertex on the constraint `cid`.
\attention This function only works for constraints that were inserted
as a pair of points or pair of vertex handles.
*/
Vertices_in_constraint_iterator
vertices_in_constraint_end(Vertex_handle va, Vertex_handle vb) const; 

/// @}


/*! \name Polyline Simplification
\cgalAdvancedBegin 
The polyline simplification algorithm described in Chapter
\ref Chapter_2D_Polyline_simplification 
operates on polyline constraints. The algorithm removes
in each simplification step
a vertex of a constraint and at the same time from the triangulation.
The class `Constrained_triangulation_plus_2` stores
for each constraint not only the sequence of vertices but
also the original sequence of points at those vertices.
As the `Vertices_in_constraint_iterator` enables the traversal of
the current set of vertices, the `Points_in_constraint_iterator` 
enables the traversal of the points that were in the constraint 
before the simplification algorithm started.

It enables the simplification algorithm to compute the error introduced by 
each simplification step:
it is the distance of the current sequence (vertices) to the original
sequence (points).

Those stored points which do not correspond to a vertex can be removed 
afterward either for a single constraint or for all constraints.

The simplification algorithm uses the following types and functions.
\cgalAdvancedEnd
*/

/// @{

/*!
\cgalAdvancedBegin
An iterator on the points of the of the original constraint
before simplification steps are applied. The value type of this iterator is `Point`.
A \link Constrained_triangulation_plus_2::Vertices_in_constraint_iterator `Vertices_in_constraint_iterator`\endlink can be converted into
a `Points_in_constraint_iterator`, but not the other way around.
\cgalAdvancedEnd
*/
typedef unspecified_type Points_in_constraint_iterator;


/*!
\cgalAdvancedBegin
Returns an iterator to the first point on the constraint before any simplification step.
\cgalAdvancedEnd
*/
Points_in_constraint_iterator points_in_constraint_begin(Constraint_id cid) const;
  
/*!
\cgalAdvancedBegin
Returns an iterator past the last point on the constraint before any simplification step.
\cgalAdvancedEnd
*/
Points_in_constraint_iterator points_in_constraint_end(Constraint_id cid) const ;

/*!
\cgalAdvancedBegin
Removes the vertex at `vicq` from the constraint and the triangulation.
The point of that vertex remains stored in the sequence of original points
of the constraint until `remove_points_without_corresponding_vertex(Constraint_id)`
or `remove_points_without_corresponding_vertex()` is called.

The polyline simplification algorithm described in Chapter
\ref Chapter_2D_Polyline_simplification 
operates on polyline constraints and applies `simplify()` to vertices in
constraints based on a cost and stop function.

\pre Each vertex of the triangulation must be either a vertex of a constraint or a vertex at the intersection of constraints.
\pre `vicq` must neither be the first nor the last vertex on a constraint.
\pre The vertex referred by vicq is not contained in any other constraint.
\pre Let `vip` and `vir` be defined as `vip = std::prev(vicq)` and `vir = std::next(vicr)`.
\pre The line segment between `*vicp->point()` and `*vicr->point()` must not intersect any constraint.
\cgalAdvancedEnd
 */
void
simplify(Vertices_in_constraint_iterator vicq);

/*!
\cgalAdvancedBegin
Removes the original points that correspond to vertices in the constraint `cid` which have
been removed by the `simplify()` function.
\cgalAdvancedEnd
*/
size_type
remove_points_without_corresponding_vertex(Constraint_id cid);


/*!
\cgalAdvancedBegin
Removes all original points that correspond to vertices in the constraints which have
been removed by the `simplify()` function.
\cgalAdvancedEnd
 */
void
remove_points_without_corresponding_vertex();




/// @}


}; /* end Constrained_triangulation_plus_2 */
} /* end namespace CGAL */
