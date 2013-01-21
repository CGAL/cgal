
namespace CGAL {

/*!
\ingroup PkgTriangulation2TriangulationClasses

The class `Polyline_constrained_triangulation_2` 
implements a constrained triangulation 
with an additional data 
structure called the constraint hierarchy 
that keeps track of the input constraints and of their refinement 
in the triangulation. 

The class `Polyline_constrained_triangulation_2<Tr>` 
inherits from its template parameter Tr, which has to be instantiated 
by a constrained or constrained Delaunay triangulation. 

According to its intersection tag, the base class 
will support intersecting input constraints or not. 
When intersections of input constraints are supported, 
the base class constructs a triangulation of the arrangement 
of the constraints, 
introducing new vertices at each proper intersection 
point and refining the input constraints into sub-constraints 
which are edges (more precisely constrained edges) of the 
triangulation. 
In this context, the constraint hierarchy 
keeps track of the input constraints and of their refinement 
in the triangulation. This data structure 
maintains for each input constraints 
the sequence of intersection vertices added on this constraint. 
The constraint hierarchy also allows the user to retrieve the set 
of constrained edges of the triangulation, and for each 
constrained edge, the set of input constraints that overlap it. 

\tparam Tr must be either a CGAL::Constrained_triangulation_2 or a CGAL::Constrained_Delaunay_triangulation_2

\sa `CGAL::Constrained_triangulation_2<Traits,Tds>` 
\sa `CGAL::Constrained_Delaunay_triangulation_2<Traits,Tds>` 
\sa `ConstrainedTriangulationTraits_2` 
\sa `ConstrainedDelaunayTriangulationTraits_2` 

*/
template< typename Tr >
class Polyline_constrained_triangulation_2 : public Tr {
public:

/// \name Types 
/// @{

/*! 
The triangulation base class. 
*/ 
typedef Tr Triangulation; 

/*! 
The intersection tag.
*/ 
typedef Itag Intersection_tag; 
  
/*!
The identifier of a polyline constraint.
*/
  typedef Hidden_type Constraint_id;

/*! 
An iterator to visit 
all the input constraints. The order of visit is arbitrary. 
The value type of this iterator is `Constraint_id` corresponding to the 
endpoints of the constraint. 
*/ 
typedef Hidden_type Constraint_iterator; 

/*! 
An iterator 
to visit all the sub-constraints of the triangulation. 
The order of visit is arbitrary. 
The value type of this iterator is a pair 
`std::pair<Vertex_handle, Vertex_handle>` 
corresponding to the vertices of the 
sub-constraint. 
*/ 
typedef Hidden_type Subconstraint_iterator; 

/*! 
An iterator on the 
vertices of the chain of triangulation edges representing a 
constraint. The value type of this iterator is `Vertex_handle`. 
*/ 
typedef Hidden_type Vertices_in_constraint_iterator; 

/*! 
This type is intended to describe 
a constraint enclosing a sub-constraint and the position of the 
sub-constraint in this constraint. 
It provides three member functions `vertices_begin()`, `vertices_end()` 
and `current()` returning 
iterators of the type `Vertices_in_constraint_iterator` 
on the sequence of vertices of the enclosing constraint. 
These iterators point 
respectively on the first vertex of the enclosing constraint, 
past the last vertex 
and on the first vertex of the sub-constraint. 
*/ 
  class Context {
    /*!
      returns the constraint id.
     */
    Constraint_id id() const;
    /*!
      returns the first vertex of the enclosing constraint.
     */
    Vertices_in_constraint_iterator vertices_begin() const;
  /*!
      returns the past the end of the vertices of the enclosing constraint.
     */
    Vertices_in_constraint_iterator vertices_end() const;
  /*!
    returns the current vertex of the enclosing constraint.
     */
    Vertices_in_constraint_iterator current() const;
  }; 

/*! 
An iterator on 
constraints enclosing a given sub-constraint. The value type of this 
iterator 
is `Context`. 
*/ 
typedef Context_iterator; 

/// @} 

/// \name Creation 
/// @{

/*! 
Introduces an empty triangulation.
*/ 
Polyline_constrained_triangulation_2(const Geom_traits& gt=Geom_traits()); 

/*! 
Copy constructor. 
*/ 
Polyline_constrained_triangulation_2(const 
Polyline_constrained_triangulation_2& ct); 


/*! 
Introduces a constrained triangulation 
from the constraints in the range `[first,last)`.
\tparam InputIterator must be an input iterator with the value type `std::pair<Point,Point>`, `Polygon_2`, or range of points.
\todo Formalize range of points. 
*/ 
template<class InputIterator> 
Polyline_constrained_triangulation_2( 
InputIterator first, 
InputIterator last,
const Geom_traits& gt= Geom_traits()); 

/// @} 

/// \name Assignment 
/// @{

/*! 
Assignment. All the vertices and faces are duplicated. 
The constraint hierarchy is also duplicated. 
*/ 
Polyline_constrained_triangulation_2 operator=(const 
Polyline_constrained_triangulation_2& tr); 

/*! 
The triangulations `tr` and this triangulation are swapped. 
This operation should be preferred to the assignment or to 
the copy constructor if `tr` is deleted after that. 
*/ 
void swap(Polyline_constrained_triangulation_2 tr); 

/// @} 

/// \name Insertion and Removal 
/// The class `Polyline_constrained_triangulation_2` overwrites the
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
\tparam InputIterator must be an input iterator with the value type `Point`. 
*/ 
template < class InputIterator > 
size_type 
insert(InputIterator first, InputIterator last); 

/*! 
Inserts the constraint segment `ab` in the triangulation. 
*/ 
void insert_constraint(Point a, Point b); 

/*! 
Inserts the constraint `c`. 
*/ 
  void push_back(const std::pair<Point,Point>& c); 

/*! 
Inserts a constraint whose endpoints are the vertices 
pointed by `va` and `vb` in the triangulation. 
*/ 
void insert_constraint(Vertex_handle va, Vertex_handle vb); 

/*! 
Removes the constraint `cid`, without removing the points from the triangulation.
*/ 
void remove_constraint(Constraint_id cid); 

/// @} 

/// \name Queries 
/// @{

/*! 
Returns a `Constraint_iterator` pointing on the first 
constraint. 
*/ 
Constraint_iterator constraints_begin() const; 

/*! 
Returns a `Constraint_iterator` pointing past the last 
constraint. 
*/ 
Constraint_iterator constraints_end() const; 

/*! 
Returns a `Subconstraint_iterator` pointing on the first 
sub-constraint. 
*/ 
Subconstraint_iterator subconstraints_begin() const; 

/*! 
Returns a `Subconstraint_iterator` pointing past the last 
sub-constraint. 
*/ 
Subconstraint_iterator subconstraints_end() const; 

/*! 
Returns the number of constraints enclosing the sub-constraint 
`(va,vb)`. 
\pre `va` and `vb` refer to the vertices of a constrained edge of the triangulation. 
*/ 
int number_of_enclosing_constraints(Vertex_handle va, 
Vertex_handle vb) const; 

/*! 
Returns the `Context` relative to one of the constraints 
enclosing the sub-constraint `(va,vb)`. 
\pre `va` and `vb` refer to the vertices of a constrained edge of the triangulation. 
*/ 
Context context(Vertex_handle va, Vertex_handle vb) const; 

/*! 
Returns an iterator pointing on the first `Context` 
of the sequence of `Contexts` 
corresponding to the constraints enclosing the sub-constraint`(va,vb)`. 
\pre `va` and `vb` refer to the vertices of a constrained edge of the triangulation. 
*/ 
Context_iterator contexts_begin(Vertex_handle va, 
Vertex_handle vb) const; 

/*! 
Returns an iterator past the last `Context` 
of the sequence of `Contexts` 
corresponding to the constraints enclosing the `(va,vb)`. 
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

/// @}

/// \name Polyline Simplification
/// The polyline simplification algorithm described in Chapter
/// \ref{chapter-polylinesimplification} uses the following types and 
/// functions.

/// @{

/*!
An iterator on the points of the chain of triangulation edges representing a
constraint. The value type of this iterator is `Point`.
A `Vertices_in_constraint_iterator` can be converted into
a `Points_in_constraint_iterator`, but not the other way round.
*/
typedef Hidden_type Points_in_constraint_iterator;

/*! 
Removes vertex at `viq` from the constraint and the triangulation.
Only the vertex but not the point is removed from the constraint `cid`.
\pre The vertex at `viq` must not be the first or last vertex, 
must not be fixed, and the line segment between the predecessor `vip` and 
the successor `vir` of `viq` must not intersect any constraint.
 */

void
  simplify(Vertices_in_constraint_iterator vip,
           Vertices_in_constraint_iterator viq,
           Vertices_in_constraint_iterator vir);

/*!
Removes the points that were kept in the constraint `cid`.
*/
size_type
remove_points_from_constraint(Constraint_id cid);}


/*!
Removes the points that were kept in the constraints.
 */
void
remove_points_from_constraints();




/// @}


}; /* end Polyline_constrained_triangulation_2 */
} /* end namespace CGAL */
