
namespace CGAL {


/*!
\ingroup PkgTriangulation2TriangulationClasses

Intersection tag for constrained triangulations, when input constraints do not intersect.
*/
struct No_intersection_tag{};


/*!
\ingroup PkgTriangulation2TriangulationClasses

Intersection tag for constrained triangulations, if exact predicates and exact constructions are provided.
*/
struct Exact_intersections_tag{};


/*!
\ingroup PkgTriangulation2TriangulationClasses

Intersection tag for constrained triangulations, if the geometric traits provides exact predicates but approximate constructions.
*/

struct Exact_predicates_tag{};

/*!
\ingroup PkgTriangulation2TriangulationClasses

A constrained triangulation is a triangulation of a set of points 
which has to include among its edges 
a given set of polylines joining the points.

The given polylines are 
called *constraints* and the corresponding 
edges in the triangulation are called *constrained edges*. 


The endpoints of constrained edges are of course vertices of the 
triangulation. However the triangulation may include 
other vertices as well. 
There are three versions of constrained triangulations 
<UL> 
<LI>In the basic version, the constrained triangulation 
does not handle intersecting constraints, and the set of input 
constraints is required to be a set of polylines that do not intersect 
except possibly at their endpoints. Any number of constrained edges 
may share the same endpoint. Constrained edges may be vertical or
have zero length. 
<LI>The two other versions support intersecting input constraints. 
In those versions, input constraints may consist of 
intersecting, overlapping or partially 
overlapping segments. 
The triangulation introduces additional vertices at each point which 
is a proper intersection point of two 
constraints. A single constraint intersecting other 
constraints will then appear as the union of several 
constrained edges of the triangulation. 
There are two ways to deal with intersecting constraints. 
<UL> 
<LI>The first one is robust when predicates are evaluated exactly but 
constructions (i. e. intersection computations) are 
approximate. 
<LI>The second one should be used with exact arithmetic (meaning exact
evaluation of predicates and exact computation of intersections.)
</UL> 
</UL> 
In order to retrieve the constrained edges of a constraint, or
the constraints overlapping with a constrained edge, we provide
the class `Constrained_triangulation_plus_2`. This class maintains
a constraint hierarchy data structure. See
Section \ref Section_2D_Triangulations_Constrained_Plus for details.
This class should also be used when doing exact intersection computations
as it avoids the cascading of intersection computations.

\image html constraints.png
\image latex constraints.png


\tparam Traits is a geometric traits class and must be a model 
of the concept `TriangulationTraits_2`. 
When intersection of input constraints are supported, 
the geometric traits class 
is required to provide additional function object types 
to compute the intersection of two segments. 
It has then to be a model of the concept 
`ConstrainedTriangulationTraits_2`. 


\tparam Tds must be a model 
of the concept `TriangulationDataStructure_2`. 

\tparam Itag is the intersection tag 
which serves to choose between the different 
strategies to deal with constraints intersections. 
\cgal provides three valid types for this parameter: 
- `No_intersection_tag` disallows intersections of 
input constraints, 
- `Exact_predicates_tag` is to be used when the traits 
class 
provides exact predicates but approximate constructions of the 
intersection points. 
- `Exact_intersections_tag` is to be used in conjunction 
with an exact arithmetic type. 

The information about constrained edges is stored in the 
faces of the triangulation. Thus the nested `Face` 
type of a constrained triangulation offers 
additional functionalities to deal with this information. 
These additional functionalities 
induce additional requirements on the face base class 
plugged into the triangulation data structure of 
a constrained Delaunay triangulation. 
The face base of a constrained Delaunay triangulation 
has to be a model of the concept 
`ConstrainedTriangulationFaceBase_2`. 

\cgal provides default instantiations for the template parameters 
`Tds` and `Itag`, and for the `ConstrainedTriangulationFaceBase_2`. 
If `Gt` is the geometric traits class 
parameter, 
the default for 
`ConstrainedTriangulationFaceBase_2` is the class 
`Constrained_triangulation_face_base_2<Gt>` 
and the default for the 
triangulation data structure parameter is the class 
`Triangulation_data_structure_2 < Triangulation_vertex_base_2<Gt>, Constrained_triangulation_face_base_2<Gt> >`. 
The default intersection tag is `No_intersection_tag`. 

\sa `CGAL::Triangulation_2<Traits,Tds>`
\sa `TriangulationDataStructure_2`
\sa `TriangulationTraits_2` 
\sa `ConstrainedTriangulationTraits_2` 
\sa `ConstrainedTriangulationFaceBase_2` 

\cgalHeading{Implementation}

The insertion of a constrained edge runs in time 
proportional to the number of triangles intersected by this edge. 

*/
template< typename Traits, typename Tds, typename Itag >
class Constrained_triangulation_2 : public Triangulation_2<Traits,Tds> {
public:

/// \name Types 
/// @{


/*! 
\deprecated The type of constraints.
*/ 
typedef std::pair<Point,Point> Constraint; 

/*! 
A bidirectional iterator to visit 
all the edges `e` of the triangulation which are constrained. 
The order of visit is undefined. 
The value type of this iterator is `Edge`. 
*/ 
typedef unspecified_type Constrained_edges_iterator; 


/*!
The intersection tag which decides how 
intersections between input constraints are dealt with. 
*/ 
typedef Itag Intersection_tag; 

/// @} 

/// \name Creation 
/// @{

/*!
%Default constructor. 
*/ 
Constrained_triangulation_2(); 

/*!
Copy constructor: All faces and vertices 
are duplicated and the constrained status of edges 
is copied. 
*/ 
Constrained_triangulation_2(const 
Constrained_triangulation_2& ct1); 


/// @} 

/// \name Queries 
/// @{

/*!
Returns `true` if edge `e` is a constrained edge. 
*/ 
bool is_constrained(Edge e) const; 

/*!
Returns `true` if at least one of the edges incident to vertex `v` 
is constrained. 
*/ 
bool are_there_incident_constraints(Vertex_handle v) const; 

/*!
Outputs the constrained edges incident to `v` 
into the output iterator `out` and returns the resulting 
output iterator. 
\tparam OutputItEdges is an `OutputIterator` with `Edge` as value 
type. 
*/ 
template<class OutputItEdges> 
OutputItEdges incident_constraints(Vertex_handle v, 
OutputItEdges out) const; 


/*!
returns an iterator that enumerates the constrained edges. 
*/
Constrained_edges_iterator constrained_edges_begin() const;

/*!
returns the past-the-end iterator. 
*/ 
Constrained_edges_iterator constrained_edges_end() const;

/// @} 

/// \name Insertion and Removal 
///
/// @{

/*!
Inserts point `p` and restores the status (constrained or not) of all 
the touched edges. If present, `f` is used as an hint 
for the location of `p`. 
*/ 
Vertex_handle insert(Point p, Face_handle f = Face_handle() ); 

/*!
Inserts point `p` in the triangulation at the location given by `(lt,loc,i)`. 
\sa `Triangulation_2::locate()`
*/ 
Vertex_handle 
insert(const Point& p, 
Locate_type& lt, 
Face_handle loc, int li ); 

/*!
Equivalent to `insert(p)`. 
*/ 
Vertex_handle push_back(const Point& p); 


/*!
Inserts points `a` and `b` in this order, and inserts the line segment `ab` as a 
constraint. Removes the faces crossed by segment `ab` and creates new 
faces instead. If a vertex `c` lies on segment `ab`, constraint `ab` is 
replaced by the two constraints `ac` and `cb`. Apart from the insertion of 
`a` and `b`, the algorithm runs in time proportional to the number of 
removed triangles. 
*/ 
void insert_constraint(Point a, Point b); 

/*!
Equivalent to `insert(c.first, c.second)`. 
*/ 
  void push_back(const std::pair<Point,Point>& c); 

/*!
Inserts the line segment `s` whose endpoints are the vertices 
`va` and 
`vb` as a constraint. The triangles intersected by `s` 
are removed and new ones are created. 
*/ 
void insert_constraint(const Vertex_handle & va, const Vertex_handle & vb); 

/*!
Inserts a polyline defined by the points in the range `[first,last)`.
The polyline is considered as a polygon if the first and last point are equal or if  `close = true`. This enables for example passing the vertex range of a `Polygon_2`.
\tparam PointIterator must be an `InputIterator` with the value type `Point`. 
*/
template < class PointIterator>
void insert_constraint(PointIterator first, PointIterator last, bool close=false);



/*! 
Removes a vertex `v`. 
\pre Vertex `v` is not incident to a constrained edge. 
*/ 
void remove(Vertex_handle v); 

/*!
Make the edges incident to vertex `v` unconstrained edges. 
*/ 
void remove_incident_constraints(Vertex_handle v); 

/*!
Make edge `(f,i)` unconstrained. 
*/ 
void remove_constrained_edge(Face_handle f, int i); 

/*!
Checks the validity of the triangulation and the consistency
of the constrained marks in edges.
*/ 
bool 
is_valid(bool verbose = false, int level = 0) const; 

/// @}

}; /* end Constrained_triangulation_2 */

/*!
Writes the triangulation as for `Triangulation_2<Traits,Tds>` and, for each face `f`, and integers `i=0,1,2`, 
writes "C" or "N" depending whether edge 
`(f,i)` is constrained or not. 
\relates Constrained_triangulation_2 
*/ 
ostream & operator<<(ostream& os, const Constrained_triangulation_2<Traits,Tds> &Ct); 

/*!
Reads a triangulation from stream `is` and assigns it to `t`. Data in the stream must have the same format `operator<<` uses. 
Note that `t` is first cleared. 
\relates Constrained_triangulation_2 
*/ 
istream& operator>>(istream& is,Constrained_triangulation_2<Traits,Tds> Ct& t); 
} /* end namespace CGAL */
