
namespace CGAL {

/*!
\ingroup PkgTriangulation2TriangulationClasses

A constrained triangulation is a triangulation of a set of points 
which has to include among its edges 
a given set of segments 
joining the points. The given segments are 
called <I>constraints</I> and the corresponding 
edges in the triangulation are called <I>constrained edges</I>. 

The endpoints of constrained edges are of course vertices of the 
triangulation. However the triangulation may include 
other vertices as well. 
There are three versions of constrained triangulations 
<UL> 
<LI>In the basic version, the constrained triangulation 
does not handle intersecting constraints, and the set of input 
constraints is required to be a set of segments that do not intersect 
except possibly at their endpoints. Any number of constrained edges 
are allowed to share the same endpoint. Vertical constrained edges 
are allowed as well as 
constrained edges with null length. 
<LI>The two other versions support intersecting input constraints. 
In those versions, input constraints are allowed to be 
intersecting, overlapping or partially 
overlapping segments. 
The triangulation introduce additional vertices at each point which 
is a proper intersection point of two 
constraints. A single constraint intersecting other 
constraints will then appear as the union of several 
constrained edges of the triangulation. 
The two versions dealing with intersecting constraints, slightly differ 
in the way intersecting constraints are dealt with. 
<UL> 
<LI>One of them is 
designed to be robust when predicates are evaluated exactly but 
constructions (i. e. intersection computations) are 
approximate. 
<LI>The other one is designed to be used 
with an exact arithmetic (meaning exact 
evaluation of predicates and exact computation of intersections.) 
This last version finds its full efficiency when used in conjunction 
with a constraint hierarchy data structure 
as provided in the class 
`Constrained_triangulation_plus_2`. See 
section \ref Section_2D_Triangulations_Constrained_Plus. 
</UL> 
</UL> 

\image html constraints.gif

The class `Constrained_triangulation_2` of the \cgal library 
implements constrained triangulations. 
The template parameter `Traits` 
stands for a geometric traits class. It has to be a model 
of the concept `TriangulationTraits_2`. 
When intersection of input constraints are supported, 
the geometric traits class 
is required to provide additional function object types 
to compute the intersection of two segments. 
It has then to be a model of the concept 
`ConstrainedTriangulationTraits_2`. 
The template parameter `Tds` 
stands for 
a triangulation data structure class that has to be a model 
of the concept `TriangulationDataStructure_2`. 
The third parameter `Itag` is the intersection tag 
which serves to choose between the different 
strategies to deal with constraints intersections. 
\cgal provides three valid types for this parameter : 

`CGAL::No_intersection_tag` disallows intersections of 
input constraints, 

`CGAL::Exact_predicates_tag` is to be used when the traits 
class 
provides exact predicates but approximate constructions of the 
intersection points. 

`CGAL::Exact_intersections_tag` is to be used in conjunction 
with an exact arithmetic type. 

The information about constrained edges is stored in the 
faces of the triangulation. Thus the nested `Face` 
type of a constrained triangulation offers 
additional functionalities to deal with this information. 
These additional functionalities 
induce additional requirements on the base face class 
plugged into the triangulation data structure of 
a constrained Delaunay triangulation. 
The base face of a constrained Delaunay triangulation 
has to be a model of the concept 
`ConstrainedTriangulationFaceBase_2`. 

\cgal provides default instantiations for the template parameters 
`Tds` and `Itag`, and for the `ConstrainedTriangulationFaceBase_2`. 
If `Gt` is the geometric traits 
parameter, 
the default for 
`ConstrainedTriangulationFaceBase_2` is the class 
`CGAL::Constrained_triangulation_face_base_2<Gt>` 
and the default for the 
triangulation data structure parameter is the class 
`CGAL::Triangulation_data_structure_2 < CGAL::Triangulation_vertex_base_2<Gt>, 		 CGAL::Constrained_triangulation_face_base_2<Gt> >`. 
The default intersection tag is `CGAL::No_intersection_tag`. 

\sa `CGAL::Triangulation_2<Traits,Tds>`
\sa `TriangulationDataStructure_2`
\sa `TriangulationTraits_2` 
\sa `ConstrainedTriangulationTraits_2` 
\sa `ConstrainedTriangulationFaceBase_2` 

### Implementation ###

The insertion of a constrained edge runs in time 
proportional to the number of triangles intersected by this edge. 

*/
template< typename Traits, typename Tds, typename Itag >
class Constrained_triangulation_2 : public Triangulation_2<Traits,Tds> {
public:

/// \name Types 
/// @{

/*! 
The type of input 
constraints 
*/ 
typedef std::pair<Point,Point> Constraint; 

/*! 
The intersection tag which decides how 
intersections between input constraints are dealt with. 
*/ 
typedef Itag Intersection_tag; 

/// @} 

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
Constrained_triangulation_2(); 

/*! 
Copy constructor, all faces and vertices 
are duplicated and the constrained status of edges 
is copied. 
*/ 
Constrained_triangulation_2(const 
Constrained_triangulation_2& ct1); 

/*! 
A templated constructor which introduces and builds 
a constrained triangulation with constrained edges in the range 
\f$ \left[\right.\f$`first`, `last`\f$\left.\right)\f$. 
\pre The `value_type` of `first` and `last` is `Constraint`. 
*/ 
template<class InputIterator> Constrained_triangulation_2( 
InputIterator first, 
InputIterator last, 
const Traits& t=Traits()); 

/// @} 

/// \name Queries 
/// @{

/*! 
Returns true if edge `e` is a constrained edge. 
*/ 
bool is_constrained(Edge e); 

/*! 
Returns true if at least one of the edges incident to vertex `v` 
is constrained. 
*/ 
bool are_there_incident_constraints(Vertex_handle v); 

/*! 
`OutputItEdges` is an output iterator with `Edge` as value 
type. 
Outputs the constrained edges incident to `v` 
in the sequence pointed to by `out` and returns the resulting 
output iterator. 
*/ 
template<class OutputItEdges> 
OutputItEdges incident_constraints(Vertex_handle v, 
OutputItEdges out) const; 

/// @} 

/// \name Insertion and Removal 
/// @{

/*! 
Inserts point `p` and restores the status (constrained or not) of all 
the touched edges. If present `f` is used as an hint 
for the location of `p`. 
*/ 
Vertex_handle insert(Point p, Face_handle f = Face_handle() ); 

/*! 
Same as above except that the location of the point 
`p` to be inserted is assumed to be given by 
`(lt,loc,i)`. 
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
Inserts the points in the range 
\f$ \left[\right.\f$`first`, `last`\f$\left.\right)\f$. 
Returns the number of inserted points. 
\pre The `value_type` of `first` and `last` is `Point`. 
*/ 
template < class InputIterator > 
std::ptrdiff_t 
insert(InputIterator first, InputIterator last); 

/*! 
Inserts points `a` and `b`, and inserts segment `ab` as a 
constraint. Removes the faces crossed by segment `ab` and creates new 
faces instead. If a vertex `c` lies on segment `ab`, constraint `ab` is 
replaced by the two constraints `ac` and `cb`. Apart from the insertion of 
`a` and `b`, the algorithm runs in time proportional to the number of 
removed triangles. 
\pre The relative interior of segment `ab` does not intersect the relative interior of another constrained edge. 
*/ 
void insert_constraint(Point a, Point b); 

/*! 
Inserts constraints `c` as above. 
*/ 
void push_back(const Constraint& c); 

/*! 
Inserts the line segment `s` whose endpoints are the vertices 
`va` and 
`vb` as a constrained edge `e`. The triangles intersected by `s` 
are removed and new ones are created. 
*/ 
void insert_constraint(const Vertex_handle & va, const Vertex_handle & vb); 

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
Make edge `(f,i)` no longer constrained. 
*/ 
void remove_constrained_edge(Face_handle f, int i); 

/*! 
\advanced Checks the validity of the triangulation and the consistency
of the constrained marks in edges.
*/ 
bool 
is_valid(bool verbose = false, int level = 0) const; 

/// @}

}; /* end Constrained_triangulation_2 */

/*! 
Writes the triangulation as for `CGAL::Triangulation_2<Traits,Tds>` and, for each face f, and integers i=0,1,2, 
write "C" or "N" depending whether edge 
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
