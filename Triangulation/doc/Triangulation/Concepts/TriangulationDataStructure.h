
/*!
\ingroup PkgTriangulationsConcepts
\cgalConcept

The `TriangulationDataStructure` concept describes objects responsible for
storing and maintaining the combinatorial part of a
\f$ d\f$-dimensional pure simplicial complex that has the topology
of the \f$ d\f$-dimensional sphere \f$ \mathbb{S}^d\f$ with \f$ d\in[-2,D]\f$.
Since the simplicial \f$ d\f$-complex is pure, all
faces are sub-faces of some \f$ d\f$-simplex. And since it has the
topology of the sphere \f$ \mathbb{S}^d\f$, it is manifold, thus
any \f$ d-1\f$-face belongs to exactly two \f$ d\f$-dimensional full cells.

The concept `TriangulationDataStructure` includes two sub-concepts
`TriangulationDataStructure::Vertex` and
`TriangulationDataStructure::FullCell`.

Possible values for the current dimension \f$ d\f$ include

<DL>
<DT><B>-2</B><DD> This corresponds to the non-existence of any object in
the triangulation.
<DT><B>-1</B><DD> This corresponds to a single vertex and a single full cell,
which is also the unique vertex and the unique full cell in the
`TriangulationDataStructure`.
In a
geometric realization of the `TriangulationDataStructure` (<I>e.g.</I>, in a
`Triangulation<TriangulationTraits_, TriangulationDataStructure_>` or a
`Delaunay_triangulation<DelaunayTriangulationTraits_, TriangulationDataStructure_>`), this vertex
corresponds to <I>the vertex at infinity</I>.

<DT><B>0</B><DD> This corresponds to two vertices, each incident to one \f$ 0\f$-face;
the two full cells being neighbors of each other. This is the unique
triangulation of the \f$ 0\f$-sphere.

<DT><B>\f$ d>0\f$</B><DD> This corresponds to a triangulation of the sphere
\f$ \mathbb{S}^d\f$.
</DL>

An \f$ i\f$-simplex is a simplex with \f$ i+1\f$ vertices. An \f$ i\f$-simplex \f$ \sigma\f$ is
incident to a \f$ j\f$-simplex \f$ \sigma'\f$, \f$ j<i\f$, if and only if \f$ \sigma'\f$
is a proper face of \f$ \sigma\f$.

We call a \f$ 0\f$-simplex a <I>vertex</I>, a \f$ (d-1)\f$-simplex a <I>facet</I> and a
\f$ d\f$-simplex a <I>full cell</I>. A <I>face</I> can have any dimension.
Two full cells are <I>neighbors</I> if they share a facet. Two faces are
<I>incident</I> if one is included in the other.

\cgalHeading{Input/Output}

The information stored in the `iostream` is:

- the current dimension (which must be smaller than or equal to `tds.maximal_dimension()`),

- the number of vertices,

- for each vertex the information of that vertex,

- the number of full cells,

- for each full cell the indices of its vertices and extra information for that full cell,

- for each full cell the indices of its neighbors.

The indices of vertices and full cells correspond to the order in the
file; the user cannot control it.
The classes `Vertex` and
`Full_cell` have to provide the relevant I/O operators
(possibly empty).

\cgalHasModelsBegin
\cgalHasModels{CGAL::Triangulation_data_structure<Dimensionality, TriangulationDSVertex_, TriangulationDSFullCell_>}
\cgalHasModelsEnd

\sa `TriangulationDataStructure::Vertex`
\sa `TriangulationDataStructure::FullCell`
\sa `TriangulationDSVertex`
\sa `TriangulationDSFullCell`
\sa `TriangulationDSFace`

*/

class TriangulationDataStructure {
public:

/// \name Types
/// @{

/*!
The vertex type, requirements for this type are described
in the concept `TriangulationDataStructure::Vertex`.
*/
typedef unspecified_type Vertex;

/*!
The full cell type, requirements for this type are described
in the concept `TriangulationDataStructure::FullCell`.
*/
typedef unspecified_type Full_cell;

/*!
A model of the concept `FullCellData`.
*/
typedef unspecified_type Full_cell_data;

/*!
The facet type, for describing faces of the triangulation with codimension 1.
The constructor `Facet(c,i)` constructs a `Facet` representing the facet of
full cell `c` opposite to its `i`-th vertex. Its dimension is
`current_dimension()-1`.
*/
typedef unspecified_type Facet;

/*!
The face type, which must be a model of the concept `TriangulationDSFace`.
*/
typedef unspecified_type Face;

/// @}


/// \name Handles
/// Vertices and full cells are manipulated via
/// <I>handles</I>. Handles support the usual two dereference
/// operators and `operator->`.
/// @{

/*!
Handle to a `Vertex`.
*/
typedef unspecified_type Vertex_handle;

/*!
Const handle to a `Vertex`.
*/
typedef unspecified_type Vertex_const_handle;

/*!
Handle to a `Full_cell`.
*/
typedef unspecified_type Full_cell_handle;

/*!
Const handle to a `Full_cell`.
*/
typedef unspecified_type Full_cell_const_handle;

/// @}

/// \name Iterators
/// Vertices, facets and full cells can be iterated over using
/// <I>iterators</I>. Iterators support the usual two dereference
/// operators and `operator->`.
/// @{

/*!
Iterator over the list of vertices.
*/
typedef unspecified_type Vertex_iterator;

/*!
Iterator over the list of full cells.

*/
typedef unspecified_type Full_cell_iterator;

/*!
Iterator over the facets of the complex.

*/
typedef unspecified_type Facet_iterator;

/*!
Size type (an unsigned integral type).
*/
typedef unspecified_type size_type;

/*!
Difference type (a signed integral type).
*/
typedef unspecified_type difference_type;

/// @}

/// \name Creation
/// @{

/*!
Creates an instance `tds` of
type `TriangulationDataStructure`. The maximal dimension of its full cells is `dim` and
`tds` is initialized to the empty triangulation. Thus,
`tds`.`current_dimension()` equals `-2`.
The parameter `dim` can be ignored by the constructor if it is already
known at compile-time. Otherwise, the following precondition holds:
\pre `dim>0`.
*/
TriangulationDataStructure(int dim = 0);

/// @}

/// \name Queries
/// @{

/*!
Returns the maximal dimension of
the full dimensional cells that can be stored in the triangulation `tds`.
\post the
returned value is positive.
*/
int maximal_dimension() const;

/*!
Returns the dimension of the
full dimensional cells stored in the triangulation. It holds that
`tds`.`current_dimension()=-2` if and only if `tds`.`empty()` is
`true`. \post the returned value `d` satisfies
\f$ -2\leq d \leq\f$`tds`.`maximal_dimension()`.
*/
int current_dimension() const;

/*!
Returns `true` if the triangulation
contains nothing, returns `false` otherwise.
*/
bool empty() const;

/*!
Returns the number of vertices in the triangulation.
*/
size_type number_of_vertices() const;

/*!
Returns the number of full cells in the triangulation.
*/
size_type number_of_full_cells() const;

/*!
Tests whether `v` is a vertex of the triangulation.
*/
bool is_vertex(const Vertex_handle & v) const;

/*!
Tests whether `c` is a full cell of the triangulation.
*/
bool is_full_cell(const Full_cell_handle & c) const;

/*!
This function computes (<I>gathers</I>) a connected set of full cells
satisfying a common criterion. Call them <I>good</I> full cells. It is assumed
that the argument `start` is a good full cell. The full cells are then
recursively explored by examining if, from a given good full cell, its adjacent
full cells are also good.

The argument `tp` is a predicate, i.e.\ a function or a functor providing
`operator()`, that takes as argument a `Facet`
whose `Full_cell` is good.
The predicate must return `true`
if the traversal of that `Facet` leads to a good full cell.

All the good full cells are output into the last argument `out`.

Returns a facet on the boundary of the set of cells.

\pre `start != Full_cell_handle()` and `start` is a good cell.
*/
template< typename TraversalPredicate, typename OutputIterator >
Facet gather_full_cells(Full_cell_handle start, TraversalPredicate & tp,
OutputIterator & out) const;

/*!
Inserts in `out` all the full cells that are incident to the vertex
`v`, i.e., the full cells that have the `Vertex v` as a vertex.
Returns the output iterator.
\pre `v!=Vertex_handle()`.

*/
template< typename OutputIterator > OutputIterator
incident_full_cells(Vertex_handle v, OutputIterator out) const;

/*!
Inserts in `out` all the full cells that are incident to the face `f`,
i.e., the full cells that have the `Face f` as a subface.
Returns the output iterator.
\pre `f.full_cell()!=Full_cell_handle()`.

*/
template< typename OutputIterator > OutputIterator
incident_full_cells(const Face & f, OutputIterator out) const;

/*!
Inserts in `out` all the full cells that share at least one vertex with the `Face f`. Returns the output iterator.
\pre `f.full_cell()!=Full_cell_handle()`.

*/
template< typename OutputIterator > OutputIterator
star(const Face & f, OutputIterator out) const;

/*!
Constructs all the `Face`s of dimension `dim` incident to
`Vertex` v and inserts them in the `OutputIterator out`. If `dim >=` `tds`.`current_dimension()`, then no `Face` is
constructed.
\pre `0 < dim` and `v!=Vertex_handle()`.

*/
template< typename OutputIterator > OutputIterator
incident_faces(Vertex_handle v, const int dim, OutputIterator
out);

/// @}

/// \name Accessing the Vertices
/// @{

/*!
Returns a handle to the `i`-th `Vertex` of the `Full_cell` `c`.
\pre \f$0 \leq i \leq \f$`tds`.`current_dimension()` and `c!=Full_cell_handle()`.
*/
Vertex_handle vertex(Full_cell_handle c, const int i) const;

/*!
Returns the index of `c` as a neighbor of its `i`-th neighbor.
\pre \f$0 \leq i \leq \f$`tds`.`current_dimension`() and `c!=Full_cell_handle()`.
*/
int mirror_index(Full_cell_handle c, int i) const;

/*!
Returns the vertex of the `i`-th neighbor of `c`
that is not vertex of `c`.
\pre \f$0 \leq i \leq \f$`tds`.`current_dimension`() and `c!=Full_cell_handle()`.
*/
Vertex_handle mirror_vertex(Full_cell_handle c, int i) const;

/*!
Iterator to the first vertex of `tds`. User has no control on the order.
*/
Vertex_iterator vertices_begin();

/*!
Iterator referring beyond the last vertex of `tds`.
*/
Vertex_iterator vertices_end();

/// @}

/// \name Accessing the Full Cells
/// @{

/*!
Returns a full cell incident to `Vertex` `v`. Note that this
full cell is
not unique (`v` is typically incident to more than one full cell).
\pre `v != Vertex_handle`
*/
Full_cell_handle full_cell(Vertex_handle v) const;

/*!
Returns a `Full_cell_handle` pointing to the `Full_cell`
opposite to the `i`-th vertex of `c`.
\pre \f$0 \leq i \leq \f$`tds`.`current_dimension()` and `c != Full_cell_handle()`
*/
Full_cell_handle neighbor(Full_cell_handle c, int i) const;

/*!
Iterator to the first full cell of `tds`. User has no control on the order.
*/
Full_cell_iterator full_cells_begin();

/*!
Iterator referring beyond the last full cell of `tds`.
*/
Full_cell_iterator full_cells_end();

/// @}

/// \name Faces and Facets
/// @{

/*!
Iterator to the first facet of the triangulation.
*/
Facet_iterator facets_begin();

/*!
Iterator referring beyond the last facet of the triangulation.
*/
Facet_iterator facets_end();

/*!
Returns a full cell containing the facet `f`
*/
Full_cell_handle full_cell(const Facet & f) const;

/*!
Returns the index of vertex of the full cell `c=tds.full_cell(f)`
which does not belong to `c`.
*/
int index_of_covertex(const Facet & f) const;

/// @}

/// \name Vertex Insertion
/// @{

/*!
Inserts a new
vertex `v` in the full cell `c` and returns a handle to
it. The full cell
`c` is subdivided into `tds`.`current_dimension()`+1 full cells which
share the vertex `v`.

\cgalFigureBegin{triangulationfiginsertfullcell,insert-in-cell.png}
Insertion in a full cell, \f$ d=2\f$
\cgalFigureEnd

\pre Current dimension is positive and `c` is a full cell of `tds`.
*/
Vertex_handle insert_in_full_cell(Full_cell_handle c);

/*!
Inserts a vertex in the triangulation data structure by subdividing the
`Face f`. Returns a handle to the newly created `Vertex`.

\cgalFigureBegin{triangulationfiginsertface,insert-in-face.png}
Insertion in face, \f$ d=3\f$
\cgalFigureEnd
*/
Vertex_handle insert_in_face(const Face & f);

/*!
Inserts a vertex in the triangulation data structure by subdividing the
`Facet ft`. Returns a handle to the newly created `Vertex`.
*/
Vertex_handle insert_in_facet(const Facet & ft);

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Removes the full cells in the range \f$ C=\f$`[s, e)`, inserts a vertex
at position `p` and fills the hole by connecting
each face of the boundary to `p`.
A `Vertex_handle` to the new `Vertex` is
returned. The facet `ft` must lie on the boundary of \f$ C\f$ and its
defining full cell, `tr`.`full_cell(ft)` must lie inside \f$ C\f$. Handles
to the newly created full cells are output in the `out` output iterator.
\pre `c` belongs to \f$ C\f$ and `c->neighbor(i)`
does not, with `f=(c,i)`.
\f$ H\f$, the union of full cells in \f$ C\f$, is simply connected and its
boundary \f$ \partial H\f$ is a
combinatorial triangulation of the sphere \f$ \mathbb{S}^{d-1}\f$.
All vertices of cells of \f$ C\f$ are on \f$ \partial H\f$.

\cgalFigureBegin{triangulationfiginserthole,insert-in-hole.png}
Insertion in a hole, \f$ d=2\f$
\cgalFigureEnd
\cgalAdvancedEnd
*/
template< class ForwardIterator > Vertex_handle
insert_in_hole(ForwardIterator start, ForwardIterator end, Facet f);

/*!
Same as above, but handles to the new full cells are
appended to the `out` output iterator.
*/
template< class ForwardIterator, class OutputIterator >
Vertex_handle insert_in_hole(ForwardIterator start, ForwardIterator end, Facet
f, OutputIterator out);

/*!
Transforms a triangulation of the sphere \f$ \mathbb{S}^d\f$ into the
triangulation of the sphere \f$ \mathbb{S}^{d+1}\f$ by adding a new vertex
`v`.
`v` is used to triangulate one of the two half-spheres of
\f$ \mathbb{S}^{d+1}\f$ (\f$ v\f$ is added as \f$ (d+2)^{th}\f$ vertex to all
full cells)
and `star` is used to triangulate the other half-sphere
(all full cells that do not already have star as vertex are duplicated,
and `star` replaces `v` in these full cells).
The indexing of the vertices in the
full cell is such that, if `f` was a full cell of maximal dimension in the
initial complex, then `(f,v)`, in this order, is the corresponding full cell
in the updated triangulation. A handle to `v` is returned
(see Figure \cgalFigureRef{triangulationfiginsertincreasedim}).
\pre
If the current dimension is -2 (empty triangulation), then `star`
can be omitted (it is ignored), otherwise
the current dimension must be strictly less than the maximal dimension
and `star` must be a vertex of `tds`.

\cgalFigureBegin{triangulationfiginsertincreasedim,insert-increase-dim.png}
Insertion, increasing the dimension from \f$ d=1\f$ to \f$ d=2\f$
\cgalFigureEnd

*/
Vertex_handle insert_increase_dimension(Vertex_handle star);

/*!
Adds a new full cell to `tds` and
returns a handle to it. The new full cell has no vertices and no neighbors yet.
*/
Full_cell_handle new_full_cell();

/*!
Adds a new vertex to `tds` and returns a handle to it. The new vertex has
no associated full cell nor index yet.
*/
Vertex_handle new_vertex();

/*!
Sets the `i`-th vertex of `c` to `v` and, if `v != Vertex_handle()`,
sets `c` as the incident full cell of `v`.
*/
void associate_vertex_with_full_cell(Full_cell_handle c, int i,
Vertex_handle v);

/*!
Sets the neighbor opposite to vertex `i` of `Full_cell` `ci` to
`cj`. Sets the neighbor opposite to vertex `j` of `Full_cell`
`cj` to `ci`.
*/
void set_neighbors(Full_cell_handle ci, int i, Full_cell_handle cj, int
j);

/*!
Forces the current dimension
of the complex to `d`.
\pre \f$-2 \leq d \leq \f$`maximal_dimension()`.
*/
void set_current_dimension(int d);

/// @}

/// \name Vertex Removal
/// @{

/*!
Reinitializes `tds` to the empty complex.
*/
void clear();

/*!
Contracts the
`Face f` to a single vertex. Returns a handle to that vertex
(see Figure \cgalFigureRef{triangulationfigcollapseface}).
\pre The boundary of the full cells incident to `f`
is a topological sphere of dimension
`tds`.`current_dimension()`-1).

\cgalFigureBegin{triangulationfigcollapseface,collapse-face.png}
Collapsing an edge in dimension \f$ d=3\f$, `v` is returned.
\cgalFigureEnd


*/
Vertex_handle collapse_face(const Face & f);

/*!
This method does exactly the opposite of
`insert_increase_dimension()`:
`v` is removed,
full cells not containing `star` are removed,
full cells containing `star` but not `v` lose vertex `star`,
full cells containing `star` and `v` lose vertex `v`
(see Figure \cgalFigureRef{triangulationfiginsertincreasedim}).
\pre All cells contain either `star` or `v`.
Edge `star-v` exists in the triangulation
and `current_dimension() != -2`.

*/
void remove_decrease_dimension(Vertex_handle v, Vertex_handle star);

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Remove the vertex `v` from the triangulation.
\cgalAdvancedEnd
\pre `v` is a vertex of `tds`.
*/
void delete_vertex(Vertex_handle v);

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Remove the full cell `c` from the triangulation.
\cgalAdvancedEnd
\pre `c` is a full cell of `tds`.
*/
void delete_full_cell(Full_cell_handle c);

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Calls `delete_full_cell` over an iterator range of value type `Full_cell_handle`.
\cgalAdvancedEnd
*/
template< typename ForwardIterator > void
delete_full_cells(ForwardIterator start, ForwardIterator end);

/// @}

/// \name Validity Check
/// @{

/*!
\cgalDebugFunction
\cgalDebugBegin
Partially checks whether `tds` is indeed a triangulation.

It must <I>at least</I>
<UL>
<LI>check the validity of the vertices and full cells of `tds` by calling
their respective `is_valid` method.
<LI>check that each full cell has no duplicate vertices and has as many
neighbors as its number of facets (`current_dimension()+1`).
<LI>check that each full cell share exactly `tds`.`current_dimension()`
vertices with each of its neighbor.
</UL>

When `verbose` is set to `true`, messages are printed to give
a precise indication on the kind of invalidity encountered.

Returns `true` if all the tests pass, `false` if any test fails. See
the documentation for the models of this concept to see the additional (if
any) validity checks that they implement.
\cgalDebugEnd
*/
bool is_valid(bool verbose=false) const;

/*!
Reads a combinatorial triangulation from `is` and assigns it to
`tds`. \pre The dimension of the input complex must be less than or
equal to `tds`.`maximal_dimension()`.
*/
std::istream & operator>>(std::istream & is, TriangulationDataStructure &
tds);

/*!
Writes `tds` into the output stream `os`.
*/
std::ostream & operator<<(std::ostream & os, const TriangulationDataStructure
& tds);

/// @}

}; /* end TriangulationDataStructure */


//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================


/*!
\ingroup PkgTriangulationsConcepts
\cgalConcept

The concept `TriangulationDataStructure::Vertex` describes the type used by a
`TriangulationDataStructure` to store the vertices.

It sets requirements of combinatorial nature
only, as geometry is not concerned here. In particular, we only require that
the vertex holds a handle to a full cell incident to it in the triangulation.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Triangulation_ds_vertex<TriangulationDataStructure_>}
\cgalHasModels{CGAL::Triangulation_vertex<TriangulationTraits_, Data, TriangulationDSVertex_>}
\cgalHasModelsEnd

\sa `TriangulationDataStructure::FullCell`
\sa `TriangulationDataStructure::Face`
\sa `TriangulationDataStructure`
*/
class TriangulationDataStructure::Vertex {
public:

/// \name Types
/// @{

/*!
A handle to a cell, which must be the same as the
nested type `TriangulationDataStructure::Full_cell_handle`.
*/
typedef unspecified_type Full_cell_handle;

/// @}

/// \name Operations
/// @{

/*!
Set `c` as the vertex's
incident full cell. \pre `c` must not be the default-constructed
`Full_cell_handle`.
*/
void set_full_cell(Full_cell_handle c);

/*!
Returns a handle to a
full cell incident to the vertex.
*/
Full_cell_handle full_cell() const;

/// @}

/// \name Validity Check
/// @{

/*!
\cgalDebugFunction
\cgalDebugBegin
Performs some validity checks on the vertex `v`.

It must <I>at least</I> check that `v` has an incident full cell, which in
turn must contain `v` as one of its vertices.

Returns `true` if all the tests pass, `false` if any test fails. See
the documentation for the models of this concept to see the additional (if
any) validity checks that they implement.
\cgalDebugEnd
*/
bool is_valid(bool verbose=false) const;

/// @}

/// \name Input/Output
/// These operators can be used directly and are called by the I/O
/// operator of class `TriangulationDataStructure`.
///
/// @{
/*!
Writes (possibly) non-combinatorial information about vertex `v` to the stream
`os`.
*/
template<class TriangulationDataStructure>
std::ostream& operator<<(std::ostream & os, const Triangulation_ds_vertex<TriangulationDataStructure> & v);

/*!
Reads from stream `is` the vertex information written by `operator<<`.
*/
template<class TriangulationDataStructure>
std::istream& operator>>(std::istream & is, Triangulation_ds_vertex<TriangulationDataStructure> & v);

/// @}

}; /* end TriangulationDataStructure::Vertex */


//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================


/*!
\ingroup PkgTriangulationsConcepts
\cgalConcept

The concept `TriangulationDataStructure::FullCell` describes the type used by a
`TriangulationDataStructure` to store the full cells.

It sets requirements of combinatorial nature
only, as geometry is not concerned here.
In the context of triangulation, the term full cell refers to a face of
<I>maximal</I> dimension. This maximality characteristic is emphasized by using
the adjective <I>full</I>.

A `TriangulationDataStructure::FullCell` is responsible for
storing handles to the vertices of the
full cell as well as handles to the adjacent full cells. Two full cells
are said to be adjacent when they share a facet. Adjacent full cells are
called hereafter neighbors.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Triangulation_ds_full_cell<TriangulationDataStructure_, DSFullCellStoragePolicy>}
\cgalHasModels{CGAL::Triangulation_full_cell<TriangulationTraits_, Data, TriangulationDSFullCell_>}
\cgalHasModelsEnd

\sa `TriangulationDataStructure::FullCell`
\sa `TriangulationDataStructure::Face`
\sa `TriangulationDataStructure`
*/
class TriangulationDataStructure::FullCell {
public:

/// \name Types
/// @{

/*!
A handle to a vertex, which must be the same as the
nested type `TriangulationDataStructure::Vertex_handle`.
*/
typedef unspecified_type Vertex_handle;

/*!
An iterator over the handles to
the vertices of the full cell.
*/
typedef unspecified_type Vertex_handle_iterator;

/*!
A handle to a full cell, which must be the same as the
nested type `TriangulationDataStructure::Full_cell_handle`.
*/
typedef unspecified_type Full_cell_handle;

/*!
A data member of this type has to be stored and accessible through
access function below.
*/
typedef TriangulationDataStructure::Full_cell_data
TDS_data;

/// @}

/// \name Access Functions
/// @{

/*!
Returns one less than the maximum
number of vertices that the full cell can store. This does not return
the dimension of the actual full cell stored in `c`.
*/
int maximal_dimension() const;

/*!
Returns an iterator to the first `Vertex_handle` stored in the
full cell.
*/
Vertex_handle_iterator vertices_begin() const;

/*!
Returns an iterator pointing beyond the last `Vertex_handle` stored in
the full cell.
*/
Vertex_handle_iterator vertices_end() const;

/*!
Returns the `i`-th vertex
of the full cell. \pre \f$0 \leq i \leq \f$ `maximal_dimension()`.
*/
Vertex_handle vertex(const int i) const;

/*!
Returns the
full cell opposite to the `i`-th vertex of the full cell `c`. \pre \f$0 \leq i \leq \f$`maximal_dimension()`.
*/
Full_cell_handle neighbor(const int i) const;

/*!
Returns the index of `c` in its \f$ i^{th}\f$ neighbor (`c.neighbor(i)`).
If the returned integer is not negative,
it holds that `c.neighbor(i)->neighbor(j) == c`. Returns
`-1` if `c` has no adjacent full cell of index `i`.
\pre \f$0 \leq i \leq \f$ `maximal_dimension()`.
*/
int mirror_index(const int i) const;

/*!
Returns the index `i`
such that `c.neighbor(i)==n`.
\pre `n` must be a neighbor of `c`.
*/
int index(Full_cell_handle n) const;

/*!
Returns the index `i` of
the vertex `v` such that `c.vertex(i)==v`. \pre `v` must be
a vertex of the `c`.
*/
int index(Vertex_handle v) const;


/// \name Internal
/// \cgalAdvancedBegin
/// These functions are used internally by the triangulation data
/// structure. The user is not encouraged to use them directly as they
/// may change in the future.
/// \cgalAdvancedEnd
/// @{

/*!
Returns the data member of
type `TDS_data`. It is typically used to mark the full cell as <I>visited</I>
during operations on a `TriangulationDataStructure`.
*/
const TDS_data & tds_data() const;

/*!
Same as above, but returns a reference to
a non-`const` object.
*/
TDS_data & tds_data();

/*!
Returns a handle to the mirror vertex of the `i`-th vertex of full cell `c`.
`cur_dim` is the current dimension of the triangulation data structure.
\cgalAdvancedBegin
This function works even if the adjacency information stored in the
neighbor full cell `*c.neighbor(i)` is corrupted. This is useful
when temporary corruption is necessary during surgical operations on a
triangulation.
\cgalAdvancedEnd

\pre \f$0 \leq i,\f$ `cur_dim` \f$ \leq \f$ `maximal_dimension()`.
*/
Vertex_handle mirror_vertex(const int i, const int cur_dim) const;

/// @}

/// \name Update Functions
/// @{

/*!
Sets the \f$ i\f$-th
vertex of the full cell.
\pre \f$0 \leq i \leq \f$ `maximal_dimension()`.
*/
void set_vertex(const int i, Vertex_handle v);

/*!
Sets the
`i`-th neighbor of `c` to `n`. Full cell `n` is
opposite to the \f$ i\f$-th vertex of `c`.
\pre \f$0 \leq i \leq \f$`maximal_dimension()`.
*/
void set_neighbor(const int i, Full_cell_handle n);

/*!
Sets the
mirror index of the \f$ i\f$-th vertex of `c` to `index`. This corresponds
to the index, in `c->neighbor(i)`, of the full cell `c`.

Note: a model of this concept may choose not to store mirror
indices, in which case this function should do nothing.
\pre \f$0 \leq i \leq \f$`maximal_dimension()`.
*/
void set_mirror_index(const int i, const int index);

/*!
Switches the orientation of the
full cell `c` by swapping its vertices with index `d1` and `d2`.
\pre \f$0 \leq d1,d2 \leq \f$`maximal_dimension()`.
*/
void swap_vertices(int d1, int d2);

/// @}

/// \name Queries
/// @{

/*!
Returns `true`
if the vertex `v` is a vertex of the full cell `c`. Returns `false`
otherwise.
*/
bool has_vertex(Vertex_handle v) const;

/*!
Returns `true` and sets the value of `ret` to the index of `v` in
`c` if the vertex `v` is a vertex of the full cell `c`. Returns
`false` otherwise.
*/
bool has_vertex(Vertex_handle v, int & ret) const;

/*!
Returns `true`
if the full cell `n` is a neighbor of the full cell `c`. Returns
`false` otherwise.
*/
bool has_neighbor(Full_cell_handle n) const;

/*!
Returns `true` and sets the value of `ret` to the index of `n` as
a neighbor of `c` if the full cell `n` is a neighbor of the full cell
`c`. Returns `false` otherwise.
*/
bool has_neighbor(Full_cell_handle n, int & ret) const;

/// @}

/// \name Validity Check
/// @{

/*!
\cgalDebugFunction
\cgalDebugBegin
Performs some validity checks on the full cell `c`.

It must <I>at least</I> check that for each <I>existing</I> neighbor `n`,
`c` is also a neighbor of `n`.

Returns `true` if all the tests pass, `false` if any test fails. See
the documentation for the models of this concept to see the additional (if
any) validity checks that they implement.
\cgalDebugEnd
*/
bool is_valid(bool verbose=false) const;

/// @}


/// \name Input/Output
/// These operators can be used directly and are called by the I/O
/// operator of class `TriangulationDataStructure`.
///
/// @{

/*!
Writes (possibly) non-combinatorial information about full cell `c` to the stream
`os`.
*/
template<class TriangulationDataStructure>
std::ostream& operator<<(std::ostream & os, const Triangulation_ds_full_cell<TriangulationDataStructure> & c);

/*!
Reads from stream `is` the full cell information written
by `operator<<`.
*/
template<class TriangulationDataStructure>
std::istream& operator>>(std::istream & is, Triangulation_ds_full_cell<TriangulationDataStructure> & c);

/// @}

}; /* end TriangulationDataStructure::Vertex */
