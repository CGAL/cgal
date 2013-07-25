
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Triangulation_simplex_3` stores a simplex of any dimension 
defined by the `Triangulation_3` class. It also defines the 
operator less such that simplices can be stored in a `map` or a 
`set` of simplices. The simplex is invalidated by any change in 
the triangulation. 

\cgalHeading{Parameters}

It is parameterized by the triangulation it derives the simplices 
from. 

\sa `CGAL::Triangulation_3<TriangulationTraits_3,TriangulationDataStructure_3>`

*/
template< typename Triangulation_3 >
class Triangulation_simplex_3 {
public:

/// \name Types 
/// @{

/*!

The simplex class itself. 
*/ 
typedef Triangulation_simplex_3<Triangulation_3> Simplex; 

/*!

*/ 
typedef Triangulation_3::Vertex_handle Vertex_handle; 

/*!

*/ 
typedef Triangulation_3::Edge Edge; 

/*!

*/ 
typedef Triangulation_3::Facet Facet; 

/*!

*/ 
typedef Triangulation_3::Cell_handle Cell_handle; 

/*!

*/ 
typedef Triangulation_3::Cell_circulator Cell_circulator; 

/*!

*/ 
typedef Triangulation_3::Facet_circulator Facet_circulator; 

/*!

*/ 
typedef Triangulation_3::Edge_iterator Edge_iterator; 

/*!

*/ 
typedef Triangulation_3::Facet_iterator Facet_iterator; 

/*!

*/ 
typedef Triangulation_3::Finite_vertices_iterator Finite_vertices_iterator; 

/*!

*/ 
typedef Triangulation_3::Finite_edges_iterator Finite_edges_iterator; 

/*!

*/ 
typedef Triangulation_3::Finite_facets_iterator Finite_facets_iterator; 

/*!

*/ 
typedef Triangulation_3::Finite_cells_iterator Finite_cells_iterator; 

/// @} 

/// \name Creation 
/// @{

/*!
Initializes the simplex to 
an invalid simplex. 
*/ 
Triangulation_simplex_3(); 

/*!

*/ 
Triangulation_simplex_3(Vertex_handle vh); 

/*!

*/ 
Triangulation_simplex_3(Edge e); 

/*!

*/ 
Triangulation_simplex_3(Facet f); 

/*!

*/ 
Triangulation_simplex_3(Cell_handle ch); 

/*!

*/ 
Triangulation_simplex_3(Cell_circulator ccir); 

/*!

*/ 
Triangulation_simplex_3(Facet_circulator fcir); 

/*!

*/ 
Triangulation_simplex_3(Edge_iterator eit); 

/*!

*/ 
Triangulation_simplex_3(Facet_iterator fit); 

/// @} 

/// \name Operations 
/// @{

/*!
returns the dimension of the 
simplex. 
*/ 
int dimension () const; 

/*!
Returns the `Vertex_handle` 
stored in the simplex. \pre dimension() == 0 
*/ 
operator Vertex_handle () const; 

/*!
Returns the `Edge` 
stored in the simplex. \pre dimension() == 1 
*/ 
operator Edge () const; 

/*!
Returns the `Facet` 
stored in the simplex. \pre dimension() == 2 
*/ 
operator Facet () const; 

/*!
Returns the `Cell_handle` 
stored in the simplex. \pre dimension() == 3 
*/ 
operator Cell_handle () const; 

/*!
Returns a cell incident 
to the simplex. 
*/ 
Cell_handle incident_cell () const; 

/*!
Test whether two 
simplices are equal. 
*/ 
bool operator==(const 
Triangulation_simplex_3<Triangulation_3> &s1); 

/*!
Defines a ordering 
on the simplices. This ordering depends on the memory layout and is 
independent of the geometry. Therefore, the ordering is not intrinsic 
*/ 
bool operator< (const 
Triangulation_simplex_3<Triangulation_3> &s1); 

/// @}

}; /* end Triangulation_simplex_3 */
} /* end namespace CGAL */
