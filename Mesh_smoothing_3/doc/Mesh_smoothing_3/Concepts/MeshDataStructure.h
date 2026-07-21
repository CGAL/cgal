/*!
\ingroup PkgMeshSmoothing3Concepts
\cgalConcept

The concept `MeshDataStructure` describes the way the tetrahedral mesh will be accessed and modified.

\sa `CGAL::Mesh_optimization::Mesh_smoother`

*/
class MeshDataStructure {
public:

/// \name Types
/// @{

/*!
Descriptor used to access a cell (tetrahedron) information
*/
using Cell_descriptor = unspecified_type;


/*!
Descriptor used to access a vertex information
*/
using Vertex_descriptor = unspecified_type;


/*!
Point type.
*/
using Point_3 = unspecified_type;

/// @}

/// \name Operations
/// The following functions are used to access and modify the mesh data:
/// @{

/*!
Used only to reserve memory. std::size_t is optional but will avoid warnings.
*/
std::size_t nb_cells() const;

/*!

*/
std::size_t nb_vertices() const;

/*!
Access the coordinates of the given point. 
*/
Point_3 vertex_coordinates(Vertex_descriptor vertex) const;

/*!
Change the coordinates of the given point. 
*/
void set_new_vertex_coordinates(Vertex_descriptor vertex, Point_3 coord);

/*!
Provide an iterable range over the Cell_descriptors of the mesh
*/
unspecified_type cell_range() const;

/*!
Access the 4 vertices of a cell. 
Returns container behaving like std::array<Vertex_descriptor, 4> 
*/
unspecified_type cell_vertices(Cell_descriptor cell) const;

/*!
Optimal shape of the given cell. 
Returns container behaving like std::array<Point_3, 4> 
*/
unspecified_type cell_reference_shape(Cell_descriptor cell) const;

/// @}



}; /* end MeshDataStructure */
