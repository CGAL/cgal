/*!
\ingroup PkgIsosurfacing3Concepts
\cgalConcept

The concept `IsosurfacingDomain` describes the set of requirements to be
fulfilled by any class used as input data for any isosurfacing algorithms.

\cgalHasModel `CGAL::Isosurfacing::Explicit_cartesian_grid_domain`
\cgalHasModel `CGAL::Isosurfacing::Implicit_cartesian_grid_domain`
\cgalHasModel `CGAL::Isosurfacing::Implicit_octree_domain`

*/

class IsosurfacingDomain {
public:
    /// \name Types
    /// @{

    /*!
    Traits type model of \cgal %Kernel
    */
    typedef unspecified_type Geom_traits;

    /*!
    The scalar type.
    */
    typedef unspecified_type FT;

    /*!
    The point type.
    */
    typedef unspecified_type Point;

    /*!
    A descriptor to uniquely identify a vertex.
    */
    typedef unspecified_type Vertex_descriptor;

    /*!
    A descriptor to uniquely identify an edge.
    */
    typedef unspecified_type Edge_descriptor;

    /*!
    A descriptor to uniquely identify a cell.
    */
    typedef unspecified_type Cell_descriptor;

    /*!
    A container for the two vertices of an edge.
    */
    typedef unspecified_type Vertices_incident_to_edge;

    /*!
    A container for the cells incident to an edge.
    */
    typedef unspecified_type Cells_incident_to_edge;

    /*!
    A container for the vertices of a cell.
    */
    typedef unspecified_type Cell_vertices;

    /*!
    A container for the edges of a cell.
    */
    typedef unspecified_type Cell_edges;


    /// @}

    /// \name Operations
    /// The following member functions must exist.
    /// @{

    /*!
    Returns the position of vertex v in 3D space
    */
    Point position(const Vertex_descriptor& v) const;

    /*!
    Returns the value of vertex v
    */
    FT value(const Vertex_descriptor& v) const;

    /*!
    Returns the two vertices incident to edge e
    */
    Vertices_incident_to_edge edge_vertices(const Edge_descriptor& e) const;

    /*!
    Returns all voxels incident to edge e
    */
    Cells_incident_to_edge cells_incident_to_edge(const Edge_descriptor& e) const;

    /*!
    Returns all vertices of the cell c
    */
    Cell_vertices cell_vertices(const Cell_descriptor& c) const;

    /*!
    Returns all edges of the cell c
    */
    Cell_edges cell_edges(const Cell_descriptor& c) const;

    /*!
    Iterate sequentially over all vertices and call the functor f on each one
    */
    template <typename Functor>
    void iterate_vertices(Functor& f, Sequential_tag) const;

    /*!
    Iterate sequentially over all edges and call the functor f on each one
    */
    template <typename Functor>
    void iterate_edges(Functor& f, Sequential_tag) const;

    /*!
    Iterate sequentially over all cells and call the functor f on each one
    */
    template <typename Functor>
    void iterate_cells(Functor& f, Sequential_tag) const;

    /*!
    (Optional) Iterate in parallel over all vertices and call the functor f on each one
    */
    template <typename Functor>
    void iterate_vertices(Functor& f, Parallel_tag) const;

    /*!
    (Optional) Iterate in parallel over all edges and call the functor f on each one
    */
    template <typename Functor>
    void iterate_edges(Functor& f, Parallel_tag) const;

    /*!
    (Optional) Iterate in parallel over all cells and call the functor f on each one
    */
    template <typename Functor>
    void iterate_cells(Functor& f, Parallel_tag) const;

    /// @}
};
