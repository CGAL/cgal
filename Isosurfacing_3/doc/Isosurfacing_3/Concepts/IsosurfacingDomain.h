/*!
\ingroup PkgIsosurfacing3Concepts
\cgalConcept

The concept `IsosurfacingDomain` describes the set of requirements to be
fulfilled by any class used as input data for any isosurfacing algorithms.

\cgalHasModel `CGAL::Isosurfacing::Explicit_cartesian_grid_domain`
\cgalHasModel `CGAL::Isosurfacing::Implicit_cartesian_grid_domain`

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
    Must be a model of the concept `RandomAccessContainer` with size 2 whose value type is `Vertex_descriptor`.
    */
    typedef unspecified_type Vertices_incident_to_edge;

    /*!
    A container for the cells incident to an edge.
    Must be a model of the concept `RandomAccessContainer` whose value type is `Cell_descriptor`.
    */
    typedef unspecified_type Cells_incident_to_edge;

    /*!
    A container for the vertices of a cell.
    Must be a model of the concept `RandomAccessContainer` whose value type is `Vertex_descriptor`.
    */
    typedef unspecified_type Cell_vertices;

    /*!
    A container for the edges of a cell.
    Must be a model of the concept `RandomAccessContainer` whose value type is `Edge_descriptor`.
    */
    typedef unspecified_type Cell_edges;


    /// @}

    /// \name Operations
    /// The following member functions must exist.
    /// @{

    /*!
    Get the position of a vertex in 3D space

    \param v the descriptor of the vertex

    \return the position of the vertex as a point
    */
    Point position(const Vertex_descriptor& v) const;

    /*!
    Get the value of the function at a vertex

    \param v the descriptor of the vertex

    \return the value of the function
    */
    FT value(const Vertex_descriptor& v) const;

    /*!
    Get the two vertices incident to an edge

    \param e the descriptor of the edge

    \return a collection of the two vertex descriptors
    */
    Vertices_incident_to_edge edge_vertices(const Edge_descriptor& e) const;

    /*!
    Get all cells incident to an edge

    \param e the descriptor of the edge

    \return a collection of cell descriptors
    */
    Cells_incident_to_edge cells_incident_to_edge(const Edge_descriptor& e) const;

    /*!
    Get all vertices of the a cell

    \param c the descriptor of the cell

    \return a collection of vertex descriptors
    */
    Cell_vertices cell_vertices(const Cell_descriptor& c) const;

    /*!
    Get all edges of the cell c

    \param c the descriptor of the cell

    \return a collection of edge descriptors
    */
    Cell_edges cell_edges(const Cell_descriptor& c) const;

    /*!
    Iterate sequentially over all vertices and call a functor on each one

    \param f the functor called with every vertex. Must implement `void operator()(const Vertex_descriptor& vertex)`.
    */
    template <typename Functor>
    void iterate_vertices(Functor& f, Sequential_tag) const;

    /*!
    Iterate sequentially over all edges and call the functor f on each one

    \param f the functor called with every edge. Must implement `void operator()(const Edge_descriptor& edge)`.
    */
    template <typename Functor>
    void iterate_edges(Functor& f, Sequential_tag) const;

    /*!
    Iterate sequentially over all cells and call the functor f on each one

    \param f the functor called with every cell. Must implement `void operator()(const Cell_descriptor& cell)`.
    */
    template <typename Functor>
    void iterate_cells(Functor& f, Sequential_tag) const;

    /*!
    (Optional) Iterate in parallel over all vertices and call a functor on each one

    \param f the functor called with every vertex. Must implement `void operator()(const Vertex_descriptor& vertex)`.
    */
    template <typename Functor>
    void iterate_vertices(Functor& f, Parallel_tag) const;

    /*!
    (Optional) Iterate in parallel over all edges and call the functor f on each one

    \param f the functor called with every edge. Must implement `void operator()(const Edge_descriptor& edge)`.
    */
    template <typename Functor>
    void iterate_edges(Functor& f, Parallel_tag) const;

    /*!
    (Optional) Iterate in parallel over all cells and call the functor f on each one

    \param f the functor called with every cell. Must implement `void operator()(const Cell_descriptor& cell)`.
    */
    template <typename Functor>
    void iterate_cells(Functor& f, Parallel_tag) const;

    /// @}
};
