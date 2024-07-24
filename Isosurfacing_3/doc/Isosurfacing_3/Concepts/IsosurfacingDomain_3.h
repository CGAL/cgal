/*!
\ingroup PkgIsosurfacing3Concepts

\cgalConcept

\brief The concept `IsosurfacingDomain_3` describes the set of requirements to be
fulfilled by any class used as input data for isosurfacing algorithms.

A model of the concept `IsosurfacingDomain_3` provides a partition of the Euclidean space in cells,
and a scalar field defined over the whole partition.
The isosurfacing algorithms traverse these cells and query the domain class
at the vertices of each cell, using the functions `point()` and `value()`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Isosurfacing::Marching_cubes_domain_3}
\cgalHasModels{CGAL::Isosurfacing::Dual_contouring_domain_3}
\cgalHasModelsEnd

\sa `IsosurfacingDomainWithGradient_3`
*/
class IsosurfacingDomain_3
{
public:
  /// \name Types
  /// @{

  /*!
  The geometric traits type.
  Must be a model of `IsosurfacingTraits_3`.
  */
  typedef unspecified_type Geom_traits;

  /*!
  The scalar type.
  */
  typedef Geom_traits::FT FT;

  /*!
  The 3D point type.
  */
  typedef Geom_traits::Point_3 Point_3;

  /*!
  A descriptor that uniquely identifies a vertex.
  Must be a model of the concepts `Descriptor` and `Hashable`.
  */
  typedef unspecified_type vertex_descriptor;

  /*!
  A descriptor that uniquely identifies an edge.
  Must be a model of the concept `Descriptor` and `Hashable`.
  */
  typedef unspecified_type edge_descriptor;

  /*!
  A descriptor that uniquely identifies a cell.
  Must be a model of the concepts `Descriptor` and `Hashable`.
  */
  typedef unspecified_type cell_descriptor;

  /*!
  A container for the two vertices of an edge.
  Must be a model of the concept `RandomAccessContainer` of size `2` whose value type is `vertex_descriptor`.
  */
  typedef unspecified_type Edge_vertices;

  /*!
  A container for the cells incident to an edge.
  Must be a model of the concept `ForwardRange` whose value type is `cell_descriptor`.
  */
  typedef unspecified_type Cells_incident_to_edge;

  /*!
  A container for the vertices of a cell.
  Must be a model of the concept `ForwardRange` whose value type is `vertex_descriptor`.
  */
  typedef unspecified_type Cell_vertices;

  /*!
  A container for the edges of a cell.
  Must be a model of the concept `ForwardRange` whose value type is `edge_descriptor`.
  */
  typedef unspecified_type Cell_edges;


  /// @}

  /// \name Operations
  /// @{

  /*!
  returns the geometric traits.
  */
  Geom_traits geom_traits();

  /*!
  returns the 3D location of the vertex `v`.
  */
  Point_3 point(vertex_descriptor v) const;

  /*!
  returns the value of the value field at the point `p`.
  */
  FT value(Point_3 p) const;

  /*!
  returns the value of the value field at the vertex `v`.
  */
  FT value(vertex_descriptor v) const;

  /*!
  returns the two vertices incident to the edge `e`.
  */
  Edge_vertices incident_vertices(edge_descriptor e) const;

  /*!
  returns all the cells incident to the edge `e`, in a clockwise or counterclockwise order.
  */
  Cells_incident_to_edge incident_cells(edge_descriptor e) const;

  /*!
  returns all the vertices of the cell `c`.
  */
  Cell_vertices cell_vertices(cell_descriptor c) const;

  /*!
  returns all the edges of the cell `c`.
  */
  Cell_edges cell_edges(cell_descriptor c) const;

  /*!
  iterates over all vertices, and calls the functor `f` on each one.

  \tparam ConcurrencyTag decides if the vertices are iterated sequentially or in parallel.
  Can be either `CGAL::Sequential_tag`, `CGAL::Parallel_if_available_tag`, or `CGAL::Parallel_tag`.
  \tparam Functor must implement `void operator()(vertex_descriptor vertex)`

  \param f the functor called on every vertex
  */
  template <typename ConcurrencyTag, typename Functor>
  void for_each_vertex(Functor& f) const;

  /*!
  iterates over all edges, and calls the functor `f` on each one.

  \tparam ConcurrencyTag decides if the edges are iterated sequentially or in parallel.
  Can be either `CGAL::Sequential_tag`, `CGAL::Parallel_if_available_tag`, or `CGAL::Parallel_tag`.
  \tparam Functor must implement `void operator()(edge_descriptor edge)`.

  \param f the functor called on every edge
  */
  template <typename ConcurrencyTag, typename Functor>
  void for_each_edge(Functor& f) const;

  /*!
  iterates over all cells, and calls the functor `f` on each one.

  \tparam ConcurrencyTag decides if the cells are iterated sequentially or in parallel.
  Can be either `CGAL::Sequential_tag`, `CGAL::Parallel_if_available_tag`, or `CGAL::Parallel_tag`.
  \tparam Functor must implement `void operator()(cell_descriptor cell)`.

  \param f the functor called on every cell
  */
  template <typename ConcurrencyTag, typename Functor>
  void for_each_cell(Functor& f) const;

  /*!
  Constructs the intersection - if it exists - between an edge and an isosurface.

  \param p_0 the location of the first vertex of the edge
  \param p_1 the location of the second vertex of the edge
  \param val_0 the value at the first vertex of the edge
  \param val_1 the value at the second vertex of the edge
  \param isovalue the isovalue defining the isosurface with which we seek an intersection
  \param p the intersection point, if it exists

  \returns `true` if the intersection point exists, `false` otherwise.
  */
  bool construct_intersection(Point_3 p_0, Point_3 p_1,
                              FT val_0, FT val_1,
                              FT isovalue,
                              Point_3& p) const;

  /// @}
};
