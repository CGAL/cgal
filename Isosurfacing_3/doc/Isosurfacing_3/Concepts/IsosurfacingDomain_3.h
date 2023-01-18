/*!
\ingroup PkgIsosurfacing3Concepts

\cgalConcept

\brief The concept `IsosurfacingDomain_3` describes the set of requirements to be
fulfilled by any class used as input data for isosurfacing algorithms.

A model of the concept `IsosurfacingDomain_3` provides a discrete representation
of an implicit function through a partition of the Euclidean space in cells.
The isosurfacing algorithms will traverse these cells and query the domain class
at the vertices of each cell, using the functions `point()` and `value()`.

\cgalHasModel `CGAL::Isosurfacing::Explicit_Cartesian_grid_domain_3`
\cgalHasModel `CGAL::Isosurfacing::Implicit_Cartesian_grid_domain_3`
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
  A descriptor to uniquely identify a vertex.
  Must be a model of the concepts `DefaultConstructible`, `CopyConstructible`, and `Assignable`.
  */
  typedef unspecified_type Vertex_descriptor;

  /*!
  A descriptor to uniquely identify an edge.
  Must be a model of the concepts `DefaultConstructible`, `CopyConstructible`, and `Assignable`.
  */
  typedef unspecified_type Edge_descriptor;

  /*!
  A descriptor to uniquely identify a cell.
  Must be a model of the concepts `DefaultConstructible`, `CopyConstructible`, and `Assignable`.
  */
  typedef unspecified_type Cell_descriptor;

  /*!
  A container for the two vertices of an edge.
  Must be a model of the concept `RandomAccessContainer` of size `2` whose value type is `Vertex_descriptor`.
  */
  typedef unspecified_type Vertices_incident_to_edge;

  /*!
  A container for the cells incident to an edge.
  Must be a model of the concept `Range` whose value type is `Cell_descriptor`.
  */
  typedef unspecified_type Cells_incident_to_edge;

  /*!
  A container for the vertices of a cell.
  Must be a model of the concept `Range` whose value type is `Vertex_descriptor`.
  */
  typedef unspecified_type Cell_vertices;

  /*!
  A container for the edges of a cell.
  Must be a model of the concept `Range` whose value type is `Edge_descriptor`.
  */
  typedef unspecified_type Cell_edges;


  /// @}

  /// \name Operations
  /// The following member functions must exist.
  /// @{

  /*!
  gets the geometric traits
  */
  Geom_traits geom_traits();

  /*!
  gets the position of a vertex in 3D space
  */
  Point_3 point(const Vertex_descriptor& v) const;

  /*!
  gets the value of the implicit function at a vertex
  */
  FT value(const Vertex_descriptor& v) const;

  /*!
  gets the two vertices incident to an edge
  */
  Vertices_incident_to_edge incident_vertices(const Edge_descriptor& e) const;

  /*!
  gets all cells incident to an edge
  */
  Cells_incident_to_edge incident_cells(const Edge_descriptor& e) const;

  /*!
  gets all vertices of the a cell
  */
  Cell_vertices cell_vertices(const Cell_descriptor& c) const;

  /*!
  gets all edges of the cell c
  */
  Cell_edges cell_edges(const Cell_descriptor& c) const;

  /*!
  iterates over all vertices and call a functor on each one

  \tparam ConcurrencyTag decides if the vertices are iterated sequentially or in parallel.
  Can be either `CGAL::Sequential_tag`, `CGAL::Parallel_if_available_tag`, or `CGAL::Parallel_tag`.

  \tparam Functor must implement `void operator()(const Vertex_descriptor& vertex)`

  \param f the functor called with every vertex
  */
  template <typename ConcurrencyTag, typename Functor>
  void iterate_vertices(Functor& f) const;

  /*!
  iterates over all edges and call a functor f on each one

  \tparam ConcurrencyTag decides if the edges are iterated sequentially or in parallel.
  Can be either `CGAL::Sequential_tag`, `CGAL::Parallel_if_available_tag`, or `CGAL::Parallel_tag`.

  \tparam Functor must implement `void operator()(const Edge_descriptor& edge)`.
  */
  template <typename ConcurrencyTag, typename Functor>
  void iterate_edges(Functor& f) const;

  /*!
  iterates over all cells and call a functor f on each one

  \tparam ConcurrencyTag decides if the cells are iterated sequentially or in parallel.
  Can be either `CGAL::Sequential_tag`, `CGAL::Parallel_if_available_tag`, or `CGAL::Parallel_tag`.

  \tparam Functor must implement `void operator()(const Cell_descriptor& cell)`.
  */
  template <typename ConcurrencyTag, typename Functor>
  void iterate_cells(Functor& f) const;

  /// @}
};
