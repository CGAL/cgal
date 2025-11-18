namespace CGAL {

/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

  This item defines attributes of edges that are of type `Complex_number` reprensenting cross-ratios.

  \tparam Traits must be a model of `HyperbolicSurfaceTraits_2`.

  \cgalModels{GenericMapItems}
*/
template<class Traits>
struct Combinatorial_map_with_cross_ratios_item
{
  template <class CMap>
  struct Dart_wrapper
  {
    typedef Cell_attribute<CMap, Complex_number<typename Traits::FT> > Edge_attrib;
    typedef std::tuple<void, Edge_attrib, void> Attributes;
  };
};

/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

  represents a triangulation of a closed orientable hyperbolic surface.

  The class provides functions such as the generation of the triangulation from a convex fundamental domain,
  the Delaunay flip algorithm, and the construction of a portion of the lift of the triangulation in the hyperbolic plane.

  \tparam Traits must be a model of `HyperbolicSurfaceTraits_2`.
  \tparam Attributes must be a model of `GenericMapItems` whose edges are
  decorated with complex numbers to represent cross ratios.
*/
template<class Traits, class Attributes = Combinatorial_map_with_cross_ratios_item<Traits> >
class Triangulation_on_hyperbolic_surface_2
{
public:
  /// \name Types
  /// @{

  /*!
    Type of combinatorial map whose edges are decorated with complex numbers.
  */
 typedef Combinatorial_map<2, Attributes> Combinatorial_map_with_cross_ratios;

  /*!
    Combinatorial map dart descriptor type.
  */
  typedef typename Combinatorial_map_with_cross_ratios::Dart_descriptor Dart_descriptor;

  /*!
      Combinatorial map dart const descriptor type.
  */
  typedef typename Combinatorial_map_with_cross_ratios::Dart_const_descriptor Dart_const_descriptor;

  /*!
    Range of one dart for each vertex (that is 0-cell) of the combinatorial map.
  */
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_range<0> Vertex_range;


  /*!
    Range of one dart for each edge (that is 1-cell) of the combinatorial map.
  */
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_range<1> Edge_range;

  /*!
    Range of one dart for each face (that is 2-cell) of the combinatorial map.
  */
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_range<2> Face_range;

  /*!
      Range of one dart for each vertex (that is 0-cell) of the combinatorial map.
  */
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_const_range<0> Vertex_const_range;

  /*!
    Range of one dart for each edge (that is 1-cell) of the combinatorial map.
  */
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_const_range<1> Edge_const_range;

  /*!
    Range of one dart for each face (that is 2-cell) of the combinatorial map.
  */
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_const_range<2> Face_const_range;

  /*!
    Point type.
  */
  typedef typename Traits::Hyperbolic_point_2 Point;

  /*!
    stores a dart \f$ d \f$ of the combinatorial map, belonging to a triangle \f$ t \f$,
    and stores the three vertices of a lift of \f$ t \f$ in the hyperbolic plane.
  */
  struct Anchor
  {
    typename Combinatorial_map_with_cross_ratios::Dart_descriptor dart;
    typename Traits::Hyperbolic_point_2 vertices[3];
  };

  /// @}

  /// \name Creation
  /// @{

  /*!
    Default constructor.
  */
  Triangulation_on_hyperbolic_surface_2() {};

  /*!
    Constructor from a convex fundamental domain: triangulates the polygon of
    the domain. (The triangulation is defined by adding an internal edge
    between domain.vertex(size-1) and the other vertices. The anchor has the
    three vertices of indices size-1, 0 and 1 and the dart is the one between
    the vertices of indices size-1 and 0.)
  */
  Triangulation_on_hyperbolic_surface_2(const Hyperbolic_fundamental_domain_2<Traits>& domain);

  /*!
    Constructor from a decorated combinatorial map and an anchor.
  */
  Triangulation_on_hyperbolic_surface_2(Combinatorial_map_with_cross_ratios& cmap, Anchor& an_anchor);

  /// @}

  /// @}

  /// \name Access Functions
  /// @{

  /*!
    returns the decorated combinatorial map.
  */
  Combinatorial_map_with_cross_ratios& combinatorial_map();

  /*!
    returns whether the triangulation has an anchor or not.
    \pre <code> is_valid() </code>
  */
  bool has_anchor() const;

  /*!
    returns the anchor.
    \pre <code> is_valid() && has_anchor() </code>
  */
  Anchor& anchor();

  /*!
    returns the anchor.
    \pre <code> is_valid() && has_anchor() </code>
  */
  const Anchor& anchor() const;

  /*!
    returns the range of vertices.
  */
Vertex_range vertices_range();

  /*!
    returns the range of edges.
  */
  Edge_range edges_range();

  /*!
    returns the range of faces.
  */
  Face_range faces_range();

  /*!
    returns the range of vertices.
  */
  Vertex_const_range vertices_const_range() const;

  /*!
    returns the range of edges.
  */
  Edge_const_range edges_const_range() const;

  /*!
    returns the range of faces.
  */
  Face_const_range faces_const_range() const;

  /// @}

  /// \name Delaunay Flip Algorithm
  /// @{

  /*!
    returns whether the edge supported by the dart is Delaunay flippable or not. An edge \f$ e \f$
    is Delaunay flippable if the imaginary part of its cross ratio is positive.

    \pre <code> is_valid() </code>
  */
  bool is_Delaunay_flippable(Dart_descriptor dart) const;

  /*!
    flips the edge supported by the dart.

    \pre <code> is_valid() </code>
  */
  void flip(Dart_descriptor dart);

  /*!
    determines if the triangulation is a valid Delaunay triangulation.
  */
  bool is_Delaunay() const;

  /*!
    applies the Delaunay flip algorithm: flips Delaunay-flippable edges until there is no such edge anymore.
    \pre <code> is_valid() </code>
  */
  int make_Delaunay();

  /// @}

  /// \name Lifting
  /// @{

  /*!
    lifts the triangulation in the hyperbolic plane.

    Returns, for every triangle \f$ t \f$ of the triangulation, one of the darts of \f$ t \f$ in the combinatorial map of the triangulation, together with a triple \f$ p,q,r \f$ of points in the hyperbolic plane.
    The points \f$ p,q,r \f$ are the vertices of a lift of \f$ t \f$ in the hyperbolic plane.
    If the center parameter is set to true, then one of the vertices of the anchor is translated to the origin \f$ 0 \f$.

    \pre <code> is_valid() && has_anchor() </code>
  */
  std::vector<std::tuple<Dart_const_descriptor, Point, Point, Point>> lift(bool center=true) const;

  /// @}

  /// \name Validity
  /// @{

  /*!
    checks that the underlying combinatorial map \f$ M \f$ has no boundary and calls the is_valid method of \f$ M \f$.

    If there is an anchor, then checks that the dart descriptor of the anchor does indeed point to a dart of \f$ M \f$,
    and checks that the three vertices of the anchor lie within the open unit disk.
  */
  bool is_valid() const;

  /// @}
};

} // namespace CGAL
