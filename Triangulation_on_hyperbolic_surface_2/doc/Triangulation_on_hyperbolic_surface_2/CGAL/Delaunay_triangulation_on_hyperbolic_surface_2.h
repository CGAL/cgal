namespace CGAL {

/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

  This item defines attributes for edges and faces.
  Edge attributes are of type `Complex_number` reprensenting cross-ratios.
  Face attributes attributes are of type `Anchor` representing a lift of a triangle.

  \tparam Traits must be a model of `HyperbolicSurfaceTraits_2`.

  \cgalModels{GenericMapItems}
*/
template<class Traits>
struct Delaunay_triangulation_attributes {
  template<class CMap>
  struct Dart_wrapper{
    typedef Cell_attribute<CMap, Complex_number<typename Traits::FT>> Edge_attrib;
    typedef Cell_attribute<CMap, typename Triangulation_on_hyperbolic_surface_2<Traits,Delaunay_triangulation_attributes<Traits>>::Anchor> Face_attrib;
    typedef std::tuple<void,Edge_attrib,Face_attrib> Attributes;
  };
};

/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

  represents a Delaunay triangulation of a closed orientable hyperbolic surface.

  The class provides functions such as the generation of the Delaunay triangulation from a convex fundamental domain,
  the location of a point in the lifted triangulation, the insertion of a point in the Delaunay triangulation, and the \f$ \varepsilon \f$-net algorithm.

  \tparam Traits must be a model of `HyperbolicSurfaceTraits_2`.

  \sa `Triangulation_on_hyperbolic_surface_2<Traits, Attributes>`
*/
template<class Traits>
class Delaunay_triangulation_on_hyperbolic_surface_2: public Triangulation_on_hyperbolic_surface_2<Traits, Delaunay_triangulation_attributes<Traits>> {
public:
  /// \name Types
  /// @{

  typedef Triangulation_on_hyperbolic_surface_2<Traits,Delaunay_triangulation_attributes<Traits>>                     Base;
  /*!
    Type of combinatorial map whose edges are decorated with complex numbers and faces are decorated with anchors.
  */
  typedef Combinatorial_map<2,Delaunay_triangulation_attributes<Traits>>                                              CMap;
  typedef typename Triangulation_on_hyperbolic_surface_2<Traits,Delaunay_triangulation_attributes<Traits>>::Anchor    Anchor;

  typedef typename Traits::FT                                         Number;
  typedef typename Traits::Complex                                    Complex_number;


  /*!
    Point in the Poincaré disk.
  */
  typedef typename Traits::Hyperbolic_point_2                         Point;

  typedef typename CMap::Dart_descriptor                              Dart_descriptor;
  typedef typename CMap::Dart_range                                   Dart_range;
  typedef typename CMap::template One_dart_per_cell_range<0>          Vertex_range;
  typedef typename CMap::template One_dart_per_cell_range<1>          Edge_range;
  typedef typename CMap::template One_dart_per_cell_range<2>          Face_range;
  typedef typename CMap::Dart_const_descriptor                        Dart_const_descriptor;
  typedef typename CMap::Dart_const_range                             Dart_const_range;
  typedef typename CMap::template One_dart_per_cell_const_range<1>    Edge_const_range;
  typedef typename CMap::template One_dart_per_cell_const_range<2>    Face_const_range;
  /// @}

  /// \name Enums
  /// @{

  /*
    The enumeration `Locate_type` is defined to specify which case occurs when locating a point in the lifted triangulation.
  */
  enum Locate_type {
    VERTEX = 0,
    EDGE,
    FACE,
    OUTSIDE
  };
  /// @}

  /// \name Creation
  /// Calls the corresponding constructor from #Base and sets an #Anchor for each face.
  /// @{
  Delaunay_triangulation_on_hyperbolic_surface_2() {};
  Delaunay_triangulation_on_hyperbolic_surface_2(CMap & cmap, Anchor & anch);
  Delaunay_triangulation_on_hyperbolic_surface_2(Hyperbolic_fundamental_domain_2<Traits> const & domain);
  Delaunay_triangulation_on_hyperbolic_surface_2(Base & triangulation);
  /// @}

  /// \name Access Functions
  /// @{
  /*!
    \return the anchor associated with the given dart.
  */
  Anchor & anchor(Dart_descriptor const dart);
  /*!
    \return the anchor associated with the given dart.
  */
  Anchor const & anchor(Dart_const_descriptor const dart) const;
  /*!
    \return an anchor of the triangulation.
  */
  Anchor & anchor();
  /*!
    \return an anchor of the triangulation.
  */
  Anchor const & anchor() const;
  /*!
    \return the index of the given dart in its face. The index is `0` if `dart` is the dart of its associated anchor.
  */
  unsigned index_in_anchor(Dart_const_descriptor const dart) const;
  /*!
    \return the `i`-th dart of the face associated with `anch`, where `anch.dart` is the `0`-th dart of the face.
  */
  Dart_descriptor ith_dart(unsigned i, Anchor const & anch);
  /// @}

  /// \name Validity
  /*!
    \return a Boolean that indicates whether the triangulation is a valid Delaunay triangulation.
  */
  bool is_valid() const;

  /// \name Point location and insertion
  /// @{
  /*!
    \return a `Locate_type` that indicates whether `query` lies on a vertex, an edge, inside or outside the lifted triangle described by the given anchor.

    If `query` lies on a vertex or an edge, the index `li` is set to the index of the vertex or edge on which `query` lies.

    If `query` lies inside the triangle, `li` is set to `-1`.

    If `query` lies outside the triangle, `li` is set to the index of the first edge such that `query` and the third point of the triangle lies on different sides.
  */
  Locate_type relative_position(Point const & query, unsigned & li, Anchor const & anch) const;

  /*!
    \return the anchor representing the lift of the triangle in which `query` lies.

    Locates `query` in the Delaunay triangulation.

    The Boolean `use_visibility` indicates whether the visibility walk algorithm is used for the point location. When `false`, the straight walk algorithm is used.
  */
  Anchor locate(Point const & query, bool use_visibility = false); // const ?
  /*!
    Same as above.

    Additionally, the point location algorithm starts from `hint` and `ld` is set to the number or triangle travarsed by the walk used for the point location algorithm.

    The variables `lt` and `li` contains information about the triangle in which `query` has been found.

    \sa relative_position
  */
  Anchor locate(Point const & query, Locate_type & lt, unsigned & li, unsigned & ld, Anchor const & hint, bool use_visibility = false); // const ?

  /*!
    Inserts `query` in the Delaunay triangulation after having located it, starting from `hint`.

    \pre <code>is_valid()</code> and <code>norm(Complex_number(query.x(), query.y())) < Number(1)</code>
  */
  void insert(Point const & query, Anchor & hint);

  /*!
    Inserts `query` in the Delaunay triangulation.
  */
  void insert(Point const & query);
  /// @}

  /// \name \f$ \varepsilon \f$-net
  /// @{
  /*!
    \return a Boolean that indicates whether the vertices of the Delaunay triangulation form a certified `epsilon`-net of the surface.

    Tries to compute an `epsilon`-net of the surface.

    When `Number` is `CGAL::Gmpq`, the algorithm rounds the coordinate of circumcenters to `CGAL::Gmpfr` with precision `p*53`.
  */
  bool epsilon_net(double const epsilon, unsigned const p = 1);
  /*!
    \return a Boolean that indicates whether the vertices of the Delaunay triangulation form a certified `epsilon`-covering of the surface.

    \pre <code>is_epsilon_packing(epsilon)</code> and <code>p > 0</code>
  */
  bool is_epsilon_covering(const double epsilon) const;
  /*!
    \return a Boolean that indicates whether the vertices of the Delaunay triangulation form a certified `epsilon`-packing of the surface.
  */
  bool is_epsilon_packing(const double epsilon) const;
  /*!
    \return a Boolean that indicates whether the vertices of the Delaunay triangulation form a certified `epsilon`-net of the surface.
  */
  bool is_epsilon_net(const double epsilon) const;
  /// @}

  /// \name Other functions
  /// @{
  /*!
    \return a `double` approximation of the length of the shortest geodesic loop in the Delaunay triangulation.
  */
  double shortest_loop() const;
  /*!
    \return a `double` approximation of the length of the shortest non-loop edge in the Delaunay triangulation.
  */
  double shortest_non_loop_edge() const;
  /// @}
}

} // namespace CGAL
