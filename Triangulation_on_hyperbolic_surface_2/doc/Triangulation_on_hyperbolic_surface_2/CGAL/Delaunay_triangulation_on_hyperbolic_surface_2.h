namespace CGAL {

/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

  This item defines attributes for edges and faces.
  Edge attributes are of type `Complex_number` representing cross-ratios encoding the metric of the surface.
  Face attributes are of type `Anchor`, representing a lift of a triangle.

  \tparam Traits must be a model of `HyperbolicSurfaceTraits_2` and the same as
  the one of the associated triangulation.

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

  \tparam Traits must be a model of `HyperbolicSurfaceDelaunayTraits_2`.

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

  /*!
    must be compatible with the method  `CGAL::to_double()` and must support a cosh function to compute the hyperbolic cosine.
   */
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

  /*!
    The enumeration `Locate_type` is defined to specify which case occurs when locating a point in the lifted triangulation.
  */
  enum Locate_type {
    VERTEX = 0,
    EDGE,
    FACE,
    OUTSIDE
  };

    /*!
    The enumeration `Locate_walk` is defined to specify which walk algorithm is
    used when locating a point in the lifted triangulation.

    The straight walk visits all triangles along the geodesic segment pq. The
    algorithm starts by rotating around p to find a triangle incident to p
    that has an edge intersecting the geodesic segment pq. The visibility walk
    consists, for each visited triangle not containing q, of moving to a
    neighbor through an edge e if q and the third vertex of the visited
    triangle are on different sides of the geodesic line supporting e.
  */
  enum Locate_walk {
    STRAIGHT = 0,
    VISIBILITY,
  };
  /// @}

  /// \name Creation
  /// The main Delaunay class can be constructed from a general triangulation or a hyperbolic fundamental domain. It constructs the Delaunay
  /// triangulation and sets an #Anchor for each face.
  /// @{
  Delaunay_triangulation_on_hyperbolic_surface_2(Traits & gt) {};
  Delaunay_triangulation_on_hyperbolic_surface_2(Traits & gt, Hyperbolic_fundamental_domain_2<Traits> const & domain);
  Delaunay_triangulation_on_hyperbolic_surface_2(Traits & gt, Base & triangulation);
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

    If `query` lies outside the triangle, `li` is set to the index of the first edge such that `query` and the third point of the triangle lie on different sides.
  */
  Locate_type relative_position(Point const & query, unsigned & li, Anchor const & anch) const;

  /*!
    \return the anchor representing the lift of the triangle in which `query` lies.

    Locates `query` in the Delaunay triangulation.
    The enumeration `walk` indicates which walk algorithm is used when locating a point in the lifted triangulation.
    The point location algorithm starts from `hint` and `ld` is set to the number of triangles traversed by the walk used for the point location algorithm.
    The variable `lt` indicates the location of `query` with respect to  the
    lifted triangle described by the returned anchor,  `li` is an index
    precising its location as described in the method `relative_position`.

    \sa relative_position
  */
  Anchor locate(Point const & query, Locate_type & lt, unsigned & li, unsigned & ld, Anchor const & hint, Locate_walk walk = STRAIGHT); // const ?

  /*!
    \return the anchor representing the lift of the triangle in which `query` lies.

    Locates `query` in the Delaunay triangulation.

    The enumeration `walk` indicates which walk algorithm is used when locating a point in the lifted triangulation.
  */
  Anchor locate(Point const & query, Locate_walk walk = STRAIGHT); // const ?

  /*!
    Inserts `p` in the Delaunay triangulation after having located it, starting from `hint`.

    \pre <code>is_valid()</code> and <code>norm(Complex_number(p.x(), p.y())) < Number(1)</code>
  */
  void insert(Point const & p, Anchor & hint);

  /*!
    Inserts `p` in the Delaunay triangulation.
  */
  void insert(Point const & p);
  /// @}

  /// \name epsilon-net
  /// @{
  /*!
    \return a Boolean that indicates whether the vertices of the Delaunay triangulation form a certified `epsilon`-net of the surface.

    Tries to compute an `epsilon`-net of the surface.

    When `Number` is `CGAL::Gmpq` and Kernel is `CGAL::Simple_cartesian<Number>` or `CGAL::Cartesian<Number>`,
    the traits class
    CGAL::Hyperbolic_surface_Delaunay_traits_2<
    CGAL::Hyperbolic_surface_traits_2<
    CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<
    CGAL::Circular_kernel_2<Kernel,CGAL::Algebraic_kernel_for_circles_2_2<Number>>
    >>>
    enables to choose the approximation precision via the
    function `set_precision()`.

    If another number type is used, the approximation uses the function
    CGAL::to_double() and there are no general guarantees whatsoever.

    Note that a Delaunay triangulation with a single vertex is always an
    espsilon-packing for any non-negative epsilon.

    \pre <code>is_epsilon_packing(epsilon)</code> and <code>p > 0</code>
  */
  bool construct_epsilon_net(double const epsilon);
  /*!
    \return a Boolean that indicates whether the vertices of the Delaunay triangulation form a certified `epsilon`-covering of the surface.
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

  /*!  \return the optimal packing value up to double precision. Returns
    inf (+infinity) if all edges are loops.
  */
  double packing_value() const;

  /*! \return the optimal covering value up to double precision.   */
  double covering_value() const;

  /*! \return a shortest loop edge of the triangulation if any, otherwise
      returns a null pointer.
   */
  Dart_const_descriptor shortest_loop_edge() const;

  /*!  Sets `p` such that when FT is Gmpq, the precision for the function
    `HyperbolicSurfaceDelaunayTraitsClass::Construct_approximate_hyperbolic_circumcenter_2`
    is `p` x 53 bits, where 53 is the precision of a double.
  */
  void set_circumcenter_approximation_precision(unsigned p);

  /*!  \return `p` such that when FT is Gmpq, the precision for the function
    `HyperbolicSurfaceDelaunayTraitsClass::Construct_approximate_hyperbolic_circumcenter_2`
    is `p` x 53 bits, where 53 is the precision of a double.
  */
  unsigned get_circumcenter_approximation_precision();

    /// @}

  /// \name Other functions
  /// @{
  /*!
    \return the hyperbolic cosine of the length of the edge associated to a dart.
   */
  Number dart_cosh_length(Dart_const_descriptor dart) const;

  /// @}
}

} // namespace CGAL
