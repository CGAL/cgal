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


/// \name Enums
/// @{

/*
The enumeration `Locate_type` is defined to specify which case occurs when locating a point in the lifted triangulation.
*/
enum Locate_type {
  VERTEX,
  EDGE,
  FACE,
  OUTSIDE = 4
};
/// @}

/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

  represents a Delaunay triangulation of a closed orientable hyperbolic surface.

  The class provides functions such as the generation of the Delaunay triangulation from a convex fundamental domain,
  the location of a point in the lifted triangulation, the insertion of a point in the Delaunay triangulation, and the \f$ \varepsilon \f$-net algorithm.

  \tparam Traits must be a model of `HyperbolicSurfaceTraits_2`.
*/
template<class Traits>
class Delaunay_triangulation_on_hyperbolic_surface_2: public Triangulation_on_hyperbolic_surface_2<Traits, Delaunay_triangulation_attributes<Traits>> {
public:
  /// \name Types
  /// @{

  /*!
    Parent type.
  */
  typedef Triangulation_on_hyperbolic_surface_2<Traits,Delaunay_triangulation_attributes<Traits>>                     Base;

  /*!
    Type of combinatorial map whose edges are decorated with complex numbers and faces are decorated with anchors.
  */
  typedef Combinatorial_map<2,Delaunay_triangulation_attributes<Traits>>                                              CMap;

  /*!
    \link Anchor \endlink type.
  */
  typedef typename Triangulation_on_hyperbolic_surface_2<Traits,Delaunay_triangulation_attributes<Traits>>::Anchor    Anchor;

  typedef typename Traits::FT                                         Number;
  typedef typename Root_of_traits<Number>::Root_of_2                  Algebraic_number;
  typedef typename Traits::Complex                                    Complex_number;
  typedef typename Traits::Hyperbolic_point_2                         Point;
  typedef typename Traits::Hyperbolic_Voronoi_point_2                 Voronoi_point;
  typedef typename CMap::Dart_descriptor                              Dart_descriptor;
  typedef typename CMap::Dart_range                                   Dart_range;
  typedef typename CMap::template One_dart_per_cell_range<0>          Vertex_range;
  typedef typename CMap::template One_dart_per_cell_range<1>          Edge_range;
  typedef typename CMap::template One_dart_per_cell_range<2>          Face_range;
  typedef typename CMap::Dart_const_descriptor                        Dart_const_descriptor;
  typedef typename CMap::Dart_const_range                             Dart_const_range;
  typedef typename CMap::template One_dart_per_cell_const_range<1>    Edge_const_range;
  typedef typename CMap::template One_dart_per_cell_const_range<2>    Face_const_range;
  typedef Hyperbolic_isometry_2<Traits>                               Isometry;
  typedef Hyperbolic_fundamental_domain_2<Traits>                     Domain;
  /// @}

  //---------- CONSTRUCTORS
  Delaunay_triangulation_on_hyperbolic_surface_2() {};
  Delaunay_triangulation_on_hyperbolic_surface_2(CMap & cmap, Anchor & anch);
  Delaunay_triangulation_on_hyperbolic_surface_2(Domain const & domain);
  Delaunay_triangulation_on_hyperbolic_surface_2(Base & triangulation);

  //---------- UTILITIES
  Anchor & anchor(Dart_descriptor const dart);
  Anchor const & anchor(Dart_const_descriptor const dart) const;
  Anchor & anchor();
  Anchor const & anchor() const;
  unsigned index_in_anchor(Dart_const_descriptor const dart) const;
  Dart_descriptor ith_dart(unsigned i, Anchor const & anch);
  void display_vertices(Anchor const & anch, bool round = true) const;
  bool is_valid() const;
  void to_stream(std::ostream & s) const;
  void from_stream(std::istream & s);

  //---------- location and insertion
  Locate_type relative_position(Point const & query, unsigned & li, Anchor const & anch) const;
  Anchor locate(Point const & query, bool use_visibility = false); // const ?
  Anchor locate(Point const & query, Locate_type & lt, unsigned & li, unsigned & ld, Anchor const & hint, bool use_visibility = false); // const ?

  //---------- Delaunay related methods
  void insert(Point const & query, Anchor & hint);
  void insert(Point const & query);

  //---------- eps-net methods
  bool epsilon_net(double const epsilon, unsigned const p = 0);
  bool is_epsilon_covering(const double epsilon) const;
  bool is_epsilon_packing(const double epsilon) const;
  bool is_epsilon_net(const double epsilon) const;

  double shortest_loop() const;
  double shortest_edge() const;
}

} // namespace CGAL
