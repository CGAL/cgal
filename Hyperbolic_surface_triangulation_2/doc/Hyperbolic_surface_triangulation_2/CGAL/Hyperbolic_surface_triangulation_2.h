namespace CGAL{

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses


This item defines attributes of edges that are
`Complex_without_sqrt` reprensenting cross-ratios.

\tparam Traits is the traits class and must be a model of `HyperbolicSurfacesTraits_2` (default model: `Hyperbolic_surface_traits_2`).

\cgalModels{GenericMapItems}
*/
template<class Traits>
struct Combinatorial_map_with_cross_ratios_item{
    template <class CMap>
    struct Dart_wrapper{
        typedef Cell_attribute<CMap, Complex_without_sqrt<typename Traits::FT>> Edge_attrib;
        typedef std::tuple<void,Edge_attrib,void>   Attributes;
    };
  };


/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Represents a triangulation of a closed orientable hyperbolic surface.

Offers facilities such as the generation of the triangulation from a convex fundamental domain,
the Delaunay flip algorithm, and the construction of a portion of the lift of the triangulation in the hyperbolic plane.

\tparam Traits is the traits class and must be a model of `HyperbolicSurfacesTraits_2` (default model: `Hyperbolic_surface_traits_2`).

\tparam Attributes must be a model of `GenericMapItems` (default model: `Combinatorial_map_with_cross_ratios_item<Traits>`).
*/
template<class Traits, class Attributes = Combinatorial_map_with_cross_ratios_item<Traits>>
class Hyperbolic_surface_triangulation_2{
public:
  /// \name Types
  /// @{
  /*!
      Type of combinatorial map whose edges are decorated with complex numbers.
  */
  typedef Combinatorial_map<2, Attributes>                                 Combinatorial_map_with_cross_ratios;
  /*!
      Combinatorial map dart handle type.
  */
  typedef typename Combinatorial_map_with_cross_ratios::Dart_handle                                     Dart_handle;
  /*!
      Combinatorial map dart const handle type.
  */
  typedef typename Combinatorial_map_with_cross_ratios::Dart_const_handle                               Dart_const_handle;
  /*!
      Point type.
  */
  typedef typename Traits::Hyperbolic_point_2                                                           Point;

  /*!
      stores a dart \f$ d \f$ of the combinatorial map, belonging to a triangle \f$ t \f$, and stores the three vertices of a lift of \f$ t \f$ in the hyperbolic plane.
  */
  struct Anchor{
    typename Combinatorial_map_with_cross_ratios::Dart_handle dart;
    typename Traits::Hyperbolic_point_2 vertices[3];
  };
  /// @}

  /// \name Creation
  /// @{
  /*!
      Default constructor.
  */
  Hyperbolic_surface_triangulation_2() {};

  /*!
      Constructor from a convex fundamental domain: triangulates the polygon of the domain, and identifies the paired sides of the domain.
  */
  Hyperbolic_surface_triangulation_2(const Hyperbolic_fundamental_domain_2<Traits>& domain);

  /*!
      Constructor from a decorated combinatorial map.
  */
  Hyperbolic_surface_triangulation_2(Combinatorial_map_with_cross_ratios& cmap);

  /*!
      Constructor from a decorated combinatorial map and an anchor.
  */
  Hyperbolic_surface_triangulation_2(Combinatorial_map_with_cross_ratios& cmap,  Anchor& anchor);
  /// @}

  /// \name Assignment
  /// @{
  /*!
      \pre <code> other.is_valid() </code>
  */
  Hyperbolic_surface_triangulation_2& operator=(Hyperbolic_surface_triangulation_2 other);
  /// @}

  /// \name Access functions
  /// @{
  /*!
      returns the decorated combinatorial map.

      \pre <code> is_valid() </code>
  */
  Combinatorial_map_with_cross_ratios& combinatorial_map();

  /*!
      tells if the triangulation has an anchor.

      \pre <code> is_valid() </code>
  */
  bool has_anchor() const;

  /*!
      returns the anchor.

      \pre <code> is_valid() && has_anchor() </code>
  */
  Anchor& anchor();
  /// @}

  /// \name Delaunay flip algorithm
  /// @{
  /*!
      tells if if the edge supported by the dart is Delaunay flippable.

      \pre <code> is_valid() </code>
  */
  bool is_Delaunay_flippable(Dart_handle dart) const;

  /*!
      flips the edge supported by the dart.

      \pre <code> is_valid() </code>
  */
  void flip(Dart_handle dart);

  /*!
      determines if the triangulation is a valid Delaunay triangulation.
  */
  bool is_Delaunay() const;

  /*!
      applies the Delaunay flip algorithm: flips Delaunay flippable edges until there is no such edge anymore.

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
      If the center parameter is set to true, then one of the lifted triangles admits the origin \f$ 0 \f$ as a vertex.

      \pre <code> is_valid() && has_anchor() </code>
  */
  std::vector<std::tuple<Dart_const_handle, Point, Point, Point>> lift(bool center=true) const;
  /// @}

  /// \name Validity
  /// @{
  /*!
      Validity test.
      Checks that the underlying combinatorial map \f$ M \f$ has no boundary and calls the is_valid method of \f$ M \f$.
      If there is an anchor, then checks that the dart handle of the anchor does indeed point to a dart of \f$ M \f$, and checks that the three vertices of the anchor lie within the open unit disk.
  */
  bool is_valid() const;
  /// @}

  /// \name Input/output
  /// @{
  /*!
      writes the triangulation in a stream.
      The format of the output is the following.
      Each dart of the triangulation is given an index between \f$ 0 \f$ and \f$ n-1 \f$, where \f$ n \f$ is the number of darts of the triangulation.
      The first line contains the number \f$ n \f$ of darts.
      The next line contains either 'yes' or 'no' and tells whether the triangulation has an anchor.
      If the triangulation has an anchor, then the four next lines print the index of the dart of the anchor, and the three vertices of the anchor.
      Then, for every triangle \f$ t \f$, the indices of the three darts of \f$ t \f$ are printed on three distinct lines.
      Finally, for every edge \f$ e \f$, the indices of the two darts of \f$ e \f$ are printed on two distinct lines, followed by a third line on which the cross ratio of \f$ e \f$ is printed.

      \pre <code> is_valid() </code>
  */
  std::ostream& operator<<(std::ostream& s, Hyperbolic_surface_triangulation_2<Traits>& Hyperbolic_surface_triangulation_2);

  /*!
      reads the triangulation from a stream.
      The format of the input should be the same as the format of the output of the '<<' operator.

      \pre <code> is_valid() </code>
  */
  std::istream& operator>>(std::istream& s, Hyperbolic_surface_triangulation_2<Traits>& triangulation);
  /// @}
};

}; // namespace CGAL
