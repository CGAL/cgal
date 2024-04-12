// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Represents a triangulation of a closed orientable hyperbolic surface.

Offers facilities such as the generation of the triangulation from a convex fundamental domain,
the Delaunay flip algorithm, and the construction of a portion of the lift of the triangulation in the hyperbolic plane.

The triangulation T is represented internally by an instance of CGAL::Combinatorial_map whose edges are attributed complex numbers. The complex number attributed to an edge e of T is the cross ratio of e in T, it can be defined as follows. Consider a lift f of e in the Poincaré disk model of the hyperbolic plane.
Orient f, and let z0 and z2 be the two complex numbers representing the first and second vertices of f.
Consider the lifted triangle of T that lies on the right of f, and let z1 be the complex number representing the third vertex of this triangle (the vertex distinct from z0 and z2).
Similarly, consider the lifted triangle of T that lies on the left of f, and let z3 be the complex number representing the third vertex of this triangle.
The cross ratio of e in T is then the complex number  (z3-z1)*(z2-z0) / ((z3-z0)*(z2-z1)). This definition does not depend on the choice of the lift f, nor on the orientation of f.

The triangulation T can optionnally contain some additional data: the anchor.
The anchor contains a lift t of a triangle of T in the hyperbolic plane: t is represented by its three vertices in the Poincaré disk model, and by a dart of the corresponding triangle in the combinatorial map of T.
While the combinatorial map and its cross ratios uniquely determine T, the anchor is used when building a portion of the lift of T in the hyperbolic plane.


\tparam Traits_2 is the traits class and must be a model of `HyperbolicSurfacesTraits_2`.
*/
template<class Traits_2>
class Hyperbolic_surface_triangulation_2{
public:
  /// \name Types
  /// @{
  struct Combinatorial_map_with_cross_ratios_item{
    template <class CMap>
    struct Dart_wrapper{
        typedef Cell_attribute<CMap, Complex_without_sqrt<typename Traits_2::FT>> Edge_attrib;
        typedef std::tuple<void,Edge_attrib,void>   Attributes;
    };
  };
  /*!
      Type of combinatorial map whose edges are attributed complex numbers.
  */
  typedef Combinatorial_map<2,Combinatorial_map_with_cross_ratios_item>                                 Combinatorial_map_with_cross_ratios;

  /*!
      Stores a dart d of the combinatorial map, belonging to a triangle t, and stores the three vertices of a lift of t in the hyperbolic plane.
  */
  struct Anchor{
    typename Combinatorial_map_with_cross_ratios::Dart_handle dart;
    typename Traits_2::Hyperbolic_point_2 vertices[3];
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
  Hyperbolic_surface_triangulation_2(const Hyperbolic_fundamental_domain_2<Traits_2>& domain);

  /*!
      Constructor from a decorated combinatorial map.
  */
  Hyperbolic_surface_triangulation_2(const Combinatorial_map_with_cross_ratios& cmap);

  /*!
      Constructor from a decorated combinatorial map and an anchor.
  */
  Hyperbolic_surface_triangulation_2(const Combinatorial_map_with_cross_ratios& cmap, const Anchor& anchor);

  /*!
      Copy constructor.
  */
  Hyperbolic_surface_triangulation_2& operator=(Hyperbolic_surface_triangulation_2 other);
  /// @}

  /// \name Get and set
  /// @{
  /*!
      Returns the decorated combinatorial map.
  */
  Combinatorial_map_with_cross_ratios& get_combinatorial_map_ref();

  /*!
      Tells if the triangulation has an anchor.
  */
  bool has_anchor() const;

  /*!
      Returns the anchor (if there is one).
  */
  Anchor& get_anchor_ref();
  /// @}

  /// \name Delaunay flip algorithm
  /// @{
  /*!
      Tells if if the edge supported by the dart is Delaunay flippable.
  */
  bool is_delaunay_flippable(Combinatorial_map_with_cross_ratios::Dart_handle dart) const;

  /*!
      Flips the edge supported by the dart.
  */
  void flip(Combinatorial_map_with_cross_ratios::Dart_handle dart);

  /*!
      Applies the Delaunay flip algorithm: flips Delaunay flippable edges until there is no such edge anymore.
  */
  int make_delaunay();
  /// @}

  /// \name Lifting
  /// @{
  /*!
      Lifts the triangulation in the hyperbolic plane. Returns, for every triangle T of the triangulation, one of the darts of T in the combinatorial map of the triangulation, together with a triple A,B,C of points in the hyperbolic plane.
      The points A,B,C are the vertices of a lift of T in the hyperbolic plane.
      This method is to be used only if the triangulation has an anchor.
  */
  std::vector<std::tuple<Combinatorial_map_with_cross_ratios::Dart_const_handle, Traits_2::Hyperbolic_point_2, Traits_2::Hyperbolic_point_2, Traits_2::Hyperbolic_point_2>> lift(bool center=true) const;
  /// @}

  /// \name Validity
  /// @{
  /*!
      Partial validity test.
  */
  bool is_valid() const;
  /// @}

  /// \name Stream input output
  /// @{
  /*!
      Writes the triangulation in a stream.

      The format of the output is the following.
      Each dart of the triangulation is given an index between 0 and n-1, where n is the number of darts of the triangulation.
      The first line contains the number n of darts.
      For every triangle t, the indices of the three darts of t are printed on three distinct lines.
      For every edge e, the indices of the two darts of e are printed on two distinct lines, followed by a third line on which the cross ratio of e is printed.
      The next line contains either 'yes' or 'no' and tells whether the triangulation has an anchor.
      If the triangulation has anchor A, then the two next lines print respectively the index of the dart of A, and the three complex numbers of A.
  */
  void to_stream(std::ostream& s) const;

  /*!
      Reads the triangulation from a stream.

      The format of the input should be the same as the format of the output of the 'from_stream' method.
  */
  void from_stream(std::istream& s);
  /// @}
};

} //namespace CGAL
