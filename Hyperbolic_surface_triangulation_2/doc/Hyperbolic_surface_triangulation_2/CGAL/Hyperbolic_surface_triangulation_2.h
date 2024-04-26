// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Represents a triangulation of a closed orientable hyperbolic surface.

Offers facilities such as the generation of the triangulation from a convex fundamental domain,
the Delaunay flip algorithm, and the construction of a portion of the lift of the triangulation in the hyperbolic plane.

\tparam Traits is the traits class and must be a model of `HyperbolicSurfacesTraits`.
*/
template<class Traits>
class Hyperbolic_surface_triangulation_2{
public:
  /// \name Types
  /// @{
  /*!
      Type of combinatorial map whose edges are attributed complex numbers.
  */
  typedef Combinatorial_map<2,unspecified_type>                                 Combinatorial_map_with_cross_ratios;
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
      Stores a dart d of the combinatorial map, belonging to a triangle t, and stores the three vertices of a lift of t in the hyperbolic plane.
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
  Hyperbolic_surface_triangulation_2(const Combinatorial_map_with_cross_ratios& cmap);

  /*!
      Constructor from a decorated combinatorial map and an anchor.
  */
  Hyperbolic_surface_triangulation_2(const Combinatorial_map_with_cross_ratios& cmap, const Anchor& anchor);
  /// @}

  /// \name Assignment
  /// @{
  /*!
      Copy constructor.
  */
  Hyperbolic_surface_triangulation_2& operator=(Hyperbolic_surface_triangulation_2 other);
  /// @}

  /// \name Access functions
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
  bool is_delaunay_flippable(Dart_handle dart) const;

  /*!
      Flips the edge supported by the dart.
  */
  void flip(Dart_handle dart);

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
      If the center parameter is set to true, then one of the lifted triangles admits the point 0 as a vertex.

      This method is to be used only if the triangulation has an anchor.
  */
  std::vector<std::tuple<Dart_const_handle, Point, Point, Point>> lift(bool center=true) const;
  /// @}

  /// \name Validity
  /// @{
  /*!
      Validity test.

      Checks that the combinatorial map has no 1,2-boundary and calls the is_valid method of the combinatorial map.
      If there is an anchor, then checks that the dart handle of the anchor does indeed point to a dart of the combinatorial map, and checks that the three vertices of the anchor lie within the open unit disk.
  */
  bool is_valid() const;
  /// @}

  /// \name Input/output
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
  template<class Traits> std::ostream& operator<<(std::ostream& s, Hyperbolic_surface_triangulation_2<Traits>& Hyperbolic_surface_triangulation_2);

  /*!
      Reads the triangulation from a stream.

      The format of the input should be the same as the format of the output of the 'from_stream' method.
  */
  template<class Traits> void operator>>(std::istream& s, Hyperbolic_surface_triangulation_2<Traits>& triangulation);
  /// @}
};

} //namespace CGAL
