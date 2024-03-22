// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

The class `Hyperbolic_surface_triangulation_2` is the main class of the Hyperbolic Surface Triangulation package.
It represents a triangulation of a closed orientable hyperbolic surface.
It offers factilities such as the generation of the triangulation from a convex fundamental domain,
the Delaunay flip algorithm, and the construction of a portion of the lift of the triangulation in the Poincaré disk model.

\tparam GeometricTraits_2 is the geometric traits class and must be a model of `HyperbolicSurfacesTraits_2`.
*/
template<class GeometricTraits_2>
class Hyperbolic_surface_triangulation_2{
public:
  /// \name Types
  /// @{
  typedef GeometricTraits_2                                                                             Geometric_traits_2;

  struct Combinatorial_map_with_cross_ratios_item{
    template <class CMap>
    struct Dart_wrapper{
        typedef Cell_attribute<CMap, Complex_without_sqrt<typename GeometricTraits_2::FT>> Edge_attrib;
        typedef std::tuple<void,Edge_attrib,void>   Attributes;
    };
  };
  typedef Combinatorial_map<2,Combinatorial_map_with_cross_ratios_item>                                 Combinatorial_map_with_cross_ratios;

  struct Anchor{
    typename Combinatorial_map_with_cross_ratios::Dart_handle dart;
    typename GeometricTraits_2::Point_2 vertices[3];
  };
  /// @}

private:
  typedef typename Combinatorial_map_with_cross_ratios::Dart_handle                                     _Dart_handle;
  typedef typename Combinatorial_map_with_cross_ratios::Dart_range                                      _Dart_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_range<0>             _Vertex_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_range<1>             _Edge_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_range<2>             _Face_range;

  typedef typename Combinatorial_map_with_cross_ratios::Dart_const_handle                               _Dart_const_handle;
  typedef typename Combinatorial_map_with_cross_ratios::Dart_const_range                                _Dart_const_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_const_range<1>       _Edge_const_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_const_range<2>       _Face_const_range;

  typedef typename GeometricTraits_2::FT                                                                _Number;
  typedef Complex_without_sqrt<_Number>                                                                 _Complex_number;
  typedef typename GeometricTraits_2::Point_2                                                           _Point;
  typedef Hyperbolic_isometry_2<GeometricTraits_2>                                                      _Isometry;
  typedef Hyperbolic_fundamental_domain_2<GeometricTraits_2>                                            _Domain;

public:
  /// \name Creation
  /// @{
  /*!
      Default constructor
  */
  Hyperbolic_surface_triangulation_2() {};

  /*!
      Constructor from a convex fundamental domain
  */
  Hyperbolic_surface_triangulation_2(const Hyperbolic_fundamental_domain_2<GeometricTraits_2>& domain);

  /*!
      Constructor from a decorated combinatorial map
  */
  Hyperbolic_surface_triangulation_2(const Combinatorial_map_with_cross_ratios& cmap);

  /*!
      Constructor from a decorated combinatorial map and an anchor
  */
  Hyperbolic_surface_triangulation_2(const Combinatorial_map_with_cross_ratios& cmap, const Anchor& anchor);

  /*!
      Copy constructor
  */
  Hyperbolic_surface_triangulation_2& operator=(Hyperbolic_surface_triangulation_2 other);
  /// @}

  /// \name Get and set
  /// @{
  /*!
      Returns the decorated combinatorial map
  */
  Combinatorial_map_with_cross_ratios& get_combinatorial_map_ref();

  /*!
      Tells whether the triangulation is anchored
  */
  bool has_anchor() const;

  /*!
      Returns the anchor
  */
  Anchor& get_anchor_ref();
  /// @}

  /// \name Stream input output
  /// @{
  /*!
      Writes the triangulation in a stream
  */
  void to_stream(std::ostream& s) const;

  /*!
      Reads the triangulation from a stream
  */
  void from_stream(std::istream& s);
  /// @}

  /// \name Delaunay flips algorithm
  /// @{
  /*!
      Tells whether an edge is Delaunay flippable
  */
  bool is_delaunay_flippable(_Dart_handle dart) const;

  /*!
      Flips a Delaunay flippable edge
  */
  void flip(_Dart_handle dart);

  /*!
      Applies the Delaunay flips algorithm
  */
  int make_delaunay();
  /// @}

  /// \name Lifting
  /// @{
  /*!
       Returns, for every triangle T of the triangulation, one of the darts of T together with a triple A,B,C of points in the hyperbolic plane.
       The points A,B,C are the vertices of a a lift of T in the hyperbolic plane.
       This method is to be used only if the triangulation has an anchor.
  */
  std::vector<std::tuple<typename Combinatorial_map_with_cross_ratios::Dart_const_handle,typename GeometricTraits_2::Point_2,typename GeometricTraits_2::Point_2,typename GeometricTraits_2::Point_2>> lift(bool center=true) const;
  /// @}

  /// \name Validity
  /// @{
  /*!
      Partial validity test
  */
  bool is_valid() const;
  /// @}
};

} //namespace CGAL
