// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Represents a convex geodesic fundamental domain of a closed orientable hyperbolic surface.
The domain is given as a polygon P given by the list of its vertices in the Poincaré disk model,
together with a pairing of the sides of P.
The n-th side of P is the side between the n-th and the (n+1)-th vertex, where indices are modulo the number of vertices of P.
The sides pairing are represented by a list of integers, such that if the n-th integer of the list is m, then the n-th side is paired to the m-th side

\tparam GeometricTraits_2 is the geometric traits class and must be a model of `HyperbolicSurfacesTraits_2`.
*/
template<class GeometricTraits_2>
class Hyperbolic_fundamental_domain_2 {
private:
  typedef typename GeometricTraits_2::Point_2                       _Point;
  typedef Hyperbolic_isometry_2<GeometricTraits_2>                  _Isometry;

  std::vector<_Point>                                               _vertices;
  std::vector<int>                                                  _pairings;

public:
  /// \name Types
  /// @{
  typedef GeometricTraits_2                                         Geometric_traits_2;
  /// @}

  /// \name Creation
  /// @{
  /*!
    Default constructor
  */
  Hyperbolic_fundamental_domain_2();

  /*!
    Constructor from the list of vertices and the sides pairing.
  */
  Hyperbolic_fundamental_domain_2(const std::vector<typename GeometricTraits_2::Point_2>& vertices, const std::vector<int>& pairings);
  /// @}

  /// \name Get and set
  /// @{
  /*!
    Sets the vertices and the sides pairings of the domain.
  */
  void set(const std::vector<typename GeometricTraits_2::Point_2>& vertices, const std::vector<int>& pairings);

  /*!
      Returns the number of vertices (equivalently, the number of sides) of the domain.
  */
  int size() const;

  /*!
      Returns the index-th vertex
  */
  typename GeometricTraits_2::Point_2 vertex(int index) const;

  /*!
      Returns the index of the side paired to side A, where A is the index-th side
  */
  int paired_side(int index) const;

  /*!
      // Returns the isometry that maps side A to side B, where B is the index-th side, and A is the side paired to B
  */
  Hyperbolic_isometry_2<GeometricTraits_2> side_pairing(int index) const;
  /// @}

  /// \name Stream input output
  /// @{
  /*!
      Reads the domain from a stream
  */
  void from_stream(std::istream& s);

  /*!
      Writes the domain in a stream
  */
  void to_stream(std::ostream& s) const;
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
