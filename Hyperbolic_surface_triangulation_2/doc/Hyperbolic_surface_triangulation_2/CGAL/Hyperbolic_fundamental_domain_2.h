// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Represents a fundamental domain of a closed orientable hyperbolic surface, the domain is convex and its edges are geodesic segments.
The domain is given as a polygon P represented by the list of its vertices in the Poincaré disk model,
together with a pairing of the sides of P.
The n-th side of P is the side between the n-th and the (n+1)-th vertex, where indices are modulo the number of vertices of P.
The sides pairing are represented by a list of integers, such that if the n-th integer of the list is m, then the n-th side is paired to the m-th side.

\tparam Traits_2 is the traits class and must be a model of `HyperbolicSurfacesTraits_2`.
*/
template<class Traits_2>
class Hyperbolic_fundamental_domain_2 {
public:
  /// \name Creation
  /// @{
  /*!
    Default constructor
  */
  Hyperbolic_fundamental_domain_2();

  /*!
    Constructor from the list of vertices and the sides pairing.
  */
  Hyperbolic_fundamental_domain_2(const std::vector<Traits_2::Hyperbolic_point_2>& vertices, const std::vector<int>& pairings);
  /// @}

  /// \name Get and set
  /// @{
  /*!
    Sets the vertices and the sides pairings of the domain.
  */
  void set(const std::vector<Traits_2::Hyperbolic_point_2>& vertices, const std::vector<int>& pairings);

  /*!
      Returns the number of vertices (equivalently, the number of sides) of the domain.
  */
  int size() const;

  /*!
      Returns the index-th vertex
  */
  typename Traits_2::Hyperbolic_point_2 vertex(int index) const;

  /*!
      Returns the index of the side paired to side A, where A is the index-th side
  */
  int paired_side(int index) const;

  /*!
       Returns the isometry that maps side A to side B, where B is the index-th side, and A is the side paired to B
  */
  Hyperbolic_isometry_2<Traits_2> side_pairing(int index) const;
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
