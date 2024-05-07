// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Represents a fundamental domain of a closed orientable hyperbolic surface.
The domain is given as a polygon P represented by the list of its vertices in the Poincaré disk model,
together with a pairing of the sides of P.
The n-th side of P is the side between the n-th and the (n+1)-th vertex, where indices are modulo the number of vertices of P.
The sides pairing are represented by a list of integers, such that if the n-th integer of the list is m, then the n-th side is paired to the m-th side.

\tparam Traits is the traits class and must be a model of `HyperbolicSurfacesTraits`.
*/
template<class Traits>
class Hyperbolic_fundamental_domain_2 {
public:
  /// \name Types
  /// @{
  /*!
  Point type.
  */
  typedef typename Traits::Hyperbolic_point_2                    Point;
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
  Hyperbolic_fundamental_domain_2(const std::vector<Point>& vertices, const std::vector<int>& pairings);
  /// @}

  /// \name Access functions
  /// @{
  /*!
    Sets the vertices and the sides pairings of the domain.
  */
  void set(const std::vector<Point>& vertices, const std::vector<int>& pairings);

  /*!
      Returns the number of vertices (equivalently, the number of sides) of the domain.
  */
  int size() const;

  /*!
      Returns the index-th vertex
  */
  Point vertex(int index) const;

  /*!
      Returns the index of the side paired to side A, where A is the index-th side
  */
  int paired_side(int index) const;

  /*!
       Returns the isometry that maps side A to side B, where B is the index-th side, and A is the side paired to B
  */
  Hyperbolic_isometry_2<Traits> side_pairing(int index) const;
  /// @}

  /// \name Input/output
  /// @{
  /*!
      Reads the domain from a stream.

      The format of the input should be the same as the format of the output of the 'from_stream' method.
  */
  template<class Traits> void operator>>(std::istream& s, Hyperbolic_fundamental_domain_2<Traits>& domain);

  /*!
      Writes the domain in a stream.

      The format of the output is the following.
      The first line prints the number n of vertices of the domain.
      For i=0 to n-1 the index of the side paired to side i is printed on a separate line.
      For i=0 to n-1 the i-th vertex is printed on a separate line.
  */
  template<class Traits> std::ostream& operator<<(std::ostream& s, const Hyperbolic_fundamental_domain_2<Traits>& domain);
  /// @}

  /// @{
  /// \name Validity
  /*!
      Validity test.

      Checks that the number of vertices is even, that there are as many side pairings as vertices, and that the vertices all lie within the open unit disk.
  */
  bool is_valid() const;
  /// @}
};

} //namespace CGAL
