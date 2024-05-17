// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

namespace CGAL{

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Represents a fundamental domain of a closed orientable hyperbolic surface.
The domain is given as a polygon \f$ P  \f$ represented by the list of its vertices in the Poincaré disk model,
together with a pairing of the sides of \f$ P \f$.
The \f$ n \f$-th side of \f$ P \f$ is the side between the \f$ n \f$-th and the \f$ (n+1) \f$-th vertex, where indices are modulo the number of vertices of \f$ P \f$.
The sides pairing are represented by a list of integers, such that if the \f$ n \f$-th integer of the list is \f$ m \f$, then the \f$ n \f$-th side is paired to the \f$ m \f$-th side.

\tparam Traits is the traits class and must be a model of `HyperbolicSurfacesTraits_2`.
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
      Returns the index-th vertex.
  */
  Point vertex(int index) const;

  /*!
      Returns the index of the side paired to the index-th side.
  */
  int paired_side(int index) const;

  /*!
       Returns the isometry that maps side \f$ A \f$ to side \f$ B \f$, where \f$ B \f$ is the index-th side, and \f$ A \f$ is the side paired to \f$ B \f$.
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
      For \f$ i=0 \f$ to \f$ n-1 \f$ the index of the side paired to side \f$ i \f$ is printed on a separate line.
      For \f$ i=0 \f$ to \f$ n-1 \f$ the \f$ i \f$-th vertex is printed on a separate line.
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

}; // namespace CGAL
