namespace CGAL{

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2MainClasses

Represents a fundamental domain of a closed orientable hyperbolic surface.
The domain is given as a polygon \f$ P  \f$ represented by the list of its vertices in the Poincaré disk model,
together with a pairing of the sides of \f$ P \f$.
The \f$ n \f$-th side of \f$ P \f$ is the side between the \f$ n \f$-th and the \f$ (n+1) \f$-th vertex, where indices are modulo the number of vertices of \f$ P \f$.
The sides pairing are represented by a list of integers, such that if the \f$ n \f$-th integer of the list is \f$ m \f$, then the \f$ n \f$-th side is paired to the \f$ m \f$-th side.

\tparam Traits is the traits class and must be a model of `HyperbolicSurfacesTraits_2` (default model: `Hyperbolic_surface_traits_2`).
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
    Constructor from vertices and side pairings interators.
  */
   Hyperbolic_fundamental_domain_2(typename std::vector<Point>::iterator vfirst,
				  typename std::vector<Point>::iterator vlast,
				  typename std::vector<int>::iterator pfirst,
				  typename std::vector<int>::iterator plast);
  /// @}

  /// \name Access functions
  /// @{
  /*!
      returns the number of vertices (equivalently, the number of sides) of the domain.

      \pre <code> is_valid() </code>
  */
  int size() const;

  /*!
      returns the \f$ i \f$-th vertex.

      \pre <code> is_valid() </code>
  */
  const Point& vertex(int i) const;

  /*!
      returns the index of the side paired to the \f$ i \f$-th side.

      \pre <code> is_valid() </code>
  */
  int paired_side(int i) const;

  /*!
       returns the isometry that maps side \f$ \overline A \f$ to side \f$ A
       \f$, where \f$ A \f$ is the \f$ i \f$-th side, and \f$ \overline A \f$ is the side paired to \f$ A \f$.

       \pre <code> is_valid() </code>
  */
  Hyperbolic_isometry_2<Traits> side_pairing(int i) const;
  /// @}

  /// \name Input/output
  /// @{
  /*!
      Reads the domain from a stream.

      The format of the input should be the same as the format of the output of the 'from_stream' method.

      \pre <code> is_valid() </code>
  */
  std::istream& operator>>(std::istream& s, Hyperbolic_fundamental_domain_2<Traits>& domain);

  /*!
      writes the domain in a stream.

      The format of the output is the following.
      The first line prints the number n of vertices of the domain.
      For \f$ i=0 \f$ to \f$ n-1 \f$ the index of the side paired to side \f$ i \f$ is printed on a separate line.
      For \f$ i=0 \f$ to \f$ n-1 \f$ the \f$ i \f$-th vertex is printed on a separate line.

      \pre <code> is_valid() </code>
  */
  std::ostream& operator<<(std::ostream& s, const Hyperbolic_fundamental_domain_2<Traits>& domain);
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
