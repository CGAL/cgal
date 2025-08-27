namespace CGAL {

/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2InputOutput

  inserts the triangulation in a stream.

  The format of the output is the following.
  Each dart of the triangulation is given an index between \f$ 0 \f$ and \f$ n-1 \f$, where \f$ n \f$
  is the number of darts of the triangulation.
  The first line contains the number \f$ n \f$ of darts.
  The next line contains either 'yes' or 'no' and tells whether the triangulation has an anchor.
  If the triangulation has an anchor, then the four next lines print the index of the dart of the anchor,
  and the three vertices of the anchor.
  Then, for every triangle \f$ t \f$, the indices of the three darts of \f$ t \f$ are printed on three distinct lines.
  Finally, for every edge \f$ e \f$, the indices of the two darts of \f$ e \f$ are printed on two distinct lines, followed by a third line on which the cross ratio of \f$ e \f$ is printed.

  \pre <code>   Triangulation_on_hyperbolic_surface_2<Traits>::is_valid() </code>
*/
std::ostream& operator<<(std::ostream& s, const Triangulation_on_hyperbolic_surface_2<Traits>& triangulation);

/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2InputOutput

  extracts the triangulation from a stream.

  The format of the input should be the same as the format of the output of
  the '<<' operator for Triangulation_on_hyperbolic_surface_2.
*/
std::istream& operator>>(std::istream& s, Triangulation_on_hyperbolic_surface_2<Traits>& triangulation);

/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2InputOutput

  inserts the domain in a stream.

  The format of the output is the following.
  The first line prints the number \f$n\f$ of vertices of the domain.
  For \f$ i=0 \f$ to \f$ n-1 \f$ the index of the side paired to side \f$ i \f$ is printed on a separate line.
  For \f$ i=0 \f$ to \f$ n-1 \f$ the i-th vertex is printed on a separate line.

  \pre <code> Hyperbolic_fundamental_domain_2< Traits >::is_valid() </code>
*/
std::ostream& operator<<(std::ostream& s, const Hyperbolic_fundamental_domain_2<Traits>& domain);

/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2InputOutput

  extracts the domain from a stream.

  The format of the input must be the same as the format of the output of
  the '<<' operator for Hyperbolic_fundamental_domain_2.
*/
std::istream& operator>>(std::istream& s, Hyperbolic_fundamental_domain_2<Traits>& domain);

/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2InputOutput

  inserts the isometry in a stream.
*/
std::ostream& operator<<(std::ostream& s, const Hyperbolic_isometry_2<Traits>& isometry);

} // namespace CGAL
