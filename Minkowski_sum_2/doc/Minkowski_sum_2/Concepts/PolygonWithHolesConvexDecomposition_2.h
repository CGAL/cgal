/*!
\ingroup PkgMinkowskiSum2Concepts
\cgalConcept

A model of the `PolygonWithHolesConvexDecomposition_2` concept is capable of
decomposing an input polygon \f$ P\f$, which may have holes, into a set of
convex sub-polygons \f$ P_1, \ldots, P_k\f$, such that
\f$ \cup_{i=1}^{k}{P_k} = P\f$.

\cgalRefines `PolygonConvexDecomposition_2`

\cgalHasModel `CGAL::Polygon_vertical_decomposition_2<Kernel,Container>`
\cgalHasModel `CGAL::Polygon_triangulation_decomposition_2<Kernel,Container>`

*/

class PolygonWithHolesConvexDecomposition_2 {
public:

  /// \name Types
  /// @{

  /*! the polygon with holes type, defined as
   * `Polygon_with_holes_2<Kernel,Container>`.
   */
  typedef unspecified_type Polygon_with_holes_2;

  /// @}

  /// \name Creation
  /// @{

  /*! default constructor. */
  PolygonWithHolesConvexDecomposition_2 ();

  /// @}

  /// \name Operations
  /// @{

  /*! decomposes the input polygon with holes `P` into convex sub-polygons,
   * and writes them to the output iterator `oi`. The value-type of the
   * output iterator must be `Polygon_2`.
   * The function returns a past-the-end iterator for the convex sub-polygons.
   */
  template <class OutputIterator>
  OutputIterator operator()(const Polygon_with_holes_2& P,
                            OutputIterator oi) const;

  /// @}

}; /* end PolygonWithHolesConvexDecomposition_2 */
