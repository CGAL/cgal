/*!
\ingroup PkgTSMAConcepts
\cgalConcept

The concept `ErrorMetric` computes the fitting error from
 a face to a `Proxy`, used in CGAL::VSA_approximation.

\cgalHasModel `CGAL::L21_metric`
\cgalHasModel `CGAL::L2_metric`
*/

class ErrorMetric {
public:
  /// Number type model of `Field` and `RealEmbeddable`
  typedef unspecified_type FT;
  /// Triangle mesh face descriptor
  typedef unspecified_type face_descriptor;
  /// Parameterized shape proxy
  typedef unspecified_type Proxy;

  /// @name Operations
  /// A model of this concept must provide:
  /// @{

  /// returns the fitting error from face f to Proxy px.
  FT compute_error(const face_descriptor &f, const Proxy &px) const;

  /// returns the fitted proxy for a range of facets.
  template <typename FacetIterator>
  Proxy fit_proxy(const FacetIterator &begin, const FacetIterator &end) const;

  /// }

};
