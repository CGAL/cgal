/*!
\ingroup PkgTSMAConcepts
\cgalConcept

The concept `FittingProxyErrorMetric` defines the notion of 
proxy and computes the fitting error from a face to a `Proxy`,
used in CGAL::VSA_approximation.

\cgalHasModel `CGAL::L2_metric_plane_proxy`
\cgalHasModel `CGAL::L21_metric_plane_proxy`
*/

class FittingProxyErrorMetric {
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

  /// returns fitting error from face f to proxy.
  FT compute_error(const face_descriptor &f, const Proxy &proxy) const;

  /// returns fitted proxy for a range of facets.
  template <typename FacetIterator>
  Proxy fit_proxy(const FacetIterator &begin, const FacetIterator &end) const;

  /// }

};
