/*!
\ingroup PkgTSMAConcepts
\cgalConcept

The concept `ProxyFitting` fits a shape primitive `Proxy` from a range of facets, used in CGAL::VSA::Mesh_approximation.

\cgalHasModel `CGAL::VSA::L21_proxy_fitting`
\cgalHasModel `CGAL::VSA::L2_proxy_fitting`
*/

class ProxyFitting {
public:
  /// Parametrized shape proxy.
  typedef unspecified_type Proxy;

  /// @name Operations
  /// A model of this concept must provide:
  /// @{

  /// returns the fitted proxy from a range of facets.
  template <typename FacetIterator>
  Proxy operator()(const FacetIterator &beg, const FacetIterator &end) const;

  /// }
};
