/*!
\ingroup PkgTSMAConcepts
\cgalConcept

The concept `ProxyFitting` is required to fit a shape `Proxy` from a range of facets.
It is used in the functions `vsa_mesh_approximation()`.

\cgalHasModel `L21ProxyFitting`
\cgalHasModel `L2ProxyFitting`
*/

class ProxyFitting {
public:
  /// Parametrized shape proxy.
  typedef unspecified_type Proxy;

  /// @name Operations
  /// A model of this concept must provide:
  /// @{

  /// returns the fitting proxy from a range of facets.
  template<FacetsIterator>
  Proxy operator()(const FacetsIterator &beg, const FacetsIterator &end) const;

  /// }
};
