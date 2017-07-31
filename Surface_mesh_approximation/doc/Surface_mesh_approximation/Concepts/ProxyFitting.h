/*!
\ingroup PkgTSMAConcepts
\cgalConcept

The function object class to fit a proxy from a range of facets.

\cgalHasModel `L21ProxyFitting`
\cgalHasModel `L2ProxyFitting`

*/
class ProxyFitting {
  typedef unspecified_type Proxy;

  template<FacetsIterator>
  Proxy operator()(const FacetsIterator &beg, const FacetsIterator &end) const;
};
