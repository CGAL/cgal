/*!
\ingroup PkgTSMAConcepts
\cgalConcept

The function object class to fit a plane from a range of facets.

\cgalHasModel `PlaneFitting`
\cgalHasModel `PCAPlaneFitting`

*/
class PlaneFitting {
  typedef unspecified_type Proxy;

  template<FacetsIterator>
  Plane operator()(const FacetsIterator &beg, const FacetsIterator &end) const;
};

