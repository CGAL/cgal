/*!
\ingroup PkgTSMAConcepts
\cgalConcept

The concept of `PlaneFitting` is a function object class that fit a plane from a range of facets.
It is used in the functions `vsa_extract()`, `vsa_approximate_and_extract()`.

\cgalHasModel `PlaneFitting`
\cgalHasModel `PCAPlaneFitting`
*/

class PlaneFitting {
public:
  /// 3D plane type
  typedef unspecified_type Plane;

  /// @name Operations
  /// A model of this concept must provide:
  /// @{

  /// returns the fitting plane from a range of facets.
  template<FacetsIterator>
  Plane operator()(const FacetsIterator &beg, const FacetsIterator &end) const;

  /// }
};
