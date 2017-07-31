/*!
\ingroup PkgTSMAConcepts
\cgalConcept

The concept `ErrorMetric` is requied to operate on a facet and a proxy, returns the fitting error.
It is used in the functions `vsa_approximate()`, `vsa_extract()`, `vsa_approximate_and_extract()`.

\cgalHasModel `L21Metric`
\cgalHasModel `L2Metric`
*/

class ErrorMetric {
public:
  /// A number type model of `Field` and `RealEmbeddable`
  typedef unspecified_type FT;
  /// Triangle mesh facet descriptor.
  typedef unspecified_type facet_descriptor;
  /// Parametrized shape proxy.
  typedef unspecified_type Proxy;

  /// @name Operations
  /// A model of this concept must provide:
  /// @{

  /// returns the fitting error of a facet f to the proxy px.
  FT operator()(const facet_descriptor &f, const Proxy &px);

  /// }
};
