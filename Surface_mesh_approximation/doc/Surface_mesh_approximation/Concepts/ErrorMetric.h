
/*!
\ingroup PkgTSMAConcepts
\cgalConcept

This is a functor to compute the fitting error of a facet to a proxy.

\cgalHasModel `L21Metric`
\cgalHasModel `L2Metric`

*/
class ErrorMetric {
public:
  /*!
   * A number type model of `Field` and `RealEmbeddable`
   */
  typedef unspecified_type FT;
  /// Triangle mesh facet descriptor.
  typedef unspecified_type facet_descriptor;
  /// Parametrized shape proxy.
  typedef unspecified_type Proxy;

  /// @name Functions
  /// @{
  /// data member to describe the proxy seed.
  FT operator()(const facet_descriptor &f, const Proxy &px);
  /// }
};

