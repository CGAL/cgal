/*!
\ingroup PkgTSMAConcepts
\cgalConcept

The concept `Proxy` describes the parameterized shape used in the Variational Shape Approximation algorithm.
It is nexsted within the `ErrorMetric` and `ProxyFitting` concepts.

\cgalHasModel `PlaneProxy`
*/

class Proxy {
public:
  /// Triangle mesh facet descriptor.
  typedef unspecified_type facet_descriptor;

  /// @name Data members
  /// @{

  /// Data member to store the proxy seed.
  facet_descriptor seed;

  /// }
};
