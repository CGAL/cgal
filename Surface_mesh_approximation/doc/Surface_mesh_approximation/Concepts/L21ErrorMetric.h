
/*!
\ingroup PkgTSMAConcepts
\cgalConcept

This is a functor to compute the L21 fitting error of a facet to a proxy.
Is this concept needed?

\cgalRefines `ErrorMetric`
\cgalHasModel `L21Metric`

*/
class L21ErrorMetric {
public:
  /*!
   * A number type model of `Field` and `RealEmbeddable`
   */
  typedef unspecified_type FT;
  /// Triangle mesh facet descriptor.
  typedef unspecified_type facet_descriptor;
  /// Parametrized shape proxy.
  typedef unspecified_type Proxy;
  /// 3D vector type
  typedef unspecified_type Vector_3;

  /// @name Functors
  /// @{
  /// Functor for computing scaled vector. Is this required?
  typedef unspecified_type Construct_scaled_vector_3;

  /// }

  /// @name Functions
  /// @{
  /// data member to describe the proxy seed.
  FT operator()(const facet_descriptor &f, const Proxy &px);
  /// }
};

