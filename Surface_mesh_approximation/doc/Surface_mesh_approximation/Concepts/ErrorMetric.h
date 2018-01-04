/*!
\ingroup PkgTSMAConcepts
\cgalConcept

The concept `ErrorMetric` computes the fitting error of a face to a `Proxy`, used in CGAL::VSA::Mesh_approximation.

\cgalHasModel `CGAL::VSA::L21_metric`
\cgalHasModel `CGAL::VSA::L2_metric`
*/

class ErrorMetric {
public:
  /// A number type model of `Field` and `RealEmbeddable`
  typedef unspecified_type FT;
  /// Triangle mesh face descriptor.
  typedef unspecified_type face_descriptor;
  /// Parametrized shape proxy.
  typedef unspecified_type Proxy;

  /// @name Operations
  /// A model of this concept must provide:
  /// @{

  /// returns the fitting error of a face f to the proxy px.
  FT operator()(const face_descriptor &f, const Proxy &px) const;

  /// }
};
