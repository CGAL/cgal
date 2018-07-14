/*!
\ingroup PkgTSMAConcepts
\cgalConcept

The concept `ErrorMetricProxy` defines the notion of proxy,
computes the fitting error from a face to a proxy,
and fits a proxy from a range of faces.

\cgalHasModel `CGAL::VSA::L21_metric_plane_proxy`
\cgalHasModel `CGAL::VSA::L2_metric_plane_proxy`
*/

class ErrorMetricProxy {
public:
  /// Number type model of `Field` and `RealEmbeddable`
  typedef unspecified_type FT;
  /// Triangle mesh
  typedef unspecified_type TriangleMesh;
  /// Triangle mesh face descriptor
  typedef unspecified_type face_descriptor;
  /// Parameterized shape proxy
  typedef unspecified_type Proxy;

  /// @name Operations
  /// A model of this concept must provide:
  /// @{

  /// Computes and returns fitting error from face f to proxy.
  FT compute_error(const TriangleMesh &tm, const face_descriptor &f, const Proxy &proxy) const;

  /// Computes and returns fitted proxy from a range of faces.
  /// @tparam FaceRange a range of
  ///   `boost::graph_traits<TriangleMesh>::%face_descriptor` model of `ConstRange`
  ///   with iterator type being model of `InputIterator`.
  template <typename FaceRange>
  Proxy fit_proxy(const FaceRange &faces, const TriangleMesh &tm) const;

  /// }

};
