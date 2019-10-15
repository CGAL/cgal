/*!
\ingroup PkgTSMAConcepts
\cgalConcept

The concept `ErrorMetricProxy` defines the notion of proxy,
computes the fitting error from a face to a proxy,
and fits a proxy from a range of faces.

\cgalHasModel `CGAL::Surface_mesh_approximation::L21_metric_plane_proxy`
\cgalHasModel `CGAL::Surface_mesh_approximation::L2_metric_plane_proxy`
*/

class ErrorMetricProxy {
public:
  /// Number type model of `Field` and `RealEmbeddable`
  typedef unspecified_type FT;
  /// Triangle mesh
  typedef unspecified_type Triangle_mesh;
  /// Triangle mesh face descriptor
  typedef unspecified_type face_descriptor;
  /// Shape proxy type
  typedef unspecified_type Proxy;

  /// @name Operations
  /// A model of this concept must provide:
  /// @{

  /// returns fitting error from face f to proxy.
  FT compute_error(const face_descriptor f, const Triangle_mesh &tm, const Proxy &proxy) const;

  /// returns a fitted proxy to a range of faces.
  /// @tparam FaceRange a range of
  ///   `boost::graph_traits<TriangleMesh>::%face_descriptor` model of `ConstRange`
  ///   with iterator type being model of `InputIterator`.
  template <typename FaceRange>
  Proxy fit_proxy(const FaceRange &faces, const Triangle_mesh &tm) const;

  /// @}

};
