
namespace CGAL {

/*!
 \ingroup PkgBasicViewerClasses

The class `Graphic_storage` store points, segments, triangles, rays and lines. Elements can be added, possibly with associated colors. Non triangular faces can be directly added, and are triangulated internally.

\tparam BufferType the type used for point coordinates: `float` by default.

*/
template <typename BufferType=float>
class Graphic_storage {
public:
  /// add the given point in the storage.
  template <typename KPoint>
  void add_point(const KPoint &p);

  /// add the given colored point in the storage.
  template <typename KPoint>
  void add_point(const KPoint &p, const CGAL::IO::Color &acolor);

  /// add the given segment in the storage.
  template <typename KPoint>
  void add_segment(const KPoint &p1, const KPoint &p2);

  /// add the given colored segment in the storage.
  template <typename KPoint>
  void add_segment(const KPoint &p1, const KPoint &p2,
                   const CGAL::IO::Color &acolor);

  /// add the given ray in the storage: an half line starting from `p` and having `v` as direction.
  template <typename KPoint, typename KVector>
  void add_ray(const KPoint &p, const KVector &v);

  /// add the given colored ray in the storage: an half line starting from `p` and having `v` as direction.
  template <typename KPoint, typename KVector>
  void add_ray(const KPoint &p, const KVector &v,
               const CGAL::IO::Color &acolor);

  /// add the given line in the storage, defined by `p` and `v` as direction.
  template <typename KPoint, typename KVector>
  void add_line(const KPoint &p, const KVector &v);

  /// add the given colored line in the storage, defined by `p` and `v` as direction.
  template <typename KPoint, typename KVector>
  void add_line(const KPoint &p, const KVector &v,
                const CGAL::IO::Color &acolor);

  /// start a new face.
  void face_begin();

  /// start a new colored face.
  void face_begin(const CGAL::IO::Color &acolor);

  /// @return `true` iff a face is started.
  bool is_a_face_started() const;

  /// add the given point in the current face.
  /// @pre `is_a_face_started()`
  template <typename KPoint> bool add_point_in_face(const KPoint &kp);

  /// add the given point in the current face, having the vertex normal.
  /// @pre `is_a_face_started()`
  template <typename KPoint, typename KVector>
  bool add_point_in_face(const KPoint &kp, const KVector &p_normal);

  /// end the current face.
  /// @pre `is_a_face_started()`
  void face_end();

  /// add the given text at the given position in the storage.
  template <typename KPoint>
  void add_text(const KPoint &kp, const char *txt);

  /// add the given text at the given position in the storage.
  template <typename KPoint>
  void add_text(const KPoint &kp, const std::string &txt);

  /// @return `true` iff the storage has no element.
  bool is_empty() const;

  /// clear the storage, i.e. remove all points, segments, triangles and text.
  void clear();
};

} // namespace CGAL
