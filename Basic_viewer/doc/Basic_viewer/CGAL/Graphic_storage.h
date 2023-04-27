
namespace CGAL {

/*!
 \ingroup PkgBasicViewerClasses

The class `Graphic_storage` store points, segments and triangles. Elements can be added, possibly with associated colors. Non triangular faces can be directly added, and triangulated internally.

\tparam BufferType the type used for point coordinates: float by default.

*/   
template <typename BufferType=float>
class Graphic_storage
{
  template <typename KPoint> void add_point(const KPoint &p);

  template <typename KPoint>
  void add_point(const KPoint &p, const CGAL::IO::Color &acolor);

  template <typename KPoint>
  void add_segment(const KPoint &p1, const KPoint &p2);

  template <typename KPoint>
  void add_segment(const KPoint &p1, const KPoint &p2,
                   const CGAL::IO::Color &acolor);

  template <typename KPoint, typename KVector>
  void add_ray(const KPoint &p, const KVector &v);

  template <typename KPoint, typename KVector>
  void add_ray(const KPoint &p, const KVector &v,
               const CGAL::IO::Color &acolor);

  template <typename KPoint, typename KVector>
  void add_line(const KPoint &p, const KVector &v);

  template <typename KPoint, typename KVector>
  void add_line(const KPoint &p, const KVector &v,
                const CGAL::IO::Color &acolor);

  template <typename KPoint> bool add_point_in_face(const KPoint &kp);

  template <typename KPoint, typename KVector>
  bool add_point_in_face(const KPoint &kp, const KVector &p_normal);

  bool is_a_face_started() const;

  void face_begin();

  void face_begin(const CGAL::IO::Color &acolor);

  void face_end();

  bool is_empty() const;

  void clear();

  template <typename KPoint>
  void add_text(const KPoint &kp, const char *txt);

  template <typename KPoint>
  void add_text(const KPoint &kp, const std::string &txt);

  void m_texts_clear();
};

} // namespace CGAL
