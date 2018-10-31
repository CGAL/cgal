#ifndef VSA_APPROXIMAITON_WRAPPER_H
#define VSA_APPROXIMAITON_WRAPPER_H

#include <CGAL/Variational_shape_approximation.h>
#include <CGAL/Surface_mesh_approximation/L2_metric_plane_proxy.h>
#include <CGAL/property_map.h>

#include <QColor>

#include "SMesh_type.h"
#include "Color_cheat_sheet.h"

namespace VSA = CGAL::Surface_mesh_approximation;

class VSA_wrapper {
  typedef typename EPICK::FT FT;
  typedef typename EPICK::Point_3 Point_3;
  typedef typename EPICK::Vector_3 Vector_3;

  typedef typename boost::property_map<SMesh, boost::vertex_point_t>::type Vertex_point_map;
  typedef boost::associative_property_map<std::map<face_descriptor, FT> > Face_area_map;
  typedef boost::associative_property_map<std::map<face_descriptor, Point_3> > Face_center_map;

#ifdef CGAL_LINKED_WITH_TBB
  typedef CGAL::Variational_shape_approximation<SMesh, Vertex_point_map,
    CGAL::Default, EPICK, CGAL::Parallel_tag> L21_approx;
#else
  typedef CGAL::Variational_shape_approximation<SMesh, Vertex_point_map,
    CGAL::Default, EPICK> L21_approx;
#endif
  typedef typename L21_approx::Error_metric L21_metric;

  typedef VSA::L2_metric_plane_proxy<SMesh> L2_metric;
#ifdef CGAL_LINKED_WITH_TBB
  typedef CGAL::Variational_shape_approximation<SMesh, Vertex_point_map,
    L2_metric, EPICK, CGAL::Parallel_tag> L2_approx;
#else
  typedef CGAL::Variational_shape_approximation<SMesh, Vertex_point_map,
    L2_metric, EPICK> L2_approx;
#endif

  // user defined point-wise compact metric
  struct Compact_metric_point_proxy {
    typedef Point_3 Proxy;

    Compact_metric_point_proxy(const Face_center_map &_center_pmap,
      const Face_area_map &_area_pmap)
      : center_pmap(_center_pmap), area_pmap(_area_pmap) {}

    FT compute_error(const face_descriptor f, const SMesh &, const Proxy &px) const {
      return FT(std::sqrt(CGAL::to_double(
        CGAL::squared_distance(center_pmap[f], px))));
    }

    template <typename FaceRange>
    Proxy fit_proxy(const FaceRange &faces, const SMesh &) const {
      CGAL_assertion(!faces.empty());

      // fitting center
      Vector_3 center = CGAL::NULL_VECTOR;
      FT area(0.0);
      BOOST_FOREACH(const face_descriptor f, faces) {
        center = center + (center_pmap[f] - CGAL::ORIGIN) * area_pmap[f];
        area += area_pmap[f];
      }
      center = center / area;
      return CGAL::ORIGIN + center;
    }

    const Face_center_map center_pmap;
    const Face_area_map area_pmap;
  };
  typedef Compact_metric_point_proxy Compact_metric;

#ifdef CGAL_LINKED_WITH_TBB
  typedef CGAL::Variational_shape_approximation<SMesh, Vertex_point_map,
    Compact_metric, EPICK, CGAL::Parallel_tag> Compact_approx;
#else
  typedef CGAL::Variational_shape_approximation<SMesh, Vertex_point_map,
    Compact_metric, EPICK> Compact_approx;
#endif

  std::size_t rand_0_255() {
    return static_cast<std::size_t>(std::rand() % 255);
  }

public:
  enum Metric { L21, L2, Compact };

  typedef CGAL::cpp11::array<std::size_t, 3> Indexed_triangle;

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  typedef typename L21_approx::Proxy_wrapper L21_proxy_wrapper;
#endif

  VSA_wrapper(const SMesh &mesh);

  ~VSA_wrapper() {
    delete m_l21_approx;
    delete m_pl21_metric;
    delete m_l2_approx;
    delete m_pl2_metric;
    delete m_compact_approx;
    delete m_pcompact_metric;
  }

  std::vector<QColor> &proxy_colors() { return m_proxy_colors; }
  const std::vector<QColor> &proxy_colors() const  { return m_proxy_colors; }

  bool initialized() { return m_initialized; }

  void set_metric(const Metric &m);

  std::size_t initialize_seeds(const VSA::Seeding_method method,
    const boost::optional<std::size_t> max_nb_of_proxies,
    const boost::optional<FT> min_error_drop,
    const std::size_t nb_relaxations);

  void run(const std::size_t nb_iterations);

  std::size_t add_one_proxy();

  std::size_t teleport_one_proxy();

  bool split(const std::size_t px_idx, const std::size_t n, const std::size_t nb_relaxations);

  bool extract_mesh(const FT subdivision_ratio,
    const bool relative_to_chord,
    const bool with_dihedral_angle,
    const bool optimize_anchor_location,
    const bool pca_plane);

  std::size_t number_of_proxies();

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  template <typename OutputIterator>
  void get_l21_proxies(OutputIterator out_itr) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->wrapped_proxies(out_itr);
      default:
        return;
    }
  }
#endif

  template <typename FaceProxyMap>
  void proxy_map(FaceProxyMap &fpmap) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->proxy_map(fpmap);
      case L2:
        return m_l2_approx->proxy_map(fpmap);
      case Compact:
        return m_compact_approx->proxy_map(fpmap);
    }
  }

  template <typename OutputIterator>
  void indexed_triangles(OutputIterator out_itr) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->indexed_triangles(out_itr);
      case L2:
        return m_l2_approx->indexed_triangles(out_itr);
      case Compact:
        return m_compact_approx->indexed_triangles(out_itr);
    }
  }

  template <typename OutputIterator>
  void anchor_points(OutputIterator out_itr) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->anchor_points(out_itr);
      case L2:
        return m_l2_approx->anchor_points(out_itr);
      case Compact:
        return m_compact_approx->anchor_points(out_itr);
    }
  }

  template <typename OutputIterator>
  void anchor_vertices(OutputIterator out_itr) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->anchor_vertices(out_itr);
      case L2:
        return m_l2_approx->anchor_vertices(out_itr);
      case Compact:
        return m_compact_approx->anchor_vertices(out_itr);
    }
  }

  template <typename OutputIterator>
  void indexed_boundary_polygons(OutputIterator out_itr) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->indexed_boundary_polygons(out_itr);
      case L2:
        return m_l2_approx->indexed_boundary_polygons(out_itr);
      case Compact:
        return m_compact_approx->indexed_boundary_polygons(out_itr);
    }
  }

private:
  Metric m_metric; // current metric

  // face property maps
  std::map<face_descriptor, Point_3> m_face_centers;
  Face_center_map m_center_pmap;
  std::map<face_descriptor, FT> m_face_areas;
  Face_area_map m_area_pmap;

  // patch color
  std::vector<QColor> m_proxy_colors;

  bool m_initialized;

  L21_metric *m_pl21_metric;
  L21_approx *m_l21_approx;

  L2_metric *m_pl2_metric;
  L2_approx *m_l2_approx;

  Compact_metric *m_pcompact_metric;
  Compact_approx *m_compact_approx;
};

#endif // VSA_APPROXIMAITON_WRAPPER_H
