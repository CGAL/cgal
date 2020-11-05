#ifndef VSA_APPROXIMAITON_WRAPPER_H
#define VSA_APPROXIMAITON_WRAPPER_H

#include <CGAL/Variational_shape_approximation.h>
#include <CGAL/Surface_mesh_approximation/L2_metric_plane_proxy.h>
#include <CGAL/Dynamic_property_map.h>

#include <QColor>

#include "SMesh_type.h"
//#include "Color_cheat_sheet.h"
#include "Color_map.h"

#ifdef surface_mesh_approximation_plugin_EXPORTS
#define VSA_WRAPPER_EXPORT Q_DECL_EXPORT
#else
#define VSA_WRAPPER_EXPORT Q_DECL_IMPORT
#endif

namespace VSA = CGAL::Surface_mesh_approximation;

namespace CGAL {
namespace Three {
class Scene_group_item;
}
}
class Scene_polygon_soup_item;
class Scene_polylines_item;

class VSA_WRAPPER_EXPORT VSA_wrapper {
  typedef EPICK::FT FT;
  typedef EPICK::Point_3 Point_3;
  typedef EPICK::Vector_3 Vector_3;

  typedef boost::property_map<SMesh, boost::vertex_point_t>::type Vertex_point_map;
  typedef CGAL::dynamic_face_property_t<Point_3> Face_center_tag;
  typedef CGAL::dynamic_face_property_t<FT> Face_area_tag;
  typedef boost::property_map<SMesh, Face_center_tag>::type Face_center_map;
  typedef boost::property_map<SMesh, Face_area_tag>::type Face_area_map;

  typedef CGAL::Variational_shape_approximation<SMesh, Vertex_point_map,
            CGAL::Default, EPICK, CGAL::Parallel_if_available_tag> L21_approx;
  typedef L21_approx::Error_metric L21_metric;

  typedef VSA::L2_metric_plane_proxy<SMesh> L2_metric;
  typedef CGAL::Variational_shape_approximation<SMesh, Vertex_point_map,
            L2_metric, EPICK, CGAL::Parallel_if_available_tag> L2_approx;

  // user defined point-wise compact metric
  struct Compact_metric_point_proxy {
    typedef Point_3 Proxy;

    Compact_metric_point_proxy(
      const Face_center_map &center_pmap_,
      const Face_area_map &area_pmap_)
      : center_pmap(center_pmap_), area_pmap(area_pmap_) {}

    FT compute_error(const face_descriptor f, const SMesh &, const Proxy &px) const {
      return CGAL::sqrt(CGAL::squared_distance(get(center_pmap, f), px));
    }

    template <typename FaceRange>
    Proxy fit_proxy(const FaceRange &faces, const SMesh &) const {
      CGAL_assertion(!faces.empty());

      // fitting center
      Vector_3 center = CGAL::NULL_VECTOR;
      FT area(0.0);
      for(const face_descriptor f : faces) {
        center = center + (get(center_pmap, f) - CGAL::ORIGIN) * get(area_pmap, f);
        area += get(area_pmap, f);
      }
      center = center / area;
      return CGAL::ORIGIN + center;
    }

    const Face_center_map center_pmap;
    const Face_area_map area_pmap;
  };
  typedef Compact_metric_point_proxy Compact_metric;

  typedef CGAL::Variational_shape_approximation<SMesh, Vertex_point_map,
            Compact_metric, EPICK, CGAL::Parallel_if_available_tag> Compact_approx;

public:
  enum Metric { L21, L2, Compact };

  typedef std::array<std::size_t, 3> Indexed_triangle;

  // visual items
  struct Visual_items {
    Visual_items() :
      group(NULL),
      seeds(NULL),
      has_meshing_items(false),
      triangles(NULL),
      polygons(NULL),
      anchors(NULL),
      planes(NULL) {}

    CGAL::Three::Scene_group_item *group;
    Scene_polylines_item *seeds;
    bool has_meshing_items;
    Scene_polygon_soup_item *triangles;
    Scene_polylines_item *polygons;
    Scene_polylines_item *anchors;
    Scene_polygon_soup_item *planes;
  };

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

  bool &initialized() { return m_initialized; }
  const bool &initialized() const { return m_initialized; }

  void set_metric(const Metric &m);

  template <typename NamedParameters>
  std::size_t initialize_seeds(const NamedParameters &np) {
    std::size_t nb_initialized = 0;
    switch (m_metric) {
      case L21: nb_initialized = m_l21_approx->initialize_seeds(np); break;
      case L2: nb_initialized = m_l2_approx->initialize_seeds(np); break;
      case Compact: nb_initialized = m_compact_approx->initialize_seeds(np); break;
    }

    // generate proxy colors
    m_proxy_colors.clear();
    m_proxy_colors.reserve(number_of_proxies());

    for (std::size_t i = 0; i < number_of_proxies(); ++i) {
      m_proxy_colors.push_back(generate_random_color());
    }

    m_initialized = true;
    return nb_initialized;
  }

  FT run(const std::size_t nb_iterations);

  std::size_t add_one_proxy();

  std::size_t teleport_one_proxy();

  bool split(const std::size_t px_idx, const std::size_t n, const std::size_t nb_relaxations);

  template <typename NamedParameters>
  bool extract_mesh(const NamedParameters &np) {
    switch (m_metric) {
      case L21: return m_l21_approx->extract_mesh(np);
      case L2: return m_l2_approx->extract_mesh(np);
      case Compact: return m_compact_approx->extract_mesh(np);
    }
    return false;
  }

  std::size_t number_of_proxies();

  template <typename NamedParameters>
  void output(const NamedParameters &np) {
    switch (m_metric) {
      case L21: return m_l21_approx->output(np);
      case L2: return m_l2_approx->output(np);
      case Compact: return m_compact_approx->output(np);
    }
  }

  template <typename OutputIterator>
  void anchor_vertices(OutputIterator out_itr) {
    switch (m_metric) {
      case L21: return m_l21_approx->anchor_vertices(out_itr);
      case L2: return m_l2_approx->anchor_vertices(out_itr);
      case Compact: return m_compact_approx->anchor_vertices(out_itr);
    }
  }

  template <typename OutputIterator>
  void indexed_boundary_polygons(OutputIterator out_itr) {
    switch (m_metric) {
      case L21: return m_l21_approx->indexed_boundary_polygons(out_itr);
      case L2: return m_l2_approx->indexed_boundary_polygons(out_itr);
      case Compact: return m_compact_approx->indexed_boundary_polygons(out_itr);
    }
  }

  template <typename OutputIterator>
  void proxy_seeds(OutputIterator out_itr) {
    switch (m_metric) {
      case L21: return proxy_seeds(*m_l21_approx, out_itr);
      case L2: return proxy_seeds(*m_l2_approx, out_itr);
      case Compact: return proxy_seeds(*m_compact_approx, out_itr);
    }
  }

  template <typename Approx, typename OutputIterator>
  void proxy_seeds(Approx &approx, OutputIterator out_itr) {
    typedef typename Approx::Proxy_wrapper Proxy_wrapper;
    std::vector<Proxy_wrapper> pxwrapper;
    approx.wrapped_proxies(std::back_inserter(pxwrapper));
    for(const Proxy_wrapper& pxw : pxwrapper)
      *out_itr++ = pxw.seed;
  }

  Visual_items &visual_items() { return m_visual_items; }
  const Visual_items &visual_items() const { return m_visual_items; }

private:
  Metric m_metric; // current metric

  // face property maps
  Face_center_map m_center_pmap;
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

  // visualization group item
  Visual_items m_visual_items;
};

#endif // VSA_APPROXIMAITON_WRAPPER_H
