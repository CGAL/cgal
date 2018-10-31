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

  VSA_wrapper()
    : m_metric(L21),
    m_center_pmap(m_face_centers),
    m_area_pmap(m_face_areas),
    m_initialized(false),
    m_pl21_metric(NULL),
    m_l21_approx(NULL),
    m_pl2_metric(NULL),
    m_l2_approx(NULL),
    m_pcompact_metric(NULL),
    m_compact_approx(NULL) {}

  ~VSA_wrapper() {
    if (m_l21_approx)
      delete m_l21_approx;
    if (m_pl21_metric)
      delete m_pl21_metric;
    if (m_l2_approx)
      delete m_l2_approx;
    if (m_pl2_metric)
      delete m_pl2_metric;
    if (m_compact_approx)
      delete m_compact_approx;
    if (m_pcompact_metric)
      delete m_pcompact_metric;
  }

  std::vector<QColor> &proxy_colors() { return m_proxy_colors; }
  const std::vector<QColor> &proxy_colors() const  { return m_proxy_colors; }

  bool initialized() { return m_initialized; }

  void set_mesh(const SMesh &mesh) {
    Vertex_point_map vpm = get(boost::vertex_point, const_cast<SMesh &>(mesh));

    m_face_centers.clear();
    m_face_areas.clear();
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      const halfedge_descriptor he = halfedge(f, mesh);
      const Point_3 &p0 = vpm[source(he, mesh)];
      const Point_3 &p1 = vpm[target(he, mesh)];
      const Point_3 &p2 = vpm[target(next(he, mesh), mesh)];

      m_face_centers.insert(std::pair<face_descriptor, Point_3>(
        f, CGAL::centroid(p0, p1, p2)));

      FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
      m_face_areas.insert(std::pair<face_descriptor, FT>(f, area));
    }

    if (m_l21_approx)
      delete m_l21_approx;
    if (m_pl21_metric)
      delete m_pl21_metric;
    if (m_l2_approx)
      delete m_l2_approx;
    if (m_pl2_metric)
      delete m_pl2_metric;
    if (m_compact_approx)
      delete m_compact_approx;
    if (m_pcompact_metric)
      delete m_pcompact_metric;

    m_pl21_metric = new L21_metric(mesh, vpm);
    m_l21_approx = new L21_approx(mesh, vpm, *m_pl21_metric);

    m_pl2_metric = new L2_metric(mesh, vpm);
    m_l2_approx = new L2_approx(mesh, vpm, *m_pl2_metric);

    m_pcompact_metric = new Compact_metric(m_center_pmap, m_area_pmap);
    m_compact_approx = new Compact_approx(mesh, vpm, *m_pcompact_metric);

    m_proxy_colors.clear();
    m_initialized = false;
  }

  void set_metric(const Metric &m) {
    m_metric = m;
    m_proxy_colors.clear();
    m_initialized = false;
  }

  std::size_t initialize_seeds(const VSA::Seeding_method method,
    const boost::optional<std::size_t> max_nb_of_proxies,
    const boost::optional<FT> min_error_drop,
    const std::size_t nb_relaxations) {
    std::size_t nb_initialized = 0;
    switch (m_metric) {
      case L21:
        nb_initialized = m_l21_approx->initialize_seeds(
          CGAL::parameters::seeding_method(method)
            .max_number_of_proxies(max_nb_of_proxies)
            .min_error_drop(min_error_drop)
            .number_of_relaxations(nb_relaxations));
        break;
      case L2:
        nb_initialized = m_l2_approx->initialize_seeds(
          CGAL::parameters::seeding_method(method)
            .max_number_of_proxies(max_nb_of_proxies)
            .min_error_drop(min_error_drop)
            .number_of_relaxations(nb_relaxations));
        break;
      case Compact:
        nb_initialized = m_compact_approx->initialize_seeds(
          CGAL::parameters::seeding_method(method)
            .max_number_of_proxies(max_nb_of_proxies)
            .min_error_drop(min_error_drop)
            .number_of_relaxations(nb_relaxations));
        break;
    }

    // generate proxy colors
    m_proxy_colors.clear();
    for (std::size_t i = 0; i < number_of_proxies(); ++i) {
      const std::size_t c = rand_0_255();
      m_proxy_colors.push_back(QColor::fromRgb(
        Color_cheat_sheet::r(c), Color_cheat_sheet::g(c), Color_cheat_sheet::b(c)));
    }

    m_initialized = true;
    return nb_initialized;
  }

  void run(const std::size_t nb_iterations) {
    FT err(0.0);
    switch (m_metric) {
      case L21:
        m_l21_approx->run(nb_iterations);
        break;
      case L2:
        m_l2_approx->run(nb_iterations);
        break;
      case Compact:
        m_compact_approx->run(nb_iterations);
        break;
    }
  }

  std::size_t add_one_proxy() {
    std::size_t nb_added = 0;
    switch (m_metric) {
      case L21:
        nb_added = m_l21_approx->add_to_furthest_proxies(1, 0);
        break;
      case L2:
        nb_added = m_l2_approx->add_to_furthest_proxies(1, 0);
        break;
      case Compact:
        nb_added = m_compact_approx->add_to_furthest_proxies(1, 0);
        break;
    }
    if (nb_added == 1) {
      const std::size_t c = rand_0_255();
      m_proxy_colors.push_back(QColor::fromRgb(
        Color_cheat_sheet::r(c), Color_cheat_sheet::g(c), Color_cheat_sheet::b(c)));
    }

    return nb_added;
  }

  std::size_t teleport_one_proxy() {
    switch (m_metric) {
      case L21:
        return m_l21_approx->teleport_proxies(1, 0, true);
      case L2:
        return m_l2_approx->teleport_proxies(1, 0, true);
      case Compact:
        return m_compact_approx->teleport_proxies(1, 0, true);
    }
    return 0;
  }

  bool split(const std::size_t px_idx, const std::size_t n, const std::size_t nb_relaxations) {
    bool splitted = false;
    switch (m_metric) {
      case L21:
        splitted = m_l21_approx->split(px_idx, n, nb_relaxations);
        break;
      case L2:
        splitted = m_l2_approx->split(px_idx, n, nb_relaxations);
        break;
      case Compact:
        splitted = m_compact_approx->split(px_idx, n, nb_relaxations);
        break;
    }
    if (splitted) {
      for (std::size_t i = m_proxy_colors.size(); i < number_of_proxies(); ++i) {
        const std::size_t c = rand_0_255();
        m_proxy_colors.push_back(QColor::fromRgb(
          Color_cheat_sheet::r(c), Color_cheat_sheet::g(c), Color_cheat_sheet::b(c)));
      }
    }

    return splitted;
  }

  bool extract_mesh(const FT subdivision_ratio,
    const bool relative_to_chord,
    const bool with_dihedral_angle,
    const bool optimize_anchor_location,
    const bool pca_plane) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->extract_mesh(
          CGAL::parameters::subdivision_ratio(subdivision_ratio).
          relative_to_chord(relative_to_chord).
          with_dihedral_angle(with_dihedral_angle).
          optimize_anchor_location(optimize_anchor_location).
          pca_plane(pca_plane));
      case L2:
        return m_l2_approx->extract_mesh(
          CGAL::parameters::subdivision_ratio(subdivision_ratio).
          relative_to_chord(relative_to_chord).
          with_dihedral_angle(with_dihedral_angle).
          optimize_anchor_location(optimize_anchor_location).
          pca_plane(pca_plane));
      case Compact:
        return m_compact_approx->extract_mesh(
          CGAL::parameters::subdivision_ratio(subdivision_ratio).
          relative_to_chord(relative_to_chord).
          with_dihedral_angle(with_dihedral_angle).
          optimize_anchor_location(optimize_anchor_location).
          pca_plane(pca_plane));
    }
    return false;
  }

  std::size_t number_of_proxies() {
    switch (m_metric) {
      case L21:
        return m_l21_approx->number_of_proxies();
      case L2:
        return m_l2_approx->number_of_proxies();
      case Compact:
        return m_compact_approx->number_of_proxies();
    }
    return 0;
  }

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  template <typename OutputIterator>
  void get_l21_proxies(OutputIterator outitr) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->wrapped_proxies(outitr);
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
  void indexed_triangles(OutputIterator outitr) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->indexed_triangles(outitr);
      case L2:
        return m_l2_approx->indexed_triangles(outitr);
      case Compact:
        return m_compact_approx->indexed_triangles(outitr);
    }
  }

  template <typename OutputIterator>
  void anchor_points(OutputIterator outitr) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->anchor_points(outitr);
      case L2:
        return m_l2_approx->anchor_points(outitr);
      case Compact:
        return m_compact_approx->anchor_points(outitr);
    }
  }

  template <typename OutputIterator>
  void anchor_vertices(OutputIterator outitr) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->anchor_vertices(outitr);
      case L2:
        return m_l2_approx->anchor_vertices(outitr);
      case Compact:
        return m_compact_approx->anchor_vertices(outitr);
    }
  }

  template <typename OutputIterator>
  void indexed_boundary_polygons(OutputIterator outitr) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->indexed_boundary_polygons(outitr);
      case L2:
        return m_l2_approx->indexed_boundary_polygons(outitr);
      case Compact:
        return m_compact_approx->indexed_boundary_polygons(outitr);
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
