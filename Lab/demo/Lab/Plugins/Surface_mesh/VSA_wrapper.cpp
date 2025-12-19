#include "VSA_wrapper.h"

VSA_wrapper::VSA_wrapper(const SMesh &mesh) :
  m_metric(L21),
  m_center_pmap(get(Face_center_tag(), const_cast<SMesh &>(mesh))),
  m_area_pmap(get(Face_area_tag(), const_cast<SMesh &>(mesh))),
  m_initialized(false)
{
  Vertex_point_map vpm = get(boost::vertex_point, const_cast<SMesh &>(mesh));

  for(face_descriptor f : faces(mesh)) {
    const halfedge_descriptor he = halfedge(f, mesh);
    const Point_3 &p0 = vpm[source(he, mesh)];
    const Point_3 &p1 = vpm[target(he, mesh)];
    const Point_3 &p2 = vpm[target(next(he, mesh), mesh)];

    put(m_center_pmap, f, CGAL::centroid(p0, p1, p2));
    put(m_area_pmap, f, CGAL::sqrt(CGAL::squared_area(p0, p1, p2)));
  }

  m_pl21_metric = new L21_metric(mesh, vpm);
  m_l21_approx = new L21_approx(mesh, vpm, *m_pl21_metric);

  m_pl2_metric = new L2_metric(mesh, vpm);
  m_l2_approx = new L2_approx(mesh, vpm, *m_pl2_metric);

  m_pcompact_metric = new Compact_metric(m_center_pmap, m_area_pmap);
  m_compact_approx = new Compact_approx(mesh, vpm, *m_pcompact_metric);

  m_proxy_colors.clear();
  m_initialized = false;
}

void VSA_wrapper::set_metric(const Metric &m) {
  m_metric = m;
  m_proxy_colors.clear();
  m_initialized = false;
}

EPICK::FT VSA_wrapper::run(const std::size_t nb_iterations) {
  switch (m_metric) {
    case L21: return m_l21_approx->run(nb_iterations);
    case L2: return m_l2_approx->run(nb_iterations);
    case Compact: return m_compact_approx->run(nb_iterations);
  }
  return 0;
}

std::size_t VSA_wrapper::add_one_proxy() {
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
    m_proxy_colors.push_back(generate_random_color());
  }

  return nb_added;
}

std::size_t VSA_wrapper::teleport_one_proxy() {
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

bool VSA_wrapper::split(const std::size_t px_idx, const std::size_t n, const std::size_t nb_relaxations) {
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
      m_proxy_colors.push_back(generate_random_color());
    }
  }

  return splitted;
}

std::size_t VSA_wrapper::number_of_proxies() {
  switch (m_metric) {
    case L21: return m_l21_approx->number_of_proxies();
    case L2: return m_l2_approx->number_of_proxies();
    case Compact: return m_compact_approx->number_of_proxies();
  }
  return 0;
}
