#include <fstream>

// Qt
#include <QtOpenGL>

// local
#include "Otr2_kerneled.h"
#include "glviewer.h"
#include <CGAL/Optimal_transportation_reconstruction_2.h>


typedef Optimal_transportation_reconstruction_kerneled_2::Rec_edge_2 PEdge;
typedef Optimal_transportation_reconstruction_kerneled_2 R_s_k_2;

void R_s_k_2::print_stats() const
{
  int nb_solid = 0;
  int nb_ghost = 0;
  for (Finite_edges_iterator ei = m_dt.finite_edges_begin();
      ei != m_dt.finite_edges_end(); ++ei)
  {
    Edge edge = *ei;
    if (m_dt.is_ghost(edge)) nb_ghost++;
    else nb_solid++;
  }

  std::cerr << "STATS" << std::endl;
  std::cerr << "# vertices : " << m_dt.number_of_vertices()-4 << std::endl;
  std::cerr << "# triangles: " << m_dt.number_of_faces() << std::endl;
  std::cerr << "# edges: " << m_dt.tds().number_of_edges() << std::endl;
  std::cerr << "# solid: " << nb_solid << std::endl;
  std::cerr << "# ghost: " << nb_ghost << std::endl;
}

QColor R_s_k_2::get_color(float value) const
{
  float hue = 240.0*(1.0 - value);
  QColor color;
  color.setHsv(hue, 255, 255);
  return color;
}

void R_s_k_2::draw_point(const Point& point)
{
  viewer->glBegin(GL_POINTS);
  viewer->glVertex2f(point.x(), point.y());
  viewer->glEnd();
}

void R_s_k_2::draw_segment(const Point& s, const Point& t)
{
  viewer->glBegin(GL_LINES);
  viewer->glVertex2d(s.x(), s.y());
  viewer->glVertex2d(t.x(), t.y());
  viewer->glEnd();
}

void R_s_k_2::draw_edge(const Edge& edge)
{
  int i = edge.second;
  Face_handle face = edge.first;
  Point a = face->vertex((i+1)%3)->point();
  Point b = face->vertex((i+2)%3)->point();
  draw_segment(a, b);
}

void R_s_k_2::draw_face(Face_handle face)
{
  viewer->glBegin(GL_TRIANGLES);
  for (int i = 0; i < 3; ++i)
  {
    Point p = face->vertex(i)->point();
    viewer->glVertex2f(p.x(), p.y());
  }
  viewer->glEnd();
}

void R_s_k_2::draw_edge_with_arrow(const Point& s, const Point& t)
{
  Vector vec = t - s;
  Vector vec90(-vec.y(),vec.x());

  // draw edge
  draw_segment(s, t);

  // draw an arrow toward merged vertex
  Point a = t - 0.4 * vec;
  Point b = a - 0.2 * vec - 0.1 * vec90;
  Point c = a - 0.2 * vec + 0.1 * vec90;
  viewer->glBegin(GL_TRIANGLES);
  viewer->glVertex2d(a.x(), a.y());
  viewer->glVertex2d(b.x(), b.y());
  viewer->glVertex2d(c.x(), c.y());
  viewer->glEnd();
}

void R_s_k_2::draw_vertices(const float point_size,
    const float red,
    const float green,
    const float blue)
{
  for (Finite_vertices_iterator vi = m_dt.finite_vertices_begin();
      vi != m_dt.finite_vertices_end(); vi++)
  {
    Vertex_handle vertex = vi;
    if (vertex->pinned())
    {
      viewer->glPointSize(point_size);
      viewer->glColor3f(0.0f, 0.0f, 0.0f);
    }
    else
    {
      viewer->glPointSize(3*point_size);
      viewer->glColor3f(red,green,blue);
    }
    draw_point(vertex->point());
  }
}

void R_s_k_2::draw_edges(const float line_width,
    const float red,
    const float green,
    const float blue)
{
  viewer->glLineWidth(line_width);
  for (Finite_edges_iterator ei = m_dt.finite_edges_begin();
      ei != m_dt.finite_edges_end(); ei++)
  {
    Edge edge = *ei;
    Edge twin = m_dt.twin_edge(edge);
    if (m_dt.is_pinned(edge) && m_dt.is_pinned(twin))
      viewer->glColor3f(0.9f,0.9f,0.75f);
    else
      viewer->glColor3f(red,green,blue);
    draw_edge(edge);
  }
}

void R_s_k_2::draw_footpoints(const float line_width,
    const float red,
    const float green,
    const float blue)
{
  draw_mesh_footpoints(m_dt, line_width, red, green, blue);
}

void R_s_k_2::draw_mesh_footpoints(const Triangulation& mesh,
    const float line_width,
    const float red,
    const float green,
    const float blue)
{
  viewer->glLineWidth(line_width);
  for (Finite_edges_iterator ei = mesh.finite_edges_begin(); ei != mesh.finite_edges_end(); ei++)
  {
    Edge edge = *ei;
    draw_edge_footpoints(mesh, edge, red, green, blue);
    draw_edge_footpoints(mesh, mesh.twin_edge(edge), red, green, blue);
  }
}

void R_s_k_2::draw_edge_footpoints(const Triangulation& mesh,
    const Edge& edge,
    const float red,
    const float green,
    const float blue)
{
  const Point& a = mesh.source_vertex(edge)->point();
  const Point& b = mesh.target_vertex(edge)->point();
  const Sample_vector& samples = edge.first->samples(edge.second);

  Sample_vector::const_iterator it;
  for (it = samples.begin(); it != samples.end(); ++it)
  {
    Sample_* sample = *it;
    Point p = sample->point();
    FT m = 0.5*(1.0 - sample->mass());

    Point q;
    if (mesh.get_plan(edge) == 0)
    {
      viewer->glColor3f(0.8f + m, m, m);
      FT Da = CGAL::squared_distance(p, a);
      FT Db = CGAL::squared_distance(p, b);
      if (Da < Db) q = a;
      else         q = b;
    }
    else
    {
      viewer->glColor3f(red + m, green + m, blue + m);
      FT t = sample->coordinate();
      q = CGAL::ORIGIN + (1.0 - t)*(a - CGAL::ORIGIN) + t*(b - CGAL::ORIGIN);
    }
    draw_segment(p, q);
  }
}

void R_s_k_2::draw_pedges(const float line_width)
{
  int nb_edges = 0;
  int nb_pinned = 0;
  int nb_cyclic = 0;
  int nb_discart = 0;
  FT min_value = (std::numeric_limits<FT>::max)();
  FT max_value = -(std::numeric_limits<FT>::max)();
  std::vector<FT>  values;
  std::vector<Edge> edges;
  for (Finite_edges_iterator ei = m_dt.finite_edges_begin(); ei != m_dt.finite_edges_end(); ++ei)
  {
    Edge edge = *ei;
    for (unsigned int i = 0; i < 2; ++i)
    {
      if (m_dt.is_pinned(edge))
      {
        nb_pinned++;
        continue;
      }

      if (m_dt.is_target_cyclic(edge))
      {
        nb_cyclic++;
        continue;
      }

      PEdge pedge;
      bool ok = create_pedge(edge, pedge);
      if (ok)
      {
        edges.push_back(edge);
        values.push_back(pedge.priority());
        min_value = (std::min)(min_value, values.back());
        max_value = (std::max)(max_value, values.back());
      }
      else
      {
        nb_discart++;
        viewer->glColor3f(1.0, 0.0, 1.0);
        draw_edge(edge);
      }
      edge = m_dt.twin_edge(edge);
      nb_edges++;
    }
  }
  if (min_value == max_value) max_value += 1.0;

  std::size_t N = values.size();
  for (unsigned int i = 0; i < N; ++i)
    draw_one_pedge(edges[i], values[i], min_value, max_value, line_width);

  std::cout << "There are: " << N << " pedges"
      << " x " << nb_discart << " discarted"
      << " x " << nb_pinned << " pinned"
      << " x " << nb_cyclic << " cyclic"
      << " = " << nb_edges << " edges"
      << std::endl;
}

void R_s_k_2::draw_one_pedge(const Edge& edge,
    const FT value,
    const FT min_value,
    const FT max_value,
    const float line_width)
{
  if (value == min_value)
  {
    viewer->glLineWidth(2*line_width);
    viewer->glColor3f(0.0f, 0.0f, 0.0f);
  }
  else
  {
    viewer->glLineWidth(line_width);
    FT color = (value - min_value) / (max_value - min_value);
    QColor qcolor = get_color(color);
    viewer->glColor3f(qcolor.redF(), qcolor.greenF(), qcolor.blueF());
  }

  Point s = m_dt.source_vertex(edge)->point();
  Point t = m_dt.target_vertex(edge)->point();
  Point c = CGAL::midpoint(s, t);
  draw_edge_with_arrow(c, t);
}

void R_s_k_2::draw_costs(const float line_width, const bool view_ghost)
{
  FT min_value = (std::numeric_limits<FT>::max)();
  FT max_value = -(std::numeric_limits<FT>::max)();
  for (Finite_edges_iterator ei = m_dt.finite_edges_begin(); ei != m_dt.finite_edges_end(); ++ei)
  {
    Edge edge = *ei;
    if (m_dt.is_ghost(edge)) continue;
    FT value = m_dt.get_cost(edge).finalize(m_alpha);
    min_value = (std::min)(min_value, value);
    max_value = (std::max)(max_value, value);
  }
  if (min_value == max_value) max_value += 1.0;

  viewer->glLineWidth(line_width);
  for (Finite_edges_iterator ei = m_dt.finite_edges_begin(); ei != m_dt.finite_edges_end(); ++ei)
  {
    Edge edge = *ei;
    draw_one_cost(edge, min_value, max_value, view_ghost);
  }
}

void R_s_k_2::draw_one_cost(const Edge& edge,
    const FT min_value,
    const FT max_value,
    const bool view_ghost)
{
  FT mass = m_dt.get_mass(edge);
  if (mass == 0.0)
  {
    if (!view_ghost) return;
    viewer->glColor3f(0.5f, 0.5f, 0.5f);
    draw_edge(edge);
    return;
  }

  if (m_dt.is_ghost(edge))
  {
    if (!view_ghost) return;
    viewer->glColor3f(0.5f, 0.5f, 0.5f);
    draw_edge(edge);
    return;
  }

  FT value = m_dt.get_cost(edge).finalize(m_alpha);
  FT color = (value - min_value) / (max_value - min_value);
  viewer->glColor3d(0.0, 1.0-color, color); // [green, blue]
                                          draw_edge(edge);
}

void R_s_k_2::draw_relevance(const float line_width, const int nb)
{
  MultiIndex mindex;
  FT min_value = (std::numeric_limits<FT>::max)();
  FT max_value = -(std::numeric_limits<FT>::max)();
  unsigned int nb_initial = 0;
  for (Finite_edges_iterator ei = m_dt.finite_edges_begin(); ei != m_dt.finite_edges_end(); ++ei)
  {
    Edge edge = *ei;
    if (m_dt.is_ghost(edge)) continue;
    FT value = m_dt.get_edge_relevance(edge); // >= 0

    nb_initial++;
    min_value = (std::min)(min_value, value);
    max_value = (std::max)(max_value, value);
    mindex.insert(PEdge(edge, value));
  }
  if (min_value == max_value) max_value += 1.0;

  viewer->glLineWidth(line_width);
  int nb_remove = (std::min)(nb, int(mindex.size()));

  viewer->glColor3d(0.5, 0.1, 0.1);
  for (int i = 0; i < nb_remove; ++i)
  {

    PEdge pedge = *(mindex.get<1>()).begin();
    (mindex.get<0>()).erase(pedge);
  }

  viewer->glColor3d(0.0, 0.5, 0.0);
  while (!mindex.empty())
  {
    PEdge pedge = *(mindex.get<1>()).begin();
    (mindex.get<0>()).erase(pedge);
    draw_edge(pedge.edge());
  }
}

void R_s_k_2::draw_bins(const float thickness)
{
  viewer->glLineWidth(thickness);
  for (Finite_edges_iterator ei = m_dt.finite_edges_begin(); ei != m_dt.finite_edges_end(); ++ei)
  {
    Edge edge = *ei;
    if (m_dt.get_plan(edge) == 0)
      draw_bins_plan0(edge);
    else
      draw_bins_plan1(edge);
  }
}

void R_s_k_2::draw_bins_plan0(const Edge& edge)
{
  Edge twin = m_dt.twin_edge(edge);
  const Point& pa = m_dt.source_vertex(edge)->point();
  const Point& pb = m_dt.target_vertex(edge)->point();

  Sample_vector samples;
  m_dt.collect_samples_from_edge(edge, samples);
  m_dt.collect_samples_from_edge(twin, samples);

  viewer->glColor3f(0.0f, 1.0f, 0.0f);
  Sample_vector_const_iterator it;
  for (it = samples.begin(); it != samples.end(); ++it)
  {
    Sample_* sample = *it;
    const Point& ps = sample->point();

    Point q = pa;
    FT Da = CGAL::squared_distance(ps, pa);
    FT Db = CGAL::squared_distance(ps, pb);
    if (Da > Db) q = pb;

    draw_segment(ps, q);
  }
}

void R_s_k_2::draw_bins_plan1(const Edge& edge)
{
  FT M = m_dt.get_mass(edge);
  Vector va = m_dt.source_vertex(edge)->point() - CGAL::ORIGIN;
  Vector vb = m_dt.target_vertex(edge)->point() - CGAL::ORIGIN;

  viewer->glColor3f(1.0f, 0.0f, 0.0f);
  SQueue queue;
  FT start = 0.0;
  m_dt.sort_samples_from_edge(edge, queue);
  while (!queue.empty())
  {
    PSample psample = queue.top();
    queue.pop();

    const FT m = psample.sample()->mass();
    const Point& ps = psample.sample()->point();

    FT bin = m/M;
    FT alpha = start + 0.5*bin;
    Point q = CGAL::ORIGIN + (1.0-alpha)*va + alpha*vb;
    start += bin;

    draw_segment(ps, q);
  }
}

void R_s_k_2::draw_relocation()
{
  for (Finite_vertices_iterator v = m_dt.finite_vertices_begin(); v != m_dt.finite_vertices_end(); ++v)
  {
    Vertex_handle vertex = v;
    if (vertex->pinned()) continue;

    const Point& pv = vertex->point();
    v->relocated() = compute_relocation(vertex);

    Vector move(0.0, 0.0);
    Edge_circulator ecirc = m_dt.incident_edges(vertex);
    Edge_circulator eend  = ecirc;
    CGAL_For_all(ecirc, eend)
    {
      Edge edge = *ecirc;
      if (m_dt.source_vertex(edge) != vertex)
        edge = m_dt.twin_edge(edge);

      Vector grad(0.0, 0.0);
      if (m_dt.get_plan(edge) == 0)
        grad = compute_gradient_for_plan0(edge);
      else
        grad = compute_gradient_for_plan1(edge);

      move = move + grad;
      viewer->glLineWidth(2.0f);
      viewer->glColor3f(1.0f, 1.0f, 0.0f);
      draw_edge_with_arrow(pv, pv-grad);
    }

    viewer->glLineWidth(1.0f);
    viewer->glColor3f(1.0f, 0.0f, 0.0f);
    draw_edge_with_arrow(pv, pv-move);
  }

  viewer->glBegin(GL_LINES);
  viewer->glLineWidth(3.0f);
  viewer->glColor3f(0.1f, 1.0f, 1.0f);
  for (Finite_edges_iterator ei = m_dt.finite_edges_begin(); ei != m_dt.finite_edges_end(); ++ei)
  {
    Edge edge = *ei;
    Edge twin = m_dt.twin_edge(edge);
    if (m_dt.is_pinned(edge) || m_dt.is_pinned(twin)) continue;

    const Point& pa = m_dt.source_vertex(edge)->relocated();
    const Point& pb = m_dt.target_vertex(edge)->relocated();
    viewer->glVertex2d(pa.x(), pa.y());
    viewer->glVertex2d(pb.x(), pb.y());
  }
  viewer->glEnd();
}

bool R_s_k_2::locate_edge(const Point& query, Edge& edge)
{
  if (m_dt.number_of_faces() == 0) return false;

  Face_handle face = m_dt.locate(query);
  if (face == Face_handle()) return false;
  if (m_dt.is_infinite(face)) return false;

  edge = m_dt.find_nearest_edge(query, face);
  if (edge.first == Face_handle()) return false;

  if (m_dt.is_pinned(edge) || m_dt.is_target_cyclic(edge))
    return false;

  return true;
}

void R_s_k_2::draw_one_ring(const float point_size, const float line_width, const Point& query)
{
  Edge edge;
  bool ok = locate_edge(query, edge);
  if (!ok) return;

  Triangulation copy;
  Edge copy_edge = copy_star(edge, copy);
  draw_mesh_one_ring(point_size, line_width, copy, copy_edge);
}

void R_s_k_2::draw_mesh_one_ring(const float point_size,
    const float line_width,
    const Triangulation& mesh,
    const Edge& edge)
{
  Vertex_handle s = mesh.source_vertex(edge);
  Vertex_handle t = mesh.target_vertex(edge);

  draw_bg_faces(mesh, 1.0f, 0.7f, 0.7f, 0.5f);
  draw_vertex_faces(s, mesh, 0.7f, 0.7f, 1.0f, 1.0f);

  viewer->glLineWidth(line_width);
  draw_bg_edges(mesh, 0.7f, 0.3f, 0.7f, 1.f, 0.f, 0.f);
  draw_vertex_edges(s, mesh, 0.f, 0.8f, 0.f, 0.f, 0.2f, 0.2f);

  viewer->glLineWidth(2.0f*line_width);
  viewer->glColor3f(0., 0., 1.);
  draw_edge_with_arrow(s->point(), t->point());

  viewer->glPointSize(0.5*point_size);
  draw_bg_vertices(mesh, 0.f, 0.f, 0.f);
  viewer->glPointSize(point_size);
  viewer->glColor3f(0., 1., 0.);
  draw_point(s->point());
  viewer->glColor3f(1., 1., 0.);
  draw_point(t->point());
}

void R_s_k_2::draw_blocking_edges(const float point_size, const float line_width, const Point& query)
{
  Edge edge;
  bool ok = locate_edge(query, edge);
  if (!ok) return;

  Triangulation copy;
  Edge copy_edge = copy_star(edge, copy);
  draw_mesh_blocking_edges(point_size, line_width, copy, copy_edge);
}

void R_s_k_2::draw_mesh_blocking_edges(const float point_size,
    const float line_width,
    const Triangulation& mesh,
    const Edge& edge)
{
  Vertex_handle s = mesh.source_vertex(edge);
  Vertex_handle t = mesh.target_vertex(edge);

  draw_mesh_one_ring(point_size, line_width, mesh, edge);

  viewer->glColor3f(0.0f, 0.0f, 0.0f);
  viewer->glLineWidth(2.0f*line_width);
  Face_circulator fcirc = mesh.incident_faces(s);
  Face_circulator fend = fcirc;
  CGAL_For_all(fcirc, fend)
  {
    Face_handle f = fcirc;
    Edge ab(f, f->index(s));
    Vertex_handle a = mesh.source_vertex(ab);
    Vertex_handle b = mesh.target_vertex(ab);
    if (!mesh.is_triangle_ccw(a, b, t))
    {
      draw_segment(a->point(), b->point());
    }
  }
}

void R_s_k_2::draw_collapsible_edge(const float point_size,
    const float line_width,
    const Point& query)
{
  Edge edge;
  bool ok = locate_edge(query, edge);
  if (!ok) return;

  Triangulation copy;
  Edge copy_edge = copy_star(edge, copy);
  Vertex_handle copy_src = copy.source_vertex(copy_edge);

  Edge_vector copy_hull;
  copy.get_edges_from_star_minus_link(copy_src, copy_hull, true);
  ok = copy.make_collapsible(copy_edge, copy_hull.begin(), copy_hull.end(), m_verbose);

  if (ok)
    draw_mesh_one_ring(point_size, line_width, copy, copy_edge);
  else
    draw_mesh_blocking_edges(point_size, line_width, copy, copy_edge);
}

void R_s_k_2::draw_cost_stencil(const float point_size,
    const float line_width,
    const Point& query)
{
  Edge edge;
  bool ok = locate_edge(query, edge);
  if (!ok) return;

  Triangulation copy;
  Edge copy_edge = copy_star(edge, copy);
  Vertex_handle copy_src = copy.source_vertex(copy_edge);

  Edge_vector copy_hull;
  copy.get_edges_from_star_minus_link(copy_src, copy_hull, true);
  ok = copy.make_collapsible(copy_edge, copy_hull.begin(), copy_hull.end(), m_verbose);
  if (!ok) return;
  copy.collapse(copy_edge, m_verbose);

  draw_bg_faces(copy, 0.7f, 0.7f, 1.0f, 1.0f);

  viewer->glLineWidth(line_width);
  viewer->glPointSize(point_size);

  Edge_vector stencil;
  collect_cost_stencil(copy, copy_hull.begin(), copy_hull.end(), stencil);
  for (Edge_vector::const_iterator it = stencil.begin(); it != stencil.end(); ++it)
  {
    Edge e = *it;
    viewer->glColor3f(0.7f, 0.4f, 0.0f);
    draw_edge(e);
    viewer->glColor3f(0.0f, 0.0f, 0.0f);
    draw_point(copy.source_vertex(e)->point());
    draw_point(copy.target_vertex(e)->point());
  }
}

void R_s_k_2::draw_remove_queue_stencil(const float point_size,
    const float line_width,
    const Point& query)
{
  Edge edge;
  bool ok = locate_edge(query, edge);
  if (!ok) return;

  Edge_vector hull;
  Edge_vector stencil;
  Edge_vector::const_iterator it;
  Vertex_handle src = m_dt.source_vertex(edge);
  m_dt.get_edges_from_star_minus_link(src, hull, true);
  collect_pqueue_stencil(m_dt, hull.begin(), hull.end(), stencil);

  draw_vertex_faces(src, m_dt, 0.7f, 0.7f, 1.0f, 1.0f);

  viewer->glLineWidth(0.5*line_width);
  for (it = stencil.begin(); it != stencil.end(); ++it)
  {
    Edge ab = *it;
    viewer->glColor3f(1.0f, 0.6f, 1.0f);
    draw_edge(ab);
  }

  viewer->glLineWidth(line_width);
  viewer->glPointSize(point_size);
  for (it = stencil.begin(); it != stencil.end(); ++it)
  {
    Edge ab = *it;
    Vertex_handle va = ab.first->vertex( (ab.second+1)%3 );
    Vertex_handle vb = ab.first->vertex( (ab.second+2)%3 );
    const Point& pa = va->point();
    const Point& pb = vb->point();
    Point pc = CGAL::midpoint(pa, pb);
    viewer->glColor3f(0.8f, 0.2f, 0.8f);
    draw_edge_with_arrow(pc, pb);
    viewer->glColor3f(0.0f, 0.0f, 0.0f);
    draw_point(pa);
    draw_point(pb);
  }
}

void R_s_k_2::draw_push_queue_stencil(const float point_size,
    const float line_width,
    const Point& query)
{
  Edge edge;
  bool ok = locate_edge(query, edge);
  if (!ok) return;

  Edge_vector hull;
  Edge_vector stencil;
  Vertex_handle src = m_dt.source_vertex(edge);
  m_dt.get_edges_from_star_minus_link(src, hull, true);
  collect_pqueue_stencil(m_dt, hull.begin(), hull.end(), stencil);

  Edge_vector::iterator it = stencil.begin();
  while (it != stencil.end())
  {
    Edge edge = *it;
    if (m_dt.source_vertex(edge) == src)
      it = stencil.erase(it);
    else if (m_dt.target_vertex(edge) == src)
      it = stencil.erase(it);
    else
      it++;
  }

  Triangulation copy;
  Edge_vector copy_hull;
  Edge_vector copy_stencil;
  Edge copy_edge = copy_star(edge, copy);
  Vertex_handle copy_src = copy.source_vertex(copy_edge);
  copy.get_edges_from_star_minus_link(copy_src, copy_hull, true);
  ok = copy.make_collapsible(copy_edge, copy_hull.begin(), copy_hull.end(), m_verbose);
  if (!ok) return;
  copy.collapse(copy_edge, m_verbose);
  collect_cost_stencil(copy, copy_hull.begin(), copy_hull.end(), copy_stencil);

  for (it = copy_stencil.begin(); it != copy_stencil.end(); ++it)
  {
    Edge edge = *it;
    Edge twin = copy.twin_edge(edge);
    if (!copy.is_pinned(edge)) stencil.push_back(edge);
    if (!copy.is_pinned(twin)) stencil.push_back(twin);
  }

  draw_bg_faces(copy, 0.7f, 0.7f, 1.0f, 1.0f);

  viewer->glLineWidth(0.5*line_width);
  for (it = stencil.begin(); it != stencil.end(); ++it)
  {
    Edge ab = *it;
    viewer->glColor3f(1.0f, 0.6f, 1.0f);
    draw_edge(ab);
  }

  viewer->glLineWidth(line_width);
  viewer->glPointSize(point_size);
  for (it = stencil.begin(); it != stencil.end(); ++it)
  {
    Edge ab = *it;
    Vertex_handle va = ab.first->vertex( (ab.second+1)%3 );
    Vertex_handle vb = ab.first->vertex( (ab.second+2)%3 );
    const Point& pa = va->point();
    const Point& pb = vb->point();
    Point pc = CGAL::midpoint(pa, pb);
    viewer->glColor3f(0.8f, 0.2f, 0.8f);
    draw_edge_with_arrow(pc, pb);
    viewer->glColor3f(0.0f, 0.0f, 0.0f);
    draw_point(pa);
    draw_point(pb);
  }
}

void R_s_k_2::draw_bg_faces(const Triangulation& mesh,
    const float red,
    const float green,
    const float blue,
    const float alpha)
{
  viewer->glEnable(GL_BLEND);
  viewer->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  viewer->glColor4f(red, green, blue, alpha);
  for (Finite_faces_iterator fi = mesh.finite_faces_begin(); fi != mesh.finite_faces_end(); ++fi)
  {
    Face_handle f = fi;
    draw_face(f);
  }
  viewer->glDisable(GL_BLEND);
}

void R_s_k_2::draw_bg_edges(const Triangulation& mesh,
    const float ri,
    const float gi,
    const float bi,
    const float ro,
    const float go,
    const float bo)
{
  for (Finite_faces_iterator fi = mesh.finite_faces_begin(); fi != mesh.finite_faces_end(); ++fi)
  {
    Face_handle f = fi;
    for (unsigned int i = 0; i < 3; ++i)
    {
      Edge e(f, i);
      e = mesh.twin_edge(e);
      if (mesh.is_infinite(e.first))
        viewer->glColor3f(ro, go, bo);
      else
        viewer->glColor3f(ri, gi, bi);
      draw_edge(e);
    }
  }
}

void R_s_k_2::draw_bg_vertices(const Triangulation& mesh,
    const float red,
    const float green,
    const float blue)
{
  viewer->glColor3f(red, green, blue);
  for (Finite_vertices_iterator vi = mesh.finite_vertices_begin(); vi != mesh.finite_vertices_end(); ++vi)
  {
    Vertex_handle v = vi;
    draw_point(v->point());
  }
}

void R_s_k_2::draw_vertex_faces(Vertex_handle vertex,
    const Triangulation& mesh,
    const float red,
    const float green,
    const float blue,
    const float alpha)
{
  viewer->glEnable(GL_BLEND);
  viewer->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  viewer->glColor4f(red, green, blue, alpha);
  Face_circulator fcirc = mesh.incident_faces(vertex);
  Face_circulator fend = fcirc;
  CGAL_For_all(fcirc, fend)
  {
    Face_handle f = fcirc;
    if (mesh.is_infinite(f)) continue;
    draw_face(f);
  }
  viewer->glDisable(GL_BLEND);
}

void R_s_k_2::draw_vertex_edges(Vertex_handle vertex,
    const Triangulation& mesh,
    const float ri,
    const float gi,
    const float bi,
    const float ro,
    const float go,
    const float bo)
{
  Face_circulator fcirc = mesh.incident_faces(vertex);
  Face_circulator fend = fcirc;
  CGAL_For_all(fcirc, fend)
  {
    Face_handle f = fcirc;
    int index = f->index(vertex);
    for (unsigned int i = 0; i < 3; ++i)
    {
      Edge e(f, i);
      if (mesh.is_infinite(e)) continue;
      if (static_cast<int>(i) == index) viewer->glColor3f(ro, go, bo);
      else viewer->glColor3f(ri, gi, bi);
      draw_edge(e);
    }
  }
}

void R_s_k_2::save_edges(std::ofstream& ofs, const int nb)
{
  MultiIndex mindex;
  for (Finite_edges_iterator ei = m_dt.finite_edges_begin(); ei != m_dt.finite_edges_end(); ++ei)
  {
    Edge edge = *ei;
    if (m_dt.is_ghost(edge)) continue;
    FT value = m_dt.get_edge_relevance(edge); // >= 0
    mindex.insert(PEdge(edge, value));
  }

  int nb_remove = (std::min)(nb, int(mindex.size()));
  for (int i = 0; i < nb_remove; ++i)
  {
    PEdge pedge = *(mindex.get<1>()).begin();
    (mindex.get<0>()).erase(pedge);

  }

  while (!mindex.empty())
  {
    PEdge pedge = *(mindex.get<1>()).begin();
    (mindex.get<0>()).erase(pedge);
    save_one_edge(ofs, pedge.edge());
  }
}

void R_s_k_2::save_one_edge(std::ofstream& ofs, const Edge& edge)
{
  int i = edge.second;
  Face_handle face = edge.first;
  Point const& a = face->vertex((i+1)%3)->point();
  Point const& b = face->vertex((i+2)%3)->point();
  ofs << a << " - " << b << std::endl;
}
