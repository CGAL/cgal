#ifndef _BVD_H_
#define _BVD_H_

#undef min
#undef max

// STL
#include <list>
#include <vector> 
// CGAL
#include <CGAL/basic.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/intersections.h> 
#include <CGAL/Polygon_2.h>
// local
#include "polygon_kernel.h"
#include "console_color.h"

template <class Kernel, class TDS>
class CBvd : public CGAL::Delaunay_triangulation_2<Kernel, TDS> {
public:
  typedef typename Kernel::FT FT;
  // type 2D
  typedef typename Kernel::Ray_2      Ray_2;
  typedef typename Kernel::Line_2     Line_2;
  typedef typename Kernel::Point_2    Point_2;
  typedef typename Kernel::Vector_2   Vector_2;
  typedef typename Kernel::Segment_2  Segment_2;
  typedef typename Kernel::Triangle_2 Triangle_2;
  typedef typename CGAL::Polygon_2<Kernel> Polygon_2;
  // type 3D
  typedef CGAL::Plane_3<Kernel>       Plane_3;
  typedef typename Kernel::Point_3    Point_3;
  typedef typename Kernel::Vector_3   Vector_3;
  typedef typename Kernel::Line_3     Line_3;
  typedef typename Kernel::Segment_3  Segment_3;
  typedef typename Kernel::Triangle_3 Triangle_3;
  typedef typename Kernel::Aff_transformation_3 Aff_transform_3;
  // Delaunay triangulation
  typedef CGAL::Delaunay_triangulation_2<Kernel, TDS> Dt;
  // handles
  typedef typename Dt::Face_handle   Face_handle;
  typedef typename Dt::Vertex_handle Vertex_handle;
  // iterators
  typedef typename Dt::Face_iterator Face_iterator;
  typedef typename Dt::Edge_iterator Edge_iterator;
  typedef typename Dt::Finite_vertices_iterator Finite_vertices_iterator;
  // circulators
  typedef typename Dt::Edge_circulator Edge_circulator;
  typedef typename Dt::Face_circulator Face_circulator;
  //containers
  typedef typename std::list<Point_3> Point_list;
  typedef typename Point_list::iterator Point_iter;
  typedef typename Point_list::const_iterator Point_const_iter;
  typedef typename std::list<CGAL::Color> Color_list;
  typedef typename Color_list::iterator Color_iter;
  typedef typename Color_list::const_iterator Color_const_iter;

  using Dt::incident_edges;
  using Dt::incident_faces;
  using Dt::is_infinite;
  using Dt::dual;
  using Dt::circumcenter;
public:
  //
  // life cyle
  //
  CBvd(const Triangle_3& triangle) {
    m_triangle_3d = triangle;
    m_to_3d = compute_transformation(m_triangle_3d.supporting_plane(),
      m_triangle_3d.vertex(0));
    m_to_2d = m_to_3d.inverse();
    m_triangle_2d = map_triangle_from_3d_to_2d(m_triangle_3d);
  }

  virtual ~CBvd() {
  }
  //
  // public API
  //
  void run(Point_list &point_list) {
    // point_list contains only the inner samples
    Point_iter pi;
    for (pi = point_list.begin(); pi != point_list.end(); pi++) {
      add_point(*pi);
    }
    point_list.clear();
    get_centroids(std::back_inserter(point_list));	// this is error_proven
  }

  void compute_voronoi_cells_and_boundaries(const Point_list &samples,
    const Vector_3 &normal, const Color_list &colors,
    std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
    std::vector<float> *pos_face_colors, std::vector<float> *pos_boundaries) {
    // step 1: build the map between points and colors
    std::map<Point_2, CGAL::Color> point_color_map;
    Point_const_iter pit;
    Color_const_iter cit;
    for (pit = samples.begin(), cit = colors.begin(); pit != samples.end();
      ++pit, ++cit) {
      add_point(*pit);
      Point_2 q = map_from_3d_to_2d(*pit);
      point_color_map.insert(std::make_pair(q, *cit));
    }
    // step 2: compute the voronoi cells and boundaries
    CGAL::Color color;
    for (Finite_vertices_iterator v = Dt::finite_vertices_begin();
      v != Dt::finite_vertices_end(); ++v) {
      const Point_2 &q = v->point();
      color = point_color_map[q];
      //::glColor3ub(color.r(), color.g(), color.b());
      compute_bounded_cell_and_boundaries(v, normal, color,
        pos_faces, pos_face_normals, pos_face_colors, pos_boundaries);
    }
  }

  void add_point(const Point_3& p) {
    Point_2 q = map_from_3d_to_2d(p);
    Dt::insert(q);
  }

  template <class OutputIterator> // value_type = Point_3
  bool get_centroids(OutputIterator out) const {
    if (Dt::dimension() < 2)
      return false;
    // compute centroids
    for (Finite_vertices_iterator v = Dt::finite_vertices_begin();
      v != Dt::finite_vertices_end(); v++) {
      Point_2 q = center_of_mass(v);
      Point_3 p = map_from_2d_to_3d(q);
      *out++ = p;
    }
    return true;
  }

  template <class OutputIterator> // value_type = FT
  bool get_cells_area(OutputIterator out) const {
    if (Dt::dimension() < 2)
      return false;
    FT area = m_triangle_3d.squared_area();
    for (Finite_vertices_iterator v = Dt::finite_vertices_begin();
      v != Dt::finite_vertices_end(); v++) {
      *out++ = cell_area(v) / area;
    }
    return true;
  }

private:
  //
  // 3D <-> 2D
  //
  Aff_transform_3 compute_transformation(const Plane_3& plane, const Point_3& point) const {
    Vector_3 translate = point - CGAL::ORIGIN;

    Vector_3 u = plane.base1();
    u = u / std::sqrt(u*u);
    Vector_3 v = plane.base2();
    v = v / std::sqrt(v*v);
    Vector_3 w = plane.orthogonal_vector();
    w = w / std::sqrt(w*w);

    return Aff_transform_3(u.x(), v.x(), w.x(), translate.x(),
      u.y(), v.y(), w.y(), translate.y(),
      u.z(), v.z(), w.z(), translate.z());
  }

  Triangle_2 map_triangle_from_3d_to_2d(const Triangle_3& triangle_3d) const {
    Point_2 a = map_from_3d_to_2d(triangle_3d.vertex(0));
    Point_2 b = map_from_3d_to_2d(triangle_3d.vertex(1));
    Point_2 c = map_from_3d_to_2d(triangle_3d.vertex(2));
    return Triangle_2(a, b, c);
  }

  Point_2 map_from_3d_to_2d(const Point_3& point_3d) const {
    Point_3 q = m_to_2d.transform(point_3d);
    return Point_2(q.x(), q.y());
  }

  Point_3 map_from_2d_to_3d(const Point_2& point_2d) const {
    Point_3 q(point_2d.x(), point_2d.y(), 0.0);
    return m_to_3d.transform(q);
  }
  //
  // Cell Area
  //
  FT cell_area(Vertex_handle v) const {
    if (has_inner_cell(v))
      return area_inner_cell(v);
    return area_bounded_cell(v);
  }

  bool has_inner_cell(Vertex_handle v) const {
    Face_circulator f = incident_faces(v);
    Face_circulator end = f;
    CGAL_For_all(f, end) {
      if (is_infinite(f))
        return false;
      if (is_outside(circumcenter(f)))
        return false;
    }
    return true;
  }

  bool is_inside(const Point_2& query) const {
    return m_triangle_2d.bounded_side(query) == CGAL::ON_BOUNDED_SIDE;
  }

  bool is_outside(const Point_2& query) const {
    return m_triangle_2d.bounded_side(query) == CGAL::ON_UNBOUNDED_SIDE;
  }

  FT area_inner_cell(Vertex_handle v) const {
    FT sum_areas = 0.0;
    Face_circulator f = incident_faces(v);
    Face_circulator end = f;
    Point_2 pivot = circumcenter(f); // pivot Voronoi vertex
    f++; f++;
    CGAL_For_all(f, end) {
      Point_2 p2 = circumcenter(f); f--;
      Point_2 p1 = circumcenter(f); f++;
      Triangle_2 triangle(pivot, p1, p2);
      sum_areas += triangle.area();
    }
    return sum_areas;
  }

  FT area_bounded_cell(Vertex_handle v) const {
    std::vector<Point_2> kernel;
    if (compute_bounded_cell(v, std::back_inserter(kernel)))
      return polygon_area(kernel);
    return 0.0; // never come here
  }

  FT polygon_area(std::vector<Point_2>& polygon) const {
    FT sum_areas = 0.0;
    Point_2 pivot = polygon[0]; // pivot vertex
    int nb_triangles = (int)polygon.size() - 2;
    for (int i = 0; i<nb_triangles; i++) {
      Point_2 p1 = polygon[(i + 1) % polygon.size()];
      Point_2 p2 = polygon[(i + 2) % polygon.size()];
      Triangle_2 triangle(pivot, p1, p2);
      sum_areas += std::abs(triangle.area());
    }
    return sum_areas;
  }
  //
  // Cell Center of Mass
  //
  Point_2 center_of_mass(Vertex_handle v) const {
    if (!is_inside(v->point())) {
      // std::cerr << red << "center_of_mass: point not strictly inside" << white << std::endl;
      return v->point();
    }
    if (has_inner_cell(v))
      return center_of_mass_inner_cell(v);
    else
      return center_of_mass_bounded_cell(v);
  }

  Point_2 center_of_mass_inner_cell(Vertex_handle v) const {
    FT sum_areas = 0.0;
    Vector_2 sum_vec = CGAL::NULL_VECTOR;
    Face_circulator f = incident_faces(v);
    Face_circulator end = f;
    Point_2 pivot = circumcenter(f); // pivot Voronoi vertex
    f++; f++;
    CGAL_For_all(f, end) {
      Point_2 p2 = circumcenter(f); f--;
      Point_2 p1 = circumcenter(f); f++;
      sum_areas += sum_centroid(Triangle_2(pivot, p1, p2), sum_vec);
    }
    if (sum_areas == 0.0) return pivot;
    return CGAL::ORIGIN + sum_vec / sum_areas;
  }

  Point_2 center_of_mass_bounded_cell(Vertex_handle v) const {
    std::vector<Point_2> kernel;
    if (compute_bounded_cell(v, std::back_inserter(kernel)))
      return polygon_center_of_mass(kernel);
    // never come here
    std::cerr << red << "bounded center of mass failed" << white << std::endl;
    return v->point();
  }

  template <class OutputIterator>
  bool compute_bounded_cell(Vertex_handle v, OutputIterator out) const {
    std::list<Segment_2> segments;
    Edge_circulator edge = incident_edges(v);
    Edge_circulator end = edge;
    CGAL_For_all(edge, end) {
      if (is_infinite(edge)) continue;
      CGAL::Object object = dual(edge);
      Segment_2 segment = convert_to_segment(object);
      segments.push_back(segment);
    }

    // add bounding domain
    segments.push_back(Segment_2(m_triangle_2d[0], m_triangle_2d[1]));
    segments.push_back(Segment_2(m_triangle_2d[1], m_triangle_2d[2]));
    segments.push_back(Segment_2(m_triangle_2d[2], m_triangle_2d[0]));

    // compute kernel of 2D polygon
    std::vector<Point_2> kernel;
    Polygon_kernel<Kernel> kernel_engine;
    return kernel_engine.run(segments.begin(), segments.end(), v->point(), out);
  }

  Point_2 polygon_center_of_mass(std::vector<Point_2>& polygon) const {
    if (polygon.size() < 3) {
      std::cerr << red << "polygon size: " << polygon.size() << white << std::endl;
      return polygon[0];
    }

    FT sum_areas = 0.0;
    Vector_2 sum_vec = CGAL::NULL_VECTOR;
    Point_2 pivot = polygon[0]; // pivot vertex
    int nb_triangles = (int)polygon.size() - 2;
    for (int i = 0; i<nb_triangles; i++) {
      Point_2 p1 = polygon[(i + 1) % polygon.size()];
      Point_2 p2 = polygon[(i + 2) % polygon.size()];
      Triangle_2 triangle(pivot, p1, p2);
      Point_2 centroid = CGAL::centroid(triangle);
      FT area = std::abs(triangle.area());
      sum_vec = sum_vec + area * (centroid - CGAL::ORIGIN);
      sum_areas += area;
    }
    if (sum_areas == 0.0) return pivot;
    return CGAL::ORIGIN + sum_vec / sum_areas;
  }

  FT sum_centroid(const Triangle_2& triangle, Vector_2& vec_sum) const {
    Point_2 centroid = CGAL::centroid(triangle); // assuming uniform density
    FT area = std::abs(triangle.area());
    vec_sum = vec_sum + area * (centroid - CGAL::ORIGIN);
    return area;
  }
  //
  // Visual elements
  // 
  void compute_bounded_cell_and_boundaries(Vertex_handle v,
    const Vector_3 &normal, const CGAL::Color &color,
    std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
    std::vector<float> *pos_face_colors, std::vector<float> *pos_boundaries) const {
    std::vector<Point_2> cell;
    if (compute_bounded_cell(v, std::back_inserter(cell))) {
      if (cell.size() < 3) {
        std::cerr << red << "polygon with degree < 3" << white << std::endl;
        return;
      }
      compute_polygon(cell, normal, color,
        pos_faces, pos_face_normals, pos_face_colors);
      compute_polygon_boundary(cell, pos_boundaries);
    }
  }

  void compute_polygon_boundary(const std::vector<Point_2> &polygon,
    std::vector<float> *pos_boundaries) const {
    int size = static_cast<int>(polygon.size());
    for (int i = 0; i < size; ++i) {
      int j = (i + 1) % size;
      Point_3 start = map_from_2d_to_3d(polygon[i]);
      Point_3 end = map_from_2d_to_3d(polygon[j]);
      compute_point(start, pos_boundaries);
      compute_point(end, pos_boundaries);
    }
  }

  void compute_polygon(const std::vector<Point_2> &polygon,
    const Vector_3 &normal, const CGAL::Color &color,
    std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
    std::vector<float> *pos_face_colors) const {
    std::vector<Point_2> points(3, polygon[0]);
    int nb_triangles = static_cast<int>(polygon.size() - 2);
    for (int i = 0; i < nb_triangles; i++) {
      Point_2 p1 = polygon[(i + 1) % polygon.size()];
      Point_2 p2 = polygon[(i + 2) % polygon.size()];
      points[1] = p1;
      points[2] = p2;
      compute_triangle(points, normal, color,
        pos_faces, pos_face_normals, pos_face_colors);
    }
  }

  void compute_triangle(const std::vector<Point_2> &points,
    const Vector_3 &normal, const CGAL::Color &color,
    std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
    std::vector<float> *pos_face_colors) const {
    for (const Point_2 &p_2 : points) {
      Point_3 p_3 = map_from_2d_to_3d(p_2);
      compute_point(p_3, pos_faces);
      compute_normal(normal, pos_face_normals);
      compute_color(color, pos_face_colors);
    }
  }

  void inline compute_point(const Point_3 &p,
    std::vector<float> *pos_faces) const {
    pos_faces->push_back(p.x());
    pos_faces->push_back(p.y());
    pos_faces->push_back(p.z());
  }

  void inline compute_normal(const Vector_3 &normal,
    std::vector<float> *pos_face_normals) const {
    pos_face_normals->push_back(normal.x());
    pos_face_normals->push_back(normal.y());
    pos_face_normals->push_back(normal.z());
  }

  void inline compute_color(const CGAL::Color &color,
    std::vector<float> *pos_face_colors) const {
    pos_face_colors->push_back(color.red() / 255.0f);
    pos_face_colors->push_back(color.green() / 255.0f);
    pos_face_colors->push_back(color.blue() / 255.0f);
  }
  //
  // Draw
  //
  void draw_bounded_cell(Vertex_handle v) const {
    std::vector<Point_2> cell;
    if (compute_bounded_cell(v, std::back_inserter(cell))) {
      draw_polygon(cell);
    }
  }

  void draw_bounded_cell_boundary(Vertex_handle v) const {
    std::vector<Point_2> cell;
    if (compute_bounded_cell(v, std::back_inserter(cell)))
      draw_polygon_boundary(cell);
  }

  void draw_polygon(const std::vector<Point_2>& polygon) const {
    if (polygon.size() < 3) {
      std::cerr << red << "polygon with degree < 3" << white << std::endl;
      return;
    }
    Point_2 p0 = polygon[0]; // pivot 
    int nb_triangles = int(polygon.size() - 2);
    for (int i = 0; i < nb_triangles; i++) {
      Point_2 p1 = polygon[(i + 1) % polygon.size()];
      Point_2 p2 = polygon[(i + 2) % polygon.size()];
      draw_triangle(p0, p1, p2);
    }
  }

  void draw_polygon_boundary(const std::vector<Point_2>& polygon) const {
    if (polygon.size() < 3) {
      std::cerr << red << "polygon with degree < 3" << white << std::endl;
      return;
    }
    glBegin(GL_LINE_LOOP);
    for (size_t i = 0; i < polygon.size(); ++i)
      gl_vertex(polygon[i]);
    glEnd();
  }

  void draw_inner_cell(Vertex_handle v) const {
    Face_circulator f = incident_faces(v);
    Face_circulator end = f;
    Point_2 p0 = circumcenter(f); // pivot Voronoi vertex
    f++; f++;
    CGAL_For_all(f, end) {
      Point_2 p2 = circumcenter(f); f--;
      Point_2 p1 = circumcenter(f); f++;
      draw_triangle(p0, p1, p2);
    }
  }

  void draw_triangle(const Point_2& a, const Point_2& b, const Point_2& c) const {
    ::glBegin(GL_TRIANGLES);
    gl_vertex(a);
    gl_vertex(b);
    gl_vertex(c);
    ::glEnd();
  }

  void clip_and_draw(const CGAL::Object& object) const {
    Segment_2 segment;
    if (CGAL::assign(segment, object)) {
      clip_and_draw(segment);
      return;
    }
    Ray_2 ray;
    if (CGAL::assign(ray, object)) {
      clip_and_draw(ray);
      return;
    }
    Ray_2 line;
    if (CGAL::assign(line, object)) {
      clip_and_draw(line);
      return;
    }
  }

  void clip_and_draw(const Ray_2& ray) const {
    std::vector<Point_2> intersections;
    if (is_inside(ray.source()) &&
      intersect_domain<Ray_2>(ray, intersections)) {
      gl_vertex(ray.source());
      gl_vertex(intersections[0]);
      return;
    }
    // source is outside
    if (intersect_domain<Ray_2>(ray, intersections) &&
      intersections.size() >= 2) {
      gl_vertex(intersections[0]);
      gl_vertex(intersections[1]); // FIXME: degenerate cases
    }
  }

  void clip_and_draw(const Line_2& line) const {
    std::vector<Point_2> intersections;
    if (intersect_domain<Line_2>(line, intersections) &&
      intersections.size() >= 2) {
      gl_vertex(intersections[0]);
      gl_vertex(intersections[1]); // FIXME: degenerate cases
    }
  }

  void clip_and_draw(const Segment_2& segment) const {
    // two are inside (easy)
    if (is_inside(segment.source()) &&
      is_inside(segment.target())) {
      gl_vertex(segment.source());
      gl_vertex(segment.target());
      return;
    }
    std::vector<Point_2> intersections;
    if (is_inside(segment.source()) &&
      intersect_domain<Segment_2>(segment, intersections)) {
      gl_vertex(segment.source());
      gl_vertex(intersections[0]);
      return;
    }

    if (is_inside(segment.target()) &&
      intersect_domain<Segment_2>(segment, intersections)) {
      gl_vertex(segment.target());
      gl_vertex(intersections[0]);
      return;
    }

    // two are outside
    if (intersect_domain<Segment_2>(segment, intersections) &&
      intersections.size() >= 2) {
      gl_vertex(intersections[0]);
      gl_vertex(intersections[1]); // FIXME: degenerate cases
    }
  }

  void gl_vertex(const Point_2& q) const {
    Point_3 p = map_from_2d_to_3d(q);
    ::glVertex3d(p.x(), p.y(), p.z());
  }

  template <class Query> // Segment_2, Ray_2 or Ray_2
  bool intersect_domain(const Query& query,
    std::vector<Point_2>& intersections) const {
    for (int i = 0; i<3; i++) {
      Point_2 intersection;
      Segment_2 segment(m_triangle_2d[i], m_triangle_2d[(i + 1) % 3]);
      CGAL::Object object = CGAL::intersection(query, segment);
      if (CGAL::assign(intersection, object))
        intersections.push_back(intersection);
    }
    return !intersections.empty();
  }

  Segment_2 convert_to_segment(const CGAL::Object& object) const {
    Segment_2 segment;
    if (CGAL::assign(segment, object))
      return segment;

    Ray_2 ray;
    if (CGAL::assign(ray, object))
      return Segment_2(ray.source(), ray.point(100));

    Line_2 line;
    if (CGAL::assign(line, object))
      return Segment_2(line.point(-100), line.point(100));

    // never come here
    std::cerr << red << "convert_to_segment failed" << white << std::endl;
    return Segment_2(CGAL::ORIGIN, CGAL::ORIGIN);
  }

private:
  Triangle_2 m_triangle_2d;
  Triangle_3 m_triangle_3d;
  Aff_transform_3 m_to_2d;
  Aff_transform_3 m_to_3d;
};

#endif

/*CGAL::Color get_rgb_color(double h, double s, double v) {
  // h(1~360), s(0~1), v(0~1)
  int hi = h / 60;
  double f = h / 60 - hi;
  double p = v * (1.0 - s);
  double q = v * (1.0 - f * s);
  double t = v * (1.0 - (1 - f) * s);
  v *= 255;
  t *= 255;
  p *= 255;
  q *= 255;
  switch (hi) {
  case 0:
    return CGAL::Color(v, t, p);
  case 1:
    return CGAL::Color(q, v, p);
  case 2:
    return CGAL::Color(p, v, t);
  case 3:
    return CGAL::Color(p, q, v);
  case 4:
    return CGAL::Color(t, p, v);
  case 5:
    return CGAL::Color(v, p, q);
  default:
    return CGAL::Color(0, 0, 0);
  }
}

void draw_voronoi_cells(Point_list &point_list,
  Color_list &colors) {
  // step 1: build the map between points and colors
  std::map<Point_2, CGAL::Color> point_to_color;
  Point_list_iterator pit;
  Colot_list_iterator cit;
  for (pit = point_list.begin(), cit = colors.begin();
    pit != point_list.end(); ++pit, ++cit) {
    add_point(*pit);
    Point_2 q = map_from_3d_to_2d(*pit);
    point_to_color.insert(std::make_pair(q, *cit));
  }
  // step 2: draw the voronoi cells
  CGAL::Color color;
  for (Finite_vertices_iterator v = Dt::finite_vertices_begin();
    v != Dt::finite_vertices_end(); ++v) {
    const Point_2 &q = v->point();
    //Point_3 p = map_from_2d_to_3d(q);
    color = point_to_color[q];
    ::glColor3ub(color.r(), color.g(), color.b());
    draw_bounded_cell(v);				// draw the cell
  }
}

void draw_voronoi_boundaries(Point_list &point_list) {
  Point_list_iterator pi;
  for (pi = point_list.begin(); pi != point_list.end(); pi++)
    add_point(*pi);
  for (Finite_vertices_iterator v = Dt::finite_vertices_begin();
    v != Dt::finite_vertices_end(); v++)
    draw_bounded_cell_boundary(v);
}

bool get_cells_areas(std::map<Point_3, FT, Point_Comp> &points_to_areas) {
  if (Dt::dimension() < 2)
    return false;
  FT area = m_triangle_3d.squared_area();
  for (Finite_vertices_iterator v = Dt::finite_vertices_begin();
    v != Dt::finite_vertices_end(); ++v) {
    const Point_2 q = v->point();
    Point_3 p = map_from_2d_to_3d(q);
    points_to_areas.insert(std::pair<Point_3, FT>(p, cell_area(v) / area));
  }
  return true;
}

template <class OutputIterator> // value_type = Point_3
bool get_centroids(std::set<Point_3, Point_Comp> &disturbed_border_points,
  std::map<Point_2, Point_3> &point_map,
  OutputIterator out) const {
  if (Dt::dimension() < 2)
    return false;
  // compute centroids
  for (Finite_vertices_iterator v = Dt::finite_vertices_begin();
    v != Dt::finite_vertices_end(); v++) {
    Point_2 q = v->point();
    // Point_3 p = map_from_2d_to_3d(q);
    Point_3 p = point_map[q];
    if (disturbed_border_points.find(p) == disturbed_border_points.end()) {
      q = center_of_mass(v);
      p = map_from_2d_to_3d(q);
    }
    *out++ = p;
  }
  return true;
}

void run(Point_list &all_points,
  std::set<Point_3, Point_Comp> &disturbed_border_points) {
  // relocate samples, all_points contains the inner samples,
  // disturbed_boundary samples will keep fixed during relocation
  std::map<Point_2, Point_3> point_map;
  for (Point_list_iterator pi = all_points.begin();
    pi != all_points.end(); pi++) {
    add_point(*pi);
    point_map[map_from_3d_to_2d(*pi)] = *pi;
  }
  all_points.clear();
  get_centroids(disturbed_border_points, point_map,
    std::back_inserter(all_points));	// error_proven
}

void run(Point_list &point_list, std::map<Point_3, FT,
  Point_Comp> &points_to_areas) {
  // calculate cell areas
  Point_iter pi;
  for (pi = point_list.begin(); pi != point_list.end(); pi++) {
    add_point(*pi);
  }
  points_to_areas.clear();
  get_cells_areas(points_to_areas);
}
*/
