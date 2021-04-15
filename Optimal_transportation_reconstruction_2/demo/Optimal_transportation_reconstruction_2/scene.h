#ifndef SCENE_H_
#define SCENE_H_

// STL
#include <fstream>

//Qt
#include <QtOpenGL>
#include <QWidget>

// local
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "Otr2_kerneled.h"
#include <CGAL/Optimal_transportation_reconstruction_2.h>


#ifdef CGAL_USE_CIMG
#define cimg_display 0 // To avoid X11 or Windows-GDI dependency
#include <CImg.h>
#endif
#include <utility>      // std::pair
#include <vector>

#include <CGAL/number_type_config.h>
#include <CGAL/Random.h>
#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/functional.h>
#include <CGAL/Iterator_range.h>

class GlViewer;

class Scene {

public:
  typedef std::pair<Point, FT> PointMassPair;

  typedef std::vector<PointMassPair> PointMassList;
  typedef PointMassList::const_iterator InputIterator;

  typedef CGAL::value_type_traits<InputIterator>::type MassPoint;

  typedef CGAL::First_of_pair_property_map<PointMassPair> Point_property_map;
  typedef CGAL::Second_of_pair_property_map<PointMassPair> Mass_property_map;

  typedef CGAL::Optimal_transportation_reconstruction_2<K, Point_property_map,
    Mass_property_map> R_s_2;

  typedef K::Segment_2 Segment;

  typedef R_s_2::Vector Vector;

  typedef R_s_2::Vertex Vertex;
  typedef R_s_2::Vertex_handle Vertex_handle;
  typedef R_s_2::Vertex_iterator Vertex_iterator;
  typedef R_s_2::Vertex_circulator Vertex_circulator;
  typedef R_s_2::Finite_vertices_iterator Finite_vertices_iterator;

  typedef R_s_2::Edge Edge;
  typedef R_s_2::Edge_circulator Edge_circulator;
  typedef R_s_2::Finite_edges_iterator Finite_edges_iterator;

  typedef R_s_2::Face_handle Face_handle;
  typedef R_s_2::Face_circulator Face_circulator;
  typedef R_s_2::Finite_faces_iterator Finite_faces_iterator;

  typedef R_s_2::Vertex_handle_map Vertex_handle_map;
  typedef R_s_2::Face_handle_map Face_handle_map;

  typedef R_s_2::Vertex_handle_set Vertex_handle_set;
  typedef R_s_2::Edge_set Edge_set;

  typedef R_s_2::Edge_vector Edge_vector;

  typedef R_s_2::Sample_ Sample_;
  typedef R_s_2::Sample_vector Sample_vector;
  typedef R_s_2::Sample_vector_const_iterator Sample_vector_const_iterator;

  typedef R_s_2::PSample PSample;
  typedef R_s_2::SQueue SQueue;

  typedef R_s_2::Rec_edge_2 PEdge;


private:

  struct Point_3_from_sample : public CGAL::cpp98::unary_function<Sample_, K::Point_3>
  {
    K::Point_3 operator() (const Sample_& sample) const
    {
      return K::Point_3 (sample.point().x(), sample.point().y(), 0.);
    }
  };

  // data
  std::vector<Sample_> m_samples;

  Optimal_transportation_reconstruction_kerneled_2* m_pwsrec;
  int m_ignore;
  bool m_init_done;
  bool is_viewer_set;
  double m_percentage;

  // bbox
  double m_bbox_x;
  double m_bbox_y;
  double m_bbox_size;

  //Random
  CGAL::Random random;

  template <class Vector>
  Vector random_vec(const double scale)
  {
    double dx = random.get_double(-scale, scale);
    double dy = random.get_double(-scale, scale);
    return Vector(dx, dy);
  }

public:
  Scene() {
    srand(0); // for sake of repeatability
    m_ignore = 0;
    m_init_done = false;
    m_percentage = 1.;
    m_bbox_x = 0.0;
    m_bbox_y = 0.0;
    m_bbox_size = 1.0;
    is_viewer_set = false;

    m_pwsrec = new Optimal_transportation_reconstruction_kerneled_2();
  }

  ~Scene() {
    delete m_pwsrec;
  }

  void clear() {
    delete m_pwsrec;
    m_pwsrec = new Optimal_transportation_reconstruction_kerneled_2();
    m_samples.clear();
  }



  void subdivide() {
    if (m_samples.size() < 3)
      return;

    std::vector<Sample_> new_samples;
    std::vector<Sample_>::const_iterator it = m_samples.begin();
    std::vector<Sample_>::const_iterator last = it++;
    while (it != m_samples.end()) {
      Point p = CGAL::midpoint(last->point(), it->point());
      FT m = 0.5 * (last->mass() + it->mass());
      new_samples.push_back(Sample_(p, m));
      last = it++;
    }
    it = m_samples.begin();
    Point p = CGAL::midpoint(last->point(), it->point());
    FT m = 0.5 * (last->mass() + it->mass());
    new_samples.push_back(Sample_(p, m));

    std::vector<Sample_> final_samples;
    std::vector<Sample_>::const_iterator it2 = new_samples.begin();
    while (it != m_samples.end() && it2 != new_samples.end()) {
      final_samples.push_back(*it);
      final_samples.push_back(*it2);
      it++;
      it2++;
    }

    m_samples = final_samples;
  }

  // SAMPLE //

  void add_sample(const Point& point, const FT mass = 1.0) {
    m_samples.push_back(Sample_(point, mass));
    m_init_done = false;
  }

  void add_outliers(const unsigned int nb) {
    std::cerr << "adding " << nb << " outliers...";
    for (unsigned int i = 0; i < nb; i++) {
      Point outlier = CGAL::ORIGIN + random_vec<Vector>(1.3);
      m_samples.push_back(outlier);
    }
    m_init_done = false;
    std::cerr << "done" << std::endl;
  }

  void noise(const FT scale) {
    std::cerr << "noising by " << scale << "...";
    std::vector<Sample_>::iterator it;
    for (it = m_samples.begin(); it != m_samples.end(); it++) {
      Sample_& sample = *it;
      Point& point = sample.point();
      point = point + random_vec<Vector>(scale);
    }
    std::cerr << "done" << std::endl;
  }

  void normalize_points() {
    noise(1e-5);
    compute_bbox(m_bbox_x, m_bbox_y, m_bbox_size);
    if (m_bbox_size == 0.0)
      return;

    Point center(m_bbox_x, m_bbox_y);
    std::vector<Sample_>::iterator it;
    for (it = m_samples.begin(); it != m_samples.end(); ++it) {
      Sample_& sample = *it;
      Vector vec = (sample.point() - center) / m_bbox_size;
      sample.point() = CGAL::ORIGIN + vec;
    }
    m_bbox_x = m_bbox_y = 0.0;
    m_bbox_size = 1.0;
  }

  void compute_bbox(double &x, double &y, double &scale) {
    if (m_samples.empty()) {
      x = y = 0.0;
      scale = 1.0;
      return;
    }

    FT x_min, x_max, y_min, y_max;
    std::vector<Sample_>::const_iterator it = m_samples.begin();
    Point p = it->point();
    x_min = x_max = p.x();
    y_min = y_max = p.y();
    ++it;
    for (; it != m_samples.end(); ++it) {
      p = it->point();
      x_min = (std::min)(x_min, p.x());
      x_max = (std::max)(x_max, p.x());
      y_min = (std::min)(y_min, p.y());
      y_max = (std::max)(y_max, p.y());
    }

    x = 0.5 * (x_min + x_max);
    y = 0.5 * (y_min + y_max);
    scale = (std::max)(x_max - x_min, y_max - y_min);
    if (scale == 0.0)
      scale = 1.0;
  }

  // IO SAMPLES //

  void load(const QString& filename, QWidget* qw) {

    if (filename.contains(".xyz", Qt::CaseInsensitive)) {
      load_xyz_file(filename);
      //      normalize_points();
      return;
    }
    if (filename.contains(".xyw", Qt::CaseInsensitive)) {
      load_xyw_file(filename);
      //      normalize_points();
      return;
    }
     if (filename.contains(".xy", Qt::CaseInsensitive)) {
      load_xy_file(filename);
      //      normalize_points();
      return;
    }

#ifdef CGAL_USE_CIMG
    if (filename.contains(".bmp", Qt::CaseInsensitive)) {
      bool use_gradient = false;

      QMessageBox::StandardButton reply;
      reply = QMessageBox::question(qw, QString("Open BMP"), "Use gradient?",
        QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);
      if (reply == QMessageBox::Yes)
        use_gradient = true;
      else if (reply == QMessageBox::No)
        use_gradient = false;
      else
        return;

      if (use_gradient)
        load_gradient(filename);
      else
        load_image(filename);
      normalize_points();
      return;
    }

    std::cerr << "Invalid file (try .xy, .bmp)" << std::endl;
#else
    CGAL_USE(qw);
    std::cerr << "Invalid file (try .xy)" << std::endl;
#endif


  }

  void load_xy_file(const QString& fileName) {

    std::cout << "filename: " << fileName.toUtf8().constData() << std::endl;
    std::ifstream ifs(qPrintable(fileName));
    std::cerr << "reading xy...";
    Point point;
    unsigned int nb = 0;
    while (ifs >> point) {
      add_sample(point, 1.0);
      nb++;
    }
    std::cerr << "done (" << nb << " points)" << std::endl;
    ifs.close();
  }

  void load_xyw_file(const QString& fileName) {

    std::cout << "filename: " << fileName.toUtf8().constData() << std::endl;
    std::ifstream ifs(qPrintable(fileName));
    std::cerr << "reading xy with weights...";
    Point point;
    FT weight;
    unsigned int nb = 0;
    while (ifs >> point >> weight) {
      add_sample(point, weight);
      nb++;
    }
    std::cerr << "done (" << nb << " points)" << std::endl;
    ifs.close();
    compute_average_spacing();
  }

  void load_xyz_file(const QString& fileName) {

    std::cout << "filename: " << fileName.toUtf8().constData() << std::endl;
    std::ifstream ifs(qPrintable(fileName));
    std::cerr << "reading xyz...";
    unsigned int nb = 0;
    std::string str;
    while (getline (ifs, str)) {
      std::istringstream iss (str);
      double x = 0., y = 0.;
      iss >> x >> y;
      str.clear();
      add_sample(Point (x, y), 1.0);
      nb++;
    }
    std::cerr << "done (" << nb << " points)" << std::endl;
    ifs.close();
  }


#ifdef CGAL_USE_CIMG
  void load_image(const QString& fileName) {
    std::cerr << "reading image...";
    cimg_library::CImg<float> image(qPrintable(fileName));
    std::cerr << "done" << std::endl;

    std::cerr << "computing grayscale...";
    cimg_library::CImg<float> grayscale =
      image.RGBtoHSV().get_channel(2).normalize(0.0f, 1.0f);
    std::cerr << "done" << std::endl;

    // turn pixels into weighted samples
    std::cerr << "adding samples...";
    for (int i = 0; i < grayscale.width(); i++) {
      for (int j = 0; j < grayscale.height(); j++) {
        float mass = 1.0f - grayscale.atXY(i, j);
        double x = double(i) / grayscale.width();
        double y = 1.0 - double(j) / grayscale.height();
        if (mass > 0.f)
          add_sample(Point(x, y), mass);
      }
    }
    std::cerr << "done (" << m_samples.size() << ")" << std::endl;
  }

  void load_gradient(const QString& fileName) {
    std::cerr << "reading image...";
    cimg_library::CImg<float> image(qPrintable(fileName));
    std::cerr << "done" << std::endl;

    std::cerr << "computing gradient...";
    cimg_library::CImgList<float> grad = image.get_gradient();
    cimg_library::CImg<float> normgrad = sqrt(
      grad[0].pow(2) + grad[1].pow(2)).normalize(0.0f, 1.0f);
    std::cerr << "done" << std::endl;

    // turn pixels into weighted samples
    std::cerr << "adding samples...";
    for (int i = 0; i < normgrad.width(); i++) {
      for (int j = 0; j < normgrad.height(); j++) {
        float mass = normgrad.atXY(i, j);
        double x = double(i) / normgrad.width();
        double y = 1.0 - double(j) / normgrad.height();
        add_sample(Point(x, y), mass);
      }
    }
    std::cerr << "done (" << m_samples.size() << ")" << std::endl;
  }
#endif

  void compute_average_spacing()
  {
    FT spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
      (CGAL::make_range (boost::make_transform_iterator (m_samples.begin(),
                                                         Point_3_from_sample()),
                         boost::make_transform_iterator (m_samples.end(),
                                                         Point_3_from_sample())),
       3, CGAL::parameters::point_map (CGAL::Identity_property_map_no_lvalue<K::Point_3>()));
    std::cerr << "Average spacing = " << spacing << std::endl;
  }

  void print_vertex(Vertex vertex) {
    std::cout << "vertex " << vertex << std::endl;
  }


  void print_edge(PEdge edge) {
    int i = ((edge).edge()).second;
    Point a = ((edge).edge()).first->vertex((i + 1) % 3)->point();
    Point b = ((edge).edge()).first->vertex((i + 2) % 3)->point();
    std::cout << "( " << (edge).priority() << ") ( " << a
      << " , " << b << " )" << std::endl;
  }


  void debug_print()
  {
    std::vector<Point> isolated_points;
    std::vector<Segment> edges;

    m_pwsrec->list_output(std::back_inserter(isolated_points), std::back_inserter(edges));

    int vertex_count = 0;
    for (std::vector<Point>::iterator it = isolated_points.begin();
      it != isolated_points.end(); it++) {
      vertex_count++;
      std::cout << *it << std::endl;
    }
    CGAL_assertion(vertex_count == 18);

    int edge_count = 0;
    for (std::vector<Segment>::iterator it = edges.begin();
      it != edges.end(); it++) {
      std::cout << *it << std::endl;
      edge_count++;
    }
  }

  void save(const QString& filename)
  {
    Sample_vector samples;
    for (std::vector<Sample_>::iterator it = m_samples.begin();
      it != m_samples.end(); ++it) {
      Sample_& s = *it;
      samples.push_back(&s);
    }

    if (filename.contains(".xy", Qt::CaseInsensitive)) {
      save_xy(filename, samples);
      return;
    }

    std::cerr << "Error: not an XY file." << std::endl;
  }



  void save_xy(const QString& filename, const Sample_vector& samples) {
    std::ofstream ofs(qPrintable(filename));
    for (Sample_vector_const_iterator it = samples.begin();
      it != samples.end(); ++it) {
      Sample_* sample = *it;
      ofs << sample->point() << std::endl;
    }
    ofs.close();
  }


  // RECONSTRUCTION //

  void set_options(const int verbose, const int mchoice,
    const bool use_flip, const unsigned int relocation,
    const double ghost, const double percentage) {

    m_pwsrec->set_verbose(verbose);
    m_pwsrec->set_random_sample_size(mchoice);
    m_pwsrec->set_use_flip(use_flip);
    m_pwsrec->set_relocation(relocation);
    m_pwsrec->set_relevance(ghost);
    m_percentage = percentage;
  }

  bool init_reconstruction(const double percentage) {
    std::cout << " init_reconstruction with " << (unsigned int)(100.*percentage) << "% of points" << std::endl;

    if (m_samples.empty()) {
      std::cerr << "initialization failed (empty point set)" << std::endl;
      return false;
    }

    Sample_vector vertices, samples;
    select_samples(percentage, vertices, samples);

    PointMassList vertices_mass_list;
    Sample_vector_const_iterator it;
    for (it = vertices.begin(); it != vertices.end(); it++) {
      vertices_mass_list.push_back(
        std::make_pair((*it)->point(), (*it)->mass()));
    }
    PointMassList samples_mass_list;
    for (it = samples.begin(); it != samples.end(); it++) {
      samples_mass_list.push_back(
        std::make_pair((*it)->point(), (*it)->mass()));
    }

    Point_property_map point_pmap;
    Mass_property_map mass_pmap;
    MassPoint mp;

    m_pwsrec->initialize_with_custom_vertices
      (samples_mass_list.begin(), samples_mass_list.end(),
       vertices_mass_list.begin(), vertices_mass_list.end(),
       point_pmap, mass_pmap);

    m_init_done = true;

    return true;
  }

  void decimate(const double percentage) {
    std::cout << "decimating from " << m_samples.size() << " to...";
    std::vector<Sample_> selected;

    std::vector<Sample_>::iterator it;
    for (it = m_samples.begin(); it != m_samples.end(); it++) {
      const double rd = random.get_double(0.0, 1.0);
      if (rd >= percentage)
        selected.push_back(*it);
    }

    m_samples.clear();
    std::copy(selected.begin(), selected.end(),
      std::back_inserter(m_samples));
    std::cout << m_samples.size() << std::endl;
  }


  void select_samples(const double percentage, Sample_vector& vertices,
    Sample_vector& samples) {
    std::vector<Sample_>::iterator it;
    for (it = m_samples.begin(); it != m_samples.end(); ++it) {
      Sample_& s = *it;

      samples.push_back(&s);
      FT rv = random.get_double(0.0, 1.0);
      if (rv <= percentage)
        vertices.push_back(&s);
    }
  }

  void reconstruct_until(const unsigned int nv) {
    std::cout << "reconstruct_until" << std::endl;
    if (!m_init_done)
      init_reconstruction(m_percentage);
    m_pwsrec->run_until(nv);
  }

  void reconstruct_wasserstein_tolerance (const double tolerance) {
    std::cout << "reconstruct_wasserstein_tolerance" << std::endl;
    if (!m_init_done)
      init_reconstruction(m_percentage);
    m_pwsrec->run_under_wasserstein_tolerance (tolerance);
  }

  void reconstruct(const unsigned int steps) {
    std::cout << "reconstruct" << std::endl;
    if (!m_init_done)
      init_reconstruction(m_percentage);
    m_pwsrec->run(steps);
  }

  void relocate_all_points() {
    std::cout << "relocate_all_points" << std::endl;
    m_pwsrec->relocate_all_points();
  }

  void output_console()
  {
    std::cout << std::endl;
    std::cout << "OFF OUTPUT" << std::endl;
    std::vector<Point> points;
    std::vector<std::size_t> isolated_vertices;
    std::vector<std::pair<std::size_t, std::size_t> > edges;

    m_pwsrec->indexed_output(
      std::back_inserter(points),
      std::back_inserter(isolated_vertices),
      std::back_inserter(edges));

    std::cout << "OFF " << points.size() << " 0 " << edges.size() << std::endl;

    // points
    std::vector<Point>::iterator pit;
    for (pit = points.begin(); pit != points.end(); pit++)
      std::cout << *pit << std::endl;

    // isolated vertices
    std::vector<std::size_t>::iterator vit;
    for (vit = isolated_vertices.begin(); vit != isolated_vertices.end(); vit++)
      std::cout << "1 " << *vit << std::endl;

    // edges
    std::vector<std::pair<std::size_t, std::size_t> >::iterator eit;
    for (eit = edges.begin(); eit != edges.end(); eit++)
      std::cout << "2 " << eit->first << " " << eit->second << std::endl;
  }

  // RENDER //

  void render(const bool view_points, const bool view_tolerance,
    const bool view_vertices,
    const bool view_edges, const bool view_ghost_edges,
    const bool view_edge_cost, const bool view_edge_priority,
    const bool view_bins, const bool view_foot_points,
    const bool view_relocation, const bool view_edge_relevance,
    const float point_size, const float vertex_size,
    const float line_thickness, GlViewer* viewer)
  {
    if (m_pwsrec == NULL) {
      return;
    }
    if(!is_viewer_set)
    {
      m_pwsrec->setViewer(viewer);
      is_viewer_set = true;
    }

    if (view_tolerance)
      draw_tolerance(viewer);

    if (view_edges)
      m_pwsrec->draw_edges(0.5f * line_thickness, 0.9f, 0.9f, 0.9f);

    if (view_edge_cost)
      m_pwsrec->draw_costs(line_thickness, view_ghost_edges);

    if (view_edge_priority)
      m_pwsrec->draw_pedges(line_thickness);

    if (view_edge_relevance)
      m_pwsrec->draw_relevance(line_thickness, m_ignore);

    if (view_relocation)
      m_pwsrec->draw_relocation();

    if (view_vertices)
      m_pwsrec->draw_vertices(vertex_size, 0.0f, 0.0f, 0.5f);

    if (view_bins)
      m_pwsrec->draw_bins(0.5f * line_thickness);

    if (view_foot_points)
      m_pwsrec->draw_footpoints(line_thickness, 0.2f, 0.8f, 0.2f);

    if (view_points)
      draw_samples(point_size, viewer);
  }

  void draw_samples(const float point_size, GlViewer* viewer);

  void draw_tolerance(GlViewer* viewer);

  // PREDEFINED EXAMPLES //

  void make_line(const unsigned int nb, const Point& start,
    const Point& end) {
    Point curr = start;
    Vector incr = (end - start) / nb;
    for (unsigned int i = 0; i < nb; i++) {
      add_sample(curr);
      curr = curr + incr;
    }
  }

  void make_circle_arc(const unsigned int nb, const Point& c,
    const double radius, const double min_angle = 0.0,
    const double max_angle = 360.0) {
    const double range = max_angle - min_angle;
    const double incr = range / double(nb);
    for (double angle = min_angle; angle < max_angle; angle += incr) {
      double angle_rad = (angle / 360.0) * 6.2831853;
      double x = c.x() + radius * cos(angle_rad);
      double y = c.y() + radius * sin(angle_rad);
      Point point(x, y);
      add_sample(point);
    }
  }

  void append_widely_variable_sampling(const float d1, const float d2) {
    double angle;
    double incr = d1;
    Point c = Point(0.5, 0.5);
    const double radius = 0.5;
    // 0-90 deg -> d1
    for (angle = 0.0; angle < 90.0; angle += incr) {
      double angle_rad = (angle / 360.0) * 6.2831853;
      double x = c.x() + radius * cos(angle_rad);
      double y = c.y() + radius * sin(angle_rad);
      Point point(x, y);
      add_sample(point);
    }
    // 90-180 deg -> d1 -> d2
    for (angle = 90.0; angle < 180.0; angle += incr) {
      double angle_rad = (angle / 360.0) * 6.2831853;
      double x = c.x() + radius * cos(angle_rad);
      double y = c.y() + radius * sin(angle_rad);
      Point point(x, y);
      add_sample(point);
      incr = d1 + (d2 - d1) / 90.0 * (angle - 90);
    }
    // 180-270 deg -> d2
    incr = d2;
    for (angle = 180.0; angle < 270.0; angle += incr) {
      double angle_rad = (angle / 360.0) * 6.2831853;
      double x = c.x() + radius * cos(angle_rad);
      double y = c.y() + radius * sin(angle_rad);
      Point point(x, y);
      add_sample(point);
    }
    // 270-360 deg -> d2 -> d1
    incr = d2;
    for (angle = 270.0; angle < 360.0; angle += incr) {
      double angle_rad = (angle / 360.0) * 6.2831853;
      double x = c.x() + radius * cos(angle_rad);
      double y = c.y() + radius * sin(angle_rad);
      Point point(x, y);
      add_sample(point);
      incr = d2 + (d1 - d2) / 90.0 * (angle - 270.0);
    }
  }

  void append_predefined_line(const int density) {
    std::cerr << "append line...";
    Point start(0.0, 0.5);
    Point end(1.0, 0.5);
    make_line(density, start, end);
    std::cerr << "done" << std::endl;
  }

  void append_predefined_parallel_lines(const int nb_lines, const float space,
    const int density) {
    std::cerr << "append parallel lines...";
    FT x[4];
    x[0] = 0.0;
    x[1] = 0.75;
    x[2] = 1.0;
    x[3] = 1.75;
    FT y = 0.0;
    for (int i = 0; i < nb_lines; ++i) {
      int j = i % 2;
      Point start(x[j], y);
      Point end(x[j + 2], y);
      make_line(density, start, end);
      y += space;
    }
    std::cerr << "done" << std::endl;
  }

  void append_predefined_circle(const int density, const float x,
    const float y, const float radius) {
    std::cerr << "append circle...";
    Point center(x, y);
    make_circle_arc(density, center, radius);
    std::cerr << "done" << std::endl;
  }

  void append_predefined_spiral(const int nb_loops, const int density) {
    std::cerr << "append spiral...";
    Point center(0.5, 0.5);
    const FT max_radius = 0.5;
    const FT spacing = 10. / density; // target spacing
    double radius = max_radius;
    const double max_angle = nb_loops * 360.0;
    for (double angle = max_angle; angle > 0.0; /**/) {
      double angle_rad = (angle / 360.0) * 6.2831853;
      double x = center.x() + radius * cos(angle_rad);
      double y = center.y() + radius * sin(angle_rad);
      Point point(x, y);
      add_sample(point);

      const double angle_incr = atan(spacing / radius);
      angle -= angle_incr;
      radius = max_radius * angle / max_angle;
    }
    std::cerr << "done" << std::endl;
  }

  void append_predefined_half_circle(const int density) {
    std::cerr << "append half circle...";
    Point center(0.5, 0.5);
    make_circle_arc(density, center, 0.5, 0.0, 180.0);
    std::cerr << "done" << std::endl;
  }

  void append_predefined_box(const int density, const float x, const float y,
    const float size_x, const float size_y) {
    std::cerr << "append box...";
    Point a(x - size_x / 2, y - size_y / 2);
    Point b(x + size_x / 2, y - size_y / 2);
    Point c(x + size_x / 2, y + size_y / 2);
    Point d(x - size_x / 2, y + size_y / 2);
    make_line(density, a, b);
    make_line(density, b, c);
    make_line(density, c, d);
    make_line(density, d, a);
    std::cerr << "done" << std::endl;
  }

  void append_predefined_box_with_boundaries(const int density) {
    std::cerr << "append box with boundaries...";

    Point a(0.1, 0.1);
    Point b(0.4, 0.1);
    Point c(0.6, 0.1);
    Point d(0.9, 0.1);
    Point e(0.9, 0.4);
    Point f(0.9, 0.6);
    Point g(0.9, 0.9);
    Point h(0.6, 0.9);
    Point i(0.4, 0.9);
    Point j(0.1, 0.9);
    Point k(0.1, 0.6);
    Point l(0.1, 0.4);

    make_line(density, a, b);
    make_line(density, c, d);
    make_line(density, d, e);
    make_line(density, f, g);
    make_line(density, g, h);
    make_line(density, i, j);
    make_line(density, j, k);
    make_line(density, l, a);

    std::cerr << "done" << std::endl;
  }

  void append_predefined_box_with_missing_corners(const int density) {
    std::cerr << "append box with missing corners...";
    Point a(0.12, 0.1);
    Point b(0.88, 0.1);
    Point c(0.9, 0.12);
    Point d(0.9, 0.88);
    Point e(0.88, 0.9);
    Point f(0.12, 0.9);
    Point g(0.1, 0.88);
    Point h(0.1, 0.12);
    make_line(density, a, b);
    make_line(density, c, d);
    make_line(density, e, f);
    make_line(density, g, h);
    std::cerr << "done" << std::endl;
  }

  void append_predefined_boxes(const int density) {
    std::cerr << "append two boxes...";
    Point a(0.0, 0.0);
    Point b(0.2, 0.0);
    Point c(0.2, 1.0);
    Point d(0.0, 1.0);
    make_line(2 * density, a, b);
    make_line(10 * density, b, c);
    make_line(2 * density, c, d);
    make_line(10 * density, d, a);

    Point e(0.3, 0.0);
    Point f(0.4, 0.0);
    Point g(0.4, 0.3);
    Point h(0.3, 0.3);
    make_line(1 * density, e, f);
    make_line(3 * density, f, g);
    make_line(1 * density, g, h);
    make_line(3 * density, h, e);
    std::cerr << "done" << std::endl;
  }

  void append_predefined_stair(const int density) {
    std::cerr << "append stair...";
    Point a(0.0, 0.0);
    Point b(0.1, 0.0);
    Point c(0.1, 0.1);
    Point d(0.2, 0.1);
    Point e(0.2, 0.2);
    Point f(0.3, 0.2);
    Point g(0.3, 0.3);
    Point h(0.4, 0.3);
    Point i(0.4, 0.4);

    make_line(density, a, b);
    make_line(density, b, c);
    make_line(density, c, d);
    make_line(density, d, e);
    make_line(density, e, f);
    make_line(density, f, g);
    make_line(density, g, h);
    make_line(density, h, i);
    std::cerr << "done" << std::endl;
  }

  void append_predefined_skyline(const int density) {
    std::cerr << "append skyline...";
    Point a(0.0, 0.0);
    Point b(0.1, 0.0);
    Point c(0.1, 0.5);
    Point d(0.3, 0.5);
    Point e(0.3, 0.2);
    Point f(0.4, 0.2);
    Point g(0.4, 0.4);
    Point h(0.6, 0.4);
    Point i(0.6, -0.1);
    Point j(0.7, -0.1);
    Point k(0.7, 0.2);
    Point l(0.8, 0.2);
    Point m(0.8, 0.7);
    Point n(0.9, 0.7);
    Point o(0.9, 0.5);
    Point p(1.0, 0.5);
    Point q(1.0, 0.1);
    Point r(1.1, 0.1);
    Point s(1.1, -0.1);
    Point t(1.2, -0.1);

    make_line(1 * density, a, b);
    make_line(5 * density, b, c);
    make_line(2 * density, c, d);
    make_line(3 * density, d, e);
    make_line(1 * density, e, f);
    make_line(2 * density, f, g);
    make_line(2 * density, g, h);
    make_line(5 * density, h, i);
    make_line(1 * density, i, j);
    make_line(3 * density, j, k);
    make_line(1 * density, k, l);
    make_line(5 * density, l, m);
    make_line(1 * density, m, n);
    make_line(2 * density, n, o);
    make_line(1 * density, o, p);
    make_line(4 * density, p, q);
    make_line(1 * density, q, r);
    make_line(2 * density, r, s);
    make_line(1 * density, s, t);
    std::cerr << "done" << std::endl;
  }

  void append_star(const int nb_branches, const int density) {
    std::cerr << "append star...";
    const double deg_in_rad = CGAL_PI / 180.0;
    const double incr = 180.0 / nb_branches;
    double angle = 0.0;
    const Point center(0.5, 0.5);
    for (int i = 0; i < nb_branches; i++) {
      const double angle_rad = angle * deg_in_rad;
      Vector v(sin(angle_rad), cos(angle_rad));
      Point a = center + v;
      Point b = center - v;
      make_line(density, a, b);
      angle += incr;
    }
    std::cerr << "done" << std::endl;
  }

  void append_predefined_increasingly_sharp_angles(const int density,
    const double min_angle) {
    const double deg_in_rad = CGAL_PI / 180.0;
    double prev_angle = 0.0;
    double curr_angle = min_angle;
    double incr = min_angle;
    const double r1 = 0.5;
    const Point center(0.5, 0.5);
    while (curr_angle < 360.0) {
      Vector va(r1 * cos(prev_angle * deg_in_rad),
        r1 * sin(prev_angle * deg_in_rad));
      Vector vb(r1 * cos(curr_angle * deg_in_rad),
        r1 * sin(curr_angle * deg_in_rad));
      const double average_angle = 0.5 * (prev_angle + curr_angle);
      Vector vc(r1 * cos(average_angle * deg_in_rad),
        r1 * sin(average_angle * deg_in_rad));
      Point a = center + va;
      Point b = center + vb;
      Point c = center + 2 * vc;

      make_line(density, a, c);
      make_line(density, b, c);

      prev_angle = curr_angle;
      curr_angle += incr;
      incr += 2.0;
    }
    noise(1e-5);
  }
};

#endif // SCENE_H_
