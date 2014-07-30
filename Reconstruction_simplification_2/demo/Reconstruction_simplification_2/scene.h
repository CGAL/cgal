#ifndef SCENE_H_
#define SCENE_H_

// STL
#include <fstream>

//Qt
#include <QtOpenGL>

// local
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "Reconstruction_simplification_kerneled_2.h"
#include "../../include/CGAL/Reconstruction_simplification_2.h"
#include "../../include/CGAL/List_output.h"


#include "third/CImg.h"
#include "random.h"
#include <utility>      // std::pair
#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>

class Scene {

public:
	typedef std::pair<Point, FT> PointMassPair;

	typedef std::list<PointMassPair> PointMassList;
	typedef PointMassList::const_iterator InputIterator;

	typedef CGAL::value_type_traits<InputIterator>::type MassPoint;

	typedef CGAL::First_of_pair_property_map<PointMassPair> PointPMap;
	typedef CGAL::Second_of_pair_property_map<PointMassPair> MassPMap;

	typedef CGAL::Reconstruction_simplification_2<K, PointPMap,
			MassPMap> R_s_2;

	typedef K::Segment_2 Segment;

	typedef R_s_2::FT FT;
	typedef R_s_2::Point Point;
	typedef R_s_2::Vector Vector;

	typedef R_s_2::Vertex Vertex;
	typedef R_s_2::Vertex_handle Vertex_handle;
	typedef R_s_2::Vertex_iterator Vertex_iterator;
	typedef R_s_2::Vertex_circulator Vertex_circulator;
	typedef R_s_2::Finite_vertices_iterator Finite_vertices_iterator;

	typedef R_s_2::Edge Edge;
	typedef R_s_2::Edge_iterator Edge_iterator;
	typedef R_s_2::Edge_circulator Edge_circulator;
	typedef R_s_2::Finite_edges_iterator Finite_edges_iterator;

	typedef R_s_2::Face Face;
	typedef R_s_2::Face_handle Face_handle;
	typedef R_s_2::Face_iterator Face_iterator;
	typedef R_s_2::Face_circulator Face_circulator;
	typedef R_s_2::Finite_faces_iterator Finite_faces_iterator;

	typedef R_s_2::Vertex_handle_map Vertex_handle_map;
	typedef R_s_2::Face_handle_map Face_handle_map;

	typedef R_s_2::Vertex_handle_set Vertex_handle_set;
	typedef R_s_2::Edge_set Edge_set;

	typedef R_s_2::Edge_list Edge_list;

	typedef R_s_2::Cost Cost;
	typedef R_s_2::Sample Sample;
	typedef R_s_2::Sample_list Sample_list;
	typedef R_s_2::Sample_list_const_iterator Sample_list_const_iterator;

	typedef R_s_2::Point_list Point_list;
	typedef R_s_2::Point_list_const_iterator Point_list_const_iterator;

	typedef R_s_2::PSample PSample;
	typedef R_s_2::SQueue SQueue;

	typedef R_s_2::Reconstruction_edge_2 PEdge;


private:
	// data
	std::list<Sample> m_samples;
	double m_min_mass;
	double m_max_mass;

	Reconstruction_simplification_kerneled_2* m_pwsrec;
	int m_ignore;

	// bbox
	double m_bbox_x;
	double m_bbox_y;
	double m_bbox_size;

public:
	Scene() {
		srand(0); // for sake of repeatability
		m_min_mass = 0.0;
		m_max_mass = 0.0;
		m_ignore = 0;
		m_bbox_x = 0.0;
		m_bbox_y = 0.0;
		m_bbox_size = 1.0;

		m_pwsrec = new Reconstruction_simplification_kerneled_2();
	}

	~Scene() {
		clear();
	}

	void clear() {
		m_pwsrec->clear();
		m_samples.clear();
		m_max_mass = 0.0;
	}

	void set_min_mass(const double min_value) {
		m_min_mass = min_value;
	}

	void set_nb_edges_to_ignore(int nb) {
		m_ignore = nb;
	}

	void invert_mass() {
		double new_max_mass = 0.0;
		std::list<Sample>::iterator it;
		for (it = m_samples.begin(); it != m_samples.end(); ++it) {
			Sample& sample = *it;
			sample.mass() = m_max_mass - sample.mass();
			new_max_mass = (std::max)(new_max_mass, sample.mass());
		}
		m_max_mass = new_max_mass;
	}

	void clamp_mass() {
		std::list<Sample> new_samples;
		std::list<Sample>::iterator it;
		for (it = m_samples.begin(); it != m_samples.end(); ++it) {
			Sample& sample = *it;
			if (sample.mass() > m_min_mass) {
				sample.mass() = 1.0;
				new_samples.push_back(sample);
			}
		}
		m_samples = new_samples;
		m_max_mass = 1.0;
	}

	void subdivide() {
		if (m_samples.size() < 3)
			return;

		std::list<Sample> new_samples;
		std::list<Sample>::const_iterator it = m_samples.begin();
		std::list<Sample>::const_iterator last = it++;
		while (it != m_samples.end()) {
			Point p = CGAL::midpoint(last->point(), it->point());
			FT m = 0.5 * (last->mass() + it->mass());
			new_samples.push_back(Sample(p, m));
			last = it++;
		}
		it = m_samples.begin();
		Point p = CGAL::midpoint(last->point(), it->point());
		FT m = 0.5 * (last->mass() + it->mass());
		new_samples.push_back(Sample(p, m));

		std::list<Sample> final_samples;
		std::list<Sample>::const_iterator it2 = new_samples.begin();
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
		m_samples.push_back(Sample(point, mass));
		m_max_mass = (std::max)(m_max_mass, mass);
	}

	void add_outliers(const unsigned int nb) {
		std::cerr << "add " << nb << " outliers...";
		for (unsigned int i = 0; i < nb; i++) {
			Point outlier = CGAL::ORIGIN + random_vec<Vector>(1.3);
			m_samples.push_back(outlier);
		}
		std::cerr << "done" << std::endl;
	}

	void noise(const FT scale) {
		std::cerr << "noise by " << scale << "...";
		std::list<Sample>::iterator it;
		for (it = m_samples.begin(); it != m_samples.end(); it++) {
			Sample& sample = *it;
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
		std::list<Sample>::iterator it;
		for (it = m_samples.begin(); it != m_samples.end(); ++it) {
			Sample& sample = *it;
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
		std::list<Sample>::const_iterator it = m_samples.begin();
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

	void load(const QString& filename, const bool gradient) {
		// TODO: load xml

		if (filename.contains(".xy", Qt::CaseInsensitive)) {
			load_xy_file(filename);
			normalize_points();
			return;
		}

		if (filename.contains(".txt", Qt::CaseInsensitive)) {
			load_xy_file(filename);
			normalize_points();
			return;
		}

		if (filename.contains(".gtn", Qt::CaseInsensitive)) {
			load_gathan_file(filename);
			normalize_points();
			return;
		}

		if (filename.contains(".dat", Qt::CaseInsensitive)) {
			load_dat_file(filename);
			normalize_points();
			return;
		}

		if (filename.contains(".bmp", Qt::CaseInsensitive)) {
			if (gradient)
				load_gradient(filename);
			else
				load_image(filename);
			normalize_points();
			return;
		}

		std::cerr << "Invalid file (try .xy, .dat, .bmp)" << std::endl;
	}

	void load_xy_file(const QString& fileName) {

		std::cout << "FILENAME " << fileName.toUtf8().constData() << std::endl;
		std::ifstream ifs(qPrintable(fileName));
		std::cerr << "read xy...";
		Point point;
		unsigned int nb = 0;
		while (ifs >> point) {
			add_sample(point, 1.0);
			nb++;
		}
		std::cerr << "done (" << nb << " points)" << std::endl;
		ifs.close();
	}

	void load_gathan_file(const QString& fileName) {
		std::ifstream ifs(qPrintable(fileName));
		std::cerr << "read gathan...";
		Point point;
		unsigned nb, x, y, z;
		ifs >> nb >> x >> y >> z;
		for (unsigned i = 0; i < nb; ++i) {
			ifs >> x >> point;
			add_sample(point, 1.0);
		}
		std::cerr << "done (" << nb << " points)" << std::endl;
		ifs.close();
	}

	void load_dat_file(const QString& fileName) {
		std::ifstream ifs(qPrintable(fileName));
		std::cerr << "read dat...";
		Point point;
		unsigned int n, m, nb = 0;
		ifs >> n;
		for (unsigned int i = 0; i < n; ++i) {
			ifs >> m;
			for (unsigned int j = 0; j < m; ++j) {
				ifs >> point;
				add_sample(point, 1.0);
			}
			nb += m;
		}
		std::cerr << "done (" << nb << " points)" << std::endl;
		ifs.close();
	}

	void load_image(const QString& fileName) {
		std::cerr << "read image...";
		cimg_library::CImg<float> image(qPrintable(fileName));
		std::cerr << "done" << std::endl;

		std::cerr << "compute grayscale...";
		cimg_library::CImg<float> grayscale =
				image.RGBtoHSV().get_channel(2).normalize(0.0f, 1.0f);
		std::cerr << "done" << std::endl;

		// turn pixels into weighted samples
		std::cerr << "add samples...";
		for (int i = 0; i < grayscale.width(); i++) {
			for (int j = 0; j < grayscale.height(); j++) {
				float mass = 1.0f - grayscale.atXY(i, j);
				double x = double(i) / grayscale.width();
				double y = 1.0 - double(j) / grayscale.height();
				add_sample(Point(x, y), mass);
			}
		}
		std::cerr << "done (" << m_samples.size() << ")" << std::endl;
	}

	void load_gradient(const QString& fileName) {
		std::cerr << "read image...";
		cimg_library::CImg<float> image(qPrintable(fileName));
		std::cerr << "done" << std::endl;

		std::cerr << "compute gradient...";
		cimg_library::CImgList<float> grad = image.get_gradient();
		cimg_library::CImg<float> normgrad = sqrt(
				grad[0].pow(2) + grad[1].pow(2)).normalize(0.0f, 1.0f);
		std::cerr << "done" << std::endl;

		// turn pixels into weighted samples
		std::cerr << "add samples...";
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

	void print_vertex(Vertex vertex) {
		std::cout <<"vertex " <<  vertex << std::endl;
	}


	void print_edge(PEdge edge) {
		int i = ((edge).edge()).second;
		Point a = ((edge).edge()).first->vertex((i+1)%3)->point();
		Point b = ((edge).edge()).first->vertex((i+2)%3)->point();
		std::cout <<"( " << (edge).priority()  <<  ") ( " << a
								<< " , " << b << " )" << std::endl;
	}


	void debug_print() {


	    std::vector<Point> isolated_points;
		std::vector<Segment> edges;

		typedef std::back_insert_iterator<std::vector<Point> > Point_it;
		typedef std::back_insert_iterator<std::vector<Segment> >  Edge_it;

		Point_it point_it(isolated_points);
		Edge_it  edge_it(edges);

		CGAL::List_output<K, Point_it, Edge_it> list_output(point_it, edge_it);

		m_pwsrec->extract_solid_elements(list_output);

		int vertex_count = 0;
		for (std::vector<Point>::iterator it = isolated_points.begin();
				it != isolated_points.end(); it++) {
			vertex_count++;
			std::cout  <<  *it << std::endl;
		}
		assert(vertex_count == 18);

		int edge_count = 0;
		for (std::vector<Segment>::iterator it = edges.begin();
				it != edges.end(); it++) {
			std::cout << *it << std::endl;
			edge_count++;
		}
	}

	void save(const QString& filename) {
		std::cout << "SAVE-------------" << std::endl;
		debug_print();

		if (filename.contains(".edges", Qt::CaseInsensitive)) {
			std::ofstream ofs(qPrintable(filename));
			m_pwsrec->save_edges(ofs, m_ignore);
			ofs.close();
			return;
		}

		Sample_list samples;
		for (std::list<Sample>::iterator it = m_samples.begin();
				it != m_samples.end(); ++it) {
			Sample& s = *it;
			if (s.mass() < m_min_mass)
				continue;
			samples.push_back(&s);
		}

		if (filename.contains(".xy", Qt::CaseInsensitive)) {
			save_xy(filename, samples);
			return;
		}

		if (filename.contains(".poff", Qt::CaseInsensitive)) {
			save_poff(filename, samples);
			return;
		}

		if (filename.contains(".gtn", Qt::CaseInsensitive)) {
			save_gtn(filename, samples);
			return;
		}

		if (filename.contains(".pwn", Qt::CaseInsensitive)) {
			save_pwn(filename, samples);
			return;
		}
	}

	void save_pwn(const QString& filename, const Sample_list& samples) {
		std::list<Vector> normals;
		compute_normals(samples, normals);
		std::ofstream ofs(qPrintable(filename));
		std::list<Vector>::const_iterator ni = normals.begin();
		for (Sample_list_const_iterator it = samples.begin();
				it != samples.end(); ++it) {
			Sample* sample = *it;
			ofs << sample->point() << " " << *ni << std::endl;
			ni++;
		}
		ofs.close();
	}

	void compute_normals(const Sample_list& samples,
			std::list<Vector>& normals) {
		normals.clear();
		Point last = samples.back()->point();
		Sample_list_const_iterator si = samples.begin();
		while (si != samples.end()) {
			Point p = (*si)->point();
			si++;
			Point next = samples.front()->point();
			if (si != samples.end())
				next = (*si)->point();

			Vector ab = p - last;
			Vector bc = next - p;
			Vector ab90(ab.y(), -ab.x());
			Vector bc90(bc.y(), -bc.x());

			Vector ni = ab90 + bc90;
			FT norm = std::sqrt(ni * ni);
			if (norm != 0.0)
				ni = ni / norm;
			normals.push_back(ni);

			last = p;
		}
	}

	void save_xy(const QString& filename, const Sample_list& samples) {
		std::ofstream ofs(qPrintable(filename));
		for (Sample_list_const_iterator it = samples.begin();
				it != samples.end(); ++it) {
			Sample* sample = *it;
			ofs << sample->point() << std::endl;
		}
		ofs.close();
	}

	void save_poff(const QString& filename, const Sample_list& samples) {
		std::ofstream ofs(qPrintable(filename));
		ofs << "POFF " << samples.size() << " 0 0" << std::endl;
		for (Sample_list_const_iterator it = samples.begin();
				it != samples.end(); ++it) {
			Sample* sample = *it;
			ofs << sample->point() << std::endl;
		}
		ofs.close();
	}

	void save_gtn(const QString& filename, const Sample_list& samples) {
		std::ofstream ofs(qPrintable(filename));
		ofs << samples.size() << " 2 0 0" << std::endl;
		unsigned i = 0;
		for (Sample_list_const_iterator it = samples.begin();
				it != samples.end(); ++it, ++i) {
			Sample* sample = *it;
			ofs << i << " " << sample->point() << std::endl;
		}
		ofs.close();
	}

	// RECONSTRUCTION //

	void set_parameters(const int verbose, const int mchoice,
			const bool use_flip, const double alpha, const double norm_tol,
			const double tang_tol, const unsigned relocation,
			const double ghost) {

		m_pwsrec->set_verbose(verbose);
		m_pwsrec->set_mchoice(mchoice);
		m_pwsrec->set_use_flip(use_flip);
		m_pwsrec->set_alpha(alpha);
		m_pwsrec->set_tang_tol(tang_tol * m_bbox_size);
		m_pwsrec->set_norm_tol(norm_tol * m_bbox_size);
		m_pwsrec->set_relocation(relocation);
		m_pwsrec->set_ghost(ghost);
	}

	bool init_reconstruction(const double percentage) {
		std::cout << " init_reconstruction " << std::endl;

		if (m_samples.empty()) {
			std::cerr << "initialization failed (empty point set)" << std::endl;
			return false;
		}

		Sample_list vertices, samples;
		select_samples(percentage, vertices, samples);

		PointMassList point_mass_list;
		Sample_list_const_iterator it;
		for (it = vertices.begin(); it != vertices.end(); it++) {
			std::cout << "Sample in Scene " << (*it)->point() << " : "
					<< (*it)->mass() << std::endl;
			point_mass_list.push_back(
					std::make_pair((*it)->point(), (*it)->mass()));
		}

		PointPMap point_pmap;
		MassPMap mass_pmap;
		MassPoint mp;

		m_pwsrec->initialize(point_mass_list.begin(), point_mass_list.end(),
				point_pmap, mass_pmap);

		return true;
	}

	void decimate(const double percentage) {
		std::cout << "Decimate from " << m_samples.size() << " to...";
		std::list<Sample> selected;

		std::list<Sample>::iterator it;
		for (it = m_samples.begin(); it != m_samples.end(); it++) {
			const double rd = random_double(0.0, 1.0);
			if (rd >= percentage)
				selected.push_back(*it);
		}

		m_samples.clear();
		std::copy(selected.begin(), selected.end(),
				std::back_inserter(m_samples));
		std::cout << m_samples.size() << std::endl;
	}

	void keep_one_point_out_of(const int n) {
		std::cout << "Decimate from " << m_samples.size() << " to...";
		std::list<Sample> selected;

		int index = 0;
		std::list<Sample>::iterator it;
		for (it = m_samples.begin(); it != m_samples.end(); it++, index++) {
			if (index % n == 0)
				selected.push_back(*it);
		}

		m_samples.clear();
		std::copy(selected.begin(), selected.end(),
				std::back_inserter(m_samples));
		std::cout << m_samples.size() << std::endl;
	}

	void select_samples(const double percentage, Sample_list& vertices,
			Sample_list& samples) {
		std::list<Sample>::iterator it;
		for (it = m_samples.begin(); it != m_samples.end(); ++it) {
			Sample& s = *it;
			if (s.mass() <= m_min_mass)
				continue;

			samples.push_back(&s);
			FT rv = random_double(0.0, 1.0);
			if (rv <= percentage)
				vertices.push_back(&s);
		}
	}

	void reconstruct_until(const unsigned nv) {
		std::cout << "reconstruct_until" << std::endl;
		m_pwsrec->reconstruct_until(nv);
	}

	void reconstruct(const unsigned steps) {
		std::cout << "reconstruct" << std::endl;
		m_pwsrec->reconstruct(steps);
	}

	void relocate_all_vertices() {
		std::cout << "relocate_all_vertices" << std::endl;
		m_pwsrec->relocate_all_vertices();
	}

	void print_stats() {
		std::cout << "print_stats" << std::endl;
		m_pwsrec->print_stats();
	}

	// RENDER //

	void render(const bool view_points, const bool view_vertices,
			const bool view_edges, const bool view_ghost_edges,
			const bool view_edge_cost, const bool view_edge_priority,
			const bool view_bins, const bool view_foot_points,
			const bool view_relocation, const bool view_tolerance,
			const bool view_incolors, const bool view_edge_relevance,
			const float point_size, const float vertex_size,
			const float line_thickness) {
		if (m_pwsrec == NULL) {
			return;
		}

		if (view_edges)
			m_pwsrec->draw_edges(0.5f * line_thickness, 0.4f, 0.4f, 0.4f);

		if (view_edge_cost)
			m_pwsrec->draw_costs(line_thickness, view_ghost_edges);

		if (view_edge_priority)
			m_pwsrec->draw_pedges(line_thickness);

		if (view_edge_relevance)
			m_pwsrec->draw_relevance(line_thickness, m_ignore, view_incolors);

		if (view_relocation)
			m_pwsrec->draw_relocation();

		if (view_vertices)
			m_pwsrec->draw_vertices(vertex_size, 0.0f, 0.0f, 0.5f);

		if (view_bins)
			m_pwsrec->draw_bins(0.5f * line_thickness);

		if (view_foot_points)
			m_pwsrec->draw_footpoints(line_thickness, 0.2f, 0.8f, 0.2f);

		if (view_tolerance)
			draw_circles();

		if (view_points)
			draw_samples(point_size);
	}

	void render_simulation(const Point& point, int option,
			const float vertex_size, const float line_width) {

		std::cout << "render_simulation" << std::endl;

		switch (option) {
		case 0:
			m_pwsrec->draw_one_ring(vertex_size, line_width, point);
			break;
		case 1:
			m_pwsrec->draw_blocking_edges(vertex_size, line_width, point);
			break;
		case 2:
			m_pwsrec->draw_collapsible_edge(vertex_size, line_width, point);
			break;
		case 3:
			m_pwsrec->draw_simulation(vertex_size, line_width, point);
			break;
		case 4:
			m_pwsrec->draw_remove_queue_stencil(vertex_size, line_width, point);
			break;
		case 5:
			m_pwsrec->draw_cost_stencil(vertex_size, line_width, point);
			break;
		case 6:
			m_pwsrec->draw_push_queue_stencil(vertex_size, line_width, point);
			break;
		default:
			break;
		}
	}

	void draw_samples(const float point_size) {
		double max_value = m_max_mass;
		if (max_value == 0.0)
			max_value = 1.0;

		::glPointSize(point_size);
		::glBegin(GL_POINTS);
		std::list<Sample>::const_iterator it;
		for (it = m_samples.begin(); it != m_samples.end(); it++) {
			double mass = it->mass();
			if (mass <= m_min_mass)
				continue;

			float value = mass / m_max_mass;
			float grey = 0.9 * (1.0f - value);
			::glColor3f(grey, grey, grey);
			const Point& p = it->point();
			::glVertex2d(p.x(), p.y());
		}
		::glEnd();
	}

	void draw_circles() {
		Sample_list vertices, samples;
		select_samples(1.0, vertices, samples);
		double percentage = 500.0 / double(vertices.size());
		percentage = (std::min)(percentage, 1.0);

		samples.clear();
		vertices.clear();
		select_samples(percentage, vertices, samples);

		::glEnable(GL_BLEND);
		::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		::glColor4f(1.0f, 1.0f, 0.2f, 0.25f);
		for (Sample_list_const_iterator it = samples.begin();
				it != samples.end(); it++) {
			Sample* sample = *it;
			draw_one_circle(sample->point());
		}
		::glDisable(GL_BLEND);
	}

	void draw_one_circle(const Point& center) {
		unsigned N = 10;
		const double r = std::sqrt(m_pwsrec->get_norm_tol());
		const double x = center.x();
		const double y = center.y();
		::glBegin (GL_POLYGON);
		for (unsigned i = 0; i < N; ++i) {
			double angle = 2.0 * M_PI * (double(i) / double(N));
			double u = r * std::cos(angle) + x;
			double v = r * std::sin(angle) + y;
			::glVertex2d(u, v);
		}
		::glEnd();
	}

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
		const double deg_in_rad = 3.1415926535897932384626 / 180.0;
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
		const double deg_in_rad = 3.1415926535897932384626 / 180.0;
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
