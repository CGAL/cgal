

#include "../../../../../Periodic_4_hyperbolic_triangulation_2/include/CGAL/Square_root_2_field.h"
#include "../../../../../Periodic_4_hyperbolic_triangulation_2/include/CGAL/Hyperbolic_octagon_translation_matrix.h"



template<class Kernel>
typename Kernel::Point_2
translate_global(typename Kernel::Circle_2 c, typename Kernel::Point_2 p) {
	double factor = sqrt(to_double(c.squared_radius()));
	double px = to_double( factor*p.x() + c.center().x() );
	double py = to_double( factor*p.y() + c.center().y() );
	return typename Kernel::Point_2(px, py);
}


template<class Kernel>
typename Kernel::Point_2
translate_unity(typename Kernel::Circle_2 c, typename Kernel::Point_2 p) {
	double factor = sqrt(to_double(c.squared_radius()));
	double px = to_double( (p.x() - c.center().x()) / factor );
	double py = to_double( (p.y() - c.center().y()) / factor );
	return typename Kernel::Point_2(px, py);
}


template<class Kernel>
std::list<typename Kernel::Point_2> 
compute_images(typename Kernel::Circle_2 circ, std::list<typename Kernel::Point_2> pts) {
	typedef Hyperbolic_octagon_translation_matrix<Kernel> Hyperbolic_matrix;
	std::vector<Hyperbolic_matrix> gens;
	get_generators(gens);
	std::list<typename Kernel::Point_2> images;

	for (typename std::list<typename Kernel::Point_2>::iterator it = pts.begin(); it != pts.end(); ++it) {
		for (int j = 0; j < gens.size(); j++) { 
			images.push_back( translate_global<Kernel>( circ, gens[j].apply( translate_unity<Kernel>(circ, *it) ) ) );
		}
	}

	return images;
}


template<class Kernel>
std::list<typename Kernel::Point_2> 
compute_images_one_direction(typename Kernel::Circle_2 circ, std::list<typename Kernel::Point_2> pts, int dir) {
	typedef Hyperbolic_octagon_translation_matrix<Kernel> Hyperbolic_matrix;
	std::vector<Hyperbolic_matrix> gens;
	get_generators(gens);
	std::list<typename Kernel::Point_2> images;

	for (typename std::list<typename Kernel::Point_2>::iterator it = pts.begin(); it != pts.end(); ++it) {
		images.push_back( translate_global<Kernel>( circ, gens[dir].apply( translate_unity<Kernel>(circ, *it) ) ) );
	}

	return images;
}

template<class Kernel>
std::list<typename Kernel::Point_2> 
compute_images_plane_model(typename Kernel::Segment_2 seg) {
	
	double sx = CGAL::to_double(seg.source().x());
	double sy = CGAL::to_double(seg.source().y());
	double tx = CGAL::to_double(seg.target().x());
	double ty = CGAL::to_double(seg.target().y());

	double mx = (sx + tx)/2.;
	double my = (sy + ty)/2.;

	// Consider the whole thing to be 10 units long, so we need to find a unit length
	double seg_length = sqrt(CGAL::to_double(seg.squared_length()));
	double unit = seg_length / 10.;

	typedef typename Kernel::Point_2 Point;
	std::list<Point> images;
	images.push_back( Point( mx, unit + my ) );
	images.push_back( Point( mx, my + unit*( 1.+sqrt(2.)-sqrt(2.*(sqrt(2.)+1.)) ) ) );
	images.push_back( Point( mx, my + unit*( 1.+sqrt(2.)+sqrt(2.*(sqrt(2.)+1.)) ) ) );

	return images;
}
