#ifndef CGAL_DELAUNAY_TRIANGULATION_ON_HYPERBOLIC_SURFACE_2_H
#define CGAL_DELAUNAY_TRIANGULATION_ON_HYPERBOLIC_SURFACE_2_H

#include <CGAL/Triangulation_on_hyperbolic_surface_2.h>
#include <boost/numeric/interval.hpp>


namespace CGAL{

template<class Traits>
struct Delaunay_triangulation_attributes {
	template<class CMap>
	struct Dart_wrapper{
		typedef Cell_attribute<CMap, Complex_number<typename Traits::FT>> Edge_attrib;
		typedef Cell_attribute<CMap, typename Triangulation_on_hyperbolic_surface_2<Traits,Delaunay_triangulation_attributes<Traits>>::Anchor> Face_attrib;
		typedef std::tuple<void,Edge_attrib,Face_attrib> Attributes;
	};
};

enum Locate_type {
	VERTEX,
	EDGE,
	FACE,
	OUTSIDE = 4
};

template<class Traits>
class Delaunay_triangulation_on_hyperbolic_surface_2: public Triangulation_on_hyperbolic_surface_2<Traits, Delaunay_triangulation_attributes<Traits>> {
public:
	typedef Triangulation_on_hyperbolic_surface_2<Traits,Delaunay_triangulation_attributes<Traits>>                     Base;
	typedef Combinatorial_map<2,Delaunay_triangulation_attributes<Traits>>                                              CMap;
	typedef typename Triangulation_on_hyperbolic_surface_2<Traits,Delaunay_triangulation_attributes<Traits>>::Anchor    Anchor;

	typedef typename Traits::FT 										Number;
	typedef typename Traits::Complex                                    Complex;
	typedef typename Traits::Hyperbolic_point_2                         Point;
	typedef typename Traits::Hyperbolic_Voronoi_point_2                 Voronoi_point;
	typedef typename CMap::Dart_descriptor                              Dart_descriptor;
	typedef typename CMap::Dart_range                                   Dart_range;
	typedef typename CMap::template One_dart_per_cell_range<0>          Vertex_range;
	typedef typename CMap::template One_dart_per_cell_range<1>          Edge_range;
	typedef typename CMap::template One_dart_per_cell_range<2>          Face_range;
	typedef typename CMap::Dart_const_descriptor                        Dart_const_descriptor;
	typedef typename CMap::Dart_const_range                             Dart_const_range;
	typedef typename CMap::template One_dart_per_cell_const_range<1>    Edge_const_range;
	typedef typename CMap::template One_dart_per_cell_const_range<2>    Face_const_range;
	typedef Hyperbolic_isometry_2<Traits>                               Isometry;
	typedef Hyperbolic_fundamental_domain_2<Traits>                     Domain;
	typedef boost::numeric::interval<double>                            Interval;

	unsigned const NB_SIDES = 3;
	unsigned const NULL_INDEX = 4;

	//---------- CONSTRUCTORS
	Delaunay_triangulation_on_hyperbolic_surface_2() {};
	Delaunay_triangulation_on_hyperbolic_surface_2(CMap & cmap, Anchor & anch);
	Delaunay_triangulation_on_hyperbolic_surface_2(Domain const & domain);
	Delaunay_triangulation_on_hyperbolic_surface_2(Base & triangulation);

	//---------- UTILITIES
	Anchor & anchor(Dart_descriptor const dart);
	Anchor const & anchor(Dart_const_descriptor const dart) const;
	Anchor & anchor();
	Anchor const & anchor() const;
	unsigned index_in_anchor(Dart_const_descriptor const dart) const;
	Dart_descriptor ith_dart(unsigned i, Anchor const & anch);
	void display_vertices(Anchor const & anch, bool round = true) const;
	bool are_triangles_equal(Anchor const & anchor1, Anchor const & anchor2);  // supprimer ?
	bool is_valid() const;
	void to_stream(std::ostream & s) const;
 	void from_stream(std::istream & s);

	//---------- location and insertion
	std::tuple<Locate_type, unsigned> relative_position(Point const & query, Anchor const & anch) const;
	Orientation hyperbolic_orientation_2(Point const & p, Point const & q, Point const & r) const;  // à mettre dans le code de hyperbolic_traits
	std::tuple<Anchor, unsigned> locate_visibility_walk(Point const & query, Anchor const & anch);  // return réf ? + const + faire 1 locate avec un param bool qui change l'algo utilisé
	std::tuple<Anchor, unsigned> locate_visibility_walk(Point const & query);
	std::tuple<Anchor, unsigned> locate_straight_walk(Point const & query, Anchor const & anch);
/*       
	SYNTAXE locate de CGAL (cf doc) :
	Face_descriptor CGAL::Triangulation_2< Traits, Tds >::locate (
			const Point &   query,
			Locate_type &   lt,
			int &   li,
			Face_descriptor     h = Face_descriptor() 
	) const */

	std::vector<Anchor> insert(Point const & query, Anchor & anch); // return réf ?
	std::vector<Anchor> insert(Point const & query);

	//---------- Delaunay related methods
	void flip(Dart_descriptor dart);
	unsigned make_Delaunay();
	std::tuple<unsigned, std::vector<Dart_descriptor>> Delaunay_insert(Point const & query, Anchor & anch); // return réf ?
	std::tuple<unsigned, std::vector<Dart_descriptor>> Delaunay_insert(Point const & query);

	//---------- eps-net methods
	bool epsilon_net(double epsilon);
	bool is_epsilon_covering(const double epsilon);
	bool is_epsilon_packing(const double epsilon);
	bool is_epsilon_net(const double epsilon);
	
	double shortest_loop() const;

private:
	//---------- CONSTRUCTORS
	void set_anchors();

	//---------- UTILITIES
	unsigned ccw(unsigned i) const;
	unsigned cw(unsigned i) const;
	void set_attribute(Dart_descriptor dart, Anchor const & anchor);
	void set_attribute(Dart_descriptor dart, Complex const & cross_ratio);
	Anchor & update_infos(Dart_descriptor dart, Point const & r, Point const & s, Point const & t, Point const & query);
	
	//---------- location and insertion
	std::vector<Anchor> insert_in_face(Point const & query, Anchor& anch);  // return réf ?
	std::vector<Anchor> insert_in_edge(Point const & query, Anchor& anch, unsigned index);

	//---------- Delaunay related methods
	void push_flippable_edge(Dart_descriptor const dart, std::list<Dart_descriptor> & darts_to_flip);
	std::tuple<unsigned, std::vector<Dart_descriptor>> restore_Delaunay(std::list<Dart_descriptor> & darts_to_flip);  // return réf ?

	//---------- eps-net methods
	Number delta(Point const & u, Point const & v) const;
	Point approx_circumcenter(Anchor const & anch) const;
	Number delta_min(Anchor const & anch, Point const & approx_center) const;
	Number delta_max(Anchor const & anch, Point const & approx_center) const;
	void push_triangle(Dart_descriptor const dart, std::list<Dart_descriptor> & triangles, size_t & triangles_list_mark);
};


//---------- CONSTRUCTORS

template<class Traits>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
Delaunay_triangulation_on_hyperbolic_surface_2(CMap & cmap, Anchor & anch)
: Base(cmap, anch)
{       
	set_anchors();
}

template<class Traits>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
Delaunay_triangulation_on_hyperbolic_surface_2(Domain const & domain)
: Base(domain)
{       
	set_anchors();
}

template<class Traits>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
Delaunay_triangulation_on_hyperbolic_surface_2(Base & triangulation)
: Base(triangulation.combinatorial_map(), triangulation.anchor())
{       
	set_anchors();
}

// only used in constructors
template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
set_anchors()
{
	CMap & cmap = this->combinatorial_map_;

	std::queue<Anchor> bfs_queue;
	size_t in_queue = cmap.get_new_mark();
	cmap.unmark_all(in_queue);
	bfs_queue.push(this->anchor_);
	cmap.mark(this->anchor_.dart, in_queue);
	cmap.mark(Base::ccw(this->anchor_.dart), in_queue);
	cmap.mark(Base::cw(this->anchor_.dart), in_queue);

	while (!bfs_queue.empty()) {
		Anchor & current = bfs_queue.front();
		set_attribute(current.dart, current);
		Dart_descriptor invader = current.dart;
		for (int i = 0; i < 3; i++) {
			Dart_descriptor invaded = Base::opposite(invader);
			if (!cmap.is_marked(invaded, in_queue)) {
				Complex cross_ratio = Base::get_cross_ratio(invader);
				Point & c = current.vertices[i % 3];
				Point & a = current.vertices[(i + 1) % 3];
				Point & b = current.vertices[(i + 2) % 3];
				Point d = Base::fourth_point_from_cross_ratio(a, b, c, cross_ratio);
				bfs_queue.push(Anchor(invaded, a, c, d));
				cmap.mark(invaded, in_queue);
				cmap.mark(Base::ccw(invaded), in_queue);
				cmap.mark(Base::cw(invaded), in_queue);
			}
			invader = Base::ccw(invader);
		}
		bfs_queue.pop();
	}
	cmap.free_mark(in_queue);

	this->anchor_ = Anchor();
	this->has_anchor_ = false;
}


//---------- UTILITIES


template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor&
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
anchor(Dart_descriptor const dart)
{
	return this->combinatorial_map_.template info<2>(dart);
}

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor const &
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
anchor(Dart_const_descriptor const dart) const
{
	return this->combinatorial_map_.template info<2>(dart);
}

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor &
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
anchor()
{
	return anchor(this->combinatorial_map_.darts().begin());
}

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor const &
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
anchor() const 
{
	return anchor(this->combinatorial_map_.darts().begin());
}

template<class Traits>
unsigned
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
index_in_anchor(Dart_const_descriptor const dart) const
{       
	unsigned index = NULL_INDEX;
	Anchor const & anch = anchor(dart);
	Dart_const_descriptor current_dart = anch.dart;
	for (unsigned i = 0; i < NB_SIDES; ++i) {
		if (current_dart == dart) {
			index = i;
			break;
		}
		current_dart = Base::const_ccw(current_dart);
	}
	CGAL_assertion(index != NULL_INDEX);
	return index;
}

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Dart_descriptor
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
ith_dart(unsigned i, Anchor const & anch)
{       
	Dart_descriptor dart = anch.dart;
	for (unsigned j = 0; j < i; ++j) {
		dart = Base::ccw(dart);
	}
	return dart;
}

template<class Traits>
unsigned
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
ccw(unsigned i) const
{
	CGAL_precondition( i <= NB_SIDES);
	return (i + 1) % NB_SIDES;
}

template<class Traits>
unsigned
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
cw(unsigned i) const
{
	CGAL_precondition( i <= NB_SIDES);
	return (i + 2) % NB_SIDES;
}

template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
set_attribute(Dart_descriptor dart, Anchor const & anchor)
{
	this->combinatorial_map_.template set_attribute<2>(dart,
		this->combinatorial_map_.template create_attribute<2>(anchor));
}

template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
set_attribute(Dart_descriptor dart, Complex const & cross_ratio)
{
	this->combinatorial_map_.template set_attribute<1>(dart,
		this->combinatorial_map_.template create_attribute<1>(cross_ratio));
}

template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
display_vertices(Anchor const & anch, bool round) const
{       
	if (round) {
		for (unsigned i = 0; i < NB_SIDES; ++i) {
			std::cout << "vertex " << i << " : ";
			std::cout << "(" << to_double(anch.vertices[i].x()) << "," << to_double(anch.vertices[i].y()) <<")" << std::endl;
		}
	} else {
		for (unsigned i = 0; i < NB_SIDES; ++i) {
			std::cout << "vertex " << i << " : ";
			std::cout << "(" << anch.vertices[i].x() << "," << anch.vertices[i].y() <<")" << std::endl;
		}
	}
}

// Output: true iff the triangles described by the anchors are equal (same triangle in the cmap and same vertices in H^2)
template<class Traits>
bool
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
are_triangles_equal(Anchor const & anchor1, Anchor const & anchor2)
{
	Anchor& triangle1 = anchor(anchor1.dart);
	Anchor& triangle2 = anchor(anchor2.dart);
	
	// check if anchors correspond to same face in the cmap
	if (triangle1.dart != triangle2.dart) {
		return false;
	}

	// check vertices
	bool res = true;
	for (unsigned i = 0; i < NB_SIDES; ++i) {
		bool found = false;
		for (unsigned j = 0; j < NB_SIDES; ++j) {
			if (anchor1.vertices[i] ==  anchor2.vertices[j]) {
				found = true;
			}
		}
		res = res && found;
	}
	return res;
}

template<class Traits>
bool
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
is_valid() const
{
	CGAL_assertion(Base::is_valid());
	if (!Base::is_valid()) {
		return false;
	}

	for (typename Face_const_range::const_iterator it = Base::faces_const_range().begin(); it != Base::faces_const_range().end(); ++it) {
		Anchor const & current = anchor(it);
		Dart_const_descriptor current_dart = current.dart;

		for (unsigned i = 0; i < NB_SIDES; ++i) {
			Dart_const_descriptor opposite_dart = Base::const_opposite(current_dart);
			Point & c1 = current.vertices[i];
			Point & a1 = current.vertices[ccw(i)];
			Point & b1 = current.vertices[cw(i)];
			Complex cross_ratio = Base::get_cross_ratio(current_dart);
			Point d1 = Base::fourth_point_from_cross_ratio(a1, b1, c1, cross_ratio);

			unsigned j = index_in_anchor(opposite_dart);
			Anchor const & neighbor = anchor(opposite_dart);
			Point & a2 = neighbor.vertices[j];
			Point & c2 = neighbor.vertices[ccw(j)];
			Isometry pair_sides = isometry_pairing_the_sides<Traits>(a2, c2, a1, c1);
			CGAL_assertion(pair_sides.evaluate(a2) == a1);
			CGAL_assertion(pair_sides.evaluate(c2) == c1);

			Point d2 = pair_sides.evaluate(anchor(opposite_dart).vertices[ccw(j)]);
			CGAL_assertion(d2 == d1);
			if (d2 != d1) {
				return false;
			}

			current_dart = Base::const_ccw(current_dart);
		}
	}
	return true;
}

template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
to_stream(std::ostream & s) const
{
	CGAL_precondition(is_valid());

	// Give indices to the darts
	std::map<Dart_const_descriptor, unsigned> darts_ids;
	unsigned current_dart_id = 0;
	for (typename Dart_const_range::const_iterator it=this->combinatorial_map_.darts().begin(); it!=this->combinatorial_map_.darts().end(); ++it) {
		darts_ids[it] = current_dart_id;
		current_dart_id++;
	}

	// Store the number of darts
	s << current_dart_id << std::endl;

	// Store the triangles and their anchor
	for (typename Face_const_range::const_iterator it = Base::faces_const_range().begin(); it != Base::faces_const_range().end(); ++it) {
		s << darts_ids[it] << std::endl;
		s << darts_ids[Base::const_cw(it)] << std::endl;
		s << darts_ids[Base::const_ccw(it)] << std::endl;
		Anchor anch = anchor(it);
		s << darts_ids[anch.dart] << std::endl;
		for(int i = 0; i < NB_SIDES ; ++i) {
			s << anch.vertices[i].x() << std::endl;
			s << anch.vertices[i].y() << std::endl;
		}
	}

	// Store the edges and their cross-ratio
	for (typename Edge_const_range::const_iterator it = Base::edges_const_range().begin(); it != Base::edges_const_range().end(); ++it) {
		s << darts_ids[it] << std::endl;
		s << darts_ids[Base::const_opposite(it)] << std::endl;
		s << Base::get_cross_ratio(it);
	}
}


template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
from_stream(std::istream & s)
{
	this->combinatorial_map_.clear();

	// Load the number of darts
	std::string line;
	s >> line;
	int nb_darts = std::stoi(line);

	// Load the triangles and the anchors
	Dart_descriptor darts_by_id[nb_darts];
	for (int k = 0; k < nb_darts / NB_SIDES; ++k) {
		Dart_descriptor triangle_dart = this->combinatorial_map_.make_combinatorial_polygon(3);
		s >> line;
		int id_0 = std::stoi(line);
		s >> line;
		int id_1 = std::stoi(line);
		s >> line;
		int id_2 = std::stoi(line);
		darts_by_id[id_0] = triangle_dart;
		darts_by_id[id_1] = Base::cw(triangle_dart);
		darts_by_id[id_2] = Base::ccw(triangle_dart);
		
		Anchor anch = Anchor();
		s >> line;
		anch.dart = darts_by_id[std::stoi(line)];
		for(int i = 0; i < NB_SIDES; ++i) {
			Number x; 
			s >> x;
			// Number x >> line;
			Number y;
			s >> y;
			// Number y >> line;
			anch.vertices[i] = Point(x, y);
		}
		set_attribute(anch.dart, anch);
  }

	// Load the edges
	for (int k = 0; k < nb_darts / 2; ++k) {
		s >> line;
		int id_0 = std::stoi(line);
		s >> line;
		int id_1 = std::stoi(line);
		Dart_descriptor dart_0 = darts_by_id[id_0];
		Dart_descriptor dart_1 = darts_by_id[id_1];
    	this->combinatorial_map_.template sew<2>(dart_0, dart_1);
    	Complex cross_ratio;
    	s >> cross_ratio;
    	set_attribute(dart_0, cross_ratio);
  }

  CGAL_assertion(is_valid());
}

//---------- location and insertion

//TODO put this out of this class
template<class Traits>
Orientation
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
hyperbolic_orientation_2(Point const & p, Point const & q, Point const & r) const
{
	Point origin = Point(0, 0);
	Orientation orientation_to_origin = orientation(p, origin, q);
	if (orientation_to_origin == COLLINEAR) {
		return orientation(p, q, r);
	}

	Traits gt;
	CGAL::internal::Side_of_oriented_hyperbolic_segment_2 orientation_test = gt.side_of_oriented_hyperbolic_segment_2_object();
	Oriented_side orientation_to_disk = orientation_test(p, q, r);
	if (orientation_to_disk == ON_POSITIVE_SIDE) {
		return orientation_to_origin;
	} else if (orientation_to_disk == ON_NEGATIVE_SIDE) {
		if (orientation_to_origin == LEFT_TURN) {
			return RIGHT_TURN;
		} else{
			return LEFT_TURN;
		}
	} else{
		return COLLINEAR;
	}
}

// Output: The locate type lt of query relative to the anchor, and an index corresponding to:
// - if lt == FACE: NULL_INDEX (= -1),
// - if lt == EDGE: index of the edge on which query lies,
// - if lt == VERTEX: index of the vertex on which query lies,
// - if lt == OUTSIDE: index of the first edge such that query and the third point of the triangle lies on different sides.
template<class Traits>
std::tuple<Locate_type, unsigned>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
relative_position(Point const & query, Anchor const & anch) const
{
	Locate_type lt = FACE;
	unsigned index = NULL_INDEX;
	for (unsigned i = 0; i < NB_SIDES; ++i) {
		Orientation ori_query = hyperbolic_orientation_2(anch.vertices[i], anch.vertices[ccw(i)], query);
		if (ori_query == RIGHT_TURN) {
			lt = OUTSIDE;
			index = i;
			break;
		}
		if (ori_query == COLLINEAR) {
			lt = EDGE;
			index = i;
			if (hyperbolic_orientation_2(anch.vertices[ccw(i)], anch.vertices[cw(i)], query) == COLLINEAR) {
				lt = VERTEX;
				index = ccw(i);
				break;
			}
		}
	}
	return std::make_tuple(lt, index);
}

// Output: an anchor of the triangle in which query lies,
// and an int corresponding to the number of traversed triangles to find it.
template<class Traits>
std::tuple<typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor, unsigned>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
locate_visibility_walk(Point const & query, Anchor const & anch)
{
	Complex z_query (query.x(), query.y());
	CGAL_precondition(norm(z_query) < Number(1));
	CGAL_expensive_precondition(Base::is_Delaunay()); 

	// initialisation
	auto [lt, index] = relative_position(query, anch);
	if (lt != OUTSIDE){
		return std::tuple(anch, 0);
	}

	bool found = false;
	unsigned count = 0;

	Dart_descriptor dart = anch.dart;
	for (unsigned i = 0; i < index; ++i) {
		dart = Base::ccw(dart);
	}
	Point c = anch.vertices[index];
	Point a = anch.vertices[ccw(index)];
	Point b = anch.vertices[cw(index)];
	Point d;

	// visibility walk
	while (!found) {
		dart = Base::opposite(dart);
		Complex cross_ratio = Base::get_cross_ratio(dart);
		d = Base::fourth_point_from_cross_ratio(a, b, c, cross_ratio);
		if (hyperbolic_orientation_2(c, d, query) == RIGHT_TURN){
			b = a;
			a = d;
			dart = Base::ccw(dart);
		} else if (hyperbolic_orientation_2(d, a, query) == RIGHT_TURN) {
			b = c;
			c = d;
			dart = Base::cw(dart);
		} else {
			found = true;
		}
		++count;
	}
	Anchor res = Anchor(dart, a, c, d);
	CGAL_assertion(std::get<0>(relative_position(query, res)) != OUTSIDE);
	return std::tuple(res, count);
}

template<class Traits>
std::tuple<typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor, unsigned>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
locate_visibility_walk(Point const & query)
{
	return locate_visibility_walk(query, anchor());
}

template<class Traits>
std::tuple<typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor, unsigned>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
locate_straight_walk(Point const & query, Anchor const & anch)
{
	Complex z_query (query.x(), query.y());
	CGAL_precondition(norm(z_query) < Number(1));

	Dart_descriptor dart = anch.dart;
	Point p = anch.vertices[0];
	Point r = anch.vertices[1];
	Point l = anch.vertices[2];

	// initialisation
	unsigned count = 0;
	if (hyperbolic_orientation_2(r, p, query) < 0) {
		while (hyperbolic_orientation_2(l, p, query) < 0) {
			dart = Base::const_opposite(Base::cw(dart));
			Point old_l = l;
			l = Base::fourth_point_from_cross_ratio(p, r, l, Base::get_cross_ratio(dart));
			r = old_l;
			++count;
		}
	} else {
		while (hyperbolic_orientation_2(r, p, query) >= 0) {
			Point old_r = r;
			r = Base::fourth_point_from_cross_ratio(r, l, p, Base::get_cross_ratio(dart));
			l = old_r;
			dart = Base::ccw(Base::const_opposite(dart));
			++count;
		}
	}

	// straight walk
	dart = Base::ccw(dart);
	Point p0 = p;
	while (hyperbolic_orientation_2(query, r, l) < 0) {
		Complex cross_ratio = Base::get_cross_ratio(dart);
		Point s = Base::fourth_point_from_cross_ratio(l, p0, r, cross_ratio);
		if (hyperbolic_orientation_2(s, p, query) < 0) {
			p0 = r;
			r = s;
			dart = Base::cw(Base::opposite(dart));
		} else {
			p0 = l;
			l = s;
			dart = Base::ccw(Base::opposite(dart));
		}
		count++;
	}
	Anchor res = Anchor(dart, r, l, p0);
	CGAL_assertion(std::get<0>(relative_position(query, res)) != OUTSIDE);
	return std::tuple(res, count);
}

// Input: dart whose cross-ratio is those of the edge [t, r] and is computed with s and a fourth vertex.
// Output: modified anch such that its vertices are {r, query, t} and its dart represents the edge between r and t.
// The cross-ratio of the edge [t, r] is updated with its value where query is replaced by s.
// The cross-ratio of the edge [r, q] is set.
// Warning: the other cross-ratios are not updated, they will be in insert_in_edge and insert_in_face.
template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor &
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
update_infos(Dart_descriptor dart, Point const & r, Point const & s, Point const & t, Point const & query)
{       
	// compute new cross-ratio
	Complex old_cr = Base::get_cross_ratio(dart);
	Point u = Base::fourth_point_from_cross_ratio(r, s, t, old_cr);
	Complex new_cr = Base::cross_ratio(r, query, t, u);

	// modify anch and cross-ratio
	Anchor & anch = anchor(dart);
	anch.dart = dart;
	this->combinatorial_map_.template info<1>(dart) = new_cr;
	anch.vertices[0] = t;
	anch.vertices[1] = r;
	anch.vertices[2] = query;
	this->combinatorial_map_.template info<2>(dart) = anch;

	set_attribute(Base::ccw(dart), Base::cross_ratio(query, t, r, s));

	return anch;
}

// Output: the three anchors incident to query after its insertion inside anch.
// the darts of the new anchors correspond to the edges of the triangle the point was inserted in
// Note to self: No need to manage triangulation's anchor as it is done in update_infos.
template<class Traits>
std::vector<typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
insert_in_face(Point const & query, Anchor& anch)
{       
	CGAL_precondition(std::get<0>(relative_position(query, anch)) == FACE);
	
	Dart_descriptor current_dart = anch.dart;
	this->combinatorial_map_.insert_cell_0_in_cell_2(anch.dart);
	std::vector<Anchor> new_anchors;
	for (unsigned i = 0; i < NB_SIDES; ++i) {
		Point & c = anch.vertices[i];
		Point & a = anch.vertices[ccw(i)];
		Point & b = anch.vertices[cw(i)];
		new_anchors.push_back(update_infos(current_dart, a, b, c, query));
		current_dart = Base::ccw(Base::opposite(Base::ccw(current_dart)));
	}

	return new_anchors;
}

// Output: the four anchors incident to query after its insertion on one edge of anch.
// the darts of the new anchors correspond to the edges of the triangle the point was inserted in
template<class Traits>
std::vector<typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
insert_in_edge(Point const & query, Anchor& anch, unsigned index)
{
	CGAL_precondition(std::get<0>(relative_position(query, anch)) == EDGE);

	// find dart on which we insert query
	Dart_descriptor insertion_dart = ith_dart(index, anch);

	// gather information
	std::vector<Anchor> new_anchors;
	Point & c = anch.vertices[index];
	Point & a = anch.vertices[ccw(index)];
	Point & b = anch.vertices[cw(index)];
	Complex cross_ratio = Base::get_cross_ratio(insertion_dart);
	Point d = Base::fourth_point_from_cross_ratio(a, b, c, cross_ratio);

	Dart_descriptor dart_ab = Base::ccw(insertion_dart);
	Dart_descriptor dart_bc = Base::ccw(dart_ab);
	Dart_descriptor dart_cd = Base::ccw(Base::opposite(insertion_dart));
	Dart_descriptor dart_da = Base::ccw(dart_cd);

	// insert vertex on edge in the cmap and create new triangles
	this->combinatorial_map_.insert_cell_0_in_cell_1(insertion_dart);
	this->combinatorial_map_.insert_cell_1_in_cell_2(dart_bc, Base::cw(dart_ab));
	this->combinatorial_map_.insert_cell_1_in_cell_2(dart_da, Base::cw(dart_cd));
	CGAL_assertion(Base::ccw(dart_ab) == Base::opposite(Base::cw(dart_bc)));
	CGAL_assertion(Base::ccw(dart_bc) == Base::opposite(Base::cw(dart_cd)));
	CGAL_assertion(Base::ccw(dart_cd) == Base::opposite(Base::cw(dart_da)));
	CGAL_assertion(Base::ccw(dart_da) == Base::opposite(Base::cw(dart_ab)));

	// set new anchors and cross-ratios
	new_anchors.push_back(update_infos(dart_ab, b, c, a, query));
	new_anchors.push_back(update_infos(dart_bc, c, a, b, query));
	new_anchors.push_back(update_infos(dart_cd, d, a, c, query));
	new_anchors.push_back(update_infos(dart_da, a, c, d, query));

	return new_anchors;
}

// the darts of the new anchors correspond to the edges of the triangle the point was inserted in 
template<class Traits>
std::vector<typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
insert(Point const & query, Anchor& anch)
{
	Anchor locate_anchor = std::get<0>(locate_visibility_walk(query, anch));
	// Anchor locate_anchor = std::get<0>(locate_straight_walk(query, anch));
	auto [lt, index] = relative_position(query, locate_anchor);
	CGAL_precondition(lt != OUTSIDE);

	std::vector<Anchor> new_anchors;
	if (lt == FACE) {
		new_anchors = insert_in_face(query, locate_anchor);
	} else if (lt == EDGE) {
		new_anchors = insert_in_edge(query, locate_anchor, index);
	}
	return new_anchors;
}

// the darts of the new anchors correspond to the edges of the triangle the point was inserted in
template<class Traits>
std::vector<typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
insert(Point const & query)
{
	return insert(query, anchor());
}



//---------- Delaunay related methods


template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
flip(Dart_descriptor dart)
{
	CGAL_expensive_precondition(is_valid());
	CGAL_precondition(!Base::has_anchor());

	// first gather all the information needed
	Dart_descriptor opposite = Base::opposite(dart);
	Anchor & anch = anchor(dart);

	unsigned index = index_in_anchor(dart);
	Point & C = anch.vertices[index];
	Point & A = anch.vertices[ccw(index)];
	Point & B = anch.vertices[cw(index)];
	Complex cross_ratio = Base::get_cross_ratio(dart);
	Point D = Base::fourth_point_from_cross_ratio(A, B, C, cross_ratio);

	// flip the edge in the cmap and update the cross-ratios
	Base::flip(dart);  // AC -> BD

	// update the two anchors
	// WARNING: do not update dart's anchor first,
	// else it will change vertices' values (due to refs) and cause errors
	this->combinatorial_map_.template info<2>(opposite) = Anchor(opposite, D, B, C);
	this->combinatorial_map_.template info<2>(dart) = Anchor(dart, B, D, A);

	CGAL_expensive_assertion(is_valid());
}

// Pushes dart in the list darts_to_flip if dart is flippable and is not already in the list
template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
push_flippable_edge(Dart_descriptor const dart, std::list<Dart_descriptor>& darts_to_flip)
{       
	if (Base::is_Delaunay_flippable(dart)) {
		bool already_there = false;
		for (Dart_descriptor const& dart_to_flip : darts_to_flip) {
			if (dart_to_flip == dart || dart_to_flip == Base::opposite(dart)) {
				already_there = true;
				break;
			}
		}
		if (!already_there) {
			darts_to_flip.push_back(dart);
		}
	}
}

// Output: number of flips done to make the triangulation Delaunay given a list of darts to flip
template<class Traits>
std::tuple<unsigned, std::vector<typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Dart_descriptor>>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
restore_Delaunay(std::list<Dart_descriptor>& darts_to_flip)
{       
	unsigned number_of_flips_done = 0;
	std::vector<Dart_descriptor> flipped_darts;  // not useless: used in epsilon-net algo
	while (!darts_to_flip.empty()) {
		Dart_descriptor current_dart = darts_to_flip.front();
		if (Base::is_Delaunay_flippable(current_dart)) {
			flip(current_dart);
			flipped_darts.push_back(current_dart);
			number_of_flips_done++;
			Dart_descriptor maybe_flippable[4] = {Base::ccw(current_dart), Base::cw(current_dart),
											  Base::ccw(Base::opposite(current_dart)), Base::cw(Base::opposite(current_dart))};
			for (int i = 0; i < 4; ++i) {
				push_flippable_edge(maybe_flippable[i], darts_to_flip);
			}
		}
		darts_to_flip.pop_front();
	}

	CGAL_expensive_assertion(is_valid());
	CGAL_expensive_assertion(Base::is_Delaunay());
	return std::make_tuple(number_of_flips_done, flipped_darts);
}

template<class Traits>
unsigned
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
make_Delaunay()
{
	unsigned number_of_flips_done = 0;

	Dart_descriptor edge_to_flip = Base::pick_edge_to_flip();
	while (edge_to_flip != nullptr) {
		flip(edge_to_flip);
		edge_to_flip = Base::pick_edge_to_flip();
		number_of_flips_done++;
	}

	CGAL_expensive_assertion(is_valid());
	CGAL_expensive_assertion(Base::is_Delaunay());
	return number_of_flips_done;
}

// Inserts query in the triangulation, with a search starting from the given anchor, and makes the triangulation Delaunay again
// Output: number of flips done to make the triangulation Delaunay after the insertion
template<class Traits>
std::tuple<unsigned, std::vector<typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Dart_descriptor>>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
Delaunay_insert(Point const & query, Anchor& anch)
{
	CGAL_expensive_precondition(Base::is_Delaunay());

	std::vector<Anchor> new_anchors = insert(query, anch);
	std::list<Dart_descriptor> darts_to_flip;
	for (int i = 0; i < new_anchors.size(); ++i) {
		push_flippable_edge(new_anchors[i].dart, darts_to_flip);
	}

	return restore_Delaunay(darts_to_flip);
}

// Same but the search starts from main anchor
template<class Traits>
std::tuple<unsigned, std::vector<typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Dart_descriptor>>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
Delaunay_insert(Point const & query)
{
	return Delaunay_insert(query, anchor());
}



//---------- eps-net methods


template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Number
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
delta(Point const & u, Point const & v) const
{
	Number num = (u.x() - v.x()) * (u.x() - v.x()) + (u.y() - v.y()) * (u.y() - v.y());
	Number den = (1 - (u.x() * u.x() + u.y() * u.y())) * (1 - (v.x() * v.x() + v.y() * v.y()));
	return 2 * num / den;
}

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Point
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
approx_circumcenter(Anchor const & anch) const
{
	Traits gt;
	CGAL::internal::Construct_hyperbolic_circumcenter_CK_2<Traits> chc(gt);
	// CGAL::internal::Construct_hyperbolic_circumcenter_2<Traits> chc(gt);
	Voronoi_point exact_center = chc(anch.vertices[0], anch.vertices[1], anch.vertices[2]);
	Number x = to_double(exact_center.x());
	Number y = to_double(exact_center.y());
	// Number x = to_double(exact_center.x().exact());
	// Number y = to_double(exact_center.y().exact());
	// std::cout << std::setprecision(17) << to_interval(exact_center.x()).first << " ; " << to_interval(exact_center.x()).second << std::endl;
	return Point(x, y);
}

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Number
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
delta_min(Anchor const & anch, Point const & approx_center) const
{
	Number res = Number(999);
	for (unsigned i = 0; i < NB_SIDES; ++i) {
		res = std::min(delta(anch.vertices[i], approx_center), res);
	}
	return res;
}

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Number
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
delta_max(Anchor const & anch, Point const & approx_center) const
{
	Number res = Number(0);
	for (unsigned i = 0; i < NB_SIDES; ++i) {
		res = std::max(delta(anch.vertices[i], approx_center), res);
	}
	return res;
}

template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
push_triangle(Dart_descriptor const dart, std::list<Dart_descriptor> & triangles, size_t & triangles_list_mark)
{
	Dart_descriptor current_dart = dart;
	for (unsigned i = 0; i < NB_SIDES; ++i){
		this->combinatorial_map_.unmark(current_dart, triangles_list_mark);
		current_dart = Base::ccw(current_dart);
	}

	this->combinatorial_map_.mark(dart, triangles_list_mark);
	triangles.push_back(dart);
}

template<class Traits>
bool
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
epsilon_net(const double epsilon)
{       
	const double BOUND = std::cosh(epsilon)-1;
	CGAL_assertion(is_epsilon_packing(epsilon));
	size_t triangles_list_mark = this->combinatorial_map_.get_new_mark();
		
	std::list<Dart_descriptor> triangles;
	for (typename Face_range::iterator it = this->combinatorial_map_.template one_dart_per_cell<2>().begin();
		it != this->combinatorial_map_.template one_dart_per_cell<2>().end(); ++it) {
		Anchor& current_anchor = anchor(it);
		Point current_center = approx_circumcenter(current_anchor);
		if (delta_min(current_anchor, current_center) > BOUND) {
			push_triangle(it, triangles, triangles_list_mark);
		}
	}

	while (!triangles.empty()) {
		Dart_descriptor current_dart = triangles.front();
		Anchor& current_anchor = anchor(current_dart); 
		triangles.pop_front();
		if(this->combinatorial_map_.is_marked(current_dart, triangles_list_mark)){
			this->combinatorial_map_.unmark(current_dart, triangles_list_mark);
			Point current_center = approx_circumcenter(current_anchor);
			if (delta_min(current_anchor, current_center) > BOUND) {
				std::vector<Anchor> new_anchors = insert(current_center, current_anchor);
				std::list<Dart_descriptor> darts_to_flip;
				for (Anchor const & new_anchor : new_anchors) {
					push_triangle(new_anchor.dart, triangles, triangles_list_mark);
					push_flippable_edge(new_anchor.dart, darts_to_flip);  //the darts of the new anchors correspond to the edges of the triangle the point was inserted in 
				}

				std::vector<Dart_descriptor> flipped_darts = std::get<1>(restore_Delaunay(darts_to_flip));
				for (Dart_descriptor const& dart : flipped_darts) {
					push_triangle(dart, triangles, triangles_list_mark);
					push_triangle(Base::opposite(dart), triangles, triangles_list_mark);
				}
			}
		}
	}
	this->combinatorial_map_.free_mark(triangles_list_mark);

	bool is_covering = is_epsilon_covering(epsilon);
	bool is_packing = is_epsilon_packing(epsilon);
	CGAL_assertion(is_covering);
	CGAL_assertion(is_packing);
	return is_covering&&is_packing;
}

template<class Traits>
bool
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
is_epsilon_covering(const double epsilon)
{
	Interval delta_epsilon = cosh(epsilon)-1;
	bool is_covering = true;
	for (typename Face_range::iterator it = this->combinatorial_map_.template one_dart_per_cell<2>().begin();
		it != this->combinatorial_map_.template one_dart_per_cell<2>().end(); ++it) {
		Anchor& current_anchor = anchor(it);

		Traits gt;
		CGAL::internal::Construct_hyperbolic_circumcenter_CK_2<Traits> chc(gt);
		// CGAL::internal::Construct_hyperbolic_circumcenter_2<Traits> chc(gt);
		Voronoi_point v = chc(current_anchor.vertices[0], current_anchor.vertices[1], current_anchor.vertices[2]);

		Point u = current_anchor.vertices[0];
		auto num = (u.x() - v.x()) * (u.x() - v.x()) + (u.y() - v.y()) * (u.y() - v.y());
		auto den = (1 - (u.x() * u.x() + u.y() * u.y())) * (1 - (v.x() * v.x() + v.y() * v.y()));
		auto d = 2 * num / den;

		auto upper_bound = Number(upper(delta_epsilon));
		if (d > upper_bound) {
			is_covering = false;
			break;
		}
	}
	return is_covering;
}

template<class Traits>
bool
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
is_epsilon_packing(const double epsilon)
{
	Interval delta_epsilon = cosh(epsilon)-1;
	bool is_packing = true;
	for (typename Edge_range::iterator it = this->combinatorial_map_.template one_dart_per_cell<1>().begin();
									   it != this->combinatorial_map_.template one_dart_per_cell<1>().end(); ++it) {
		Anchor& current_anchor = anchor(it);
		unsigned index = index_in_anchor(it);
		Point & a = current_anchor.vertices[index];
		Point & b = current_anchor.vertices[ccw(index)];
		
		if (delta(a, b) < lower(delta_epsilon)) {
			is_packing = false;  // consider that it's not a packing
			Dart_descriptor next = Base::ccw(it);
			auto doc = this->combinatorial_map_.template darts_of_cell<0>(it);
			// check if the edge is a loop:
			// if next and it are darts of the same vertex, then a = b
			for (auto dart = doc.begin(); dart != doc.end(); ++dart) {
				if (next == dart) {
					is_packing = true;  // the edge is a loop, so actually it's ok
				}
			}
		}
	}
	return is_packing;
}

template<class Traits>
bool
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
is_epsilon_net(const double epsilon)
{
	return is_epsilon_covering(epsilon) && is_epsilon_packing(epsilon);
}

template<class Traits>
double
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
shortest_loop() const
{
	Number min_delta_length = 999;
	for (typename Edge_range::iterator it = this->combinatorial_map_.template one_dart_per_cell<1>().begin();
									   it != this->combinatorial_map_.template one_dart_per_cell<1>().end(); ++it) {
		Anchor& current_anchor = anchor(it);
		unsigned index = index_in_anchor(it);
		Point & a = current_anchor.vertices[index];
		Point & b = current_anchor.vertices[ccw(index)];
		
		Dart_descriptor next = Base::ccw(it);
		auto doc = this->combinatorial_map_.template darts_of_cell<0>(it);
		// check if the edge is a loop:
		// if next and it are darts of the same vertex, then a = b
		for (auto dart = doc.begin(); dart != doc.end(); ++dart) {
			if (next == dart) {
				Number delta_length = delta(a, b);  // the edge is a loop
				min_delta_length = min(min_delta_length,delta_length);
			}
		}
	}
	return std::acosh(1 + to_double(min_delta_length));
}

//---------- GLOBAL FUNCTIONS

template<class Traits>
std::ostream&
operator<<(std::ostream& s, const Delaunay_triangulation_on_hyperbolic_surface_2<Traits>& triangulation)
{
  triangulation.to_stream(s);
  return s;
}

template<class Traits>
void
operator>>(std::istream& s, Delaunay_triangulation_on_hyperbolic_surface_2<Traits>& triangulation)
{
  triangulation.from_stream(s);
}

}  // namespace CGAL

#endif  //CGAL_DELAUNAY_TRIANGULATION_ON_HYPERBOLIC_SURFACE_2