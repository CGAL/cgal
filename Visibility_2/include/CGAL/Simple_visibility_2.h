#ifndef CGAL_SIMPLE_VISIBILITY_2_H
#define CGAL_SIMPLE_VISIBILITY_2_H

#include <CGAL/Arrangement_2.h>
#include <stack>

namespace CGAL {

namespace Visibility_2 {

template<class Arrangement_2> 
class Simple_visibility_2 {

public:
	typedef typename Arrangement_2::Geometry_traits_2				Geometry_traits_2;
	// Currently only consider with same type for both
	typedef Arrangement_2											Input_Arrangement_2;
	typedef Arrangement_2											Output_Arrangement_2;

	typedef typename Arrangement_2::Halfedge_const_handle			Halfedge_const_handle;
	typedef typename Arrangement_2::Ccb_halfedge_const_circulator 	Ccb_halfedge_const_circulator;
	typedef typename Arrangement_2::Face_const_handle				Face_const_handle;

	typedef typename Geometry_traits_2::Point_2						Point_2;
	typedef typename Geometry_traits_2::Ray_2						Ray_2;
	typedef typename Geometry_traits_2::Segment_2					Segment_2;

	Simple_visibility_2() : p_arr(NULL) {};

	/*! Constructor given an arrangement and the Regularization tag. */
	Simple_visibility_2(const Input_Arrangement_2 &arr/*, Regularization_tag r_t*/): p_arr(&arr) {};

	bool is_attached() {
		return (p_arr != NULL);
	}

	void attach(const Input_Arrangement_2 &arr) {
		p_arr = &arr;
	}

	void detach() {
		p_arr = NULL;
	}

	Input_Arrangement_2 arr() {
		return *p_arr;
	}

	void visibility_region(Point_2 &q, 
						   const Face_const_handle face,
						   Output_Arrangement_2 &out_arr
						   ) {

		int i = 0;
		bool ccw = false;

		typename Input_Arrangement_2::Ccb_halfedge_const_circulator circ = face->outer_ccb();
		typename Input_Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
  		typename Input_Arrangement_2::Halfedge_const_handle he = curr;

  		int parity = 0;

  		std::vector<Point_2> temp_vertices;
  		int index_v0 = 0;
  		int index = 0;
  		Point_2 curr_min = he->source()->point();

  		// Push all vertices
  		do {
			he = curr;  		
			Point_2 curr_vertex = he->target()->point();
			if (curr_vertex.x() < curr_min.x() && (curr_vertex.x() > q.x())) {
				curr_min = curr_vertex;
				index_v0 = index;
			}
			temp_vertices.push_back(curr_vertex);
			index++;
  		} while (++curr != circ);

  		// Now create vector so that first vertex v0 has the smallest positive x-coordinate (thus - visible from the query point)
  		for (unsigned int k = index_v0 ; k < temp_vertices.size() ; k++) {
  			vertices.push_back(temp_vertices[k]);
  		}
  		for (unsigned int k = 0 ; k < index_v0 ; k++) {
  			vertices.push_back(temp_vertices[k]);
  		}
  		// Push first vertex again to fulfill algo precondition
  		vertices.push_back(vertices[0]);

  		std::cout << "Vertices:\n";
  		for (unsigned int k = 0 ; k < vertices.size() ; k++) {
  			std::cout << vertices[k] << std::endl;
  		}

		Point_2 w;

		if (orientation(q, vertices[0], vertices[1]) == CGAL::LEFT_TURN) {
			std::cout << "left" << std::endl;
			upcase = LEFT;
			i = 1;
			w = vertices[1];
			s.push(vertices[0]);
			s.push(vertices[1]);
		}
		else {
			std::cout << "scana" << std::endl;
			upcase = SCANA;
			i = 1;
			w = vertices[1];
			s.push(vertices[0]);
		}
		int counter = 0;
		do {
			std::cout << "CASE: " << upcase << std::endl;
			switch(upcase) {
				case LEFT: 
					std::cout << "before entering left" << std::endl;
					left(i, w, q);
					std::cout << "after exiting left" << std::endl;
					break;
				case RIGHT:
					right(i, w, q);
					break;
				case SCANA:
					scana(i, w, q);
					break;
				case SCANB:
					scanb(i, w, q);
					break;
				case SCANC:
					scanc(i, w, q);
					break;
				case SCAND:
					scand(i, w, q);
					break;
			}
			if (upcase == LEFT) {
				// Check if (s_t-1, s_t) intersects (q, vn) 
				Point_2 s_t = s.top();
				std::cout << "s_t= " << s_t << std::endl;
				s.pop();
				Point_2 s_t_prev = s.top();
				std::cout << "s_t-1= " << s_t_prev << std::endl;
				Segment_2 s1(s_t_prev, s_t);
				Segment_2 s2(q, vertices[vertices.size()-1]);
				CGAL::Object result = CGAL::intersection(s1, s2);

				if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) { 
					Segment_2 s3(s_t_prev, vertices[i]);
					CGAL::Object result2 = CGAL::intersection(s3, s2);
					if (const Point_2 *vertex_new = CGAL::object_cast<Point_2>(&result2)) {
						if ((*vertex_new) != (s_t_prev) && (*vertex_new != s_t)) {
							std::cout << "switch to scanb" << std::endl;
							upcase = SCANB;
							s.push(*vertex_new);
						}
						else { // Do not alter stack if it doesn't intersect - push back s_t
							std::cout << "skipping because it is endpoint" << std::endl;
							s.push(s_t);
						}
					}
					else {
						s.push(s_t);
					}
				}
				else {
					s.push(s_t);
				}
			}
			std::cout << "gets out" << std::endl;
			if (counter == 2) {
			//	exit(0);	
			}
			counter++;
		} while(upcase != FINISH);

		std::cout << "RESULT: " << std::endl;
		typename std::list<Segment_2> segments;
		if (!s.empty()) {
			Point_2 prev_pt = s.top();
			s.pop();
			while(!s.empty()) {
				Point_2 curr_pt = s.top();
				segments.push_front(Segment_2(curr_pt, prev_pt));
				prev_pt = curr_pt;
				s.pop();
			}
		}
		for (typename std::list<Segment_2>::iterator it = segments.begin(); it != segments.end(); it++) {
		    std::cout << it->source() << " " << it->target() << std::endl;
		}
		CGAL::insert(out_arr, segments.begin(), segments.end());
	}

	void visibility_region(const Point_2 &q, 
						   const Halfedge_const_handle he,
						   Output_Arrangement_2 &out_arr
						   ) {

	}
protected:
	const Input_Arrangement_2 *p_arr;
	std::stack<Point_2> s;
	std::vector<Point_2> vertices;
	enum {LEFT, RIGHT, SCANA, SCANB, SCANC, SCAND, FINISH} upcase;

	bool do_overlap(const Point_2 &a, const Point_2 &b, const Point_2 &c) {
		if (collinear(a, b, c)) {
			std::cout << a << " " << b << " " << c << " are collinear" << std::endl;
			Segment_2 s1(a, b);
			Segment_2 s2(a, c);
			const Segment_2 *seg_overlap;
			CGAL::Object result = CGAL::intersection(s1, s2);
			if (seg_overlap = CGAL::object_cast<Segment_2>(&result)) { 
					return true;
			}
		}
		return false;
	}

	void left(int &i, Point_2 &w, const Point_2 &query_pt) {
		std::cout << "begin with i = " << i << std::endl;
		if (i == vertices.size() - 1) {
			std::cout << "done" << std::endl;
			upcase = FINISH;
		}
		else if (orientation(query_pt, vertices[i], vertices[i+1]) == CGAL::LEFT_TURN) {
			std::cout << "left::left turn with i =" << i << std::endl;
			upcase = LEFT;
			s.push(vertices[i+1]);
			w = vertices[i+1];
			i++;
		}
		else if (orientation(query_pt, vertices[i], vertices[i+1]) == CGAL::RIGHT_TURN) {
			Point_2 s_t = s.top();
			s.pop();
			Point_2 s_t_prev = s.top();
			s.pop();
			if (orientation(s_t_prev, vertices[i], vertices[i+1]) == CGAL::RIGHT_TURN) {
				upcase = SCANA;
				i++;
			}
			s.push(s_t_prev);
			s.push(s_t);
			w = vertices[i+1];
		}
		else {
			upcase = RIGHT;
			i++;
			w = vertices[i];
		}
	}

	void right(int &i, Point_2 &w, const Point_2 &query_pt) {
		// Scan s_t, s_t-1, ..., s_1, s_0 for the first edge (s_j, s_j-1) such that
		// (a) (z, s_j, v_i) is a right turn and (z, s_j-1, v_i) is a left turn, or
		// (b) (z, s_j-1, s_j) is a forward move and (v_i-1, v_i) intersects (s_j-1, s_j)
		bool found = false;

		while(!found) {
			Point_2 s_j = s.top();
			if (!s.empty()) {
				s.pop();
				Point_2 s_j_prev = s.top();

				// Check condition (a)
				if ((orientation(query_pt, s_j, vertices[i]) == CGAL::RIGHT_TURN)
					&& (orientation(query_pt, s_j_prev, vertices[i]) == CGAL::LEFT_TURN)) {
					found = true;
					Segment_2 s1(s_j_prev, s_j);
					Segment_2 s2(query_pt, vertices[i]);
					CGAL::Object result = intersection(s1, s2);
					if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
						s_j = *ipoint;
					}

					if (orientation(query_pt, vertices[i], vertices[i+1]) == CGAL::RIGHT_TURN) {
						upcase = RIGHT;
						s.push(s_j);
						i++;
						w = vertices[i];
					}
					else if ((orientation(query_pt, vertices[i], vertices[i+1]) == CGAL::LEFT_TURN)
							&& (orientation(vertices[i-1], vertices[i], vertices[i+1]) == CGAL::RIGHT_TURN)) {

						upcase = LEFT;
						i++;
						s.push(s_j);
						s.push(vertices[i]);
						s.push(vertices[i+1]);
						w = vertices[i+1];
					}
					else {
						upcase = SCANC;
						i++;
						s.push(s_j);
					}
				}
				else { // Case (b)
					Segment_2 s1(s_j_prev, s_j);
					Segment_2 s2(vertices[i-1], vertices[i]);
					CGAL::Object result = intersection(s1, s2);
					if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
						// Keep s_j off the stack
						upcase = SCAND;
						w = *ipoint;
					}
				}
			}
		}
	}

	void scana(int &i, Point_2 &w, const Point_2 &query_pt) {
		// Scan v_i, v_i+1, ..., v_n for the first edge to intersect (z, s_t)
		bool found = false;
		int k = i;
		Point_2 intersection_pt;
		while (k+1 < vertices.size()-1) {
			std::cout << "Considering edge: " << vertices[k] << " " << vertices[k+1] << std::endl;
			Segment_2 s1(vertices[k], vertices[k+1]);
			std::cout << "With edge: " << query_pt << " " << s.top() << std::endl;
			Ray_2 s2(query_pt, s.top());
			CGAL::Object result = intersection(s1, s2);
			if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) { 
				found = true;
				std::cout << " found edge!" << std::endl;
				intersection_pt = *ipoint;
				std::cout << intersection_pt << std::endl;
				break;
			}
			k++;
		}
		if (found) {
			if ((orientation(query_pt, vertices[k], vertices[k+1]) == CGAL::RIGHT_TURN)
				&& (!do_overlap(query_pt, s.top(), intersection_pt))) {
				std::cout << "if1" << std::endl;
				upcase = RIGHT;
				i = k+1;
				w = intersection_pt;
			}
			else if ((orientation(query_pt, vertices[k], vertices[k+1]) == CGAL::RIGHT_TURN)
				&& (do_overlap(query_pt, s.top(), intersection_pt))) {
				std::cout << "elseif1" << std::endl;
				upcase = SCAND;
				i = k+1;
				w = intersection_pt;
			}
			else if ((orientation(query_pt, vertices[k], vertices[k+1]) == CGAL::LEFT_TURN)
				&& (do_overlap(query_pt, s.top(), intersection_pt))) {
				std::cout << "elseif2" << std::endl;
				upcase = LEFT;
				i = k+1;
				s.push(intersection_pt);
				s.push(vertices[k+1]);
				w = vertices[k+1];
			}
			else {
				std::cout << "should never occur" << std::endl;
				// This case never occurs
			}
		}
	}

	void scanb(int &i, Point_2 &w, const Point_2 &query_pt) {
		// Scan v_i, v_i+1, ..., v_n-1, v_n for the first edge to intersect (s_t, v_n]
		Point_2 s_t = s.top();
		int k = i;
		bool found = false;
		Point_2 intersection_pt;
		while (k+1 < vertices.size()-1) {
			Segment_2 s1(vertices[k], vertices[k+1]);
			Segment_2 s2(s_t, vertices[vertices.size()-1]);
			CGAL::Object result = intersection(s1, s2);
			if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) { 
				intersection_pt = *ipoint;
				found = true;
				break;
			}
			k++;
		}
		if (found) {
			if ((intersection_pt == vertices[k+1]) && (intersection_pt == vertices[vertices.size()-1])) {
				upcase = FINISH;
				w = vertices[vertices.size()-1];
				s.push(vertices[vertices.size()-1]);
			}
			else {
				upcase = RIGHT;
				i = k+1;
				w = intersection_pt;
			}
		}
	}

	void scanc(int &i,Point_2 &w, const Point_2 &query_pt) {
		// Scan v_i, v_i+1, ..., v_n-1, v_n for the first edge to intersect (s_t, w)
		Point_2 s_t = s.top();
		int k = i;
		bool found = false;
		Point_2 intersection_pt;
		while (k+1 < vertices.size()-1) {
			Segment_2 s1(vertices[k], vertices[k+1]);
			Segment_2 s2(s_t, w);
			CGAL::Object result = intersection(s1, s2);
			if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
				found = true;
				intersection_pt = *ipoint;
				break;
			}
			k++;
		}
		if (found) {
			upcase = RIGHT;
			i = k+1;
			w = intersection_pt;
		}
	}

	void scand(int &i, Point_2 &w, const Point_2 &query_pt) {
		// Scan v_i, v_i+1, v_n-1, v_n for the fist edge to intersect (s_t, w)
		Point_2 s_t = s.top();
		int k = i;
		bool found = false;
		Point_2 intersection_pt;
		while (k+1 < vertices.size()-1) {
			Segment_2 s1(vertices[k], vertices[k+1]);
			Segment_2 s2(s_t, w);
			CGAL::Object result = intersection(s1, s2);
			if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
				found = true;
				intersection_pt = *ipoint;
				break;
			}
			k++;
		}
		if (found) {
			upcase = LEFT;
			i = k+1;
			s.push(intersection_pt);
			s.push(vertices[k+1]);
			w = vertices[k+1];
		}
	}
};

} // namespace Visibility_2
} // namespace CGAL

#endif