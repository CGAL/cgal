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
	typedef typename Geometry_traits_2::Vector_2					Vector_2;
	typedef typename Geometry_traits_2::Segment_2					Segment_2;

	/*! Constructor given an arrangement and the Regularization tag. */
	Simple_visibility_2(const Input_Arrangement_2 &arr/*, Regularization_tag r_t*/): p_arr(arr) {};

	bool is_attached() {
		return (p_arr != NULL);
	}

	void attach(const Input_Arrangement_2 &arr) {
		p_arr = arr;
	}

	void detach() {
		p_arr = NULL;
	}

	Input_Arrangement_2 arr() {
		return p_arr;
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

  		std::vector<Point_2> translated_vertices;
  		int index_v0 = 0;
  		int index = 0;
  		double curr_min_x = 1000000;
  		// Push all vertices and translate coordinate system such that query point is at the origin
  		do {
			he = curr;  			
			Point_2 curr_vertex = he->target()->point();
			Point_2 translated_vertex = Point_2(curr_vertex.x() - q.x(), curr_vertex.y() - q.y());
			if (CGAL::to_double(translated_vertex.x()) < curr_min_x && CGAL::to_double(translated_vertex.x()) > 0) {
				curr_min_x = CGAL::to_double(translated_vertex.x());
				index_v0 = index;
			}
			translated_vertices.push_back(translated_vertex);
			index++;
  		} while (++curr != circ);

  		// Now create vector so that first vertex v0 has the smallest positive x-coordinate (thus - visible from the query point)
  		for (unsigned int k = index_v0 ; k < translated_vertices.size() ; k++) {
  			vertices.push_back(translated_vertices[k]);
  		}
  		for (unsigned int k = 0 ; k < index_v0 ; k++) {
  			vertices.push_back(translated_vertices[k]);
  		}
  		// Push first vertex again to fulfill algo precondition
  		vertices.push_back(vertices[0]);

  		std::cout << "Vertices:\n";
  		for (unsigned int k = 0 ; k < vertices.size() ; k++) {
  			std::cout << vertices[k] << std::endl;
  		}

		Point_2 w;

		if (CGAL::orientation(q, vertices[0], vertices[1]) == CGAL::LEFT_TURN) {
			i = 1;
			w = vertices[1];
			s.push(vertices[0]);
			s.push(vertices[1]);
			left(i, w, q);
		}
		else {
			i = 1;
			w = vertices[1];
			s.push(vertices[0]);
			scana(i, w, q);
		}
		std::cout << "gets here\n";
		Point_2 stored_q(q);
		q = Point_2(0, 0);
		do {
			switch(upcase) {
				case LEFT: 
					left(i, w, q);
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
				s.pop();
				Point_2 s_t_prev = s.top();
				Segment_2 s1(s_t_prev, s_t);
				Segment_2 s2(q, vertices[vertices.size()-1]);
				CGAL::Object result = CGAL::intersection(s1, s2);

				if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) { 
					upcase = SCANB;
					Segment_2 s3(s_t_prev, vertices[i]);
					CGAL::Object result2 = CGAL::intersection(s3, s2);
					if (const Point_2 *vertex_new = CGAL::object_cast<Point_2>(&result2)) {
						s.push(*vertex_new);
					}
				}
				else { // Do not alter stack if it doesn't intersect - push back s_t
					s.push(s_t);
				}
			}
		} while(upcase != FINISH);

		std::cout << "RESULT: " << std::endl;
		
		while(!s.empty()) {
			Point_2 curr_pt = s.top();
			Point_2 final_pt(curr_pt.x() + stored_q.x(), curr_pt.y() + stored_q.y());
			std::cout << final_pt << std::endl;
			s.pop();
		}
	}

	void visibility_region(const Point_2 &q, 
						   const Halfedge_const_handle he,
						   Output_Arrangement_2 &out_arr
						   ) {

	}
protected:
	Input_Arrangement_2 p_arr;
	std::stack<Point_2> s;
	std::vector<Point_2> vertices;
	enum {LEFT, RIGHT, SCANA, SCANB, SCANC, SCAND, FINISH} upcase;

	bool do_overlap(const Point_2 &a, const Point_2 &b, const Point_2 &c) {

		if (CGAL::collinear(a, b, c)) {
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

		if (i == vertices.size()) {
			upcase = FINISH;
		}
		else if (CGAL::orientation(query_pt, vertices[i], vertices[i+1]) == CGAL::LEFT_TURN) {
			upcase = LEFT;
			s.push(vertices[i+1]);
			w = vertices[i+1];
			i++;
		}
		else if (CGAL::orientation(query_pt, vertices[i], vertices[i+1]) == CGAL::RIGHT_TURN) {
			Point_2 s_t = s.top();
			s.pop();
			Point_2 s_t_prev = s.top();
			s.pop();
			if (CGAL::orientation(s_t_prev, vertices[i], vertices[i+1]) == CGAL::RIGHT_TURN) {
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
				if ((CGAL::orientation(query_pt, s_j, vertices[i]) == CGAL::RIGHT_TURN)
					&& (CGAL::orientation(query_pt, s_j_prev, vertices[i]) == CGAL::LEFT_TURN)) {
					found = true;
					Segment_2 s1(s_j_prev, s_j);
					Segment_2 s2(query_pt, vertices[i]);
					CGAL::Object result = CGAL::intersection(s1, s2);
					if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
						s_j = *ipoint;
					}

					if (CGAL::orientation(query_pt, vertices[i], vertices[i+1]) == CGAL::RIGHT_TURN) {
						upcase = RIGHT;
						s.push(s_j);
						i++;
						w = vertices[i];
					}
					else if ((CGAL::orientation(query_pt, vertices[i], vertices[i+1]) == CGAL::LEFT_TURN)
							&& (CGAL::orientation(vertices[i-1], vertices[i], vertices[i+1]) == CGAL::RIGHT_TURN)) {

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
					CGAL::Object result = CGAL::intersection(s1, s2);
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
		const Point_2 *ipoint;
		while (k+1 < vertices.size()-1) {
			Segment_2 s1(vertices[k], vertices[k+1]);
			Segment_2 s2(query_pt, s.top());
			CGAL::Object result = CGAL::intersection(s1, s2);
			if (ipoint = CGAL::object_cast<Point_2>(&result)) { 
				found = true;
				break;
			}
			k++;
		}
		if (found) {
			if ((CGAL::orientation(query_pt, vertices[k], vertices[k+1]) == CGAL::RIGHT_TURN)
				&& (do_overlap(query_pt, s.top(), *ipoint))) {

				upcase = RIGHT;
				i = k+1;
				w = *ipoint;
			}
			else if ((CGAL::orientation(query_pt, vertices[k], vertices[k+1]) == CGAL::RIGHT_TURN)
				&& (!do_overlap(query_pt, s.top(), *ipoint))) {

				upcase = SCAND;
				i = k+1;
				w = *ipoint;
			}
			else if ((CGAL::orientation(query_pt, vertices[k], vertices[k+1]) == CGAL::LEFT_TURN)
				&& (!do_overlap(query_pt, s.top(), *ipoint))) {

				upcase = LEFT;
				i = k+1;
				s.push(*ipoint);
				s.push(vertices[k+1]);
				w = vertices[k+1];
			}
			else {
				// This case never occurs
			}
		}
	}

	void scanb(int &i, Point_2 &w, const Point_2 &query_pt) {
		// Scan v_i, v_i+1, ..., v_n-1, v_n for the first edge to intersect (s_t, v_n]
		Point_2 s_t = s.top();
		int k = i;
		bool found = false;
		const Point_2 *ipoint;
		while (k+1 < vertices.size()-1) {
			Segment_2 s1(vertices[k], vertices[k+1]);
			Segment_2 s2(s_t, vertices[vertices.size()-1]);
			CGAL::Object result = CGAL::intersection(s1, s2);
			if (ipoint = CGAL::object_cast<Point_2>(&result)) { 
				found = true;
				break;
			}
			k++;
		}
		if (found) {
			if ((*ipoint == vertices[k+1]) && (*ipoint == vertices[vertices.size()-1])) {
				upcase = FINISH;
				w = vertices[vertices.size()-1];
				s.push(vertices[vertices.size()-1]);
			}
			else {
				upcase = RIGHT;
				i = k+1;
				w = *ipoint;
			}
		}
	}

	void scanc(int &i,Point_2 &w, const Point_2 &query_pt) {
		// Scan v_i, v_i+1, ..., v_n-1, v_n for the first edge to intersect (s_t, w)
		Point_2 s_t = s.top();
		int k = i;
		bool found = false;
		const Point_2 *ipoint;
		while (k+1 < vertices.size()-1) {
			Segment_2 s1(vertices[k], vertices[k+1]);
			Segment_2 s2(s_t, w);
			CGAL::Object result = CGAL::intersection(s1, s2);
			if (ipoint = CGAL::object_cast<Point_2>(&result)) {
				found = true;
				break;
			}
			k++;
		}
		if (found) {
			upcase = RIGHT;
			i = k+1;
			w = *ipoint;
		}
	}

	void scand(int &i, Point_2 &w, const Point_2 &query_pt) {
		// Scan v_i, v_i+1, v_n-1, v_n for the fist edge to intersect (s_t, w)
		Point_2 s_t = s.top();
		int k = i;
		bool found = false;
		const Point_2 *ipoint;
		while (k+1 < vertices.size()-1) {
			Segment_2 s1(vertices[k], vertices[k+1]);
			Segment_2 s2(s_t, w);
			CGAL::Object result = CGAL::intersection(s1, s2);
			if (ipoint = CGAL::object_cast<Point_2>(&result)) {
				found = true;
				break;
			}
			k++;
		}
		if (found) {
			upcase = LEFT;
			i = k+1;
			s.push(*ipoint);
			s.push(vertices[k+1]);
			w = vertices[k+1];
		}
	}
};

} // namespace Visibility_2
} // namespace CGAL

#endif