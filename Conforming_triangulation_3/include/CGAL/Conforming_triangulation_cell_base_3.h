//A tetrahedral cell in a conforming Delaunay triangulation.
//Copyright (C) 2012  Utrecht University
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Author(s): Thijs van Lankveld

// Based on CGAL/Triangulation_cell_base_3.h
//

// Note that a constrained 3D Delaunay triangulation is partially conforming and guided.
// 1-dimensional constraints are specified by segments.
// 2-demensional constraints are specified by closure complete polygons.
// The constrained 3D DT must contain a collection of edges and faces that exactly cover
// each constrained segment or polygon. For example, a constrained edge will not necessarily
// be maintained as vertices are inserted, but after updating the triangulation to the vertex
// insertion, the constrained edge must be the union of one or more constrained edges.

#ifndef CGAL_CONFORMING_TRIANGULATION_CELL_BASE_3_H
#define CGAL_CONFORMING_TRIANGULATION_CELL_BASE_3_H

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_cell_base_3.h>

CGAL_BEGIN_NAMESPACE

template < typename Gt, typename Cb = Triangulation_cell_base_3<Gt> >
class Conforming_triangulation_cell_base_3: public Cb {
public:
	typedef typename Cb::Vertex_handle		Vertex_handle;
	typedef typename Cb::Cell_handle		Cell_handle;

	typedef Gt								Geom_traits;
	typedef typename Geom_traits::Point_3	Point;

	typedef Point*							Point_container;
	typedef Point*							Point_iterator;
	typedef const Point*					Point_const_iterator;

	template < typename TDS2 >
	struct Rebind_TDS {
		typedef typename Cb::template Rebind_TDS<TDS2>::Other	Cb2;
		typedef Conforming_triangulation_cell_base_3<Gt, Cb2>	Other;
	}; // Rebind_TDS
	
protected:
	// The edge states.
	bool EF[6];
	bool mm[6];

public:
	Conforming_triangulation_cell_base_3(): Cb() {clear_conforming();}
	Conforming_triangulation_cell_base_3(const Vertex_handle& v0, const Vertex_handle& v1,
										  const Vertex_handle& v2, const Vertex_handle& v3)
										  : Cb(v0, v1, v2, v3) {clear_conforming();}
	Conforming_triangulation_cell_base_3(const Vertex_handle& v0, const Vertex_handle& v1,
										  const Vertex_handle& v2, const Vertex_handle& v3,
										  const Cell_handle&   n0, const Cell_handle&   n1,
										  const Cell_handle&   n2, const Cell_handle&   n3)
										  : Cb(v0, v1, v2, v3, n0, n1, n2, n3) {clear_conforming();}

	bool is_conforming(int i, int j) const {return EF[getEdgeIndex(i,j)];}
	bool is_marked(int i, int j) const {return mm[getEdgeIndex(i,j)];}

	void set_conforming(bool c0, bool c1, bool c2, bool c3, bool c4, bool c5) {
		EF[0] = c0;
		EF[1] = c1;
		EF[2] = c2;
		EF[3] = c3;
		EF[4] = c4;
		EF[5] = c5;
		clear_marked();
	}
	void set_conforming(int i, int j, bool c) {
		int index = getEdgeIndex(i,j);
		EF[index] = c;
		mm[index] = false;
	}
	void mark(int i, int j) {
		int index = getEdgeIndex(i,j);
		EF[index] = true;
		mm[index] = true;
	}

	void clear_conforming() {set_conforming(false, false, false, false, false, false);}
	void clear_marked() {
		mm[0] = false;
		mm[1] = false;
		mm[2] = false;
		mm[3] = false;
		mm[4] = false;
		mm[5] = false;
	}
	virtual void clear() {clear_conforming();}

	bool has_conforming() const {return EF[0] || EF[1] || EF[2] || EF[3] || EF[4] || EF[5];}

	virtual void reorient() {
		bool tmp = EF[1];
		EF[1] = EF[3];
		EF[3] = tmp;
		tmp = mm[1];
		mm[1] = mm[3];
		mm[3] = tmp;
	}

	virtual std::istream& read_cell(std::istream& is);
	virtual std::ostream& write_cell(std::ostream& os) const;

private:
	int getEdgeIndex(int i, int j) const {
		CGAL_triangulation_precondition(i>=0 && i<=3);
		CGAL_triangulation_precondition(j>=0 && j<=3);
		CGAL_triangulation_precondition(i != j);
		return (i==0 || j==0) ? i+j-1 : i+j;
	}
}; // Conforming_triangulation_cell_base_3

template < class Gt, class Cb >
inline std::istream& operator>>(std::istream& is, Conforming_triangulation_cell_base_3<Gt, Cb>& c) {
	is >> static_cast<Cb&>(c);
	return c.read_cell(is);
}

template < class Gt, class Cb >
inline std::ostream& operator<<(std::ostream &os, const Conforming_triangulation_cell_base_3<Gt, Cb>& c) {
	os << static_cast<const Cb&>(c);
	return c.write_cell(os);
}

template < class Gt, class Cb >
std::istream& Conforming_triangulation_cell_base_3<Gt, Cb>::read_cell(std::istream& is) {
	char s;
	for (int i = 0; i < 4; ++i) {
		for (int j = i+1; j < 4; ++j) {
			if (is_ascii(is))
				is >> s;
			else
				read(is, s);
			if (s == 'C')
				set_conforming(i, j, true);
		}
	}

	return is;
}

template < class Gt, class Cb >
std::ostream& Conforming_triangulation_cell_base_3<Gt, Cb>::write_cell(std::ostream& os) const {
	for (int i = 0; i < 4; ++i) {
		for (int j = i+1; j < 4; ++j) {
			if (is_conforming(i, j)) {
				os << "C";
			}
			else 
				os << "N";
			if (is_ascii(os))
				os << ' ';
		}
	}

	return os;
}

CGAL_END_NAMESPACE

#endif // CGAL_CONFORMING_TRIANGULATION_CELL_BASE_3_H