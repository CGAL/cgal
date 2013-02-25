//A vertex in the conforming Delaunay triangulation structure.
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

#ifndef CGAL_CONFORMING_TRIANGULATION_VERTEX_BASE_3_H
#define CGAL_CONFORMING_TRIANGULATION_VERTEX_BASE_3_H

#include <CGAL/Triangulation_vertex_base_3.h>

CGAL_BEGIN_NAMESPACE

template < typename GT, typename Vb = Triangulation_vertex_base_3<GT> >
class Conforming_triangulation_vertex_base_3: public Vb {
	bool _s;

public:
	typedef typename Vb::Cell_handle							Cell_handle;
	typedef typename Vb::Point									Point;

	template < typename TDS2 > struct Rebind_TDS {
		typedef typename Vb::template Rebind_TDS<TDS2>::Other	Vb2;
		typedef Conforming_triangulation_vertex_base_3<GT, Vb2>	Other;
	};

	Conforming_triangulation_vertex_base_3(): Vb(), _s(false) {}
	Conforming_triangulation_vertex_base_3(const Point& p): Vb(p), _s(false) {}
	Conforming_triangulation_vertex_base_3(const Point& p, bool s): Vb(p), _s(s) {}
	Conforming_triangulation_vertex_base_3(const Point& p, Cell_handle c): Vb(p, c), _s(false) {}
	Conforming_triangulation_vertex_base_3(const Point& p, bool s, Cell_handle c): Vb(p, c), _s(s) {}
	Conforming_triangulation_vertex_base_3(Cell_handle c): Vb(c), _s(false) {}

	const bool& steiner() const {return _s;}
	bool& steiner() {return _s;}
};

template < class GT, class Vb >
std::istream&
operator>>(std::istream &is, Conforming_triangulation_vertex_base_3<GT, Vb> &v) {
	is >> static_cast<Vb&>(v);
	if (is_ascii(is)) is >> v.steiner();
	else read(is, v.steiner());
	return is;
}

template < class GT, class Vb >
std::ostream&
operator<<(std::ostream &os, const Conforming_triangulation_vertex_base_3<GT, Vb> &v) {
	os << static_cast<const Vb&>(v);
	if (is_ascii(os)) os << ' ';
	return os << v.steiner();
}

CGAL_END_NAMESPACE

#endif // CGAL_CONFORMING_TRIANGULATION_VERTEX_BASE_3_H
