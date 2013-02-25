//A traverser that collects points encroaching upon its path.
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

// A class that traverses the cells of a triangulation by following a segment from source to target,
// while collecting points that encroach upon the segment.

#ifndef CGAL_ENCROACHING_COLLECTER_3_H
#define CGAL_ENCROACHING_COLLECTER_3_H

#include "Triangulation_segment_traverser_3.h"
//#include "Conforming_Delaunay_triangulation_3.h"

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds, class Itag >
class Conforming_Delaunay_triangulation_3;

template < class Gt, class Tds, class Itag = No_intersection_tag, class Out = std::back_insert_iterator<std::list<Point_3<Gt>>>>
class Encroaching_collecter_3: public Triangulation_segment_traverser_3<Gt,Tds> {
	typedef Encroaching_collecter_3<Gt,Tds,Itag,Out>			Self;
	typedef Triangulation_segment_traverser_3<Gt,Tds>			Base;

	typedef typename Gt::Plane_3								Plane;

public:
	typedef Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>	CDT;

	typedef typename Tds::Vertex								Vertex;
	typedef typename Tds::Cell									Cell;
	typedef typename Tds::Edge									Edge;
	typedef typename Tds::Facet									Facet;
	typedef typename Tds::Vertex_handle							Vertex_handle;
	typedef typename Tds::Cell_handle							Cell_handle;
	typedef typename Tds::Cell_iterator							Cell_iterator;
	typedef typename Tds::Cell_circulator						Cell_circulator;
	typedef typename Tds::Facet_iterator						Facet_iterator;
	typedef typename Tds::Facet_circulator						Facet_circulator;

	typedef typename Gt::Point_3								Point;
	typedef typename CDT::Locate_type							Locate_type;
	typedef typename CDT::Bi_vertex								Bi_vertex;

protected:
	Out	_out;

	Bi_vertex se; // These are used when the source lies on an edge.
            
public:
	Encroaching_collecter_3(): Base() {}
	Encroaching_collecter_3(Vertex_handle v, const Point& t, const CDT* tr, Out output)
		: Base(v, t, tr), _out(output) {if (_lt == CDT::EDGE) se = _tr->sort_bi_vertex(_pos, _li, _lj);}
	Encroaching_collecter_3(const Point& s, const Point& t, const CDT* tr, Out output, Cell_handle hint = Cell_handle())
		: Base(s, t, tr, hint), _out(output) {if (_lt == CDT::EDGE) se = _tr->sort_bi_vertex(_pos, _li, _lj);}

	virtual bool	barrier_hit() const {return _lt == CDT::EDGE && _pos->is_conforming(_li, _lj) && _tr->sort_bi_vertex(_pos, _li, _lj) != se;}

protected:
	// Traverse to the next cell along the line.
	virtual void	increment();
};

template < class Gt, class Tds, class Itag, class Out >
void Encroaching_collecter_3<Gt,Tds,Itag,Out>::increment() {
	// Walk to the next cell.
	Base::increment();
	
	// Check for encroaching points.
	// Which points are checked is based on the type of simplex traversed.
	const CDT* _cdt = dynamic_cast<const CDT*>(_tr);
	switch (_lt) {
		case CDT::FACET:
			for (int i = 0; i < 4; ++i)
				if (i != _li && !_cdt->is_infinite(_pos->vertex(i)) && _cdt->is_encroaching(_source, _target, _pos->vertex(i)->point()))
					*_out++ = _pos->vertex(i)->point();
			break;
		case CDT::EDGE: {
			Vertex_handle vi = _pos->vertex(_li), vj = _pos->vertex(_lj), vk;
			Cell_handle c;

			// Check the vertices in the star around the edge.
			Facet_circulator fit = _cdt->incident_facets(_pos, _li, _lj), start(fit);
			do {
				c = fit->first;
				vk = c->vertex(6 - c->index(vi) - c->index(vj) - fit->second);
				if (!_cdt->is_infinite(vk) && !collinear(_source, _target, vk->point()) && _cdt->is_encroaching(_source, _target, vk->point()))
					*_out++ = vk->point();
				++fit;
			} while (fit != start);

			// Check the vertices of the edge.
			if (_cdt->is_encroaching(_source, _target, vi->point()))
				*_out++ = vi->point();
			if (_cdt->is_encroaching(_source, _target, vj->point()))
				*_out++ = vj->point();
			break;
		}
		case CDT::VERTEX:
			CGAL_triangulation_assertion(collinear(_source, _target, _pos->vertex(_li)->point()));
			break;
		default:
			CGAL_triangulation_assertion(false);
			break;
	}
}

CGAL_END_NAMESPACE

#endif // CGAL_ENCROACHING_COLLECTER_3_H