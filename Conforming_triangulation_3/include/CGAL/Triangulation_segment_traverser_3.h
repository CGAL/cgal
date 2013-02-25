//A class that follows a straight line through a Delaunay triangulation structure.
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

// A class that traverses the cells of a triangulation by following a segment from source to target.

#ifndef CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_H
#define CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_H

#include "Delaunay_triangulation_utils_3.h"

#include <CGAL/Random.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds>
class Delaunay_triangulation_utils_3;

template < class Gt, class Tds >
class Triangulation_segment_traverser_3: public Forward_circulator_base< typename Tds::Cell, std::ptrdiff_t, std::size_t > {
	typedef Triangulation_segment_traverser_3<Gt,Tds>	Self;

public:
	typedef Delaunay_triangulation_utils_3<Gt,Tds>		Tr;

	typedef typename Tds::Vertex						Vertex;
	typedef typename Tds::Cell							Cell;
	typedef typename Tds::Edge							Edge;
	typedef typename Tds::Facet							Facet;
	typedef typename Tds::Vertex_handle					Vertex_handle;
	typedef typename Tds::Cell_handle					Cell_handle;
	typedef typename Tds::Cell_iterator					Cell_iterator;
	typedef typename Tds::Cell_circulator				Cell_circulator;
	typedef typename Tds::Facet_iterator				Facet_iterator;
	typedef typename Tds::Facet_circulator				Facet_circulator;

	typedef typename Gt::Point_3						Point;
	typedef typename Tr::Locate_type					Locate_type;

protected:
	// The triangulation being traversed.
	const Tr* _tr;

	// The source and target points of the traversal.
	Point _source, _target;

	// The cell currently traversed and the previous one.
	Cell_handle _pos, _prev;

	// Otherwise, they indicate the simplex through which this cell was entered,
	// or the location of the source if it is in this cell.
	int _li, _lj;
	Locate_type _lt;

	// This bit signifies when a cell containing the target is found.
	bool _done;

	// Where possible, facets are checked in a random order.
	mutable Random rng;
            
public:
	Triangulation_segment_traverser_3(): _tr(NULL), _pos(), _prev(), _li(-1), _lj(-1), _lt(Tr::CELL), _done(true) {}
	Triangulation_segment_traverser_3(Vertex_handle s, const Point& t, const Tr* tr);
	Triangulation_segment_traverser_3(const Point& s, const Point& t, const Tr* tr, Cell_handle hint = Cell_handle());

	Self&			operator++();
	Self			operator++(int);
	Cell*			operator->()										{return &*_pos;}
	Cell&			operator*()											{return *_pos;}
	Cell_handle		handle()											{return _pos;}
	operator const	Cell_handle() const									{return _pos;}
	bool			operator==(const Self& ct) const;
	bool			operator!=(const Self& ct) const;

	bool			operator==(const Cell_handle& ch) const				{return ch == _pos;}
	bool			operator!=(const Cell_handle& ch) const				{return ch != _pos;}

	bool			operator==(Nullptr_t CGAL_triangulation_assertion_code(n)) const;
	bool			operator!=(Nullptr_t n) const;

	// Check if a cell containing the target has been reached.
	bool			has_next() const {return !_done;}

	// Access the current and previous cells.
	Cell_handle		previous() const {return _prev;}

	// Traverse to a cell containing the target.
	Cell_handle		traverse();

	// Get the type of simplex traversed last.
	void			traversed(Locate_type& lt, int& li, int& lj) const {lt = _lt; li = _li; lj = _lj;}

protected:
	// Traverse to the next cell along the line.
	virtual void	increment();

private:
	inline int edgeIndex(int i, int j) const {
		CGAL_triangulation_precondition(i>=0 && i<=3);
		CGAL_triangulation_precondition(j>=0 && j<=3);
		CGAL_triangulation_precondition(i != j);
		return (i==0 || j==0) ? i+j-1 : i+j;
	}

	inline Orientation coplanar_orientation(const Point& p, const Point& q, const Point& r, const Point& s) {
		return Gt().coplanar_orientation_3_object()(p, q, r, s);
	}
}; // class Triangulation_segment_traverser_3

template < class Gt, class Tds >
inline bool operator==(typename Tds::Cell_handle ch, Triangulation_segment_traverser_3<Gt, Tds> tci) {return tci == ch;}

template < class Gt, class Tds >
inline bool operator!=(typename Tds::Cell_handle ch, Triangulation_segment_traverser_3<Gt, Tds> tci) {return tci != ch;}

template < class Gt, class Tds >
Triangulation_segment_traverser_3<Gt,Tds>::Triangulation_segment_traverser_3(Vertex_handle s, const Point& t, const Tr* tr)
: _tr(tr), _pos(), _prev(), _lj(-1), _lt(Tr::VERTEX), _done(false) {
	CGAL_triangulation_precondition(!_tr->is_infinite(s));
	CGAL_triangulation_precondition(s->point() != t);
	CGAL_triangulation_precondition(_tr->dimension() >= 2);
	CGAL_triangulation_assertion(_tr->dimension() == 3 || _tr->plane(*_tr->finite_facets_begin()).has_on(_target));

	_source = s->point();
	_target = t;

	_pos = s->cell();
	int inf;
	if (_pos->has_vertex(_tr->infinite_vertex(), inf))
		_pos = _pos->neighbor(inf);
	_li = _pos->index(s);
	
	CGAL_triangulation_postcondition(_pos != Cell_handle());
}

template < class Gt, class Tds >
Triangulation_segment_traverser_3<Gt,Tds>::Triangulation_segment_traverser_3(const Point& s, const Point& t, const Tr* tr, Cell_handle hint)
: _tr(tr), _pos(), _prev(), _li(-1), _lj(-1), _done(false) {
	CGAL_triangulation_precondition(s != t);
	CGAL_triangulation_precondition(_tr->dimension() >= 2);
	CGAL_triangulation_assertion(_tr->dimension() == 3 || _tr->plane(*_tr->finite_facets_begin()).has_on(_target));

	_source = s;
	_target = t;

	_pos = _tr->locate(s, _lt, _li, _lj, hint);

	CGAL_triangulation_postcondition(_pos != Cell_handle());
}

template < class Gt, class Tds >
inline Triangulation_segment_traverser_3<Gt,Tds>&
Triangulation_segment_traverser_3<Gt,Tds>::operator++() {
	CGAL_triangulation_precondition(_pos != Cell_handle());
	increment();
	return *this;
}

template < class Gt, class Tds >
inline Triangulation_segment_traverser_3<Gt,Tds>
Triangulation_segment_traverser_3<Gt,Tds>::operator++(int) {
	Self tmp(*this);
	++(*this);
	return tmp;
}

template < class Gt, class Tds >
inline bool Triangulation_segment_traverser_3<Gt,Tds>::operator==(const Self& ct) const {
	CGAL_triangulation_precondition(_pos != Cell_handle() && ct._pos != Cell_handle());
	return (_pos == ct._pos &&
			_prev == ct._prev &&
			_tr == ct._tr &&
			_li == ct._li &&
			_lj == ct._lj &&
			_lt == ct._lt &&
			_source == ct._source &&
			_target == ct._target &&
			_done == ct._done);
}

template < class Gt, class Tds >
inline bool Triangulation_segment_traverser_3<Gt,Tds>::operator!=(const Self& ct) const {
	return !(*this == ct);
}

template < class Gt, class Tds >
inline bool Triangulation_segment_traverser_3<Gt,Tds>::operator==(Nullptr_t CGAL_triangulation_assertion_code(n)) const {
	CGAL_triangulation_assertion(n == NULL);
	return _pos == Cell_handle();
}

template < class Gt, class Tds >
inline bool Triangulation_segment_traverser_3<Gt,Tds>::operator!=(Nullptr_t n) const {
	CGAL_triangulation_assertion(n == NULL);
	return !(*this == n);
}

template < class Gt, class Tds >
inline typename Triangulation_segment_traverser_3<Gt,Tds>::Cell_handle
Triangulation_segment_traverser_3<Gt,Tds>::traverse() {
	while (!_done)
		increment();
	return _pos;
}

template < class Gt, class Tds >
void Triangulation_segment_traverser_3<Gt,Tds>::increment() {
	CGAL_triangulation_precondition(!_done);

	// Walks to the next cell over a facet intersected by the line from source to target.
	// This method is based on Triangulation_3::locate().
	switch (_tr->dimension()) {
		case 3: {
			// Infinite cells should be handled differently.
			int inf;
			if (_pos->has_vertex(_tr->infinite_vertex(), inf)) {
				// If this cell was reached by traversal from a finite one, it must be the final cell.
				Cell_handle fin = _pos->neighbor(inf);
				if (fin == _prev) {
					_done = true;
					// Here, I do not change the _lt etc. because they were already set.
					return;
				}

				const Point* pts[4];
				for (int i = 0; i != 4; ++i)
					if (i != inf)
						pts[i] = &(_pos->vertex(i)->point());
				pts[inf] = &_target;
				Orientation o[4];

				// For the remembering stochastic walk, we start trying with a random index:
				int li = rng.template get_bits<2>();

				// Check if the target lies outside the convex hull.
				if (orientation(*pts[0], *pts[1], *pts[2], *pts[3]) == POSITIVE) {
					// The target lies in an infinite cell.
					_done = true;
					return;

					// Note that we do not traverse to other infinite cells.
				}

				pts[inf] = &_source;
				CGAL_triangulation_assertion(orientation(*pts[0], *pts[1], *pts[2], *pts[3]) == POSITIVE);

				// Check if the line enters an adjacent infinite cell.
				// This occurs if the target lies on the other side of
				// a plane through one of the finite edges and the source point.
				for (int j = 0; j != 4; ++j, li = (li+1)&3) {
					if (li == inf) {
						o[li] = COPLANAR;
						continue;
					}

					Cell_handle next = _pos->neighbor(li);
					if (next == _prev) {
						o[li] = POSITIVE;
						continue;
					}

					const Point* backup = pts[li];
					pts[li] = &_target;
					o[li] = orientation(*pts[0], *pts[1], *pts[2], *pts[3]);

					if (o[li] != NEGATIVE) {
						pts[li] = backup;
						continue;
					}

					// The target lies behind the plane through the source and two finite vertices.
					// Traverse to the incident infinite cell.
					_prev = _pos;
					_pos = next;
					CGAL_triangulation_assertion(_tr->is_infinite(_pos));

					_lt = Tr::FACET;
					_li = _pos->index(_prev);
					return;
				}

				// The line enters the convex hull here (or lies on the finite facet).
				_prev = _pos;
				_pos = fin;

				// Check through which simplex the line traverses.
				switch (o[0]+o[1]+o[2]+o[3]) {
					case 3:
						_lt = Tr::FACET;
						_li = _pos->index(_prev);
						return;
					case 2:
						_lt = Tr::EDGE;
						for (int i = 0; i < 4; ++i) {
							if (o[i] == COPLANAR && i != inf) {
								Edge opp = _tr->opposite_edge(_prev, inf, i);
								_li = _pos->index(_prev->vertex(opp.second));
								_lj = _pos->index(_prev->vertex(opp.third));
								return;
							}
						}
						CGAL_triangulation_assertion(false);
						return;
					case 1:
						_lt = Tr::VERTEX;
						for (int i = 0; i < 4; ++i) {
							if (o[i] == POSITIVE) {
								_li = _pos->index(_prev->vertex(i));
								return;
							}
						}
						CGAL_triangulation_assertion(false);
						return;
					default:
						CGAL_triangulation_assertion(false);
						return;
				}
			}

			const Point* pts[4] = {&(_pos->vertex(0)->point()),
								   &(_pos->vertex(1)->point()),
								   &(_pos->vertex(2)->point()),
								   &(_pos->vertex(3)->point())};
			
			// We check in which direction the target lies
			// by comparing its position relative to the planes through the
			// source and the edges of the cell.
			Orientation o[6];
			bool calc[6] = {false, false, false, false, false, false};

			if (_lt == Tr::VERTEX) {
				// The three planes through the vertex are set to coplanar.
				for (int j = 0; j < 4; ++j) {
					if (_li != j) {
						int ij = edgeIndex(_li, j);
						o[ij] = COPLANAR;
						calc[ij] = true;
					}
				}
			}
			else if (_lt == Tr::EDGE) {
				// The plane through the edge is set to coplanar.
				int ij = edgeIndex(_li, _lj);
				o[ij] = COPLANAR;
				calc[ij] = true;
			}

			// For the remembering stochastic walk, we start trying with a random facet:
			int li = rng.template get_bits<2>();

			CGAL_triangulation_assertion_code(bool incell = true;)
			for (int k = 0; k < 4; ++k, li = (li+1)&3) {
				Cell_handle next = _pos->neighbor(li);
				if (next == _prev)
					continue;

				// Check if the target is outside the cell.
				const Point* backup = pts[li];
				pts[li] = &_target;
				if (orientation(*pts[0], *pts[1], *pts[2], *pts[3]) != NEGATIVE) {
					pts[li] = backup;
					continue;
				}

				CGAL_triangulation_assertion_code(incell = false;)

				// Check if the target is inside the pyramid.
				int lj = rng.template get_bits<2>();
				int or = 0;
				for (int l = 0; l < 4; ++l, lj = (lj+1)&3) {
					if (li == lj)
						continue;

					// We check the orientation of the target compared to the plane
					// Through the source and the edge opposite of ij.
					int oij = 5 - edgeIndex(li, lj);
					if (!calc[oij]) {
						const Point* backup2 = pts[lj];
						pts[lj] = &_source;
						o[oij] = orientation(*pts[0], *pts[1], *pts[2], *pts[3]);
						pts[lj] = backup2;
						calc[oij] = true;
					}

					or -= o[oij];
					if (o[oij] == POSITIVE) {
						// The target is not inside the pyramid.
						// Invert the planes.
						// This can be safely done because either
						// they were not calculated yet,
						// or they will no longer be used.
						for (int j = 0; j < 4; ++j) {
							if (li == j) continue;
							int oij = 5 - edgeIndex(li, j);
							o[oij] = -o[oij];
						}
						or = 0;
						break;
					}
				}

				if (or == 0) {
					// Either the target is not inside the pyramid,
					// or the pyramid is degenerate.
					pts[li] = backup;
					continue;
				}

				// The target is inside the pyramid.
				_prev = _pos;
				_pos = next;

				switch (or) {
					case 3:
						_lt = Tr::FACET;
						_li = _pos->index(_prev);
						return;
					case 2:
						_lt = Tr::EDGE;
						for (int j = 0; j < 4; ++j) {
							if (li != j && o[5-edgeIndex(li, j)] == COPLANAR) {
								Edge opp = _tr->opposite_edge(_prev, li, j);
								_li = _pos->index(_prev->vertex(opp.second));
								_lj = _pos->index(_prev->vertex(opp.third));
								return;
							}
						}
						CGAL_triangulation_assertion(false);
						return;
					case 1:
						_lt = Tr::VERTEX;
						for (int j = 0; j < 4; ++j) {
							if (li != j && o[5-edgeIndex(li, j)] == NEGATIVE) {
								_li = _pos->index(_prev->vertex(j));
								return;
							}
						}
						CGAL_triangulation_assertion(false);
						return;
					default:
						CGAL_triangulation_assertion(false);
						return;
				}
			}

			// The target lies inside this cell.
			CGAL_triangulation_assertion(incell);
			_done = true;
			return;
		}
		case 2: {
			int inf;
			if (_pos->has_vertex(_tr->infinite_vertex(), inf) && inf < 3) {
				// If this cell was reached by traversal from a finite one, it must be the final cell.
				Cell_handle fin = _pos->neighbor(inf);
				if (fin == _prev) {
					_done = true;
					return;
				}

				// Check the neighboring cells.
				if (coplanar_orientation(_source, _pos->vertex(_tr->ccw(inf))->point(), _pos->vertex(_tr->cw(inf))->point(), _target) == NEGATIVE) {
					Cell_handle tmp = _pos->neighbor(_tr->cw(inf));
					_prev = _pos;
					_pos = tmp;
					return;
				}
				if (coplanar_orientation(_source, _pos->vertex(_tr->cw(inf))->point(), _pos->vertex(_tr->ccw(inf))->point(), _target) == NEGATIVE) {
					Cell_handle tmp = _pos->neighbor(_tr->ccw(inf));
					_prev = _pos;
					_pos = tmp;
					return;
				}
				if (coplanar_orientation(_pos->vertex(_tr->ccw(inf))->point(), _pos->vertex(_tr->cw(inf))->point(), _source, _target) != POSITIVE) {
					_prev = _pos;
					_pos = fin;
					return;
				}

				// The target lies in this infinite cell.
				_done = true;
				return;
			}

			const Point* pts[3] = {&(_pos->vertex(0)->point()),
								   &(_pos->vertex(1)->point()),
								   &(_pos->vertex(2)->point())};
			
			switch (_lt) {
				case Tr::VERTEX: {
					// First we try the incident edges.
					Orientation ocw = coplanar_orientation(*pts[_li], *pts[(_li+2)%3], *pts[(_li+1)%3], _target);
					if (_pos->neighbor((_li+1)%3) != _prev &&
						ocw == NEGATIVE) {
							Cell_handle tmp = _pos->neighbor((_li+1)%3);
							_prev = _pos;
							_pos = tmp;
							_li = _pos->index(_prev->vertex(_li));
							return;
					}
					Orientation occw = coplanar_orientation(*pts[_li], *pts[(_li+1)%3], *pts[(_li+2)%3], _target);
					if (_pos->neighbor((_li+2)%3) != _prev &&
						occw == NEGATIVE) {
							Cell_handle tmp = _pos->neighbor((_li+2)%3);
							_prev = _pos;
							_pos = tmp;
							_li = _pos->index(_prev->vertex(_li));
							return;
					}

					// Then we try the opposite edge.
					if (coplanar_orientation(*pts[(_li+1)%3], *pts[(_li+2)%3], *pts[_li], _target) == NEGATIVE) {
						Cell_handle tmp = _pos->neighbor(_li);
						_prev = _pos;
						_pos = tmp;

						switch (ocw+occw) {
							case 2: {
								_lt = Tr::EDGE;
								_li = _pos->index(_prev->vertex((_li+1)%3));
								_lj = _pos->index(_prev->vertex((_li+2)%3));
								return;
							}
							case 1: {
								_lt = Tr::VERTEX;
								if (ocw == COLLINEAR) _li = _pos->index(_prev->vertex((_li+2)%3));
								else _li = _pos->index(_prev->vertex((_li+1)%3));
								return;
							}
							default:
								CGAL_triangulation_assertion(false);
								return;
						}
					}

					// The target lies in this cell.
					_done = true;
					return;
				}
				case Tr::EDGE: {
					int lk = 3-_li-_lj;

					if (_pos->neighbor(lk) != _prev) {
						// Check the edge itself
						switch (coplanar_orientation(*pts[_li], *pts[_lj], *pts[lk], _target)) {
							case COLLINEAR: {
								// The target lies in this cell.
								_done = true;
								return;
							}
							case NEGATIVE: {
								// The target lies opposite of the edge.
								Cell_handle tmp = _pos->neighbor(lk);
								_prev = _pos;
								_pos = tmp;
								_li = _pos->index(_prev->vertex(_li));
								_lj = _pos->index(_prev->vertex(_lj));
								return;
							}
							default:
								break;
						}
					}

					Orientation o = coplanar_orientation(_source, *pts[lk], *pts[_li], _target);
					switch (o) {
						case POSITIVE: {
							// The ray passes through the edge ik.
							if (coplanar_orientation(*pts[lk], *pts[_li], _source, _target) == NEGATIVE) {
								Cell_handle tmp = _pos->neighbor(_lj);
								_prev = _pos;
								_pos = tmp;

								if (collinear(_source, *pts[_li], _target)) {
									_lt = Tr::VERTEX;
									_li = _pos->index(_prev->vertex(_li));
								}
								else {
									_lt = Tr::EDGE;
									_li = _pos->index(_prev->vertex(_li));
									_lj = _pos->index(_prev->vertex(lk));
								}
								return;
							}
							break;
						}
						default: {
							// The ray passes through the edge jk.
							if (coplanar_orientation(*pts[lk], *pts[_lj], _source, _target) == NEGATIVE) {
								Cell_handle tmp = _pos->neighbor(_li);
								_prev = _pos;
								_pos = tmp;

								if (collinear(_source, *pts[_lj], _target)) {
									_lt = Tr::VERTEX;
									_li = _pos->index(_prev->vertex(_lj));
								}
								else if (o == COLLINEAR) {
									_lt = Tr::VERTEX;
									_li = _pos->index(_prev->vertex(lk));
								}
								else {
									_lt = Tr::EDGE;
									_li = _pos->index(_prev->vertex(lk));
									_lj = _pos->index(_prev->vertex(_lj));
								}
								return;
							}
							break;
						}
					}

					// The target lies in this cell.
					_done = true;
					return;
				}
				case Tr::FACET: {
					// We test its edges in a random order until we find a neighbor to go further
					int li = rng.get_int(0, 3);

					Orientation o[3];
					bool calc[3] = {false, false, false};

					for (int j = 0; j != 3; ++j, li = (li+1)%3) {
						Cell_handle next = _pos->neighbor(li);
						if (next == _prev)
							continue;

						// The target should lie on the other side of the edge.
						if (coplanar_orientation(*pts[(li+1)%3], *pts[(li+2)%3], *pts[li], _target) != NEGATIVE)
							continue;

						// The target should lie inside the wedge.
						if (!calc[(li+1)%3]) {
							o[(li+1)%3] = coplanar_orientation(_source, *pts[(li+1)%3], *pts[(li+2)%3], _target);
							calc[(li+1)%3] = true;
						}
						if (o[(li+1)%3] == NEGATIVE) continue;

						if (!calc[(li+2)%3]) {
							o[(li+2)%3] = coplanar_orientation(_source, *pts[(li+2)%3], *pts[(li+1)%3], _target);
							calc[(li+2)%3] = true;
						}
						if (o[(li+2)%3] == POSITIVE) continue;

						_prev = _pos;
						_pos = next;

						switch (o[(li+1)%3]+o[(li+2)%3]) {
							case 2: {
								_lt = Tr::EDGE;
								_li = _pos->index(_prev->vertex((li+1)%3));
								_lj = _pos->index(_prev->vertex((li+2)%3));
								return;
							}
							case 1: {
								_lt = Tr::VERTEX;
								if (o[(li+1)%3] == COLLINEAR) _li = _pos->index(_prev->vertex((li+1)%3));
								else _li = _pos->index(_prev->vertex((li+2)%3));
								return;
							}
							default:
								CGAL_triangulation_assertion(false);
								return;
						}
					}

					// The target lies in this cell.
					_done = true;
					return;
				}
				default:
					CGAL_triangulation_assertion(false);
			}
		}
	}
}

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_H