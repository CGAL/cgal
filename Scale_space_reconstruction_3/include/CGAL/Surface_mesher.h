//The ball pivoting algorithm implemented using the alpha-shape.
//Copyright (C) 2013  INRIA - Sophia Antipolis
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
// Author(s):      Thijs van Lankveld


#ifndef BALL_PIVOTING
#define BALL_PIVOTING

#include <CGAL/Delaunay_triangulation_3.h>

#define BALL_PIVOTING_FIXED
#define BALL_PIVOTING_CONNECTED

#ifdef BALL_PIVOTING_FIXED
#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>
#else
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#endif

#include "Flagged_triangulation_cell_base_3.h"

// Collect the vertices within a ball of fixed radius.
template < typename Kernel,
		   typename Vb = CGAL::Triangulation_vertex_base_3<Kernel>,
		   typename Cb = CGAL::Triangulation_cell_base_3<Kernel> >
class Surface_mesher {
	typedef typename Kernel::FT											Scalar;

	// Convert the DT types to be in an alpha-shape.
	// Vertex.
#ifdef BALL_PIVOTING_FIXED
	typedef CGAL::Fixed_alpha_shape_vertex_base_3<Kernel, Vb>			AVb;
#else
	typedef CGAL::Alpha_shape_vertex_base_3<Kernel, Vb>					AVb;
#endif

	// Cell.
	typedef CGAL::Flagged_triangulation_cell_base_3<Kernel, Cb>			FCb;
#ifdef BALL_PIVOTING_FIXED
	typedef CGAL::Fixed_alpha_shape_cell_base_3<Kernel, FCb>			ACb;
#else
	typedef CGAL::Alpha_shape_cell_base_3<Kernel, FCb>					ACb;
#endif

	// Triangulation.
	typedef CGAL::Triangulation_data_structure_3<AVb, ACb>				Tds;
	typedef CGAL::Delaunay_triangulation_3<Kernel, Tds>					ADT;
#ifdef BALL_PIVOTING_FIXED
	typedef CGAL::Fixed_alpha_shape_3<ADT>								AS;
#else
	typedef CGAL::Alpha_shape_3<ADT>									AS;
#endif

public:
	typedef ADT															Triangulation;
	typedef AS															Alpha_shape;
	
	typedef typename AS::Vertex_handle									Vertex_handle;
	typedef typename AS::Cell_handle									Cell_handle;
	typedef typename AS::Facet											Facet;

private:
	typedef typename AS::All_cells_iterator								All_cells_iterator;
	typedef typename AS::Finite_facets_iterator							Finite_facets_iterator;
	typedef typename AS::Classification_type							Classification_type;
	
	AS* _shape;

	Scalar radius2;
	
	// Once a facet is outputted, it is marked as handled.
	inline bool isHandled(Cell_handle c, unsigned int li) const {
		switch (li) {
		case 0: return (c->flag()&1) != 0;
		case 1: return (c->flag()&2) != 0;
		case 2: return (c->flag()&4) != 0;
		case 3: return (c->flag()&8) != 0;
		}
		return false;
	}
	inline bool isHandled(const Facet& f) const {return isHandled(f.first, f.second);}

	inline void setHandled(Cell_handle c, unsigned int li) {
		switch (li) {
		case 0: c->flag() |= 1; return;
		case 1: c->flag() |= 2; return;
		case 2: c->flag() |= 4; return;
		case 3: c->flag() |= 8; return;
		}
	}
	inline void setHandled(Facet f) {setHandled(f.first, f.second);}

public:
	typedef typename Kernel::Point_3									Point;
	typedef typename Kernel::Triangle_3									Triangle;

	Surface_mesher(const Scalar& r2): _shape(0), radius2(r2) {}
	Surface_mesher(const Triangulation& tr, const Scalar& r2): radius2(r2) {
#ifdef BALL_PIVOTING_FIXED
		_shape = new AS(tr, radius2);
#else
		_shape = new AS(tr, radius2, AS::GENERAL);
#endif
	}
	~Surface_mesher() {if (_shape != 0) {delete _shape; _shape = 0;}}

	template < class InputIterator >
	void construct_shape(InputIterator start, InputIterator end) {
		std::cout << "BP: Construct shape";
		ADT tr(start, end);

		if (_shape != 0)
			delete _shape;
#ifdef BALL_PIVOTING_FIXED
		_shape = new AS(tr, radius2);
#else
		_shape = new AS(tr, radius2, AS::GENERAL);
#endif
		std::cout << ": " << _shape->number_of_vertices() << " points" << std::endl;
	}
	bool is_constructed() const {return _shape != 0;}
	void clear() {if (_shape != 0) {delete _shape; _shape = 0;}}

	AS* shape() {return _shape;}

	Scalar get_radius2() const {return radius2;}
	void set_radius2(const Scalar& r2) {
		radius2 = r2;
		if (_shape != 0) {
			ADT tr;
			tr.swap(*_shape);

			delete _shape;
#ifdef BALL_PIVOTING_FIXED
			_shape = new AS(tr, radius2);
#else
			_shape = new AS(tr, radius2, AS::GENERAL);
#endif
		}
	}

	template < class OutputIterator >
	OutputIterator collect_shell(Cell_handle c, unsigned int li, OutputIterator out) {
		// Collect one surface mesh from the alpha-shape in a fashion similar to ball-pivoting.
		// Invariant: the facet is regular or singular.

		// To stop stack overflows: use own stack.
		std::stack<Facet> stack;
		stack.push(Facet(c, li));

		Facet f;
		Cell_handle n, p; int ni, pi;
		Vertex_handle a;
		Classification_type cl;
		while (!stack.empty()) {
			f = stack.top();
			stack.pop();

			// Check if the cell was already handled.
			// Note that this is an extra check that in many cases is not necessary.
			if (isHandled(f))
				continue;

			// The facet is part of the surface.
			CGAL_triangulation_assertion(!_shape->is_infinite(f));
			*out++ = f;
			setHandled(f);
		
#ifdef BALL_PIVOTING_CONNECTED
			// Pivot over each of the facet's edges and continue the surface at the next regular or singular facet.
			for (unsigned int i = 0; i < 4; ++i) {
				// Skip the current facet.
				if (i == f.second || isHandled(f.first, i))
					continue;

				// Rotate around the edge (starting from the shared facet in the current cell) until a regular or singular facet is found.
				n = f.first;
				ni = i;
				a = f.first->vertex(f.second);
				cl = _shape->classify(Facet(n, ni));
				while (cl != AS::REGULAR && cl != AS::SINGULAR) {
					p = n;
					n = n->neighbor(ni);
					ni = n->index(a);
					pi = n->index(p);
					a = n->vertex(pi);
					cl = _shape->classify(Facet(n, ni));
				}

				// Continue the surface at the next regular or singular facet.
				stack.push(Facet(n, ni));
			}
#endif
		}

		return out;
	}

	template < class OutputIterator >
	OutputIterator collect_shell(const Facet& f, OutputIterator out) {
		return collect_shell(f.first, f.second, out);
	}

	template < class OutputIterator >
	OutputIterator collect_surface(OutputIterator out) {
		std::cout << "BP: Collect surfaces" << std::endl;
		// Collect all surface meshes from the alpha-shape in a fashion similar to ball-pivoting.
		// Reset the facet handled markers.
		for (All_cells_iterator cit = _shape->all_cells_begin(); cit != _shape->all_cells_end(); ++cit)
			cit->flag() = 0;

		// We check each of the facets: if it is not handled and either regular or singular,
		// we start collecting the next surface from here.
		Facet m;
		int ns = 0;
		for (Finite_facets_iterator fit = _shape->finite_facets_begin(); fit != _shape->finite_facets_end(); ++fit) {
			m = _shape->mirror_facet(*fit);
			switch (_shape->classify(*fit)) {
			case AS::REGULAR:
				if (!isHandled(*fit) && !isHandled(m))
					++ns;
				// Build a surface from the outer cell.
				if (_shape->classify(fit->first) == AS::EXTERIOR)
					collect_shell(*fit, out);
				else
					collect_shell(m, out);
				break;
			case AS::SINGULAR:
				if (!isHandled(*fit))
					++ns;
				// Build a surface from both incident cells.
				collect_shell(*fit, out);
				if (!isHandled(m))
					++ns;
				collect_shell(m, out);
				break;
			}
		}

		std::cout << ": " << ns << std::flush;

		return out;
	}

	template < class InputIterator, class OutputIterator >
	OutputIterator operator()(InputIterator start, InputIterator end, OutputIterator out) {
		construct_shape(start, end);
		return collect_surface(out);
	}

	bool locate_vertex(const Point& p, Vertex_handle& v, Cell_handle hint = Cell_handle()) const {
		typename AS::Locate_type lt; int li, lj;
		Cell_handle c = _shape->locate(p, lt, li, lj, hint);
		if (lt != AS::VERTEX)
			return false;
		v = c->vertex(li);
		return true;
	}

	void ordered_vertices(const Facet& f, Vertex_handle& v0, Vertex_handle& v1, Vertex_handle& v2) {
		if ((f.second&1) == 0) {
			v0 = f.first->vertex((f.second+2)&3);
			v1 = f.first->vertex((f.second+1)&3);
			v2 = f.first->vertex((f.second+3)&3);
		}
		else {
			v0 = f.first->vertex((f.second+1)&3);
			v1 = f.first->vertex((f.second+2)&3);
			v2 = f.first->vertex((f.second+3)&3);
		}
	}

	Triangle oriented_triangle(Cell_handle c, unsigned int li) const {return _shape->triangle(c, li);}
	Triangle oriented_triangle(const Facet& f) const {return _shape->triangle(f);}
}; // class Surface_mesher

#endif // BALL_PIVOTING
