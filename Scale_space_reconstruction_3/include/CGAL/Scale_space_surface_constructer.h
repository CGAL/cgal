//A method to construct a surface.
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


#ifndef SCALE_SPACE_SURFACE_CONSTRUCTER
#define SCALE_SPACE_SURFACE_CONSTRUCTER

#include <iostream>
#include <list>
#include <map>

#include "Mean_neighborhood_ball.h"
#include "Scale_space_transform.h"
#include "Surface_mesher.h"

// Smooth a set of points and then construct a surface mesh on the smoothed set;
// Finally transport the surface back to the original point locations.
template < typename Kernel,
		   typename Vb = CGAL::Triangulation_vertex_base_3<Kernel>,
		   typename Cb = CGAL::Triangulation_cell_base_3<Kernel> >
class Scale_space_surface_constructer {
	typedef typename Mean_neighborhood_ball<Kernel>			ComputeRadius;
	typedef typename Scale_space_transform<Kernel>			Scale_space_transform;
	typedef typename Surface_mesher<Kernel, Vb, Cb>			Surface_mesher;

public:
	typedef typename Kernel::FT								Scalar;

	typedef typename Kernel::Point_3						Point;
	typedef typename Kernel::Triangle_3						Triangle;
	
	typedef typename Surface_mesher::Vertex_handle			Vertex_handle;
	typedef typename Surface_mesher::Facet					Facet;

	typedef typename std::vector<Point>						PointCollection;

private:
	PointCollection _moved;

	ComputeRadius mnb;
	unsigned int iterations;
	Surface_mesher bp;

public:
	Scale_space_surface_constructer(unsigned int iter = 1): mnb(), iterations(iter), bp(0) {}

	void set_iterations(unsigned int iter) {iterations = iter;}

	const PointCollection& moved() const {return _moved;}

	// Input: iterators over Point, output: collection of Triangles.
	template < class InputIterator, class OutputIterator >
	OutputIterator operator()(InputIterator start, InputIterator end, OutputIterator out) {
		typedef std::map<Point, Point>						Map;
		typedef std::list<Facet>							List;

		// Compute the radius for which the mean ball would contain 30 points.
		Scalar radius2 = mnb(start, end);
		radius2 *= radius2;

		// Perturb the points.
		Scale_space_transform ssp(radius2);
		ssp(start, end, iterations);
		// Note that the transform has constructed a novel vector to hold the tranformed points.
		_moved = PointCollection(ssp.points().begin(), ssp.points().end());

		// Perform ball-pivoting on the perturbed points.
		std::list<Facet> facets;
		bp = Surface_mesher(radius2);
		bp(_moved.begin(), _moved.end(), std::back_inserter(facets));

		// Transport the points back to their original position.
		// Note that if you don't need the vertex base to have info,
		// you can replace the vertices by CGAL::Triangulation_vertex_base_with_info_3<size_t, Vb>
		// to reduce the computation time for retrieving the original point locations
		// (just set the info of each vertex to its index in the original collection).
		Map map;
		InputIterator oit = start;
		for (PointCollection::iterator pit = _moved.begin(); pit != _moved.end(); ++pit, ++oit)
			map[*pit] = *oit;

		Triangle t;
		for (List::const_iterator fit = facets.begin(); fit != facets.end(); ++fit) {
			t = bp.oriented_triangle(*fit);
			*out++ = Triangle(map[t.vertex(0)],
							  map[t.vertex(1)],
							  map[t.vertex(2)]);
		}

		return out;
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
}; // class Scale_space_surface_constructer

#endif // SCALE_SPACE_SURFACE_CONSTRUCTER