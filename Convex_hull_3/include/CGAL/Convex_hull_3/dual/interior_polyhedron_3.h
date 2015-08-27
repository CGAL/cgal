// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Jocelyn Meyron
//                 Pierre Alliez
//

#ifndef INTERIOR_POLYHEDRON_3_H
#define INTERIOR_POLYHEDRON_3_H

// LP solver to compute an interior point of a polyhedron
#include <CGAL/QP_functions.h>
#include <CGAL/QP_models.h>
#include <limits>
#include <CGAL/number_utils.h>
#include <CGAL/assertions.h>

// Description taken from http://www.qhull.org/html/qhalf.htm

// If you do not know an interior point for the halfspaces, use linear programming
// to find one. Assume, n halfspaces defined by: aj*x1+bj*x2+cj*x3+dj>=0, j=1..n.
// Perform the following linear program:
//		max(x5) aj*x1+bj*x2+cj*x3+dj*x4-x5>=0, j=1..n

// Then, if [x1,x2,x3,x4,x5] is an optimal m_solution with x4,x5>0 we get:
//		aj*(x1/x4)+bj*(x2/x4)+cj*(x3/x4)+dj>=(x5/x4)>0, j=1..n
// and conclude that the point [x1/x4,x2/x4,x3/x4] is in the interior of all
// the halfspaces. Note that x5 is optimal, so this point is "way in" the
// interior (good for precision errors).

// After finding an interior point, the rest of the intersection algorithm is
// from Preparata & Shamos ['85, p. 316, "A simple case ..."]. Translate the
// halfspaces so that the interior point is the origin. Calculate the dual
// polytope. The dual polytope is the convex hull of the vertices dual to the
// original faces in regard to the unit sphere (i.e., halfspaces at distance
// d from the origin are dual to vertices at distance 1/d). Then calculate
// the resulting polytope, which is the dual of the dual polytope, and
// translate the origin back to the interior point [S. Spitz and S. Teller].

// NOTE here we change this to max(x4) under constraints aj*x1 + bj*x2 + cj*x3 + dj - x4 >= 0, j=1..n
//  i.e. aj*x1 + bj*x2 + cj*x3 - x4 >= -dj, j=1..n
// Then, if [x1,x2,x3,x4] is an optimal m_solution with x3 > 0 we pick
// the point [x1,x2,x3] as inside point.

template <class Kernel, class ET>
class Interior_polyhedron_3 {
        // 3D
        typedef typename Kernel::FT FT;
        typedef typename Kernel::Plane_3 Plane;
        typedef typename Kernel::Point_3 Point;

        // program and solution types
        typedef CGAL::Quadratic_program<double> LP;
        typedef typename CGAL::Quadratic_program_solution<ET> Solution;
        typedef typename Solution::Variable_value_iterator Variable_value_iterator;
        typedef CGAL::Real_embeddable_traits<typename Variable_value_iterator::value_type> RE_traits;
        typename RE_traits::To_double to_double;
        Solution m_solution;
        Point m_inside_point;
        Point m_optimal_point;

    public:
        Point& inside_point() { return m_inside_point; }
        const Point& inside_point() const { return m_inside_point; }
        Point& optimal_point() { return m_optimal_point; }
        const Point& optimal_point() const { return m_optimal_point; }

        // Determines if a value is infinite or not
        template<typename T>
        inline bool isinf(T value) {
            return value == std::numeric_limits<T>::infinity();
        }

        // Find a point inside the polyhedron defined by a list of planes
        // InputIterator::value_type = Plane
        template < class InputIterator >
        bool find(InputIterator begin, InputIterator end) {
            // solve linear program
            LP lp(CGAL::LARGER,false); // with constraints Ax >= b

            // column indices
            const int index_x1 = 0;
            const int index_x2 = 1;
            const int index_x3 = 2;
            const int index_x4 = 3;

            // assemble linear program
            int j = 0; // row index
            // iterate over segments
            InputIterator it;
            for(it = begin; it != end; ++it, j++) {
                const Plane& plane = *it;
                const double aj = CGAL::to_double(plane.a());
                const double bj = CGAL::to_double(plane.b());
                const double cj = CGAL::to_double(plane.c());
                const double dj = CGAL::to_double(plane.d());

                CGAL_assertion(!isinf(aj));
                CGAL_assertion(!isinf(bj));
                CGAL_assertion(!isinf(cj));
                CGAL_assertion(!isinf(dj));

                // plane defined the halfspace: aj * x1 + bj * x2 + cj * x3 + dj <= 0
                // <=> - (aj * x1 + bj * x2 + cj * x3 + dj) >= 0
                // j^th constraint: -(aj * x1 + bj * x2 + cj * x3 + x4) >= dj
                lp.set_a(index_x1, j,   -aj);
                lp.set_a(index_x2, j,   -bj);
                lp.set_a(index_x3, j,   -cj);
                lp.set_a(index_x4, j, -1.0);

                // right hand side
                lp.set_b(j, dj);
            }

            // objective function -> max x4 (negative sign set because
            // the lp solver always minimizes an objective function)
            lp.set_c(index_x4,-1.0);

            // solve the linear program
            m_solution = CGAL::solve_linear_program(lp, ET());

            if(m_solution.is_infeasible())
                return false;

            if(!m_solution.is_optimal())
                return false;

            // get variables
            Variable_value_iterator X = m_solution.variable_values_begin();

            // solution if x4 > 0
            double x4 = to_double(X[index_x4]);
            if(x4 <= 0.0)
                return false;

            // define inside point as (x1;x2;x3)
            double x1 = to_double(X[index_x1]);
            double x2 = to_double(X[index_x2]);
            double x3 = to_double(X[index_x3]);
            m_inside_point = Point(x1,x2,x3);

            return true;
        }
};

#endif // INTERIOR_POLYHEDRON_3_H

