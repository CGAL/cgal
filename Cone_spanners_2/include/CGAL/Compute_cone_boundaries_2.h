// Copyright (c) 2013-2015  The University of Western Sydney, Australia.
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
// Authors: Weisheng Si, Quincy Tse

/*! \file Compute_cone_boundaries_2.h
 *
 * This header implements the functor for computing the directions of cone boundaries with a given
 *  cone number and a given initial direction either exactly or inexactly.
 */

#ifndef CGAL_COMPUTE_CONE_BOUNDARIES_2_H
#define CGAL_COMPUTE_CONE_BOUNDARIES_2_H

// if leda::real is used, pls modify the following definition
#define CGAL_USE_CORE 1

#include <iostream>
#include <cstdlib>
#include <utility>
#include <CGAL/Polynomial.h>
#include <CGAL/number_type_config.h>    // defining CGAL_PI
#include <CGAL/enum.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Aff_transformation_2.h>

namespace CGAL {

/*! \ingroup PkgConeBasedSpanners
 *
 *  \brief The functor for computing the directions of cone boundaries with a given
 *  cone number and a given initial direction. 
 *
 *  This computation can be either inexact by simply dividing an approximate Pi by the cone number
 *  (which is quick), or exact by using roots of polynomials (requiring number types such as `CORE::Expr` or `LEDA::Real`,
 *  which are slow). The inexact computation is done by the general functor definition,
 *  while the exact computation is done by a specialization of this functor.
 *   
 *  \tparam Kernel_   If this parameter is `Exact_predicates_exact_constructions_kernel_with_sqrt`,
 *                    the specialization functor will be called; otherwise, the general functor will
 *                    be called.
 *
 *  In the construction of Yao graph and Theta graph implemented by this package,
 *  all predicates and construction functions are from \cgal.
 *  Therefore, if the kernel `Exact_predicates_exact_constructions_kernel_with_sqrt` is used,
 *  the Yao or Theta graph will be constructed exactly, otherwise inexactly.
 *
 *  Of course, this functor can also be used in other applications where the plane needs to be divided
 *  into equally-angled cones.
 *
 *  \cgalModels `ComputeConeBoundaries_2`
 *
 */
template <typename Kernel_>
class Compute_cone_boundaries_2 {

public:
	/*! Indicate the type of the \cgal kernel. */
    typedef  Kernel_                    Kernel_type;

private:
    typedef  typename Kernel_::FT                FT;
    typedef  typename Kernel_::Direction_2       Direction_2;
	typedef  typename Kernel_::Aff_transformation_2    Transformation;

public:
	/* No member variables in this class, so a custom constructor is not needed. */
	// Compute_cone_boundaries_2() {};

	/*! \brief The operator(). 
     *
	 * \details The direction of the first ray can be specified by the parameter `initial_direction`, 
	 * which allows the first ray to start at any direction. The remaining rays are calculated in
	 * counter-clockwise order.
	 *
	 * \param[in] cone_number The number of cones
	 * \param[in] initial_direction The direction of the first ray
	 * \param[out] rays  The results, a vector of directions
	 */
    void operator()(const unsigned int cone_number,
                    Direction_2& initial_direction,
                    std::vector<Direction_2>& rays)  {
        if (cone_number<2) {
            std::cout << "The number of cones should be larger than 1!" << std::endl;
            std::exit(1);
        }

        if (rays.size() > 0) {
            std::cout << "Initially, the vector rays must contain no elements!" << std::endl;
            std::exit(1);
        }

        rays.push_back(initial_direction);

        const double cone_angle = 2*CGAL_PI/cone_number;
        double sin_value, cos_value;
        for (unsigned int i = 1; i < cone_number; i++) {
            sin_value = std::sin(i*cone_angle);
            cos_value = std::cos(i*cone_angle);
            Direction_2 ray_i = Transformation(cos_value, -sin_value, sin_value, cos_value)(initial_direction);
            rays.push_back(ray_i);
        }
    }

};


/* 
 The specialised functor for computing the directions of cone boundaries exactly
 with a given cone number and a given initial direction. 
*/
template <>
class Compute_cone_boundaries_2<Exact_predicates_exact_constructions_kernel_with_sqrt> {

public:
	/* Indicate the type of the cgal kernel. */
    typedef  Exact_predicates_exact_constructions_kernel_with_sqrt                kernel_type;

private:
    typedef  typename Exact_predicates_exact_constructions_kernel_with_sqrt::FT            FT;
    typedef  typename Exact_predicates_exact_constructions_kernel_with_sqrt::Direction_2   Direction_2;
	typedef  typename Exact_predicates_exact_constructions_kernel_with_sqrt::Aff_transformation_2   Transformation;

public:
	/* No member variables in this class, so a Constructor is not needed. */
	// Compute_cone_boundaries_2() {};

	/* The operator().  

      The direction of the first ray can be specified by the parameter
      initial_direction, which allows the first ray to start at any direction. 
      The remaining rays are calculated in counter-clockwise order.

      \param[in] cone_number The number of cones
      \param[in] initial_direction The direction of the first ray
      \param[out] rays  The results, a vector of directions
    */
    void operator()(const unsigned int cone_number,
                    Direction_2& initial_direction,
                    std::vector< Direction_2 >& rays) {

        if (cone_number<2) {
            std::cout << "The number of cones should be larger than 1!" << std::endl;
            std::exit(1);
        }

        if (rays.size() > 0) {
            std::cout << "Initially, the vector rays must contain no elements!" << std::endl;
            std::exit(1);
        }

        //std::cout << "Specialization is called!" << std::endl;

        // here -x is actually used instead of x, since CGAL::root_of() gives the 
		// k-th smallest root, but we want the second largest one with no need to count.
        Polynomial<FT> x(CGAL::shift(Polynomial<FT>(-1), 1));
        Polynomial<FT> double_x(2*x);
        Polynomial<FT> a(1), b(x);
        for (unsigned int i = 2; i <= cone_number; ++i) {
            Polynomial<FT> c = double_x*b - a;
            a = b;
            b = c;
        }
        a = b - 1;

        unsigned int m, i;
        if (cone_number % 2 == 0)
            m = cone_number/2;       // for even number of cones
        else
            m= cone_number/2 + 1;    // for odd number of cones

        FT cos_value, sin_value;
        Direction_2 ray_i;
        // add the first half number of rays in counter clockwise order
        for (i = 1; i <= m; i++) {
            cos_value = - root_of(i, a.begin(), a.end());
            sin_value = sqrt(FT(1) - cos_value*cos_value);
            ray_i = Transformation(cos_value, -sin_value, sin_value, cos_value)(initial_direction);
            rays.push_back(ray_i);
        }

        // add the remaining half number of rays in ccw order
        if (cone_number % 2 == 0) {
            for (i = 0; i < m; i++) {
                rays.push_back(-rays[i]);
            }
        } else {
            for (i = 0; i < m-1; i++) {
                cos_value = - root_of(m-i, a.begin(), a.end());
                sin_value = - sqrt(FT(1) - cos_value*cos_value);
                ray_i = Transformation(cos_value, -sin_value, sin_value, cos_value)(initial_direction);
                rays.push_back(ray_i);
            }
        }

    };      // end of operator()
};      // end of functor specialization: Compute_cone_..._2

}  // namespace CGAL

#endif
