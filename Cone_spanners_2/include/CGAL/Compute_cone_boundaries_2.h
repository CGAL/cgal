// Copyright (c) 2013-2015  The University of Western Sydney, Australia.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

#include <CGAL/license/Cone_spanners_2.h>


#include <iostream>
#include <cstdlib>
#include <vector>
#include <utility>
#include <CGAL/config.h>      // included compiler_config.h, defining CGAL_USE_CORE, etc.
#include <CGAL/Polynomial.h>
#include <CGAL/number_type_config.h>    // CGAL_PI is defined there
#include <CGAL/enum.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_root_of.h>
#include <CGAL/Aff_transformation_2.h>

namespace CGAL {

/*! \ingroup PkgConeSpanners2Ref
 *
 *  \brief The functor for computing the directions of cone boundaries with a given
 *  cone number and a given initial direction.
 *
 *  This computation can be either inexact by simply dividing an approximate \f$ \pi \f$ by the cone number
 *  (which is quick), or exact by using roots of polynomials (requiring number types such as `CORE::Expr` or `leda_real`,
 *  which are slow). The inexact computation is done by the general functor definition,
 *  while the exact computation is done by a specialization of this functor.
 *
 *  In the construction of cone-based spanners such as Yao graph and Theta graph implemented by this package,
 *  this functor is called first to compute the cone boundaries.
 *  Of course, this functor can also be used in other applications where the plane needs to be divided
 *  into equally-angled cones.
 *
 * \tparam Traits_  Must be either `CGAL::Exact_predicates_exact_constructions_kernel_with_root_of`
 *                  or `CGAL::Exact_predicates_inexact_constructions_kernel`.
 *
 */
template <typename Traits_>
class Compute_cone_boundaries_2 {

public:
    /*! the geometric traits class. */
    typedef  Traits_      Traits;

    /*! the direction type. */
    typedef  typename Traits::Direction_2       Direction_2;

private:

    typedef  typename Traits::Aff_transformation_2    Transformation;

public:
    /* Note: No member variables in this class, so a custom constructor is not needed. */

    /*! \brief The operator().
     *
     * \details The direction of the first ray can be specified by the parameter `initial_direction`,
     * which allows the first ray to start at any direction.
     * This operator first places the `initial_direction` at the
     * position pointed by `result`. Then, it calculates the remaining directions (cone boundaries)
     * and output them to `result` in the counterclockwise order.
     * Finally, the past-the-end iterator for the resulting directions is returned.
     *
         * \tparam DirectionOutputIterator  an `OutputIterator` with value type `Direction_2`.
     * \param cone_number The number of cones
     * \param initial_direction The direction of the first ray
     * \param result  The output iterator
     */
    template<class DirectionOutputIterator>
    DirectionOutputIterator operator()(const unsigned int cone_number,
                                       const Direction_2& initial_direction,
                                       DirectionOutputIterator result)  {
        if (cone_number<2) {
            std::cout << "The number of cones must be larger than 1!" << std::endl;
            CGAL_assertion(false);
        }

        *result++ = initial_direction;

        const double cone_angle = 2*CGAL_PI/cone_number;
        double sin_value, cos_value;
        Direction_2 ray;
        for (unsigned int i = 1; i < cone_number; i++) {
            sin_value = std::sin(i*cone_angle);
            cos_value = std::cos(i*cone_angle);
            ray = Transformation(cos_value, -sin_value, sin_value, cos_value)(initial_direction);
            *result++ = ray;
        }

        return result;
    } // end of operator

};


/*
 The specialised functor for computing the directions of cone boundaries exactly
 with a given cone number and a given initial direction.
*/
template <>
class Compute_cone_boundaries_2<Exact_predicates_exact_constructions_kernel_with_root_of> {

public:
    /* Indicate the type of the cgal kernel. */
    typedef  Exact_predicates_exact_constructions_kernel_with_root_of  Kernel_type;

private:
    typedef Kernel_type::FT            FT;
    typedef Kernel_type::Direction_2   Direction_2;
    typedef Kernel_type::Aff_transformation_2   Transformation;

public:
    /* Note: No member variables in this class, so a Constructor is not needed. */

    /* The operator().

      The direction of the first ray can be specified by the parameter
      initial_direction, which allows the first ray to start at any direction.
      The remaining directions are calculated in counter-clockwise order.

      \param cone_number The number of cones
      \param initial_direction The direction of the first ray
      \param result  The output iterator
    */
    template<typename DirectionOutputIterator>
    DirectionOutputIterator operator()(const unsigned int cone_number,
                                       const Direction_2& initial_direction,
                                       DirectionOutputIterator result)  {

        if (cone_number<2) {
            std::cout << "The number of cones must be larger than 1!" << std::endl;
            std::exit(1);
        }

        // Since CGAL::root_of() gives the k-th smallest root,
        // here -x is actually used instead of x.
        // But we want the second largest one with no need to count.
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
        bool is_even;
        if (cone_number % 2 == 0) {
            is_even = true;
            m = cone_number/2;       // for even number of cones
        }
        else {
            m= cone_number/2 + 1;    // for odd number of cones
            is_even = false;
        }

        FT cos_value, sin_value;
        // for storing the intermediate result
        Direction_2 ray;
        // For saving the first half number of rays when cone_number is even
        std::vector<Direction_2> ray_store;

        // add the first half number of rays in counter clockwise order
        for (i = 1; i <= m; i++) {
            cos_value = - root_of(i, a.begin(), a.end());
            sin_value = sqrt(FT(1) - cos_value*cos_value);
            ray = Transformation(cos_value, -sin_value, sin_value, cos_value)(initial_direction);
            *result++ = ray;
            if (is_even)
                ray_store.push_back(ray);
        }

        // add the remaining half number of rays in ccw order
        if (is_even) {
            for (i = 0; i < m; i++) {
                *result++ = -ray_store[i];
            }
        } else {
            for (i = 0; i < m-1; i++) {
                cos_value = - root_of(m-i, a.begin(), a.end());
                sin_value = - sqrt(FT(1) - cos_value*cos_value);
                ray = Transformation(cos_value, -sin_value, sin_value, cos_value)(initial_direction);
                *result++ = ray;
            }
        }

        return result;

    };      // end of operator()
};      // end of functor specialization: Compute_cone_..._2

}  // namespace CGAL

#endif
