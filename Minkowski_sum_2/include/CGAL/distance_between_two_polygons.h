//
// Created by kabir on 11/04/21.
//

#pragma once

#include <CGAL/minkowski_sum_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Aff_transformation_2.h>

namespace CGAL {


    template<typename Kernel_ , typename Container_>
    typename Kernel_::FT
    squared_distance_from_origin(const Polygon_with_holes_2<Kernel_, Container_>& sum) {
        typedef Kernel_                                         Kernel;
        typedef Container_                                      Container;
        typedef CGAL::Polygon_2<Kernel, Container>              Polygon_2;
        typedef CGAL::Polygon_with_holes_2<Kernel, Container>   Polygon_with_holes_2;
        typedef typename Kernel::Point_2                        Point_2;

        Kernel kernel;

        CGAL::Origin origin;
        auto comp_sqr_distance = kernel.compute_squared_distance_2_object();
        auto n = sum.outer_boundary().size();
        const auto &edge = sum.outer_boundary().edge(1);
        auto sqr_distance = comp_sqr_distance(origin, edge);
        for (auto i = 1; i < n; ++i) {
            const auto &edge = sum.outer_boundary().edge(i);
            auto tmp = comp_sqr_distance(origin, edge);
            if (tmp < sqr_distance) sqr_distance = tmp;
        }
        return sqr_distance;
    }

    template<typename Kernel_, typename Container_ >
    typename Kernel_::FT
    distance_between_two_polygons(const Polygon_2<Kernel_, Container_>& P,
                                  const Polygon_2<Kernel_, Container_>& Q) {
        typedef Kernel_                                       Kernel;
        typedef Container_                                    Container;
        typedef CGAL::Polygon_2<Kernel, Container>            Polygon_2;
        typedef CGAL::Polygon_with_holes_2<Kernel, Container> Polygon_with_holes_2;
        typedef typename Kernel::Aff_transformation_2         Transformation;
        typedef typename Kernel::Point_2                      Point_2;
        typedef typename Kernel_::FT                          FT;
        Kernel kernel;

        Transformation rotate(CGAL::ROTATION, sin(M_PI), cos(M_PI));

        Polygon_2 Q_Reflection;
        for (auto &i : Q) {
            //Reflection
            Q_Reflection.push_back(rotate(i));
        }
        Polygon_with_holes_2 sum = CGAL::minkowski_sum_2(P, Q_Reflection);
        Kernel traits;
        FT x;
        switch (CGAL::bounded_side_2(sum.outer_boundary().begin(), sum.outer_boundary().end(), Point_2(0.0, 0.0), traits)) {
            case CGAL::ON_BOUNDED_SIDE :
                x = 0;
                break;
            case CGAL::ON_BOUNDARY:
            case CGAL::ON_UNBOUNDED_SIDE:
                x = squared_distance_from_origin(sum);
                break;
        }
        return x;
    }

}


