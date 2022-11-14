// Copyright (c) 2013  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Olivier Devillers
//                 Pedro Machado Manhaes de Castro

#ifndef CGAL_HILBERT_SORT_ON_SPHERE_3_H
#define CGAL_HILBERT_SORT_ON_SPHERE_3_H

#include <CGAL/Hilbert_sort_2.h>
#include <CGAL/Spatial_sorting/internal/Transform_coordinates_traits_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/double.h>
#include <algorithm>
#include <vector>

namespace CGAL {

template <class K,
          class Hilbert_policy,
          class P = typename K::Point_3>
class Hilbert_sort_on_sphere_3
{
        typedef P Point_3;
        typedef typename K::FT FT;

        // Face 1, x > sqrt(1/3)
        // Face 2, y > sqrt(1/3)
        // Face 3, x < -sqrt(1/3)
        // Face 4, z > sqrt(1/3)
        // Face 5, y < -sqrt(1/3)
        // Face 6, z < -sqrt(1/3)

        typedef internal::Transform_coordinates_traits_3<K,0,1,1,0> Face_1_traits_3;  // +y +z
        typedef internal::Transform_coordinates_traits_3<K,-1,0,1,0> Face_2_traits_3; // -x +z
        typedef internal::Transform_coordinates_traits_3<K,0,1,-1,1> Face_3_traits_3; // +z -y
        typedef internal::Transform_coordinates_traits_3<K,-1,1,0,1> Face_4_traits_3;  // -y +x
        typedef internal::Transform_coordinates_traits_3<K,-1,0,1,1> Face_5_traits_3;  // -z +x
        typedef internal::Transform_coordinates_traits_3<K,1,1,0,0> Face_6_traits_3;  // +x +y



        Hilbert_sort_2<Face_1_traits_3, Hilbert_policy > _hs_1_object;
        Hilbert_sort_2<Face_2_traits_3, Hilbert_policy > _hs_2_object;
        Hilbert_sort_2<Face_3_traits_3, Hilbert_policy > _hs_3_object;
        Hilbert_sort_2<Face_4_traits_3, Hilbert_policy > _hs_4_object;
        Hilbert_sort_2<Face_5_traits_3, Hilbert_policy > _hs_5_object;
        Hilbert_sort_2<Face_6_traits_3, Hilbert_policy > _hs_6_object;

        K _k;
        Point_3 _p;
        FT _sq_r;
        const FT _sqrt_of_one_over_three;

public:
        Hilbert_sort_on_sphere_3 (const K &k = K(),
                                  FT sq_r = 1.0,
                                  const Point_3 &p = Point_3(0,0,0),
                                  std::ptrdiff_t limit=1)
        : _hs_1_object(Face_1_traits_3(),limit),
          _hs_2_object(Face_2_traits_3(),limit),
          _hs_3_object(Face_3_traits_3(),limit),
          _hs_4_object(Face_4_traits_3(),limit),
          _hs_5_object(Face_5_traits_3(),limit),
          _hs_6_object(Face_6_traits_3(),limit),
          _k(k), _p(p), _sq_r(sq_r),
          _sqrt_of_one_over_three(CGAL_NTS sqrt(FT(1)/FT(3)))
        {
                CGAL_precondition( sq_r > 0 );
        }


        template <class RandomAccessIterator>
        void operator()(RandomAccessIterator begin, RandomAccessIterator end) const {
                typedef typename std::iterator_traits<RandomAccessIterator>::value_type Point;
                std::vector< Point > vec[6];

                const FT mulcte = _sqrt_of_one_over_three * CGAL_NTS sqrt(_sq_r);
                const FT lxi = _p.x() - mulcte, lxs = _p.x() + mulcte;
                const FT lyi = _p.y() - mulcte, lys = _p.y() + mulcte;
                const FT lzs = _p.z() + mulcte;

                for(RandomAccessIterator i = begin; i != end; ++i) {
                        const Point &p = *i;
                        const typename K::FT x = _k.compute_x_3_object()(p);
                        const typename K::FT y = _k.compute_y_3_object()(p);
                        const typename K::FT z = _k.compute_z_3_object()(p); // for unit sphere
                        if(x > lxs) vec[0].push_back(p);             // Face 1, x > sqrt(1/3)
                        else if(y > lys) vec[1].push_back(p);        // Face 2, y > sqrt(1/3)
                        else if(x < lxi) vec[2].push_back(p);        // Face 3, x < -sqrt(1/3)
                        else if(z > lzs) vec[3].push_back(p);        // Face 4, z > sqrt(1/3)
                        else if(y < lyi) vec[4].push_back(p);        // Face 5, y < -sqrt(1/3)
                        else vec[5].push_back(p);                    // Face 6, z < -sqrt(1/3)
                }
                if(vec[0].size()) _hs_1_object(vec[0].begin(), vec[0].end());
                if(vec[1].size()) _hs_2_object(vec[1].begin(), vec[1].end());
                if(vec[2].size()) _hs_3_object(vec[2].begin(), vec[2].end());
                if(vec[3].size()) _hs_4_object(vec[3].begin(), vec[3].end());
                if(vec[4].size()) _hs_5_object(vec[4].begin(), vec[4].end());
                if(vec[5].size()) _hs_6_object(vec[5].begin(), vec[5].end());

                // this is the order that set of points in a face should appear
                // after sorting points wrt each face
                for(int i=0; i<6; i++)
                        for(std::size_t j=0; j<vec[i].size(); j++)
                                *begin++ = vec[i][j];
        }
};

} // namespace CGAL

#endif//CGAL_HILBERT_SORT_ON_SPHERE_3_H
