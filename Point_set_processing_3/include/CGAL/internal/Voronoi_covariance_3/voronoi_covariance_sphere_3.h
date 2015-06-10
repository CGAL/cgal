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
// Author(s) : Jocelyn Meyron and Quentin Mérigot
//

#ifndef CGAL_INTERNAL_VCM_VORONOI_COVARIANCE_SPHERE_3_HPP
#define CGAL_INTERNAL_VCM_VORONOI_COVARIANCE_SPHERE_3_HPP

#include <CGAL/point_generators_3.h>
#include <cmath>

namespace CGAL {
    namespace Voronoi_covariance_3 {
        template <class K>
            class Sphere_discretization
            {
                typedef typename K::FT FT;
                FT _R;
                size_t _N;

                public:
                Sphere_discretization (FT R, size_t N = 20) :
                    _R(R), _N(N)
                {
                    if (_N != 8 && _N != 20)
                        _N = 20;
                }

                template <class OutputIterator>
                    void
                    operator ()(OutputIterator out) const
                    {
                        typedef typename K::Plane_3 Plane;

                        if (_N == 8)
                        {
                            static const FT phi = (FT(1) + std::sqrt(5.0))/FT(2);
                            static const FT s = FT(1) / std::sqrt(phi + FT(2));

                            *out ++ = Plane(0, +s, +s*phi, -_R);
                            *out ++ = Plane(0, -s, +s*phi, -_R);
                            *out ++ = Plane(0, +s, -s*phi, -_R);
                            *out ++ = Plane(0, -s, -s*phi, -_R);

                            *out ++ = Plane(+s, +s*phi, 0, -_R);
                            *out ++ = Plane(+s, -s*phi, 0, -_R);
                            *out ++ = Plane(-s, +s*phi, 0, -_R);
                            *out ++ = Plane(-s, -s*phi, 0, -_R);

                            *out ++ = Plane(+s*phi, 0, +s, -_R);
                            *out ++ = Plane(-s*phi, 0, +s, -_R);
                            *out ++ = Plane(+s*phi, 0, -s, -_R);
                            *out ++ = Plane(-s*phi, 0, -s, -_R);
                        }
                        else if (_N == 20)
                        {
                            const FT phi = (FT(1) + std::sqrt(5.0))/FT(2);
                            const FT one_phi = FT(1)/phi;
                            const FT s = FT(1) / std::sqrt(3.0);

                            *out ++ = Plane(+s, +s, +s, -_R);
                            *out ++ = Plane(-s, +s, +s, -_R);
                            *out ++ = Plane(+s, -s, +s, -_R);
                            *out ++ = Plane(-s, -s, +s, -_R);
                            *out ++ = Plane(+s, +s, -s, -_R);
                            *out ++ = Plane(-s, +s, -s, -_R);
                            *out ++ = Plane(+s, -s, -s, -_R);
                            *out ++ = Plane(-s, -s, -s, -_R);

                            *out ++ = Plane(0, +s*one_phi, +s*phi, -_R);
                            *out ++ = Plane(0, -s*one_phi, +s*phi, -_R);
                            *out ++ = Plane(0, +s*one_phi, -s*phi, -_R);
                            *out ++ = Plane(0, -s*one_phi, -s*phi, -_R);

                            *out ++ = Plane(+s*one_phi, +s*phi, 0, -_R);
                            *out ++ = Plane(-s*one_phi, +s*phi, 0, -_R);
                            *out ++ = Plane(+s*one_phi, -s*phi, 0, -_R);
                            *out ++ = Plane(-s*one_phi, -s*phi, 0, -_R);

                            *out ++ = Plane(+s*phi, 0, +s*one_phi, -_R);
                            *out ++ = Plane(-s*phi, 0, +s*one_phi, -_R);
                            *out ++ = Plane(+s*phi, 0, -s*one_phi, -_R);
                            *out ++ = Plane(-s*phi, 0, -s*one_phi, -_R);
                        }
                    }
            };

    } // namespace Voronoi_covariance_3
} // namespace CGAL

#endif // CGAL_INTERNAL_VCM_VORONOI_COVARIANCE_SPHERE_3_HPP

