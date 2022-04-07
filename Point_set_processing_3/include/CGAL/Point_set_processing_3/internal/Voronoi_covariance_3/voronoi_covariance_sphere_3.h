// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Jocelyn Meyron and Quentin Mérigot
//

#ifndef CGAL_INTERNAL_VCM_VORONOI_COVARIANCE_SPHERE_3_H
#define CGAL_INTERNAL_VCM_VORONOI_COVARIANCE_SPHERE_3_H

#include <CGAL/license/Point_set_processing_3.h>


#include <CGAL/point_generators_3.h>
#include <cmath>
#include <cstddef>

/// \cond SKIP_IN_MANUAL

namespace CGAL {
    namespace Voronoi_covariance_3 {
        template <class K>
            class Sphere_discretization
            {
                typedef typename K::FT FT;
                FT R;
                std::size_t  N;

                public:
              Sphere_discretization (FT R, std::size_t N = 20) :
                    R(R), N(N)
                {
                    if (N != 8 && N != 20)
                        N = 20;
                }

                template <class OutputIterator>
                    void
                    operator ()(OutputIterator out) const
                    {
                        typedef typename K::Plane_3 Plane;

                        if (N == 8)
                        {
                            static const FT phi = (FT(1) + std::sqrt(5.0))/FT(2);
                            static const FT s = FT(1) / std::sqrt(phi + FT(2));

                            *out ++ = Plane(0, +s, +s*phi, -R);
                            *out ++ = Plane(0, -s, +s*phi, -R);
                            *out ++ = Plane(0, +s, -s*phi, -R);
                            *out ++ = Plane(0, -s, -s*phi, -R);

                            *out ++ = Plane(+s, +s*phi, 0, -R);
                            *out ++ = Plane(+s, -s*phi, 0, -R);
                            *out ++ = Plane(-s, +s*phi, 0, -R);
                            *out ++ = Plane(-s, -s*phi, 0, -R);

                            *out ++ = Plane(+s*phi, 0, +s, -R);
                            *out ++ = Plane(-s*phi, 0, +s, -R);
                            *out ++ = Plane(+s*phi, 0, -s, -R);
                            *out ++ = Plane(-s*phi, 0, -s, -R);
                        }
                        else if (N == 20)
                        {
                            const FT phi = (FT(1) + std::sqrt(5.0))/FT(2);
                            const FT one_phi = FT(1)/phi;
                            const FT s = FT(1) / std::sqrt(3.0);

                            *out ++ = Plane(+s, +s, +s, -R);
                            *out ++ = Plane(-s, +s, +s, -R);
                            *out ++ = Plane(+s, -s, +s, -R);
                            *out ++ = Plane(-s, -s, +s, -R);
                            *out ++ = Plane(+s, +s, -s, -R);
                            *out ++ = Plane(-s, +s, -s, -R);
                            *out ++ = Plane(+s, -s, -s, -R);
                            *out ++ = Plane(-s, -s, -s, -R);

                            *out ++ = Plane(0, +s*one_phi, +s*phi, -R);
                            *out ++ = Plane(0, -s*one_phi, +s*phi, -R);
                            *out ++ = Plane(0, +s*one_phi, -s*phi, -R);
                            *out ++ = Plane(0, -s*one_phi, -s*phi, -R);

                            *out ++ = Plane(+s*one_phi, +s*phi, 0, -R);
                            *out ++ = Plane(-s*one_phi, +s*phi, 0, -R);
                            *out ++ = Plane(+s*one_phi, -s*phi, 0, -R);
                            *out ++ = Plane(-s*one_phi, -s*phi, 0, -R);

                            *out ++ = Plane(+s*phi, 0, +s*one_phi, -R);
                            *out ++ = Plane(-s*phi, 0, +s*one_phi, -R);
                            *out ++ = Plane(+s*phi, 0, -s*one_phi, -R);
                            *out ++ = Plane(-s*phi, 0, -s*one_phi, -R);
                        }
                    }
            };

    } // namespace Voronoi_covariance_3
} // namespace CGAL

/// \endcond

#endif // CGAL_INTERNAL_VCM_VORONOI_COVARIANCE_SPHERE_3_H

