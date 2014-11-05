#ifndef CGAL_VORONOI_COVARIANCE_SPHERE_3_HPP
#define CGAL_VORONOI_COVARIANCE_SPHERE_3_HPP

#include <CGAL/point_generators_3.h>

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
                            static const FT phi = (FT(1) + ::sqrt(5))/FT(2);
                            static const FT s = FT(1) / ::sqrt(phi + FT(2));

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
                            const FT phi = (FT(1) + ::sqrt(5))/FT(2);
                            const FT one_phi = FT(1)/phi;
                            const FT s = FT(1) / ::sqrt(3.0);

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

    }; // namespace Voronoi_covariance_3
}; // namespace CGAL

#endif // CGAL_VORONOI_COVARIANCE_SPHERE_3_HPP

