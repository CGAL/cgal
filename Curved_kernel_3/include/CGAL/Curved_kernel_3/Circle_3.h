// Copyright (c) 2005  INRIA Sophia-Antipolis (France) 
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Julien Hazebrouck
//           Damien Leroy
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_CIRCLE_3_H
#define CGAL_SPHERICAL_KERNEL_CIRCLE_3_H

#include <CGAL/Curved_kernel_3/internal_functions_on_sphere_3.h>

namespace CGAL {
  namespace CGALi{
    template <class SK> class Circle_3 {

      typedef typename SK::Plane_3      Plane_3;
      typedef typename SK::Sphere_3     Sphere_3;
      typedef typename SK::Point_3      Point_3;
      typedef typename SK::FT           FT;

    private:
      typedef std::pair<Sphere_3, Plane_3>  Rep;
      typedef typename SK::template Handle<Rep>::type  Base;

      Base base;

    public:
      Circle_3() {}

      Circle_3(const Circle_3 &c) 
      {
        base = Rep(c.diametral_sphere(), c.supporting_plane());
      }

      Circle_3(const Point_3& center, const FT& squared_r, const Plane_3& p)
      {
        // the plane contains the center
	CGAL_kernel_assertion((p.a() * center.x() +
                               p.b() * center.y() +
                               p.c() * center.z() +
                               p.d()) == ZERO);
        // It is not allowed non-positive radius 
        // Should we keep this pre-condition?
        CGAL_kernel_assertion(squared_r > 0);
        base = Rep(Sphere_3(center,squared_r), p);
      }

      Circle_3(const Sphere_3 &s1, 
               const Sphere_3 &s2) {
         std::vector<Object> sols;
         SK().intersect_3_object()(s1, s2, std::back_inserter(sols));
         // s1,s2 must intersect
         CGAL_kernel_precondition(sols.size() != 0);
         typename SK::Circle_3 circle;
         // the intersection must be a circle (no point allowed)
         CGAL_kernel_precondition(assign(circle,sols[0]));
         assign(circle,sols[0]);
         *this = circle.rep();
      }

      Circle_3(const Plane_3 &p, 
               const Sphere_3 &s) {
         std::vector<Object> sols;
         SK().intersect_3_object()(p, s, std::back_inserter(sols));
         // s1,s2 must intersect
         CGAL_kernel_precondition(sols.size() != 0);
         typename SK::Circle_3 circle;
         // the intersection must be a circle (no point allowed)
         CGAL_kernel_precondition(assign(circle,sols[0]));
         assign(circle,sols[0]);
         *this = circle.rep();
      }

      const Plane_3& supporting_plane () const {
        return get(base).second;
      }

      const Point_3& center () const {
        return get(base).first.center();
      }

      const FT& squared_radius () const {
        return get(base).first.squared_radius();
      }

      const Sphere_3& diametral_sphere() const {
        return get(base).first;
      }

      bool operator==(const Circle_3 &) const;
      bool operator!=(const Circle_3 &) const;

    };

    template < class SK >
    CGAL_KERNEL_INLINE
    bool
    Circle_3<SK>::operator==(const Circle_3<SK> &t) const
    {
      if (CGAL::identical(base, t.base))
        return true;
      return CGAL::SphericalFunctors::non_oriented_equal<SK>(*this, t);
    } 

    template < class SK >
    CGAL_KERNEL_INLINE
    bool
    Circle_3<SK>::operator!=(const Circle_3<SK> &t) const
    {
      return !(*this == t);
    }

  }
}

#endif

