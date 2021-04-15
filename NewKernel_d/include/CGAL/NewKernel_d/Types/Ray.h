// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_KERNELD_TYPES_RAY_H
#define CGAL_KERNELD_TYPES_RAY_H
#include <utility>
#include <CGAL/NewKernel_d/functor_tags.h>
#include <CGAL/Kernel/mpl.h>
namespace CGAL {
template <class R_> class Ray {
        typedef typename Get_type<R_, FT_tag>::type FT_;
        typedef typename Get_type<R_, Point_tag>::type        Point_;
        typedef typename Get_type<R_, Vector_tag>::type        Vector_;
        typedef std::pair<Point_,Vector_> Data_;
        Data_ data;
        public:
        Ray(){}
        Ray(Point_ const&a, Vector_ const&b): data(a,b) {}
        Point_ source()const{
          return data.first;
        }
        // FIXME: return a R_::Direction?
        Vector_ direction()const{
          return data.second;
        }
};
namespace CartesianDKernelFunctors {
  template <class R_> struct Construct_ray : Store_kernel<R_> {
    CGAL_FUNCTOR_INIT_STORE(Construct_ray)
    typedef typename Get_type<R_, Ray_tag>::type        result_type;
    typedef typename Get_type<R_, Point_tag>::type        Point;
    typedef typename Get_type<R_, Vector_tag>::type        Vector;
    typedef typename Get_functor<R_, Difference_of_points_tag>::type Dp_;
    //typedef typename Get_functor<R_, Translated_point_tag>::type Tp_;
    //typedef typename Get_functor<R_, Scaled_vector_tag>::type Sv_;
    result_type operator()(Point const&a, Vector const&b)const{
      return result_type(a,b);
    }
    result_type operator()(Point const&a, typename First_if_different<Point,Vector>::Type const&b)const{
      Dp_ dp(this->kernel());
      return result_type(a,dp(b,a));
    }
  };
}
CGAL_KD_DEFAULT_TYPE(Ray_tag,(CGAL::Ray<K>),(Point_tag,Vector_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Construct_ttag<Ray_tag>,(CartesianDKernelFunctors::Construct_ray<K>),(Point_tag,Ray_tag,Vector_tag),(Difference_of_points_tag));

} // namespace CGAL

#endif // CGAL_KERNELD_TYPES_RAY_H
