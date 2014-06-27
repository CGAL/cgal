// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Marc Glisse

#ifndef CGAL_KERNELD_TYPES_LINE_H
#define CGAL_KERNELD_TYPES_LINE_H
#include <utility>
#include <CGAL/NewKernel_d/functor_tags.h>
#include <CGAL/Kernel/mpl.h>
namespace CGAL {
template <class R_> class Line {
	typedef typename Get_type<R_, FT_tag>::type FT_;
	typedef typename Get_type<R_, Point_tag>::type	Point_;
	typedef std::pair<Point_,Point_> Data_;
	Data_ data;
	public:
	Line(){}
	Line(Point_ const&a, Point_ const&b): data(a,b) {}
	Point_ point(int i)const{
	  if(i==0) return data.first;
	  if(i==1) return data.second;
	  throw "not implemented";
	}
	Line opposite()const{
		return Line(data.second,data.first);
	}
};
namespace CartesianDKernelFunctors {
  template <class R_> struct Construct_line : Store_kernel<R_> {
    CGAL_FUNCTOR_INIT_STORE(Construct_line)
    typedef typename Get_type<R_, Line_tag>::type	result_type;
    typedef typename Get_type<R_, Point_tag>::type	Point;
    typedef typename Get_type<R_, Vector_tag>::type	Vector;
    typedef typename Get_functor<R_, Translated_point_tag>::type Tp_;
    //typedef typename Get_functor<R_, Difference_of_points_tag>::type Dp_;
    //typedef typename Get_functor<R_, Scaled_vector_tag>::type Sv_;
    result_type operator()(Point const&a, Point const&b)const{
      return result_type(a,b);
    }
    result_type operator()(Point const&a, typename First_if_different<Vector,Point>::Type const&b)const{
      Tp_ tp(this->kernel());
      return result_type(a,tp(a,b));
    }
  };
}
CGAL_KD_DEFAULT_TYPE(Line_tag,(CGAL::Line<K>),(Point_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Construct_ttag<Line_tag>,(CartesianDKernelFunctors::Construct_line<K>),(Line_tag,Point_tag,Vector_tag),(Translated_point_tag));

} // namespace CGAL

#endif // CGAL_KERNELD_TYPES_LINE_H
