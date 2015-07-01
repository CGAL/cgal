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

#ifndef CGAL_KERNELD_TYPES_ISO_BOX_H
#define CGAL_KERNELD_TYPES_ISO_BOX_H
#include <utility>
#include <CGAL/basic.h>
#include <CGAL/NewKernel_d/functor_tags.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/transforming_pair_iterator.h>
namespace CGAL {
template <class R_> class Iso_box {
	typedef typename Get_type<R_, FT_tag>::type FT_;
	typedef typename Get_type<R_, Point_tag>::type	Point_;
	typedef std::pair<Point_,Point_> Data_;
	Data_ data;
	public:
	Iso_box(){}
	Iso_box(Point_ const&a, Point_ const&b): data(a,b) {}
	Point_ min BOOST_PREVENT_MACRO_SUBSTITUTION ()const{
	  return data.first;
	}
	Point_ max BOOST_PREVENT_MACRO_SUBSTITUTION ()const{
	  return data.second;
	}
};
namespace CartesianDKernelFunctors {
  template <class R_> struct Construct_iso_box : Store_kernel<R_> {
    CGAL_FUNCTOR_INIT_STORE(Construct_iso_box)
    typedef typename Get_type<R_, Iso_box_tag>::type	result_type;
    typedef typename Get_type<R_, RT_tag>::type RT;
    typedef typename Get_type<R_, Point_tag>::type	Point;
    typedef typename Get_functor<R_, Construct_ttag<Point_tag> >::type Cp_;
    typedef typename Get_functor<R_, Construct_ttag<Point_cartesian_const_iterator_tag> >::type Ci_;
    result_type operator()(Point const&a, Point const&b)const{
      Cp_ cp(this->kernel());
      Ci_ ci(this->kernel());
      return result_type(cp(
	  make_transforming_pair_iterator(ci(a,Begin_tag()), ci(b,Begin_tag()), Min<RT>()),
	  make_transforming_pair_iterator(ci(a,End_tag()), ci(b,End_tag()), Min<RT>())),
      cp(
	  make_transforming_pair_iterator(ci(a,Begin_tag()), ci(b,Begin_tag()), Max<RT>()),
	  make_transforming_pair_iterator(ci(a,End_tag()), ci(b,End_tag()), Max<RT>())));
    }
  };

  template <class R_> struct Construct_min_vertex {
    CGAL_FUNCTOR_INIT_IGNORE(Construct_min_vertex)
    typedef typename Get_type<R_, Iso_box_tag>::type argument_type;
    //TODO: make result_type a reference
    typedef typename Get_type<R_, Point_tag>::type result_type;
    result_type operator()(argument_type const&b)const{
      return b.min BOOST_PREVENT_MACRO_SUBSTITUTION ();
    }
  };
  template <class R_> struct Construct_max_vertex {
    CGAL_FUNCTOR_INIT_IGNORE(Construct_max_vertex)
    typedef typename Get_type<R_, Iso_box_tag>::type argument_type;
    typedef typename Get_type<R_, Point_tag>::type result_type;
    result_type operator()(argument_type const&b)const{
      return b.max BOOST_PREVENT_MACRO_SUBSTITUTION ();
    }
  };
}
//TODO (other types as well) only enable these functors if the Iso_box type is the one defined in this file...
CGAL_KD_DEFAULT_TYPE(Iso_box_tag,(CGAL::Iso_box<K>),(Point_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Construct_ttag<Iso_box_tag>,(CartesianDKernelFunctors::Construct_iso_box<K>),(Iso_box_tag,Point_tag),(Construct_ttag<Point_cartesian_const_iterator_tag>,Construct_ttag<Point_tag>));
CGAL_KD_DEFAULT_FUNCTOR(Construct_min_vertex_tag,(CartesianDKernelFunctors::Construct_min_vertex<K>),(Iso_box_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Construct_max_vertex_tag,(CartesianDKernelFunctors::Construct_max_vertex<K>),(Iso_box_tag),());
} // namespace CGAL

#endif // CGAL_KERNELD_TYPES_ISO_BOX_H
