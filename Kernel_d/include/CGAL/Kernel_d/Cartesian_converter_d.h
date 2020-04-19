// Copyright (c) 2001-2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Menelaos Karavelas <mkaravel@cse.nd.edu>
//                 Samuel Hornus <Samuel.Hornus@sophia.inria.fr>

#ifndef CGAL_CARTESIAN_CONVERTER_D_H
#define CGAL_CARTESIAN_CONVERTER_D_H

// This file contains the definition of a dD kernel converter, based on Cartesian
// representation.  It should work between *Cartesian_d<A> and *Cartesian_d<B>,
// provided you give a NT converter from A to B.
// TODO: There's NO Homogeneous counterpart at the moment.
// TODO: Maybe we can derive Cartesian_converter_d from Cartesian_converter ?
//       But I don't know if it is interesting or not.

#include <CGAL/basic.h>
#include <CGAL/NT_converter.h>
#include <CGAL/Enum_converter.h>
#include <CGAL/Iterator_transform.h>
#include <vector>

#include <CGAL/Referenced_argument.h>

namespace CGAL {

// Guess which compiler needs this work around ?
namespace internal {
template < typename K1, typename K2 >
struct Default_converter_d {
  typedef typename K1::FT FT1;
  typedef typename K2::FT FT2;
  typedef ::CGAL::NT_converter<FT1, FT2> Type;
};
} // namespace internal

template < class K1, class K2,
          //class Converter = NT_converter<typename K1::FT, typename K2::FT> >
           class Converter = typename internal::Default_converter_d<K1, K2>::Type >
class Cartesian_converter_d : public Enum_converter
{
    typedef Enum_converter   Base;
    typedef Cartesian_converter_d<K1, K2, Converter>        Self;

public:
    typedef K1         Source_kernel;
    typedef K2         Target_kernel;
    typedef Converter  Number_type_converter;

#ifdef CGAL_CFG_USING_BASE_MEMBER_BUG
    bool operator()(bool b) const { return Base::operator()(b); }
    Sign operator()(Sign s) const { return Base::operator()(s); }

    Oriented_side operator()(Oriented_side os) const {
      return Base::operator()(os);
    }

    Bounded_side operator()(Bounded_side bs) const {
      return Base::operator()(bs);
    }

    Comparison_result operator()(Comparison_result cr) const {
      return Base::operator()(cr);
    }

    Angle operator()(Angle a) const { return Base::operator()(a); }
#else
    using Base::operator();
#endif

    Cartesian_converter_d() // To shut up a warning with SunPRO.
        : c(), k(), result_point_(20) {}

    Origin
    operator()(const Origin& o) const
    {
        return o;
    }

    Null_vector
    operator()(const Null_vector& n) const
    {
        return n;
    }

    typename K2::FT
    operator()(const typename K1::FT &a) const
    {
        return c(a);
    }

    std::vector<Object>
    operator()(const std::vector<Object>& v) const
    {
      std::vector<Object> res;
      res.reserve(v.size());
      for(unsigned int i = 0; i < v.size(); i++) {
        res.push_back(operator()(v[i]));
      }
      return res;
    }

    int operator()(const int &a)
    {
        return a;
    }

    typename K2::Point_d
    operator()(const typename K1::Point_d &a) const
    {
        typedef typename K2::Point_d Point_d;
        typedef typename K1::Point_d::Cartesian_const_iterator Coord_iter;
        typedef Iterator_transform<Coord_iter, Converter> It;
        return Point_d(a.dimension(), It(a.cartesian_begin()), It(a.cartesian_end()));
    }

    typename K2::Vector_d
    operator()(const typename K1::Vector_d &a) const
    {
        typedef typename K2::Vector_d  Vector_d;
        typedef typename K1::Point_d::Cartesian_const_iterator Coord_iter;
        typedef Iterator_transform<Coord_iter, Converter> It;
        return Vector_d(a.dimension(), It(a.cartesian_begin()), It(a.cartesian_end()));
    }

    typename K2::Direction_d
    operator()(const typename K1::Direction_d &a) const
    {
        typedef typename K2::Direction_d  Direction_d;
        return Direction_d(operator()(a.vector()));
    }

    typename K2::Segment_d
    operator()(const typename K1::Segment_d &a) const
    {
        typedef typename K2::Segment_d  Segment_d;
        return Segment_d(operator()(a.source()), operator()(a.target()));
    }

    typename K2::Line_d
    operator()(const typename K1::Line_d &a) const
    {
        typedef typename K2::Line_d Line_d;
        return Line_d(operator()(a.point(0)), operator()(a.direction()));
    }

    typename K2::Ray_d
    operator()(const typename K1::Ray_d &a) const
    {
        typedef typename K2::Ray_d  Ray_d;
        return Ray_d(operator()(a.source()), operator()(a.second_point()));
    }

    typename K2::Sphere_d
    operator()(const typename K1::Sphere_d &a) const
    {
        typedef typename K2::Sphere_d        Sphere_d;
        typedef typename K1::Sphere_d::point_iterator Coord_iter;
        // TODO: Check that the use of Iterator_transform is correct
        typedef Iterator_transform<Coord_iter, Converter> It;
        return Sphere_d(a.dimension(), It(a.points_begin()), It(a.points_end()));
    }

    typename K2::Hyperplane_d
    operator()(const typename K1::Hyperplane_d &a) const
    {
        typedef typename K2::Hyperplane_d        Hyperplane_d;
        typedef typename K1::Hyperplane_d::Coefficient_const_iterator        Coord_iter;
        // TODO: Check that the use of Iterator_transform is correct
        typedef Iterator_transform<Coord_iter, Converter> It;
        return Hyperplane_d(a.dimension(), It(a.coefficients_begin()), It(a.coefficients_end()));
    }

    typename K2::Iso_box_d
    operator()(const typename K1::Iso_box_d &a) const
    {
        typedef typename K2::Iso_box_d  Iso_box_d;
        return Iso_box_d(operator()((a.min)()), operator()((a.max)()));
    }

    /*std::vector<int> &
    operator()(const std::vector<int> & a) const
    {
        return const_cast<std::vector<int> &>(a);
    }*/
    std::vector<int>::iterator
    operator()(const std::vector<int>::iterator & a) const
    {
        return a;
    }

    template<typename T>
    Referenced_argument<T> & operator()(const Referenced_argument<T> & a) const
    {
        return const_cast<Referenced_argument<T> &> (a);
    }

    // ---- The following allows to convert iterators to Point_d, Vector_d,
    // etc...  by adapting the iterator with an 'on-the-fly' converter. This is
    // not the most efficient because, if the iterator is dereferenced N times,
    // its value will be converted N times... TODO: Convert once and for all:
    // this should probably be done in Filtered_kernel.h (or a new
    // Filtered_kernel_d.h file).
    template<typename IterBase>
    struct Specialized_converter : public Self
    {
        typedef typename IterBase::value_type                argument_type;
        typedef typename argument_type::
                    template WithAnotherKernel<K2>::Type        result_type;
    };

    template<typename IterBase>
    Iterator_transform<IterBase, Specialized_converter<IterBase> >
    operator()(const IterBase & a) const
    {
        return Iterator_transform<IterBase, Specialized_converter<IterBase> >(a);
    }
    // ---- ---- ---- ---- ---- ---- //

private:
    Converter c;
    K2 k;
        typename K2::Point_d result_point_;
};

// Specialization when converting to the same kernel,
// to avoid making copies.
template < class K, class C >
class Cartesian_converter_d <K, K, C>
{
public:
  typedef K Source_kernel;
  typedef K Target_kernel;
  typedef C Number_type_converter;

  template < typename T >
  const T& operator()(const T&t) const { return t; }
};

} //namespace CGAL

#endif // CGAL_CARTESIAN_CONVERTER_D_H
