// Copyright (c) 2022 Institut Géographique National - IGN (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mathieu Brédif

#ifndef CGAL_BBOX_D_H
#define CGAL_BBOX_D_H

#include <boost/config.hpp> // defines BOOST_PREVENT_MACRO_SUBSTITUTION
#include <stddef.h>
#include <limits>
#include <array>
#include <iostream>
#include <iterator>
#include <CGAL/use.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/assertions.h>
#include <CGAL/Dimension.h>
#include <CGAL/Concatenate_iterator.h>
#include <boost/math/special_functions/next.hpp>

namespace CGAL {
namespace Impl {

template<typename Container, typename Derived>
class Bbox
{
protected:
    typedef typename Container::value_type T;
    Container min_values;
    Container max_values;

public:
    Bbox& operator+=(const Bbox& bbox)
    {
        CGAL_assertion(min_values.size() == 0 || min_values.size() == bbox.min_values.size());
        if(min_values.size() == 0){
            *this =  bbox;
        }
        int dim = bbox.min_values.size();
        for(int i=0; i<dim; ++i)
        {
            if(min_values[i] > bbox.min_values[i]) min_values[i] = bbox.min_values[i];
            if(max_values[i] < bbox.max_values[i]) max_values[i] = bbox.max_values[i];
        }
        return *this;
    }

    inline int dimension() const
    {
        return static_cast<const Derived*>(this)->dimension();
    }

    inline T min BOOST_PREVENT_MACRO_SUBSTITUTION (int i) const
    {
        return min_values[i];
    }

    inline T max BOOST_PREVENT_MACRO_SUBSTITUTION (int i) const
    {
        return max_values[i];
    }

    inline T& min BOOST_PREVENT_MACRO_SUBSTITUTION (int i)
    {
        return min_values[i];
    }

    inline T& max BOOST_PREVENT_MACRO_SUBSTITUTION (int i)
    {
        return max_values[i];
    }

    inline T measure() const {
        T result = max_values[0] - min_values[0];
        if (result <= 0) return 0;
        for(int i=1; i<dimension(); ++i)
            result *= max_values[i] - min_values[i];
        return result;
    }

    inline T intersection_measure(const Bbox& bbox) const {
        CGAL_assertion(dimension() == bbox.dimension());
        T result = 1;
        for(int i=0; i<dimension(); ++i) {
            result *= (std::min)((max)(i), (bbox.max)(i)) -
                      (std::max)((min)(i), (bbox.min)(i));
            if (result <= 0) return 0;
        }
        return result;
    }

    inline
    void
    scale(double factor)
    {
      CGAL_precondition(factor > 0);

      if (factor == 1.)
        return;

      int d = dimension();
      for(int i=0; i<d; ++i){
        auto half_width = ((max)(i) - (min)(i)) * 0.5;
        auto center = (min)(i) + half_width;
        min_values[i] = center - factor * half_width;
        max_values[i] = center + factor * half_width;
      }
    }

    void dilate(int dist)
    {
      using boost::math::float_advance;
      int d = dimension();
      for(int i=0; i<d; ++i){
        min_values[i] = float_advance(min_values[i],-dist);
        max_values[i] = float_advance(max_values[i], dist);
      }
    }



    bool operator==(const Bbox& bbox) const {
        for(int i=0; i<dimension(); ++i)
            if (min_values[i] != bbox.min_values[i] || max_values[i] != bbox.max_values[i])
                return false;
        return true;
    }

    bool operator!=(const Bbox& bbox) const { return !operator==(bbox); }

protected:
    void init(int d, T range = -std::numeric_limits<T>::infinity()) {
        for(int i=0; i<d; ++i)
        {
            min_values[i] = -range;
            max_values[i] =  range;
        }
    }

    template <typename I>
    void init(int d, I b, I e) {
        CGAL_USE(e);
        CGAL_assertion(d == std::distance(b,e));
        for(int i=0; i<d; ++i,++b)
        {
            min_values[i] = (*b).first;
            max_values[i] = (*b).second;
        }
    }

    void init(const Bbox_2& bb){
        min_values[0] = bb.xmin();
        min_values[1] = bb.ymin();
        max_values[0] = bb.xmax();
        max_values[1] = bb.ymax();
    }

    void init(const Bbox_3& bb){
        min_values[0] = bb.xmin();
        min_values[1] = bb.ymin();
        min_values[2] = bb.zmin();
        max_values[0] = bb.xmax();
        max_values[1] = bb.ymax();
        max_values[2] = bb.zmax();
    }
};

}

template <typename Di>
class Bbox_d;

// A fixed D-dimensional axis aligned box
template<int N>
class Bbox_d<Dimension_tag<N>> : public Impl::Bbox<std::array<double, N>, Bbox_d<Dimension_tag<N>>>
{
    enum { D = N };
    using array_const_iterator = typename std::array<double, N>::const_iterator;
public:
    using Cartesian_const_iterator = Concatenate_iterator<array_const_iterator,array_const_iterator>;

    inline constexpr int dimension() const { return D; }
    Bbox_d() { this->init(N      ); }
    Bbox_d(double range) { this->init(N, range); }
    template <typename I>
    Bbox_d(I b, I e) { this->init(N, b, e); }

    Bbox_d(const Bbox_2& bb2){ this->init(bb2);}
    Bbox_d(const Bbox_3& bb3){ this->init(bb3);}

    Cartesian_const_iterator cartesian_begin() const
    {
        return Cartesian_const_iterator(this->min_values.end(), this->max_values.begin(), this->min_values.begin());
    }

    Cartesian_const_iterator cartesian_end() const
    {
        return Cartesian_const_iterator(this->min_values.end(), this->max_values.begin(), this->max_values.end(),0);
    }
};

// A dynamic D-dimensional axis aligned box
template<>
class Bbox_d<Dynamic_dimension_tag> : public Impl::Bbox<std::vector<double>, Bbox_d<Dynamic_dimension_tag>>
{
public:
    inline int dimension() const { return this->min_values.size(); }
    Bbox_d(int d = 0          ) { init_values(d); this->init(d       ); }
    Bbox_d(int d, double range) { init_values(d); this->init(d, range); }
    template <typename I>
    Bbox_d(int d, I b, I e) { init_values(d); this->init(d, b, e); }

protected:
    void init_values(int d) {
        this->min_values.resize(d);
        this->max_values.resize(d);
    }
};

template<typename Container, typename Derived>
std::ostream& operator<<(std::ostream& out, const Impl::Bbox<Container, Derived>& bbox)
{
    int d = bbox.dimension();
    for(int i=0; i<d; ++i)
        out  << (bbox.min)(i) << "  " << (bbox.max)(i) << " ";
    return out;
}

template<typename Container, typename Derived>
std::istream& operator>>(std::istream& in, Impl::Bbox<Container, Derived>& bbox)
{
    int d = bbox.dimension();
    for(int i=0; i<d; ++i)
        in  >> (bbox.min)(i) >> (bbox.max)(i);
    return in;
}

template<int N>
Bbox_d<Dimension_tag<N>> operator+(Bbox_d<Dimension_tag<N>> bbox, const Bbox_d<Dimension_tag<N>>& other)
{
    bbox += other;
    return bbox;
}

template<typename Di>
inline
bool
do_overlap(const Bbox_d<Di>& bb1, const Bbox_d<Di>& bb2)
{
    // check for emptiness ??
    int d = bb1.dimension();
    for(int i=0; i<d; ++i){
    if ((bb1.max)(i) < (bb2.min)(i) || (bb2.max)(i) < (bb1.min)(i))
        return false;
    }
    return true;
}



} // namespace CGAL

#endif // CGAL_BBOX_D_H
