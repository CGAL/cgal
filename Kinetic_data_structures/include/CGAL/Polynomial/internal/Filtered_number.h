// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
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
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_INTERNAL_FILTERED_NUMBER_H
#define CGAL_POLYNOMIAL_INTERNAL_FILTERED_NUMBER_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/FPU.h>

namespace CGAL { namespace POLYNOMIAL {

template <class NT_t>
class Filtered_number
{
    typedef Filtered_number<NT_t> This;
    public:
        typedef NT_t NT;
        Filtered_number(NT_t nt):nt_(nt) {
            double_=false;
        }
        Filtered_number(double d):d_(d) {
            double_=true;
        }
        Filtered_number(){}
        void set(NT_t nt) {
            nt_=nt;
            double_=false;
        }
        void set(double d) {
            d_=d;
            if (!double_) {
                nt_=NT_t();
                double_=true;
            }
        }
        NT_t nt() const
        {
            if (double_) return NT_t(d_);
            return nt_;
        }
        double d() const
        {
            if (!double_) return CGAL_POLYNOMIAL_TO_DOUBLE(nt_);
            return d_;
        }
        const std::pair<double, double> it() const
        {
            if (!double_) return CGAL_POLYNOMIAL_TO_INTERVAL(nt_);
            return std::pair<double,double>(d_, d_);
        }
        This operator-() const
        {
            if (double_) {
                return This(-d_);
            }
            else {
                return This(-nt_);
            }
        }

        This inf() const
        {
            if ( double_ ) {
                double dbl = d_ - 1.0;
                if ( dbl < d_ && CGAL::is_finite(dbl) ) {
                    return This(dbl);
                }
                else {
                    return This( NT(d_) - NT(1) );
                }
            }

            return This(nt_ - NT(1));
        }

        This sup() const
        {
            if ( double_ ) {
                double dbl = d_ + 1.0;
                if ( dbl > d_ && CGAL::is_finite(dbl) ) {
                    return This(dbl);
                }
                else {
                    return This( NT(d_) + NT(1) );
                }
            }

            return This(nt_ + NT(1));
        }

        bool is_double() const { return double_; }

        bool operator<(const This &o) const
        {
            return compare(o) == CGAL::SMALLER;
        }
        bool operator>(const This &o) const
        {
            return compare(o) == CGAL::LARGER;
        }
        bool operator<=(const This &o) const
        {
            return compare(o) != CGAL::LARGER;
        }
        bool operator>=(const This &o) const
        {
            return compare(o) != CGAL::SMALLER;
        }
        bool operator==(const This &o) const
        {
            return compare(o) == CGAL::EQUAL;
        }
        bool operator!=(const This &o) const
        {
            return compare(o) != CGAL::EQUAL;
        }

        static This midpoint(const This &a, const This &b) {
            CGAL_precondition( a != b );

            if (a.double_ && b.double_) {
                double mid = (a.d_ + b.d_) * 0.5;
                double dmin, dmax;

// gcc in x86 has register with bigger size than memory
// as a result although we have underflow this is not
// detected; remedy: copy variable to memory and retrieve it;
// this is what the following macro does;
// many thanks to Sylvain for pointing this out and providing a
// solution
                mid = CGAL_IA_FORCE_TO_DOUBLE(mid);

                if ( a.d_ < b.d_ ) {
                    dmin = a.d_;
                    dmax = b.d_;
                }
                else {
                    dmin = b.d_;
                    dmax = a.d_;
                }
                if ( dmin < mid && mid < dmax ) {
                    return This(mid);
                }
                else {
                    return This(NT(.5)*(a.nt() + b.nt()));
                }
            }
            else {
                return This(NT(.5)*(a.nt() + b.nt()));
            }
        }
        double compute_double() const
        {
            return d();
        }
        std::pair<double, double> compute_interval() const
        {
            if (double_) {
                return std::pair<double, double>(d(), d());
            }
            else {
                return CGAL_POLYNOMIAL_TO_INTERVAL(nt());
            }
        }
        void write(std::ostream &out) const
        {
            if (double_) {
                out << d_;
            }
            else {
                out << nt_;
            }
        }
        template <class Functor>
        static typename Functor::result_type apply(const Functor &f, const This &a) {
            if (a.double_) {
                return f(a.d_);
            }
            else {
                return f(a.nt_);
            }
        }
        template <class Functor>
        static typename Functor::result_type apply(const Functor &f, const This &a,  const This &b) {
            if (a.double_ && b.double_) {
                return f(a.d_, b.d_);
            }
            else {
                return f(a.nt(), b.nt());
            }
        }

        template <class Functor, class Data>
            static typename Functor::result_type apply(const Functor &f, const This &a,  const This &b,
        const Data &da, const Data &db) {
            if (a.double_ && b.double_) {
                return f(a.d_, b.d_, da, db);
            }
            else {
                return f(a.nt(), b.nt(), da, db);
            }
        }

        template <class Functor, class Data>
            static typename Functor::result_type apply(const Functor &f,
            const This &a,
            const This &b,
        const Data& data) {
            if (a.double_ && b.double_) {
                return f(a.d_, b.d_, data);
            }
            else {
                return f(a.nt(), b.nt(), data);
            }
        }

    protected:
        NT nt_;
        double d_;
        bool double_;

        template <class NT>
        static CGAL::Comparison_result compare(const NT &a, const NT &b) {
            if (a < b) return CGAL::SMALLER;
            else if (a==b) return CGAL::EQUAL;
            else return CGAL::LARGER;
        }

        CGAL::Comparison_result compare(const This &o) const
        {
            if (double_ && o.double_) {
                return compare(d_, o.d_);
            }
            else {
                return compare(nt(), o.nt());
            }
        }
};

template <class NT>
Filtered_number<NT> midpoint(const Filtered_number<NT> &a, const Filtered_number<NT> &b)
{
    return Filtered_number<NT>::midpoint(a,b);
}


template <class NT, class Functor>
typename Functor::result_type apply(const Functor &f,  const Filtered_number<NT> &a)
{
    return Filtered_number<NT>::apply(f, a);
}


template <class NT, class Functor>
typename Functor::result_type apply(const Functor &f,  const Filtered_number<NT> &a,
const Filtered_number<NT> &b)
{
    return Filtered_number<NT>::apply(f, a, b);
}


template <class NT, class Functor, class Data>
typename Functor::result_type apply(const Functor &f,  const Filtered_number<NT> &a,
const Filtered_number<NT> &b, const Data &da, const Data &db)
{
    return Filtered_number<NT>::apply(f, a, b, da, db);
}


template <class NT>
double to_double(const Filtered_number<NT> &a)
{
    return a.to_double();
}


template <class NT>
std::pair<double, double> to_interval(const Filtered_number<NT> &a)
{
    return a.compute_interval();
}


template <class NT>
std::ostream &operator<<(std::ostream &out, const Filtered_number<NT> &a)
{
    a.write(out);
    return out;
}


} } //namespace CGAL::POLYNOMIAL

namespace CGAL {
template <class NT>
double to_double(const CGAL_POLYNOMIAL_NS::Filtered_number<NT> &a)
{
    return a.compute_double();
}


template <class NT>
std::pair<double, double> to_interval(const CGAL_POLYNOMIAL_NS::Filtered_number<NT> &a)
{
    return a.to_interval();
}


} //namespace CGAL

template <class NT>
std::ostream &operator<<(std::ostream &out, const CGAL_POLYNOMIAL_NS::Filtered_number<NT> &a)
{
    a.write(out);
    return out;
}
#endif
