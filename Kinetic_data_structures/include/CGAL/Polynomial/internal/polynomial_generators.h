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

#ifndef CGAL_POLYNOMIAL_GENERATORS_H
#define CGAL_POLYNOMIAL_GENERATORS_H

#include <CGAL/Polynomial/basic.h>

namespace CGAL { namespace POLYNOMIAL {

template <class K>
class Chebychev_generator
{
    public:
        Chebychev_generator(const K &k): k_(k){}
        typedef typename K::Function result_type;
        typedef unsigned int argument_type;
        result_type operator()(argument_type n) const
        {

            typename K::Construct_function cf= k_.construct_function_object();
            if ( n == 0 ) { return result_type(0); }

            if ( n == 1 ) { return cf(0, 1); }

            if ( n == 2 ) { return cf(-1, 0, 2); }

            int counter = static_cast<int>(n) - 2;

            result_type C2 = cf(0, 1);
            result_type C1 = cf(-1, 0, 2);
            result_type TWOX = cf(0,2);

            result_type Cnext;
            while ( counter > 0 ) {
                Cnext = TWOX * C1 - C2;
                C2 = C1;
                C1 = Cnext;
                counter--;
            }

            return C1;
        }
    protected:
        K k_;
};

//---------------------------------------------------------------

template <class K>
class Laguerre_generator
{
    public:
        Laguerre_generator(const K &k): k_(k){}
        typedef unsigned int argument_type;
        typedef typename K::Function result_type;
        result_type operator()(argument_type n) const
        {
            typedef typename result_type::NT   NT;
            typename K::Construct_function cf= k_.construct_function_object();

            if ( n == 0 ) { return result_type(1); }
            if ( n == 1 ) { return cf(1, -1); }

            int counter = static_cast<int>(n) - 1;
            int deg = 2;

            result_type L2(1);
            result_type L1 = cf(1, -1);

            result_type Lnext;
            while ( counter > 0 ) {
                result_type C1 = cf(NT(2 * deg - 1),-1);
                result_type C2 = cf(deg - 1);
                NT coef = NT(1)/NT(deg);

                Lnext = (C1 * L1 - C2 * L2) * coef;
                L2 = L1;
                L1 = Lnext;

                counter--;
                deg++;
            }

            return L1;
        }
    protected:
        K k_;
};

//---------------------------------------------------------------

template <class T>
class Wilkinson_generator
{
    public:
        Wilkinson_generator(const T &k): k_(k){}
        typedef unsigned int argument_type;
        typedef typename T::Function result_type;
        result_type operator()(unsigned int n) const
        {
            if ( n == 0 ) {
	      return result_type(typename result_type::NT(0));
            }

            result_type w(typename result_type::NT(1));
            typename  T::Construct_function cf = k_.construct_function_object();

            for (unsigned int i = 1; i <= n; i++) {
                int r = static_cast<int>(i);
                w = w * cf(-r, 1);
            }

            return w;
        }
    protected:
        T k_;
};

//---------------------------------------------------------------

template <class T>
class Mignotte_generator
{
    public:
        Mignotte_generator(const T &k): k_(k){}
        typedef unsigned int argument_type;
        typedef typename T::Function result_type;
        result_type operator()(argument_type n) const
        {
            typename T::Construct_function cf = k_.construct_function_object();

            if ( n < 3 ) {
                result_type m = cf(-2, 20, -50);
                return m;
            }
            else {
                std::vector<typename result_type::NT> v(n+1);
                v[0] = -2;
                v[1] = 20;
                v[2] = -50;
                for (unsigned int i = 3; i < n; i++) {
                    v[i] = 0;
                }
                v[n] = 1;
                result_type m(v.begin(), v.end());

                return m;
            }
        }
    protected:
        T k_;
};

} } //namespace CGAL::POLYNOMIAL
#endif                                            // CGAL_POLYNOMIAL_GENERATORS_H
