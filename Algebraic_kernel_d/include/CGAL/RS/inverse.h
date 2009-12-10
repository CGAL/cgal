// Copyright (c) 2007-2008 Inria Lorraine (France). All rights reserved.
// 
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author: Luis Peñaranda <luis.penaranda@loria.fr>

#ifndef CGAL_RS__INVERSE_H
#define CGAL_RS__INVERSE_H

#include <stdint.h>

#define n_(A)   (A<0?-A:A)
//#define u_(A) (A<0?-1:1)
#define u_(A,C) (A<0?(C<0?1:-1):(C<0?-1:1))

namespace CGAL{
namespace RS_MGCD{

class Inverse{

    protected:
        // given a and b, returns s such that gcd(a,b)=s*a+t*b (GCL, page 36)
        // s*a+t*q=1 => s is the inverse of a, mod q (pafe 173)
        static int64_t eea_s(uint32_t a,uint32_t b){
            int64_t c1,d1,r1;//,c2,d2,r2,t,s;
            uint32_t r,c,d;//,q;
            // a and b are positive, so a=n_(a) and b=n_(b)
            //c=n_(a);      d=n_(b);
            c=a;    d=b;
            c1=1;   d1=0;
            //c2=0; d2=1;
            while(d){
                    //q=c/d;
                    r=c%d;
                    r1=c1-d1*(c/d); //r2=c2-d2*q;
                    c=d;    c1=d1;  //c2=d2;
                    d=r;    d1=r1;  //d2=r2;
            }
            // gcd(a,b) is n_(c)
            //t=c2/(_u(b)*u_(c));
            // a and c are always positive, so s=c1/u_(a,c) equals c1
            //s=c1/u_(a,c);
            //return s;
            return c1;
        };
}; // class Inverse

} // namespace RS_MGCD
} // namespace CGAL

#endif  // CGAL_RS__INVERSE_H

// vim: tabstop=4: softtabstop=4: smarttab: shiftwidth=4: expandtab
