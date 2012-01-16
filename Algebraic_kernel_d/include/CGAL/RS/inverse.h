// Copyright (c) 2007-2008 Inria Lorraine (France). All rights reserved.
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
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS__INVERSE_H
#define CGAL_RS__INVERSE_H

#ifdef _MSC_VER
#  define CGALRS_S64 __int64
#  define CGALRS_U64 unsigned __int64
#  define CGALRS_U32 unsigned __int32
#else
#  include <stdint.h>
#  define CGALRS_S64 int64_t
#  define CGALRS_U64 uint64_t
#  define CGALRS_U32 uint32_t
#endif

namespace CGAL{
namespace RS_MGCD{

#define CGALRS_N(A)   (A<0?-A:A)
//#define CGALRS_U(A) (A<0?-1:1)
#define CGALRS_U(A,C) (A<0?(C<0?1:-1):(C<0?-1:1))

class Inverse{

    protected:
        // given a and b, returns s such that gcd(a,b)=s*a+t*b (GCL, page 36)
        // s*a+t*q=1 => s is the inverse of a, mod q (pafe 173)
        static CGALRS_S64 eea_s(CGALRS_U32 a,CGALRS_U32 b){
            CGALRS_S64 c1,d1,r1;//,c2,d2,r2,t,s;
            CGALRS_U32 r,c,d;//,q;
            // a and b are positive, so a=CGALRS_N(a) and b=CGALRS_N(b)
            //c=CGALRS_N(a);      d=CGALRS_N(b);
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
            // gcd(a,b) is CGALRS_N(c)
            //t=c2/(CGALRS_U(b)*CGALRS_U(c));
            // a and c are always positive, so s=c1/CGALRS_U(a,c) equals c1
            //s=c1/CGALRS_U(a,c);
            //return s;
            return c1;
        };
}; // class Inverse

#undef CGALRS_N
#undef CGALRS_U

} // namespace RS_MGCD
} // namespace CGAL

#endif  // CGAL_RS__INVERSE_H
