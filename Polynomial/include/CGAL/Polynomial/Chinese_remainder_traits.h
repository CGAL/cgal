// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
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
// Author(s)     :  Michael Hemmer <hemmer@mpi-inf.mpg.de>

#ifndef CGAL_POLYNOMIAL_CHINESE_REMAINDER_TRAITS
#define CGAL_POLYNOMIAL_CHINESE_REMAINDER_TRAITS

#include <CGAL/basic.h>
#include <CGAL/Chinese_remainder_traits.h>
#include <CGAL/Polynomial/Polynomial_type.h>

namespace CGAL {

template <class NT> 
class Chinese_remainder_traits<Polynomial<NT> >{
public:
    typedef Polynomial<NT> Type;
    typedef Chinese_remainder_traits<NT> CRT_NT;
   
    typedef typename CRT_NT::Scalar_type Scalar_type;
    
    struct Chinese_remainder{
        void operator()(
                const Scalar_type& m1, const Scalar_type& m2, 
                const Scalar_type& m, 
                const Scalar_type& s,  const Scalar_type& t,  
                const Type& u1, const Type& u2, 
                Type& u) const {
            
            typename CRT_NT::Chinese_remainder chinese_remainder_nt;
            
            CGAL_precondition(u1.degree() == u2.degree());
            
            std::vector<NT> coeffs(u1.degree()+1);
            for(int i = 0; i <= u1.degree(); i++){
                NT c; 
                chinese_remainder_nt(m1,m2,m,s,t,u1[i],u2[i],c);
                coeffs[i] = c;  
            }
            u = Polynomial<NT>(coeffs.begin(),coeffs.end());
        }

        void operator()(
                const Scalar_type& m1, const Type& u1,
                const Scalar_type& m2, const Type& u2,
                Scalar_type& m, Type& u) const {
            Scalar_type s,t; 
            
            CGAL::extended_euclidean_algorithm(m1,m2,s,t);
            m = m1 * m2;
            this->operator()(m1,m2,m,s,t,u1,u2,u);
        }
    };
};


} // namespace CGAL 

#endif // CGAL_POLYNOMIAL_CHINESE_REMAINDER_TRAITS
