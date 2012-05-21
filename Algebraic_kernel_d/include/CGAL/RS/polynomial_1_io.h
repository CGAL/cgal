// Copyright (c) 2006-2009 Inria Lorraine (France). All rights reserved.
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

#ifndef CGAL_RS_POLYNOMIAL_1_IO_H
#define CGAL_RS_POLYNOMIAL_1_IO_H

namespace CGAL{

inline
std::ostream& operator<<(std::ostream &os,const RS_polynomial_1 &p){
        if(is_pretty(os)){
                bool printed = false;
                int degree=p.get_degree();
                mpz_t *coef=p.get_coefs();
                if(!degree)
                        return(os<<Gmpz(coef[0]));
                for(int i=degree;i>=0;--i){
                        if(mpz_sgn(coef[i])){
                                if(printed&&(mpz_sgn(coef[i])==1))
                                        os<<"+";
                                printed = true;
                                bool flag = false;
                                if((!mpz_cmp_si(coef[i],-1))&&i)
                                        os<<"-";
                                else
                                        if((mpz_cmp_ui(coef[i],1))||(!i)){
                                                flag = true;
                                                os<<Gmpz(coef[i]);
                                        }
                                if(i){
                                        if (flag)
                                                os<<"*";
                                        os<<"x";
                                        if(i-1)
                                                os<<"^"<<i;
                                }
                        }
                }
                if (!printed)
                        os<<'0';
                return os;
        }else{
                int d=p.get_degree();
                os<<"[ d="<<d<<" ";
                os<<"[ ";
                for(int i=0;i<d+1;++i)
                        os<<Gmpz(p.coef(i))<<" ";
                os<<"] ]";
                return os;
        }
}

inline
std::istream& operator>>(std::istream &is,RS_polynomial_1 &pol){
        std::istream::int_type c;
        std::ios::fmtflags old_flags=is.flags();
        is.unsetf(std::ios::skipws);
        gmpz_eat_white_space(is);
        if(is_pretty(is)){
                std::string s;
                while(((c=is.get())>='0'&&c<='9')||
                      c=='+'||c=='-'||c=='*'||c==' '||c=='x'||c=='^')
                        s.push_back(c);
                is.putback(c);
                pol=RS_polynomial_1(s);
                goto ret_is;
        }else{
                int deg;
                Gmpz z;
                c=is.get();     // eat a '['
                if(c!='[')
                        goto ret_fail_is;
                gmpz_eat_white_space(is);
                c=is.get();     // c='d'
                if(c!='d')
                        goto ret_fail_is;
                gmpz_eat_white_space(is);
                c=is.get();     // c='='
                if(c!='=')
                        goto ret_fail_is;
                gmpz_eat_white_space(is);
                c=is.peek();
                deg=0;
                while(c>='0'&&c<='9'){
                        c=is.get();
                        deg=10*deg+(c-'0');
                        c=is.peek();
                }
                pol=RS_polynomial_1(0);
                pol.set_degree(deg);
                gmpz_eat_white_space(is);
                c=is.get();
                if(c!='[')
                        goto ret_fail_is;
                gmpz_eat_white_space(is);
                for(int k=0;k<deg+1;++k){
                        is>>z;
                        pol.set_coef(k,z);
                        gmpz_eat_white_space(is);
                }
                c=is.get();     // ']'
                if(c!=']')
                        goto ret_fail_is;
                gmpz_eat_white_space(is);
                c=is.get();     // ']'
                if(c!=']')
                        goto ret_fail_is;
                goto ret_is;
        }
ret_fail_is:
        is.setstate(std::ios_base::failbit);
ret_is:
        is.flags(old_flags);
        return is;
}

} // namespace CGAL

#endif  // CGAL_RS_POLYNOMIAL_1_IO_H
