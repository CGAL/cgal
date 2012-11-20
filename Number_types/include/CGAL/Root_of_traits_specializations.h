// Copyright (c) 2010 GeometryFactory
// Copyright (c) 1999-2007  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Michael Hemmer
//                 Sebastien Loriot
//                 Sylvain Pion

#ifndef CGAL_ROOT_OF_TRAITS_SPECIALIZATIONS_H
#define CGAL_ROOT_OF_TRAITS_SPECIALIZATIONS_H

#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Root_of_traits.h>

namespace CGAL {

// We create a type of new node in Lazy_exact_nt's DAG
// for the make_root_of_2() and solve_1(of degree 2) operation.

template <typename ET >
struct Lazy_exact_ro2
  : public Lazy_exact_nt_rep< typename Root_of_traits<ET>::Root_of_2 >
{
    typedef typename Root_of_traits<ET>::Root_of_2   RO2;
    typedef Lazy_exact_nt_rep<RO2>                  Base;
    typedef typename Base::AT::Protector            P;


    mutable Lazy_exact_nt<ET> op1, op2, op3;
    bool smaller;
    bool old_rep;//if the rep=true then representation with polynomial coeff, else alpha, beta, gamma


    Lazy_exact_ro2 (const Lazy_exact_nt<ET> &a,
                    const Lazy_exact_nt<ET> &b,
                    const Lazy_exact_nt<ET> &c, bool s)
      : Base((P(), make_root_of_2(a.approx(), b.approx(), c.approx(), s))),
        op1(a), op2(b), op3(c), smaller(s), old_rep(true) {}

    Lazy_exact_ro2 (const Lazy_exact_nt<ET> &a,
                    const Lazy_exact_nt<ET> &b,
                    const Lazy_exact_nt<ET> &c)
      : Base((P(), make_root_of_2(a.approx(), b.approx(), c.approx()))),
        op1(a), op2(b), op3(c), smaller(true), old_rep(false) {}

    void update_exact() const
    {
        if (old_rep)
          this->et = new RO2(make_root_of_2(op1.exact(), op2.exact(),
                                            op3.exact(), smaller));
        else
          this->et = new RO2(make_root_of_2(op1.exact(), op2.exact(),
                                            op3.exact()));
        if (!this->approx().is_point())
            this->at = to_interval(*(this->et));
        this->prune_dag();

    }

    void prune_dag() const
    {
        op1 = op2 = op3 = Lazy_exact_nt<ET>();
    }
};

template <typename NT >
struct Root_of_traits< Lazy_exact_nt < NT > >
{
private:
    typedef Root_of_traits<NT> T;
public:
    typedef Root_of_traits< Lazy_exact_nt < NT > > Base;
    typedef Lazy_exact_nt< typename T::Root_of_1 > Root_of_1;
    typedef Lazy_exact_nt< typename T::Root_of_2 > Root_of_2;
    typedef Root_of_2 RootOf_2;
    typedef Root_of_1 RootOf_1;

    struct Make_root_of_2{
        typedef Root_of_2 result_type;
        Root_of_2
        operator()(const Lazy_exact_nt<NT>& a, const Lazy_exact_nt<NT>& b, const Lazy_exact_nt<NT>& c) const {
            return new Lazy_exact_ro2<NT>(a, b, c);
        };
        Root_of_2
        operator()(const Lazy_exact_nt<NT>& a, const Lazy_exact_nt<NT>& b, const Lazy_exact_nt<NT>& c, bool smaller) const{
          return new Lazy_exact_ro2<NT>(a, b, c, smaller);
        };
    };
  
private:
  typedef CGAL::Algebraic_structure_traits<Root_of_2> AST;
public:
  typedef typename AST::Square  Square; 
  typedef typename AST::Inverse Inverse;

  struct Make_sqrt{
    typedef Root_of_2 result_type;
    Root_of_2 operator()(const Lazy_exact_nt<NT>& x) const {
      return new Lazy_exact_ro2<NT>( Lazy_exact_nt<NT>(0), Lazy_exact_nt<NT>(1) , x);
    }
  };
};


//these two functions for test suite requirement
template < typename RT >
typename CGAL::Root_of_traits<CGAL::Lazy_exact_nt<RT> >::Root_of_2 make_sqrt(const CGAL::Lazy_exact_nt< RT> & r)
{
  typedef Lazy_exact_nt< RT> TT;
  CGAL_assertion(r >= 0);
  if(CGAL_NTS is_zero(r)) return make_root_of_2((TT) 1,(TT) 0,(TT) 0);
  return make_root_of_2((TT) 1,(TT) 0,-r,false);
}

} //namespace CGAL

#endif // CGAL_ROOT_OF_TRAITS_SPECIALIZATIONS_H
