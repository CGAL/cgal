// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Monique Teillaud, Sylvain Pion, Athanasios Kakargias

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_ROOT_OF_CGAL_LAZY_EXACT_NT_H
#define CGAL_ROOT_OF_CGAL_LAZY_EXACT_NT_H

#include <CGAL/Binary_operator_result.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/NT_extensions_Root_of/CGAL_Interval_nt.h>
#include <CGAL/Binary_operator_result.h>


#if 1

CGAL_BEGIN_NAMESPACE



// We create a type of new node in Lazy_exact_nt's DAG
// for the make_root_of_2() operation.

template <typename ET >
struct Lazy_exact_ro2
  : public Lazy_exact_rep< typename Root_of_traits<ET>::RootOf_2 >
{
    typedef typename Root_of_traits<ET>::RootOf_2   RO2;
    typedef Lazy_exact_rep<RO2>                     Base;
    typedef typename Base::AT::Protector            P;


    mutable Lazy_exact_nt<ET> op1, op2, op3;
    bool smaller;

    Lazy_exact_ro2 (const Lazy_exact_nt<ET> &a,
		    const Lazy_exact_nt<ET> &b,
		    const Lazy_exact_nt<ET> &c, bool s)
#ifndef CGAL_CFG_COMMA_BUG
      : Base((P(), make_root_of_2(a.approx(), b.approx(), c.approx(), s))),
        op1(a), op2(b), op3(c), smaller(s) {}
#else
      : Base(a.approx() /* dummy value */, a),
	op1(a), op2(b), op3(c), smaller(s)
  {
    P p;
    this->approx() = make_root_of_2(a.approx(), b.approx(),
				    c.approx(), s);
  } 
#endif

    void update_exact()
    {
	this->et = new RO2(make_root_of_2(op1.exact(), op2.exact(),
					  op3.exact(), smaller));

	if (!this->approx().is_point())
            this->at = CGAL::to_interval(*(this->et));
	this->prune_dag();

    }

    void prune_dag() const
    {
	op1 = op2 = op3 = Lazy_exact_nt<ET>::zero();
    }
};

template < typename ET >
inline
Lazy_exact_nt< typename Root_of_traits<ET>::RootOf_2 >
make_root_of_2( const Lazy_exact_nt<ET> &a,
                const Lazy_exact_nt<ET> &b,
                const Lazy_exact_nt<ET> &c, bool d)
{
    return new Lazy_exact_ro2<ET>(a, b, c, d);
}

template <typename NT >
struct Root_of_traits< Lazy_exact_nt < NT > >
{
    typedef Root_of_traits<NT> T;
    typedef Lazy_exact_nt< typename T::RootOf_1 > RootOf_1;
    typedef Lazy_exact_nt< typename T::RootOf_2 > RootOf_2;
    typedef Lazy_exact_nt< typename T::RootOf_3 > RootOf_3;
    typedef Lazy_exact_nt< typename T::RootOf_4 > RootOf_4;
};

CGAL_END_NAMESPACE

#endif

#endif // CGAL_ROOT_OF_CGAL_LAZY_EXACT_NT_H
