// Copyright (c) 2010  INRIA Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial

#include <CGAL/basic.h>

namespace CGAL{
namespace Test{

/// Force the double (non-extended) precision, even when the x87 Intel FPU
/// is used.
///
/// Applied to another number type, does nothing.
///@{
template <typename T> T ftd(const T x) { return x; }
double ftd(const double x) { return CGAL_IA_FORCE_TO_DOUBLE(x); }
///@}

template <class T, class RootOf1, class RootOf2>
void test_root_of_traits(){
    // pure type checking
    typedef CGAL::Root_of_traits<T> RoT;
    typedef typename RoT::Root_of_1 Root_of_1;
    typedef typename RoT::Root_of_2 Root_of_2;

    CGAL_static_assertion((::boost::is_same<RootOf1,Root_of_1>::value));
    CGAL_static_assertion((::boost::is_same<RootOf2,Root_of_2>::value));

    typedef typename RoT::Make_root_of_2 Make_root_of_2;
    typedef typename RoT::Make_sqrt      Make_sqrt;
    typedef typename RoT::Inverse        Inverse;
    typedef typename RoT::Square         Square;

    const Make_root_of_2& make_root_of_2 = Make_root_of_2();
    const Make_sqrt&      make_sqrt      = Make_sqrt();
    const Inverse&        inverse        = Inverse();
    const Square&         square         = Square();

    CGAL_static_assertion((::boost::is_same<Root_of_2,typename Make_root_of_2::result_type>::value));
    CGAL_static_assertion((::boost::is_same<Root_of_2,typename Make_sqrt::result_type>::value));
    CGAL_static_assertion((::boost::is_same<Root_of_2,typename Inverse::result_type>::value));
    CGAL_static_assertion((::boost::is_same<Root_of_2,typename Square::result_type>::value));


    {
      Root_of_2 r  = ftd(make_root_of_2(T(0),T(-1),T(2))); //-sqrt(2)
      Root_of_2 rl = ftd(make_root_of_2(T(1),T(0),T(-2),true)); //-sqrt(2);
      Root_of_2 rr = ftd(make_root_of_2(T(1),T(0),T(-2),false)); //+sqrt(2)
      assert(r == rl);
      assert(rl != rr);

      assert( ftd(r * Root_of_1(2)) == ftd(make_root_of_2(T(0),T(-2),T(2))));
      assert( ftd(r * T(2)) == ftd(make_root_of_2(T(0),T(-2),T(2))));
    }{
      Root_of_2 r  = ftd(CGAL::make_root_of_2(T(0),T(-1),T(2))); //-sqrt(2)
      Root_of_2 rl = ftd(CGAL::make_root_of_2(T(1),T(0),T(-2),true)); //-sqrt(2);
      Root_of_2 rr = ftd(CGAL::make_root_of_2(T(1),T(0),T(-2),false)); //+sqrt(2)
      assert(r == rl);
      assert(rl != rr);

      assert( ftd(r * Root_of_1(2)) == ftd(CGAL::make_root_of_2(T(0),T(-2),T(2))));
      assert( ftd(r * T(2)) == ftd(CGAL::make_root_of_2(T(0),T(-2),T(2))));
    }


    {
      Root_of_2 r  = ftd(make_sqrt(T(2))); //sqrt(2)
      Root_of_2 rr = ftd(make_root_of_2(T(1),T(0),T(-2),false)); //+sqrt(2)
      assert(r == rr);
    }{
      Root_of_2 r  = ftd(CGAL::make_sqrt(T(2))); //sqrt(2)
      Root_of_2 rr = ftd(CGAL::make_root_of_2(T(1),T(0),T(-2),false)); //+sqrt(2)
      assert(r == rr);
    }

    {
      Root_of_2 r  = ftd(inverse(ftd(CGAL::make_sqrt(T(2)))));
      Root_of_2 rr = ftd(1/ftd(CGAL::make_sqrt(T(2))));
      assert(r == rr);
    }{
        Root_of_2 r  = ftd(CGAL::inverse(ftd(CGAL::make_sqrt(T(2)))));
        Root_of_2 rr = ftd(1/ftd(CGAL::make_sqrt(T(2))));
        assert(r == rr);
    }

    {
      Root_of_2 r  = ftd(square(ftd(CGAL::make_sqrt(T(2)))));
      Root_of_2 rr = ftd(ftd(CGAL::make_sqrt(T(2)))*ftd(CGAL::make_sqrt(T(2))));
      assert(r == rr);
    }{
      Root_of_2 r  = ftd(CGAL::square(ftd(CGAL::make_sqrt(T(2)))));
      Root_of_2 rr = ftd(ftd(CGAL::make_sqrt(T(2)))*ftd(CGAL::make_sqrt(T(2))));
      assert(r == rr);
    }

    bool is_not_exact = !CGAL::Algebraic_structure_traits<T>::Is_exact::value;
    {
      std::vector<Root_of_2> roots;
      CGAL::compute_roots_of_2(T(1),T(0),T(-2),std::back_inserter(roots));
      assert(roots.size()==2);
      assert(roots[0]==-CGAL::make_sqrt(T(2)) || is_not_exact );
      assert(roots[1]== CGAL::make_sqrt(T(2)) || is_not_exact );
    }
    {
      Root_of_2 roots[2]= {Root_of_2(1),Root_of_2(1)};
      CGAL::compute_roots_of_2(T(13),T(4),T(-23),roots);
      assert(roots[0]==CGAL::make_root_of_2(T(13),T(4),T(-23),true)  || is_not_exact );
      assert(roots[1]==CGAL::make_root_of_2(T(13),T(4),T(-23),false) || is_not_exact );
    }
    {
      std::vector<Root_of_2> roots;
      CGAL::compute_roots_of_2(T(1),T(-6),T(9),std::back_inserter(roots));
      assert(roots.size()==1);
      assert(roots[0]==Root_of_2(3) || is_not_exact );
    }
    {
      std::vector<Root_of_2> roots;
      CGAL::compute_roots_of_2(T(1),T(0),T(2),std::back_inserter(roots));
      assert(roots.size()==0);
    }
    {
      std::vector<Root_of_2> roots;
      CGAL::compute_roots_of_2(T(0),T(2),T(3),std::back_inserter(roots));
      assert(roots.size()==1);
      assert(roots[0]==-Root_of_2(3)/Root_of_2(2) || is_not_exact );
    }

}

} // namespace Test
} // namespace CGAL
