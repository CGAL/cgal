// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
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
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/Polynomial/include/CGAL/Polynomial.h $
// $Id: Polynomial.h 47254 2008-12-06 21:18:27Z afabri $
// 
//
// Author(s)     :   Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//                   Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#include <CGAL/basic.h>
#include <cassert>

//#include <CGAL/_test_basic.h>

#ifndef CGAL_TEST_ALGEBRAIC_CURVE_KERNEL_2_H
#define CGAL_TEST_ALGEBRAIC_CURVE_KERNEL_2_H

namespace CGAL {

namespace internal {

static const char *ACK_2_ascii_polys[] = {

    "P[12(0,P[12(0,-403)(1,445)(2,-244)(3,-261)(4,182)(5,-284)(6,321)(7,-150)(8,-312)(9,-395)(10,451)(11,166)(12,-268)])(1,P[11(0,-149)(1,-194)(2,190)(3,397)(4,-423)(5,507)(6,-321)(7,215)(8,-459)(9,-112)(10,35)(11,-296)])(2,P[10(0,326)(1,296)(2,-461)(3,-175)(4,-123)(5,432)(6,195)(7,447)(8,-135)(9,-87)(10,8)])(3,P[9(0,299)(1,10)(2,-494)(3,192)(4,269)(5,-177)(6,-418)(7,461)(8,352)(9,23)])(4,P[8(0,193)(1,-357)(2,-187)(3,-140)(4,-126)(5,195)(6,-178)(7,-156)(8,-118)])(5,P[7(0,-11)(1,-205)(2,-357)(3,-157)(4,-321)(5,295)(6,-120)(7,-72)])(6,P[6(0,322)(1,443)(2,-26)(3,-29)(4,196)(5,-188)(6,132)])(7,P[5(0,372)(1,-15)(2,-184)(3,-4)(4,453)(5,-134)])(8,P[4(0,26)(1,131)(2,-170)(3,256)(4,282)])(9,P[3(0,198)(1,-424)(2,134)(3,350)])(10,P[2(0,-49)(1,-445)(2,343)])(11,P[1(0,86)(1,417)])(12,P[0(0,-295)])]", // 0
    
    "P[6(0,P[6(0,132371318621897968)(1,62063386199148599600)(2,-151639931559080760922)(3,93583844480566235832)(4,-83154766543076124512)(5,-35432091262507130904)(6,227271934222282896800)])(1,P[5(0,-71827641993430936816)(1,368348281512476672028)(2,-286334917753887437776)(3,204520624662741644980)(4,-153426196026247710960)(5,-806483254283204758968)])(2,P[4(0,-207823824331267897258)(1,407437416975079236352)(2,-109980596108783730996)(3,560473260289176309328)(4,707425975812414448406)])(3,P[3(0,-229793885496425988888)(1,-16258161769731113484)(2,-225314637374518033384)(3,390228626765477166248)])(4,P[2(0,30511269636646463220)(1,-323492766503557473768)(2,-917252466360862051372)])(5,P[1(0,159815596286749757976)(1,567591812420861278704)])(6,P[0(0,-148967343686356735666)])]", // 1
    
    "P[2(0,P[2(0,35204645504740)(1,18690431343280)(2,-33659097374560)])(1,P[1(0,-5878113292740)(1,25167512604160)])(2,P[0(0,41597705375085)])]", // 2

    "P[6(0,P[6(0,-4041376884475638882353)(1,7050387367979191249572)(2,13657295334099123447165)(3,25540501499116672333266)(4,35818873525049083535243)(5,5246125427237612693902)(6,14604340604186011880149)])(1,P[5(0,6901581474556480012157)(1,-7529175981573515299061)(2,-11055581988553126878879)(3,12006542139536462689689)(4,14452481220936855448204)(5,13877670792372320318222)])(2,P[4(0,-11020699216485230516718)(1,-17311543433158205744767)(2,-31286655398777064659338)(3,-31056714506153479404652)(4,4356717231425141453473)])(3,P[3(0,3592444310543685259781)(1,6258214164280712933460)(2,-21871202396326925470570)(3,-7534603152074886821296)])(4,P[2(0,9063109123142130614750)(1,7920859532025064873268)(2,-6296939225750115900496)])(5,P[1(0,5169051546475000208976)(1,893576149924111317664)])(6,P[0(0,973644041162690858383)])]", // 3
    
    "P[3(0,P[3(0,1)(1,3)(2,3)(3,1)])(1,P[2(0,-3)(1,-6)(2,-3)])(2,P[1(0,3)(1,3)])(3,P[0(0,-1)])]", // (x-y+1)^3 // 4
    
    "P[1(0,P[1(1,1)])(1,P[0(0,-1)])]", // (x-y) // 5
    
    "P[2(0,P[2(2,1)])(2,P[0(0,-1)])]", // x^2-y^2 // 6

    "P[1(0,P[2(2,1)])(1,P[0(0,1)])]", // x^2+y // 7 

    "P[2(0,P[2(0,1)(2,2)])(2,P[0(0,-1)])]", // 2x^2-y^2+1 // 8

    "P[7(1,P[7(1,-500)(3,300)(5,-60)(7,4)])(3,P[5(1,300)(3,-147)(5,12)])(5,P[3(1,-60)(3,12)])(7,P[1(1,4)])]", // (4*x)*y^7 + (12*x^3 + (-60)*x)*y^5 + (12*x^5 + (-147)*x^3 + 300*x)*y^3 + (4*x^7 + (-60)*x^5 + 300*x^3 + (-500)*x)*y // 9

    "P[2(0,P[2(0,148)(1,20)(2,1)])(1,P[0(0,16)])(2,P[0(0,1)])]", // y^2 + 16*y + (x^2 + 20*x + 148) //10
    
    "P[0(0,P[1(0,10)(1,1)])]", // x+10 // 11
};



static const int ACK_2_n_polys = 12;

template< class AlgebraicCurveKernel_2  >
void test_algebraic_curve_kernel_2() {

    typedef AlgebraicCurveKernel_2 AK_2;
    
  /*  CGAL_static_assertion( (::boost::is_same< 
            Algebraic_real_1, typename AK::Algebraic_real_1 >::value) );

    CGAL_static_assertion((::boost::is_same<
            Isolator,
            typename AK::Isolator >::value) );
            
    CGAL_static_assertion((::boost::is_same< 
            Coefficient, 
            typename AK::Coefficient >::value));
            
    CGAL_static_assertion((::boost::is_same<
            Polynomial_1,
            typename AK::Polynomial_1 >::value));*/

    typedef typename AK_2::Polynomial_2 Poly_2;
    typedef typename AK_2::Curve_analysis_2 Curve_analysis_2;
    typedef typename Curve_analysis_2::Status_line_1
        Status_line_1;

    typedef typename AK_2::Algebraic_real_1 Algebraic_real_1;
    typedef typename AK_2::Algebraic_real_2 Algebraic_real_2;

    typedef typename AK_2::Coordinate_1 Coordinate_1;
    CGAL_USE_TYPE(Coordinate_1);
    typedef typename AK_2::Coordinate_2 Coordinate_2;
    CGAL_USE_TYPE(Coordinate_2);

    Poly_2 polys[ACK_2_n_polys];

    ::CGAL::set_mode(std::cerr, ::CGAL::IO::PRETTY);
    
    //std::cerr << "constructing curves..\n";
    for(int i = 0; i < ACK_2_n_polys; i++) {
      std::istringstream in(ACK_2_ascii_polys[i]);
        in >> polys[i];    
    }



        

    ///////// testing curve construction //////////

    
    AK_2 kernel_2;
    Curve_analysis_2 c0 = kernel_2.construct_curve_2_object()(polys[0]),
            // make it decomposable
            c1 = kernel_2.construct_curve_2_object()(polys[1]*polys[2]),
            c2 = kernel_2.construct_curve_2_object()(polys[2]),
            c3 = kernel_2.construct_curve_2_object()(polys[3]),
            c5 = kernel_2.construct_curve_2_object()(polys[5]),
            c6 = kernel_2.construct_curve_2_object()(polys[6]);
    //std::cerr << "done..\n";
            
    Status_line_1 line1, line2;
    Algebraic_real_2 xy1, xy2, xy3, xy4;

    ///////////// testing sign_at_2 for non-coprime case /////////////

    {    
        Curve_analysis_2 c7_c6 =
            kernel_2.construct_curve_2_object()(polys[7]*polys[6]);
        assert(c7_c6.number_of_status_lines_with_event() > 0);
        Status_line_1 line = c7_c6.status_line_at_event(0);
        assert(line.number_of_events() > 0);
        Algebraic_real_2 xy = line.algebraic_real_2(0);

        //std::cerr << "done..1.5\n";
        assert(kernel_2.sign_at_2_object()(c7_c6, xy) == CGAL::ZERO);
    }
    {    
        Curve_analysis_2 c7_c6 =
                   kernel_2.construct_curve_2_object()(polys[7]*polys[6]),
           c8_c6 = kernel_2.construct_curve_2_object()(polys[8]*polys[6]);
        
        assert(c7_c6.number_of_status_lines_with_event() > 0);
        Status_line_1 line = c7_c6.status_line_at_event(0);
        assert(line.number_of_events() > 0);
        Algebraic_real_2 xy = line.algebraic_real_2(0);

        //std::cerr << "done..1.6\n";
        assert(kernel_2.sign_at_2_object()(c8_c6, xy) == CGAL::ZERO);
    }

    //std::cerr << "done..2\n";
    

    ///////// test buggy y() //////////////

    { 
        Curve_analysis_2 c9 = kernel_2.construct_curve_2_object()(polys[9]);

        typename Algebraic_real_1::Polynomial_1 f(25,0,-11,0,1);
        
        typedef typename Algebraic_real_1::Rational Rational;
        
        Rational left(-3), right(-2);
        
        Algebraic_real_1 x(f,left,right);

        Algebraic_real_2 xy(x,c9,0);

        xy.y();
    }

    ///////// testing comparison predicates //////////
    
    line1 = c0.status_line_of_interval(0);
    xy1 = line1.algebraic_real_2(1);
    line2 = c1.status_line_at_event(1);
    xy2 = line2.algebraic_real_2(2);
    xy3 = line2.algebraic_real_2(1);

    assert(kernel_2.compare_x_2_object()(xy1, xy2) == CGAL::SMALLER);
    assert(kernel_2.compare_x_2_object()(xy2, xy3) == CGAL::EQUAL);

    assert(kernel_2.compare_xy_2_object()(xy1, xy2) ==
        CGAL::SMALLER);
    assert(kernel_2.compare_xy_2_object()(xy2, xy3) ==
        CGAL::LARGER);

    xy1 = c2.status_line_at_event(0).algebraic_real_2(0);
    line2 = c3.status_line_at_event(0);
    xy2 = line2.algebraic_real_2(2);
    xy3 = line2.algebraic_real_2(3);
    xy4 = line2.algebraic_real_2(4);

    //std::cerr << "y_comp 1" << std::flush;
    assert(kernel_2.compare_y_2_object()(xy1, xy2) == CGAL::LARGER);
    //std::cerr << " 2" << std::flush;
    assert(kernel_2.compare_y_2_object()(xy1, xy3) == CGAL::SMALLER);
    //std::cerr << " 3" << std::flush;
    assert(kernel_2.compare_y_2_object()(xy1, xy4) == CGAL::SMALLER);
    //std::cerr << " 4" << std::flush;
    assert(kernel_2.compare_y_2_object()(xy2, xy3) == CGAL::SMALLER);
    //std::cerr << " done" << std::endl;

    /////// testing squarefreeness and coprimality /////////

    /*assert(
        kernel_2.has_finite_number_of_self_intersections_2_object()
            (polys[0]));
    assert(
        !kernel_2.has_finite_number_of_self_intersections_2_object()
            (polys[4])); // non-squarefree*/
    
    assert(
        kernel_2.has_finite_number_of_intersections_2_object()
        (c2.polynomial_2(), c3.polynomial_2())); // coprime
            
    assert(
        !kernel_2.has_finite_number_of_intersections_2_object()
            (c1.polynomial_2(), c2.polynomial_2())); // non-coprime
            
    assert(
        !kernel_2.has_finite_number_of_intersections_2_object()
            (c5.polynomial_2(), c6.polynomial_2())); // non-coprime

    //////// testing decompose ///////////
            
   assert((kernel_2.decompose_2_object()(polys[4])) !=
        polys[4]); // non-squarefree

   assert((kernel_2.decompose_2_object()(polys[3])) ==
        polys[3]);

    typedef std::vector<Curve_analysis_2> Curves_2;
    typedef std::vector<int> Int_vector;
    Curves_2 parts;
    Int_vector mults;
    typename Curves_2::const_iterator cit;
    typename Int_vector::const_iterator iit;

    //std::cerr << "n_factors: " << kernel_2.decompose_2_object()(c1,
    //    std::back_inserter(parts),  std::back_inserter(mults));
    
    for(cit = parts.begin(), iit = mults.begin(); cit != parts.end();
            cit++, iit++) {
        //std::cerr << "part: " << cit->f() << "; mult: " << *iit << "\n";
    }

    parts.clear();
    mults.clear();


    //std::cerr << "n_factors :" <<
    //    kernel_2.decompose_2_object()(c3, std::back_inserter(parts),
    //        std::back_inserter(mults));
    for(cit = parts.begin(), iit = mults.begin(); cit != parts.end();
            cit++, iit++) {
        //std::cerr << "part: " << cit->f() << "; mult: " << *iit << "\n";
    }

    Curves_2 fgs, fs, gs;

    assert(kernel_2.decompose_2_object()(c5, c6,
        std::back_inserter(fs), std::back_inserter(gs),
            std::back_inserter(fgs))); // must have a common part

    assert(fgs.size() == 1);
    assert(fs.size() == 0);
    assert(gs.size() == 1);
    
    /*std::cerr << "common: " << fgs.size() << "; " <<
        fs.size() << "; " << gs.size() << "\n";*/

    fgs.clear();
    fs.clear();
    gs.clear();
    assert(!kernel_2.decompose_2_object()(c2, c3,
        std::back_inserter(fs), std::back_inserter(gs),
            std::back_inserter(fgs))); // must be coprime

    assert(fgs.size() == 0);
    assert(fs.size() == 1);
    assert(gs.size() == 1);
            
    /*std::cerr << "coprime: " << fgs.size() << "; " <<
        fs.size() << "; " << gs.size() << "\n";*/
    
    //////// testing x/y-critical points ////////
        
    typedef std::vector<Algebraic_real_2> Xy_coords;

    Xy_coords points;
    typename Xy_coords::const_iterator xyit;
    kernel_2.x_critical_points_2_object()(c0, std::back_inserter(points));

    assert(points.size() == 10);
    //std::cerr << "testing x-critical points: " << points.size() << "\n";

    points.clear();
    kernel_2.y_critical_points_2_object()(c0, std::back_inserter(points));

    assert(points.size() == 10);

    ///////// testing sign_2 //////////////

    Status_line_1 line;
    Algebraic_real_2 xy;
        
    int n_lines = c0.number_of_status_lines_with_event(), ii, jj,
        n_events;
    for(ii = 0; ii < n_lines; ii++) {
        
        line = c0.status_line_at_event(ii);
        n_events = line.number_of_events();
        //std::cout << ii << "pts at event: \n";
        for(jj = 0; jj < n_events; jj++) {
            xy = line.algebraic_real_2(jj);
            //std::cout << "sign 2: " <<
            //    kernel_2.sign_at_2_object()(c1, xy) << "\n\n";
        }
            
        //std::cout << ii << "pts over interval: \n";
        line = c0.status_line_of_interval(ii);
        n_events = line.number_of_events();
        for(jj = 0; jj < n_events; jj++) {
            xy = line.algebraic_real_2(jj);
            //std::cout << " sign 2: " <<
            //    kernel_2.sign_at_2_object()(c1, xy) << "\n\n";
        }
    }
    
    
    ///////////// testing solve_2 /////////////

    points.clear();
    mults.clear();
    kernel_2.solve_2_object()(c2, c3, std::back_inserter(points),
        std::back_inserter(mults));

    for(xyit = points.begin(), iit = mults.begin(); xyit != points.end();
            xyit++, iit++) {
        //std::cerr << "pt: " << *xyit << "; mult: " << *iit << "\n";
    }

    points.clear();
    mults.clear();
    Curve_analysis_2 c10 = kernel_2.construct_curve_2_object()(polys[10]),
        c11 = kernel_2.construct_curve_2_object()(polys[11]);
    kernel_2.solve_2_object()(c10, c11, std::back_inserter(points),
        std::back_inserter(mults));
    assert(points.size()==2);
    points.clear();
    mults.clear();
    kernel_2.solve_2_object()(c11, c10, std::back_inserter(points),
        std::back_inserter(mults));
    assert(points.size()==2);



    ///////////// testing Swap_x_and_y /////////////
    {
        Curve_analysis_2 c7_c6 =
            kernel_2.construct_curve_2_object()(polys[7]*polys[6]);
        
        Curve_analysis_2 swapped = kernel_2.swap_x_and_y_2_object()(c7_c6);
        (void) swapped;

    }
    ///////////// Overloaded sign_at_1 functor
    {
      typedef typename AK_2::Polynomial_1 Polynomial_1;
      Polynomial_1 p(5,-4,3,-2,-1,1);
      AK_2 kernel_2;
      Algebraic_real_1 ar(p,-3,3);
      assert(kernel_2.sign_at_1_object()(p,ar,1000)==CGAL::ZERO);
      Polynomial_1 q(1,0,0,-1,1,0,0,1);
      Algebraic_real_1 ar2(q,-2,0);
      assert(kernel_2.sign_at_1_object()(p,ar2)==CGAL::POSITIVE);
      Polynomial_1 r(500001,-400000,300000,-200000,-100000,100000);
      Algebraic_real_1 ar3(r,-5,5);
      assert(kernel_2.sign_at_1_object()(p,ar3)==CGAL::NEGATIVE);
    }
    
}

} //namespace internal

} //namespace CGAL

#endif // CGAL_TEST_ALGEBRAIC_CURVE_KERNEL_2_H
