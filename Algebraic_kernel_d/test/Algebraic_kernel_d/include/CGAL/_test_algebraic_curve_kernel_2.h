// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     :   Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#include <CGAL/basic.h>
#include <CGAL/Testsuite/assert.h>

//#include <CGAL/Algebraic_kernel_1.h>

//#include <CGAL/_test_basic.h>

#ifndef CGAL_TEST_ALGEBRAIC_CURVE_KERNEL_2_H
#define CGAL_TEST_ALGEBRAIC_CURVE_KERNEL_2_H

CGAL_BEGIN_NAMESPACE

namespace CGALi {

static const char *ACK_2_ascii_polys[] = {

    "P[12(0,P[12(0,-403)(1,445)(2,-244)(3,-261)(4,182)(5,-284)(6,321)(7,-150)(8,-312)(9,-395)(10,451)(11,166)(12,-268)])(1,P[11(0,-149)(1,-194)(2,190)(3,397)(4,-423)(5,507)(6,-321)(7,215)(8,-459)(9,-112)(10,35)(11,-296)])(2,P[10(0,326)(1,296)(2,-461)(3,-175)(4,-123)(5,432)(6,195)(7,447)(8,-135)(9,-87)(10,8)])(3,P[9(0,299)(1,10)(2,-494)(3,192)(4,269)(5,-177)(6,-418)(7,461)(8,352)(9,23)])(4,P[8(0,193)(1,-357)(2,-187)(3,-140)(4,-126)(5,195)(6,-178)(7,-156)(8,-118)])(5,P[7(0,-11)(1,-205)(2,-357)(3,-157)(4,-321)(5,295)(6,-120)(7,-72)])(6,P[6(0,322)(1,443)(2,-26)(3,-29)(4,196)(5,-188)(6,132)])(7,P[5(0,372)(1,-15)(2,-184)(3,-4)(4,453)(5,-134)])(8,P[4(0,26)(1,131)(2,-170)(3,256)(4,282)])(9,P[3(0,198)(1,-424)(2,134)(3,350)])(10,P[2(0,-49)(1,-445)(2,343)])(11,P[1(0,86)(1,417)])(12,P[0(0,-295)])]", // 0
    
    "P[6(0,P[6(0,132371318621897968)(1,62063386199148599600)(2,-151639931559080760922)(3,93583844480566235832)(4,-83154766543076124512)(5,-35432091262507130904)(6,227271934222282896800)])(1,P[5(0,-71827641993430936816)(1,368348281512476672028)(2,-286334917753887437776)(3,204520624662741644980)(4,-153426196026247710960)(5,-806483254283204758968)])(2,P[4(0,-207823824331267897258)(1,407437416975079236352)(2,-109980596108783730996)(3,560473260289176309328)(4,707425975812414448406)])(3,P[3(0,-229793885496425988888)(1,-16258161769731113484)(2,-225314637374518033384)(3,390228626765477166248)])(4,P[2(0,30511269636646463220)(1,-323492766503557473768)(2,-917252466360862051372)])(5,P[1(0,159815596286749757976)(1,567591812420861278704)])(6,P[0(0,-148967343686356735666)])]", // 1
    
    "P[2(0,P[2(0,35204645504740)(1,18690431343280)(2,-33659097374560)])(1,P[1(0,-5878113292740)(1,25167512604160)])(2,P[0(0,41597705375085)])]", // 2

    "P[6(0,P[6(0,-4041376884475638882353)(1,7050387367979191249572)(2,13657295334099123447165)(3,25540501499116672333266)(4,35818873525049083535243)(5,5246125427237612693902)(6,14604340604186011880149)])(1,P[5(0,6901581474556480012157)(1,-7529175981573515299061)(2,-11055581988553126878879)(3,12006542139536462689689)(4,14452481220936855448204)(5,13877670792372320318222)])(2,P[4(0,-11020699216485230516718)(1,-17311543433158205744767)(2,-31286655398777064659338)(3,-31056714506153479404652)(4,4356717231425141453473)])(3,P[3(0,3592444310543685259781)(1,6258214164280712933460)(2,-21871202396326925470570)(3,-7534603152074886821296)])(4,P[2(0,9063109123142130614750)(1,7920859532025064873268)(2,-6296939225750115900496)])(5,P[1(0,5169051546475000208976)(1,893576149924111317664)])(6,P[0(0,973644041162690858383)])]", // 3
    
    "P[3(0,P[3(0,1)(1,3)(2,3)(3,1)])(1,P[2(0,-3)(1,-6)(2,-3)])(2,P[1(0,3)(1,3)])(3,P[0(0,-1)])]", // (x-y+1)^3 // 4
    
    "P[1(0,P[1(1,1)])(1,P[0(0,-1)])]", // (x-y) // 5
    
    "P[2(0,P[2(2,1)])(2,P[0(0,-1)])]", // x^2-y^2 // 6

    "P[1(0,P[2(2,1)])(1,P[0(0,1)])]", // x^2+y // 7 

    "P[2(0,P[2(0,1)(2,2)])(2,P[0(0,-1)])]" // 2x^2-y^2+1 // 8

};



static const int ACK_2_n_polys = 9;

template< class AlgebraicCurveKernel_2  >
void test_algebraic_curve_kernel_2() {

    typedef AlgebraicCurveKernel_2 AK_2;
    
  /*  BOOST_STATIC_ASSERT( (::boost::is_same< 
            Algebraic_real_1, typename AK::Algebraic_real_1 >::value) );

    BOOST_STATIC_ASSERT((::boost::is_same<
            Isolator,
            typename AK::Isolator >::value) );
            
    BOOST_STATIC_ASSERT((::boost::is_same< 
            Coefficient, 
            typename AK::Coefficient >::value));
            
    BOOST_STATIC_ASSERT((::boost::is_same<
            Polynomial_1,
            typename AK::Polynomial_1 >::value));*/

    typedef typename AK_2::Curve_2 Curve_2;
    typedef typename Curve_2::Poly_d Internal_poly_2;
    typedef typename AK_2::Curve_analysis_2 Curve_analysis_2;
    typedef typename Curve_analysis_2::Status_line_1
        Status_line_1;

    typedef typename AK_2::X_coordinate_1 X_coordinate_1;
    typedef typename AK_2::Xy_coordinate_2 Xy_coordinate_2;

    Internal_poly_2 polys[ACK_2_n_polys];

    ::CGAL::set_mode(std::cerr, ::CGAL::IO::PRETTY);
    
    std::cerr << "constructing curves..\n";
    for(int i = 0; i < ACK_2_n_polys; i++) {
        istringstream in(ACK_2_ascii_polys[i]);
        in >> polys[i];    
    }

    ///////// testing curve construction //////////
    
    AK_2 kernel_2;
    Curve_2 c0 = kernel_2.construct_curve_2_object()(polys[0]),
            // make it decomposable
            c1 = kernel_2.construct_curve_2_object()(polys[1]*polys[2]),
            c2 = kernel_2.construct_curve_2_object()(polys[2]),
            c3 = kernel_2.construct_curve_2_object()(polys[3]),
            c5 = kernel_2.construct_curve_2_object()(polys[5]),
            c6 = kernel_2.construct_curve_2_object()(polys[6]);
    std::cerr << "done..\n";
            
    Curve_analysis_2 ca0(c0), ca1(c1);
    Status_line_1 line1, line2;
    Xy_coordinate_2 xy1, xy2, xy3, xy4;

    ///////////// testing sign_at_2 for non-coprime case /////////////
    {    
        Curve_2 c7_c6 = kernel_2.construct_curve_2_object()(polys[7]*polys[6]);
        Curve_analysis_2 ca_c6(c6);
        CGAL_test_assert(ca_c6.number_of_status_lines_with_event() > 0);
        Status_line_1 line = ca_c6.status_line_at_event(0);
        CGAL_test_assert(line.number_of_events() > 0);
        Xy_coordinate_2 xy = line.algebraic_real_2(0);
        CGAL_test_assert(kernel_2.sign_at_2_object()(c7_c6,xy) == CGAL::ZERO);
    }
    {    
        Curve_2 c7_c6 = kernel_2.construct_curve_2_object()(polys[7]*polys[6]);
        Curve_2 c8_c6 = kernel_2.construct_curve_2_object()(polys[8]*polys[6]);
        Curve_analysis_2 ca_c7c6(c7_c6);
        CGAL_test_assert(ca_c7c6.number_of_status_lines_with_event() > 0);
        Status_line_1 line = ca_c7c6.status_line_at_event(0);
        CGAL_test_assert(line.number_of_events() > 0);
        Xy_coordinate_2 xy = line.algebraic_real_2(0);
        CGAL_test_assert(kernel_2.sign_at_2_object()(c8_c6,xy) == CGAL::ZERO);
    }
    


    ///////// testing comparison predicates //////////
    
    line1 = ca0.status_line_of_interval(0);
    xy1 = line1.algebraic_real_2(1);
    line2 = ca1.status_line_at_event(1);
    xy2 = line2.algebraic_real_2(2);
    xy3 = line2.algebraic_real_2(1);

    CGAL_test_assert(kernel_2.compare_x_2_object()(xy1, xy2) == CGAL::SMALLER);
    CGAL_test_assert(kernel_2.compare_x_2_object()(xy2, xy3) == CGAL::EQUAL);

    CGAL_test_assert(kernel_2.compare_xy_2_object()(xy1, xy2) ==
        CGAL::SMALLER);
    CGAL_test_assert(kernel_2.compare_xy_2_object()(xy2, xy3) ==
        CGAL::LARGER);

    Curve_analysis_2 ca2(c2), ca3(c3);
    xy1 = ca2.status_line_at_event(0).algebraic_real_2(0);
    line2 = ca3.status_line_at_event(0);
    xy2 = line2.algebraic_real_2(2);
    xy3 = line2.algebraic_real_2(3);
    xy4 = line2.algebraic_real_2(4);

    std::cerr << "y_comp 1" << std::flush;
    CGAL_test_assert(kernel_2.compare_y_2_object()(xy1, xy2) == CGAL::LARGER);
    std::cerr << " 2" << std::flush;
    CGAL_test_assert(kernel_2.compare_y_2_object()(xy1, xy3) == CGAL::SMALLER);
    std::cerr << " 3" << std::flush;
    CGAL_test_assert(kernel_2.compare_y_2_object()(xy1, xy4) == CGAL::SMALLER);
    std::cerr << " 4" << std::flush;
    CGAL_test_assert(kernel_2.compare_y_2_object()(xy2, xy3) == CGAL::SMALLER);
    std::cerr << " done" << std::endl;

    /////// testing squarefreeness and coprimality /////////
     
    /*CGAL_test_assert(
        kernel_2.has_finite_number_of_self_intersections_2_object()
            (polys[0]));
    CGAL_test_assert(
        !kernel_2.has_finite_number_of_self_intersections_2_object()
            (polys[4])); // non-squarefree*/

    CGAL_test_assert(
        kernel_2.has_finite_number_of_intersections_2_object()
            (c2, c3)); // coprime
            
    CGAL_test_assert(
        !kernel_2.has_finite_number_of_intersections_2_object()
            (c1, c2)); // non-coprime
            
    CGAL_test_assert(
        !kernel_2.has_finite_number_of_intersections_2_object()
            (c5, c6)); // non-coprime

    //////// testing decompose ///////////
            
   CGAL_test_assert((kernel_2.decompose_2_object()(polys[4])) !=
        polys[4]); // non-squarefree

   CGAL_test_assert((kernel_2.decompose_2_object()(polys[3])) ==
        polys[3]);

    typedef std::vector<Curve_2> Curves_2;
    typedef std::vector<int> Int_vector;
    Curves_2 parts;
    Int_vector mults;
    typename Curves_2::const_iterator cit;
    typename Int_vector::const_iterator iit;

    std::cerr << "n_factors: " << kernel_2.decompose_2_object()(c1,
        std::back_inserter(parts),  std::back_inserter(mults));
    
    for(cit = parts.begin(), iit = mults.begin(); cit != parts.end();
            cit++, iit++) {
        //std::cerr << "part: " << cit->f() << "; mult: " << *iit << "\n";
    }

    parts.clear();
    mults.clear();

    std::cerr << "n_factors :" <<
        kernel_2.decompose_2_object()(c3, std::back_inserter(parts),
            std::back_inserter(mults));
    for(cit = parts.begin(), iit = mults.begin(); cit != parts.end();
            cit++, iit++) {
        //std::cerr << "part: " << cit->f() << "; mult: " << *iit << "\n";
    }

    Curves_2 fgs, fs, gs;

    CGAL_test_assert(kernel_2.decompose_2_object()(c5, c6,
        std::back_inserter(fs), std::back_inserter(gs),
            std::back_inserter(fgs))); // must have a common part

    CGAL_test_assert(fgs.size() == 1);
    CGAL_test_assert(fs.size() == 0);
    CGAL_test_assert(gs.size() == 1);
    
    /*std::cerr << "common: " << fgs.size() << "; " <<
        fs.size() << "; " << gs.size() << "\n";*/

    fgs.clear();
    fs.clear();
    gs.clear();
    CGAL_test_assert(!kernel_2.decompose_2_object()(c2, c3,
        std::back_inserter(fs), std::back_inserter(gs),
            std::back_inserter(fgs))); // must be coprime

    CGAL_test_assert(fgs.size() == 0);
    CGAL_test_assert(fs.size() == 1);
    CGAL_test_assert(gs.size() == 1);
            
    /*std::cerr << "coprime: " << fgs.size() << "; " <<
        fs.size() << "; " << gs.size() << "\n";*/
    
    ///////// testing derivation //////////

    typename AK_2::NiX2CGAL_converter cvt;
    typename AK_2::CGAL2NiX_converter cvt_back;

    cvt_back(kernel_2.derivative_x_2_object()(cvt(polys[0])));
        
    cvt_back(kernel_2.derivative_y_2_object()(cvt(polys[0])));

    //////// testing x/y-critical points ////////
        
    typedef std::vector<Xy_coordinate_2> Xy_coords;

    Xy_coords points;
    typename Xy_coords::const_iterator xyit;
    kernel_2.x_critical_points_2_object()(c0, std::back_inserter(points));

    CGAL_test_assert(points.size() == 10);
    //std::cerr << "testing x-critical points: " << points.size() << "\n";

    points.clear();
    kernel_2.y_critical_points_2_object()(c0, std::back_inserter(points));

    CGAL_test_assert(points.size() == 10);

    ///////// testing sign_2 //////////////

    Status_line_1 line;
    Xy_coordinate_2 xy;
        
    int n_lines = ca0.number_of_status_lines_with_event(), ii, jj,
        n_events;
    for(ii = 0; ii < n_lines; ii++) {
        
        line = ca0.status_line_at_event(ii);
        n_events = line.number_of_events();
        std::cout << ii << "pts at event: \n";
        for(jj = 0; jj < n_events; jj++) {
            xy = line.algebraic_real_2(jj);
            std::cout << "sign 2: " <<
                kernel_2.sign_at_2_object()(c1, xy) << "\n\n";
        }
            
        std::cout << ii << "pts over interval: \n";
        line = ca0.status_line_of_interval(ii);
        n_events = line.number_of_events();
        for(jj = 0; jj < n_events; jj++) {
            xy = line.algebraic_real_2(jj);
            std::cout << " sign 2: " <<
                kernel_2.sign_at_2_object()(c1, xy) << "\n\n";
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

    

}

} //namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_TEST_ALGEBRAIC_CURVE_KERNEL_2_H
