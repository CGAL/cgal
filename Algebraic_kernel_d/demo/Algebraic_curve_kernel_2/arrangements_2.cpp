// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================


#include <CGAL/Algebraic_curve_kernel_2/flags.h>

// demo-speficic flags

// Allows to use the Filtered_curve_kernel_via_analysis_2
#ifndef CGAL_ACK_USE_FILTERED_CKvA_2
#define CGAL_ACK_USE_FILTERED_CKvA_2 0
#endif

// If set, only one (random) segment of each input curve is choosen
#ifndef CGAL_ACK_ONE_SEGMENT_PER_CURVE
#define CGAL_ACK_ONE_SEGMENT_PER_CURVE 0
#endif

// Defines the underlying ArithmeticKernel
#define ARITHMETIC_KERNEL CGAL::CORE_arithmetic_kernel



#include <CGAL/basic.h>

#include <CGAL/Timer.h>
CGAL::Timer overall_timer;

#include <sstream>
#include <fstream>

#include <CGAL/Timer.h>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_1.h>
#include <CGAL/Algebraic_curve_kernel_2.h>


#if CGAL_ACK_USE_FILTERED_CKvA_2
#include <CGAL/Filtered_algebraic_curve_kernel_2.h>
#include <CGAL/Filtered_curved_kernel_via_analysis_2.h>
#else
#include <CGAL/Curved_kernel_via_analysis_2.h>
#endif

#include <CGAL/Arrangement_2.h>

template<typename Poly1> Poly1 
random_dense_univariate_polynomial(int degree,int bitsize) {
    typedef typename Poly1::NT NT;
    std::vector<NT> coeffs;
    for(int i=0;i<=degree;i++) {
        NT coeff=0;
        for(int j=0;j<bitsize-1;j++) {
            coeff = 2*coeff + (lrand48()%2);
        }
        // The last bit determines the sign
        if(lrand48()%2==0) {
            coeff=-coeff;
        }    
        coeffs.push_back(coeff);
    }
    return Poly1(coeffs.begin(),coeffs.end());
} 

template<typename Poly2> Poly2
random_dense_bivariate_polynomial(int degree,int bitsize) {
    typedef typename Poly2::NT Poly1;
    std::vector<Poly1> coeffs;
    for(int i=0;i<=degree;i++) {
        coeffs.push_back(random_dense_univariate_polynomial<Poly1>(degree-i,bitsize));
    }
    return Poly2(coeffs.begin(),coeffs.end());
}

int main(int argc, char** argv) {
    if(argc<3) {
        std::cerr << "Needs method specification and Input file" << std::endl;
        return 1;
    }
    typedef ARITHMETIC_KERNEL AK;

    typedef AK::Integer Coefficient;
    typedef AK::Rational Rational;

    typedef CGAL::CGALi::Algebraic_real_quadratic_refinement_rep_bfi
        < Coefficient, Rational > Rep_class;
    typedef CGAL::CGALi::Bitstream_descartes< CGAL::Polynomial< Coefficient >, 
        Rational > Isolator;
    
    typedef CGAL::Algebraic_kernel_1<Coefficient,Rational,Rep_class, Isolator> 
        Algebraic_kernel_1;

#if !CGAL_ACK_USE_FILTERED_CKvA_2
    typedef CGAL::Algebraic_curve_kernel_2<Algebraic_kernel_1>
        Algebraic_curve_kernel_2;
#else
    typedef CGAL::Filtered_algebraic_curve_kernel_2<Algebraic_kernel_1>
        Algebraic_curve_kernel_2;
#endif
    typedef Algebraic_curve_kernel_2::Polynomial_2 Polynomial_2;

    std::vector<Polynomial_2> curves;

    int arrangement_type;
  
    std::string str(argv[1]);

    if(str=="LEDA" || str=="Leda" || str=="leda") {
        arrangement_type=1;
    } else if(str=="CGAL" || str=="Cgal" || str=="cgal") {
        arrangement_type=2;
    } else if(str=="NAIV" || str=="Naiv" || str=="naiv") {
        arrangement_type=3;
    } else {
        std::cerr << "Second argument must specify the arrangement type: " 
                  << "LEDA or CGAL possible" << std::endl;
        return 1;
    }

    if(argc<3) {
        std::cerr << "Input file needed! (or random for random generated arrangements)" 
                  << std::endl;
        std::exit(1);
    }
    std::string file(argv[2]);

    if(file=="RANDOM" || file=="random" || file=="Random") {
        if(argc<6) {
            std::cout << "Need to specify: Number of curves, degree, bitsize" 
                      << std::endl;
        }
        int no_curves = atoi(argv[3]);
        int max_degree = atoi(argv[4]);
        int max_coeff = atoi(argv[5]);
        for(int i=0;i<no_curves;i++) {
            Polynomial_2 curr_curve = random_dense_bivariate_polynomial<Polynomial_2>(max_degree,max_coeff);
            curves.push_back(curr_curve);
            std::cout << curr_curve << std::endl;
        }
    } else {
      
        std::ifstream input(file.c_str());
    
        while(!input.eof()) {
            Polynomial_2 f;
            int g = input.get();
            if(g=='P') {
                input.putback(g);
                input >> f; 
                curves.push_back(CGAL::CGALi::canonicalize_polynomial(f));
            }
        }
    }
  
    ::CGAL::set_pretty_mode(std::cout);
    std::cout << "Ready for the arrangement computation of " 
              << curves.size() << " input curves:" << std::endl;
    for(std::vector<Polynomial_2>::iterator it=curves.begin();
        it!=curves.end();it++) {
        std::cout << *it << std::endl;
    }
    CGAL::Timer overall_timer;
    overall_timer.start();
    
#if !CGAL_ACK_USE_FILTERED_CKvA_2
    typedef CGAL::Curved_kernel_via_analysis_2< Algebraic_curve_kernel_2 > 
        Curved_kernel_2; 
#else
    typedef CGAL::Curved_kernel_via_analysis_2< Algebraic_curve_kernel_2 > 
        Exact_curved_kernel_2; 
    typedef CGAL::Filtered_curved_kernel_via_analysis_2<Exact_curved_kernel_2>
        Curved_kernel_2; 
#endif

    Curved_kernel_2 curve_kernel;

    Curved_kernel_2::Make_x_monotone_2 make_x_monotone = 
        Curved_kernel_2::instance().make_x_monotone_2_object();

    std::vector<CGAL::Object> sweepable_objects;  

    for(std::vector<Polynomial_2>::iterator it=curves.begin();
        it!=curves.end();it++) {
        Curved_kernel_2::Curve_kernel_2::Construct_curve_2 construct_curve =
            Curved_kernel_2::instance().kernel().construct_curve_2_object();
        
        Curved_kernel_2::Curve_kernel_2::Curve_analysis_2 curr_curve = 
            construct_curve(*it);
#if CGAL_ACK_ONE_SEGMENT_PER_CURVE
#warning Warning, only one segment per curve is chosen
        std::vector<CGAL::Object> curr_sweepable_objects;  
        make_x_monotone(curr_curve,std::back_inserter(curr_sweepable_objects));
        int no_segments = static_cast<int>(curr_sweepable_objects.size());
        if(no_segments > 0) {
            int curr_index = lrand48()%no_segments;
            std::cout << curr_index << ", " << std::flush;
            sweepable_objects.push_back(curr_sweepable_objects[curr_index]);
        }
        
#else
        make_x_monotone(curr_curve,std::back_inserter(sweepable_objects));
#endif

    }
    std::cout << sweepable_objects.size() 
              << " sweepable segments found" << std::endl;
    
    std::cout << "Time for the one-curve analysis so far: " 
              << overall_timer.time() << std::endl;
    
    if(arrangement_type==1) {

        // TODO EXACUS-sweep!
        
        CGAL_error_msg("The sweep-line method based on LEDA is not supported at the moment, try CGAL instead");

    } else if(arrangement_type==2) {
        typedef CGAL::Arrangement_2<Curved_kernel_2> CGAL_Arrangement_2;
        CGAL_Arrangement_2 cgal_arrangement;
        
/*
        std::vector<CGAL::Object> new_sweepable_objects;
        if(argc==3) {
            new_sweepable_objects=sweepable_objects;
        } else {
            for(int i=3;i<argc;i++) {
                new_sweepable_objects.push_back(sweepable_objects[atoi(argv[i])]);
            }
        }
*/
        std::cout << "Start sweep with " << sweepable_objects.size() << " segments" << std::endl;
        std::vector<Curved_kernel_2::X_monotone_curve_2> segments;
        std::vector<Curved_kernel_2::Point_2> isol_points;

        for( std::vector<CGAL::Object>::iterator it 
                 = sweepable_objects.begin();
             it != sweepable_objects.end();
             it++ ) {

            Curved_kernel_2::X_monotone_curve_2 curr_segment;
            Curved_kernel_2::Point_2 curr_point;

            if(CGAL::assign(curr_segment,*it)) {
                segments.push_back(curr_segment);
            } else {
                CGAL_assertion_code(bool check = )
                    CGAL::assign(curr_point,*it);
                CGAL_assertion(check);
                isol_points.push_back(curr_point);
            }

        }

/*
        CGAL::make_x_monotone(curves.begin(), curves.end(),
                              std::back_inserter(segments),
                              std::back_inserter(isol_points),
                              &curve_kernel);
*/
      

        CGAL::insert_empty(cgal_arrangement,
                           segments.begin(),
                           segments.end(),
                           isol_points.begin(),
                           isol_points.end());

        overall_timer.stop();
        std::cout << "****************** RESULTS ***************** " << std::endl;
        std::cout << cgal_arrangement.number_of_vertices() << " nodes" << std::endl;
        std::cout << cgal_arrangement.number_of_edges() << " edges" << std::endl;
        std::cout << cgal_arrangement.number_of_faces() << " faces" << std::endl;
        std::cout << overall_timer.time() << " time elpased in total" << std::endl;
    }
    
    return 0;
}
