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

#ifndef AcX_DEBUG_PRINT
#define AcX_DEBUG_PRINT 0
#endif

#if AcX_DEBUG_PRINT
#define AcX_DSTREAM(str) std::cout << str;
#else
#define AcX_DSTREAM(str) 
#endif


#ifndef AcX_STATIC_SEED
#define AcX_STATIC_SEED 0
#endif

#ifndef BITSTREAM_USES_E08_TREE
#define BITSTREAM_USES_E08_TREE 0
#endif

#ifndef AcX_USE_CGAL_FILTERED_CKvA_2
#define AcX_USE_CGAL_FILTERED_CKvA_2 0
#endif

#ifndef AcX_ONE_SEGMENT_PER_CURVE
#define AcX_ONE_SEGMENT_PER_CURVE 0
#endif

#ifndef CGAL_ACK_USE_CORE
#define CGAL_ACK_USE_CORE 1
#endif
#ifndef CGAL_ACK_USE_LEDA
#define CGAL_ACK_USE_LEDA 0
#endif

#ifndef AcX_USE_NO_BFI_APPROX_IN_BITSTREAM_TRAITS
#define AcX_USE_NO_BFI_APPROX_IN_BITSTREAM_TRAITS 1
#endif

#include <CGAL/basic.h>

#ifndef NiX_USE_QUADRATIC_REFINEMENT 
#ifndef NiX_USE_QUADRATIC_REFINEMENT_BFI
#define NiX_REFINEMENTS_BEFORE_GCD 16
#endif
#endif

#ifndef AcX_CHECK_POLYNOMIALS_FOR_COPRIMABILITY
#define AcX_CHECK_POLYNOMIALS_FOR_COPRIMABILITY 0
#endif

#ifndef AcX_SECOND_RUN_OF_SWEEP
#define AcX_SECOND_RUN_OF_SWEEP 0
#endif

#ifndef AcX_USE_MAPLE_FOR_MODULAR_RESULTANT
#define AcX_USE_MAPLE_FOR_MODULAR_RESULTANT 0
#endif
#ifndef AcX_SPEED_UP_FOR_REGULAR_CURVES
#define AcX_SPEED_UP_FOR_REGULAR_CURVES 0
#endif

#ifndef AcX_SPEED_UP_FOR_DEGREE_GREATER_EQUAL
#define AcX_SPEED_UP_FOR_DEGREE_GREATER_EQUAL 999
#endif

#ifndef AcX_NO_ARC_FLIP
#define AcX_NO_ARC_FLIP 0
#endif

#include <CGAL/Timer.h>
CGAL::Timer overall_timer;

#if AcX_USE_GLOBAL_TIMERS_ARR
#include<AcX/macros.h>
#ifndef AcX_USE_GLOBAL_TIMERS
#define AcX_USE_GLOBAL_TIMERS 1
#endif
AcX_TIMERS
#endif

#ifndef AcX_COEFFICIENT
#define AcX_COEFFICIENT Integer
#endif

#include <sstream>
#include <fstream>

#include <CGAL/Timer.h>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_1.h>
#include <CGAL/Algebraic_curve_kernel_2.h>


#if AcX_USE_CGAL_FILTERED_CKvA_2
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
#if CGAL_ACK_USE_CORE
    typedef CGAL::CORE_arithmetic_kernel AK;
#elif defined CGAL_ACK_USE_LEDA
    typedef CGAL::LEDA_arithmetic_kernel AK;
#else
#error
#endif

    typedef AK::Integer Coefficient;
    typedef AK::Rational Rational;

    typedef CGAL::CGALi::Algebraic_real_quadratic_refinement_rep_bfi
        < Coefficient, Rational > Rep_class;
    typedef CGAL::CGALi::Bitstream_descartes< CGAL::Polynomial< Coefficient >, 
        Rational > Isolator;
    
    typedef CGAL::Algebraic_kernel_1<Coefficient,Rational,Rep_class, Isolator> 
        Algebraic_kernel_1;

#if !AcX_USE_CGAL_FILTERED_CKvA_2
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
    
#if !AcX_USE_CGAL_FILTERED_CKvA_2
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

#if AcX_ONE_SEGMENT_PER_CURVE
#if AcX_STATIC_SEED
    srand48(AcX_STATIC_SEED);
#else
    srand48(time(NULL));
#endif
#endif
    for(std::vector<Polynomial_2>::iterator it=curves.begin();
        it!=curves.end();it++) {
        Curved_kernel_2::Curve_kernel_2::Construct_curve_2 construct_curve =
            Curved_kernel_2::instance().kernel().construct_curve_2_object();
        
        Curved_kernel_2::Curve_kernel_2::Curve_analysis_2 curr_curve = 
            construct_curve(*it);
#if AcX_ONE_SEGMENT_PER_CURVE
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
