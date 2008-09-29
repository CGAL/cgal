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

#ifndef CGAL_ACK_DEBUG_FLAG
#define CGAL_ACK_DEBUG_FLAG 1
#endif

#ifndef CGAL_ACK_DEBUG_PRINT
#define CGAL_ACK_DEBUG_PRINT std::cout
#endif

#include<CGAL/Algebraic_curve_kernel_2/flags.h>

// demo-specific flags

// Is the input curve made square free beforehand, 
// or does the program rely on a square free input?
#ifndef CGAL_ACK_MAKE_SQUARE_FREE
#define CGAL_ACK_MAKE_SQUARE_FREE 0
#endif

#ifndef CGAL_ACK_USE_APPROXIMATE_ROTATION
#define CGAL_ACK_USE_APPROXIMATE_ROTATION 0
#endif

#if !CGAL_ACK_USE_APPROXIMATE_ROTATION
#ifndef CGAL_ACK_BASE_ANGLE
#define CGAL_ACK_BASE_ANGLE 30
#endif
#endif

#if CGAL_ACK_USE_APPROXIMATE_ROTATION
#ifndef CGAL_ACK_ANGLE_PRECISION
#define CGAL_ACK_ANGLE_PRECISION 8
#endif
#endif


// What is the coefficient type of the underlying basic kernel?
#ifndef CGAL_ACK_COEFFICIENT
#define CGAL_ACK_COEFFICIENT CGAL::CORE_arithmetic_kernel::Integer
#define CGAL_ACK_COEFFICIENT_IS_INTEGER 1
#endif

#include <CGAL/Timer.h>
CGAL::Timer overall_timer;

// Bitsize for random coefficients, if created
#define COEFFICIENT_BITSIZE 32

#include <sstream>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_curve_kernel_2_generator.h>

#if !CGAL_ACK_USE_APPROXIMATE_ROTATION
#include <CGAL/Rotated_algebraic_curve_kernel_2.h>
#endif

#include <CGAL/Polynomial_parser_2.h>

void print_parse_error(std::string str) {
    std::cout << "Interpreting " << str << " as a polynomial in MAPLE format, "
              << "parser reports an error" << std::endl;
}

void print_help(char* execname) {

    // TODO: Update!
    std::cout 
        << "Usage: " << execname 
        << " (RANDOM n)|(FILE filename)|(FILESEQ filename)|(INPUT input) angle"
        << std::endl << std::endl
        << "Takes a bivariate polynomial with integer coefficients, "
        << "interprets it as algebraic plane curve, rotates it by some angle, "
        << "and anlyses the rotated curve." << std::endl << std::endl
        << "if RANDOM is used, a random bivariate integer "
        << "polynomial of bidegree "
        << "n is generated" << std::endl
        << "if FILE is used, the file filename is taken as input. It is "
        << "expected that the file contains a curve as its head, the "
        << "remainder of the file is ignored"
        << std::endl
        << "if FILESEQ is used, the file filename is taken as input. it is"
        << " expected that each character 'P' marks the begin of a curve"
        << " which ends with the next endline symbol."
        << std::endl
        << "if INPUT is used, the next argument is taken as definition of "
        << "the polynomial" << std::endl << std::endl
        << "For the options FILE, FILESEQ and INPUT, the program expects the "
        << "polynomials to be in MAPLE format, " 
        << "with variables x and y, or in the EXACUS format P[...]"
        << " (for the FILESEQ option, the delimiter is '\\n'"
        << std::endl << std::endl
        << "angle defines the angle that the curve is rotated. It must be "
        << "a multiple of the compiler flag CGAL_ACK_BASE_ANGLE that can be "
        << "passed for the compilation." << std::endl << std::endl
        << "Example: " << std::endl << "\t" <<  execname << " RANDOM 4  "
        << "(analyses random polynomial of bidegree 4)" << std::endl
        << "\t" << execname 
        << " INPUT \"x^2+y^2-1\" "
        << "(analyses the unit circle)" << std::endl
        << "\t" << execname << " FILE poly1  (analyses the polynomial "
        << "defined in the textfile poly1)" << std::endl <<std::endl;
}


template<class NT>
std::vector<NT> randvector(int degree,long bitsize) {
    std::vector<NT> coeffs(degree+1);
    for(int i=0;i<=degree;i++) {
        // Creating the coefficients VERY elementary...
        NT coeff=0;
        for(int j=0;j<bitsize-1;j++) {
            coeff = 2*coeff + (lrand48()%2);
        }
        // The last bit determines the sign
        if(lrand48()%2==0) {
            coeff=-coeff;
        }    
        coeffs[i]=coeff;
    }
    return coeffs;
}

int main(int argc,char** argv) {
  
    if(argc<4) {
        print_help(argv[0]);
        exit(-1);
    }

    ::CGAL::set_pretty_mode(std::cout);

    typedef CGAL_ACK_COEFFICIENT Coefficient;

    typedef CGAL::Algebraic_curve_kernel_2_generator<Coefficient>
        ::Algebraic_curve_kernel_with_qir_and_bitstream_2
        Basic_algebraic_curve_kernel_2;
    typedef Basic_algebraic_curve_kernel_2::Polynomial_2 Input_polynomial_2;
    typedef Basic_algebraic_curve_kernel_2::Polynomial_1 Input_polynomial_1;
    typedef Basic_algebraic_curve_kernel_2::Boundary Rational;
    

#if CGAL_ACK_USE_APPROXIMATE_ROTATION
    typedef Basic_algebraic_curve_kernel_2
        Rotated_algebraic_curve_kernel_2;
#else
    typedef CGAL::Rotation_traits_for_base_angle
        <Input_polynomial_2,CGAL_ACK_BASE_ANGLE> Rotation_traits;
    
    typedef CGAL::Rotated_algebraic_curve_kernel_2<Rotation_traits>
        Rotated_algebraic_curve_kernel_2;
#endif
    typedef Rotated_algebraic_curve_kernel_2::Curve_analysis_2 
        Curve_analysis_2;
    typedef Curve_analysis_2::Integer Integer;
  

    std::vector<Input_polynomial_2> curves;
    Input_polynomial_2 f;

    if(strcmp(argv[1],"FILE")!=0 
       && strcmp(argv[1],"RANDOM")!=0  
       && strcmp(argv[1],"INPUT")!=0 
       && strcmp(argv[1],"FILESEQ")!=0) {
        print_help(argv[0]);
        exit(-1);
    }

    if(strcmp(argv[1],"FILE")==0) {
        std::ifstream input(argv[2]);
        if(input.peek()=='P') {
            input >> f;
        } else {
            std::string str;
            std::getline(input,str);
            bool check = CGAL::Polynomial_parser_2<Input_polynomial_2>() 
                (str, f);
            if(! check) {
                print_parse_error(str);
                std::exit(-1);
            }
        }
        curves.push_back(f);
    }
    if(strcmp(argv[1],"FILESEQ")==0) {
        std::ifstream input(argv[2]);
        while(!input.eof()) {
            if(input.peek()=='P') {
                input >> f; 
                curves.push_back(f);
            } else {
                std::string str;
                std::getline(input,str);
                if(str.length()>0) {
                    CGAL::Polynomial_parser_2<Input_polynomial_2>() (str, f);
                    curves.push_back(f);
                }
            }
            
        }
    }
    if(strcmp(argv[1],"RANDOM")==0) {
#if CGAL_ACK_COEFFICIENT_IS_INTEGER 
        srand48(time(NULL));
        // Create random polynomial of given degree
     
        int degree = atoi(argv[2]);
        std::vector<Input_polynomial_1> coeffs;
        for(int i=0;i<=degree;i++) {
            std::vector<Coefficient> curr_coeffs 
                = randvector<Coefficient>(degree-i,COEFFICIENT_BITSIZE);
            coeffs.push_back(Input_polynomial_1
                             (curr_coeffs.begin(),curr_coeffs.end()));
        }
        f=Input_polynomial_2(coeffs.begin(),coeffs.end());
        curves.push_back(f);
#else
        std::cerr << "Choosen coefficient type does not allow "
                  << " creation of random numbers" << std::endl;
        std::exit(1);
#endif
    }
    if(strcmp(argv[1],"INPUT")==0) {
        std::stringstream ss(argv[2]);
        if(ss.peek()=='P') {
            ss >> f;
        } else {
            std::string str = ss.str();
            bool check = CGAL::Polynomial_parser_2<Input_polynomial_2>() 
                (str, f);
            if(! check) {
                print_parse_error(str);
                std::exit(-1);
            }
        }
        curves.push_back(f);
    }
    for(int i=0;i<static_cast<int>(curves.size());i++) {

        std::cout << "Input polynomial: " << curves[i] << std::endl;
            
       
#if CGAL_ACK_MAKE_SQUARE_FREE
        std::cout << "Make it square free..." << std::flush;
        f = CGAL::CGALi::make_square_free(curves[i]);
        std::cout << "done" << std::endl;
#else
        f = curves[i];
#endif

#if CGAL_ACK_USE_APPROXIMATE_ROTATION
        double angle_double = atof(argv[3]);

        Integer num,denom(1);
        while((double)((long)angle_double)!=angle_double) {
            std::cout << "angle_double=" << angle_double << std::endl;
            denom*=10;
            angle_double*=10;
        }
        num = Integer((long)angle_double);

        Rational angle = CGAL::Fraction_traits<Rational>::Compose()(num,denom);
#else
        int angle = atoi(argv[3]);
#endif

#if CGAL_ACK_USE_APPROXIMATE_ROTATION
        int prec;
        if(argc>4) {
            prec = atoi(argv[4]);
        } else {
            prec=CGAL_ACK_ANGLE_PRECISION;
        }
#endif 


        overall_timer.start();
        std::cout << "*************** Analysis starts ***************" << std::endl;
        Curve_analysis_2 algebraic_curve
#if CGAL_ACK_USE_APPROXIMATE_ROTATION
            = Rotated_algebraic_curve_kernel_2::Construct_curve_2()
            (f,angle,prec);
#else
        = Rotated_algebraic_curve_kernel_2::Construct_curve_2()(f,angle);
#endif

        std::cout << "Now refine..." << std::flush;
     
        algebraic_curve.refine_all(Rational(1,1000));

        overall_timer.stop();        

        std::cout << "done" << std::endl;

        std::cout << "*************** Analysis finished ***************" << std::endl;
     
        std::cout << "Overall time: " << overall_timer.time() << std::endl;
        std::cout << algebraic_curve;
     
    }
    return 0;
      
}
