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
#define CGAL_ACK_DEBUG_FLAG 0
#endif

#ifndef CGAL_ACK_DEBUG_PRINT
#define CGAL_ACK_DEBUG_PRINT std::cout
#endif

#include <CGAL/basic.h>

#ifndef CGAL_ACK_USE_BEZOUT_MATRIX_FOR_SUBRESULTANTS
#define CGAL_ACK_USE_BEZOUT_MATRIX_FOR_SUBRESULTANTS 0
#endif

#ifndef CGAL_ACK_RESULTANT_FIRST_STRATEGY
#define CGAL_ACK_RESULTANT_FIRST_STRATEGY 0
#endif

#ifdef CGAL_HAVE_CORE
#ifndef CGAL_ACK_USE_CORE
#define CGAL_ACK_USE_CORE 1
#endif
#else
#ifdef CGAL_USE_LEDA
#ifndef CGAL_ACK_USE_LEDA
#define CGAL_ACK_USE_LEDA 1
#endif
#endif
#endif

#define COEFFICIENT_BITSIZE 8

#include <sstream>

#include <CGAL/Timer.h>
CGAL::Timer overall_timer;
CGAL::Timer first_curve_timer, second_curve_timer,
  pair_timer;

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_1.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_curve_kernel_2.h>

const char* name="Curve_pair_analyser";


void print_help() {
  std::cout << "Usage: " << name 
	    << " (RANDOM n)|(FILE filename)|(FILESEQ filename)|(INPUT input)"
	    << std::endl << std::endl
	    << "if RANDOM is used, a random bivariate polynomial of bidegree "
	    << "n is generated and analysed." << std::endl
	    << "if FILE is used, the file filename is taken as input. It is "
	    << "expected that the file contains a curve as its head, and the "
	    << "remainder of the file is ignored"
	    << std::endl
	    << "if FILESEQ is used, the file filename is taken as input. it is"
	    << " expected that each character 'P' marks the begin of a curve"
	    << " which ends with the next endline symbol."
	    << std::endl
	    << "if INPUT is used, the next argument is taken as definition of "
	    << "the polynomial" << std::endl << std::endl
	    << "For the options FILE and INPUT, the program expects the "
	    << "polynomials to be of the pattern P[...], according to the "
	    << "EXACUS format" << std::endl << std::endl
	    << "A file called \"gen_plot\" is created in the same "
	    << "directory giving a (very primitive) topological description "
	    << "of the polynomial. "
	    << "Try \"graph -T X gen_plot\" to display the graph "
	    << "on X-window systems" << std::endl << std::endl
	    << "Example: " << std::endl << "\t" <<  name << " RANDOM 4  "
	    << "(analyses random polynomial of bidegree 4)" << std::endl
	    << "\t" << name 
	    << " INPUT \"P[2(0,P[2(0,-1)(2,1)])(2,P[0(0,1)])]\" "
	    << "(analyses the unit circle)" << std::endl
	    << "\t" << name << " FILE poly1  (analyses the polynomial "
	    << "defined in the textfile poly1)" << std::endl <<std::endl;

}

void reset_timers() {
  overall_timer.reset();
  first_curve_timer.reset();
  second_curve_timer.reset();
  pair_timer.reset();
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
  
  if(argc<3) {
    print_help();
    exit(-1);
  }

  ::CGAL::set_pretty_mode(std::cout);

#if CGAL_ACK_USE_CORE
  //std::cout << "Use CORE library" << std::endl;
  typedef CGAL::CORE_arithmetic_kernel Arithmetic_kernel;
#elif defined(CGAL_ACK_USE_LEDA)
  //std::cout << "Use LEDA library" << std::endl;
  typedef CGAL::LEDA_arithmetic_kernel Arithmetic_kernel;
#else
  std::cerr << "CORE or LEDA required!" << std::endl;
  std::exit(1);
#endif

  typedef Arithmetic_kernel::Rational Rational;
  typedef Arithmetic_kernel::Integer Integer;
  
  typedef Integer Coefficient;
  typedef CGAL::Polynomial<Coefficient> Poly_int1;
  typedef CGAL::Polynomial<Poly_int1> Poly_int2;
  
  typedef CGAL::CGALi::Algebraic_real_quadratic_refinement_rep_bfi
      < Coefficient, Rational > Rep_class;
  typedef CGAL::CGALi::Bitstream_descartes< CGAL::Polynomial< Coefficient >, 
      Rational > Isolator;
  
  typedef CGAL::Algebraic_kernel_1<Coefficient,Rational,Rep_class, Isolator> 
      Algebraic_kernel_1;

  typedef CGAL::Algebraic_curve_kernel_2<Algebraic_kernel_1>
      Algebraic_curve_kernel_2;
  
  typedef Algebraic_curve_kernel_2::Curve_analysis_2 
      Curve_analysis_2;
  
  typedef Algebraic_curve_kernel_2::Curve_pair_analysis_2 
      Curve_pair_analysis_2;
  
  std::vector<std::pair<Poly_int2,Poly_int2> > curve_pairs;
  Poly_int2 f,g;
  
  if(strcmp(argv[1],"FILE")!=0 
     && strcmp(argv[1],"RANDOM")!=0  
     && strcmp(argv[1],"INPUT")!=0 
     && strcmp(argv[1],"FILESEQ")!=0) {
    print_help();
    exit(-1);
  }

   if(strcmp(argv[1],"FILE")==0) {
    std::ifstream input(argv[2]);
    input >> f; 
    input >> g;
    curve_pairs.push_back(std::make_pair(f,g));
   }
   /* THINK: Include this later
   if(strcmp(argv[1],"FILESEQ")==0) {
     std::ifstream input(argv[2]);
     while(!input.eof()) {
       int g = input.get();
       if(g=='P') {
	 input.putback(g);
	 input >> f; 
	 curve_pairs.push_back(f);
       }
     }
   }
   */
   if(strcmp(argv[1],"RANDOM")==0) {
     srand48(time(NULL));
     // Create random polynomial of given degree
     
     int degree1 = atoi(argv[2]);
     int degree2 = atoi(argv[3]);
     std::vector<Poly_int1> coeffs;
     for(int i=0;i<=degree1;i++) {
       std::vector<Integer> curr_coeffs 
	 = randvector<Integer>(degree1-i,COEFFICIENT_BITSIZE);
       coeffs.push_back(Poly_int1(curr_coeffs.begin(),curr_coeffs.end()));
     }
     f=Poly_int2(coeffs.begin(),coeffs.end());
     coeffs.clear();
     for(int i=0;i<=degree2;i++) {
       std::vector<Integer> curr_coeffs 
	 = randvector<Integer>(degree2-i,COEFFICIENT_BITSIZE);
       coeffs.push_back(Poly_int1(curr_coeffs.begin(),curr_coeffs.end()));
     }
     g=Poly_int2(coeffs.begin(),coeffs.end());
     curve_pairs.push_back(std::make_pair(f,g));
   }
   if(strcmp(argv[1],"INPUT")==0) {
     std::stringstream ss1(argv[2]);
     ss1 >> f;
     std::stringstream ss2(argv[3]);
     ss2 >> g;
     curve_pairs.push_back(std::make_pair(f,g));
   }
   for(int i=0;i<static_cast<int>(curve_pairs.size());i++) {
     f = curve_pairs[i].first;
     g = curve_pairs[i].second;
     std::cout << "Input polynomials:\n" << "f=" << f << "\ng=" << g 
               << std::endl;
     reset_timers();
     overall_timer.start();
     std::cout << "+++++++++++++++ One curve analysis starts +++++++++++" << std::endl;
     first_curve_timer.start();
     Curve_analysis_2 F(f);
     first_curve_timer.stop();
     second_curve_timer.start();
     Curve_analysis_2 G(g);
     second_curve_timer.stop();
     std::cout << "+++++++++++++++ One curve analysis ends +++++++++++" << std::endl;
     std::cout << "*************** Two  curve analysis starts ***************" << std::endl;
     pair_timer.start();
     Curve_pair_analysis_2 algebraic_curve_pair(F,G);
     pair_timer.stop();
     overall_timer.stop();   
     std::cout << "*************** Two curve analysis ends ***************" << std::endl;

     std::cout << algebraic_curve_pair << std::endl;

     std::cout << "TIMINGS: " << std::endl;
     std::cout << "Analysing 1st curve:  " << first_curve_timer.time() 
               << std::endl; 
     std::cout << "Analysing 2nd curve:  " << second_curve_timer.time() 
               << std::endl;
     std::cout << "Analysing curve pair: " << pair_timer.time() << std::endl; 
     std::cout << "----" << std::endl;
     std::cout << "Overall timer:        " << overall_timer.time() << std::endl; 

   }
   return 0;
      
}

