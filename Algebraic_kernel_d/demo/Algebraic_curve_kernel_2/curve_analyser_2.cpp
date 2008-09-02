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

// Does the demo create a file "gen_plot" 
// that displays the input curves topology?
#ifndef CGAL_ACK_WRITE_GRAPH_TO_FILE
#define CGAL_ACK_WRITE_GRAPH_TO_FILE 1
#endif

// Is the input curve made square free beforehand, 
// or does the program rely on a square free input?
#ifndef CGAL_ACK_MAKE_SQUARE_FREE
#define CGAL_ACK_MAKE_SQUARE_FREE 0
#endif

// What is the coefficient type of the input?
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

#include <CGAL/Polynomial_parser_2.h>

void print_parse_error(std::string str) {
    std::cout << "Interpreting " << str << " as a polynomial in MAPLE format, "
              << "parser reports an error" << std::endl;
}

void print_help(char* execname) {
    std::cout 
        << "Usage: " << execname 
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
        << "For the options FILE, FILESEQ and INPUT, the program expects the "
        << "polynomials to be in MAPLE format, " 
        << "with variables x and y, or in the EXACUS format P[...]"
        << " (for the FILESEQ option, the delimiter is '\\n'"
        << std::endl << std::endl
        << "A file called \"gen_plot\" is created in the same "
        << "directory giving a (very primitive) topological description "
        << "of the polynomial. "
        << "Try \"graph -T X gen_plot\" to display the graph "
        << "on X-window systems" << std::endl << std::endl
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

// Finds good y-values for placing the vertical asymptotes in the drawing
template<typename Rational,
	 typename InputIterator1,
	 typename InputIterator2>
void find_asymptote_places(double& minus,double& plus,
			   InputIterator1 event_begin,
			   InputIterator1 event_end,
			   InputIterator2 intermediate_begin,
			   InputIterator2 intermediate_end) {
    plus=0;
    minus=0;
    for(InputIterator1 curr=event_begin;curr!=event_end;curr++) {
        int n = curr->number_of_events();
        if(n>0) {
            if(minus > CGAL::to_double(curr->lower_boundary(0))) {
                minus= CGAL::to_double(curr->lower_boundary(0));
            }
            if(plus < CGAL::to_double(curr->upper_boundary(n-1))) {
                plus=CGAL::to_double(curr->upper_boundary(n-1));
            }
        }
    }
    for(InputIterator2 curr=intermediate_begin;curr!=intermediate_end;curr++) {
        int n = curr->number_of_events();
        if(n>0) {
            if(minus > CGAL::to_double(curr->lower_boundary(0))) {
                minus=CGAL::to_double(curr->lower_boundary(0));
            }
            if(plus < CGAL::to_double(curr->upper_boundary(n-1))) {
                plus=CGAL::to_double(curr->upper_boundary(n-1));
            }
        }
    }
    minus-=1;
    plus+=1;
    return;
}


template<typename Algebraic_real,
	 typename Rational,
	 typename LineType1,
	 typename LineType2> 
void connect (std::ofstream& out,
	      LineType1 line1,
	      LineType2 line2,
	      double asym_y_value_minus,
	      double asym_y_value_plus) {
    // Quick and dirty...
    std::vector<int> indices1,indices2;
  
    int lower_asym,upper_asym;
  
    typedef typename LineType1::Arc_pair Arc_pair1;
  
    Arc_pair1 apair11 = line1.number_of_branches_approaching_minus_infinity();
    Arc_pair1 apair12 = line1.number_of_branches_approaching_plus_infinity();
  
    lower_asym = apair11.second;

    upper_asym = apair12.second;

    // Asymptotes from -infty:
    for(int i=0;i<lower_asym;i++) {
        indices1.push_back(-1);
    }
    for(int i=0;i<line1.number_of_events();i++) {
        if(line1.is_event(i)) {
            for(int j = 0;j<line1.number_of_incident_branches(i).second;j++) {
                indices1.push_back(i);
            }
        }
        else {
            indices1.push_back(i);
        }
    }
    // Asymptotes from +infty:
    for(int i=0;i<upper_asym;i++) {
        indices1.push_back(-2);
    }

    typedef typename LineType2::Arc_pair Arc_pair2;

    Arc_pair2 apair21 = line2.number_of_branches_approaching_minus_infinity();
    Arc_pair2 apair22 = line2.number_of_branches_approaching_plus_infinity();

    lower_asym = apair21.first;
    upper_asym = apair22.first;

    // Asymptotes from -infty:
    for(int i=0;i<lower_asym;i++) {
        indices2.push_back(-1);
    }
    for(int i=0;i<line2.number_of_events();i++) {
        if(line2.is_event(i)) {
            for(int j = 0;j<line2.number_of_incident_branches(i).first;j++) {
                indices2.push_back(i);
            }
        }
        else {
            indices2.push_back(i);
        }
    }
    // Asymptotes from +infty:
    for(int i=0;i<upper_asym;i++) {
        indices2.push_back(-2);
    }
    //  std::cout << indices1.size() << " " << indices2.size() << std::endl;
    CGAL_assertion(indices1.size()==indices2.size());
    double xl = CGAL::to_double(line1.x());
    double xr = CGAL::to_double(line2.x());
    for(int i = 0;i<(int)indices1.size();i++) {
        double yl;
        int ind1=indices1[i];
        if(ind1==-1) {
            yl=asym_y_value_minus;
        }
        else {
            if(ind1==-2) {
                yl=asym_y_value_plus;
            }
            else {
                yl = (CGAL::to_double(line1.upper_boundary(indices1[i]))
                      + CGAL::to_double(line1.upper_boundary(indices1[i])))/2.0;
            }
        }
    
        double yr;
        int ind2=indices2[i];
        if(ind2==-1) {
            yr=asym_y_value_minus;
        }
        else {
            if(ind2==-2) {
                yr=asym_y_value_plus;
            }
            else {
                yr = (CGAL::to_double(line2.upper_boundary(indices2[i]))
                      + CGAL::to_double(line2.upper_boundary(indices2[i])))/2.0;
            }
        }
        out << "#m=1,S=0\n" << xl << " " << yl << "\n" << xr << " " << yr << "\n";
    
    }
}


int main(int argc,char** argv) {
  
    if(argc<3) {
        print_help(argv[0]);
        exit(-1);
    }

    ::CGAL::set_pretty_mode(std::cout);

    typedef CGAL_ACK_COEFFICIENT Coefficient;

    typedef CGAL::Algebraic_curve_kernel_2_generator<Coefficient>
        ::Algebraic_curve_kernel_with_qir_and_bitstream_2
        Algebraic_curve_kernel_2;

    typedef Algebraic_curve_kernel_2::Curve_analysis_2 Curve_analysis_2;
    typedef Algebraic_curve_kernel_2::Polynomial_2 Polynomial_2;
    typedef Algebraic_curve_kernel_2::Polynomial_1 Polynomial_1;
    typedef Algebraic_curve_kernel_2::Boundary Rational;
  

    std::vector<Polynomial_2> curves;
    Polynomial_2 f;

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
            bool check = CGAL::Polynomial_parser_2<Polynomial_2>() (str, f);
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
                    CGAL::Polynomial_parser_2<Polynomial_2>() (str, f);
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
        std::vector<Polynomial_1> coeffs;
        for(int i=0;i<=degree;i++) {
            std::vector<Coefficient> curr_coeffs 
                = randvector<Coefficient>(degree-i,COEFFICIENT_BITSIZE);
            coeffs.push_back(Polynomial_1(curr_coeffs.begin(),curr_coeffs.end()));
        }
        f=Polynomial_2(coeffs.begin(),coeffs.end());
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
            bool check = CGAL::Polynomial_parser_2<Polynomial_2>() (str, f);
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

        overall_timer.start();
        std::cout << "*************** Analysis starts ***************" << std::endl;
        Curve_analysis_2 algebraic_curve(f);

        std::cout << "Now refine..." << std::flush;
     
        algebraic_curve.refine_all(Rational(1,1000));

        overall_timer.stop();        

        std::cout << "done" << std::endl;

        std::cout << "*************** Analysis finished ***************" << std::endl;
     
        std::cout << algebraic_curve;

        std::cout << "Overall time: " << overall_timer.time() << std::endl;
     
#if CGAL_ACK_WRITE_GRAPH_TO_FILE
        // Now, create the plot:
        typedef Curve_analysis_2::Status_line_1 Event_line;
        typedef Curve_analysis_2::Status_line_1 Intermediate_line;
        typedef Curve_analysis_2::Event_line_iterator Event_iterator;
        typedef Curve_analysis_2::Intermediate_line_iterator 
            Intermediate_iterator;
        typedef Curve_analysis_2::Algebraic_real_1 Algebraic_real;
        std::ofstream out("gen_plot");
     
        // Event points
        out << "#m=0,S=3\n";
     
        for(Event_iterator ev_it = algebraic_curve.event_begin();
            ev_it!=algebraic_curve.event_end();ev_it++) {
            for(int i = 0; i<ev_it->number_of_events();i++) {
                if(ev_it->is_event(i)) {
                    double x=CGAL::to_double(ev_it->x());
                    double y=(CGAL::to_double(ev_it->upper_boundary(i))
                              +CGAL::to_double(ev_it->lower_boundary(i)))/2.0;
                    out <<  x << " " << y << "\n"; 
                }
            }
        }
        //Edges:
        double asym_y_value_plus,asym_y_value_minus;
   
        find_asymptote_places<Rational>
            (asym_y_value_minus,asym_y_value_plus,
             algebraic_curve.event_begin(),
             algebraic_curve.event_end(),
             algebraic_curve.intermediate_begin(),
             algebraic_curve.intermediate_end());
   
        Intermediate_iterator inter_it = algebraic_curve.intermediate_begin();
        Event_iterator ev_it = algebraic_curve.event_begin();
        while(inter_it!=algebraic_curve.intermediate_end() 
              && ev_it !=algebraic_curve.event_end()) {
            connect<Algebraic_real,Rational>(out,*inter_it,*ev_it,
                                             asym_y_value_minus,asym_y_value_plus);
            ++inter_it;
            connect<Algebraic_real,Rational>(out,*ev_it,*inter_it,
                                             asym_y_value_minus,asym_y_value_plus);
            ++ev_it;
        }
        out.close();
#endif
    }
    return 0;
      
}

