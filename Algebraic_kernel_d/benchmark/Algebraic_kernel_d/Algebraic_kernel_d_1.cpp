// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Sebastian Limbach <slimbach@mpi-inf.mpg.de> 
//
// ============================================================================

// Benchmark of Algebraic_kernel

#define NiX_POLYNOMIAL_GCD_AVOID_CANONICALIZE 1
//#define CGAL_TEST_ONLY_SOLVE 1

#include <CGAL/basic.h>

#include <CGAL/Benchmark/Benchmark.hpp>
#include <CGAL/Benchmark/Option_parser.hpp>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Algebraic_kernel_1.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>

// TODO: Copied from benchmark_helper
template <class Integer> 
int get_max_bit_size(Integer x){
    return CGAL::CGALi::ceil_log2_abs(x);
}
   
template <class NT, class ROOT>
int get_max_bit_size(CGAL::Sqrt_extension<NT,ROOT> ext){
    int max = 0; 
    max = std::max(max,get_max_bit_size(ext.a0()));
    max = std::max(max,get_max_bit_size(ext.a1()));
    max = std::max(max,get_max_bit_size(ext.root()));
    return max; 
}

template <class NT>
int get_max_bit_size(CGAL::Polynomial<NT> poly){
    typedef CGAL::Polynomial<NT> Poly; 
    typedef CGAL::Polynomial_traits_d<Poly> PT;
    typename PT::Innermost_coefficient_begin begin;
    typename PT::Innermost_coefficient_end end;
    typedef typename PT::Innermost_coefficient_iterator IT;
    
    int max = 0; 
    for(IT it = begin(poly); it != end(poly); it++){
        max = std::max(max,get_max_bit_size(*it));
    }
    return max; 
}

template <class T>
int get_max_bit_size(std::vector<T> v) {
    int max = 0; 
    for(unsigned int i = 0; i < v.size(); i++){
        max = std::max(max,get_max_bit_size(v[i]));
    }       
    return max; 
}







struct Benchmark_result {
    int number_of_polys;
    double bits;
    int number_of_real_roots_found;
    int degree_of_polynomials;
    
    float solve_time;
    float solve_time_no_mult;
    float sort_time;
    float to_double_time;
    float total_time;
};

std::ostream& operator<<( std::ostream& os, const Benchmark_result& br ) {
    os << br.number_of_polys << "\t" << br.bits << "\t" << br.number_of_real_roots_found 
       << "\t" << br.degree_of_polynomials << "\t" << br.solve_time
       << "\t" << br.solve_time_no_mult << "\t" << br.sort_time
       << "\t" << br.to_double_time << "\t" << br.total_time;
    return os;
}

template< class AlgebraicKernel >
class Bench_solve_1 {
    private:
        typedef AlgebraicKernel AK;
        typedef std::vector< typename AK::Polynomial_1 >     Poly_vec;
        typedef std::vector< typename AK::Algebraic_real_1 > Root_vec;
        typedef std::vector< int >                           Mult_vec;
    
        typename Poly_vec::iterator polys_begin;
        typename Poly_vec::iterator polys_end;
        Root_vec* pRoot_vec;
        Mult_vec* pMult_vec;
        
    public:
        void prepare_op( typename Poly_vec::iterator polys_begin,
                         typename Poly_vec::iterator polys_end,
                         Root_vec* pRoot_vec,
                         Mult_vec* pMult_vec ) {
            this->polys_begin = polys_begin;
            this->polys_end = polys_end;
            this->pRoot_vec = pRoot_vec;
            this->pMult_vec = pMult_vec;
        }
        
        int init() { return 0; }
        void clean() {}
        void sync() {}
        void op() {
           // Clear for the case of multiple op calls
           pRoot_vec->clear();
           pMult_vec->clear();
           
           // Calculate roots and multiplicities of all polynomials
            typename AK::Solve_1 solve_1;
                
            for( typename Poly_vec::iterator it = polys_begin; it != polys_end; ++it ) {
	      typename AK::Polynomial_1 poly = (*it);
	      CGAL::remove_scalar_factor( poly );
	      solve_1( (*it), std::back_inserter( *pRoot_vec ), std::back_inserter( *pMult_vec ) );    
            }            
        }
};

template< class AlgebraicKernel >
class Bench_solve_1_no_mult {
    private:
        typedef AlgebraicKernel AK;
        typedef std::vector< typename AK::Polynomial_1 >     Poly_vec;
        typedef std::vector< typename AK::Algebraic_real_1 > Root_vec;
    
        typename Poly_vec::iterator polys_begin;
        typename Poly_vec::iterator polys_end;
        Root_vec* pRoot_vec;
          
    public:
        void prepare_op( typename Poly_vec::iterator polys_begin,
                         typename Poly_vec::iterator polys_end,
                         Root_vec* pRoot_vec ) {
            this->polys_begin = polys_begin;
            this->polys_end = polys_end;
            this->pRoot_vec = pRoot_vec;
        }
        
        int init() { return 0; }
        void clean() {}
        void sync() {}
        void op() {
           // Clear for the case of multiple op calls
           pRoot_vec->clear();
           
           // Calculate roots and multiplicities of all polynomials
            typename AK::Solve_1 solve_1;
                
            for( typename Poly_vec::iterator it = polys_begin; it != polys_end; ++it ) {
	      typename AK::Polynomial_1 poly = (*it);
	      CGAL::remove_scalar_factor( poly );
	      solve_1( (*it), std::back_inserter( *pRoot_vec ) );    
            }            
        }
};

template< class AlgebraicKernel >
class Bench_sort {
    private:
        typedef AlgebraicKernel AK;
        typedef std::vector< typename AK::Algebraic_real_1 > Root_vec;
    
        Root_vec* pRoot_vec;
        
    public:
        void prepare_op( Root_vec* pRoot_vec ) {
            this->pRoot_vec = pRoot_vec;
        }
        
        int init() { return 0; }
        void clean() {}
        void sync() {}
        void op() {
            // Copy array of roots for the case of multiple op calls
//            Root_vec root_vec_copy = *pRoot_vec;
            // Sorting the roots
            std::sort( pRoot_vec->begin(), pRoot_vec->end() );
        }
};

template< class AlgebraicKernel >
class Bench_to_double {
    private:
        typedef AlgebraicKernel AK;
        typedef std::vector< typename AK::Algebraic_real_1 > Root_vec;
    
        Root_vec roots;
        
    public:
        void prepare_op( Root_vec roots ) {
            this->roots = roots;
        }
        
        int init() { return 0; }
        void clean() {}
        void sync() {}
        void op() {
            // Copy array of roots for the case of multiple op calls
//            Root_vec root_vec_copy = roots;
            // Calling to_double to force the refinement of the intervals to 53 bit precission
            for( typename Root_vec::iterator rit = roots.begin();
                 rit != roots.end(); ++rit ) {
                CGAL::to_double( (*rit) );
            }
        }
};

template< class AlgebraicKernel >
Benchmark_result do_benchmark( std::string filename, int samples = 5 ) {
    typedef AlgebraicKernel AK;
    typedef std::vector< typename AK::Polynomial_1 >     Poly_vec;
    typedef std::vector< typename AK::Algebraic_real_1 > Root_vec;
    typedef std::vector< int >                           Mult_vec;
    
    Poly_vec   polys;
    Root_vec   roots;
    Mult_vec   mults;
    
    Benchmark_result result;
    
    CGAL::benchmark::Benchmark< Bench_solve_1< AK > > bench_solve_1( filename, 0, true );
    CGAL::benchmark::Benchmark< Bench_solve_1_no_mult< AK > > bench_solve_1_no_mult( filename, 0, true );
    CGAL::benchmark::Benchmark< Bench_sort< AK > > bench_sort( filename, 0, true );
    CGAL::benchmark::Benchmark< Bench_to_double< AK > > bench_to_double( filename, 0, true );
    bench_solve_1.set_samples( samples );
    bench_solve_1_no_mult.set_samples( samples );
    bench_sort.set_samples( samples );
    bench_to_double.set_samples( samples );
    
    // Read in the polynomials for testing
    std::ifstream file( filename.c_str() );
     
    CGAL_assertion_msg( file, (std::string("File not found: ") + filename).c_str() );
     
    int numPolys;
    file >> numPolys; 
          
    for( int i = 0; i < numPolys; ++i ) {
       typename AK::Polynomial_1 poly;
       file >> poly;
       CGAL::remove_scalar_factor( poly );
       polys.push_back( poly );
       CGAL_assertion( !file.eof() );
    }
          
    file.close();       
    
    result.number_of_polys = numPolys;    
    result.degree_of_polynomials = polys.begin()->degree();
    
    double total_bits = 0.0;
    for( typename Poly_vec::iterator poly_it = polys.begin(); poly_it != polys.end(); ++poly_it )
        total_bits += get_max_bit_size( (*poly_it) );
        
    result.bits = total_bits / numPolys;
        
    // Bench Solve_1
    bench_solve_1.get_benchable().prepare_op( polys.begin(), polys.end(), &roots, &mults );
    bench_solve_1();
    result.solve_time = bench_solve_1.get_period() / samples;
    
    result.number_of_real_roots_found = roots.size();

    // Bench Solve_1_no_mult
    bench_solve_1_no_mult.get_benchable().prepare_op( polys.begin(), polys.end(), &roots );
    bench_solve_1_no_mult();
    result.solve_time_no_mult = bench_solve_1_no_mult.get_period() / samples;
    
    result.number_of_real_roots_found = roots.size();
    
    // Bench Sort
#ifndef CGAL_TEST_ONLY_SOLVE    
    bench_sort.get_benchable().prepare_op( &roots );
    bench_sort();
    result.sort_time = bench_sort.get_period() / samples;
    
    // Bench to_double
    bench_to_double.get_benchable().prepare_op( roots );
    bench_to_double();
    result.to_double_time = bench_to_double.get_period() / samples;
#else
    result.sort_time = 0;
    result.to_double_time = 0;
#endif

    result.total_time = result.solve_time + result.solve_time_no_mult + result.sort_time + result.to_double_time;
    
    return result;
}


template< class Coeff_, class Boundary_, class RepClass, class Isolator_ >
void single_benchmark( std::string filename, int samples = 5 ) {
    typedef Coeff_      Coeff;
    typedef Boundary_   Boundary;
    typedef RepClass    Rep_class;
    typedef Isolator_   Isolator;
    typedef CGAL::Algebraic_kernel_1< Coeff, Boundary, Rep_class, Isolator > AK;

    // Output result to cerr
    // I'm using cerr because the benchmark results are written to cout.
    // One can stream them into a result file e.g. with 2>>
    std::cerr << do_benchmark< AK >( filename, samples ) << std::endl;
}

template< class Coeff_, class Boundary_, class RepClass >
void single_benchmark( std::string filename, std::string isolator, int samples = 5 ) {
    typedef Coeff_      Coeff;
    typedef Boundary_   Boundary;
    typedef RepClass    Rep_class;

    // Select isolator
    if( isolator == "Descartes" )
        single_benchmark< Coeff, Boundary, Rep_class,
            CGAL::CGALi::Descartes< typename CGAL::Polynomial< Coeff >, Boundary > >( filename, samples );
    else if( isolator == "Bitstream_descartes" )
        single_benchmark< Coeff, Boundary, Rep_class,
            CGAL::CGALi::Bitstream_descartes< typename CGAL::Polynomial< Coeff >, Boundary > >( filename, samples );
    else
        CGAL_error_msg( "Unknown isolator class" );
}

template< class Coeff_ >
void single_benchmark( std::string filename, std::string rep_class, std::string isolator, int samples = 5 ) {
    typedef Coeff_ Coeff;
    typedef typename CGAL::Get_arithmetic_kernel< Coeff >::Arithmetic_kernel::Rational Boundary;

    // Select rep class    
    if( rep_class == "Algebraic_real_rep" )
        single_benchmark< Coeff, Boundary,
             CGAL::CGALi::Algebraic_real_rep< Coeff, Boundary > >( filename, isolator, samples );
    else if( rep_class == "Algebraic_real_rep_bfi" )
        single_benchmark< Coeff, Boundary,
             CGAL::CGALi::Algebraic_real_rep_bfi< Coeff, Boundary > >( filename, isolator, samples );
    else if( rep_class == "Algebraic_real_quadratic_refinement_rep_bfi" )
        single_benchmark< Coeff, Boundary,
             CGAL::CGALi::Algebraic_real_quadratic_refinement_rep_bfi< Coeff, Boundary > >( filename, isolator, samples );
    else
         CGAL_error_msg( "Unknown rep class" );
}

int main( int argc, char** argv ) {

    typedef CGAL::CORE_arithmetic_kernel Arithmetic_kernel;
    typedef Arithmetic_kernel::Integer Integer;

    if( argc > 1 ) {
    
        CGAL_assertion( argc >= 5 );
        int samples = (argc == 6 ) ? std::atoi( argv[5] ) : 5;
        if( std::string( argv[1] ) == "INT" )
            single_benchmark< Integer >( argv[4], argv[2], argv[3], samples );
        else if( std::string( argv[1] ) == "EXT" ) 
            single_benchmark< CGAL::Sqrt_extension< Integer, Integer > >
                ( argv[4], argv[2], argv[3], samples );
        else
            CGAL_error_msg( "Unknown coefficient type" );
        
    } else {
        std::cerr << "No parameters found" << std::endl;    
    }
            
    return 0;    
}
