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
#include <CGAL/basic.h>

#include <CGAL/Benchmark/Benchmark.hpp>
#include <CGAL/Benchmark/Option_parser.hpp>

#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>


template< class AlgebraicKernel >
class Bench_algebraic_kernel_d_1 {
private:
    typedef AlgebraicKernel AK;
    typedef std::vector< typename AK::Polynomial_1 >     Poly_vec;
    typedef std::vector< typename AK::Algebraic_real_1 > Root_vec;
    typedef std::vector< int >                           Mult_vec;
    
    std::string filename;
    
    Poly_vec   polys;

public:
    int num_roots;

  int init(void) { 
     // Read in the polynomials for testing
     std::ifstream file( filename.c_str() );
     
     CGAL_assertion_msg( file, (std::string("File not found: ") + filename).c_str() );
     
     int numPolys;
     file >> numPolys; 
          
     for( int i = 0; i < numPolys; ++i ) {
        typename AK::Polynomial_1 poly;
        file >> poly;
        polys.push_back( poly );
        CGAL_assertion( !file.eof() );
     }
          
     file.close();       
    
     return 0; 
  }
  
  void clean(void) { 
    polys.clear();
  }
  
  void sync(void) {
  }
  
  void op(void) {
   // Calculate roots and multiplicities of all polynomials
    typename AK::Solve_1 solve_1;

    Root_vec   roots;
    Mult_vec   mults;

    for( typename Poly_vec::iterator it = polys.begin(); it != polys.end(); ++it ) {
        typename AK::Polynomial_1 poly = (*it);
        CGAL::remove_scalar_factor( poly );
        solve_1( poly, std::back_inserter( roots ), std::back_inserter( mults ) );    
    }

    // Sorting the roots
    std::sort( roots.begin(), roots.end() );
    
    // Calling to_double to force the refinement of the intervals to 53 bit precission
    for( typename Root_vec::iterator rit = roots.begin(); rit != roots.end(); ++rit ) {
        CGAL::to_double( (*rit) );
    }
    
    // Save the number of roots to correct the results
    num_roots = roots.size(); 
    
    return; 
  }
  
  void set_filename( std::string filename ) {
    this->filename = filename;
  }

};

void do_benchmark( int numQuadrics, int bitsFrom, int bitsTo, int samples = 5 ) {
    typedef CGAL::benchmark::Benchmark< Bench_algebraic_kernel_d_1< CGAL::Algebraic_kernel_d_1< CGAL::Sqrt_extension< CORE::BigInt, CORE::BigInt > > > > Bench;

    for( int bits = bitsFrom; bits <= bitsTo; bits += 10 ) {
        std::stringstream filename;
        filename.fill( '0' );
        filename << "data/resultants_from_" 
                 << std::setw(4) << numQuadrics
                 << "_quadrics_with_" 
                 << std::setw(4) << bits << "_bits.nix";
        
        std::stringstream benchmarkname;
        benchmarkname << std::setw(4) << numQuadrics << " quadrics / " << bits << " bits"; 
        
        Bench bench( benchmarkname.str(), 0, bits == bitsFrom );
        bench.set_samples( samples );     
        bench.get_benchable().set_filename( filename.str() );
        bench();
    }    
};

template< class Coeff_, class Boundary_, class RepClass, class Isolator_ >
void single_benchmark( std::string filename, int samples = 5 ) {
    typedef Coeff_      Coeff;
    typedef Boundary_   Boundary;
    typedef RepClass    Rep_class;
    typedef Isolator_   Isolator;
    typedef CGAL::Algebraic_kernel_d_1< Coeff, Boundary, Rep_class, Isolator > AK;
    typedef CGAL::benchmark::Benchmark< Bench_algebraic_kernel_d_1< AK > > Bench;

    Bench bench( filename, 0, true );
    bench.set_samples( samples );
    bench.get_benchable().set_filename( filename );
    bench();

    // Output result to cerr
    // I'm using cerr because the benchmark results are written to cout.
    // One can stream them into a result file e.g. with 2>>
    std::cerr << (bench.get_period() / samples / bench.get_benchable().num_roots) << std::endl;
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
        CGAL_error( "Unknown isolator class" );
}

template< class Coeff_ >
void single_benchmark( std::string filename, std::string rep_class, std::string isolator, int samples = 5 ) {
    typedef Coeff_ Coeff;
    typedef typename CGAL::CGALi::Get_arithmetic_kernel< Coeff >::Arithmetic_kernel::Rational Boundary;

    // Select rep class    
    if( rep_class == "Algebraic_real_rep" )
        single_benchmark< Coeff, Boundary,
             CGAL::CGALi::Algebraic_real_rep< Coeff, Boundary > >( filename, isolator, samples );
    else if( rep_class == "Algebraic_real_rep_bfi" )
        single_benchmark< Coeff, Boundary,
             CGAL::CGALi::Algebraic_real_rep_bfi< Coeff, Boundary > >( filename, isolator, samples );
    else
         CGAL_error( "Unknown rep class" );
}

int main( int argc, char** argv ) {
    if( argc > 1 ) {
    
        CGAL_assertion( argc >= 5 );
        int samples = (argc == 6 ) ? std::atoi( argv[5] ) : 5;
            
        if( std::string( argv[1] ) == "INT" )
            single_benchmark< leda_integer >( argv[4], argv[2], argv[3], samples );
//            single_benchmark< CORE::BigInt >( argv[4], argv[2], argv[3], samples );
        else if( std::string( argv[1] ) == "EXT" ) 
            single_benchmark< CGAL::Sqrt_extension< leda_integer, leda_integer > >( argv[4], argv[2], argv[3], samples );
//            single_benchmark< CGAL::Sqrt_extension< CORE::BigInt, CORE::BigInt > >( argv[4], argv[2], argv[3], samples );
        else
            CGAL_error( "Unknown coefficient type" );
    
    } else {
    
        for( int num_quadrics = 10; num_quadrics <= 90; num_quadrics += 10 ) {
            do_benchmark( num_quadrics, 10, 90 );
        }

    }
            
    return 0;    
}
