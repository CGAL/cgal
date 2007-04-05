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

#include <CGAL/basic.h>

#include <CGAL/Benchmark/Benchmark.hpp>
#include <CGAL/Benchmark/Option_parser.hpp>

#include <CGAL/Algebraic_kernel_d_1.h>

template< class AlgebraicKernel >
class Bench_algebraic_kernel_d_1 {
private:
    typedef AlgebraicKernel AK;
    typedef std::vector< typename AK::Polynomial_1 >     Poly_vec;
    typedef std::vector< typename AK::Algebraic_real_1 > Root_vec;
    typedef std::vector< int >                           Mult_vec;
    
    std::string filename;
    
    Poly_vec   polys;
    Root_vec   roots;
    Mult_vec   mults;
public:
  int init(void) { 
     // Read in the polynomials for testing
     std::ifstream file( filename.c_str() );
     
     CGAL_assertion_msg( file, "File not found" );
     
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
    roots.clear();
    mults.clear(); 
  }
  
  void sync(void) {
  }
  
  void op(void) {    
    // Calculate roots and multiplicities of all polynomials
    typename AK::Solve_1 solve_1;

    for( typename Poly_vec::iterator it = polys.begin(); it != polys.end(); ++it ) {
        solve_1( (*it), std::back_inserter( roots ), std::back_inserter( mults ) );    
    }
    
    // Sorting the roots
    std::sort( roots.begin(), roots.end() );
    
    return; 
  }
  
  void set_filename( std::string filename ) {
    this->filename = filename;
  }

};

void do_benchmark( int numQuadrics, int bitsFrom, int bitsTo, int seconds = 2 ) {
    typedef CGAL::benchmark::Benchmark< Bench_algebraic_kernel_d_1< CGAL::Algebraic_kernel_d_1< CGAL::Sqrt_extension< leda_integer, leda_integer > > > > Bench;

    for( int bits = bitsFrom; bits <= bitsTo; bits += 10 ) {
        std::stringstream filename;
        filename.fill( '0' );
        filename << "data/resultants_from_" 
                 << std::setw(4) << numQuadrics
                 << "_quadrics_with_" 
                 << std::setw(4) << bits << "_bits.nix";
        
        std::stringstream benchmarkname;
        benchmarkname << std::setw(4) << numQuadrics << " quadrics / " << bits << " bits"; 
        
        Bench bench( benchmarkname.str(), seconds, bits == bitsFrom );
        bench.set_iterations( 1 );     
        bench.get_benchable().set_filename( filename.str() );
        bench();
    }
    
};

int main( int argc, char** argv ) {
    do_benchmark( 10, 10, 90 );
        
    return 0;    
}
