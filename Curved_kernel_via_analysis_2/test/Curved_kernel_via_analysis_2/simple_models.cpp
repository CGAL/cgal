// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#include <vector>
#include <sstream>

using std::output_iterator_tag; // compiler complains

// not compiled with preconditions (NumeriX library conflicts)
#define CGAL_NO_PRECONDITIONS

#include <CGAL/basic.h>

#ifdef CGAL_USE_LEDA
#undef CGAL_USE_LEDA // doesn't compile with leda ?
#endif

// temporarily required while CKvA depends on NumeriX library
#include <NiX/resultant.h>

#include <CGAL/Curved_kernel_via_analysis/test/simple_models.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Curved_kernel_via_analysis_2.h>

int main(int argc, char** argv) {

    typedef CGAL::Simple_algebraic_kernel_2 Kernel_2;

     // GPA arrangement traits
    typedef CGAL::Curved_kernel_via_analysis_2<Kernel_2> Traits;
    
    typedef Traits::Point_2 Point_2;
    typedef Traits::X_monotone_curve_2 Arc_2;
    Kernel_2 kernel_2;

    std::vector<Kernel_2::Curve_2> curves(2);

    typedef CGAL::Arrangement_2<Traits> Arrangement;
    Arrangement arrangement;    
    
    CGAL::insert(arrangement, curves.begin(), curves.end());     
 
    return 0;
}
