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

#ifndef CGAL_ACK_DEGENERACY_STRATEGY
#define CGAL_ACK_DEGENERACY_STRATEGY 1

CGAL_BEGIN_NAMESPACE

/*!
 * \brief Represents different strategies how to handle 
 * degenerate cases during the analysis
 *
 * Currently, there are two possible strategies implemented. See the 
 * constructor of \c Curve_analysis_2 for more details.
 */
enum Degeneracy_strategy {
    
    SHEAR_STRATEGY = 0,
    EXCEPTION_STRATEGY = 1,
    SHEAR_ONLY_AT_IRRATIONAL_STRATEGY = 2
};

CGAL_END_NAMESPACE

#endif
