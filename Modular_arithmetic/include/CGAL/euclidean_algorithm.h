
// Author(s)     : Lutz Kettner   <kettner@mpi-inf.mpg.de>
// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>


/*! \file CGAL/euclidean_algorithm.h
    \brief Defines funciton related to euclids algorithm. 
*/

#ifndef CGAL_EUCLIDEAN_ALGORITHM_H
#define CGAL_EUCLIDEAN_ALGORITHM_H 1

// This forward declaration is required to resolve the circular dependency
// between euclidean_algorithm and the partial specializations of NT_Traits
// for the built-in number types.

namespace CGAL {
    template <class NT>
        NT euclidean_algorithm( const NT& a, const NT& b);
}

#include <CGAL/basic.h>
#include <CGAL/Algebraic_structure_traits.h>


namespace CGAL {

// We have a circular header file inclusion dependency with 
// CGAL/Algebraic_structure_traits.h.
// As a consequence, we might not get the declaration for 
// CGAL::Algebraic_structure_traits
// although we include the header file above. We repeat the declaration
// here. We still include the header file to hide this dependency from users
// such that they get the full Algebraic_structure_traits declaration after 
// including euclid_algorithm.h.

template <class NT> class Algebraic_structure_traits;


/*! \brief generic Euclids algorithm, returns the unit 
    normal greatest common  devisor (gcd) of \a a and \a b.

    Requires the number type \c NT to be a model of the concepts
    \c EuclideanRing, however, it uses only
    the functors \c Mod and \c Unit_part from the \c 
    Algebraic_structure_traits<NT>, and the equality comparison operator. 
    The implementation uses loop unrolling to avoid swapping the local 
    variables all the time.
*/
template <class NT>
NT euclidean_algorithm( const NT& a, const NT& b) {
    typedef Algebraic_structure_traits<NT> AST;
    typename AST::Mod mod;
    typename AST::Unit_part unit_part;
    typename AST::Integral_division idiv; 
    // First: the extreme cases and negative sign corrections.
    if (a == NT(0)) {
        if (b == NT(0))  
            return NT(0);
        return idiv(b,unit_part(b));
    }
    if (b == NT(0))
        return idiv(a,unit_part(a));
    NT u = idiv(a,unit_part(a));
    NT v = idiv(b,unit_part(b));
    // Second: assuming mod is the most expensive op here, we don't compute it
    // unnecessarily if u < v
    if (u < v) {
        v = mod(v,u);
        // maintain invariant of v > 0 for the loop below
        if ( v == 0)
            return idiv(u,unit_part(u));
    }
    // Third: generic case of two positive integer values and u >= v.
    // The standard loop would be:
    //      while ( v != 0) {
    //          int tmp = mod(u,v);
    //          u = v;
    //          v = tmp;
    //      }
    //      return u;
    //
    // But we want to save us all the variable assignments and unroll
    // the loop. Before that, we transform it into a do {...} while()
    // loop to reduce branching statements.
    NT w;
    do {
        w = mod(u,v);
        if ( w == 0)
            return idiv(v,unit_part(v));
        u = mod(v,w);
        if ( u == 0)
            return idiv(w,unit_part(w));
        v = mod(w,u);
    } while (v != 0);
    return idiv(u,unit_part(u));;
}

// TODO: do we need a variant for unit normal inputs?

} // namespace CGAL

#endif // CGAL_EUCLIDEAN_ALGORITHM_H //
// EOF
