//Author(s) : Michael Hemmer <mhemmer@uni-mainz.de>

/*! \file CGAL/Modular.h
    \brief Defines the class CGAL::Modular and CGAL::Modular_traits.
 
    Provides the \c CGAL::Modular_traits specialization for the build in number 
    types. 
*/


#ifndef CGAL_MODULAR_H
#define CGAL_MODULAR_H 1

#include <CGAL/basic.h>
#include <CGAL/Modular_type.h>
#include <CGAL/Modular_traits.h>
#include <cfloat>

namespace CGAL {

// I/O 

inline std::ostream& operator << (std::ostream& os, const Modular& p) {   
    typedef Modular MOD;
    os <<"("<< p.x()<<"%"<<MOD::get_current_prime()<<")";
    return os;
}


inline std::istream& operator >> (std::istream& is, Modular& p) {
    typedef Modular MOD;
    char ch;
    int prime;

    is >> p.x();
    is >> ch;    // read the %
    is >> prime; // read the prime
    CGAL_precondition(prime==MOD::get_current_prime());
    return is;
}

/*! \brief Specialization of CGAL::NT_traits for \c Modular, which is a model
 * of the \c Field concept. 
 * \ingroup CGAL_NT_traits_spec
 */
template <>
class Algebraic_structure_traits<Modular>
    : public Algebraic_structure_traits_base< Modular ,Field_tag >{};



}///namespace CGAL

#endif //#ifnedef CGAL_MODULAR_H 1
 
