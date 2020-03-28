#include <CGAL/IO/io.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/number_utils.h>
#include <CGAL/int.h>

template< typename NT > NT unit_part(const NT& x);
template< typename NT >
NT unit_part_(const NT& x, CGAL::Field_tag);
template< typename NT >
NT unit_part_(const NT& x, CGAL::Integral_domain_without_division_tag);

template< typename NT >
NT unit_part(const NT& x){
    // the unit part of 0 is defined as 1.
    if (x == 0 ) return NT(1);

    typedef CGAL::Algebraic_structure_traits<NT> AST;
    typedef typename AST::Algebraic_category Algebraic_category;
    return unit_part_(x,Algebraic_category());
}

template< typename NT >
NT unit_part_(const NT& x, CGAL::Integral_domain_without_division_tag){
    // For many other types the only units are just -1 and +1.
    return NT(int(CGAL::sign(x)));
}

template< typename NT >
NT unit_part_(const NT& x, CGAL::Field_tag){
    // For Fields every x != 0 is a unit.
    // Therefore, every x != 0 is its own unit part.
    return x;
}

int main(){
    // Function call for a model of EuclideanRing, i.e. int.
    std::cout<< "int:    unit_part(-3  ): " << unit_part(-3  ) << std::endl;
    // Function call for a model of FieldWithSqrt, i.e. double
    std::cout<< "double: unit_part(-3.0): " << unit_part(-3.0) << std::endl;
    return 0;
}

// Note that this is just an example
// This implementation for unit part won't work for some types, e.g.,
// types that are not RealEmbeddable or types representing structures that have
// more units than just -1 and +1. (e.g. MP_Float representing Z[1/2])
// From there Algebraic_structure_traits provides the functor Unit_part.

