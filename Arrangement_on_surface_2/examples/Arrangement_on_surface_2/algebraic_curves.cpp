#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE
#include <iostream>
int main ()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl; 
  return 0;
}
#else

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>

#include <CGAL/CORE_BigInt.h>

typedef CORE::BigInt Integer;
typedef CGAL::Arr_algebraic_segment_traits_2<Integer> Arr_traits_2;
typedef CGAL::Arrangement_2<Arr_traits_2> Arrangement_2;
typedef Arr_traits_2::Curve_2 Curve_2;
typedef Arr_traits_2::Polynomial_2 Polynomial_2;

int main() {

    // For nice printouts
    CGAL::set_pretty_mode(std::cout);

    Arr_traits_2 arr_traits;

    // Traits class for Polynomial type
    typedef CGAL::Polynomial_traits_d<Polynomial_2> Polynomial_traits_2;

    // Functor to create Polynomial_2 objects
    Polynomial_traits_2::Construct_polynomial construct_polynomial;

    // Functor to create a curve from a Polynomial_2
    Arr_traits_2::Construct_curve_2 construct_curve
        = arr_traits.construct_curve_2_object();
  
    Arrangement_2 arr(&arr_traits);

    // Construct an (unbounded line) with equation 3x-5y+2=0
    std::vector<std::pair<CGAL::Exponent_vector,Integer> > coeffs1;
    coeffs1.push_back(std::make_pair(CGAL::Exponent_vector(1,0),3));
    coeffs1.push_back(std::make_pair(CGAL::Exponent_vector(0,1),-5));
    coeffs1.push_back(std::make_pair(CGAL::Exponent_vector(0,0),2));
    Polynomial_2 f1 = construct_polynomial(coeffs1.begin(),coeffs1.end());

    Curve_2 cv1 = construct_curve(f1);

    std::cout << "Adding curve " << f1 << " to the arrangement" << std::endl;
    CGAL::insert(arr,cv1);

    // Construct the ellipse x^2+3*y^2-10=0
    std::vector<std::pair<CGAL::Exponent_vector,Integer> > coeffs2;
    coeffs2.push_back(std::make_pair(CGAL::Exponent_vector(2,0),1));
    coeffs2.push_back(std::make_pair(CGAL::Exponent_vector(0,2),3));
    coeffs2.push_back(std::make_pair(CGAL::Exponent_vector(0,0),-10));
    Polynomial_2 f2 = construct_polynomial(coeffs2.begin(),coeffs2.end());

    Curve_2 cv2 = construct_curve(f2);

    std::cout << "Adding curve " << f2 << " to the arrangement" << std::endl;
    CGAL::insert(arr,cv2);

    // Construct a cubic curve with isoated point, and vertical asymptote
    // x^2+y^2+xy^2
    std::vector<std::pair<CGAL::Exponent_vector,Integer> > coeffs3;
    coeffs3.push_back(std::make_pair(CGAL::Exponent_vector(2,0),1));
    coeffs3.push_back(std::make_pair(CGAL::Exponent_vector(0,2),1));
    coeffs3.push_back(std::make_pair(CGAL::Exponent_vector(1,2),1));
    Polynomial_2 f3 = construct_polynomial(coeffs3.begin(),coeffs3.end());

    Curve_2 cv3 = construct_curve(f3);

    std::cout << "Adding curve " << f3 << " to the arrangement" << std::endl;
    CGAL::insert(arr,cv3);

    // Construct a curve of degree 6 with equation x^6+y^6-x^3y^3-12
    std::vector<std::pair<CGAL::Exponent_vector,Integer> > coeffs4;
    coeffs4.push_back(std::make_pair(CGAL::Exponent_vector(6,0),1));
    coeffs4.push_back(std::make_pair(CGAL::Exponent_vector(0,6),1));
    coeffs4.push_back(std::make_pair(CGAL::Exponent_vector(3,3),-1));
    coeffs4.push_back(std::make_pair(CGAL::Exponent_vector(0,0),-12));
    Polynomial_2 f4 = construct_polynomial(coeffs4.begin(),coeffs4.end());

    Curve_2 cv4 = construct_curve(f4);

    std::cout << "Adding curve " << f4 << " to the arrangement" << std::endl;
    CGAL::insert(arr,cv4);

    // Print the arrangement size.
    std::cout << "The arrangement size:" << std::endl
              << "   V = " << arr.number_of_vertices()
              << ",  E = " << arr.number_of_edges() 
              << ",  F = " << arr.number_of_faces() << std::endl;
    
    return 0;
}

#endif
