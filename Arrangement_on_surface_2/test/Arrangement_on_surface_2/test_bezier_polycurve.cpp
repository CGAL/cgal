#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE
#include <iostream>
int main ()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl; 
  return 0;
}
#else

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             NT;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef Rat_kernel::Point_2                             Rat_point_2;

//bezier traits
typedef CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits> Bezier_traits_2;
typedef Bezier_traits_2::Curve_2                               Bezier_curve;
typedef Bezier_traits_2::X_monotone_curve_2                    Bezier_x_monotone_curve;

//polyline traits
typedef CGAL::Arr_polyline_traits_2<Bezier_traits_2>            Bezier_polycurve_traits;
typedef Bezier_polycurve_traits::Curve_2                        Polycurve_bezier;
typedef Bezier_polycurve_traits::X_monotone_curve_2             X_polycurve_bezier;


int main (int argc, char *argv[])
{
  Bezier_traits_2                         bezier_traits;
  Bezier_polycurve_traits                 poly_traits;

  //vectors required for polylines 
  std::vector<Rat_point_2>                points_vector;
  std::vector<Bezier_curve>               curves_vector;
  std::vector<Bezier_x_monotone_curve>    x_curves_vector;
  
  //creating a bezier curve
  points_vector.push_back( Rat_point_2(0,0) );
  //points_vector.push_back( Rat_point_2(500,200) );
  points_vector.push_back( Rat_point_2(100,200) );
  points_vector.push_back( Rat_point_2(500,200) );
  points_vector.push_back( Rat_point_2(900,0) );
  Bezier_curve curve_1 (points_vector.begin(), points_vector.end());
  
  std::vector<CGAL::Object> obj_vector;

  //creating x-mono bezier
  bezier_traits.make_x_monotone_2_object()( curve_1, std::back_inserter(obj_vector));
  //std::cout << "number of x_curves: " << obj_vector.size() << std::endl;
  Bezier_x_monotone_curve x_curve_1 = CGAL::object_cast<Bezier_x_monotone_curve>( (obj_vector[0]) );

  //std::cout << x_curve << std::endl;

  points_vector.clear();
  points_vector.push_back( Rat_point_2(900,0) );
  points_vector.push_back( Rat_point_2(1000,400) );
  points_vector.push_back( Rat_point_2(1500,200) );
  points_vector.push_back( Rat_point_2(2000,0) );
  Bezier_curve curve_2 (points_vector.begin(), points_vector.end());
  //creating x-monotne
  obj_vector.clear();
  bezier_traits.make_x_monotone_2_object()( curve_2, std::back_inserter(obj_vector));
  Bezier_x_monotone_curve x_curve_2 = CGAL::object_cast<Bezier_x_monotone_curve>( (obj_vector[0]) );

  //push curves into polyline vectors
  curves_vector.push_back(curve_1);
  curves_vector.push_back(curve_2);
  x_curves_vector.push_back(x_curve_1);
  x_curves_vector.push_back(x_curve_2);
  
  //create polycurves
  Polycurve_bezier polycurve_1 = poly_traits.construct_curve_2_object()(curves_vector.begin(), curves_vector.end());
  X_polycurve_bezier x_polycurve_1 = poly_traits.construct_x_monotone_curve_2_object()(x_curves_vector.begin(), x_curves_vector.end());

  std::cout << "polycurve : " << polycurve_1 << std::endl;
  std::cout << "x_polycurve : " << x_polycurve_1 << std::endl;

  return 0;
}

#endif
