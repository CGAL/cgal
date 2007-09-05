#ifndef CGAL_ARRANGEMENT_OF_SPHERES_3_BASIC_H
#define CGAL_ARRANGEMENT_OF_SPHERES_3_BASIC_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/CORE_BigRat.h>
//#define CGAL_ARRANGEMENT_OF_SPHERES_3_USE_TEMPLATES





#define CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE CGAL_BEGIN_NAMESPACE namespace CGAL_AOS3_internal {
#define CGAL_AOS3_END_INTERNAL_NAMESPACE CGAL_END_NAMESPACE }

#define CGAL_AOS3_INTERNAL_NS CGAL::CGAL_AOS3_internal

#include <CGAL/Arrangement_of_spheres_3/Coordinate_index.h>

/*template <class NT>
typename CGAL::Root_of_traits<NT>::RootOf_2 translate_quadratic(const CGAL::Sqrt_extension<NT,NT> &q) {
  typedef typename CGAL::Root_of_traits<NT>::RootOf_2 R;
  
  make_root_of_2(q.a1()/2,	
		 -q.a0()*q.a1(),	
		 (q.a0()*q.a0() * q.a1()	
		  -q.root()*q.a1()*q.a1()*q.a1()),
		 false);
		 }*/
#if 1
template <class NT>
typename CGAL::Sqrt_extension<NT, NT> construct_quadratic(NT a, NT b, NT c, bool smaller) {
  NT root=b*b-4*a*c;
  
  typename CGAL::Sqrt_extension<NT, NT> ret;
  if (root != 0) {
    ret=   CGAL::Sqrt_extension<NT, NT>(-b/(2*a),
						(smaller? -1: 1) /(2*a),
						b*b-4*a*c);
  } else {
    ret= -b/(2*a);
  }
  //typename CGAL::Root_of_traits<NT>::RootOf_2 check= CGAL::make_root_of_2(a,b,c, smaller);
  /*std::cout << CGAL::Interval_nt<true>(CGAL::to_interval(ret)) << "==" 
    << CGAL::Interval_nt<true>(CGAL::to_interval(check)) << std::endl;*/
  return ret;
}


#define CGAL_AOS3_DEFINE_QUADRATIC typedef typename CGAL::Sqrt_extension<NT, NT> Quadratic_NT
#define CGAL_AOS3_CONSTRUCT_QUADRATIC(a,b,c,smaller) construct_quadratic(a,b,c,smaller)

//a = (1/2)*y, c = (1/2)*x^2*y-(1/2)*y^3*z, b = -x*y

#else
#define CGAL_AOS3_DEFINE_QUADRATIC  typedef typename CGAL::Root_of_traits<NT>::RootOf_2 Quadratic_NT
#define CGAL_AOS3_CONSTRUCT_QUADRATIC(a,b,c,smaller) make_root_of_2(a,b,c,smaller)
#endif
// from x+ y sqrt(z) = (-b +-sqrt(b^2-4ac))/2a
// take a=y/2, b=-xy, z=b^2-4ac=x^2y^2-2yc c=(-z+x^2 y^2)/2y



//




#ifndef CGAL_AOS3_USE_TEMPLATES
#define CGAL_AOS3_TYPENAME

#define CGAL_AOS3_TEMPLATE
#define CGAL_AOS3_TARG
#define CGAL_AOS3_TRAITS public: typedef Arrangement_of_spheres_traits_3 Traits; private: friend class Semicolon_eater
 
CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE
//
struct Arrangement_of_spheres_3_geom_traits: Cartesian<Gmpq>{};
CGAL_AOS3_END_INTERNAL_NAMESPACE
#include <CGAL/Arrangement_of_spheres_3/Rule_direction.h>
#include <CGAL/Arrangement_of_spheres_traits_3.h>





#else
#define CGAL_AOS3_TYPENAME typename
#define CGAL_AOS3_TEMPLATE template <class Traits_t>
#define CGAL_AOS3_TARG <Traits_t>
#define CGAL_AOS3_TRAITS public: typedef Traits_t Traits; private: friend class Semicolon_eater


#endif



#endif
