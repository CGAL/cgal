
#include <CGAL/basic.h>
#include <iostream>

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

#define CGAL_ACK_DEBUG_FLAG 0

//#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 1


static const char *ACK_2_ascii_polys[] = {
    
    "P[2(0,P[2(2,1)])(2,P[0(0,-1)])]", // x^2-y^2 

    "P[1(0,P[2(2,1)])(1,P[0(0,1)])]", // x^2+y7 
};

static const int ACK_2_n_polys = 2;

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_curve_kernel_2.h>


int main()
{
  typedef CGAL::CORE_arithmetic_kernel AT;
    typedef AT::Integer Coefficient;
    typedef AT::Rational Rational;
      
    typedef CGAL::internal::Algebraic_real_quadratic_refinement_rep_bfi
        < Coefficient, Rational > Rep_class;
    typedef CGAL::internal::Bitstream_descartes< 
        CGAL::internal::Bitstream_descartes_rndl_tree_traits<
        CGAL::internal::Bitstream_coefficient_kernel<Coefficient > > > 
        Isolator;
    
    typedef CGAL::Algebraic_kernel_d_1<Coefficient,Rational,Rep_class, Isolator> 
        Algebraic_kernel_d_1;

  
    typedef CGAL::Algebraic_curve_kernel_2<Algebraic_kernel_d_1>  AK_2;
 
    typedef  AK_2::Polynomial_2 Poly_2;
    typedef  AK_2::Curve_analysis_2 Curve_analysis_2;
    typedef  Curve_analysis_2::Status_line_1 Status_line_1;

    Poly_2 polys[ACK_2_n_polys];
    
    for(int i = 0; i < ACK_2_n_polys; i++) {
      std::istringstream in(ACK_2_ascii_polys[i]);
        in >> polys[i];    
    }
    
    AK_2 kernel_2;
     
    Curve_analysis_2 c7_c6 =
      kernel_2.construct_curve_2_object()(polys[0]*polys[1]); 
    
    Status_line_1 line = c7_c6.status_line_at_event(0);

    return 0;
}

#else

int main()
{
  std::cerr << "Needs CGAL_HAS_CORE_ARITHMETIC_KERNEL" << std::endl;
  return 0;
}

#endif
