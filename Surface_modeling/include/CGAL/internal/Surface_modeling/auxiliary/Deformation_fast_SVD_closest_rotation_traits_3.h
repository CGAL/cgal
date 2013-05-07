#ifndef CGAL_DEFORMATION_FAST_SVD_CLOSEST_ROTATION_TRAITS_3
#define CGAL_DEFORMATION_FAST_SVD_CLOSEST_ROTATION_TRAITS_3

#include <CGAL/Deformation_Eigen_closest_rotation_traits_3.h>

namespace CGAL {

class Deformation_fast_SVD_closest_rotation_traits_3 :
  public Deformation_Eigen_closest_rotation_traits_3{
public:

  // #define USE_SCALAR_IMPLEMENTATION
  #define USE_SSE_IMPLEMENTATION
  #define COMPUTE_V_AS_MATRIX
  #define COMPUTE_U_AS_MATRIX

  #include "Singular_Value_Decomposition_Preamble.hpp"

  void compute_close_rotation(const Matrix& m, Matrix& R)
  {

    #include "Singular_Value_Decomposition_Kernel_Declarations.hpp"

    ENABLE_SCALAR_IMPLEMENTATION(Sa11.f=m(0,0);) 
    ENABLE_SCALAR_IMPLEMENTATION(Sa21.f=m(1,0);) 
    ENABLE_SCALAR_IMPLEMENTATION(Sa31.f=m(2,0);) 
    ENABLE_SCALAR_IMPLEMENTATION(Sa12.f=m(0,1);) 
    ENABLE_SCALAR_IMPLEMENTATION(Sa22.f=m(1,1);) 
    ENABLE_SCALAR_IMPLEMENTATION(Sa32.f=m(2,1);) 
    ENABLE_SCALAR_IMPLEMENTATION(Sa13.f=m(0,2);) 
    ENABLE_SCALAR_IMPLEMENTATION(Sa23.f=m(1,2);) 
    ENABLE_SCALAR_IMPLEMENTATION(Sa33.f=m(2,2);) 

    ENABLE_SSE_IMPLEMENTATION(Va11=_mm_set1_ps(m(0,0));)
    ENABLE_SSE_IMPLEMENTATION(Va21=_mm_set1_ps(m(1,0));)
    ENABLE_SSE_IMPLEMENTATION(Va31=_mm_set1_ps(m(2,0));)
    ENABLE_SSE_IMPLEMENTATION(Va12=_mm_set1_ps(m(0,1));)
    ENABLE_SSE_IMPLEMENTATION(Va22=_mm_set1_ps(m(1,1));)
    ENABLE_SSE_IMPLEMENTATION(Va32=_mm_set1_ps(m(2,1));)
    ENABLE_SSE_IMPLEMENTATION(Va13=_mm_set1_ps(m(0,2));)
    ENABLE_SSE_IMPLEMENTATION(Va23=_mm_set1_ps(m(1,2));)
    ENABLE_SSE_IMPLEMENTATION(Va33=_mm_set1_ps(m(2,2));)

    #include "Singular_Value_Decomposition_Main_Kernel_Body.hpp"  
    
    Matrix U, V;

#ifdef USE_SCALAR_IMPLEMENTATION
    U(0,0) = Su11.f;
    U(1,0) = Su21.f;
    U(2,0) = Su31.f;
    U(0,1) = Su12.f;
    U(1,1) = Su22.f;
    U(2,1) = Su32.f;
    U(0,2) = Su13.f;
    U(1,2) = Su23.f;
    U(2,2) = Su33.f;

    V(0,0) = Sv11.f;
    V(1,0) = Sv21.f;
    V(2,0) = Sv31.f;
    V(0,1) = Sv12.f;
    V(1,1) = Sv22.f;
    V(2,1) = Sv32.f;
    V(0,2) = Sv13.f;
    V(1,2) = Sv23.f;
    V(2,2) = Sv33.f;
#endif

#ifdef USE_SSE_IMPLEMENTATION
    float buf[4];
    _mm_storeu_ps(buf,Vu11);U(0,0)=buf[0];
    _mm_storeu_ps(buf,Vu21);U(1,0)=buf[0];
    _mm_storeu_ps(buf,Vu31);U(2,0)=buf[0];
    _mm_storeu_ps(buf,Vu12);U(0,1)=buf[0];
    _mm_storeu_ps(buf,Vu22);U(1,1)=buf[0];
    _mm_storeu_ps(buf,Vu32);U(2,1)=buf[0];
    _mm_storeu_ps(buf,Vu13);U(0,2)=buf[0];
    _mm_storeu_ps(buf,Vu23);U(1,2)=buf[0];
    _mm_storeu_ps(buf,Vu33);U(2,2)=buf[0];

    _mm_storeu_ps(buf,Vv11);V(0,0)=buf[0];
    _mm_storeu_ps(buf,Vv21);V(1,0)=buf[0];
    _mm_storeu_ps(buf,Vv31);V(2,0)=buf[0];
    _mm_storeu_ps(buf,Vv12);V(0,1)=buf[0];
    _mm_storeu_ps(buf,Vv22);V(1,1)=buf[0];
    _mm_storeu_ps(buf,Vv32);V(2,1)=buf[0];
    _mm_storeu_ps(buf,Vv13);V(0,2)=buf[0];
    _mm_storeu_ps(buf,Vv23);V(1,2)=buf[0];
    _mm_storeu_ps(buf,Vv33);V(2,2)=buf[0];
#endif

    R = V * U.transpose(); // no need to check determinant (according to the paper their algo automatically fix reflections)
  }
};

}//namespace CGAL
#endif // CGAL_DEFORMATION_FAST_SVD_CLOSEST_ROTATION_TRAITS_3