// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de  
//                 Ralf Schindlbeck <rschindl@mpi-inf.mpg.de>
//                 Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

/*! \file Rotated_algebraic_curve_kernel_2.h
 *  \brief refines Algebraic_curve_kernel_2 to support fixed angle degree
 *  rotations
 */

#ifndef CGAL_ROTATED_ALGEBRAIC_CURVE_KERNEL_2_H
#define CGAL_ROTATED_ALGEBRAIC_CURVE_KERNEL_2_H

#include <CGAL/basic.h>
#include <CGAL/Algebraic_curve_kernel_2.h>

#include <CGAL/Algebraic_curve_kernel_2/trigonometric_approximation.h>

CGAL_BEGIN_NAMESPACE

namespace { 

// quick hack to rebind curve pairs & ak_1
template <class Coefficient, class Rational>
struct Rebind_helper {
    
    typedef CGAL::CGALi::Algebraic_real_quadratic_refinement_rep_bfi
        < Coefficient, Rational > Rep_class;
    typedef CGAL::CGALi::Bitstream_descartes< CGAL::Polynomial< Coefficient >,
        Rational > Isolator;
    
    typedef CGAL::Algebraic_kernel_1<Coefficient, Rational, Rep_class,
         Isolator>  Kernel_1;
   
};

#ifdef CnX_USE_whatever_blabla

// put some CnX-specific code

#endif     

// a bunch of rationals for convenience
template <class Rational>
struct _Packed {

    _Packed() : a1_sine(0), a2_sine(0), b1_sine(0), b2_sine(0),
        a1_cosine(0), a2_cosine(0), b1_cosine(0), b2_cosine(0),
        a3_sine(0), a3_cosine(0) {
    }  

    Rational a1_sine, a2_sine, b1_sine, b2_sine;
    Rational a1_cosine, a2_cosine, b1_cosine, b2_cosine;
    Rational a3_sine, a3_cosine; 
};

// a lengthy switch-case to select appropriate angle coefficients
template <class Rational>
void _15_degree_selector(_Packed<Rational>& t, int angle) {
               
    switch(angle) {
        case 0:
            //sin(0)=0
            //, also equals sin(180)=(-1)*0
            //cos(0)=1
            //, also equals cos(180)=(-1)*1
            t.a1_cosine = Rational(1);
            break;
        case 15:
            //sin(15)=-1/4*sqrt(2)+1/4*sqrt(6)
            //, also equals sin(195)=(-1)*(-1/4*sqrt(2)+1/4*sqrt(6))
            //cos(15)=1/4*sqrt(2)+1/4*sqrt(6)
            //, also equals cos(195)=(-1)*(1/4*sqrt(2)+1/4*sqrt(6))
            t.b1_sine = Rational(-1,4); t.b2_sine = Rational(1,4);
            t.b1_cosine = Rational(1,4); t.b2_cosine = Rational(1,4);
            break;
        case 30:
            //sin(30)=1/2
            //, also equals sin(210)=(-1)*(1/2)
            //cos(30)=1/2*sqrt(3)
            //, also equals cos(210)=(-1)*(1/2*sqrt(3))
            t.a1_sine = Rational(1,2);
            t.a2_cosine = Rational(1,2);
            break;
        case 45:
            //sin(45)=1/2*srt(2)
            //, also equals sin(225)=(-1)*(1/2*srt(2))
            //cos(45)=1/2*sqrt(2)
            //, also equals cos(225)=(-1)*(1/2*sqrt(2))
            t.b1_sine = Rational(1,2);
            t.b1_cosine = Rational(1,2);
            break;
        case 60:
            //sin(60)=1/2*sqrt(3)
            //, also equals sin(240)=(-1)*(1/2*sqrt(3))
            //cos(60)=1/2
            //, also equals cos(240)=(-1)*(1/2)
            t.a2_sine = Rational(1,2);
            t.a1_cosine = Rational(1,2);
            break;
        case 75:
            //sin(75)=1/4*sqrt(2)+1/4*sqrt(6)
            //, also equals sin(255)=(-1)*(1/4*sqrt(2)+1/4*sqrt(6))
            //cos(75)=-1/4*sqrt(2)+1/4*sqrt(6)
            //, also equals cos(255)=(-1)*(-1/4*sqrt(2)+1/4*sqrt(6))
            t.b1_sine = Rational(1,4); t.b2_sine = Rational(1,4);
            t.b1_cosine = Rational(-1,4); t.b2_cosine = Rational(1,4);
            break;
        case 90:
            //sin(90)=1
            //, also equals sin(270)=(-1)*1
            //cos(90)=0
            //, also equals cos(270)=(-1)*0
            t.a1_sine = Rational(1);
            break;
        case 105:
            //sin(105)=1/4*sqrt(2)+1/4*sqrt(6)
            //, also equals sin(285)=(-1)*(1/4*sqrt(2)+1/4*sqrt(6))
            //cos(105)=1/4*sqrt(2)-1/4*sqrt(6)
            //, also equals cos(285)=(-1)*(1/4*sqrt(2)-1/4*sqrt(6))
            t.b1_sine = Rational(1,4); t.b2_sine = Rational(1,4);
            t.b1_cosine = Rational(1,4); t.b2_cosine = Rational(-1,4);
            break;
        case 120:
            //sin(120)=1/2*sqrt(3)
            //, also equals sin(300)=(-1)*(1/2*sqrt(3))
            //cos(120)=-1/2
            //, also equals cos(300)=(-1)*(-1/2)
            t.a2_sine = Rational(1,2);
            t.a1_cosine = Rational(-1,2);
            break;
        case 135:
            //sin(135)=1/2*srt(2)
            //, also equals sin(315)=(-1)*(1/2*srt(2))
            //cos(135)=-1/2*sqrt(2)
            //, also equals cos(315)=(-1)*(-1/2*sqrt(2))
            t.b1_sine = Rational(1,2);
            t.b1_cosine = Rational(-1,2);
            break;
        case 150:
            //sin(150)=1/2
            //, also equals sin(330)=(-1)*1/2
            //cos(150)=-1/2*sqrt(3)
            //, also equals cos(330)=(-1)*(-1/2*sqrt(3))
            t.a1_sine = Rational(1,2);
            t.a2_cosine = Rational(-1,2);
            break;
        case 165:
            //sin(165)=-1/4*sqrt(2)+1/4*sqrt(6)
            //, also equals sin(345)=(-1)*(-1/4*sqrt(2)+1/4*sqrt(6))
            //cos(165)=-1/4*sqrt(2)-1/4*sqrt(6)
            //, also equals cos(345)=(-1)*(-1/4*sqrt(2)-1/4*sqrt(6))
            t.b1_sine = Rational(-1,4); t.b2_sine = Rational(1,4);
            t.b1_cosine = Rational(-1,4); t.b2_cosine = Rational(-1,4);
            break;
        default:
            CGAL_error_msg("The angle must be multiple of 15 degrees !!\n");
    }
}

// the same for 18 degree
template <class Rational>
void _18_degree_selector(_Packed<Rational>& t, int angle) {

    switch(angle) {
        case 0:
            //sin(0)=0
            //, also equals sin(180)=(-1)*0
            //cos(0)=1
            //, also equals cos(180)=(-1)*1
            t.a1_cosine = Rational(1);
            break;
        case 18:
            //sin(18)=-1/4+1/4*sqrt(5)
            //, also equals sin(198)=(-1)*(-1/4+1/4*sqrt(5))
            //cos(18)=1/4*sqrt(10+2*sqrt(5))
            //, also equals cos(198)=(-1)*(1/4*sqrt(10+2*sqrt(5)))
            t.a1_sine = Rational(-1,4), t.b1_sine = Rational(1,4);
            t.a2_cosine = Rational(1,4);
            break;
        case 36:
            //sin(36)=1/4*sqrt(10-2*sqrt(5))
            //, also equals sin(216)=(-1)*(1/4*sqrt(10-2*sqrt(5)))
            //cos(36)=1/4+1/4*sqrt(5)
            //, also equals cos(216)=(-1)*(1/4+1/4*sqrt(5))
            t.a2_sine = Rational(-1,8), t.b2_sine = Rational(1,8);
            t.a1_cosine = Rational(1,4), t.b1_cosine = Rational(1,4);
            break;
        case 54:
            //sin(54)=1/4+1/4*sqrt(5)
            //, also equals sin(234)=(-1)*(1/4+1/4*sqrt(5))
            //cos(54)=1/4*sqrt(10-2*sqrt(5))
            //, also equals cos(234)=(-1)*(1/4*sqrt(10-2*sqrt(5)))
            t.a1_sine = Rational(1,4), t.b1_sine = Rational(1,4);
            t.a2_cosine = Rational(-1,8), t.b2_cosine = Rational(1,8);
            break;
        case 72:
            //sin(72)=1/4*sqrt(10+2*sqrt(5))
            //, also equals sin(252)=(-1)*(1/4*sqrt(10+2*sqrt(5)))
            //cos(72)=-1/4+1/4*sqrt(5)
            //, also equals cos(252)=(-1)*(-1/4+1/4*sqrt(5))
            t.a2_sine = Rational(1,4);
            t.a1_cosine = Rational(-1,4), t.b1_cosine = Rational(1,4);
            break;
        case 90:
            //sin(90)=1
            //, also equals sin(270)=(-1)*1
            //cos(90)=0
            //, also equals cos(270)=(-1)*0
            t.a1_sine = Rational(1);
            break;
        case 108:
            //sin(108)=1/4*sqrt(10+2*sqrt(5))
            //, also equals sin(288)=(-1)*(1/4*sqrt(10+2*sqrt(5)))
            //cos(108)=1/4-1/4*sqrt(5)
            //, also equals cos(288)=(-1)*(1/4-1/4*sqrt(5))
            t.a2_sine = Rational(1,4);
            t.a1_cosine = Rational(1,4), t.b1_cosine = Rational(-1,4);
            break;
        case 126:
            //sin(126)=1/4+1/4*sqrt(5)
            //, also equals sin(306)=(-1)*(1/4+1/4*sqrt(5))
            //cos(126)=1/8*sqrt(10+2*sqrt(5))-1/8*sqrt(10+2*sqrt(5))*sqrt(5)
            //, also equals cos(306)=(-1)*(-1/4*sqrt(10-2*sqrt(5)))
            t.a1_sine = Rational(1,4); t.b1_sine = Rational(1,4);
            t.a2_cosine = Rational(1,8), t.b2_cosine=Rational(-1,8);
            break;
        case 144:
            //sin(144)=-1/8*sqrt(10+2*sqrt(5))+1/8*sqrt(10+2*sqrt(5))*sqrt(5)
            //, also equals sin(324)=(-1)*(1/4*sqrt(10-2*sqrt(5)))
            //cos(144)=-1/4-1/4*sqrt(5)
            //, also equals cos(324)=(-1)*(-1/4-1/4*sqrt(5))
            t.a2_sine = Rational(-1,8), t.b2_sine = Rational(1,8);
            t.a1_cosine = Rational(-1,4), t.b1_cosine = Rational(-1,4);
            break;
        case 162:
            //sin(162)=-1/4+1/4*sqrt(5)
            //, also equals sin(342)=(-1)*(-1/4+1/4*sqrt(5))
            //cos(162)=-1/4*sqrt(10+2*sqrt(5))
            //, also equals cos(342)=(-1)*(-1/4*sqrt(10+2*sqrt(5)))
            t.a1_sine = Rational(-1,4), t.b1_sine = Rational(1,4);
            t.a2_cosine = Rational(-1,4);
            break;

        default:
            CGAL_error_msg("The angle must be multiple of 18 degrees !!\n");
    }
}

/*!\brief 
 * maps \c BaseAngle to the one of supported base angles or -1 otherwise 
 * (to prevent spurious template instantiations)
 */
template <int BaseAngle>
struct Normalized_angle {
    enum { 
        // wrap around the base angle to [0..360) range 
        tmp = BaseAngle % 360,
        modulo = tmp + (tmp < 0 ? 360 : 0),
        angle = ((modulo % 45) == 0 ? 45 :
                 (modulo % 30) == 0 ? 30 :
                 (modulo % 15) == 0 ? 15 : (modulo % 18) == 0 ? 18 : 
                 (modulo % 3) == 0 ? 3 : -1)
    };    
};

template <class ArithmeticKernel, int BaseAngle>
struct Rotation_traits_base {

    typedef void Extended_coefficient;
};

/*!\brief
 * Rotation kernel for angles which are multiple of 3 degrees
 */
template <class AlgebraicCurveKernel_2>
struct Rotation_traits_base<AlgebraicCurveKernel_2, 3> {

    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    typedef typename CGAL::Get_arithmetic_kernel<
        typename Algebraic_curve_kernel_2::Boundary>::Arithmetic_kernel
            Arithmetic_kernel;    

    typedef typename Arithmetic_kernel::Integer Integer;
    typedef typename Arithmetic_kernel::Rational Rational;

private:
    typedef CGAL::Sqrt_extension<Rational, Integer> EXT2; //root(2)
    typedef CGAL::Sqrt_extension<Integer, Integer> EXT2_int; //root(5)

    typedef CGAL::Sqrt_extension<EXT2, Integer> EXT3; //root(3)
    typedef CGAL::Sqrt_extension<EXT3, Integer> EXT4; //root(5)
    typedef CGAL::Sqrt_extension<EXT2_int, Integer> EXT3_int; //root(5)
    typedef CGAL::Sqrt_extension<EXT3_int, Integer> EXT4_int; //root(5)
    typedef CGAL::Sqrt_extension<EXT4, EXT4_int> EXT5; //root(5+root(5))

public:
    typedef EXT5 Extended_rational;

    typedef typename CGAL::Fraction_traits<Extended_rational>::Numerator_type
        Extended_coefficient;

    //! bivaritate polynomial over integers 
    typedef typename Algebraic_curve_kernel_2::Polynomial_2 Poly_int_2;
    //! univariate polynomial over sqrts
    typedef CGAL::Polynomial<Extended_coefficient> Poly_ext_1;
    //! bivaritate polynomial over sqrts
    typedef CGAL::Polynomial<Poly_ext_1> Poly_ext_2;

    //!\name functor invokation 
    //!@{

    void operator()(Extended_rational& esine, Extended_rational& ecosine,
         Extended_rational& ezero, int& angle) {
        
        if(angle % 3) 
            CGAL_error_msg("angle is not a multiple of 3");

        int angle_help = angle, sign = 1;
        
        // indicates that addition theorems are necessary
        bool use_addition = false;
        bool new_start_point = false;
        bool angle_360_deg = false;

        // the nearest rotation-angle which is multiple of 15 degrees
        if(angle % 15 != 0) {
            angle_help = (angle / 15) * 15;
            
            if((angle - angle_help) > 7) {
                angle_help += 15;
                new_start_point = true;
            }

            if(angle_help == 360) {
                angle_help = 0;
                angle_360_deg = true;
            }
            use_addition = true;
        }

        if(angle_help >= 180) {
            sign = -1;
            angle_help -= 180;
        }
        
        _Packed<Rational> t;
        _15_degree_selector<Rational>(t, angle_help);

        Extended_rational esine3, ecosine3;

        esine3 =
            EXT5(EXT4(EXT3(EXT2(Rational(0),Rational(-1,16),Integer(2)),
                           EXT2(Rational(0),Rational(-1,16),Integer(2)),
                           Integer(3)),
                      EXT3(EXT2(Rational(0),Rational(1,16),Integer(2)),
                           EXT2(Rational(0),Rational(1,16),Integer(2)),
                           Integer(3)),
                      Integer(5)),
                 EXT4(EXT3(EXT2(Rational(1,8),Rational(0),Integer(2)),
                           EXT2(Rational(-1,8),Rational(0),Integer(2)),
                           Integer(3)),
                      EXT3(EXT2(Rational(0),Rational(0),Integer(2)),
                           EXT2(Rational(0),Rational(0),Integer(2)),
                           Integer(3)),
                      Integer(5)),
                 EXT4_int(EXT3_int(EXT2_int(Integer(5),Integer(0),Integer(2)),
                           EXT2_int(Integer(0),Integer(0),Integer(2)),
                           Integer(3)),
                      EXT3_int(EXT2_int(Integer(1),Integer(0),Integer(2)),
                           EXT2_int(Integer(0),Integer(0),Integer(2)),
                           Integer(3)),
                      Integer(5)));
        ecosine3 =
            EXT5(EXT4(EXT3(EXT2(Rational(0),Rational(1,16),Integer(2)),
                           EXT2(Rational(0),Rational(-1,16),Integer(2)),
                           Integer(3)),
                      EXT3(EXT2(Rational(0),Rational(-1,16),Integer(2)),
                           EXT2(Rational(0),Rational(1,16),Integer(2)),
                           Integer(3)),
                      Integer(5)),
                 EXT4(EXT3(EXT2(Rational(1,8),Rational(0),Integer(2)),
                           EXT2(Rational(1,8),Rational(0),Integer(2)),
                           Integer(3)),
                      EXT3(EXT2(Rational(0),Rational(0),Integer(2)),
                           EXT2(Rational(0),Rational(0),Integer(2)),
                           Integer(3)),
                      Integer(5)),
                 EXT4_int(EXT3_int(EXT2_int(Integer(5),Integer(0),Integer(2)),
                           EXT2_int(Integer(0),Integer(0),Integer(2)),
                           Integer(3)),
                      EXT3_int(EXT2_int(Integer(1),Integer(0),Integer(2)),
                           EXT2_int(Integer(0),Integer(0),Integer(2)),
                           Integer(3)),
                      Integer(5)));

         esine = EXT5(EXT4(EXT3(EXT2(t.a1_sine,t.b1_sine,Integer(2)),
                               EXT2(t.a2_sine,t.b2_sine,Integer(2)),
                               Integer(3)),
                          EXT3(EXT2(Rational(0),Rational(0),Integer(2)),
                               EXT2(Rational(0),Rational(0),Integer(2)),
                               Integer(3)),
                          Integer(5)),
                     EXT4(EXT3(EXT2(Rational(0),Rational(0),Integer(2)),
                               EXT2(Rational(0),Rational(0),Integer(2)),
                               Integer(3)),
                          EXT3(EXT2(Rational(0),Rational(0),Integer(2)),
                               EXT2(Rational(0),Rational(0),Integer(2)),
                               Integer(3)),
                          Integer(5)),
                     EXT4_int(EXT3_int(EXT2_int(Integer(5),Integer(0),
                               Integer(2)),
                               EXT2_int(Integer(0),Integer(0),Integer(2)),
                               Integer(3)),
                          EXT3_int(EXT2_int(Integer(1),Integer(0),Integer(2)),
                               EXT2_int(Integer(0),Integer(0),Integer(2)),
                               Integer(3)),
                          Integer(5)));
        ecosine = EXT5(EXT4(EXT3(EXT2(t.a1_cosine,t.b1_cosine,Integer(2)),
                               EXT2(t.a2_cosine,t.b2_cosine,Integer(2)),
                               Integer(3)),
                          EXT3(EXT2(Rational(0),Rational(0),Integer(2)),
                               EXT2(Rational(0),Rational(0),Integer(2)),
                               Integer(3)),
                          Integer(5)),
                     EXT4(EXT3(EXT2(Rational(0),Rational(0),Integer(2)),
                               EXT2(Rational(0),Rational(0),Integer(2)),
                               Integer(3)),
                          EXT3(EXT2(Rational(0),Rational(0),Integer(2)),
                               EXT2(Rational(0),Rational(0),Integer(2)),
                               Integer(3)),
                          Integer(5)),
                     EXT4_int(EXT3_int(EXT2_int(Integer(5),Integer(0),
                          Integer(2)),
                               EXT2_int(Integer(0),Integer(0),Integer(2)),
                               Integer(3)),
                          EXT3_int(EXT2_int(Integer(1),Integer(0),Integer(2)),
                               EXT2_int(Integer(0),Integer(0),Integer(2)),
                               Integer(3)),
                          Integer(5)));
        
        if(sign == -1) {
            esine = -esine;
            ecosine = -ecosine;
        }
    
        if(use_addition)
        {
            if(sign == -1)
                angle_help += 180;
            
            if(angle_360_deg)
                angle_help = 360;
            
            int inc = 3;
            if(new_start_point) {
                esine3 = -esine3;
                inc = -3;
            }
            while((angle - angle_help) != 0) {
                //sin(x+y) = sin(x)*cos(y) + sin(y)*cos(x)
                Extended_rational tmp = esine*ecosine3 + esine3*ecosine;
                //cos(x+y) = cos(x)*cos(y) - sin(x)*sin(y)
                ecosine = ecosine*ecosine3 - esine*esine3;
                esine = tmp;                    
                angle_help += inc;
            }
        }

        //defines the zero, needed for the subs_poly_* creation
        ezero =
            EXT5(EXT4(EXT3(EXT2(Rational(0),Rational(0),Integer(2)),
                           EXT2(Rational(0),Rational(0),Integer(2)),
                           Integer(3)),
                      EXT3(EXT2(Rational(0),Rational(0),Integer(2)),
                           EXT2(Rational(0),Rational(0),Integer(2)),
                           Integer(3)),
                      Integer(5)),
                 EXT4(EXT3(EXT2(Rational(0),Rational(0),Integer(2)),
                           EXT2(Rational(0),Rational(0),Integer(2)),
                           Integer(3)),
                      EXT3(EXT2(Rational(0),Rational(0),Integer(2)),
                           EXT2(Rational(0),Rational(0),Integer(2)),
                           Integer(3)),
                      Integer(5)),
                 EXT4_int(EXT3_int(EXT2_int(Integer(5),Integer(0),Integer(2)),
                           EXT2_int(Integer(0),Integer(0),Integer(2)),
                           Integer(3)),
                      EXT3_int(EXT2_int(Integer(1),Integer(0),Integer(2)),
                           EXT2_int(Integer(0),Integer(0),Integer(2)),
                           Integer(3)),
                      Integer(5)));
    }
    //!@}
};

/*!\brief
 * Rotation kernel for angles which are multiple of 15 degrees
 */
template <class AlgebraicCurveKernel_2>
struct Rotation_traits_base<AlgebraicCurveKernel_2, 15> {

    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    typedef typename CGAL::Get_arithmetic_kernel<
        typename Algebraic_curve_kernel_2::Boundary>::Arithmetic_kernel
            Arithmetic_kernel;    

    typedef typename Arithmetic_kernel::Integer Integer;
    typedef typename Arithmetic_kernel::Rational Rational;

private:
    typedef CGAL::Sqrt_extension<Rational, Integer> EXT2; //root(2)
    typedef CGAL::Sqrt_extension<Integer, Integer> EXT2_int; //root(5)
    typedef CGAL::Sqrt_extension<EXT2, Integer> EXT3; //root(3)

public:
    typedef EXT3 Extended_rational;

    typedef typename CGAL::Fraction_traits<Extended_rational>::Numerator_type
        Extended_coefficient;

    //! bivaritate polynomial over integers 
    typedef typename Algebraic_curve_kernel_2::Polynomial_2 Poly_int_2;
    //! univariate polynomial over sqrts
    typedef CGAL::Polynomial<Extended_coefficient> Poly_ext_1;
    //! bivaritate polynomial over sqrts
    typedef CGAL::Polynomial<Poly_ext_1> Poly_ext_2;

    //!\name functor invokation
    //!@{    

    void operator()(Extended_rational& esine, Extended_rational& ecosine,
         Extended_rational& ezero, int& angle) {
        
        if(angle % 15) 
            CGAL_error_msg("angle is not a multiple of 15");

        int angle_help = angle, sign = 1;
        if(angle_help >= 180) {
            sign = -1;
            angle_help -= 180;
        }

        _Packed<Rational> t;
        _15_degree_selector<Rational>(t, angle_help);

        //now, coefficients are set in the esine and ecosine
        //if rot_angle >= 180 => multiply by -1 (-> reflection)
        esine = EXT3(EXT2(t.a1_sine,t.b1_sine,Integer(2)),
                     EXT2(t.a2_sine,t.b2_sine,Integer(2)),
                     Integer(3));
        ecosine = EXT3(EXT2(t.a1_cosine,t.b1_cosine,Integer(2)),
                     EXT2(t.a2_cosine,t.b2_cosine,Integer(2)),
                     Integer(3));

        if(sign == -1) {
            esine = -esine;
            ecosine = -ecosine;
        }

        ezero =
            EXT3(EXT2(Rational(0),Rational(0),Integer(2)),
                 EXT2(Rational(0),Rational(0),Integer(2)),
                 Integer(3));
    }
    //!@}
};

/*!\brief
 * Rotation kernel for angles which are multiple of 18 degrees
 */
template <class AlgebraicCurveKernel_2>
struct Rotation_traits_base<AlgebraicCurveKernel_2, 18> {

    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    typedef typename CGAL::Get_arithmetic_kernel<
        typename Algebraic_curve_kernel_2::Boundary>::Arithmetic_kernel
            Arithmetic_kernel;    

    typedef typename Arithmetic_kernel::Integer Integer;
    typedef typename Arithmetic_kernel::Rational Rational;

private:
    typedef CGAL::Sqrt_extension<Rational, Integer> EXT2; //root(2)
    typedef CGAL::Sqrt_extension<Integer, Integer> EXT2_int; //root(5)

    typedef CGAL::Sqrt_extension<EXT2_int, EXT2_int> EXT3i_int; 
    typedef CGAL::Sqrt_extension<EXT2, EXT2_int> EXT3i; //root(10+2*root(5))
    
public:
    typedef EXT3i Extended_rational;

    typedef typename CGAL::Fraction_traits<Extended_rational>::Numerator_type
        Extended_coefficient;

    //! bivaritate polynomial over integers 
    typedef typename Algebraic_curve_kernel_2::Polynomial_2 Poly_int_2;
    //! univariate polynomial over sqrts
    typedef CGAL::Polynomial<Extended_coefficient> Poly_ext_1;
    //! bivaritate polynomial over sqrts
    typedef CGAL::Polynomial<Poly_ext_1> Poly_ext_2;

    //!\name functor invokation
    //!@{    

    // computes rotation angles
    void operator()(Extended_rational& esine, Extended_rational& ecosine,
         Extended_rational& ezero, int& angle) {
        
        if(angle % 18) 
            CGAL_error_msg("angle is not a multiple of 18");

        int angle_help = angle, sign = 1;
        if(angle_help >= 180) {
            sign = -1;
            angle_help -= 180;
        }

        _Packed<Rational> t;
        _18_degree_selector<Rational>(t, angle_help);

        esine = EXT3i(EXT2(t.a1_sine, t.b1_sine, Integer(5)),
                      EXT2(t.a2_sine, t.b2_sine, Integer(5)),
                      EXT2_int(Integer(10),Integer(2),Integer(5)));
                    
        ecosine = EXT3i(EXT2(t.a1_cosine, t.b1_cosine, Integer(5)),
                        EXT2(t.a2_cosine, t.b2_cosine, Integer(5)),
                        EXT2_int(Integer(10),Integer(2),Integer(5)));

        ezero = EXT3i(EXT2(Rational(0),Rational(0),Integer(5)),
                      EXT2(Rational(0),Rational(0),Integer(5)),
                      EXT2_int(Integer(10),Integer(2),Integer(5)));

        if(sign == -1) {
            esine = -esine;
            ecosine = -ecosine;
        }
    }
    //!@}
};

/*!\brief
 * Rotation kernel for angles which are multiple of 30 degrees
 */
template <class AlgebraicCurveKernel_2>
struct Rotation_traits_base<AlgebraicCurveKernel_2, 30> {

    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    typedef typename CGAL::Get_arithmetic_kernel<
        typename Algebraic_curve_kernel_2::Boundary>::Arithmetic_kernel
            Arithmetic_kernel;    

    typedef typename Arithmetic_kernel::Integer Integer;
    typedef typename Arithmetic_kernel::Rational Rational;

private:
    typedef CGAL::Sqrt_extension<Rational, Integer> EXT2; //root(3)
    typedef CGAL::Sqrt_extension<Integer, Integer> EXT2_int; 

public:
    typedef EXT2 Extended_rational;

    typedef typename CGAL::Fraction_traits<Extended_rational>::Numerator_type
        Extended_coefficient;

    //! bivaritate polynomial over integers 
    typedef typename Algebraic_curve_kernel_2::Polynomial_2 Poly_int_2;
    //! univariate polynomial over sqrts
    typedef CGAL::Polynomial<Extended_coefficient> Poly_ext_1;
    //! bivaritate polynomial over sqrts
    typedef CGAL::Polynomial<Poly_ext_1> Poly_ext_2;

    //!\name functor invokation
    //!@{    

    // computes rotation angles
    void operator()(Extended_rational& esine, Extended_rational& ecosine,
         Extended_rational& ezero, int& angle) {
        
        if(angle % 30) 
            CGAL_error_msg("angle is not a multiple of 30");

        int angle_help = angle, sign = 1;
        if(angle_help >= 180) {
            sign = -1;
            angle_help -= 180;
        }

        _Packed<Rational> t;
        _15_degree_selector<Rational>(t, angle_help);

        esine = EXT2(t.a1_sine,t.a2_sine,Integer(3)); 

        ecosine = EXT2(t.a1_cosine,t.a2_cosine, Integer(3));

        ezero = EXT2(Rational(0),Rational(0),Integer(3));

        if(sign == -1) {
            esine = -esine;
            ecosine = -ecosine;
        }
    }
    //!@}
};

/*!\brief
 * Rotation kernel for angles which are multiple of 45 degrees
 */
template <class AlgebraicCurveKernel_2>
struct Rotation_traits_base<AlgebraicCurveKernel_2, 45> {

    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    typedef typename CGAL::Get_arithmetic_kernel<
        typename Algebraic_curve_kernel_2::Boundary>::Arithmetic_kernel
            Arithmetic_kernel;    

    typedef typename Arithmetic_kernel::Integer Integer;
    typedef typename Arithmetic_kernel::Rational Rational;

private:
    typedef CGAL::Sqrt_extension<Rational, Integer> EXT2; //root(2)
    typedef CGAL::Sqrt_extension<Integer, Integer> EXT2_int; 

public:
    typedef EXT2 Extended_rational;

    typedef typename CGAL::Fraction_traits<Extended_rational>::Numerator_type
        Extended_coefficient;

    //! bivaritate polynomial over integers 
    typedef typename Algebraic_curve_kernel_2::Polynomial_2 Poly_int_2;
    //! univariate polynomial over sqrts
    typedef CGAL::Polynomial<Extended_coefficient> Poly_ext_1;
    //! bivaritate polynomial over sqrts
    typedef CGAL::Polynomial<Poly_ext_1> Poly_ext_2;

    //!\name functor invokation
    //!@{    

    // computes rotation angles
    void operator()(Extended_rational& esine, Extended_rational& ecosine,
         Extended_rational& ezero, int& angle) {
        
        if(angle % 45) 
            CGAL_error_msg("angle is not a multiple of 45");

        int angle_help = angle, sign = 1;
        if(angle_help >= 180) {
            sign = -1;
            angle_help -= 180;
        }

        _Packed<Rational> t;
        _15_degree_selector<Rational>(t, angle_help);

        esine = EXT2(t.a1_sine,t.b1_sine,Integer(2)); 

        ecosine = EXT2(t.a1_cosine,t.b1_cosine, Integer(2));

        ezero = EXT2(Rational(0),Rational(0),Integer(2));

        if(sign == -1) {
            esine = -esine;
            ecosine = -ecosine;
        }
    }
    //!@}
};

template<typename Poly_int_2, typename Poly_ext_2>
Poly_ext_2 _substitute_xy(const Poly_int_2& p, 
                          const Poly_ext_2& x, const Poly_ext_2& y) {
    
    typename Poly_int_2::const_iterator rit = p.end()-1;
    Poly_ext_2 r = rit->evaluate(x);
    while((rit--) != p.begin())
        r = r * y + rit->evaluate(x);
    return r;
}    

/*!\brief
 * defines coefficient number types and polynomial rotation functions for a
 * set of fixed angles
 */
template <class AlgebraicCurveKernel_2, int BaseAngle> 
class Rotation_traits {
    
public:
     //! this instance's template argument
    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;
    
    typedef Rotation_traits_base<Algebraic_curve_kernel_2, BaseAngle> Base;
    //! bivaritate polynomial over integers 
    typedef typename Base::Poly_int_2 Poly_int_2;
    //! bivaritate polynomial over sqrts
    typedef typename Base::Poly_ext_2 Poly_ext_2;

private:
    typedef typename Base::Arithmetic_kernel Arithmetic_kernel;

    typedef typename Arithmetic_kernel::Integer Integer;
    typedef typename Arithmetic_kernel::Rational Rational;

    typedef typename Base::Extended_rational Extended_rational;
    typedef CGAL::Polynomial<Extended_rational> Poly_sqrt_1;
    typedef CGAL::Polynomial<Poly_sqrt_1> Poly_sqrt_2;

    /* TODO: add rotations cache ?
     template <class KeyType_, class ValueType_,
        class Hash_ = boost::hash<KeyType_>,
        class Pred_ = std::equal_to<KeyType_>,
        class Canonicalizer_ = CGAL::Identity<KeyType_>,
        class Creator_ = CGAL::Creator_1<KeyType_, ValueType_> >
    class LRU_hashed_map
    
    typedef std::pair<Poly_int_2, int> Key_type;
    typedef CGALi::LRU_hashed_map<Key_type, Poly_ext_2,*/
    
public:
    //!\name public typedefs
    //!@{

    //! integer coefficients over sqrts
    typedef typename Base::Extended_coefficient Extended_coefficient;

    typedef Rebind_helper<Extended_coefficient, Rational> Rebind_helper;

    //! rebound kernel after substituting coefficient type
    typedef typename AlgebraicCurveKernel_2::template 
        rebind<typename Rebind_helper::Kernel_1>::Other Rebound_kernel;

    //!@}
    //!\name functor invokation 
    //!@{

    typedef Poly_int_2 argument_type;
    typedef Poly_ext_2 result_type;

    /*!\brief 
     * returns polynomial over sqrt extensions which corresponds to polynomial
     * \c poly_int rotated by \c angle degrees clockwise
     *
     * unless \c angle is multiple of \c BaseAngle , it is clamped to the 
     * nearest angle
     */
    Poly_ext_2 operator()(const Poly_int_2& poly_int, int angle) const {

        angle %= 360; 
        if(angle < 0) 
            angle += 360;

        Extended_rational esine, ecosine, ezero;
        Base engine;
        engine(esine, ecosine, ezero, angle);

        //std::cout << "nearest angle: " << angle << "\n";

        // CGAL_precondition_code macro cannot eat this
        typedef CGAL::Coercion_traits<Extended_rational,
                CGAL::Interval_nt<true> > CT;
        CGAL_precondition_code( // validity check
            typename CT::Cast cast;
            typedef typename CT::Type Type;
            Type esined = cast(esine);
            Type ecosined = cast(ecosine);
            Type ezerod = cast(ezero); 
            double sined = sin(angle * M_PI/180.0); 
            double cosined = cos(angle * M_PI/180.0);
            CGAL_precondition_msg(
                std::abs(sined - CGAL::to_double(esined)) < 1e-15 &&
                std::abs(cosined - CGAL::to_double(ecosined)) < 1e-15 && 
                std::abs(CGAL::to_double(ezerod)) < 1e-15,
                    "Wrond angle computed !!");
        );
        Poly_sqrt_2 sub_x(Poly_sqrt_1(ezero, ecosine), Poly_sqrt_1(esine)), 
            sub_y(Poly_sqrt_1(ezero, -esine), Poly_sqrt_1(ecosine)), res;

        res = _substitute_xy(poly_int, sub_x, sub_y);
        
        //std::cout << "rotated poly: " << res << std::endl;
        // integralize polynomial
        typedef CGAL::Fraction_traits<Poly_sqrt_2> FT;
        typename FT::Denominator_type dummy;
        Poly_ext_2 num;
        typename FT::Decompose()(res, num, dummy);
        
        //std::cout << "integralized poly: " << num << "\n\n";
        return num;
    }
    //!@}   

};

/*!\brief
 * required to prevent redundant template instantiations for angles
 * which are multiples of the same base angle
 */
template <class AlgebraicCurveKernel_2, int BaseAngle>
struct Rotated_algebraic_kernel_base :
    public Rotation_traits<AlgebraicCurveKernel_2, BaseAngle>::
            Rebound_kernel {

private:
    typedef Rotation_traits<AlgebraicCurveKernel_2, BaseAngle> Rotation_traits;
    
    typedef typename Rotation_traits::Rebound_kernel Base;

public:
    
    //! default constructor
    Rotated_algebraic_kernel_base() :
        Base() {
    }

    //! bivariate polynomial over integers
    typedef typename Rotation_traits::Poly_int_2 Poly_int_2;

    //! bivariate polynomial over sqrt-exts
    typedef typename Rotation_traits::Poly_ext_2 Poly_ext_2;

    typedef typename Base::Curve_analysis_2 Curve_analysis_2;
    
    struct Construct_curve_2 {
            
        Curve_analysis_2 operator()(const Poly_int_2& f, int angle=0) const {
            Rotation_traits traits;
            return Base::curve_cache_2()(traits(f, angle));
        }

        Curve_analysis_2 operator()(const Poly_ext_2& f) const {
            return Base::curve_cache_2()(f);
        };
    };
    
    Construct_curve_2 construct_curve_2_object() const {
        return Construct_curve_2();
    } 
  
    //Rotation_traits _m_traits;
};

template <class AlgebraicCurveKernel_2>
struct Approximately_rotated_algebraic_curve_kernel_base
    : AlgebraicCurveKernel_2 {

    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;
    typedef Algebraic_curve_kernel_2 Base;

    typedef typename Base::Curve_analysis_2 Curve_analysis_2;
    typedef typename Base::Polynomial_2 Polynomial_2;
    typedef typename Base::Boundary Rational;
    typedef typename Curve_analysis_2::Integer Integer;

    typedef CGAL::Polynomial<Rational> Poly_rat_1;
    typedef CGAL::Polynomial<Poly_rat_1> Poly_rat_2;

    struct Construct_curve_2 {
            
        Curve_analysis_2 operator()(const Polynomial_2& f, 
                                    Rational angle,
                                    long final_prec) {
            
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "angle=" << angle << std::endl;
            CGAL_ACK_DEBUG_PRINT << "final_prec=" << final_prec << std::endl;
#endif            

            while(abs(angle)>180) {
                angle>0 ?  angle-=360 : angle+=360;
            }
            std::cout << angle << std::endl;
            bool between_90_and_270 = abs(angle)>=90;
            bool greater_180 =(angle<0);
            // Normalize
            angle=abs(angle);
            if(between_90_and_270) {
                angle=180-angle;
            }
            CGAL_assertion(angle>=0 && angle <= 90);

            bool greater_45 = (angle>=45);
            if(greater_45) {
                angle=90-angle;
            }

            Rational sine, cosine;

            // Filter boundary case of 0 degree
            if(angle==Rational(0) || 
               Rational(1)/angle>CGAL::ipower(Integer(2),final_prec)) {
                sine = 0;
                cosine = 1;
            } else {

                typedef typename 
                    CGAL::Get_arithmetic_kernel<Integer>::Arithmetic_kernel AT;
                typedef typename AT::Bigfloat_interval Bigfloat_interval;
                typedef typename Bigfloat_interval_traits<Bigfloat_interval>
                    ::Boundary Bigfloat_boundary;


                long old_prec = CGAL::get_precision(Bigfloat_interval());

                long prec = 16;
                Rational t;
                while(true) {
                    CGAL::set_precision(Bigfloat_interval(),prec);
                    Bigfloat_interval pi = CGAL::pi<AT>(prec);
                    Bigfloat_interval s 
                        = CGAL::sin<AT>
                        (CGAL::median(
                                 pi*CGAL::convert_to_bfi(angle)/
                                 CGAL::convert_to_bfi(Integer(180))),
                         prec);
                    Bigfloat_boundary x 
                        = CGAL::median(CGAL::convert_to_bfi(Integer(1))/s + 
                                       CGAL::sqrt
                                       (CGAL::convert_to_bfi(Integer(1))/(s*s)-CGAL::convert_to_bfi(Integer(1))));

                    int n = 0;
                    typename CGAL::Fraction_traits<Rational>::Compose compose;
                    
                    bool success=false;
                    
                    Bigfloat_boundary e0=x, e1=-1, olde0;
                    Integer p0=0, q0=1, p1=1, q1=0,oldp0,oldq0;
                    while(true) {
                        Integer r = CGAL::CGALi::floor(e0/e1);
                        Integer oldp0=p0, oldq0=q0;
                        Bigfloat_boundary olde0=e0;
                        e0=e1;p0=p1;q0=q1;
                        e1=olde0-r*e1;
                        p1=oldp0-r*p1;
                        q1=oldq0-r*q1;
                        if(q1!=Integer(0)) {
                            t = compose(p1,q1);
                            CGAL::simplify(t);
                            sine = CGAL::abs(2/(t+1/t));
                            CGAL::simplify(sine);
                            Bigfloat_interval asin 
                                = CGAL::arcsin<AT>(sine,prec)
                                * CGAL::convert_to_bfi(Integer(180))/pi;
                            long bound = CGAL::CGALi::ceil_log2_abs
                                (CGAL::abs(asin-CGAL::convert_to_bfi(angle)));
                            success = (bound <= -final_prec);
                            typename 
                                CGAL::Coercion_traits<Rational, 
                                                      Bigfloat_boundary>::Cast
                                cast;
                            if((cast(t)==cast(x)) || success) {
                                break;
                            }
                        }
                        n++;
                    }
                    if(success) {
                        CGAL::set_precision(Bigfloat_interval(),old_prec);
                        break;
                    } else {
                        prec*=2;
                    }
                }
                
                cosine = (t-1/t)/(t+1/t);
                CGAL::simplify(cosine);
                
            }
            if(greater_45) {
                std::swap(sine,cosine);
            }
            
            if(greater_180) {
                sine = -sine;
            }
            if(between_90_and_270) {
                cosine=-cosine;
            }
            
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "sine=" << sine << std::endl;
            CGAL_ACK_DEBUG_PRINT << "cosine=" << cosine << std::endl;
#endif
            
            Poly_rat_2 
                sub_x(Poly_rat_1(Rational(0), cosine), Poly_rat_1(sine)), 
                    sub_y(Poly_rat_1(Rational(0), -sine), Poly_rat_1(cosine)), 
                res;
            
            res = _substitute_xy(f, sub_x, sub_y);
            
            CGAL::simplify(res);
            
            // integralize polynomial
            typedef CGAL::Fraction_traits<Poly_rat_2> FT;
            typename FT::Denominator_type dummy;
            Polynomial_2 num;
            typename FT::Decompose()(res, num, dummy);
            
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "integralized poly: " << num << std::endl;
#endif

            return Base::curve_cache_2()(num);
            }


        Curve_analysis_2 operator()(const Polynomial_2& f) const {
            return Base::curve_cache_2()(f);
        };
    
    };
    
    Construct_curve_2 construct_curve_2_object() const {
        return Construct_curve_2();
    } 
};

} // anonymous namespace

/*!\brief 
 *  defines \c Algebraic_curve_kernel_2 with rotation support 
 *
 *  \c BaseAngle (divisible by 3) specifies rotation traits which are used to
 *  compute otations by degrees multiple of \c BaseAngle
 */
template <class AlgebraicCurveKernel_2, int BaseAngle>
struct Rotated_algebraic_curve_kernel_2 :
    public Rotated_algebraic_kernel_base< AlgebraicCurveKernel_2,
        Normalized_angle< BaseAngle >::angle >
{ };

/*!\brief 
 *  defines \c Algebraic_curve_kernel_2 with rotation support for
 *  approximate rotations by arbitrary angles
 *
 * \Todo More documentation
 */
template <class AlgebraicCurveKernel_2>
struct Approximately_rotated_algebraic_curve_kernel_2   
    : public Approximately_rotated_algebraic_curve_kernel_base
        <AlgebraicCurveKernel_2>
{ 

    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;
    typedef Algebraic_curve_kernel_2 Base;
    
};


CGAL_END_NAMESPACE

#endif //CGAL_ROTATED_ALGEBRAIC_CURVE_KERNEL_2_H
