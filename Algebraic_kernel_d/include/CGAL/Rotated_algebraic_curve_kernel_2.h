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

#include <CGAL/Polynomial_type_generator.h>

CGAL_BEGIN_NAMESPACE

namespace { 

#ifdef CnX_USE_whatever_blabla

// put some CnX-specific code

#endif     

// For rotations by multiples of 3, we need 4 adjoint square roots, thus
// we work in a 16-dimensional Q-vector space. Each sin and cos is represented
// by the coefficients of the basis of this vector space
// x stands for sqrt(2)
// y stands for sqrt(3)
// z stands for sqrt(5)
// w stands for sqrt(10+2\sqrt(5))
template <class Rational>
struct Angle_coefficients {

    struct _Angle_coefficients_help {
        typedef _Angle_coefficients_help Self;
        Rational xyzw,xyz,xyw,xy,xzw,xz,xw,x,yzw,yz,yw,y,zw,z,w,c;
        Self operator-() {
            Self sn;
            sn.xyzw=-xyzw;sn.xyz=-xyz;sn.xyw=-xyw;sn.xy=-xy;sn.xzw=-xzw;
            sn.xz=-xz;sn.xw=-xw;sn.x=-x;sn.yzw=-yzw;sn.yz=-yz;sn.yw=-yw;
            sn.y=-y;sn.zw=-zw;sn.z=-z;sn.w=-w;sn.c=-c;
            return sn;
        }
        const Self& operator=(const Self& sn) {
            this->xyzw=sn.xyzw;this->xyz=sn.xyz;this->xyw=sn.xyw;
            this->xy=sn.xy;this->xzw=sn.xzw;this->xz=sn.xz;
            this->xw=sn.xw;this->x=sn.x;this->yzw=sn.yzw;this->yz=sn.yz;
            this->yw=sn.yw;this->y=sn.y;this->zw=sn.zw;
            this->z=sn.z;this->w=sn.w;this->c=sn.c;
            return (*this);
        }
    };
    
    Angle_coefficients(int angle) {
        select(angle);
    }

    void select(int angle) {
        CGAL_assertion(angle>=0 && angle<360);
        if(angle>180) {
            select(360-angle);
            this->sin = -this->sin;
        } else if(angle>90) {
            select(180-angle);
            this->cos = -this->cos;
        } else if(angle>45) {
            select(90-angle);
            _Angle_coefficients_help swap=sin;
            sin=cos;
            cos=swap;
        } else {
            CGAL_assertion(angle>=0 && angle<=45);
            typename CGAL::Fraction_traits<Rational>::Compose compose;
            switch(angle) { // these values have been computed by Maple
            case 0: {
                this->cos.xyzw=compose(0,1);
                this->cos.xyz=compose(0,1);
                this->cos.xyw=compose(0,1);
                this->cos.xy=compose(0,1);
                this->cos.xzw=compose(0,1);
                this->cos.xz=compose(0,1);
                this->cos.xw=compose(0,1);
                this->cos.x=compose(0,1);
                this->cos.yzw=compose(0,1);
                this->cos.yz=compose(0,1);
                this->cos.yw=compose(0,1);
                this->cos.y=compose(0,1);
                this->cos.zw=compose(0,1);
                this->cos.z=compose(0,1);
                this->cos.w=compose(0,1);
                this->cos.c=compose(1,1);
                this->sin.xyzw=compose(0,1);
                this->sin.xyz=compose(0,1);
                this->sin.xyw=compose(0,1);
                this->sin.xy=compose(0,1);
                this->sin.xzw=compose(0,1);
                this->sin.xz=compose(0,1);
                this->sin.xw=compose(0,1);
                this->sin.x=compose(0,1);
                this->sin.yzw=compose(0,1);
                this->sin.yz=compose(0,1);
                this->sin.yw=compose(0,1);
                this->sin.y=compose(0,1);
                this->sin.zw=compose(0,1);
                this->sin.z=compose(0,1);
                this->sin.w=compose(0,1);
                this->sin.c=compose(0,1);
                break;
            }
            case 3: {
                this->cos.xyzw=compose(0,1);
                this->cos.xyz=compose(1,16);
                this->cos.xyw=compose(1,16);
                this->cos.xy=compose(-1,16);
                this->cos.xzw=compose(0,1);
                this->cos.xz=compose(-1,16);
                this->cos.xw=compose(1,16);
                this->cos.x=compose(1,16);
                this->cos.yzw=compose(0,1);
                this->cos.yz=compose(0,1);
                this->cos.yw=compose(0,1);
                this->cos.y=compose(0,1);
                this->cos.zw=compose(0,1);
                this->cos.z=compose(0,1);
                this->cos.w=compose(0,1);
                this->cos.c=compose(0,1);
                this->sin.xyzw=compose(0,1);
                this->sin.xyz=compose(1,16);
                this->sin.xyw=compose(-1,16);
                this->sin.xy=compose(-1,16);
                this->sin.xzw=compose(0,1);
                this->sin.xz=compose(1,16);
                this->sin.xw=compose(1,16);
                this->sin.x=compose(-1,16);
                this->sin.yzw=compose(0,1);
                this->sin.yz=compose(0,1);
                this->sin.yw=compose(0,1);
                this->sin.y=compose(0,1);
                this->sin.zw=compose(0,1);
                this->sin.z=compose(0,1);
                this->sin.w=compose(0,1);
                this->sin.c=compose(0,1);
                break;
            }
            case 6: {
                this->cos.xyzw=compose(0,1);
                this->cos.xyz=compose(0,1);
                this->cos.xyw=compose(0,1);
                this->cos.xy=compose(0,1);
                this->cos.xzw=compose(0,1);
                this->cos.xz=compose(0,1);
                this->cos.xw=compose(0,1);
                this->cos.x=compose(0,1);
                this->cos.yzw=compose(0,1);
                this->cos.yz=compose(1,8);
                this->cos.yw=compose(0,1);
                this->cos.y=compose(1,8);
                this->cos.zw=compose(1,16);
                this->cos.z=compose(0,1);
                this->cos.w=compose(-1,16);
                this->cos.c=compose(0,1);
                this->sin.xyzw=compose(0,1);
                this->sin.xyz=compose(0,1);
                this->sin.xyw=compose(0,1);
                this->sin.xy=compose(0,1);
                this->sin.xzw=compose(0,1);
                this->sin.xz=compose(0,1);
                this->sin.xw=compose(0,1);
                this->sin.x=compose(0,1);
                this->sin.yzw=compose(1,16);
                this->sin.yz=compose(0,1);
                this->sin.yw=compose(-1,16);
                this->sin.y=compose(0,1);
                this->sin.zw=compose(0,1);
                this->sin.z=compose(-1,8);
                this->sin.w=compose(0,1);
                this->sin.c=compose(-1,8);
                break;
            }
            case 9: {
                this->cos.xyzw=compose(0,1);
                this->cos.xyz=compose(0,1);
                this->cos.xyw=compose(0,1);
                this->cos.xy=compose(0,1);
                this->cos.xzw=compose(1,16);
                this->cos.xz=compose(1,8);
                this->cos.xw=compose(-1,16);
                this->cos.x=compose(1,8);
                this->cos.yzw=compose(0,1);
                this->cos.yz=compose(0,1);
                this->cos.yw=compose(0,1);
                this->cos.y=compose(0,1);
                this->cos.zw=compose(0,1);
                this->cos.z=compose(0,1);
                this->cos.w=compose(0,1);
                this->cos.c=compose(0,1);
                this->sin.xyzw=compose(0,1);
                this->sin.xyz=compose(0,1);
                this->sin.xyw=compose(0,1);
                this->sin.xy=compose(0,1);
                this->sin.xzw=compose(-1,16);
                this->sin.xz=compose(1,8);
                this->sin.xw=compose(1,16);
                this->sin.x=compose(1,8);
                this->sin.yzw=compose(0,1);
                this->sin.yz=compose(0,1);
                this->sin.yw=compose(0,1);
                this->sin.y=compose(0,1);
                this->sin.zw=compose(0,1);
                this->sin.z=compose(0,1);
                this->sin.w=compose(0,1);
                this->sin.c=compose(0,1);
                break;
            }
            case 12: {
                this->cos.xyzw=compose(0,1);
                this->cos.xyz=compose(0,1);
                this->cos.xyw=compose(0,1);
                this->cos.xy=compose(0,1);
                this->cos.xzw=compose(0,1);
                this->cos.xz=compose(0,1);
                this->cos.xw=compose(0,1);
                this->cos.x=compose(0,1);
                this->cos.yzw=compose(0,1);
                this->cos.yz=compose(0,1);
                this->cos.yw=compose(1,8);
                this->cos.y=compose(0,1);
                this->cos.zw=compose(0,1);
                this->cos.z=compose(1,8);
                this->cos.w=compose(0,1);
                this->cos.c=compose(-1,8);
                this->sin.xyzw=compose(0,1);
                this->sin.xyz=compose(0,1);
                this->sin.xyw=compose(0,1);
                this->sin.xy=compose(0,1);
                this->sin.xzw=compose(0,1);
                this->sin.xz=compose(0,1);
                this->sin.xw=compose(0,1);
                this->sin.x=compose(0,1);
                this->sin.yzw=compose(0,1);
                this->sin.yz=compose(-1,8);
                this->sin.yw=compose(0,1);
                this->sin.y=compose(1,8);
                this->sin.zw=compose(0,1);
                this->sin.z=compose(0,1);
                this->sin.w=compose(1,8);
                this->sin.c=compose(0,1);
                break;
            }
            case 15: {
                this->cos.xyzw=compose(0,1);
                this->cos.xyz=compose(0,1);
                this->cos.xyw=compose(0,1);
                this->cos.xy=compose(1,4);
                this->cos.xzw=compose(0,1);
                this->cos.xz=compose(0,1);
                this->cos.xw=compose(0,1);
                this->cos.x=compose(1,4);
                this->cos.yzw=compose(0,1);
                this->cos.yz=compose(0,1);
                this->cos.yw=compose(0,1);
                this->cos.y=compose(0,1);
                this->cos.zw=compose(0,1);
                this->cos.z=compose(0,1);
                this->cos.w=compose(0,1);
                this->cos.c=compose(0,1);
                this->sin.xyzw=compose(0,1);
                this->sin.xyz=compose(0,1);
                this->sin.xyw=compose(0,1);
                this->sin.xy=compose(1,4);
                this->sin.xzw=compose(0,1);
                this->sin.xz=compose(0,1);
                this->sin.xw=compose(0,1);
                this->sin.x=compose(-1,4);
                this->sin.yzw=compose(0,1);
                this->sin.yz=compose(0,1);
                this->sin.yw=compose(0,1);
                this->sin.y=compose(0,1);
                this->sin.zw=compose(0,1);
                this->sin.z=compose(0,1);
                this->sin.w=compose(0,1);
                this->sin.c=compose(0,1);
                break;
            }
            case 18: {
                this->cos.xyzw=compose(0,1);
                this->cos.xyz=compose(0,1);
                this->cos.xyw=compose(0,1);
                this->cos.xy=compose(0,1);
                this->cos.xzw=compose(0,1);
                this->cos.xz=compose(0,1);
                this->cos.xw=compose(0,1);
                this->cos.x=compose(0,1);
                this->cos.yzw=compose(0,1);
                this->cos.yz=compose(0,1);
                this->cos.yw=compose(0,1);
                this->cos.y=compose(0,1);
                this->cos.zw=compose(0,1);
                this->cos.z=compose(0,1);
                this->cos.w=compose(1,4);
                this->cos.c=compose(0,1);
                this->sin.xyzw=compose(0,1);
                this->sin.xyz=compose(0,1);
                this->sin.xyw=compose(0,1);
                this->sin.xy=compose(0,1);
                this->sin.xzw=compose(0,1);
                this->sin.xz=compose(0,1);
                this->sin.xw=compose(0,1);
                this->sin.x=compose(0,1);
                this->sin.yzw=compose(0,1);
                this->sin.yz=compose(0,1);
                this->sin.yw=compose(0,1);
                this->sin.y=compose(0,1);
                this->sin.zw=compose(0,1);
                this->sin.z=compose(1,4);
                this->sin.w=compose(0,1);
                this->sin.c=compose(-1,4);
                break;
            }
            case 21: {
                this->cos.xyzw=compose(1,32);
                this->cos.xyz=compose(1,16);
                this->cos.xyw=compose(-1,32);
                this->cos.xy=compose(1,16);
                this->cos.xzw=compose(-1,32);
                this->cos.xz=compose(1,16);
                this->cos.xw=compose(1,32);
                this->cos.x=compose(1,16);
                this->cos.yzw=compose(0,1);
                this->cos.yz=compose(0,1);
                this->cos.yw=compose(0,1);
                this->cos.y=compose(0,1);
                this->cos.zw=compose(0,1);
                this->cos.z=compose(0,1);
                this->cos.w=compose(0,1);
                this->cos.c=compose(0,1);
                this->sin.xyzw=compose(1,32);
                this->sin.xyz=compose(-1,16);
                this->sin.xyw=compose(-1,32);
                this->sin.xy=compose(-1,16);
                this->sin.xzw=compose(1,32);
                this->sin.xz=compose(1,16);
                this->sin.xw=compose(-1,32);
                this->sin.x=compose(1,16);
                this->sin.yzw=compose(0,1);
                this->sin.yz=compose(0,1);
                this->sin.yw=compose(0,1);
                this->sin.y=compose(0,1);
                this->sin.zw=compose(0,1);
                this->sin.z=compose(0,1);
                this->sin.w=compose(0,1);
                this->sin.c=compose(0,1);
                break;
            }
            case 24: {
                this->cos.xyzw=compose(0,1);
                this->cos.xyz=compose(0,1);
                this->cos.xyw=compose(0,1);
                this->cos.xy=compose(0,1);
                this->cos.xzw=compose(0,1);
                this->cos.xz=compose(0,1);
                this->cos.xw=compose(0,1);
                this->cos.x=compose(0,1);
                this->cos.yzw=compose(1,16);
                this->cos.yz=compose(0,1);
                this->cos.yw=compose(-1,16);
                this->cos.y=compose(0,1);
                this->cos.zw=compose(0,1);
                this->cos.z=compose(1,8);
                this->cos.w=compose(0,1);
                this->cos.c=compose(1,8);
                this->sin.xyzw=compose(0,1);
                this->sin.xyz=compose(0,1);
                this->sin.xyw=compose(0,1);
                this->sin.xy=compose(0,1);
                this->sin.xzw=compose(0,1);
                this->sin.xz=compose(0,1);
                this->sin.xw=compose(0,1);
                this->sin.x=compose(0,1);
                this->sin.yzw=compose(0,1);
                this->sin.yz=compose(1,8);
                this->sin.yw=compose(0,1);
                this->sin.y=compose(1,8);
                this->sin.zw=compose(-1,16);
                this->sin.z=compose(0,1);
                this->sin.w=compose(1,16);
                this->sin.c=compose(0,1);
                break;
            }
            case 27: {
                this->cos.xyzw=compose(0,1);
                this->cos.xyz=compose(0,1);
                this->cos.xyw=compose(0,1);
                this->cos.xy=compose(0,1);
                this->cos.xzw=compose(0,1);
                this->cos.xz=compose(1,8);
                this->cos.xw=compose(1,8);
                this->cos.x=compose(-1,8);
                this->cos.yzw=compose(0,1);
                this->cos.yz=compose(0,1);
                this->cos.yw=compose(0,1);
                this->cos.y=compose(0,1);
                this->cos.zw=compose(0,1);
                this->cos.z=compose(0,1);
                this->cos.w=compose(0,1);
                this->cos.c=compose(0,1);
                this->sin.xyzw=compose(0,1);
                this->sin.xyz=compose(0,1);
                this->sin.xyw=compose(0,1);
                this->sin.xy=compose(0,1);
                this->sin.xzw=compose(0,1);
                this->sin.xz=compose(-1,8);
                this->sin.xw=compose(1,8);
                this->sin.x=compose(1,8);
                this->sin.yzw=compose(0,1);
                this->sin.yz=compose(0,1);
                this->sin.yw=compose(0,1);
                this->sin.y=compose(0,1);
                this->sin.zw=compose(0,1);
                this->sin.z=compose(0,1);
                this->sin.w=compose(0,1);
                this->sin.c=compose(0,1);
                break;
            }
            case 30: {
                this->cos.xyzw=compose(0,1);
                this->cos.xyz=compose(0,1);
                this->cos.xyw=compose(0,1);
                this->cos.xy=compose(0,1);
                this->cos.xzw=compose(0,1);
                this->cos.xz=compose(0,1);
                this->cos.xw=compose(0,1);
                this->cos.x=compose(0,1);
                this->cos.yzw=compose(0,1);
                this->cos.yz=compose(0,1);
                this->cos.yw=compose(0,1);
                this->cos.y=compose(1,2);
                this->cos.zw=compose(0,1);
                this->cos.z=compose(0,1);
                this->cos.w=compose(0,1);
                this->cos.c=compose(0,1);
                this->sin.xyzw=compose(0,1);
                this->sin.xyz=compose(0,1);
                this->sin.xyw=compose(0,1);
                this->sin.xy=compose(0,1);
                this->sin.xzw=compose(0,1);
                this->sin.xz=compose(0,1);
                this->sin.xw=compose(0,1);
                this->sin.x=compose(0,1);
                this->sin.yzw=compose(0,1);
                this->sin.yz=compose(0,1);
                this->sin.yw=compose(0,1);
                this->sin.y=compose(0,1);
                this->sin.zw=compose(0,1);
                this->sin.z=compose(0,1);
                this->sin.w=compose(0,1);
                this->sin.c=compose(1,2);
                break;
            }
            case 33: {
                this->cos.xyzw=compose(0,1);
                this->cos.xyz=compose(-1,16);
                this->cos.xyw=compose(1,16);
                this->cos.xy=compose(1,16);
                this->cos.xzw=compose(0,1);
                this->cos.xz=compose(1,16);
                this->cos.xw=compose(1,16);
                this->cos.x=compose(-1,16);
                this->cos.yzw=compose(0,1);
                this->cos.yz=compose(0,1);
                this->cos.yw=compose(0,1);
                this->cos.y=compose(0,1);
                this->cos.zw=compose(0,1);
                this->cos.z=compose(0,1);
                this->cos.w=compose(0,1);
                this->cos.c=compose(0,1);
                this->sin.xyzw=compose(0,1);
                this->sin.xyz=compose(1,16);
                this->sin.xyw=compose(1,16);
                this->sin.xy=compose(-1,16);
                this->sin.xzw=compose(0,1);
                this->sin.xz=compose(1,16);
                this->sin.xw=compose(-1,16);
                this->sin.x=compose(-1,16);
                this->sin.yzw=compose(0,1);
                this->sin.yz=compose(0,1);
                this->sin.yw=compose(0,1);
                this->sin.y=compose(0,1);
                this->sin.zw=compose(0,1);
                this->sin.z=compose(0,1);
                this->sin.w=compose(0,1);
                this->sin.c=compose(0,1);
                break;
            }
            case 36: {
                this->cos.xyzw=compose(0,1);
                this->cos.xyz=compose(0,1);
                this->cos.xyw=compose(0,1);
                this->cos.xy=compose(0,1);
                this->cos.xzw=compose(0,1);
                this->cos.xz=compose(0,1);
                this->cos.xw=compose(0,1);
                this->cos.x=compose(0,1);
                this->cos.yzw=compose(0,1);
                this->cos.yz=compose(0,1);
                this->cos.yw=compose(0,1);
                this->cos.y=compose(0,1);
                this->cos.zw=compose(0,1);
                this->cos.z=compose(1,4);
                this->cos.w=compose(0,1);
                this->cos.c=compose(1,4);
                this->sin.xyzw=compose(0,1);
                this->sin.xyz=compose(0,1);
                this->sin.xyw=compose(0,1);
                this->sin.xy=compose(0,1);
                this->sin.xzw=compose(0,1);
                this->sin.xz=compose(0,1);
                this->sin.xw=compose(0,1);
                this->sin.x=compose(0,1);
                this->sin.yzw=compose(0,1);
                this->sin.yz=compose(0,1);
                this->sin.yw=compose(0,1);
                this->sin.y=compose(0,1);
                this->sin.zw=compose(1,8);
                this->sin.z=compose(0,1);
                this->sin.w=compose(-1,8);
                this->sin.c=compose(0,1);
                break;
            }
            case 39: {
                this->cos.xyzw=compose(1,32);
                this->cos.xyz=compose(1,16);
                this->cos.xyw=compose(-1,32);
                this->cos.xy=compose(1,16);
                this->cos.xzw=compose(1,32);
                this->cos.xz=compose(-1,16);
                this->cos.xw=compose(-1,32);
                this->cos.x=compose(-1,16);
                this->cos.yzw=compose(0,1);
                this->cos.yz=compose(0,1);
                this->cos.yw=compose(0,1);
                this->cos.y=compose(0,1);
                this->cos.zw=compose(0,1);
                this->cos.z=compose(0,1);
                this->cos.w=compose(0,1);
                this->cos.c=compose(0,1);
                this->sin.xyzw=compose(-1,32);
                this->sin.xyz=compose(1,16);
                this->sin.xyw=compose(1,32);
                this->sin.xy=compose(1,16);
                this->sin.xzw=compose(1,32);
                this->sin.xz=compose(1,16);
                this->sin.xw=compose(-1,32);
                this->sin.x=compose(1,16);
                this->sin.yzw=compose(0,1);
                this->sin.yz=compose(0,1);
                this->sin.yw=compose(0,1);
                this->sin.y=compose(0,1);
                this->sin.zw=compose(0,1);
                this->sin.z=compose(0,1);
                this->sin.w=compose(0,1);
                this->sin.c=compose(0,1);
                break;
            }
            case 42: {
                this->cos.xyzw=compose(0,1);
                this->cos.xyz=compose(0,1);
                this->cos.xyw=compose(0,1);
                this->cos.xy=compose(0,1);
                this->cos.xzw=compose(0,1);
                this->cos.xz=compose(0,1);
                this->cos.xw=compose(0,1);
                this->cos.x=compose(0,1);
                this->cos.yzw=compose(0,1);
                this->cos.yz=compose(1,8);
                this->cos.yw=compose(0,1);
                this->cos.y=compose(-1,8);
                this->cos.zw=compose(0,1);
                this->cos.z=compose(0,1);
                this->cos.w=compose(1,8);
                this->cos.c=compose(0,1);
                this->sin.xyzw=compose(0,1);
                this->sin.xyz=compose(0,1);
                this->sin.xyw=compose(0,1);
                this->sin.xy=compose(0,1);
                this->sin.xzw=compose(0,1);
                this->sin.xz=compose(0,1);
                this->sin.xw=compose(0,1);
                this->sin.x=compose(0,1);
                this->sin.yzw=compose(0,1);
                this->sin.yz=compose(0,1);
                this->sin.yw=compose(1,8);
                this->sin.y=compose(0,1);
                this->sin.zw=compose(0,1);
                this->sin.z=compose(-1,8);
                this->sin.w=compose(0,1);
                this->sin.c=compose(1,8);
                break;
            }
            case 45: {
                this->cos.xyzw=compose(0,1);
                this->cos.xyz=compose(0,1);
                this->cos.xyw=compose(0,1);
                this->cos.xy=compose(0,1);
                this->cos.xzw=compose(0,1);
                this->cos.xz=compose(0,1);
                this->cos.xw=compose(0,1);
                this->cos.x=compose(1,2);
                this->cos.yzw=compose(0,1);
                this->cos.yz=compose(0,1);
                this->cos.yw=compose(0,1);
                this->cos.y=compose(0,1);
                this->cos.zw=compose(0,1);
                this->cos.z=compose(0,1);
                this->cos.w=compose(0,1);
                this->cos.c=compose(0,1);
                this->sin.xyzw=compose(0,1);
                this->sin.xyz=compose(0,1);
                this->sin.xyw=compose(0,1);
                this->sin.xy=compose(0,1);
                this->sin.xzw=compose(0,1);
                this->sin.xz=compose(0,1);
                this->sin.xw=compose(0,1);
                this->sin.x=compose(1,2);
                this->sin.yzw=compose(0,1);
                this->sin.yz=compose(0,1);
                this->sin.yw=compose(0,1);
                this->sin.y=compose(0,1);
                this->sin.zw=compose(0,1);
                this->sin.z=compose(0,1);
                this->sin.w=compose(0,1);
                this->sin.c=compose(0,1);
                break;
            }

            } // of switch

        } // of else {
    }

    _Angle_coefficients_help sin,cos;
};

template <typename ArithmeticKernel,int base_angle>
struct Rotation_traits_for_base_angle_base {
    typedef Null_functor Rotated_coefficient;
    typedef Null_functor Rotated_rational_coefficient;
    typedef Null_functor Sin;
    typedef Null_functor Cos;
    typedef Null_functor Zero;
};

// ZZ[sqrt(2)]
template <typename Polynomial_2>
struct Rotation_traits_for_base_angle_base<Polynomial_2,45> {
    

    typedef Polynomial_2 Unrotated_polynomial_2;

    typedef typename CGAL::Polynomial_traits_d<Unrotated_polynomial_2>
        ::Innermost_coefficient_type Integer;

    typedef typename CGAL::Get_arithmetic_kernel<Integer>::Arithmetic_kernel
    Arithmetic_kernel;

    typedef typename Arithmetic_kernel::Rational Rational;

    typedef CGAL::Sqrt_extension<Rational, Integer> 
        Rotated_rational_coefficient; //root(2)
    typedef CGAL::Sqrt_extension<Integer, Integer> 
        Rotated_coefficient;

    struct Sin : public std::unary_function<int,Rotated_rational_coefficient> {

        Rotated_rational_coefficient operator() (int angle) {
            
            Angle_coefficients<Rational> coeffs(angle);
            CGAL_assertion_code(
              typename 
              CGAL::Algebraic_structure_traits<Rational>::Is_zero zero;
            );
            CGAL_assertion(zero(coeffs.sin.xyzw));
            CGAL_assertion(zero(coeffs.sin.xyz));
            CGAL_assertion(zero(coeffs.sin.xy));
            CGAL_assertion(zero(coeffs.sin.xzw));
            CGAL_assertion(zero(coeffs.sin.xyw));
            CGAL_assertion(zero(coeffs.sin.xz));
            CGAL_assertion(zero(coeffs.sin.xw));
            CGAL_assertion(zero(coeffs.sin.y));
            CGAL_assertion(zero(coeffs.sin.yz));
            CGAL_assertion(zero(coeffs.sin.yw));
            CGAL_assertion(zero(coeffs.sin.yzw));
            CGAL_assertion(zero(coeffs.sin.zw));
            CGAL_assertion(zero(coeffs.sin.z));
            CGAL_assertion(zero(coeffs.sin.w));
            return Rotated_rational_coefficient(coeffs.sin.c,
                                                coeffs.sin.x,
                                                Integer(2));
        }

    };

    struct Cos : public std::unary_function<int,Rotated_rational_coefficient> {

        Rotated_rational_coefficient operator() (int angle) {
            
            Angle_coefficients<Rational> coeffs(angle);
            CGAL_assertion_code(
                    typename 
                    CGAL::Algebraic_structure_traits<Rational>::Is_zero zero;
            );
            CGAL_assertion(zero(coeffs.cos.xyzw));
            CGAL_assertion(zero(coeffs.cos.xyz));
            CGAL_assertion(zero(coeffs.cos.xy));
            CGAL_assertion(zero(coeffs.cos.xzw));
            CGAL_assertion(zero(coeffs.cos.xyw));
            CGAL_assertion(zero(coeffs.cos.xz));
            CGAL_assertion(zero(coeffs.cos.xw));
            CGAL_assertion(zero(coeffs.cos.y));
            CGAL_assertion(zero(coeffs.cos.yz));
            CGAL_assertion(zero(coeffs.cos.yw));
            CGAL_assertion(zero(coeffs.cos.yzw));
            CGAL_assertion(zero(coeffs.cos.zw));
            CGAL_assertion(zero(coeffs.cos.z));
            CGAL_assertion(zero(coeffs.cos.w));
            return Rotated_rational_coefficient(coeffs.cos.c,
                                                coeffs.cos.x,
                                                Integer(2));
        }

    };
    
    struct Zero {
        
        typedef Rotated_rational_coefficient result_type;

        Rotated_rational_coefficient operator() () {
            Rational zero(0);
            return Rotated_rational_coefficient(zero,zero,Integer(2));
        }
    };        

};

// ZZ[sqrt(3)]
template <typename Polynomial_2>
struct Rotation_traits_for_base_angle_base<Polynomial_2,30> {
    
    typedef Polynomial_2 Unrotated_polynomial_2;

    typedef typename CGAL::Polynomial_traits_d<Unrotated_polynomial_2>
        ::Innermost_coefficient_type Integer;

    typedef typename CGAL::Get_arithmetic_kernel<Integer>::Arithmetic_kernel
    Arithmetic_kernel;

    typedef typename Arithmetic_kernel::Rational Rational;

    typedef CGAL::Sqrt_extension<Rational, Integer> 
        Rotated_rational_coefficient; //root(3)
    typedef CGAL::Sqrt_extension<Integer, Integer> 
        Rotated_coefficient;

    struct Sin : public std::unary_function<int,Rotated_rational_coefficient> {

        Rotated_rational_coefficient operator() (int angle) {
            
            Angle_coefficients<Rational> coeffs(angle);
            CGAL_assertion_code(
                    typename 
                    CGAL::Algebraic_structure_traits<Rational>::Is_zero zero;
            );
            CGAL_assertion(zero(coeffs.sin.xyzw));
            CGAL_assertion(zero(coeffs.sin.xyz));
            CGAL_assertion(zero(coeffs.sin.xy));
            CGAL_assertion(zero(coeffs.sin.xzw));
            CGAL_assertion(zero(coeffs.sin.xyw));
            CGAL_assertion(zero(coeffs.sin.xz));
            CGAL_assertion(zero(coeffs.sin.xw));
            CGAL_assertion(zero(coeffs.sin.x));
            CGAL_assertion(zero(coeffs.sin.yz));
            CGAL_assertion(zero(coeffs.sin.yw));
            CGAL_assertion(zero(coeffs.sin.yzw));
            CGAL_assertion(zero(coeffs.sin.zw));
            CGAL_assertion(zero(coeffs.sin.z));
            CGAL_assertion(zero(coeffs.sin.w));
            return Rotated_rational_coefficient(coeffs.sin.c,
                                                coeffs.sin.y,
                                                Integer(3));
        }

    };

    struct Cos : public std::unary_function<int,Rotated_rational_coefficient> {

        Rotated_rational_coefficient operator() (int angle) {
            
            Angle_coefficients<Rational> coeffs(angle);
            CGAL_assertion_code(
                    typename 
                    CGAL::Algebraic_structure_traits<Rational>::Is_zero zero;
            );
            CGAL_assertion(zero(coeffs.cos.xyzw));
            CGAL_assertion(zero(coeffs.cos.xyz));
            CGAL_assertion(zero(coeffs.cos.xy));
            CGAL_assertion(zero(coeffs.cos.xzw));
            CGAL_assertion(zero(coeffs.cos.xyw));
            CGAL_assertion(zero(coeffs.cos.xz));
            CGAL_assertion(zero(coeffs.cos.xw));
            CGAL_assertion(zero(coeffs.cos.x));
            CGAL_assertion(zero(coeffs.cos.yz));
            CGAL_assertion(zero(coeffs.cos.yw));
            CGAL_assertion(zero(coeffs.cos.yzw));
            CGAL_assertion(zero(coeffs.cos.zw));
            CGAL_assertion(zero(coeffs.cos.z));
            CGAL_assertion(zero(coeffs.cos.w));
            return Rotated_rational_coefficient(coeffs.cos.c,
                                                coeffs.cos.y,
                                                Integer(3));
        }

    };
    
    struct Zero {
        
        typedef Rotated_rational_coefficient result_type;

        Rotated_rational_coefficient operator() () {
            Rational zero(0);
            return Rotated_rational_coefficient(zero,zero,Integer(3));
        }
    };        

};

// ZZ[sqrt(5),sqrt(10+2sqrt(5)))]
template <typename Polynomial_2>
struct Rotation_traits_for_base_angle_base<Polynomial_2,18> {
    
    typedef Polynomial_2 Unrotated_polynomial_2;

    typedef typename CGAL::Polynomial_traits_d<Unrotated_polynomial_2>
        ::Innermost_coefficient_type Integer;

    typedef typename CGAL::Get_arithmetic_kernel<Integer>::Arithmetic_kernel
    Arithmetic_kernel;

    typedef typename Arithmetic_kernel::Rational Rational;

private:

    typedef CGAL::Sqrt_extension<Rational,Integer> Rat_with_sqrt_5;
    typedef CGAL::Sqrt_extension<Integer,Integer> Int_with_sqrt_5;

public:

    typedef CGAL::Sqrt_extension<Rat_with_sqrt_5,Int_with_sqrt_5> 
        Rotated_rational_coefficient; 
    typedef CGAL::Sqrt_extension<Int_with_sqrt_5,Int_with_sqrt_5> 
        Rotated_coefficient;

    struct Sin : public std::unary_function<int,Rotated_rational_coefficient> {

        Rotated_rational_coefficient operator() (int angle) {
            
            Angle_coefficients<Rational> coeffs(angle);
            CGAL_assertion_code(
                    typename 
                    CGAL::Algebraic_structure_traits<Rational>::Is_zero zero;
            );
            CGAL_assertion(zero(coeffs.sin.xyzw));
            CGAL_assertion(zero(coeffs.sin.xyz));
            CGAL_assertion(zero(coeffs.sin.xzw));
            CGAL_assertion(zero(coeffs.sin.xyw));
            CGAL_assertion(zero(coeffs.sin.xz));
            CGAL_assertion(zero(coeffs.sin.xw));
            CGAL_assertion(zero(coeffs.sin.xy));
            CGAL_assertion(zero(coeffs.sin.x));
            CGAL_assertion(zero(coeffs.sin.y));
            CGAL_assertion(zero(coeffs.sin.yz));
            CGAL_assertion(zero(coeffs.sin.yw));
            CGAL_assertion(zero(coeffs.sin.yzw));
            return Rotated_rational_coefficient
                (Rat_with_sqrt_5(coeffs.sin.c,coeffs.sin.z,Integer(5)),
                 Rat_with_sqrt_5(coeffs.sin.w,coeffs.sin.zw,Integer(5)),
                 Int_with_sqrt_5(Integer(10),Integer(2),Integer(5)));
        }

    };

    struct Cos : public std::unary_function<int,Rotated_rational_coefficient> {

        Rotated_rational_coefficient operator() (int angle) {
            
            Angle_coefficients<Rational> coeffs(angle);
            CGAL_assertion_code(
                    typename 
                    CGAL::Algebraic_structure_traits<Rational>::Is_zero zero;
            );
            CGAL_assertion(zero(coeffs.cos.xyzw));
            CGAL_assertion(zero(coeffs.cos.xyz));
            CGAL_assertion(zero(coeffs.cos.xzw));
            CGAL_assertion(zero(coeffs.cos.xyw));
            CGAL_assertion(zero(coeffs.cos.xz));
            CGAL_assertion(zero(coeffs.cos.xw));
            CGAL_assertion(zero(coeffs.cos.xy));
            CGAL_assertion(zero(coeffs.cos.x));
            CGAL_assertion(zero(coeffs.cos.y));
            CGAL_assertion(zero(coeffs.cos.yz));
            CGAL_assertion(zero(coeffs.cos.yw));
            CGAL_assertion(zero(coeffs.cos.yzw));
            return Rotated_rational_coefficient
                (Rat_with_sqrt_5(coeffs.cos.c,coeffs.cos.z,Integer(5)),
                 Rat_with_sqrt_5(coeffs.cos.w,coeffs.cos.zw,Integer(5)),
                 Int_with_sqrt_5(Integer(10),Integer(2),Integer(5)));
        }

    };
    
    struct Zero {
        
        typedef Rotated_rational_coefficient result_type;

        Rotated_rational_coefficient operator() () {
            Rational zero(0);
            return Rotated_rational_coefficient
                (Rat_with_sqrt_5(zero,zero,Integer(5)),
                 Rat_with_sqrt_5(zero,zero,Integer(5)),
                 Int_with_sqrt_5(Integer(10),Integer(2),Integer(5)));
        }
    };        

};

// ZZ[sqrt(2),sqrt(3)]
template <typename Polynomial_2>
struct Rotation_traits_for_base_angle_base<Polynomial_2,15> {
    
    typedef Polynomial_2 Unrotated_polynomial_2;

    typedef typename CGAL::Polynomial_traits_d<Unrotated_polynomial_2>
        ::Innermost_coefficient_type Integer;

    typedef typename CGAL::Get_arithmetic_kernel<Integer>::Arithmetic_kernel
    Arithmetic_kernel;

    typedef typename Arithmetic_kernel::Rational Rational;

private:

    typedef CGAL::Sqrt_extension<Rational,Integer> Rat_with_sqrt_2;
    typedef CGAL::Sqrt_extension<Integer,Integer> Int_with_sqrt_2;

public:

    typedef CGAL::Sqrt_extension<Rat_with_sqrt_2,Integer> 
        Rotated_rational_coefficient; 
    typedef CGAL::Sqrt_extension<Int_with_sqrt_2,Integer> 
        Rotated_coefficient;

    struct Sin : public std::unary_function<int,Rotated_rational_coefficient> {

        Rotated_rational_coefficient operator() (int angle) {
            
            Angle_coefficients<Rational> coeffs(angle);
            CGAL_assertion_code(
                    typename 
                    CGAL::Algebraic_structure_traits<Rational>::Is_zero zero;
            );
            CGAL_assertion(zero(coeffs.sin.xyzw));
            CGAL_assertion(zero(coeffs.sin.xyz));
            CGAL_assertion(zero(coeffs.sin.xzw));
            CGAL_assertion(zero(coeffs.sin.xyw));
            CGAL_assertion(zero(coeffs.sin.xz));
            CGAL_assertion(zero(coeffs.sin.xw));
            CGAL_assertion(zero(coeffs.sin.yz));
            CGAL_assertion(zero(coeffs.sin.yw));
            CGAL_assertion(zero(coeffs.sin.yzw));
            CGAL_assertion(zero(coeffs.sin.zw));
            CGAL_assertion(zero(coeffs.sin.z));
            CGAL_assertion(zero(coeffs.sin.w));
            return Rotated_rational_coefficient
                (Rat_with_sqrt_2(coeffs.sin.c,coeffs.sin.x,Integer(2)),
                 Rat_with_sqrt_2(coeffs.sin.y,coeffs.sin.xy,Integer(2)),
                 Integer(3));
        }

    };

    struct Cos : public std::unary_function<int,Rotated_rational_coefficient> {

        Rotated_rational_coefficient operator() (int angle) {
            
            Angle_coefficients<Rational> coeffs(angle);
            CGAL_assertion_code(
                    typename 
                    CGAL::Algebraic_structure_traits<Rational>::Is_zero zero;
            );
            CGAL_assertion(zero(coeffs.cos.xyzw));
            CGAL_assertion(zero(coeffs.cos.xyz));
            CGAL_assertion(zero(coeffs.cos.xzw));
            CGAL_assertion(zero(coeffs.cos.xyw));
            CGAL_assertion(zero(coeffs.cos.xz));
            CGAL_assertion(zero(coeffs.cos.xw));
            CGAL_assertion(zero(coeffs.cos.yz));
            CGAL_assertion(zero(coeffs.cos.yw));
            CGAL_assertion(zero(coeffs.cos.yzw));
            CGAL_assertion(zero(coeffs.cos.zw));
            CGAL_assertion(zero(coeffs.cos.z));
            CGAL_assertion(zero(coeffs.cos.w));
            return Rotated_rational_coefficient
                (Rat_with_sqrt_2(coeffs.cos.c,coeffs.cos.x,Integer(2)),
                 Rat_with_sqrt_2(coeffs.cos.y,coeffs.cos.xy,Integer(2)),
                 Integer(3));
        }

    };
    
    struct Zero {
        
        typedef Rotated_rational_coefficient result_type;

        Rotated_rational_coefficient operator() () {
            Rational zero(0);
            return Rotated_rational_coefficient
                (Rat_with_sqrt_2(zero,zero,Integer(2)),
                 Rat_with_sqrt_2(zero,zero,Integer(2)),
                 Integer(3));
        }
    };        

};

// ZZ[sqrt(2),sqrt(5),sqrt(10+2sqrt(5))]
template <typename Polynomial_2>
struct Rotation_traits_for_base_angle_base<Polynomial_2,9> {
    
    typedef Polynomial_2 Unrotated_polynomial_2;

    typedef typename CGAL::Polynomial_traits_d<Unrotated_polynomial_2>
        ::Innermost_coefficient_type Integer;

    typedef typename CGAL::Get_arithmetic_kernel<Integer>::Arithmetic_kernel
    Arithmetic_kernel;

    typedef typename Arithmetic_kernel::Rational Rational;

private:

    typedef CGAL::Sqrt_extension<Rational,Integer> Rat_with_sqrt_2;
    typedef CGAL::Sqrt_extension<Integer,Integer> Int_with_sqrt_2;
    typedef CGAL::Sqrt_extension<Integer,Integer> Int_with_sqrt_5;
    typedef CGAL::Sqrt_extension<Rat_with_sqrt_2,Integer> Rat_with_sqrt_2_5;
    typedef CGAL::Sqrt_extension<Int_with_sqrt_2,Integer> Int_with_sqrt_2_5;

public:

    typedef CGAL::Sqrt_extension<Rat_with_sqrt_2_5,Int_with_sqrt_2_5> 
        Rotated_rational_coefficient; 
    typedef CGAL::Sqrt_extension<Int_with_sqrt_2_5,Int_with_sqrt_2_5> 
        Rotated_coefficient;

    struct Sin : public std::unary_function<int,Rotated_rational_coefficient> {

        Rotated_rational_coefficient operator() (int angle) {
            
            Angle_coefficients<Rational> coeffs(angle);
            CGAL_assertion_code(
                    typename 
                    CGAL::Algebraic_structure_traits<Rational>::Is_zero zero;
            );
            CGAL_assertion(zero(coeffs.sin.xyzw));
            CGAL_assertion(zero(coeffs.sin.xyz));
            CGAL_assertion(zero(coeffs.sin.xyw));
            CGAL_assertion(zero(coeffs.sin.xy));
            CGAL_assertion(zero(coeffs.sin.yzw));
            CGAL_assertion(zero(coeffs.sin.yz));
            CGAL_assertion(zero(coeffs.sin.yw));
            CGAL_assertion(zero(coeffs.sin.y));
            return Rotated_rational_coefficient
                (Rat_with_sqrt_2_5
                 (Rat_with_sqrt_2(coeffs.sin.c,coeffs.sin.x,Integer(2)),
                  Rat_with_sqrt_2(coeffs.sin.z,coeffs.sin.xz,Integer(2)),
                  Integer(5)),
                 Rat_with_sqrt_2_5
                 (Rat_with_sqrt_2(coeffs.sin.w,coeffs.sin.xw,Integer(2)),
                  Rat_with_sqrt_2(coeffs.sin.zw,coeffs.sin.xzw,Integer(2)),
                  Integer(5)),
                 Int_with_sqrt_2_5
                 (Int_with_sqrt_2(Integer(10),Integer(0),Integer(2)),
                  Int_with_sqrt_2(Integer(2),Integer(0),Integer(2)),
                  Integer(5)));
        }

    };

    struct Cos : public std::unary_function<int,Rotated_rational_coefficient> {

        Rotated_rational_coefficient operator() (int angle) {
            Angle_coefficients<Rational> coeffs(angle);
            CGAL_assertion_code(
                    typename 
                    CGAL::Algebraic_structure_traits<Rational>::Is_zero zero;
            );

            CGAL_assertion(zero(coeffs.cos.xyzw));
            CGAL_assertion(zero(coeffs.cos.xyz));
            CGAL_assertion(zero(coeffs.cos.xyw));
            CGAL_assertion(zero(coeffs.cos.xy));
            CGAL_assertion(zero(coeffs.cos.yzw));
            CGAL_assertion(zero(coeffs.cos.yz));
            CGAL_assertion(zero(coeffs.cos.yw));
            CGAL_assertion(zero(coeffs.cos.y));
            return Rotated_rational_coefficient
                (Rat_with_sqrt_2_5
                 (Rat_with_sqrt_2(coeffs.cos.c,coeffs.cos.x,Integer(2)),
                  Rat_with_sqrt_2(coeffs.cos.z,coeffs.cos.xz,Integer(2)),
                  Integer(5)),
                 Rat_with_sqrt_2_5
                 (Rat_with_sqrt_2(coeffs.cos.w,coeffs.cos.xw,Integer(2)),
                  Rat_with_sqrt_2(coeffs.cos.zw,coeffs.cos.xzw,Integer(2)),
                  Integer(5)),
                 Int_with_sqrt_2_5
                 (Int_with_sqrt_2(Integer(10),Integer(0),Integer(2)),
                  Int_with_sqrt_2(Integer(2),Integer(0),Integer(2)),
                  Integer(5)));

        }

    };
    
    struct Zero {
        
        typedef Rotated_rational_coefficient result_type;

        Rotated_rational_coefficient operator() () {
            Rational zero(0);
            return Rotated_rational_coefficient
                (Rat_with_sqrt_2_5
                 (Rat_with_sqrt_2(zero,zero,Integer(2)),
                  Rat_with_sqrt_2(zero,zero,Integer(2)),
                  Integer(5)),
                 Rat_with_sqrt_2_5
                 (Rat_with_sqrt_2(zero,zero,Integer(2)),
                  Rat_with_sqrt_2(zero,zero,Integer(2)),
                  Integer(5)),
                 Int_with_sqrt_2_5
                 (Int_with_sqrt_2(Integer(10),Integer(0),Integer(2)),
                  Int_with_sqrt_2(Integer(2),Integer(0),Integer(2)),
                  Integer(5)));
        }
    };        

};


// ZZ[sqrt(3),sqrt(5),sqrt(10+2sqrt(5))]
template <typename Polynomial_2>
struct Rotation_traits_for_base_angle_base<Polynomial_2,6> {
    
    typedef Polynomial_2 Unrotated_polynomial_2;

    typedef typename CGAL::Polynomial_traits_d<Unrotated_polynomial_2>
        ::Innermost_coefficient_type Integer;

    typedef typename CGAL::Get_arithmetic_kernel<Integer>::Arithmetic_kernel
    Arithmetic_kernel;

    typedef typename Arithmetic_kernel::Rational Rational;

private:

    typedef CGAL::Sqrt_extension<Rational,Integer> Rat_with_sqrt_3;
    typedef CGAL::Sqrt_extension<Integer,Integer> Int_with_sqrt_3;
    typedef CGAL::Sqrt_extension<Integer,Integer> Int_with_sqrt_5;
    typedef CGAL::Sqrt_extension<Rat_with_sqrt_3,Integer> Rat_with_sqrt_3_5;
    typedef CGAL::Sqrt_extension<Int_with_sqrt_3,Integer> Int_with_sqrt_3_5;

public:

    typedef CGAL::Sqrt_extension<Rat_with_sqrt_3_5,Int_with_sqrt_3_5> 
        Rotated_rational_coefficient; 
    typedef CGAL::Sqrt_extension<Int_with_sqrt_3_5,Int_with_sqrt_3_5> 
        Rotated_coefficient;

    struct Sin : public std::unary_function<int,Rotated_rational_coefficient> {

        Rotated_rational_coefficient operator() (int angle) {
            
            Angle_coefficients<Rational> coeffs(angle);
            CGAL_assertion_code(
                    typename 
                    CGAL::Algebraic_structure_traits<Rational>::Is_zero zero;
            );
            CGAL_assertion(zero(coeffs.sin.xyzw));
            CGAL_assertion(zero(coeffs.sin.xyz));
            CGAL_assertion(zero(coeffs.sin.xzw));
            CGAL_assertion(zero(coeffs.sin.xyw));
            CGAL_assertion(zero(coeffs.sin.xz));
            CGAL_assertion(zero(coeffs.sin.xw));
            CGAL_assertion(zero(coeffs.sin.xy));
            CGAL_assertion(zero(coeffs.sin.x));
            return Rotated_rational_coefficient
                (Rat_with_sqrt_3_5
                 (Rat_with_sqrt_3(coeffs.sin.c,coeffs.sin.y,Integer(3)),
                  Rat_with_sqrt_3(coeffs.sin.z,coeffs.sin.yz,Integer(3)),
                  Integer(5)),
                 Rat_with_sqrt_3_5
                 (Rat_with_sqrt_3(coeffs.sin.w,coeffs.sin.yw,Integer(3)),
                  Rat_with_sqrt_3(coeffs.sin.zw,coeffs.sin.yzw,Integer(3)),
                  Integer(5)),
                 Int_with_sqrt_3_5
                 (Int_with_sqrt_3(Integer(10),Integer(0),Integer(3)),
                  Int_with_sqrt_3(Integer(2),Integer(0),Integer(3)),
                  Integer(5)));
        }

    };

    struct Cos : public std::unary_function<int,Rotated_rational_coefficient> {

        Rotated_rational_coefficient operator() (int angle) {
            Angle_coefficients<Rational> coeffs(angle);
            CGAL_assertion_code(
                    typename 
                    CGAL::Algebraic_structure_traits<Rational>::Is_zero zero;
            );
            CGAL_assertion(zero(coeffs.cos.xyzw));
            CGAL_assertion(zero(coeffs.cos.xyz));
            CGAL_assertion(zero(coeffs.cos.xzw));
            CGAL_assertion(zero(coeffs.cos.xyw));
            CGAL_assertion(zero(coeffs.cos.xz));
            CGAL_assertion(zero(coeffs.cos.xw));
            CGAL_assertion(zero(coeffs.cos.xy));
            CGAL_assertion(zero(coeffs.cos.x));
            return Rotated_rational_coefficient
                (Rat_with_sqrt_3_5
                 (Rat_with_sqrt_3(coeffs.cos.c,coeffs.cos.y,Integer(3)),
                  Rat_with_sqrt_3(coeffs.cos.z,coeffs.cos.yz,Integer(3)),
                  Integer(5)),
                 Rat_with_sqrt_3_5
                 (Rat_with_sqrt_3(coeffs.cos.w,coeffs.cos.yw,Integer(3)),
                  Rat_with_sqrt_3(coeffs.cos.zw,coeffs.cos.yzw,Integer(3)),
                  Integer(5)),
                 Int_with_sqrt_3_5
                 (Int_with_sqrt_3(Integer(10),Integer(0),Integer(3)),
                  Int_with_sqrt_3(Integer(2),Integer(0),Integer(3)),
                  Integer(5)));
        }

    };
    
    struct Zero {
        
        typedef Rotated_rational_coefficient result_type;

        Rotated_rational_coefficient operator() () {
            Rational zero(0);
            return Rotated_rational_coefficient
                (Rat_with_sqrt_3_5
                 (Rat_with_sqrt_3(zero,zero,Integer(3)),
                  Rat_with_sqrt_3(zero,zero,Integer(3)),
                  Integer(5)),
                 Rat_with_sqrt_3_5
                 (Rat_with_sqrt_3(zero,zero,Integer(3)),
                  Rat_with_sqrt_3(zero,zero,Integer(3)),
                  Integer(5)),
                 Int_with_sqrt_3_5
                 (Int_with_sqrt_3(Integer(10),Integer(0),Integer(3)),
                  Int_with_sqrt_3(Integer(2),Integer(0),Integer(3)),
                  Integer(5)));
        }
    };        

};

// ZZ[sqrt(2),sqrt(3),sqrt(5),sqrt(10+2sqrt(5))]
template <typename Polynomial_2>
struct Rotation_traits_for_base_angle_base<Polynomial_2,3> {
    
    typedef Polynomial_2 Unrotated_polynomial_2;

    typedef typename CGAL::Polynomial_traits_d<Unrotated_polynomial_2>
        ::Innermost_coefficient_type Integer;

    typedef typename CGAL::Get_arithmetic_kernel<Integer>::Arithmetic_kernel
    Arithmetic_kernel;

    typedef typename Arithmetic_kernel::Rational Rational;

private:

    typedef CGAL::Sqrt_extension<Rational,Integer> Rat_with_sqrt_2;
    typedef CGAL::Sqrt_extension<Integer,Integer> Int_with_sqrt_2;
    typedef CGAL::Sqrt_extension<Integer,Integer> Int_with_sqrt_5;
    typedef CGAL::Sqrt_extension<Rat_with_sqrt_2,Integer> Rat_with_sqrt_2_3;
    typedef CGAL::Sqrt_extension<Int_with_sqrt_2,Integer> Int_with_sqrt_2_3;
    typedef CGAL::Sqrt_extension<Rat_with_sqrt_2_3,Integer> 
        Rat_with_sqrt_2_3_5;
    typedef CGAL::Sqrt_extension<Int_with_sqrt_2_3,Integer> 
        Int_with_sqrt_2_3_5;

public:

    typedef CGAL::Sqrt_extension<Rat_with_sqrt_2_3_5,Int_with_sqrt_2_3_5> 
        Rotated_rational_coefficient; 
    typedef CGAL::Sqrt_extension<Int_with_sqrt_2_3_5,Int_with_sqrt_2_3_5> 
        Rotated_coefficient;

    struct Sin : public std::unary_function<int,Rotated_rational_coefficient> {

        Rotated_rational_coefficient operator() (int angle) {
            
            Angle_coefficients<Rational> coeffs(angle);
            return Rotated_rational_coefficient
                (Rat_with_sqrt_2_3_5
                 (Rat_with_sqrt_2_3
                  (Rat_with_sqrt_2(coeffs.sin.c,coeffs.sin.x,Integer(2)),
                   Rat_with_sqrt_2(coeffs.sin.y,coeffs.sin.xy,Integer(2)),
                   Integer(3)),
                  Rat_with_sqrt_2_3
                   (Rat_with_sqrt_2(coeffs.sin.z,coeffs.sin.xz,Integer(2)),
                    Rat_with_sqrt_2(coeffs.sin.yz,coeffs.sin.xyz,Integer(2)),
                   Integer(3)),
                  Integer(5)),
                 Rat_with_sqrt_2_3_5
                  (Rat_with_sqrt_2_3
                   (Rat_with_sqrt_2(coeffs.sin.w,coeffs.sin.xw,Integer(2)),
                    Rat_with_sqrt_2(coeffs.sin.yw,coeffs.sin.xyw,Integer(2)),
                    Integer(3)),
                   Rat_with_sqrt_2_3
                   (Rat_with_sqrt_2(coeffs.sin.zw,coeffs.sin.xzw,Integer(2)),
                    Rat_with_sqrt_2(coeffs.sin.yzw,coeffs.sin.xyzw,Integer(2)),
                    Integer(3)),
                   Integer(5)),
                 Int_with_sqrt_2_3_5
                  (Int_with_sqrt_2_3
                   (Int_with_sqrt_2(Integer(10),Integer(0),Integer(2)),
                    Int_with_sqrt_2(Integer(0),Integer(0),Integer(2)),
                    Integer(3)),
                   Int_with_sqrt_2_3
                   (Int_with_sqrt_2(Integer(2),Integer(0),Integer(2)),
                    Int_with_sqrt_2(Integer(0),Integer(0),Integer(2)),
                    Integer(3)),
                   Integer(5)));

                 
        }

    };

    struct Cos : public std::unary_function<int,Rotated_rational_coefficient> {

        Rotated_rational_coefficient operator() (int angle) {
            
            Angle_coefficients<Rational> coeffs(angle);
            return Rotated_rational_coefficient
                (Rat_with_sqrt_2_3_5
                 (Rat_with_sqrt_2_3
                  (Rat_with_sqrt_2(coeffs.cos.c,coeffs.cos.x,Integer(2)),
                   Rat_with_sqrt_2(coeffs.cos.y,coeffs.cos.xy,Integer(2)),
                   Integer(3)),
                  Rat_with_sqrt_2_3
                   (Rat_with_sqrt_2(coeffs.cos.z,coeffs.cos.xz,Integer(2)),
                    Rat_with_sqrt_2(coeffs.cos.yz,coeffs.cos.xyz,Integer(2)),
                   Integer(3)),
                  Integer(5)),
                 Rat_with_sqrt_2_3_5
                  (Rat_with_sqrt_2_3
                   (Rat_with_sqrt_2(coeffs.cos.w,coeffs.cos.xw,Integer(2)),
                    Rat_with_sqrt_2(coeffs.cos.yw,coeffs.cos.xyw,Integer(2)),
                    Integer(3)),
                   Rat_with_sqrt_2_3
                   (Rat_with_sqrt_2(coeffs.cos.zw,coeffs.cos.xzw,Integer(2)),
                    Rat_with_sqrt_2(coeffs.cos.yzw,coeffs.cos.xyzw,Integer(2)),
                    Integer(3)),
                   Integer(5)),
                 Int_with_sqrt_2_3_5
                  (Int_with_sqrt_2_3
                   (Int_with_sqrt_2(Integer(10),Integer(0),Integer(2)),
                    Int_with_sqrt_2(Integer(0),Integer(0),Integer(2)),
                    Integer(3)),
                   Int_with_sqrt_2_3
                   (Int_with_sqrt_2(Integer(2),Integer(0),Integer(2)),
                    Int_with_sqrt_2(Integer(0),Integer(0),Integer(2)),
                    Integer(3)),
                   Integer(5)));
        }

    };
    
    struct Zero {
        
        typedef Rotated_rational_coefficient result_type;

        Rotated_rational_coefficient operator() () {
            Rational zero(0);
            
            return Rotated_rational_coefficient
                (Rat_with_sqrt_2_3_5
                 (Rat_with_sqrt_2_3
                  (Rat_with_sqrt_2(zero,zero,Integer(2)),
                   Rat_with_sqrt_2(zero,zero,Integer(2)),
                   Integer(3)),
                  Rat_with_sqrt_2_3
                  (Rat_with_sqrt_2(zero,zero,Integer(2)),
                   Rat_with_sqrt_2(zero,zero,Integer(2)),
                   Integer(3)),
                  Integer(5)),
                 Rat_with_sqrt_2_3_5
                  (Rat_with_sqrt_2_3
                   (Rat_with_sqrt_2(zero,zero,Integer(2)),
                    Rat_with_sqrt_2(zero,zero,Integer(2)),
                    Integer(3)),
                   Rat_with_sqrt_2_3
                   (Rat_with_sqrt_2(zero,zero,Integer(2)),
                    Rat_with_sqrt_2(zero,zero,Integer(2)),
                    Integer(3)),
                   Integer(5)),
                 Int_with_sqrt_2_3_5
                  (Int_with_sqrt_2_3
                   (Int_with_sqrt_2(Integer(10),Integer(0),Integer(2)),
                    Int_with_sqrt_2(Integer(0),Integer(0),Integer(2)),
                    Integer(3)),
                   Int_with_sqrt_2_3
                   (Int_with_sqrt_2(Integer(2),Integer(0),Integer(2)),
                    Int_with_sqrt_2(Integer(0),Integer(0),Integer(2)),
                    Integer(3)),
                   Integer(5)));
        }
    };        

};


} // anonymous namespace

template <typename Polynomial_2,int BaseAngle> 
class Rotation_traits_for_base_angle 
    : Rotation_traits_for_base_angle_base<Polynomial_2,BaseAngle> {

public:

    typedef int Angle;
    typedef Polynomial_2 Unrotated_polynomial_2;

    typedef typename Polynomial_traits_d<Unrotated_polynomial_2>
        ::Innermost_coefficient_type Integer;

    typedef Rotation_traits_for_base_angle_base<Integer,BaseAngle>
        Base;

    typedef typename Base::Rotated_coefficient Rotated_coefficient;
    
    typedef typename 
    CGAL::Get_arithmetic_kernel<Integer>::Arithmetic_kernel::Rational Rational;

    typedef typename 
    CGAL::Algebraic_curve_kernel_2_generator<Rotated_coefficient,Rational>
        ::Filtered_algebraic_curve_kernel_with_qir_and_bitstream_2 
        Algebraic_kernel_with_analysis_2;

    typedef typename  Algebraic_kernel_with_analysis_2::Polynomial_2 
        Rotated_polynomial_2;

    struct Rotate : public binary_function<Unrotated_polynomial_2,
                                           Angle,
                                           Rotated_polynomial_2> {

        Rotated_polynomial_2 operator() (Unrotated_polynomial_2 f,
                                         Angle angle) const {

            if(angle%BaseAngle) {
                std::cerr <<  angle << " is not a multiple of " << BaseAngle
                          << std::endl;
                CGAL_error();
            }
            
            angle %= 360; 
            if(angle < 0) {
                angle += 360;
            }

            typedef typename Base::Rotated_rational_coefficient 
                Rotated_rational_coefficient;

            Rotated_rational_coefficient esin = typename Base::Sin()(angle);
            CGAL::simplify(esin);
            Rotated_rational_coefficient ecos = typename Base::Cos()(angle);
            CGAL::simplify(ecos);
            Rotated_rational_coefficient ezero = typename Base::Zero()();
            
            typedef typename CGAL::Polynomial_traits_d<Rotated_polynomial_2>
                ::template Rebind<Rotated_rational_coefficient,1>::Other::Type 
                Rotated_rational_polynomial_1;
            typedef typename CGAL::Polynomial_traits_d<Rotated_polynomial_2>
                ::template Rebind<Rotated_rational_coefficient,2>::Other::Type 
                Rotated_rational_polynomial_2;
            
            Rotated_rational_polynomial_2 
                sub_x(Rotated_rational_polynomial_1(ezero, ecos), 
                      Rotated_rational_polynomial_1(esin)), 
                sub_y(Rotated_rational_polynomial_1(ezero, -esin), 
                      Rotated_rational_polynomial_1(ecos));
            CGAL::simplify(sub_x);
            CGAL::simplify(sub_y);
            std::vector<Rotated_rational_polynomial_2> subs;
            subs.push_back(sub_x);
            subs.push_back(sub_y);
            
            Rotated_rational_polynomial_2 result 
                = typename 
                    CGAL::Polynomial_traits_d<Unrotated_polynomial_2>
                ::Substitute() (f, subs.begin(), subs.end());
            CGAL::simplify(result);
            //std::cout << "rotated poly: " << res << std::endl;
            // integralize polynomial
            typedef CGAL::Fraction_traits<Rotated_rational_polynomial_2> FT;
            typename FT::Denominator_type dummy;
            Rotated_polynomial_2 num;
            typename FT::Decompose()(result, num, dummy);
        
            std::cout << "integralized poly: " << num << "\n\n";
            return num;
        }

    };

};

template <typename RotationTraits>
struct Rotated_algebraic_curve_kernel_2 :
    public RotationTraits::Algebraic_kernel_with_analysis_2 {

public:
    typedef RotationTraits Rotation_traits;
    
    typedef typename Rotation_traits::Algebraic_kernel_with_analysis_2 
        Rotated_kernel_2;

    typedef typename Rotation_traits::Unrotated_polynomial_2
        Poly_int_2;

    typedef typename Rotation_traits::Angle Angle;

    typedef typename Rotation_traits::Rotate Rotate;

public:
    
    //! default constructor
    Rotated_algebraic_curve_kernel_2() :
        Rotated_kernel_2() {
    }

    //! bivariate polynomial over sqrt-exts
    typedef typename Rotated_kernel_2::Polynomial_2 Poly_ext_2;

    typedef typename Rotated_kernel_2::Curve_analysis_2 Curve_analysis_2;
    
    struct Construct_curve_2 {
            
        Curve_analysis_2 operator()(const Poly_int_2& f, 
                                    Angle angle=Angle()) 
            const {
            return Rotated_kernel_2::curve_cache_2()(Rotate()(f, angle));
        }

        Curve_analysis_2 operator()(const Poly_ext_2& f) const {
            return Rotated_kernel_2::curve_cache_2()(f);
        };
    };
    
    Construct_curve_2 construct_curve_2_object() const {
        return Construct_curve_2();
    } 
  
};


CGAL_END_NAMESPACE

#endif //CGAL_ROTATED_ALGEBRAIC_CURVE_KERNEL_2_H
