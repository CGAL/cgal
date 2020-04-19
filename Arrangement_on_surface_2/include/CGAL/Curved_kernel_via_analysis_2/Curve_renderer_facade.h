// Copyright (c) 2008,2009,2010,2011 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

/*!\file CGAL/Curved_kernel_via_analysis_2/Curve_renderer_facade.h
 * \brief
 * definition of \c Curve_renderer_facade<>
 *
 * high-level interface to the curve renderer
 */

#ifndef CGAL_CKVA_CURVE_RENDERER_FACADE_H
#define CGAL_CKVA_CURVE_RENDERER_FACADE_H

// do not compile curve renderer code (for fast debugging)
//#define CGAL_CKVA_DUMMY_RENDERER

// whether to use multi-precision arithmetic
#define CGAL_CKVA_USE_MULTIPREC_ARITHMETIC
// whether to use exact rational arithmetic
#define CGAL_CKVA_USE_RATIONAL_ARITHMETIC

// this turns on a signleton curve renderer
// (not recommended for multi-threaded applications)
//#define CGAL_CKVA_USE_STATIC_RENDERER

// prints out detailed info of the rendering process
// (only for internal debugging)
#ifndef Gfx_DETAILED_OUT
//#define Gfx_USE_DETAILED_OUT
#ifdef Gfx_USE_DETAILED_OUT
#define Gfx_DETAILED_OUT(x) std::cerr << x
#else
#define Gfx_DETAILED_OUT(x) static_cast<void>(0)
#endif
#endif

// prints out only high-level debug messages
#ifndef Gfx_OUT
// #define Gfx_USE_OUT
#ifdef Gfx_USE_OUT
#define Gfx_OUT(x) std::cerr << x
#else
#define Gfx_OUT(x) static_cast<void>(0)
#endif
#endif

// I still breathe !!!
#define STILL_ALIVE std::cout << __LINE__ << "\n";

#include <CGAL/basic.h>
#include <CGAL/tss.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Arithmetic_kernel.h>

#include <boost/array.hpp>
#include <CGAL/Curved_kernel_via_analysis_2/gfx/Curve_renderer_2.h>


namespace CGAL {

template <class CurvedKernelViaAnalysis_2, class Float_>
class Curve_renderer_interface;

/*!\brief
 * represents a single curve renderer instance with usual parameter set to
 * speed up rendering of several objects supported by the same curve
 *
 * @warning not recommended to use in multi-threaded applications
 */
template <class CurvedKernelViaAnalysis_2>
class Curve_renderer_facade
{
    Curve_renderer_facade() { // private constructor

    }

public:
    typedef CGAL::Interval_nt<true> Interval_double;

    typedef Curve_renderer_interface<CurvedKernelViaAnalysis_2,
             Interval_double > Curve_renderer;

    static Curve_renderer& instance() {
        CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(Curve_renderer, _this);
        return _this;
    }

    static void setup(const CGAL::Bbox_2& bbox, int res_w, int res_h) {
        int _w, _h;
        CGAL::Bbox_2 tmp;
//         CORE::CORE_init(2);
//     CORE::setDefaultPrecision(70, CORE::extLong::getPosInfty());
        instance().get_setup_parameters(&tmp, _w, _h);
        if(bbox != tmp || res_w != _w || res_h != _h) {
            instance().setup(bbox, res_w, res_h);
        }
    }
};

/*! \brief
 * CKvA interface to the curve renderer. Provides three levels of increasing
 * arithmetic precision
 *
 * the curve renderer is instantiated with the base \c Float_ which can be
 * either integral or user-defined floating-point number type, and optionally
 * with multi-precision \c Bigfloat and exact rational number types as
 * provided by the currently used \c Arithmetic_kernel
 */
template <class CurvedKernelViaAnalysis_2, class Float_>
class Curve_renderer_interface
{
public:
    //! \name Public typedefs
    //!@{

    //! this instance's first argument
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! base floating-point number type
    typedef Float_ Float;

    //! type of a curve kernel
    typedef typename Curved_kernel_via_analysis_2::Curve_kernel_2
        Curve_kernel_2;

    //! exact rational number type
    typedef typename ::CGAL::Get_arithmetic_kernel<
            typename Curve_kernel_2::Coefficient>::Arithmetic_kernel
        Arithmetic_kernel;

    //! multi-precision float NT
    typedef typename Arithmetic_kernel::Bigfloat_interval BFI;
    typedef typename CGAL::Bigfloat_interval_traits<BFI>::Bound
        Bigfloat;

    //! rational NT
    typedef typename Arithmetic_kernel::Rational Rational;

    //! an arc of generic curve
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Arc_2;

    //! a point on curve
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;

#ifndef  CGAL_CKVA_DUMMY_RENDERER
    //! instance of a curve renderer
    typedef Curve_renderer_2<Curved_kernel_via_analysis_2, Float>
        Default_renderer_2;
#ifdef CGAL_CKVA_USE_MULTIPREC_ARITHMETIC
    //! the curve renderer instantiated with BigFloat
    typedef Curve_renderer_2<Curved_kernel_via_analysis_2, Bigfloat>
        Bigfloat_renderer_2;
#endif
#ifdef CGAL_CKVA_USE_RATIONAL_ARITHMETIC
    //! the curve renderer instantiated with Rationals
    typedef Curve_renderer_2<Curved_kernel_via_analysis_2, Rational>
        Exact_renderer_2;
#endif

#ifdef CGAL_CKVA_USE_STATIC_RENDERER
    static Default_renderer_2& renderer()
    {
        CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(Default_renderer_2, rend);
        return rend;
    }
#ifdef CGAL_CKVA_USE_MULTIPREC_ARITHMETIC
    static Bigfloat_renderer_2& bigfloat_renderer()
    {
        CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(Bigfloat_renderer_2, rend);
        return rend;
    }
#endif
#ifdef CGAL_CKVA_USE_RATIONAL_ARITHMETIC
    static Exact_renderer_2& exact_renderer()
    {
        CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(Exact_renderer_2, rend);
        return rend;
    }
#endif
#else // !CGAL_CKVA_USE_STATIC_RENDERER
    Default_renderer_2 rend;
    Default_renderer_2& renderer() {
        return rend;
    }
#ifdef CGAL_CKVA_USE_MULTIPREC_ARITHMETIC
    Bigfloat_renderer_2 bigfloat_rend;
    Bigfloat_renderer_2& bigfloat_renderer() {
        return bigfloat_rend;
    }
#endif
#ifdef CGAL_CKVA_USE_RATIONAL_ARITHMETIC
    Exact_renderer_2 exact_rend;
    Exact_renderer_2& exact_renderer() {
        return exact_rend;
    }
#endif
#endif // !CGAL_CKVA_USE_STATIC_RENDERER

#endif // !CGAL_CKVA_DUMMY_RENDERER

    //!@}
public:
    //! \name Public methods and properties
    //!@{

    //!@note pass null-pointer for \c pbox parameter if you don't need
    //! the drawing window
    inline void get_setup_parameters(CGAL::Bbox_2 *pbox, int& res_w,
                int& res_h) {
#ifndef CGAL_CKVA_DUMMY_RENDERER
        return renderer().get_setup_parameters(pbox, res_w, res_h);
#endif
    }

    /*!\brief
     * initializes renderer with drawing window dimensions
     * and pixel resolution
     *
     * <tt>[x_min; y_min]x[x_max; y_max]</tt> - drawing window
     * \c res_w and \c res_h - h/v pixel resolution
     */
    inline void setup(const Bbox_2& bbox, int res_w, int res_h)
    {
#ifndef CGAL_CKVA_DUMMY_RENDERER
        renderer().setup(bbox, res_w, res_h);
#endif
    }

    /*!\brief
     * rasterizes an x-monotone curve \c arc
     *
     * inserts the list of sequences of pixel coordinates as objects of type
     * \c Coord_2 to \c pts , \c Container must support \c push_back and \c
     * clear operations
     *
     * \c Coord_2 must be constructible from a pair of integers / doubles
     * depending on the renderer type
     *
     * computes optionaly end-point coordinates (even if they lie outside the
     * window)
     */
    template < class Coord_2, template < class, class > class Container,
        class Allocator >
    inline void draw(const Arc_2& arc,
            Container< std::vector< Coord_2 >, Allocator >& pts,
            boost::optional< Coord_2 > *end_pt1 = nullptr,
            boost::optional< Coord_2 > *end_pt2 = nullptr) {

#ifndef CGAL_CKVA_DUMMY_RENDERER
        Bbox_2 bbox;
#ifdef CGAL_CKVA_USE_MULTIPREC_ARITHMETIC
        int res_w, res_h;
#endif

        renderer().set_IA_method(0); // start with QF
Lrestart:
        try {
            renderer().draw(arc, pts, end_pt1, end_pt2);
        }
        catch(internal::Insufficient_rasterize_precision_exception) {

            int prev = renderer().set_IA_method(1);
            if(prev == 0) {
                std::cerr << "Restarting with MAA\n";
                pts.clear();
                goto Lrestart;
            }

            std::cerr << "Switching to multi-precision arithmetic" <<
                std::endl;
#ifdef CGAL_CKVA_USE_MULTIPREC_ARITHMETIC
            if(::boost::is_same<typename Algebraic_structure_traits< Float >::
                    Is_exact, CGAL::Tag_true>::value)
                goto Lexit;

            get_setup_parameters(&bbox, res_w, res_h);
            exact_renderer().setup(bbox, res_w, res_h);
            exact_renderer().set_IA_method(1); // always use MAA
            try {
                pts.clear();
                exact_renderer().draw(arc, pts, end_pt1, end_pt2);
                return;
            }
            catch(internal::Insufficient_rasterize_precision_exception) {

            // HACK HACK HACK: uses exact arithmetic
                goto Lexit;
                std::cerr << "Switching to exact arithmetic" << std::endl;
#ifdef CGAL_CKVA_USE_RATIONAL_ARITHMETIC

                if(::boost::is_same<
                    typename Algebraic_structure_traits< Float >::Is_exact,
                         CGAL::Tag_true>::value)
                    goto Lexit;

                get_setup_parameters(&bbox, res_w, res_h);
                exact_renderer().setup(bbox, res_w, res_h);
                exact_renderer().set_IA_method(1); // always use MAA
                try {
                    pts.clear();
                    exact_renderer().draw(arc, pts, end_pt1, end_pt2);
                    return;
                }
                catch(internal::Insufficient_rasterize_precision_exception) {
                    goto Lexit;
                }
                return;
#endif  // CGAL_CKVA_USE_RATIONAL_ARITHMETIC
            }
Lexit:  std::cerr << "Sorry, this does not work even with exact "
                    "arithmetic, bailing out..." << std::endl;
        std::cerr << "polynomial: " << renderer().curve().polynomial_2() <<
            std::endl;

        renderer().get_setup_parameters(&bbox, res_w, res_h);
        std::cerr << "window: " << bbox << "; resolution: " <<
            res_w << " x " << res_h << std::endl;

#endif  // CGAL_CKVA_USE_MULTIPREC_ARITHMETIC
        }
#endif // !CGAL_CKVA_DUMMY_RENDERER
    }

    /*!\brief
     * rasterizes a point on curve, returns point coordinates as objects of
     * type \c Coord_2 which are constructible from a pair of ints / doubles
     *
     * retunrs \c false if point lies outside the window or cannot be
     * rasterized due to precision problems
     */
    template < class Coord_2 >
    bool draw(const Point_2& point, Coord_2& coord) {
#ifndef CGAL_CKVA_DUMMY_RENDERER
        try {
            return renderer().draw(point, coord);
        }
        catch(internal::Insufficient_rasterize_precision_exception) {
            std::cerr << "Unable to rasterize point..\n";
            return false;
        }
#else
        return true;
#endif // !CGAL_CKVA_DUMMY_RENDERER
    }

   //!@}
}; // Curve_renderer_interface

} //namespace CGAL

#endif // CGAL_CKVA_CURVE_RENDERER_FACADE_H

