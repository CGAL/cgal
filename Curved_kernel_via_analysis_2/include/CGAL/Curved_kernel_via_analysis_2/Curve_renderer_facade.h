// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
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
//#define CGAL_CKVA_USE_MULTIPREC_ARITHMETIC
// whether to use exact rational arithmetic
//#define CGAL_CKVA_USE_RATIONAL_ARITHMETIC

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
//#define Gfx_USE_OUT
#ifdef Gfx_USE_OUT
#define Gfx_OUT(x) std::cerr << x
#else
#define Gfx_OUT(x) static_cast<void>(0)
#endif
#endif

#include <CGAL/basic.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Curved_kernel_via_analysis_2/gfx/Curve_renderer_2.h>

CGAL_BEGIN_NAMESPACE

template <class CurvedKernelViaAnalysis_2, class Float_>
class Curve_renderer_interface;

/*!\brief 
 * represents a single curve renderer instance with usual parameter set to
 * speed up rendering of several objects supported by the same curve
 * 
 * @warning not recommended to use for multi-threaded applications
 */
template <class CurvedKernelViaAnalysis_2>
class Curve_renderer_facade
{
    Curve_renderer_facade() { // private constructor
    }

public:
    typedef CGAL::Interval_nt<true> Interval_double;

    typedef Curve_renderer_interface<CurvedKernelViaAnalysis_2,
             Interval_double> Curve_renderer;

    typedef typename Curve_renderer::Coord_vec_2 Coord_vec_2;

    typedef typename Curve_renderer::Coord_2 Coord_2;

    static Curve_renderer& instance() {
        static Curve_renderer _this;
        return _this;
    }

    static void setup(const CGAL::Bbox_2& bbox, int res_w, int res_h) {
        int _w, _h;
        CGAL::Bbox_2 tmp;
        instance().get_resolution(_w, _h);
        instance().get_window(tmp);

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

    //! approximation coordinates
    typedef CGALi::Coord_2 Coord_2;
    
    //! vector of coordinates
    typedef CGALi::Coord_vec_2 Coord_vec_2;
    
    //! exact rational number type
    typedef typename ::CGAL::Get_arithmetic_kernel<
            typename Curve_kernel_2::Coefficient>::Arithmetic_kernel
        Arithmetic_kernel;

    //! multi-precision float NT
    typedef typename Arithmetic_kernel::Bigfloat_interval BFI;
    typedef typename CGAL::Bigfloat_interval_traits<BFI>::Boundary 
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
        static Default_renderer_2 rend;
        return rend;
    }
#ifdef CGAL_CKVA_USE_MULTIPREC_ARITHMETIC
    static Bigfloat_renderer_2& bigfloat_renderer()
    {
        static Bigfloat_renderer_2 rend;
        return rend;
    }
#endif
#ifdef CGAL_CKVA_USE_RATIONAL_ARITHMETIC
    static Exact_renderer_2& exact_renderer()
    {
        static Exact_renderer_2 rend;
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

    inline void get_window(CGAL::Bbox_2& bbox) {
        return renderer().get_window(bbox);
    }
    
    inline void get_resolution(int& res_w, int& res_h) {
        return renderer().get_resolution(res_w, res_h);
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
#endif //  !CGAL_CKVA_DUMMY_RENDERER
    }

    /*!\brief
     * rasterizes an x-monotone curve \c arc
     *
     * returns a list of sequences of pixel coordinates in \c points and
     * end-point coordinats in \c end_points 
     */
    inline void draw(const Arc_2& arc, std::list<Coord_vec_2>& points,
        std::pair<Coord_2, Coord_2>& end_points)
    {
#ifndef CGAL_CKVA_DUMMY_RENDERER
        try {
            renderer().draw(arc, points, end_points);
        } 
        catch(CGALi::Insufficient_rasterize_precision_exception) {
            std::cerr << "Switching to multi-precision arithmetic" << 
                std::endl;
#ifdef CGAL_CKVA_USE_MULTIPREC_ARITHMETIC
            Bbox_2 bbox;
            int res_w, res_h;
            if(::boost::is_same<typename NiX::NT_traits<Float>::Is_exact,
                    CGAL::Tag_true>::value)
                goto Lexit;
            
            get_window(bbox);
            get_resolution(res_w, res_h);
            bigfloat_renderer().setup(bbox, res_w, res_h);
            try {
                points.clear();
                bigfloat_renderer().draw(arc, points, end_points);
                return;
            }
            catch(CGALi::Insufficient_rasterize_precision_exception) {

                std::cerr << "Switching to exact arithmetic" << std::endl;
#ifdef CGAL_CKVA_USE_RATIONAL_ARITHMETIC
                Bbox_2 bbox;
                int res_w, res_h;
                if(::boost::is_same<typename NiX::NT_traits<Float>::Is_exact,
                        CGAL::Tag_true>::value)
                    goto Lexit;

                get_window(bbox);
                get_resolution(res_w, res_h);
                exact_renderer().setup(bbox, res_w, res_h);

                try {
                    points.clear();
                    exact_renderer().draw(arc, points, end_points);
                    return;
                }
                catch(CGALi::Insufficient_rasterize_precision_exception) {
                    goto Lexit;
                }
                return;
#endif  // CGAL_CKVA_USE_RATIONAL_ARITHMETIC
            }
Lexit:  std::cerr << "Sorry, this does not work even with exact "
                    "arithmetic, bailing out..." << std::endl;
        std::cerr << "polynomial: " << renderer().curve().polynomial_2() <<
            std::endl;
        
        renderer().get_window(bbox);
        std::cerr << "window: " << bbox << std::endl;
        
#endif  // CGAL_CKVA_USE_MULTIPREC_ARITHMETIC
        }
#endif // !CGAL_CKVA_DUMMY_RENDERER
    }
    
    /*!\brief
     * rasterizes a point on curve
     */  
    bool draw(const Point_2& point, CGALi::Coord_2& coord) {
#ifndef CGAL_CKVA_DUMMY_RENDERER
        try {
            return renderer().draw(point, coord);
        }    
        catch(CGALi::Insufficient_rasterize_precision_exception) {
            std::cerr << "Unable to rasterize point..\n";
            return false;
        }
#else
        return true;
#endif // !CGAL_CKVA_DUMMY_RENDERER
    }

   //!@}
}; // Curve_renderer_interface



CGAL_END_NAMESPACE

#endif // CGAL_CKVA_CURVE_RENDERER_FACADE_H

