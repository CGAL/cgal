// Copyright (c) 2004-2008 Max-Planck-Institute Saarbruecken (Germany).
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

/*!\file CGAL/Curved_kernel_via_analysis_2/gfx/Curve_renderer_internals.h
 * \brief
 * contains various low-level routines used by \c Curve_renderer_2 and
 * \c Subdivision_1
 *
 * provide caching for polynomials and polynomial evaluations; 1D range
 * analysis using First Quadratic Affine Form (QF) and Modified Affine
 * Arithmetic (MAA), both equipped with recursive derivative information
 */

#ifndef CGAL_CKVA_CURVE_RENDERER_INTERNALS_H
#define CGAL_CKVA_CURVE_RENDERER_INTERNALS_H 1

#include <vector>
#include <stack>
// #include <boost/multi_index_container.hpp>
// #include <boost/multi_index/member.hpp>
// #include <boost/multi_index/hashed_index.hpp>
// #include <boost/multi_index/sequenced_index.hpp>
#include <boost/functional.hpp>

#include <CGAL/Interval_nt.h>
#include <CGAL/Polynomial/Real_embeddable_traits.h>

#include <CGAL/Curved_kernel_via_analysis_2/gfx/Curve_renderer_traits.h>

// using boost::multi_index::multi_index_container;
// using boost::multi_index::get;
// using boost::multi_index::project;

namespace CGAL {

namespace internal {

#define CGAL_N_CACHES    2    // maximal number of cache instances

#define CGAL_X_RANGE     0       // polynomial in x-range

#define CGAL_Y_RANGE     1       // polynomial in y-range

// defines subdivision level beyond which univariate polynomials are not
// cached (if 0 cache is off)
#define CGAL_MAX_POLY_CACHE_LEVEL   8

// maximal number of entries in univariate polynomial cache containers
#define CGAL_POLY_CACHE_SIZE 1024*1024

// maximal degree of the derivative taken into account during recursive
// derivative range analysis
#define CGAL_RECURSIVE_DER_MAX_DEGREE  7

// 8-pixel neighbouthood directions
static const struct { int x; int y; } directions[] = {
    { 1, 0}, { 1, 1}, { 0, 1}, {-1, 1},
    {-1, 0}, {-1,-1}, { 0,-1}, { 1,-1}};
// map from rectangle sides to pixel directions
// 0 - left side, 1 - right side, 2 - bottom side, 3 - top side
static const int DIR_MAP[][3] =
    {{5, 4, 3}, {7, 0, 1}, {5, 6, 7}, {3, 2, 1}};
// map from 8 pixel directions to "left" (0) or "right" (1) direction of motion
// vertical directions (2 and 6) are mapped to -1
static const int DIR_TAKEN_MAP[] =
    {1, 1, -1, 0, 0, 0, -1, 1};

template <class Integer>
struct Pixel_2_
{
    int x, y;             // pixel coordinates relative to the drawing area
    unsigned level;       // subdivision level: 0 - for pixels,
                          // 1 - for 1/2 pixels, and so on)
#ifdef CGAL_CKVA_RENDER_WITH_REFINEMENT
    double xv, yv; // double approximation of curve-pixel intersection point
#endif

    Integer sub_x, sub_y; // subpixel coordinates relative to pixel's boundary
                          // (always 0 for pixels)

    Pixel_2_& operator =(const Pixel_2_& pix) = default;

    bool operator ==(const Pixel_2_& pix) const {
        return (
          x == pix.x && y == pix.y && level == pix.level &&
          sub_x == pix.sub_x && sub_y == pix.sub_y);
    }
};

// structure describing the seed point and backward direction, support for
// multiple seed points
template <class Integer>
struct Seed_point_
{
    Seed_point_()
    { }

    Seed_point_(const Pixel_2_<Integer>& start_, int dir_,
            int taken_, bool coincide_) : start(start_),
        back_dir(dir_), direction_taken(taken_),
        branches_coincide(coincide_)
    { }

    Pixel_2_<Integer> start; // starting pixel
    int back_dir; // backward direction
    int direction_taken;
    bool branches_coincide;
};

template <class Rational, class AlgebraicReal>
struct Clip_point_ // describes a bottom/top clip point
{
    Clip_point_() { }
    Clip_point_(Rational left_, Rational right_, int arcno_=-1) :
        left(left_), right(right_), arcno(arcno_) { }
    Rational left, right;       // isolating interval boundaries
    int arcno;          // arcno of a segment this point belongs to
    AlgebraicReal alpha;   // this clip point event line
};

template <class Integer>
std::ostream& operator <<(std::ostream& os, const Pixel_2_<Integer>& pix)
{
    os << " (" << pix.x << "/" << pix.sub_x << "; " << pix.y << "/" <<
        pix.sub_y << ") level = " << pix.level;
    return os;
}

/*! \brief defines class \c Curve_renderer_internals
 *
 * provides an interface to low-level range analysis and polynomials evaluation
 * methods along with caching mechanism. \c Coeff_ defines internal polynomial
 * coefficients used by the renderer.
 * \c CurveKernel_2 specifies underlying 2D curve kernel
 */
template <class CurveKernel_2, class Coeff_>
class Curve_renderer_internals
{
public:
    //! \name typedefs
    //!@{

    //! this instance's first template argument
    typedef CurveKernel_2 Curve_kernel_2;

    //! this instance's second template argument
    typedef Coeff_ Coeff;

    //! type of 1-curve analysis
    typedef typename Curve_kernel_2::Curve_analysis_2 Curve_analysis_2;

    //! type of supporting polynomial
    typedef typename Curve_analysis_2::Polynomial_2 Polynomial_2;

    //! underlying univariate polynomial type
    typedef typename Polynomial_2::NT Poly_dst_1;

    //! rational number type
    typedef typename ::CGAL::Get_arithmetic_kernel<
        typename Curve_kernel_2::Coefficient>::Arithmetic_kernel::Rational
            Rational;

    typedef Curve_renderer_traits<Coeff, Rational> Renderer_traits;

    /// polynomial traits should be used whenever possible
    typedef CGAL::Polynomial_traits_d<Polynomial_2> Polynomial_traits_2;

    typedef CGAL::Coercion_traits<Rational,
        typename Curve_kernel_2::Coefficient> Coercion;

    //! coercion between rational and polynomial coefficients number type
    typedef typename Coercion::Type Rat_coercion_type;

    //! base number type for all internal computations
    typedef typename Renderer_traits::Float NT;

    //! instance of a univariate polynomial
    typedef CGAL::Polynomial<Coeff> Poly_1;
    //! instance of a bivariate polynomial
    typedef CGAL::Polynomial<Poly_1> Poly_2;
    //! container's const iterator (random access)
    typedef typename Poly_1::const_iterator const_iterator_1;
    //! container's const iterator (random access)
    typedef typename Poly_2::const_iterator const_iterator_2;

    //! a univariate rational polynomial
    typedef CGAL::Polynomial<Rat_coercion_type> Rational_poly_1;
    //! a bivariate rational polynomial
    typedef CGAL::Polynomial<Rational_poly_1> Rational_poly_2;

    //! conversion from \c Rational type to used number type
    typename Renderer_traits::Rat_to_float rat2float;
    //! conversion from the basic NT to \c Rational
    typename Renderer_traits::Float_to_rat float2rat;
    //! makes the result exact after inexact operation (applicable only for
    //! exact number types
    typename Renderer_traits::Make_exact make_exact;

    //!@}
private:
    //!\name private typedefs
    //!@{

    //! container to store derivative coefficients
    typedef std::vector<NT> Derivative_1;
    //! instance of polynomial derivatives type
    typedef std::vector<CGAL::Polynomial<NT> > Derivative_2;
    //! const_iterator_2 for polynomial derivatives
    typedef typename Derivative_2::const_iterator der_iterator_2;
    //! const_iterator_1 for polynomial derivatives
    typedef typename std::vector<NT>::const_iterator der_iterator_1;

    //! hashed container element's type for univariate polynomials
    typedef std::pair<NT, Poly_1> Poly_entry;

    //! hash function to cache univariate polynomials
//     struct hash_func_poly {
//         typename Renderer_traits::Hash_function hash;
//         std::size_t operator()(const NT& key) const {
//             return hash(key);
//         }
//     };
//
    //! hash function to cache function evaluations
//     struct hash_func_eval {
//         typename Renderer_traits::Hash_function hash;
//         std::size_t operator()(const Eval_hash_key& key) const {
//             std::size_t h1 = hash(key.first), h2 = hash(key.second);
//             return (h1 - h2);
//         }
//     };

    //! hashed map container with LRU capabilities used to store precomputed
    //! univariate polynomials at fixed x or y coordinates
//     typedef boost::multi_index::multi_index_container<
//         Poly_entry,
//         boost::multi_index::indexed_by<
//             boost::multi_index::sequenced<>,
//             boost::multi_index::hashed_unique<
//                 BOOST_MULTI_INDEX_MEMBER(Poly_entry, NT, first),
//                 hash_func_poly > > >
//     Poly_cache;

    //! hashed map used to store precomputed polynomial evaluations
//     typedef boost::multi_index::multi_index_container<
//         Eval_entry,
//         boost::multi_index::indexed_by<
//             boost::multi_index::sequenced<>,
//             boost::multi_index::hashed_unique<
//                 BOOST_MULTI_INDEX_MEMBER(Eval_entry, Eval_hash_key, first),
//                 hash_func_eval > > >
//     Eval_cache;

    //!@}
public:
    //!\name public interface
    //!@{

    //! default constructor
    Curve_renderer_internals() : window(0.0, 0.0, 0.0, 0.0), res_w(0),
         res_h(0) {
    }

    //! \brief checks whether interval contains zero
    inline bool is_zero(const NT& low, const NT& up)
    {
        return (low*up < 0); //(low < 0&&up > 0);
    }

    //! \brief evalutates a certain polynomial derivative at x
    //!
    //! \c der_coeffs is a set of derivative coefficients,
    //! \c poly - polynomial coefficients
    NT evaluate_der(const CGAL::Polynomial<NT>& der_coeffs, const Poly_1& poly,
        const NT& x)
    {
        typename Renderer_traits::Extract_eval extract;
        const_iterator_1 poly_it = poly.end() - 1;
        der_iterator_1 der_it = der_coeffs.end() - 1;
        NT y(extract(*poly_it--) * (*der_it));
        while(der_it != der_coeffs.begin()) {
            --der_it;
            y = y * x + extract(*poly_it--) * (*der_it);
        }
        return y;
    }

    //! \brief the same as \c evaluate but arguments are passed by value
    //! (needed to substitute variables in bivariate polynomial)
    inline static NT binded_eval(Poly_1 poly, NT x)
    {
        return evaluate(poly, x);
    }

    //! \brief evalutates a polynomial at certain x-coordinate
    static NT evaluate(const Poly_1& poly, const NT& x,
        bool *error_bounds_ = nullptr)
    {
        typename Renderer_traits::Extract_eval extract;
        int n = poly.degree()+1, m = (n-1)>>1, odd = n&1;
        if(error_bounds_ != nullptr)
            *error_bounds_ = false;
        if(n == 1)
            return extract(poly.lcoeff(), error_bounds_);
        Coeff cc = static_cast<Coeff>(x);
        const_iterator_1 it1 = poly.end()-1;
        const_iterator_1 it2 = it1 - (n>>1);
        // const_iterator_1 beg = poly.begin()+odd;
        Coeff y1 = *it1, y2 = *it2, mul = cc, y;
        // unrolled loop for better instruction pairing
        while(m-- > odd) {
            it1--; it2--;
            y1 = y1*cc + (*it1);
            y2 = y2*cc + (*it2);
            mul = mul*cc;
        }
        y = y2 + y1*mul;
        if(odd)
            y = poly[0] + y*cc;
        //Gfx_OUT("eval results x = " << x << " res: [" << y.lower() <<
              //  "; " << y.upper() << "\n");

        return extract(y, error_bounds_);
    }

    //! \brief evaluates polynomial at a point (x; y)
    //!
    //! y is the outermost variable, x is the innermost
    static Rat_coercion_type substitute_xy(
        const Rational_poly_2& poly, const Rational& x, const Rational& y) {

        typename Rational_poly_2::const_iterator rit = poly.end()-1;
        Rat_coercion_type r = rit->evaluate(x),
            yc = typename Coercion::Cast()(y);

        while(rit != poly.begin()) {
            --rit;
            r = r * yc + rit->evaluate(x);
        }
        return r;
    }

    //! destructor
    ~Curve_renderer_internals()
    {
        clear_caches();
    }

    //!@}
public:
    //!\name public data (can be accessed from hosting curve renderer)
    //!@{

    CGAL::Bbox_2 window;            //! drawing window
    NT x_min, x_max, y_min, y_max;
    NT pixel_w, pixel_h;            //! pixel dimensions w.r.t. resolution

    Rational x_min_r, y_min_r, x_max_r, y_max_r; //! rational versions
    Rational pixel_w_r, pixel_h_r;

    int res_w, res_h; //! pixel resolution

    Poly_2 *coeffs_x, *coeffs_y; //! f(x(y)) / f(y(x))
    Derivative_2 *der_x, *der_y;  //! derivative coefficients df/dx (df/dy)

    //! caches of precomputed univariate polynomials in x and y directions
    //Poly_cache *cached_x, *cached_y;

    Rational_poly_2 *rational_x, *rational_y;   //! poly with rational coeffs

    Rational_poly_2 *rational_fx, *rational_fy; //! partial derivatives

    //! 0 - the 1st (2nd) derivative does not have zero over an
    //! interval; 1 - does have 0; -1 - not computed
    int first_der, second_der;
    bool sign_change;  //! detected a sign change over interval

    bool zero_bounds;  //! indicates that the result of the last range
                       //! evaluation has at least one zero boundary

   // static bool show_dump;    //! for debugging

    //!@}
private:
    //!\name private members for caching support
    //!@{

    //! data structures grouped into arrays to facilitate caching
    Poly_2 coeffs_x_[CGAL_N_CACHES], coeffs_y_[CGAL_N_CACHES];
    Derivative_2 der_x_[CGAL_N_CACHES], der_y_[CGAL_N_CACHES];
    //Poly_cache cached_x_[CGAL_N_CACHES], cached_y_[CGAL_N_CACHES];
    Rational_poly_2 rational_x_[CGAL_N_CACHES], rational_y_[CGAL_N_CACHES];
    Rational_poly_2 rational_fx_[CGAL_N_CACHES], rational_fy_[CGAL_N_CACHES];

    //!@}
public:
    //! \name public methods
    //!@{

// template <class CurveKernel_2, class Coeff_>
// bool Curve_renderer_internals<CurveKernel_2, Coeff_>::show_dump = false;

//! sets up drawing window and pixel resolution
//!
//! returns \c false if parameters are incorrect
bool setup(const CGAL::Bbox_2& box_, int res_w_, int res_h_)
{
    window = box_;
    x_min = static_cast<NT>(box_.xmin());
    y_min = static_cast<NT>(box_.ymin());
    x_max = static_cast<NT>(box_.xmax());
    y_max = static_cast<NT>(box_.ymax());

    res_w = res_w_;
    res_h = res_h_;

    if (x_min >= x_max) {
      std::cerr << "Incorrect setup parameters: "
                << "x_min = " << x_min << " >= " << x_max << " = x_max"
                << std::endl;
      return false;
    }
    if (y_min >= y_max) {
      std::cerr << "Incorrect setup parameters: "
                << "y_min = " << y_min << " >= " << y_max << " = y_max"
                << std::endl;
      return false;
    }
    if (res_w < 4) {
      std::cerr << "Incorrect setup parameters: "
                << "res_w = " << res_w << " < 4"
                << std::endl;
      return false;
    }
    if (res_h < 4) {
      std::cerr << "Incorrect setup parameters: "
                << "res_h = " << res_h << " < 4"
                << std::endl;
      return false;
    }
     if (res_w > 2048) {
      std::cerr << "Incorrect setup parameters: "
                << "res_w = " << res_w << " > 2048"
                << std::endl;
      return false;
    }
    if (res_h > 2048) {
      std::cerr << "Incorrect setup parameters: "
                << "res_h = " << res_h << " > 2048"
                << std::endl;
      return false;
    }

    pixel_w = (x_max - x_min) / res_w;
    pixel_h = (y_max - y_min) / res_h;

    //srand(time(nullptr));
    // from 0.1 to 0.5
    double rmin = 0.1, rmax = 0.5;
    NT sx = pixel_w * static_cast<NT>(rmin +
            (static_cast<double>(rand()) / RAND_MAX) * (rmax - rmin));
    NT sy = pixel_h * static_cast<NT>(rmin +
            (static_cast<double>(rand()) / RAND_MAX) * (rmax - rmin));

    Gfx_OUT("---------------------- sx: " << (sx / pixel_w) << "; sy: " <<
        (sy / pixel_h) << "\n");
    x_min += sx, x_max += sx, y_min += sy, y_max += sy;

    x_min_r = Rational(x_min), y_min_r = Rational(y_min);
    x_max_r = Rational(x_max), y_max_r = Rational(y_max);

    pixel_w_r = (x_max_r - x_min_r) / res_w;
    pixel_h_r = (y_max_r - y_min_r) / res_h;

    pixel_w = rat2float(pixel_w_r);
    pixel_h = rat2float(pixel_h_r);
    make_exact(pixel_w);
    make_exact(pixel_h);

    //show_dump = false;
    return true;
}

//! \brief computes a range of polynomial values \c f([lower,upper]) using
//! the First Affine Form (AF1) with Recursive Derivative information
bool get_range_AF1_1(int var,
    const NT& l_, const NT& r_, const NT& key, const Poly_1& poly, int check)
{
    Derivative_2 *der = (var == CGAL_X_RANGE ? der_x : der_y);
    NT l(l_), r(r_), l1, r1, low, up;
    NT v1, v2;
    int eval1, eval2;
    der_iterator_2 der_it_2 = der->end()-1;
    der_iterator_1 der_it, der_begin;
    const_iterator_1 cache_it, begin;

    first_der = false, sign_change = false;
    if(poly.degree() == 0) {
        zero_bounds = false;
        return (poly.lcoeff() == NT(0.0));
    }

    if(l_ == r_) {
        zero_bounds = false;
        return false;
    }
    if(l > r) {
        l = r_;
        r = l_;
    }

    eval1 = evaluate_generic(var, l, key, poly);
    eval2 = evaluate_generic(var, r, key, poly);
    bool sign_change = (eval1*eval2 < 0);
    zero_bounds = ((eval1&eval2) == 0);

    if((sign_change||zero_bounds)&&check==1)
        return true;
    if(var == CGAL_X_RANGE) {
        l1 = x_min + l*pixel_w;
        r1 = x_min + r*pixel_w;
    } else {
        l1 = y_min + l*pixel_h;
        r1 = y_min + r*pixel_h;
    }

    typename Renderer_traits::Extract_eval extract;
    unsigned index = CGAL_RECURSIVE_DER_MAX_DEGREE;
    if(index >= der->size()) {
        low = up = extract(poly.lcoeff()) * (*der_it_2).lcoeff();
    } else {
        der_it_2 = der->begin() + index;
        low = 1;
        up = -1;
    }

    NT x0 = (l1 + r1)/2, x1 = (r1 - l1)/2;
    make_exact(x0);
    make_exact(x1);
    NT x0_abs = CGAL_ABS(x0), x1_abs = CGAL_ABS(x1);

    while(der_it_2 != der->begin()) {
        --der_it_2;
        // iterate through derivative coefficients
        der_it = der_it_2->end() - 1;
        der_begin = der_it_2->begin();
        cache_it = poly.end() - 1; // iterate through precomputed y-values
        // if a derivative does not straddle zero we can
        // calculate the exact boundaries for f(x)
        if(low * up >= 0) {
            v1 = v2 = extract(*cache_it--) * (*der_it);
            // calculate the ith derivative at xa and xb
            while(der_it != der_begin) {
                --der_it;

                NT cc1 = extract(*cache_it--) * (*der_it);
                v1 = v1 * l1 + cc1;
                v2 = v2 * r1 + cc1;
            }
            low = v1, up = v2;
        } else { // use affine arithmetic to compute bounds

            NT y0 = extract(*cache_it) * (*der_it), y1(0), e1(0);
            cache_it--;
            while(der_it!=der_begin) {
                --der_it;

                e1 = x0_abs*e1 + x1_abs*(CGAL_ABS(y1) + CGAL_ABS(e1));
                y1 = y1*x0 + x1*y0;
                y0 = x0*y0 + extract(*cache_it) * (*der_it);
                cache_it--;
            }
            NT spread = CGAL_ABS(y1) + CGAL_ABS(e1);
            low = y0 - spread, up = y0 + spread;
        }
        if(der_it_2 == der->begin() && check == 3) {
            first_der = //(eval1*eval2 < 0);
                    is_zero(low, up);
            if(sign_change||zero_bounds)
                return true;
        }
    }
    if(low * up < 0) {
        cache_it = poly.end()-1;
        begin = poly.begin();

        NT y0 = extract(*cache_it), y1(0), e1(0);
        while(cache_it != begin) {
            --cache_it;
            e1 = x0_abs*e1 + x1_abs*(CGAL_ABS(y1) + CGAL_ABS(e1));
            y1 = y1*x0 + x1*y0;
            y0 = x0*y0 + extract(*cache_it);
        }

        NT spread = CGAL_ABS(y1) + CGAL_ABS(e1);
        low = y0 - spread, up = y0 + spread;

        eval1 = CGAL_SGN(low), eval2 = CGAL_SGN(up);
    }
    zero_bounds = ((eval1 & eval2) == 0);
    if(eval1 * eval2 < 0) {
        //std::cerr << "eval1 = " << low << "; eval2 = " << up << std::endl;
        return true;
    }
    return false;
}

//! \brief evaluates a univariate polynomial over an interval using
//! First Quadratic Affine Form
bool get_range_QF_1(int var, const NT& l_, const NT& r_, const NT& key,
    const Poly_1& poly, int check = 1)
{
    Derivative_2 *der = (var == CGAL_X_RANGE ? der_x : der_y);
    NT l(l_), r(r_), l1, r1, low, up;
    NT v1, v2;
    int eval1, eval2;
    der_iterator_2 der_it_2 = der->end()-1;
    der_iterator_1 der_it, der_begin;
    const_iterator_1 cache_it, begin;

    first_der = false, sign_change = false;
    if(poly.degree()==0) {
        zero_bounds = false;
        return (poly.lcoeff()==NT(0.0));
    }
    if(l_ == r_) {
        zero_bounds = false;
        return false;
    }
    if(l > r) {
        l = r_;
        r = l_;
    }

    eval1 = evaluate_generic(var, l, key, poly);
    eval2 = evaluate_generic(var, r, key, poly);
    sign_change = (eval1*eval2 < 0);

    zero_bounds = ((eval1 & eval2) == 0);
    if((sign_change || zero_bounds) && check==1)
        return true;

    if(var == CGAL_X_RANGE) {
        l1 = x_min + l*pixel_w;
        r1 = x_min + r*pixel_w;
    } else {
        l1 = y_min + l*pixel_h;
        r1 = y_min + r*pixel_h;
    }

    typename Renderer_traits::Extract_eval extract;
    unsigned index = CGAL_RECURSIVE_DER_MAX_DEGREE;
    if(index >= der->size()) {
        low = up = extract(poly.lcoeff()) * (*der_it_2).lcoeff();
    } else {
        der_it_2 = der->begin() + index;
        low = 1;
        up = -1;
    }

    NT x0 = (l1 + r1)/2, x1 = (r1 - l1)/2;
    make_exact(x0);
    make_exact(x1);
    NT x0_abs = CGAL_ABS(x0), xsum_abs = CGAL_ABS(x1) + x0_abs;

    while(der_it_2!=der->begin()) {
        --der_it_2;
        // iterate through derivative coefficients
        der_it = (*der_it_2).end()-1;
        der_begin = (*der_it_2).begin();
        cache_it = poly.end()-1; // iterate through precomputed y-values

        // if a derivative does not straddle zero we can obtain exact bounds
        // by evaluating a polynomial at end-points
        if(low * up >= 0) {
            v1 = v2 = extract(*cache_it--)* (*der_it);
            // calculate the ith derivative at xa and xb
            while(der_it != der_begin) {
                --der_it;

                NT cc1 = extract(*cache_it--)* (*der_it);
                v1 = v1 * l1 + cc1;
                v2 = v2 * r1 + cc1;
            }
            low = v1;
            up = v2;
        } else { // use Quadratic Form to compute the bounds

            NT y0 = extract(*cache_it) * (*der_it), y1(0), z1(0), e1(0);
            cache_it--;
            while(der_it != der_begin) {
                --der_it;

                e1 = xsum_abs * e1 + CGAL_ABS(x1 * z1);
                z1 = x0*z1 + x1*y1;
                y1 = y1*x0 + x1*y0;
                y0 = x0*y0 + extract(*cache_it)*(*der_it);
                cache_it--;
            }
            NT spread = CGAL_ABS(y1) + e1;
            low = spread;
            up = spread;
            if(z1 > 0)
                up = up + z1;
            else
                low = low - z1;
            low = y0 - low;
            up = y0 + up;
        }

        if(der_it_2 == der->begin() && check==3) {
            first_der = //(eval1*eval2 < 0);//
                    is_zero(low, up);
            if(sign_change||zero_bounds)
                return true;
        }
    }

    if(low * up < 0) {
        cache_it = poly.end()-1;
        begin = poly.begin();
        NT y0 = extract(*cache_it), y1(0), z1(0), e1(0);
        while(cache_it != begin) {
            --cache_it;
            e1 = xsum_abs * e1 + CGAL_ABS(x1*z1);
            z1 = x0*z1 + x1*y1;
            y1 = y1*x0 + x1*y0;
            y0 = x0*y0 + extract(*cache_it);
        }
        NT spread = CGAL_ABS(y1) + e1;
        low = spread;
        up = spread;
        if(z1 > 0)
            up = up + z1;
        else
            low = low - z1;
        low = y0 - low;
        up = y0 + up;
        eval1 = CGAL_SGN(low);
        eval2 = CGAL_SGN(up);
    }
    //if(eval1*eval2 < 0)
      //  std::cerr << "eval1 = " << low << "; eval2 = " << up << std::endl;
    zero_bounds = ((eval1 & eval2) == 0);
    return (eval1*eval2 < 0);
}

//! \brief evaluates a univariate polynomial over an interval using
//! Modified Affine Arithmetic
bool get_range_MAA_1(int var, const NT& l_, const NT& r_, const NT& key,
    const Poly_1& poly, int check = 1)
{
    Derivative_2 *der = (var == CGAL_X_RANGE) ? der_x : der_y;
    // stores precomputed polynomial derivatives and binominal coeffs
    Derivative_1 der_cache
        //(der->size()+1, NT(0))
        , binom;//(der->size()+1, NT(0));
    NT l(l_), r(r_), l1, r1, low = NT(1), up = NT(-1);
    NT v1, v2, v;
    int eval1, eval2;

    first_der = false, sign_change = false;
    if(poly.degree()==0) {
        zero_bounds = false;
        return (poly.lcoeff()==NT(0.0));
    }
    if(l_ == r_) {
        first_der = false;
        return false;
    }
    if(l > r) {
        l = r_;
        r = l_;
    }
    eval1 = evaluate_generic(var, l, key, poly);
    eval2 = evaluate_generic(var, r, key, poly);

    sign_change = (eval1*eval2 < 0);
    zero_bounds = ((eval1&eval2) == 0);

    if((sign_change || zero_bounds) && check == 1)
        return true;

    if(var == CGAL_X_RANGE) {
        l1 = x_min + l*pixel_w;
        r1 = x_min + r*pixel_w;
    } else {
        l1 = y_min + l*pixel_h;
        r1 = y_min + r*pixel_h;
    }
    NT x0 = (l1 + r1)/2, x1 = (r1 - l1)/2;
    make_exact(x0);
    make_exact(x1);

    int d = 0;
    v1 = evaluate(poly, x0);
    //der_cache[d] = v1;
    der_cache.push_back(v1);
    v = x1, d++;
    der_iterator_2 der_it_2;
    for(der_it_2 = der->begin(); der_it_2 != der->end(); der_it_2++) {
        v1 = evaluate_der((*der_it_2), poly, x0);
        der_cache.push_back(v1);
        //der_cache[d] = v1; // replace by push_backs ?
        //binom[d] = v;
        binom.push_back(v);
        d++;
        //v *= x1/d;
        //make_exact(v);
        v *= x1; // ONLY when derivative coefficients are normalized
    }

    typename Renderer_traits::Extract_eval extract;
    unsigned index = CGAL_RECURSIVE_DER_MAX_DEGREE;
    if(index >= der->size()) {
        low = up = extract(poly.lcoeff()) * (der->end()-1)->lcoeff();
    } else {
        der_it_2 = der->begin()+index;
        low = 1;
        up = -1;
    }
    // assume we have an array of derivatives:
    // der_cache: {f^(0); f^(1); f^(2); ...}
    // and binominal coefficients: [h; h^2/2; h^3/6; ... h^d/d!]
    der_iterator_1 eval_it = der_cache.end()-1, local_it, binom_it,
                   eval_end = der_cache.end();
    d = poly.degree();
    der_it_2 = der->end()-1;
    while(1) {//der_it!=end) {
        if(low * up < 0) {
            v2 = *eval_it;
            local_it = eval_it + 1;
            binom_it = binom.begin();
            low = v2, up = v2;
            while((local_it) != eval_end) {// calculate derivative bounds
                v1 = (*local_it++) * (*binom_it++);
                if(v1 >= 0) { // derivative index is odd
                    up += v1;
                    low -= v1;
                } else {
                    up -= v1;
                    low += v1;
                }
                if(local_it == eval_end)
                    break;
                // derivative index is even
                v1 = (*local_it++) * (*binom_it++);
                (v1 > 0) ? up += v1 : low += v1;
            }
            eval1 = CGAL_SGN(low);
            eval2 = CGAL_SGN(up);

        } else if(d > 0) {
            low = evaluate_der((*der_it_2), poly, l1);
            up = evaluate_der((*der_it_2), poly, r1);
            eval1 = CGAL_SGN(low);
            eval2 = CGAL_SGN(up);

        } else if(d == 0) {
            ;//Gfx_DETAILED_OUT("MAA bounds: sign change\n");
        }

        if(d == 1 && check == 3) {
            first_der = (eval1*eval2 < 0);
            if(sign_change||zero_bounds)
                return true;
        } else if(d == 0) {

            zero_bounds = ((eval1&eval2) == 0);
            return (eval1*eval2 < 0);
        }
        d--; der_it_2--;
        if(eval_it == der_cache.begin())
            break;
        eval_it--;
    }
    return true;
}

//! \brief fixes one coordinate of a bivariate polynomial, uses caching
//! if appropriate
void get_precached_poly(int var, const NT& key, int /* level */, Poly_1& poly)
{
    NT key1;
    //Poly_cache *cached = cached_x;
    Poly_2 *coeffs = coeffs_x;
    //typename boost::multi_index::nth_index_iterator<Poly_cache,1>::type it;
    bool not_cached = true,
            /*(level >= CGAL_MAX_POLY_CACHE_LEVEL)*/ not_found = false;

    if(var == CGAL_Y_RANGE) {
        //cached = cached_y;
        coeffs = coeffs_y;
        key1 = x_min + key*pixel_w;
    } else
        key1 = y_min + key*pixel_h;

//     if(!not_cached) {
//         typename boost::multi_index::nth_index<Poly_cache,1>::type&
//             idx = cached->get<1>();
//         it = idx.find(key1); //*4
//         not_found = (it == idx.end());
//     }

    if(not_cached||not_found) {
        poly = Poly_1(::boost::make_transform_iterator(coeffs->begin(),
                            boost::bind2nd(std::ptr_fun(binded_eval), key1)),
                      ::boost::make_transform_iterator(coeffs->end(),
                            boost::bind2nd(std::ptr_fun(binded_eval), key1)));
        if(not_cached)
            return;
    // all available space consumed: drop the least recently used entry
//         if(cached->size() >= CGAL_POLY_CACHE_SIZE)
//             cached->pop_back();
//         cached->push_front(Poly_entry(key1,poly)); //*4
//         return;
    }
//     cached->relocate(cached->begin(), cached->project<0>(it));
//     poly = (*it).second;
}

//! \brief computes the sign of univariate polynomial at (x; y).
//!
//! if computation with the current precision is not reliable, the sign is
//! recomputed with modular or exact rational arithmetic
int evaluate_generic(int var, const NT& c, const NT& key, const Poly_1& poly)
{
    NT c1;

    if(var == CGAL_X_RANGE)
        c1 = x_min + c*pixel_w;
    else
        c1 = y_min + c*pixel_h;

    bool error_bounds_;
    int sign;
    NT res = evaluate(poly, c1, &error_bounds_);

    if(error_bounds_)
        sign = evaluate_rational(var, c, key);
    else
        sign = CGAL_SGN(res);

    return sign;
}

//! \brief evaluates a polynomial at (x, y) using exact rational arithmetic
int evaluate_rational(int var, const NT& c, const NT& key)
{
    //Rational_poly_1 poly;
    Rational c_r = float2rat(c), key_r = float2rat(key);
    c_r = (var == CGAL_X_RANGE) ? (x_min_r + c_r*pixel_w_r) :
        (y_min_r + c_r*pixel_h_r);

    Rational_poly_2 *rat = rational_x;
    if(var == CGAL_Y_RANGE) {
        rat = rational_y;
        key_r = x_min_r + key_r*pixel_w_r;
    } else
        key_r = y_min_r + key_r*pixel_h_r;

    Rat_coercion_type res = substitute_xy(*rat, key_r, c_r);
    return CGAL_SGN(res);
}

//! precomputes polynomials and derivative coefficients
void precompute(const Polynomial_2& in) {

    typedef typename Polynomial_traits_2::Innermost_coefficient_type Coeff_src;

    Max_coeff<Coeff_src> max_coeff;
    Coeff_src max_c = max_coeff(in);
    /////// magic symbol ////////////
    // somehow relates to double precision fix
    // std::cerr << ' ';

    typedef Reduce_by<Rat_coercion_type, Coeff_src> Reduce_op;
    Transform<Rational_poly_2, Polynomial_2, Reduce_op> transform;
    Reduce_op op(max_c);

    typedef CGAL::Polynomial_traits_d<Rational_poly_2> RP_traits;
    *rational_y = transform(in, op); // divides by maximal coefficient
    *rational_x = typename RP_traits::Swap()(*rational_y, 0, 1);

    // rational fx and fy must have y outermost var and x innermost
    *rational_fx = typename RP_traits::Differentiate()(*rational_y, 0);
    *rational_fy = typename RP_traits::Differentiate()(*rational_y, 1);

    typename Renderer_traits::Convert_poly convert_poly;
    ////////////////////////////////////////////////////////
    /////// ATTENTION: need to call makeExact for bigfloats after conversion
    ////////////////////////////////////////////////////////
    *coeffs_y = convert_poly(*rational_y);
    *coeffs_x = typename CGAL::Polynomial_traits_d<Poly_2>::
            Swap()(*coeffs_y, 0, 1);

    int degree_x = coeffs_x->degree(),
        degree_y = coeffs_y->degree();
/*    cached_x->clear();
    cached_y->clear(); */
    der_x->clear();
    der_y->clear();
    int i, j, maxdeg = (degree_x > degree_y ? degree_x : degree_y);
    std::vector<NT> X(maxdeg, NT(0));
    der_x->reserve(degree_x);
    der_y->reserve(degree_y);

    NT det(1.0);
    for(i = 0; i < degree_x; i++) {
        if(i != 0)
            det = X[0];
        for(j = 1; j <= degree_x - i; j++) {
            // divide by the lowest coefficient ?
            X[j - 1] = (i == 0 ? NT(j) : X[j] * NT(j) / det);
            make_exact(X[j-1]);
        }
        der_x->push_back(CGAL::Polynomial<NT>(X.begin(),
            (X.begin() + degree_x - i)));
    }

    det = NT(1.0);
    for(i = 0; i < degree_y; i++) {
        if(i != 0)
            det = X[0];
        for(j = 1; j <= degree_y - i; j++) {
            // divide by the lowest coefficient ?
            X[j - 1] = (i == 0 ? NT(j) : X[j] * NT(j) / det);
            make_exact(X[j-1]);
        }
        der_y->push_back(CGAL::Polynomial<NT>(X.begin(),
                (X.begin() + degree_y - i)));
    }
}

//! \brief activates the cache entry \c cache_id
void select_cache_entry(std::size_t cache_id)
{
    coeffs_x = coeffs_x_ + cache_id;
    coeffs_y = coeffs_y_ + cache_id;
    rational_x = rational_x_ + cache_id;
    rational_y = rational_y_ + cache_id;
    rational_fx = rational_fx_ + cache_id;
    rational_fy = rational_fy_ + cache_id;
    der_x = der_x_ + cache_id;
    der_y = der_y_ + cache_id;
/*    cached_x = cached_x_ + cache_id;
    cached_y = cached_y_ + cache_id;*/
}

//! \brief empties all cache instances
void clear_caches() {

    for(unsigned i = 0; i < CGAL_N_CACHES; i++) {
/*        cached_x_[i].clear();
        cached_y_[i].clear(); */
        der_x_[i].clear();
        der_y_[i].clear();
    }
}

//!@}

}; // Curve_renderer_internals

} // namespace internal

} //namespace CGAL

#endif // CGAL_CKVA_CURVE_RENDERER_INTERNALS_H
