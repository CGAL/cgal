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

/*!\file CGAL/Curved_kernel_via_analysis_2/gfx//Subdivision_1.h
 * \brief definition of Subdivision_1<>
 * 1D space subdivision for rasterization of planar curves
 */

#ifndef CGAL_CKVA_SUBDIVISION_1_H
#define CGAL_CKVA_SUBDIVISION_1_H 1

#include <vector>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
//#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>

#include <CGAL/Polynomial.h>
#include <CGAL/Interval_nt.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>

#include <CGAL/Curved_kernel_via_analysis_2/gfx/Curve_renderer_internals.h>

using boost::multi_index::multi_index_container;
using boost::multi_index::get;
using boost::multi_index::project;

namespace CGAL {

namespace internal {

#ifndef SoX_CURVE_RENDERER_DEFS
#define SoX_CURVE_RENDERER_DEFS

#define SoX_WINDOW_ENLARGE  15 // # of pixels by which the drawing window
                                // is enlarged in y-direction: this is used
                                // to skip closely located clip-points

#define SoX_REFINE_X        1000    // refine factor for X-intervals
                                    // (in pixel size)
#define SoX_REFINE_Y        100000  // refine factor for Y-intervals
                                    // (in pixel size)

#endif // SoX_CURVE_RENDERER_DEFS

#define SoX_REFINE_ISOLATED_POINTS  2 // refine factor for clip points

template <class NT>
struct Range_templ
{
    Range_templ() { }
    Range_templ(const NT& l_, const NT& r_) : left(l_), right(r_)
    {
        sign_change = true;
    }
    NT left, right;
    bool sign_change;
};

} // namespace internal

/*!\brief
 * The class template \c Subdivision_1 and its associate functions.
 *
 * The class implements a space method to plot algebraic curves, we use Affine
 * Arithmetic with recursive derivative information
 */
template <class Coeff_, class Algebraic_curve_2_>
class Subdivision_1
{
public:
    //! \name public typedefs
    //!@{

    //! this instance's first template argument
    typedef Coeff_ Coeff;
    //! this instance's second template argument
    typedef Algebraic_curve_2_ Algebraic_curve_2;

    //!@}
private:
    //! \name private typedefs
    //!@{
    //! rational number type
    typedef typename Algebraic_curve_2::Rational Rational;

    //! special number type traits dependent on polynomial coefficient
    typedef SoX::Curve_renderer_traits<Coeff> Renderer_traits;

    //! specialized integer number type
    typedef typename Renderer_traits::Integer Integer;
    //! supporting bivariate polynomial type
    typedef typename Algebraic_curve_2::Poly_d Poly_dst_2;
    //! supporting univariate polynomial type
    typedef typename Algebraic_curve_2::Poly_d::NT Poly_dst_1;

    //! a univariate rational polynomial
    typedef CGAL::Polynomial<Rational> Rational_poly_1;
    //! a bivariate rational polynomial
    typedef CGAL::Polynomial<Rational_poly_1> Rational_poly_2;

    //! basic number type used in all computations
    typedef typename Renderer_traits::Float NT;
    //! instance of a univariate polynomial
    typedef CGAL::Polynomial<Coeff> Poly_1;
    //! instance of a bivariate polynomial
    typedef CGAL::Polynomial<Poly_1> Poly_2;

    //! conversion from the basic number type to doubles
    typename CGAL::Real_embeddable_traits<NT>::To_double to_double;
    //! conversion from the basic number type to integers
    typename Renderer_traits::To_integer to_integer;
    //! conversion from \c Integer type to built-in integer
    typename Renderer_traits::To_machine_int to_int;
    //! conversion from \c Rational type to used number type
    typename Renderer_traits::From_rational from_rat;
    //! conversion from the basic NT to \c Rational
    typename Renderer_traits::To_rational to_rat;
    //! makes the result exact after inexact operation (applicable only for
    //! exact number types
    typename Renderer_traits::Make_exact make_exact;

    //! returns \c true whenever a precision limit of used number type is
    //! reached
    typename Renderer_traits::Precision_limit limit;
    //! maximum level of subdivision, dependent on used data type
    static const unsigned MAX_SUBDIVISION_LEVEL =
            Renderer_traits::MAX_SUBDIVISION_LEVEL;
    //! isolated point
    typedef Intern::Range_templ<NT> Isolated_point;
    //! a list of isolated points
    typedef std::list<Isolated_point> Isolated_points;

    //! map container element's type for maintaining a list of cache instances
    typedef std::pair<int,int> LRU_entry;

    //! LRU list used for effective cache switching
    typedef boost::multi_index::multi_index_container<
        LRU_entry, boost::multi_index::indexed_by<
            boost::multi_index::sequenced<>,
            boost::multi_index::ordered_unique<
                BOOST_MULTI_INDEX_MEMBER(LRU_entry,int,first) > > >
    LRU_list;

    //!@}
public:
    //! \name Constructors
    //!@{

    //! default constructor
    Subdivision_1() : initialized(false), one(1), polynomial_set(false) {}
    //!@}

public:
    //! \name public methods
    //!@{

    //! sets up drawing window and pixel resolution
    void setup(const double& x_min_,const double& y_min_,
                const double& x_max_,const double& y_max_,
                int res_w_, int res_h_) {
        initialized = engine.setup(x_min_, y_min_, x_max_, y_max_, res_w_,
            res_h_);
    }

    //! sets up the underlying polynomial
    void set_polynomial(const Poly_dst_2& poly)
    {
        input_poly = poly;
        polynomial_set = true;
        select_cache_entry();
    }
    //! returns the underlying polynomial
    Poly_dst_2 get_polynomial() const
    {
        return input_poly;
    }

    //! \brief returns drawing window boundaries
    void get_window(double& x_min_, double& y_min_, double& x_max_,
        double& y_max_) const
    {
        x_min_ = to_double(engine.x_min);
        x_max_ = to_double(engine.x_max);
        y_min_ = to_double(engine.y_min);
        y_max_ = to_double(engine.y_max);
    }

    //! \brief returns pixel resolution
    void get_resolution(int& res_w_, int& res_h_) const
    {
        res_w_ = engine.res_w;
        res_h_ = engine.res_h;
    }

    //! \brief the main rendering procedure, copies a set of pixel coordinates
    //! to the output iterator \c oi
    template <class OutputIterator>
    OutputIterator draw(OutputIterator oi);

    ~Subdivision_1()
    {
        cache_list.clear();
    }

    //!@}
private:
    //! \name Private methods
    //!@{

    //! makes certain cache entry active
    void select_cache_entry();

    //! refines isolated point intervals until pixel size
    void refine_points(int var, const NT& key, const Poly_1& poly,
        Isolated_points& points);

    //! isolates curve points along a sweep-line
    bool isolate_recursive(int var, const NT& beg, const NT& end,
        const NT& clip, const Poly_1& poly, Isolated_points& points);

    //! \brief returns whether a polynomial has zero over an interval,
    //! we are not interested in concrete values
    //!
    //! a bit mask \c check indicates which boundaries are to be computed
    //! 0th bit - of a polynomial itself (0th derivative, default)
    //! 1st bit - of the first derivative (sets \c first_der flag)
    //! 2nd bit - of the second derivative
    bool get_range_1(int var, const NT& lower, const NT& upper, const NT& key,
        const Poly_1& poly, int check = 1);

    //!@}
private:
    //! \name Private properties
    //!@{

    Poly_dst_2 input_poly; //! supporting polynomial

    //! an instance of rendering engine
    Intern::Curve_renderer_internals<Coeff, Algebraic_curve_2> engine;

    int cache_id;        //! index of currently used cache instance
    LRU_list cache_list; //! list of indices of cache instances

    bool initialized;  //! indicates whether the renderer has been initialized
                       //! with correct parameters
    const Integer one; //! just "one"

    bool polynomial_set; //! true if input polynomial was set

    //!@}
}; // class Subdivision_1<>

//! \brief main rasterization procedure, copies in the the output iterator
//! \c oi a set of pixel coordinates

template <class Coeff_, class Algebraic_curve_2_>
template <class OutputIterator>
OutputIterator Subdivision_1<Coeff_, Algebraic_curve_2_>::draw(
    OutputIterator oi)
{
    if(!initialized||!polynomial_set)
        return oi;
    //std::cout << "resolution: " << res_w << " x " << res_h << std::endl;
    //std::cout << "box: [" << x_min << "; " << y_min << "]x[" << x_max << "; "
        // <<   y_max << "]" << std::endl;

    Isolated_points points;
    Poly_1 poly;
    int x, y;

    for(x = 0; x < engine.res_w; x++) {
        NT key = NT(x);
        points.clear();

        engine.get_precached_poly(SoX_Y_RANGE, key, 0, poly);
        isolate_recursive(SoX_Y_RANGE, NT(0), NT(engine.res_h), key, poly,
            points);
        refine_points(SoX_Y_RANGE, key, poly, points);

        typename Isolated_points::iterator it = points.begin();
        while(it != points.end()) {
        //          std::cout << "(" << x << "; " << y << ")  ";
            y = engine.res_h - (int)floor((*it).left);
            //if((*it).sign_change)
                //painter->drawPoint(x, y);
            //else {
            *oi++ = std::make_pair(x, y);
              //int yend = engine.res_h - (int)floor((*it).right);
              //painter->moveTo(x,y);
              //painter->lineTo(x,yend);
            it++;
        }
    }

    for(y = 0; y < engine.res_h; y++) {
        NT key = NT(y);
        points.clear();

        engine.get_precached_poly(SoX_X_RANGE, key, 0, poly);
        isolate_recursive(SoX_X_RANGE, NT(0), NT(engine.res_h), key, poly,
            points);
        refine_points(SoX_X_RANGE, key, poly, points);

        typename Isolated_points::iterator it = points.begin();
        while(it != points.end()) {
        //          std::cout << "(" << x << "; " << y << ")  ";
            x = (int)floor((*it).left);
            *oi++ = std::make_pair(x, engine.res_h - y);

                //int yend = res_h - (int)floor((*it).right);
                //painter->moveTo(x,y);
                //painter->lineTo(x,yend);

            it++;
        }
    }
    //std::cout << "exit normal" << std::endl;
    return oi;
}

template <class Coeff_, class Algebraic_curve_2_>
void Subdivision_1<Coeff_, Algebraic_curve_2_>::refine_points(int var,
        const NT& key, const Poly_1& poly, Isolated_points& points)
{
    NT dist = NT(1)/NT(SoX_REFINE_ISOLATED_POINTS);
    make_exact(dist);
    typename Isolated_points::iterator it = points.begin();
    while(it != points.end()) {
        NT l = (*it).left, r = (*it).right;
        //Gfx_OUT << l << " " << r << "\n";
        int eval_l = engine.evaluate_generic(var, l, key, poly),
            eval_r = engine.evaluate_generic(var, r, key, poly);
        if((eval_l^eval_r)==0) { // no sign change - interval is too small
            /*if(eval_l==0)
                (*it).left = (*it).right = l;
            else if(eval_r==0)
                (*it).left = (*it).right = r;*/
            (*it).sign_change = false;
            //std::cout << "no sign change\n";
            it++;
            continue;
        }
        while(r - l > dist) {
            NT mid = (l+r)/2;
            make_exact(mid);
            int eval_m = engine.evaluate_generic(var, mid, key, poly);
            if(eval_m == 0)
                l = r = mid;
            else if(eval_m == eval_l)
                l = mid;
            else
                r = mid;
        }
        (*it).left = l;
        (*it).right = r;
        it++;
    }
}

template <class Coeff_, class Algebraic_curve_2_>
bool Subdivision_1<Coeff_, Algebraic_curve_2_>::isolate_recursive
    (int var, const NT& beg, const NT& end, const NT& key, const Poly_1& poly,
            Isolated_points& points)
{
    // y_clip: 0 for top boundary, res_h-1 for bottom boundary
    // beg, end: 0 and res_w-1 respectively

    if(!get_range_1(var, beg, end, key, poly, 3))
        return false;

    if(!engine.first_der) {
        //std::cout << "interval found: " << (x_end - x_beg) << std::endl;
        points.push_back(Isolated_point(beg, end));
        return true;
    }

    NT dist = NT(1)/NT(SoX_REFINE_ISOLATED_POINTS);
    make_exact(dist);
    if(end - beg < dist) {// pixel size is reached
        //std::cout << "WARNING: roots are too close..\n";

        // in this case we should generate an exception..
        points.push_back(Isolated_point(beg, end));
        return true;
    }
    NT mid = (beg + end)/2;
    make_exact(mid);

    isolate_recursive(var, beg, mid, key, poly, points);
    isolate_recursive(var, mid, end, key, poly, points);
    return true;
}

//! \brief returns whether a polynomial has zero at a given interval,
//! we are not interested in concrete values
//!
//! if \t der_check is set a range for the first derivative is also computed
template <class Coeff_, class Algebraic_curve_2_>
inline bool Subdivision_1<Coeff_, Algebraic_curve_2_>::get_range_1(int var,
    const NT& lower, const NT& upper, const NT& key, const Poly_1& poly,
        int check)
{
    bool res = engine.get_range_QF_1(var, lower, upper, key, poly, check);
    return res;
}

//! \brief switches to a certain cache instance depending on currently used
//! algebraic curve
template <class Coeff_, class Algebraic_curve_2_>
void Subdivision_1<Coeff_, Algebraic_curve_2_>::select_cache_entry()
{
    cache_id = 0; // no cache is currently used

    engine.select_cache_entry(cache_id);
    engine.precompute(input_poly);
}

} //namespace CGAL

#endif // CGAL_CKVA_SUBDIVISION_1_H
