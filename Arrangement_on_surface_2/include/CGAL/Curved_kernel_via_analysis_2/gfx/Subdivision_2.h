// Copyright (c) 2004-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

/*!\file CGAL/Curved_kernel_via_analysis_2/gfx//Subdivision_2.h
 * \brief definition of Subdivision_2<> 
 * 2D space subdivision for rasterization of planar curves
 */

#ifndef CGAL_CKVA_SUBDIVISION_2_H
#define CGAL_CKVA_SUBDIVISION_2_H 1
#warning this file is considered obsolete

#include <vector>
#include <CGAL/Polynomial.h>
#include <CGAL/Interval_nt.h>

#include <CGAL/Curved_kernel_via_analysis_2/gfx/Curve_renderer_internals.h>

namespace CGAL {

namespace internal {

template <class NT>
struct Range_templ
{
	Range_templ() { }
	Range_templ(const NT& l_, const NT& u_) : lower(l_), upper(u_) { }
	NT lower, upper;
}; 

} // namespace internal

/*!\brief 
 * The class template \c Subdivision_2 and its associate functions.
 * 
 * The class implements a space method to plot algebraic curves, we use Affine
 * Arithmetic with recursive derivative information
 */
template <class NT_, class Algebraic_curve_2_>
class Subdivision_2
{
public: 
	//! \name public typedefs 
    //!@{ 
	//! this instance's first template argument
	typedef NT_ NT;
	//! this instance's second template argument
	typedef Algebraic_curve_2_ Algebraic_curve_2;
	//! specialized integer number type
	typedef typename SoX::Curve_renderer_traits<NT>::Integer Integer;
	//! instance of univariate polynomial
	typedef CGAL::Polynomial<NT> Poly_1;
	//! instance of bivariate polynomial
	typedef CGAL::Polynomial<Poly_1> Poly_2;
	//! rational number type instance
	typedef typename Algebraic_curve_2::Rational Rational;
	//! instance of input bivariate polynomial
	typedef typename Algebraic_curve_2::Poly_d Poly_input_2;
	//! instance of input univariate polynomuial
	typedef typename Poly_input_2::NT Poly_input_1;
	//! instance of coefficient number type
	typedef typename Poly_input_1::NT Coeff;
	//! container used to store coefficient sequence
	typedef typename Poly_1::Vector Vector_1;
	//! container's const iterator (random access)
	typedef typename Poly_1::const_iterator const_iterator_1; 
	//! container used to store univariate polynomials sequence
	typedef typename Poly_2::Vector Vector_2;
	//! container's const iterator (random access)
	typedef typename Poly_2::const_iterator const_iterator_2; 
	//@}
private:	
	//! \name private typedefs 
    //@{ 
		
	//! conversion from the basic number type to integers
	typename SoX::Curve_renderer_traits<NT>::To_integer to_integer;
	//! conversion from \c Integer type to built-in integer
	typename SoX::Curve_renderer_traits<NT>::To_machine_int to_int;
    //! conversion from \c Rational type to used number type
	typename SoX::Curve_renderer_traits<NT>::From_rational from_rat;
	//! makes the result exact after inexact operation (applicable only for
	//! exact number types
	typename SoX::Curve_renderer_traits<NT>::Make_exact make_exact;
	//! an interval
	typedef Intern::Range_templ<NT> Range;
	//! affine form instance
	typedef SoX::Affine_form<NT> Affine_form;
	//@}
public: 
	//! \name Constructors
    //@{ 
	//! default constructor
	Subdivision_2() : initialized(false), polynomial_set(false) {}
	//@}
public:
	//! \name public methods
    //@{ 
	//! specifies drawing window and pixel resolution
	void setup(const double& x_min_,const double& y_min_,
				const double& x_max_,const double& y_max_,
				int res_w_, int res_h_) 
	{ 
		x_min = static_cast<NT>(x_min_); 
		x_max = static_cast<NT>(x_max_); 
		y_min = static_cast<NT>(y_min_); 
		y_max = static_cast<NT>(y_max_);
		res_w = res_w_; 
		res_h = res_h_;
		if(x_min >= x_max||y_min >= y_max||res_w < 5||res_h < 5||res_w > 1024||
			res_h > 1024) {
			std::cout << "Incorrect setup parameters" << std::endl;
			initialized = false;
			return;
		}
		pixel_w = (x_max - x_min) / res_w;
		pixel_h = (y_max - y_min) / res_h;
		make_exact(pixel_w);
		make_exact(pixel_h);
		initialized = true;	
	}    
	//! sets up the underlying polynomial
	void set_polynomial(const Poly_input_2& poly)
	{
		input_poly = poly;
		precompute();
	}
	//! returns the underlying polynomial
	Poly_input_2 get_polynomial() const
	{
		return input_poly;
	}
	//! \brief returns drawing window boundaries
	void get_window(double& x_min_, double& y_min_, double& x_max_, 
		double& y_max_) const
	{
		x_min_ = to_double(x_min); 
		x_max_ = to_double(x_max); 
		y_min_ = to_double(y_min); 
		y_max_ = to_double(y_max);
	}
	//! \brief returns pixel resolution
	void get_resolution(int& res_w_, int& res_h_) const
	{
		res_w_ = res_w; 
		res_h_ = res_h;
	}
	//! \brief the main rendering procedure
	void draw(QPainter *painter_); 
	//@}	
private:
	//! \name Private methods
    //@{ 
	//! precomputes polynomials and derivative coefficients
	void precompute();
	//! \brief switches to another cache instance depending on the
	//! supporting curve of a segment
		//! \brief evalutates the ith derivative at certain x
	//!
	//! \c cache_it - an intetator pointing to the end of an array of 
	//! polynomial coefficients, \c der_it - an iterator for derivative 
	//! coefficients
	NT evaluate_der(const_iterator_1 der_it, const_iterator_1 begin,
			const_iterator_1 cache_it, const NT& x)
	{
		NT val((*cache_it--) * (*der_it));
		while((der_it--)!=begin) 
			val = val * x + (*cache_it--) * (*der_it);
		return val;
	}
	//! evalutates a function at a certain x
	NT evaluate(const Poly_1& poly, const NT& x)
	{
		const_iterator_1 it = poly.end() - 1, begin = poly.begin();
		NT val(*it);
		while((it--)!=begin) 
			val = val * x + (*it);
		return val;
	}
	//! \brief computes the ranges of univariate polynomial values over an interval
	void get_range_1(int var, const NT& lower, const NT& upper, 
		const Poly_1& poly, NT& l, NT& h);
	//! Affine Arithmetic with Recursive Derivative information 
	//! for univariate case
	void get_range_AARD_1(int var, const NT& lower, const NT& upper, 
		const Poly_1& poly, NT& l, NT& h);
	//! \brief Recursive Taylor, bivariate case
	//!
	//! returns a range of polynomial values as Affine_form
	void get_range_RT_2(const NT& x_low,  const NT& x_high, const NT& y_low,
		const NT& y_high, int depth, int index, Affine_form& res);
	//! checks a rectangular area with 2D range analysis: either discrads it,
	//! subdivides further or draws a pixel
	void quad_tree(const NT& x_low, const NT& x_high, const NT& y_low, 
		const NT& y_high);
	//! \brief recursive quad tree subdivision
	void subdivide(const NT& x_low, const NT& x_high, const NT& y_low, 
		const NT& y_high);
public:
	//! destructor
	~Subdivision_2()
	{
	
	}
	//@}
private:
	//! \name Private properties
    //@{ 
	NT x_min, x_max, y_min, y_max; //! drawing window boundaries
	int res_w, res_h;			   //! pixel resolution
	NT pixel_w, pixel_h; 		   //! pixel dimensions w.r.t. resolution
	Poly_2 coeffs_x, coeffs_y; //! f(x(y)) / f(y(x))
	Poly_input_2 input_poly;
	Vector_2 der_x, der_y;  //! derivative coefficients df/dx (df/dy)
	std::vector<Poly_2> mixed_derivatives;
	int max_deg; //! the maximal degree w.r.t. x and y
	bool initialized, polynomial_set;	
	QPainter *painter;
	//@} 
}; // class Subdivision_2<>

//! \brief main rasterization procedure
template <class NT, class Algebraic_curve_2_>
void Subdivision_2<NT, Algebraic_curve_2_>::draw(QPainter *painter_)
{
	if(!initialized||!polynomial_set||painter_==NULL)
		return;
	painter = painter_;
	//std::cout << " P(x(y)): " << coeffs_x << std::endl; 
	//std::cout << " P(y(x)): " << coeffs_y << std::endl; 
	std::cout << "resolution: " << res_w << " x " << res_h << std::endl;
	//std::cout << "box: [" << x_min << "; " << y_min << "]x[" << x_max << "; "
		// <<	y_max << "]" << std::endl;
	
	//quad_tree(x_min/15, x_max/15, y_min/15, y_max/15);
	quad_tree(x_min, x_max, y_min, y_max);
	//std::cout << "exit normal" << std::endl;
}

//! \brief checks a rectangular area with 2D range analysis: either discrads it,
//! subdivides further or draws a pixel
template <class NT_, class Algebraic_curve_2_>
void Subdivision_2<NT_, Algebraic_curve_2_>::quad_tree(const NT& x_low, 
	  const NT& x_high, const NT& y_low, const NT& y_high)
{
	Affine_form res;
	NT lower, upper;
	get_range_RT_2(x_low, x_high, y_low, y_high, 0, 0, res);
	res.convert(lower, upper);
	if(lower*upper > 0)
		return;
	if(x_high - x_low <= pixel_w&&y_high - y_low <= pixel_h) {
		NT x = (x_low), y = (y_low);
		int pix_x = static_cast<int>(CGAL::to_double((x - x_min) / pixel_w)),
		    pix_y = static_cast<int>(CGAL::to_double((y - y_min) / pixel_h));
		painter->drawPoint(pix_x, res_h - pix_y);
		//painter->drawEllipse(pix_x-2,res_h-pix_y-2,4,4);
	}
	else 
		subdivide(x_low, x_high, y_low, y_high);
}

//! \brief recursive quad tree subdivision
template <class NT_, class Algebraic_curve_2_>
void Subdivision_2<NT_, Algebraic_curve_2_>::subdivide(const NT& x_low, 
	  const NT& x_high, const NT& y_low, const NT& y_high)
{
	NT x_mid = (x_low + x_high)/2, y_mid = (y_low + y_high)/2;
	quad_tree(x_low,x_mid,y_low,y_mid);
	quad_tree(x_low,x_mid,y_mid,y_high);
	quad_tree(x_mid,x_high,y_mid,y_high);
	quad_tree(x_mid,x_high,y_low,y_mid);
}

//! \brief Recursive Taylor, bivariate case
template <class NT_, class Algebraic_curve_2_>
void Subdivision_2<NT_, Algebraic_curve_2_>::get_range_RT_2(
	const NT& x_low, const NT& x_high, const NT& y_low, const NT& y_high,
		int depth, int index, Affine_form& res)
{
	//std::cout << "range for [" << x_low << "; " << y_low << "]x[" << 
	//x_high <<
		//"; " << y_high << "]: (" << low << "; " << high << ")" << std::endl;
	typename std::vector<Poly_2>::const_iterator der_it =
		mixed_derivatives.begin() + index;
	if((*der_it).degree()==0) {
		NT c = (*der_it).lcoeff().lcoeff();
		res = Affine_form(c,c);
		return;
	} else if((*der_it).degree()==1) {
		Poly_1 p_x_low = NiX::substitute_x((*der_it), x_low),
			   p_x_high = NiX::substitute_x((*der_it), x_high);
		NT eval_1 = p_x_low.evaluate(y_low), 
		   eval_2 = p_x_low.evaluate(y_high),
		   eval_3 = p_x_high.evaluate(y_low), 
		   eval_4 = p_x_high.evaluate(y_high);
		NT min1 = eval_1, max1 = eval_2;
		if(min1 > eval_2) {
			min1 = eval_2;
			max1 = eval_1;
		}
		NT min2 = eval_3, max2 = eval_4;
		if(min2 > eval_4) {
			min2 = eval_4;
			max2 = eval_3;
		}
		if(min1 > min2)
			min1 = min2;
		if(max1 < max2)
			max1 = max2;
		res = Affine_form(min1,max1);
		return;
	}
	
	if(depth >= 4)
	{
	NT l, h;
	std::vector<Affine_form> x_forms;
	const_iterator_2 it_2 = (*der_it).begin(); 
	while(it_2 != (*der_it).end())
	{
		// it_2 - a poly in x_range
		if((*it_2).is_zero())
			l = h = 0;
		else
			get_range_1(X_RANGE, x_low, x_high, *it_2, l, h); 
		x_forms.push_back(Affine_form(l, h));
		it_2++;
	}
	typename std::vector<Affine_form>::const_iterator it = 
		x_forms.end()-1;
	Affine_form y_form(y_low, y_high);
	res = Affine_form(*it);
	while((it--) != x_forms.begin())
		res = res * y_form + (*it);
	//res.convert(lower, upper);
	return;
	}
	
	NT x0 = (x_low+x_high)/2, y0 = (y_low+y_high)/2, x1 = (x_high-x_low)/2,
		y1 = (y_high-y_low)/2;
	Affine_form one1(-1,1), one2(-1,1), one3(-1,1),
			zero1(0,1), zero2(0,1), fxx, fxy, fyy;
	NT eval_f = NiX::substitute_xy(*der_it, x0, y0),
	   eval_fx = NiX::substitute_xy(*(der_it+depth+1), x0, y0)*x1,
	   eval_fy = NiX::substitute_xy(*(der_it+depth+2), x0, y0)*y1;
	res = eval_f + eval_fx*one1 + eval_fy*one2;
	   
	int idx = index+2*depth+3;
	if(idx < (int)mixed_derivatives.size()) {
		get_range_RT_2(x_low, x_high, y_low, y_high, depth+2, idx,   fxx);
		get_range_RT_2(x_low, x_high, y_low, y_high, depth+2, idx+1, fxy);
		get_range_RT_2(x_low, x_high, y_low, y_high, depth+2, idx+2, fyy);
		res = res + ((x1*x1/2*zero1)*fxx) + ((y1*y1/2*zero2)*fyy) +
			((x1*y1*one3)*fxy);
	} 
		
	//res.convert(lower, upper);
	//std::cout << "range for depth = " << depth << " index = " << index << 
		//" [" << lower << "; " << upper << "]" << std::endl;
}

//! \brief computes the range of polynomial values \c f([lower,upper]) using
//! Affine Arithmetic with Recursive Derivative information
//!
//! \c var = \c X_RANGE: \c y is fixed - polynomial is given in \c x 
//! coordinates, otherwise: \c x is fixed - polynomial is given in \c y 
//! coordinates 
template <class NT_, class Algebraic_curve_2_>
void Subdivision_2<NT_, Algebraic_curve_2_>::get_range_AARD_1(int var, 
	const NT& l_, const NT& r_, const Poly_1& poly,  NT& l1, NT& h1)
{
	Vector_2 *der = &der_y;
	if(var == X_RANGE) 
		der = &der_x; 
	NT low, up, l(l_), r(r_);
	const_iterator_2 der_it_2 = der->end()-1; 
	const_iterator_1 der_it, cache_it, begin;
	if(poly.degree()==0) {
		l1 = h1 = poly.lcoeff();
		return;
	}
	low = up = poly.lcoeff() * (*der_it_2).lcoeff();
	if(l > r) {
		l = r_;
		r = l_;
	}
	Affine_form f(l,r);
	/*while((der_it_2--)!=der->begin()) {
		// iterate through derivative coefficients
		der_it = (*der_it_2).end()-1; 
		begin = (*der_it_2).begin();
		cache_it = poly.end()-1; // iterate through precomputed y-values
		// if a derivative does not straddle zero we can 
		// calculate the exact boundaries for f(x)		
		if(low * up > 0) {
			v1 = v2 = (*cache_it--) * (*der_it);
			// calculate the ith derivative at xa and xb
			while((der_it--)!=begin) {
				v1 = v1 * l + (*cache_it) * (*der_it);
				v2 = v2 * r + (*cache_it--) * (*der_it);
			}
			if(low < 0 && up < 0) { 
				v = v1; 
				v1 = v2; 
				v2 = v; 
			}
			low = v1; 
			up = v2;
		} else { // use affine arithmetic to compute bounds
			Affine_form<NT>	range(((*cache_it--) * (*der_it)));
			// calculate the ith derivative using affine arithmetic	
			while((der_it--)!=begin) 
				range = (*cache_it--) * (*der_it) + (range * f);
			range.convert(low,up);
		}
	}
	if(low * up >= 0) {
		Gfx_OUT << "+ " << std::endl;
		v1 = evaluate(poly, l);
		v2 = evaluate(poly, r); 
		//Gfx_OUT << "v1 = " << v1 << " v2 = " << v2 << std::endl;
		if(low < 0 && up < 0) { 
			v = v1; 
			v1 = v2; 
			v2 = v; 
		}
		low = v1; 
		up = v2;
	} else */{ // use affine arithmetic to compute bounds
		cache_it = poly.end()-1;
		begin = poly.begin();
		Affine_form res(*cache_it);
		while((cache_it--)!=begin) 
			res = (*cache_it) + (res * f);
		res.convert(low,up);
	}
	l1 = low;
	h1 = up;
//std::cout << "AARD bounds: [" << low << "; " << up << "]" << std::endl;	
}

//! \brief returns whether a polynomial has zero at a given interval,
//! since we are not interested in concrete values
//!
//! flag \t der_check forces to calculate boundaries for the first 
//! derivative instead of function itself. Caching is used if
//! appropriate
template <class NT_, class Algebraic_curve_2_>
void Subdivision_2<NT_, Algebraic_curve_2_>::get_range_1(int var, 
	const NT& lower, const NT& upper, const Poly_1& poly, NT& l, NT& h)
{
	//std::cout << "range for: [" << lower << "; " << upper << "] poly: " <<
		//poly << std::endl;
	get_range_AARD_1(var, lower, upper, poly, l, h);
}

//! precomputes some initial data, necessary for subsequent computations
template <class NT_, class Algebraic_curve_2_>
void Subdivision_2<NT_, Algebraic_curve_2_>::precompute()
{
	Intern::Max_coeff<Coeff> max_coeff;
	Coeff mx1 = max_coeff(input_poly);
	Intern::Transform<Poly_2, Coeff, Rational, 
		typename SoX::Curve_renderer_traits<NT>::From_rational,
		typename SoX::Curve_renderer_traits<NT>::Make_exact> tr;
	tr.max = Rational(mx1);
	coeffs_y = tr(input_poly);
	// x - outer variable, y - inner
	coeffs_x = transpose_bivariate_polynomial(coeffs_y); 
	int degree_x = coeffs_x.degree(),
		degree_y = coeffs_y.degree();
	der_x.clear();
	der_y.clear();
	int i, j; 
	max_deg = degree_x;
	if(degree_y > max_deg) 
		max_deg = degree_y;
	NT *X = new NT[max_deg];
	NT det(1.0);
	std::cout << "start" << std::endl;
	for(i = 0; i < degree_x; i++) {
		if(i != 0) 
			det = X[0];
		for(j = 1; j <= degree_x - i; j++)
			if(i == 0) 
				X[j-1] = j;
			else { 
				X[j-1] = X[j] * j / det; // divide by the lowest coefficient ?
				make_exact(X[j-1]);
			}	
		der_x.push_back(Poly_1(X,(X + degree_x - i)));
	}
	for(i = 0; i < degree_y; i++) {
		if(i != 0) 
			det = X[0];
		for(j = 1; j <= degree_y - i; j++)
			if(i == 0) 
				X[j-1] = j; // divide by the lowest coefficient ?
			else { 
				X[j-1] = X[j] * j / det; 
				make_exact(X[j-1]);
			}
		der_y.push_back(Poly_1(X,(X + degree_y - i)));
	}
	delete []X;
	mixed_derivatives.clear();
	// f, fx, fy, fxx, fxy, fyy, fxxx, fxxy, fxyy, fyyy, ...
	mixed_derivatives.push_back(coeffs_y);
	int idx = 0;
	for(i = 1; i <= max_deg; i++)
	{
		mixed_derivatives.push_back(NiX::diff_x(mixed_derivatives[idx]));
		for(j = 0; j < i; j++) 
// compute fx^(i-j)y^(j), i.e. (i-j) times derivate by x; and j times
// derivate by y
		{
			Poly_2 p = mixed_derivatives[idx+j];
			p.diff();
			mixed_derivatives.push_back(p);
		}
		idx += i;
	}
	std::cout << "finished" << std::endl;
	polynomial_set = true;
	/*typename std::vector<Poly_2>::const_iterator der_it = 
		mixed_derivatives.end()-1;
	for(i = max_deg; i >= 0; i--)
	{
		std::cout << i << "th mixed derivatives: " << std::endl;
		for(j = 0; j < i+1; j++, der_it--)
			std::cout << *der_it << std::endl;
	}*/
}

} //namespace CGAL

#endif // CGAL_CKVA_SUBDIVISION_2_H
