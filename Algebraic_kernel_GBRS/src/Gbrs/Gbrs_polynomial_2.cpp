// Copyright (c) 2006 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
//
// Author(s)     : Luis Peñaranda <penarand@loria.fr>

#include <CGAL/assertions.h>
#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gbrs_algebraic_1.h>
#include <CGAL/Gbrs_polynomial_2.h>
#include <iostream>
#include <gmp.h>

CGAL_BEGIN_NAMESPACE

// constructors
Rational_polynomial_2::Rational_polynomial_2 () {
	degree_x = 0;
	degree_y = 0;
	coef = (mpq_t**)malloc (sizeof(mpq_t*));
	coef[0] = (mpq_t*)malloc (sizeof(mpq_t));
	mpq_init (coef[0][0]);
};

Rational_polynomial_2::Rational_polynomial_2 (unsigned int dx,
		unsigned int dy) {
	degree_x = (int)dx;
	degree_y = (int)dy;
	coef = (mpq_t**)malloc (sizeof(mpq_t*)*(degree_x+1));
	for (int i=0; i<degree_x+1; ++i) {
		coef[i] = (mpq_t*)malloc (sizeof(mpq_t)*(degree_y+1));
		for (int j=0; j<degree_y+1; ++j)
			mpq_init (coef[i][j]);
	}
};

Rational_polynomial_2::Rational_polynomial_2 (int dx, int dy) {
	degree_x = dx<0?0:dx;
	degree_y = dy<0?0:dy;
	coef = (mpq_t**)malloc (sizeof(mpq_t*)*(degree_x+1));
	for (int i=0; i<degree_x+1; ++i) {
		coef[i] = (mpq_t*)malloc (sizeof(mpq_t)*(degree_y+1));
		for (int j=0; j<degree_y+1; ++j)
			mpq_init (coef[i][j]);
	}
};

Rational_polynomial_2::Rational_polynomial_2 (const Rational_polynomial_2 &p) {
	degree_x = p.degree_x;
	degree_y = p.degree_y;
	mpq_t **p_coef = p.get_coefs ();
	coef = (mpq_t**)malloc (sizeof(mpq_t*)*(degree_x+1));
	// we have to copy the contents, not just the pointer
	for (int i=0; i<degree_x+1; ++i) {
		coef[i] = (mpq_t*)malloc (sizeof(mpq_t)*(degree_y+1));
		for (int j=0; j<degree_y+1; ++j) {
			mpq_init (coef[i][j]);
			mpq_set (coef[i][j], p_coef[i][j]);
		}
	}
};

// destructor
Rational_polynomial_2::~Rational_polynomial_2 () {
	for (int i=0; i<degree_x+1; ++i) {
		for (int j=0; j<degree_y+1; ++j)
			mpq_clear (coef[i][j]);
		free (coef[i]);
	}
	free (coef);
};

Rational_polynomial_2&
Rational_polynomial_2::operator= (const Rational_polynomial_2 &p) {
	// destroy the current data
	for (int i=0; i<degree_x+1; ++i) {
		for (int j=0; j<degree_y+1; ++j)
			mpq_clear (coef[i][j]);
		free (coef[i]);
	}
	free (coef);
	// copy data from p
	degree_x = p.degree_x;
	degree_y = p.degree_y;
	mpq_t **p_coef = p.get_coefs ();
	coef = (mpq_t**)malloc (sizeof(mpq_t*)*(degree_x+1));
	// we have to copy the contents, not just the pointer
	for (int i=0; i<degree_x+1; ++i) {
		coef[i] = (mpq_t*)malloc (sizeof(mpq_t)*(degree_y+1));
		for (int j=0; j<degree_y+1; ++j) {
			mpq_init (coef[i][j]);
			mpq_set (coef[i][j], p_coef[i][j]);
		}
	}
	return *this;
};

Rational_polynomial_2 Rational_polynomial_2::operator-()const{
	Rational_polynomial_2 opposite(degree_x,degree_y);
	mpq_t **coef_o=opposite.get_coefs();
	for(int i=0;i<degree_x+1;++i)
		for(int j=0;j<degree_y+1;++j)
			mpq_neg(coef_o[i][j],coef[i][j]);
	return opposite;
};

Rational_polynomial_2& Rational_polynomial_2::operator+=(const Rational_polynomial_2 &s){
	int sx,sy,minorx,majorx,minory,majory,i,j;
	bool am_i_bigger_x,am_i_bigger_y,am_i_bigger,am_i_smaller;
	mpq_t **coef_s=s.get_coefs();
	if(am_i_bigger_x=((sx=s.get_degree_x())<degree_x)){
		minorx=sx;
		majorx=degree_x;
	}else{
		minorx=degree_x;
		majorx=sx;
	}
	if(am_i_bigger_y=((sy=s.get_degree_y())<degree_y)){
		minory=sy;
		majory=degree_y;
	}else{
		minory=degree_y;
		majory=sy;
	}
	Rational_polynomial_2 sum(majorx,majory);
	mpq_t **coef_sum=(mpq_t**)malloc(sizeof(mpq_t*)*(majorx+1));
	for(i=0;i<minorx+1;++i){
		coef_sum[i]=(mpq_t*)malloc(sizeof(mpq_t)*(majory+1));
		for(j=0;j<minory+1;++j){
			mpq_init(coef_sum[i][j]);
			mpq_add(coef_sum[i][j],coef_s[i][j],coef[i][j]);
			mpq_clear(coef[i][j]);
		}
		for(j=minory+1;j<majory+1;++j){
			mpq_init(coef_sum[i][j]);
			if(am_i_bigger_y){
				mpq_set(coef_sum[i][j],coef[i][j]);
				mpq_clear(coef[i][j]);
			}else
				mpq_set(coef_sum[i][j],coef_s[i][j]);
		}
		free(coef[i]);
	}
	am_i_bigger=am_i_bigger_x&&am_i_bigger_y;
	am_i_smaller=(!am_i_bigger_x)&&(!am_i_bigger_y);
	for(i=minorx+1;i<majorx+1;++i){
		coef_sum[i]=(mpq_t*)malloc(sizeof(mpq_t)*(majory+1));
		for(j=0;j<minory+1;++j){
			mpq_init(coef_sum[i][j]);
			if(am_i_bigger_x){
				mpq_set(coef_sum[i][j],coef[i][j]);
				mpq_clear(coef[i][j]);
			}else
				mpq_set(coef_sum[i][j],coef_s[i][j]);
		}
		if(am_i_bigger||am_i_smaller)
			for(j=minory+1;j<majory+1;++j){
				mpq_init(coef_sum[i][j]);
				if(am_i_bigger){
					mpq_set(coef_sum[i][j],coef[i][j]);
					mpq_clear(coef[i][j]);
				}else
					mpq_set(coef_sum[i][j],coef_s[i][j]);
			}
		if(am_i_bigger_x)
			free(coef[i]);
	}
	free(coef);
	coef=coef_sum;
	degree_x=majorx;
	degree_y=majory;
	return *this;
};

// XXX: is it a good idea to implement this using operator+=?
Rational_polynomial_2
Rational_polynomial_2::operator+(const Rational_polynomial_2 &s)const{
	int sx=s.get_degree_x();
	int sy=s.get_degree_y();
	int minorx,majorx,minory,majory;
	bool am_i_bigger_x,am_i_bigger_y, am_i_bigger,am_i_smaller;
	mpq_t **coef_s=s.get_coefs();
	if(sx<degree_x){
		am_i_bigger_x=true;
		minorx=sx;
		majorx=degree_x;
	}else{
		am_i_bigger_x=false;
		minorx=degree_x;
		majorx=sx;
	}
	if(sy<degree_y){
		am_i_bigger_y=true;
		minory=sy;
		majory=degree_y;
	}else{
		am_i_bigger_y=false;
		minory=degree_y;
		majory=sy;
	}
	am_i_bigger=(am_i_bigger_x&&am_i_bigger_y);
	am_i_smaller=((!am_i_bigger_x)&&(!am_i_bigger_y));
	Rational_polynomial_2 sum(majorx,majory);
	mpq_t **coef_sum=sum.get_coefs();
	int i,j;
	for(i=0;i<minorx+1;++i){
		for(j=0;j<minory+1;++j)
			mpq_add(coef_sum[i][j],coef_s[i][j],coef[i][j]);
		for(j=minory+1;j<majory+1;++j)
			mpq_set(coef_sum[i][j],(am_i_bigger_y?
						coef[i][j]:
						coef_s[i][j]));
	}
	for(i=minorx+1;i<majorx+1;++i){
		for(j=0;j<minory+1;++j)
			mpq_set(coef_sum[i][j],(am_i_bigger_x?
						coef[i][j]:
						coef_s[i][j]));
		if (am_i_bigger||am_i_smaller)
			for(j=minory+1;j<majory+1;++j)
				mpq_set(coef_sum[i][j],(am_i_bigger?
							coef[i][j]:
							coef_s[i][j]));
	}
	return sum;
};

// multiplies the polynomial by scale * x^shift_x * y^shift_y
// (preconditions: shift_[xy] >= 0)
/*Rational_polynomial_2& Rational_polynomial_2::scale_and_shift(mpq_t &scale,
		int shift_x,int shift_y){
	int i,j;
	degree_x+=shift_x;
	degree_y+=shift_y;
	mpq_t **new_coef=(mpq_t**)malloc(sizeof(mpq_t*)*(degree_x+1));
	for(i=0;i<degree_x+1;++i){
		new_coef[i]=(mpq_t*)malloc(sizeof(mpq_t)*(degree_y+1));
		for(j=0;j<degree_y+1;++j)
			mpq_init(new_coef[i][j]);
	}
	for(i=shift_x;i<degree_x+1;++i){
		for (j=shift_y; j<degree_y+1; ++j) {
			mpq_mul(new_coef[i][j],coef[i-shift_x][j-shift_y],scale);
			mpq_clear(coef[i-shift_x][j-shift_y]);
		}
		free(coef[i-shift_x]);
	}
	free(coef);
	coef=new_coef;
	return *this;
};*/

// TODO: implement a better multiplication algorithm
Rational_polynomial_2 Rational_polynomial_2::operator*
(const Rational_polynomial_2 &f)const{
	/*// how to multiply:
	// 1. create a polynomial to store the result (with zeros)
	// 2. multiply *this by every monomial of &f (with shift)
	// 3. sum all
	Rational_polynomial_2 product;
	mpq_t **f_coefs=f.get_coefs();
	int xf=f.get_degree_x();
	int yf=f.get_degree_y();
	Rational_polynomial_2 partial;
	for (int i=0; i<xf+1; ++i)
		for (int j=0; j<yf+1; ++j)
			product+=(partial=*this).scale_and_shift(f_coefs[i][j],i,j);
	return product;*/
	// TODO: finish testing this (if it doesn't work, use the infraoptimal
	// c++-flavored implementation above)
	mpq_t **coef_f=f.get_coefs();
	int x_f=f.get_degree_x();
	int y_f=f.get_degree_y();
	int x_p=degree_x+x_f;
	int y_p=degree_y+y_f;
	Rational_polynomial_2 product(x_p,y_p);
	mpq_t **coef_p=product.get_coefs();
	mpq_t temp;
	mpq_init(temp);
	for(int coef_i=0;coef_i<x_p+1;++coef_i)
	  for(int coef_j=0;coef_j<y_p+1;++coef_j){
	    int max_i=(coef_i<degree_x?coef_i:degree_x)+1;
	    int max_j=(coef_j<degree_y?coef_j:degree_y)+1;
	    for(int i=0;i<max_i;++i)
	      for(int j=0;j<max_j;++j)
		if((coef_i-i<=x_f)&&(coef_j-j<=y_f)){
		  mpq_mul(temp,coef[i][j],coef_f[coef_i-i][coef_j-j]);
		  mpq_add(coef_p[coef_i][coef_j],coef_p[coef_i][coef_j],temp);
		}
	  }
	mpq_clear(temp);
	return product;
};

// TODO: test this
bool Rational_polynomial_2::operator==(const Rational_polynomial_2 &p)const{
	mpq_t **coef_p=p.get_coefs();
	int p_x=p.get_degree_x();
	int p_y=p.get_degree_y();
	int i,j,minorx,majorx,minory,majory;
	bool am_i_bigger_x, am_i_bigger_y;
	if(degree_x<p_x){
		am_i_bigger_x=false;
		majorx=p_x;
		minorx=degree_x;
	}else{
		am_i_bigger_x=true;
		majorx=degree_x;
		minorx=p_x;
	}
	if(degree_y<p_y){
		am_i_bigger_y=false;
		majory=p_y;
		minory=degree_y;
	}else{
		am_i_bigger_y=true;
		majory=degree_y;
		minory=p_y;
	}
	bool am_i_bigger=am_i_bigger_x&&am_i_bigger_y;
	bool am_i_smaller=(!am_i_bigger_x)&&(!am_i_bigger_y);
	for(i=0;i<minorx+1;++i){
		for(j=0;j<minory+1;++j)
			if(mpq_cmp(coef[i][j],coef_p[i][j]))
				return false;
		for(j=minory+1;j<majory+1;++j)
			if(mpq_sgn(am_i_bigger_y?coef[i][j]:coef_p[i][j]))
				return false;
	}
	for(i=minorx+1;i<majorx+1;++i){
		for(j=0;j<minory+1;++j)
			if(mpq_sgn(am_i_bigger_x?coef[i][j]:coef_p[i][j]))
				return false;
		if(am_i_bigger)
			for(j=minory+1;j<majory+1;++j)
				if(mpq_sgn(coef[i][j]))
					return false;
		else
			if(am_i_smaller)
				for(j=minory+1;j<majory+1;++j)
					if(mpq_sgn(coef_p[i][j]))
						return false;
	}
	return true;
};

std::ostream& Rational_polynomial_2::show (std::ostream &s) const {
	bool printed = false;
	for (int i=degree_x;i>-1;--i)
		for (int j=degree_y;j>-1;--j) {
			if(int sgn=mpq_sgn(coef[i][j])){
				if((sgn==1)&&printed)
					s<<"+";
				if (mpq_cmp_ui(coef[i][j],1,1)) {
					s << coef[i][j];
					if (i)
						s<<"*";
				}
				printed = true;
				if (i) {
					s<<"x";
					if (i>1)
						s<<"^"<<i;
				}
				if (j) {
					if (i)
						s<<"*";
					s<<"y";
					if (j>1)
						s<<"^"<<j;
				}
			}
		}
#ifdef CGAL_RS_DEBUG
	if (!printed)
		s<<"0";
	s<<" [ dx="<<degree_x<<" dy="<<degree_y<<" ";
	int i,j;
	for (i=0; i<degree_x+1; ++i) {
		s<<"[ ";
		for (j=0; j<degree_y+1; ++j)
			s<<coef[i][j]<<" ";
		s<<"] ";
	}
	s<<"]";
	return s;
#else
	return (printed?s:(s<<"0"));
#endif
};

Rational_polynomial_2& Rational_polynomial_2::operator*=
(const Rational_polynomial_2 &f){
	Rational_polynomial_2 aux(*this);
	return (*this=aux*f);
};

Rational_polynomial_2& Rational_polynomial_2::operator*=(const mpq_t &s){
	for(int i=0;i<degree_x+1;++i)
		for(int j=0;j<degree_y+1;++j)
			mpq_mul(coef[i][j],coef[i][j],s);
	return *this;
};

Rational_polynomial_2& Rational_polynomial_2::operator*=(const CGAL::Gmpq &s){
	for(int i=0;i<degree_x+1;++i)
		for(int j=0;j<degree_y+1;++j)
			mpq_mul(coef[i][j],coef[i][j],s.mpq());
	return *this;
};

CGAL_END_NAMESPACE
