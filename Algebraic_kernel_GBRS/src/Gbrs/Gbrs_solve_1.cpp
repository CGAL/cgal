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

// this file is based on the RS original example, rs_demo_mp.c (because we
// don't know RS, this little funny suprising amazing fantastic candy-flavored
// library)

#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>
#include <CGAL/assertions.h>
#include <CGAL/Gbrs_algebraic_1.h>
#include <CGAL/Gbrs_polynomial_1.h>
#include <CGAL/Gbrs_solve_1.h>
#include <rs_exports.h>

CGAL_BEGIN_NAMESPACE

int init_solver(){
	static unsigned int rs_inits=0;
	if (!rs_inits) {
		rs_init_rs();
		++rs_inits;
	}
	return 0;
}

int affiche_sols_eqs(mpfi_ptr *x){
	int ident_sols_eqs,nb_elts,ident_node,ident_vect,ident_elt,i;
	ident_sols_eqs=rs_get_default_sols_eqs ();
	nb_elts=rs_export_list_vect_ibfr_nb(ident_sols_eqs);
	ident_node=rs_export_list_vect_ibfr_firstnode(ident_sols_eqs);
	//x=(mpfi_ptr*)malloc(nb_elts*sizeof(mpfi_ptr));
	for (i=0; i<nb_elts; ++i) {
		ident_vect=rs_export_list_vect_ibfr_monnode(ident_node);
		CGAL_assertion_msg(rs_export_dim_vect_ibfr(ident_vect)==1,
				"the dimension of vector must be 1");
		ident_elt=rs_export_elt_vect_ibfr(ident_vect,0);
		x[i]=(mpfi_ptr)rs_export_ibfr_mpfi(ident_elt);
		ident_node=rs_export_list_vect_ibfr_nextnode(ident_node);
	}
	return nb_elts;
}

void affiche_sols_constr(int nr,mpfi_ptr p){
	int ident_sols_eqs,nb_elts,ident_node,ident_vect,nb,ident_elt;
	ident_sols_eqs=rs_get_default_sols_ineqs();
	nb_elts=rs_export_list_vect_ibfr_nb(ident_sols_eqs);
	ident_node=rs_export_list_vect_ibfr_firstnode(ident_sols_eqs);
	for(int i=0;i<nb_elts;++i){
		ident_vect=rs_export_list_vect_ibfr_monnode(ident_node);
		nb=rs_export_dim_vect_ibfr(ident_vect);
		CGAL_assertion_msg((nb==1),
				"the vector must contain one element");
		ident_elt=rs_export_elt_vect_ibfr(ident_vect,0);
		if(i==nr){
			mpfi_set(p,(mpfi_ptr)rs_export_ibfr_mpfi(ident_elt));
			//break;
		}
		ident_node=rs_export_list_vect_ibfr_nextnode(ident_node);
	}
}

Sign affiche_signs_constr(const Algebraic_1 &a){
	int ident_sols_eqs,nb_elts,ident_node,ident_vect, nb, ident_elt;
	mpfi_t tmp;
	mpfi_init(tmp);
	ident_sols_eqs = rs_get_default_sols_ineqs ();
	nb_elts = rs_export_list_vect_ibfr_nb (ident_sols_eqs);
	ident_node = rs_export_list_vect_ibfr_firstnode (ident_sols_eqs);
	for (int i=1; i<nb_elts+1; ++i) {
		ident_vect = rs_export_list_vect_ibfr_monnode (ident_node);
		nb = rs_export_dim_vect_ibfr (ident_vect);
		CGAL_assertion_msg ((nb == 1),
				"the vector must contain one element");
		ident_elt = rs_export_elt_vect_ibfr (ident_vect, 0);
		if (i == a.nr()+1) {
			mpfi_set(tmp,(mpfi_ptr)rs_export_ibfr_mpfi(ident_elt));
			break;
		}
		ident_node = rs_export_list_vect_ibfr_nextnode (ident_node);
	}
	/*std::cerr << "\nreturned value: ";
	mpfi_out_str(stderr,10,0,tmp);
	std::cerr << std::endl;*/
	/* mpfi_is_zero(tmp) doesn't work. The reason is that MPFR_SIGN in
	   the mpfi code returns 1 when applied to the left and right zeros.
	   This is not surprising because zero is signed in IEEE 754, and MPFR
	   adopts it. Nevertheless, mpfr_sgn returns 0, but mpfi doesn't use
	   it to implement mpfi_is_zero.
	   Here is the difference (from MPFR source code):
	    define mpfr_sgn(_x)      (mpfr_zero_p(_x) ? 0 : MPFR_SIGN(_x))
	*/
	if(mpfr_zero_p(&(tmp->right))&&mpfr_zero_p(&(tmp->left)))
		return ZERO;
	// the same holds for mpfi_is_pos and mpfi_is_neg
	if((mpfr_sgn(&(tmp->left))>=0)&&(mpfr_sgn(&(tmp->right)))>0)
		return POSITIVE;
	if((mpfr_sgn(&(tmp->left))<0)&&(mpfr_sgn(&(tmp->right))<=0))
		return NEGATIVE;
	// if we arrive here, it is because the signs of the endpoints are -
	// and +, and (I think) RS guarantees that this never happens
	CGAL_assertion_msg(false,"error in sign calculation");
	return ZERO;
}

void create_rs_upoly (mpz_t *poly, const int deg, const int ident_pol) {
	int ident_mon, ident_coeff, i;
	rs_import_uppring ("T");	// what the heck is this?
	for (i=0; i<=deg; ++i)
		if (mpz_sgn (poly[i]))	{	// don't add if == 0
			// (this is one of the few things we know about RS)
			ident_mon = rs_export_new_mon_upp_bz ();
			ident_coeff = rs_export_new_gmp ();
			rs_import_bz_gmp
				(ident_coeff, TO_RSPTR_IN (&(poly[i])));
			rs_dset_mon_upp_bz (ident_mon, ident_coeff, i);
			rs_dappend_list_mon_upp_bz (ident_pol, ident_mon);
		}
}

void create_rs_uconstr (mpz_t ** list_constr,
		const int * list_degs, const int ident_list) {
	int ident_poly;
	ident_poly = rs_export_new_list_mon_upp_bz ();
	create_rs_upoly (*list_constr, *list_degs, ident_poly);
	rs_dappend_list_sup_bz (ident_list, ident_poly);
}

int solve_1(mpfi_ptr *x,Rational_polynomial_1 &p1,unsigned int prec){
	rs_reset_all();
	create_rs_upoly(p1.get_coefs(),p1.get_degree(),rs_get_default_up());
	set_rs_precisol(prec);
	set_rs_verbose(CGAL_RS_VERB);
	rs_run_algo("UISOLE");
	return affiche_sols_eqs(x);
}

// y=p1(a)
void eval_1(const Rational_polynomial_1 &p1,const Algebraic_1 &a,mpfi_ptr y){
// I think it's not a good idea to call RS to do this simple evaluation
#if 0
	CGAL_assertion((a.is_consistent())/*&&(a.nr()>=0)*/);
	mpz_t **constr;
	int *degs;
	rs_reset_all();
	create_rs_upoly(a.pol().get_coefs(),a.pol().get_degree(),
			rs_get_default_up());
	constr=(mpz_t**)malloc(sizeof(mpz_t*));
	*constr=p1.get_coefs();
	degs=(int*)malloc(sizeof(int));
	*degs=p1.get_degree();
	create_rs_uconstr(constr,degs,rs_get_default_ineqs_u());
	set_rs_precisol(a.rsprec());
	set_rs_verbose(CGAL_RS_VERB);
	rs_run_algo("UISOLES");
	affiche_sols_constr(a.nr(),y);
#endif
	p1.eval_mpfi(y,a.mpfi());
}

Sign sign_1_rs
(const Rational_polynomial_1 &p1,const Algebraic_1 &a,unsigned int prec){
	mpz_t **constr;
	int *degs;
	// XXX: is this always necessary? can I do it just once?
	rs_reset_all ();
	// tell RS to find the roots of this polynomial
	create_rs_upoly (a.pol().get_coefs (), a.pol().get_degree (),
			rs_get_default_up ());
	// the constraint vector will have just one element
	constr = (mpz_t**)malloc(sizeof(mpz_t*));
	*constr = p1.get_coefs ();
	degs = (int*)malloc(sizeof(int));
	*degs = p1.get_degree ();
	create_rs_uconstr (constr, degs, rs_get_default_ineqs_u ());
	set_rs_precisol (prec);
	set_rs_verbose (CGAL_RS_VERB);
	rs_run_algo ("UISOLES");
	//Sign s=affiche_signs_constr(a);
	//std::cout<<"sign of "<<p1<<" in the root of "<<a.pol()<<" = "<<s<<std::endl;
	//return s;
	return affiche_signs_constr (a);
}

// This function calculates the sign of the evaluation of a polynomial at a
// given algebraic number. If it is impossible to know the sign evaluating
// the interval, it calls sign_1_rs, which uses RS to do it.
Sign sign_1
(const Rational_polynomial_1 &p,const Algebraic_1 &x,unsigned int prec){
	mpfi_t y;
	mpfi_init2(y,mpfi_get_prec(x.mpfi()));
	eval_1(p,x,y);
	if(mpfr_zero_p(&(y->right))&&mpfr_zero_p(&(y->left))){
		mpfi_clear(y);
		return ZERO;
	}
	if((mpfr_sgn(&(y->left))<0)&&(mpfr_sgn(&(y->right))<0)){
		mpfi_clear(y);
		return NEGATIVE;
	}
	if((mpfr_sgn(&(y->left))>0)&&(mpfr_sgn(&(y->right))>0)){
		mpfi_clear(y);
		return POSITIVE;
	}
	// interval arithmetic couldn't find the sign, we call RS
	mpfi_clear(y);
	CGAL_assertion(x.is_consistent());
	return sign_1_rs(p,x);
}

// new function (implement evaluation using Horner's rule)
Sign exactsignat(const Rational_polynomial_1 &p,mpfr_srcptr xcoord){
	// we have numbers H0=h0*2^e0 and X=x*2^ex
	mpz_t h0,x,temp;
	mp_exp_t e0,ex;
	// we convert the evaluation point (mpfr) to fit our needs
	mpz_init(x);
	ex=mpfr_get_z_exp(x,xcoord);
	// data from the polynomial
	mpz_t *c=p.get_coefs();
	int d=p.get_degree();
	// we start the algorithm by setting H0=c[d]*2^0
	mpz_init_set(h0,c[d]);
	e0=0;
	// the iteration is H0=H0*X+c[d-i]
	mpz_init(temp);
	for(int i=1;i<d+1;++i){
		mpz_mul(h0,h0,x);
		e0+=ex;	// at this point, we did H0*=X
		// now, we have to add H0+c[d-i]
		if(e0>0){	// shift h0 and add
			mpz_mul_2exp(h0,h0,e0);
			mpz_add(h0,h0,c[d-i]);
			e0=0;
		}else{
			if(e0<0){	// shift c[d-i] and add
				mpz_mul_2exp(temp,c[d-i],-e0);
				mpz_add(h0,h0,temp);
				e0=0;
			}else	// we are lucky, e0=0
				mpz_add(h0,h0,c[d-i]);
		}
		// at this point, H0 is the evaluation of the polynomial
	}
	int s=mpz_sgn(h0);
	mpz_clear(h0);
	mpz_clear(x);
	mpz_clear(temp);
	switch(s){
		case 1:return POSITIVE;
		case -1:return NEGATIVE;
		default:return ZERO;
	}
}

// This function uses interval arithmetic to calculate the sign of the
// evaluation. First, it converts the given mpfr number to an mpfi interval
// with the given precision. Later, it evaluates the polynomial in the
// interval an looks for the sign. If it is not able to determine it, it
// calls the exact signat function (which is slower when the size of the
// mpfr number is big).
Sign quicksignat
(const Rational_polynomial_1 &p,mpfr_srcptr xcoord,mp_prec_t prec){
	mpfi_t x,result;
	int d=p.get_degree();
	mpz_t *c=p.get_coefs();
	mpfi_init2(x,prec);
	mpfi_init2(result,prec);
	mpfi_set_fr(x,xcoord);
	mpfi_set_z(result,c[d]);
	for(int i=1;i<d+1;++i){
		mpfi_mul(result,result,x);
		mpfi_add_z(result,result,c[d-i]);
	}
	mpfi_clear(x);
	if(mpfi_is_strictly_pos(result)){
		mpfi_clear(result);
		return POSITIVE;
	}
	if(mpfi_is_strictly_neg(result)){
		mpfi_clear(result);
		return NEGATIVE;
	}
	if(mpfi_is_zero(result)){
		mpfi_clear(result);
		return ZERO;
	}
	// interval arithmetic couldn't find the sign, we have to do it better
	mpfi_clear(result);
	return exactsignat(p,xcoord);
}

// This function is supposed to return the signat calculated in the best way.
// Obviously, it must be tested and fixed.
Sign signat(const Rational_polynomial_1 &p,mpfr_srcptr xcoord){
	mp_prec_t xprec=mpfr_get_prec(xcoord);
	if(xprec<160)
		return exactsignat(p,xcoord);
	if(xprec<640)
		return quicksignat(p,xcoord,xprec/4);
	return quicksignat(p,xcoord,xprec/8);
}

// this function returns the number of refinements made
// TODO: this function should be applied to a square-free polynomial
// TODO: study when is it better to call quicksignat and when signat
int refine(Algebraic_1 &a,unsigned n){
	Sign sl,sc;
	mp_prec_t pl,pr,pc;
	mpfr_t center;
	unsigned i;
	pl=mpfr_get_prec(&(a.mpfi()->left));
	sl=signat(a.pol(),&(a.mpfi()->left));
	if(sl==ZERO)
		return 0;
	pr=mpfr_get_prec(&(a.mpfi()->right));
	mpfr_init2(center,MPFR_PREC_MIN);
	for(i=0;i<n;++i){
		pc=pl<pr?pr+1:pl+1;
		mpfr_set_prec(center,pc);
		mpfr_add(center,&(a.mpfi()->left),&(a.mpfi()->right),GMP_RNDN);
		mpfr_div_2ui(center,center,1,GMP_RNDN);
		sc=signat(a.pol(),center);
		if(sc==ZERO){	// we have a root
			mpfr_set_prec(&(a.mpfi()->left),pc);
			mpfr_set(&(a.mpfi()->left),center,GMP_RNDN);
			// with mpfr_swap, we don't have to worry about prec.
			mpfr_swap(&(a.mpfi()->right),center);
			break;
		}else{
			if(sc==sl)
				mpfr_swap(&(a.mpfi()->left),center);
			else
				mpfr_swap(&(a.mpfi()->right),center);
		}
	}
	mpfr_clear(center);
	return i;
}

#if 0
int get_root (mpfi_ptr x, int n) {
	int ident_sols_eqs = rs_get_default_sols_eqs ();
	int ident_node = rs_export_list_vect_ibfr_firstnode (ident_sols_eqs);
	for (int i=0; i<rs_export_list_vect_ibfr_nb (ident_sols_eqs); ++i) {
		if (i == n) {
			mpfi_set (x, (mpfi_ptr)rs_export_ibfr_mpfi
					(rs_export_elt_vect_ibfr
					 (rs_export_list_vect_ibfr_monnode
					  (ident_node), 0)));
			return 1;	// true: we succeeded
		}
		ident_node = rs_export_list_vect_ibfr_nextnode (ident_node);
	}
	return 0;	// false: we couldn't copy the root
}
#endif

// TODO: rewrite using gcd
Comparison_result compare_1(Algebraic_1 &r1,Algebraic_1 &r2){
	if(r1.overlaps(r2)){
		if((r1.pol()==r2.pol())&&(r1.nr()==r2.nr()))
			return EQUAL;
		if(sign_1(r2.pol(),r1)==ZERO)
			return EQUAL;
		else{
			do{
				refine(r1,5);
				refine(r2,5);
			}while(r1.overlaps(r2));
		}
	}
	return(r1<r2?SMALLER:LARGER);
}

CGAL_END_NAMESPACE
