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
	rs_init_rs ();
	return 0;
}

int affiche_sols_eqs(mpfi_ptr *&x){
	int ident_sols_eqs,nb_elts,ident_node,ident_vect,i,ident_elt;
	ident_sols_eqs=rs_get_default_sols_eqs ();
	nb_elts=rs_export_list_vect_ibfr_nb(ident_sols_eqs);
	ident_node=rs_export_list_vect_ibfr_firstnode(ident_sols_eqs);
	x=(mpfi_ptr*)malloc(nb_elts*sizeof(mpfi_ptr));
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

Sign affiche_sols_constr(const Algebraic_1 &a){
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
	/*std::cout << "\nreturned value: ";
	mpfi_out_str(stdout,10,0,tmp);
	std::cout << std::endl;*/
	// mpfi_is_zero(tmp) doesn't work. The reason is that MPFR_SIGN in the
	// mpfi code returns 1 when applied to the left and right zeros. This
	// is not surprising, because zero is signed in IEEE 754-1985, and MPFR
	// adopts it. Nevertheless, mpfr_sgn returns 0, but mpfi doesn't use
	// it to implement mpfi_is_zero.
	// Here is the difference (from MPFR source code):
	// #define mpfr_sgn(_x)      (mpfr_zero_p(_x) ? 0 : MPFR_SIGN(_x))
	// Why? zeros are signed in IEEE-754.
	if(mpfr_zero_p(&(tmp->right))&&mpfr_zero_p(&(tmp->left)))
		return ZERO;
	// the same holds for mpfi_is_pos and mpfi_is_neg
	if((mpfr_sgn(&(tmp->left))>=0)&&(mpfr_sgn(&(tmp->right)))>0)
		return POSITIVE;
	if((mpfr_sgn(&(tmp->left))<0)&&(mpfr_sgn(&(tmp->right))<=0))
		return NEGATIVE;
	// if we arrive here, it is because the signs of the endpoints are -
	// and +, and (I think) RS guarantees that this never happen
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

int solve_1(mpfi_ptr *&x,Rational_polynomial_1 &p1,unsigned int prec){
	rs_reset_all();
	create_rs_upoly(p1.get_coefs(),p1.get_degree(),rs_get_default_up());
	set_rs_precisol(prec);
	set_rs_verbose(CGAL_RS_VERB);
	rs_run_algo("UISOLE");
	return affiche_sols_eqs(x);
}

Sign sign_1(const Rational_polynomial_1 &p1,const Algebraic_1 &a,
		unsigned int prec){
	CGAL_assertion(a.is_consistent());
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
	return affiche_sols_constr (a);
}

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

int refine_1_rs(Algebraic_1 &a){
	rs_reset_all ();
	create_rs_upoly (a.pol().get_coefs (), a.pol().get_degree (),
			rs_get_default_up ());
	int newprec = a.rsprec() * CGAL_RS_PREC_FACTOR;
	a.set_rsprec (newprec);
	set_rs_precisol (newprec);
	set_rs_verbose (CGAL_RS_VERB);
	rs_run_algo ("UISOLE");
	return get_root (a.mpfi(), a.nr());
}

// TODO: increase the precision by bigger steps, to avoid creating and erasing
// many times the mpfrs
int refine_1(Algebraic_1 &a){
	mpfr_t left,center,right,eval_l,eval_c,eval_r;
	mpfr_inits(left,right,eval_l,NULL);
	a.get_endpoints(left,right);
	mp_prec_t prec_l=mpfr_get_prec(left);
	a.pol().eval_mpfr(eval_l,left,prec_l);
	int sign_l=mpfr_sgn(eval_l);
	if(sign_l==0){
		mpfr_clears(left,right,eval_l,NULL);
		return refine_1_rs(a);
	}
	mpfr_init(eval_r);
	mp_prec_t prec_r=mpfr_get_prec(right);
	a.pol().eval_mpfr(eval_r,right,prec_r);
	int sign_r=mpfr_sgn(eval_r);
	if(sign_r==0){
		mpfr_clears(left,right,eval_l,eval_r,NULL);
		return refine_1_rs(a);
	}
	CGAL_assertion(sign_l!=sign_r);
	mp_prec_t prec_c=prec_l<prec_r?prec_r:prec_l;
	mpfr_inits2(prec_c,center,eval_c,NULL);
	mpfi_get_fr(center,a.mpfi());
	a.pol().eval_mpfr(eval_c,center,prec_c);
	int sign_c=mpfr_sgn(eval_c);
	if(sign_c==0){
		mpfr_clears(left,center,right,eval_l,eval_c,eval_r,NULL);
		return refine_1_rs(a);
	}
	if(sign_l==sign_c)
		mpfi_interv_fr(a.mpfi(),center,right);
	else
		mpfi_interv_fr(a.mpfi(),left,center);
	a.set_rsprec(1+a.rsprec());
	mpfr_clears(left,center,right,eval_l,eval_c,eval_r,NULL);
	return 1;
}

Comparison_result compare_1 (Algebraic_1 &r1, Algebraic_1 &r2) {
	try {return ((r1==r2)?EQUAL:((r1<r2)?SMALLER:LARGER));}
	catch (CGAL::comparison_overlap_exn &o) {
		int p1, p2;
		if (((p1=r1.rsprec())>CGAL_RS_MAX_PREC) &&
				((p2=r2.rsprec())>CGAL_RS_MAX_PREC))
			return EQUAL;
		refine_1 ((p1<p2)?r1:r2);
		return compare_1 (r1, r2);
	}
}

CGAL_END_NAMESPACE
