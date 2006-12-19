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
	/*std::cout << "\nreturned value: ";
	mpfi_out_str(stdout,10,0,tmp);
	std::cout << std::endl;*/
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
	Sign s=affiche_signs_constr(a);
	//std::cout<<"sign of "<<p1<<" in the root of "<<a.pol()<<" = "<<s<<std::endl;
	return s;
	//return affiche_signs_constr (a);
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

// TODO: rewrite this awful function
// TODO: test cases where RS is called
Comparison_result refine_and_compare_1(Algebraic_1 &r1,Algebraic_1 &r2){
	mpfr_t left1,right1,center1,evall1,evalc1,evalr1,
	       left2,right2,center2,evall2,evalc2,evalr2;
	mp_prec_t prec1,prec2,prec;
	if((prec1=r1.rsprec())<(prec2=r2.rsprec()))
		prec=(prec2*=2);
	else
		prec=(prec1*=2);
	mp_prec_t local_prec=prec;
	mpfr_inits2(local_prec,left1,right1,center1,evall1,evalc1,evalr1,
			left2,right2,center2,evall2,evalc2,evalr2,NULL);
	bool flag1,flag2;
	int res1,res2;
	int prec_changes=1;
	do{
		// refine r1
		flag1=true;
		r1.get_endpoints(left1,right1);
		mpfi_get_fr(center1,r1.mpfi());
		r1.pol().eval_mpfr(evall1,left1,prec);
		r1.pol().eval_mpfr(evalc1,center1,prec);
		r1.pol().eval_mpfr(evalr1,right1,prec);
		int sl1=mpfr_sgn(evall1);
		int sr1=mpfr_sgn(evalr1);
		// if the following assertion fails, it means that the
		// precision of the mpfr's was not correctly calculated
		//std::cout<<"sl1="<<sl1<<", sr1="<<sr1<<std::endl;
		CGAL_assertion(sl1&&sr1&&(sl1!=sr1));
		int sc1=mpfr_sgn(evalc1);
		while(r1.rsprec()<(int)prec){
			if(!sc1){
				refine_1_rs(r1);
				r1.get_endpoints(left1,right1);
				flag1=false;
				break;
			}else{
				if(sl1==sc1){
					mpfr_set(left1,center1,GMP_RNDN);
					mpfr_set(evall1,evalc1,GMP_RNDN);
					sl1=sc1;
				}else{
					mpfr_set(right1,center1,GMP_RNDN);
					mpfr_set(evalr1,evalc1,GMP_RNDN);
					sr1=sc1;
				}
				mpfr_add(center1,left1,right1,GMP_RNDN);
				mpfr_div_ui(center1,center1,2,GMP_RNDN);
				r1.pol().eval_mpfr(evalc1,center1,prec);
				sc1=mpfr_sgn(evalc1);
				r1.set_rsprec(r1.rsprec()+1);
			}
		}
		// refine r2
		flag2=true;
		r2.get_endpoints(left2,right2);
		mpfi_get_fr(center2,r2.mpfi());
		r2.pol().eval_mpfr(evall2,left2,prec);
		r2.pol().eval_mpfr(evalc2,center2,prec);
		r2.pol().eval_mpfr(evalr2,right2,prec);
		int sl2=mpfr_sgn(evall2);
		int sr2=mpfr_sgn(evalr2);
		CGAL_assertion(sl2&&sr2&&(sl2!=sr2));
		int sc2=mpfr_sgn(evalc2);
		while(r2.rsprec()<(int)prec){
			if(!sc2){
				refine_1_rs(r2);
				r2.get_endpoints(left2,right2);
				flag2=false;
				break;
			}else{
				if(sl2==sc2){
					mpfr_set(left2,center2,GMP_RNDN);
					mpfr_set(evall2,evalc2,GMP_RNDN);
					sl2=sc2;
				}else{
					mpfr_set(right2,center2,GMP_RNDN);
					mpfr_set(evalr2,evalc2,GMP_RNDN);
					sr2=sc2;
				}
				mpfr_add(center2,left2,right2,GMP_RNDN);
				mpfr_div_ui(center2,center2,2,GMP_RNDN);
				r2.pol().eval_mpfr(evalc2,center2,prec);
				sc2=mpfr_sgn(evalc2);
				r2.set_rsprec(r2.rsprec()+1);
			}
		}
		res1=mpfr_less_p(right1,left2);
		res2=mpfr_less_p(right2,left1);
		if((prec*=2)>local_prec){
			++prec_changes;
			local_prec*=(prec_changes*prec_changes);
			// change al the precisions
			mpfr_set_prec(left1,local_prec);
			mpfr_set_prec(center1,local_prec);
			mpfr_set_prec(right1,local_prec);
			mpfr_set_prec(evall1,local_prec);
			mpfr_set_prec(evalc1,local_prec);
			mpfr_set_prec(evalr1,local_prec);
			mpfr_set_prec(left2,local_prec);
			mpfr_set_prec(center2,local_prec);
			mpfr_set_prec(right2,local_prec);
			mpfr_set_prec(evall2,local_prec);
			mpfr_set_prec(evalc2,local_prec);
			mpfr_set_prec(evalr2,local_prec);
		}
	}while(!res1&&!res2);
	if(flag1)
		mpfi_interv_fr(r1.mpfi(),left1,right1);
	if(flag2)
		mpfi_interv_fr(r2.mpfi(),left2,right2);
	mpfr_clears(left1,right1,center1,evall1,evalc1,evalr1,
			left2,right2,center2,evall2,evalc2,evalr2,NULL);
	if(res1)
		return SMALLER;
	if(res2)
		return LARGER;
	CGAL_assertion_msg(false,"this point should never be reached");
	return EQUAL;
}

Comparison_result compare_1(Algebraic_1 &r1,Algebraic_1 &r2){
	// there is a strange RS result where (r1<r2)&&(r2<r1)&&(r1!=r2)
	try{
		if((r1==r2)/*||((r1<r2)&&(r2<r1))*/)
			return EQUAL;
	}
	catch(CGAL::comparison_overlap_exn &o){
		if((r1.pol()==r2.pol())||(sign_1(r2.pol(),r1)==ZERO))
			return EQUAL;
	}
	if(r1.overlaps(r2))
		if((r1.pol()==r2.pol())||(sign_1(r2.pol(),r1)==ZERO))
			return EQUAL;
		else
			return refine_and_compare_1(r1,r2);
	// at this point, we know the intervals aren't equal
	return(r1<r2?SMALLER:LARGER);
}

CGAL_END_NAMESPACE
