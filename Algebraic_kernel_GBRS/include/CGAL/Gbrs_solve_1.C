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
#include <mpfi.h>
#include <CGAL/assertions.h>
#include <CGAL/Gbrs_polynomial_1.h>
#include <rs_exports.h>

CGAL_BEGIN_NAMESPACE

// the default precision of gb/rs
#ifndef CGAL_RS_DEF_PREC
#define CGAL_RS_DEF_PREC 23
#endif

#ifndef CGAL_RS_MIN_PREC
#define CGAL_RS_MIN_PREC 5
#endif

#ifndef CGAL_RS_MAX_PREC
#define CGAL_RS_MAX_PREC 50
#endif

int init_rs () {
	rs_init_rs ();
	return 0;
}

int affiche_vect_ibfr (mpfi_t *&x, int index, const int ident_vect) {
	int ident_elt;
	int nb = rs_export_dim_vect_ibfr (ident_vect);
	CGAL_assertion_msg (nb == 1, "the dimension of vector must be 1");
	ident_elt = rs_export_elt_vect_ibfr (ident_vect, 0);
	mpfi_set (x[index-1], (mpfi_ptr)rs_export_ibfr_mpfi (ident_elt));
	return nb;
}

int affiche_sols_eqs (mpfi_t *&x) {
	int ident_sols_eqs, nb_elts, ident_node, ident_vect, i;
	ident_sols_eqs = rs_get_default_sols_eqs ();
	// the number of solutions
	nb_elts = rs_export_list_vect_ibfr_nb (ident_sols_eqs);
	ident_node = rs_export_list_vect_ibfr_firstnode (ident_sols_eqs);
	// allocate space
	x = (mpfi_t *) malloc (nb_elts*sizeof(mpfi_t));
	for (i=0; i<nb_elts; i++) {
		mpfi_init (x[i]);
		ident_vect = rs_export_list_vect_ibfr_monnode (ident_node);
		affiche_vect_ibfr (x, i+1, ident_vect);
		ident_node = rs_export_list_vect_ibfr_nextnode (ident_node);
	}
	return nb_elts;
}

Sign affiche_sols_constr (const Algebraic_1 &a) {
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
	/*std::cout << "returned value: ";
	mpfi_out_str(stdout,10,0,tmp);
	std::cout << std::endl;*/
	if (mpfi_is_strictly_pos (tmp))
		return POSITIVE;
	if (mpfi_is_strictly_neg (tmp))
		return NEGATIVE;
	// TODO: refine the result when needed
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

int solve_1 (mpfi_t *&x, const Rational_polynomial_1 &p1, unsigned int prec) {
	// the solver must be initialized
	rs_reset_all ();
	create_rs_upoly
		(p1.get_coefs (), p1.get_degree (), rs_get_default_up ());
	set_rs_precisol (prec);	// precision
	set_rs_verbose (0);	// we don't want any output
	rs_run_algo ("UISOLE");	// isolate roots
	// XXX remember to free the results array x and
	// clear all the mpfi_t elements
	return affiche_sols_eqs (x);	// return the number of solutions
}

inline int solve_1 (mpfi_t *&x, const Rational_polynomial_1 &p1) {
	return solve_1 (x, p1, CGAL_RS_DEF_PREC);
}

Sign sign_1 (const Rational_polynomial_1 &p1, const Algebraic_1 &a) {
	mpz_t **constr;
	int *degs;
	CGAL_assertion_msg (a.is_consistent (),
		"this number was not calculated as a root of a polynomial");
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
	set_rs_precisol (CGAL_RS_MIN_PREC);
	set_rs_verbose (0);
	rs_run_algo ("UISOLES");
	return affiche_sols_constr (a);
}

CGAL_END_NAMESPACE
