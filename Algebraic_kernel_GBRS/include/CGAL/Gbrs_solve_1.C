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

// this file is based on the RS original example, rs_demo_mp.c

#include <gmp.h>
#include <mpfi.h>
#include <CGAL/assertions.h>
#include <CGAL/Gbrs_polynomial_1.h>
#include <rs_exports.h>

CGAL_BEGIN_NAMESPACE

// the default precision of gb/rs
#ifndef DEF_PREC
#define DEF_PREC 23
#endif

int affiche_vect_ibfr (mpfi_t *&x, int index, const int ident_vect) {
	int ident_elt;
	int nb = rs_export_dim_vect_ibfr (ident_vect);
	CGAL_assertion_msg (nb == 1, "the dimension of vector must be 1");
	ident_elt = rs_export_elt_vect_ibfr (ident_vect, 0);
	mpfi_set (x[index-1], rs_export_ibfr_mpfi (ident_elt));
	return nb;
}

int affiche_sols_eqs (mpfi_t *&x) {
	int ident_sols_eqs, nb_elts, ident_node, ident_vect;
	int i;
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

void create_rs_upoly(mpq_t *poly, const int deg) {
	int i;
	int ident_pol, ident_mon, ident_coeff;
	rs_import_uppring ("T");
	ident_pol = rs_get_default_up ();
	for (i=0; i<deg+1; i++)
		if (mpq_sgn (poly[deg-i]) != 0)	{	// don't add if == 0
			ident_mon = rs_export_new_mon_upp_bz ();
			ident_coeff = rs_export_new_gmp ();
			CGAL_assertion_msg
				(mpz_cmp_ui (mpq_denref (poly[deg-i]), 1) == 0,
				 "by now, RS works only with integer coeffs");
			rs_import_bz_gmp
				(ident_coeff,
				 TO_RSPTR_IN (mpq_numref (poly[deg-i])));
			/*rs_import_bz_gmp
				(ident_coeff,
				 TO_RSPTR_IN (&(poly[deg-i])));*/
			rs_dset_mon_upp_bz (ident_mon, ident_coeff, i);
			rs_dappend_list_mon_upp_bz (ident_pol, ident_mon);
		}
}

int solve_1 (mpfi_t *&x, const Rational_polynomial_1 &p1, unsigned int prec) {
	rs_init_rs ();
	rs_reset_all ();
	create_rs_upoly (p1.get_coefs (), p1.get_degree ());	// the RS poly
	set_rs_precisol (prec);	// preciseness
	set_rs_verbose (0);	// we don't want any output
	rs_run_algo ("UISOLE");	// isolate roots
	// XXX remember to free the results array x and
	// clear all the mpfi_t elements
	return affiche_sols_eqs (x);	// return the number of solutions
}

int solve_1 (mpfi_t *&x, const Rational_polynomial_1 &p1) {
	return solve_1 (x, p1, DEF_PREC);
}

CGAL_END_NAMESPACE
