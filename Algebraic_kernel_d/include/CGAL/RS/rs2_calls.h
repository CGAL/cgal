// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.

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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_RS2_CALLS_H
#define CGAL_RS_RS2_CALLS_H

#include <CGAL/Gmpz.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/Gmpfi.h>
#include <CGAL/Polynomial.h>
#include <CGAL/tss.h>

#include <rs_exports.h>

#ifdef CGAL_RS_OLD_INCLUDES
#define CGALRS_PTR(a)   long int a
#else
#define CGALRS_PTR(a)   void *a
#endif

// RS3 does not work with MPFR 3.1.3 to 3.1.6. In case RS3 is enabled and
// the version of MPFR is one of those buggy versions, abort the compilation
// and instruct the user to update MPFR or don't use RS3.
#ifdef CGAL_USE_RS3
#include <boost/static_assert.hpp>
BOOST_STATIC_ASSERT_MSG(
        MPFR_VERSION_MAJOR!=3 ||
        MPFR_VERSION_MINOR!=1 ||
        MPFR_VERSION_PATCHLEVEL<3 || MPFR_VERSION_PATCHLEVEL>6,
        "RS3 does not work with MPFR versions 3.1.3 to 3.1.6. "
        "Please update MPFR or disable RS3.");
#endif // CGAL_USE_RS3

namespace CGAL{
namespace RS2{

struct RS2_calls{

        static void init_solver(){
                CGAL_STATIC_THREAD_LOCAL_VARIABLE(bool, first,true);
                if(first){
                        first=false;
                        rs_init_rs();
                        rs_reset_all();
                }else
                        rs_reset_all();
        }

        static void create_rs_upoly(CGAL::Polynomial<CGAL::Gmpz> poly,
                                    CGALRS_PTR(ident_pol)){
                CGALRS_PTR(ident_mon);
                CGALRS_PTR(ident_coeff);
                rs_import_uppring((char*)"T");
                for(int i=0;i<=poly.degree();++i)
                        if(mpz_sgn(poly[i].mpz())){   // don't add if == 0
                                ident_mon=rs_export_new_mon_upp_bz();
                                ident_coeff=rs_export_new_gmp();
                                rs_import_bz_gmp(ident_coeff,
                                                 TO_RSPTR_IN(&(poly[i].mpz())));
                                rs_dset_mon_upp_bz(ident_mon,ident_coeff,i);
                                rs_dappend_list_mon_upp_bz(ident_pol,
                                                           ident_mon);
                        }
        }

        static int affiche_sols_eqs(mpfi_ptr *x){
                CGALRS_PTR(ident_sols_eqs);
                CGALRS_PTR(ident_node);
                CGALRS_PTR(ident_vect);
                CGALRS_PTR(ident_elt);
                int nb_elts;
                ident_sols_eqs=rs_get_default_sols_eqs();
                nb_elts=rs_export_list_vect_ibfr_nb(ident_sols_eqs);
                ident_node=rs_export_list_vect_ibfr_firstnode(ident_sols_eqs);
                mpfi_t *roots=(mpfi_t*)malloc(nb_elts*sizeof(mpfi_t));
                for(int i=0;i<nb_elts;++i){
                        ident_vect=rs_export_list_vect_ibfr_monnode
                                (ident_node);
                        CGAL_assertion_msg(rs_export_dim_vect_ibfr
                                                (ident_vect)==1,
                                           "vector dimension must be 1");
                        ident_elt=rs_export_elt_vect_ibfr(ident_vect,0);
                        mpfi_ptr root_pointer=
                                (mpfi_ptr)rs_export_ibfr_mpfi(ident_elt);
                        mpfi_init2(roots[i],mpfi_get_prec(root_pointer));
                        mpfi_set(roots[i],root_pointer);
                        x[i]=roots[i];
                        // This doesn't work because RS relocates the
                        // mpfrs that form the mpfi. Nevertheless, the
                        // mpfi address is not changed.
                        //x[i]=(mpfi_ptr)rs_export_ibfr_mpfi(ident_elt);
                        ident_node=rs_export_list_vect_ibfr_nextnode
                                (ident_node);
                }
                return nb_elts;
        }

        template<class OutputIterator>
        static OutputIterator insert_roots(OutputIterator x){
                CGALRS_PTR(ident_sols_eqs);
                CGALRS_PTR(ident_node);
                CGALRS_PTR(ident_vect);
                CGALRS_PTR(ident_elt);
                int nb_elts;
                ident_sols_eqs=rs_get_default_sols_eqs();
                nb_elts=rs_export_list_vect_ibfr_nb(ident_sols_eqs);
                ident_node=rs_export_list_vect_ibfr_firstnode(ident_sols_eqs);
                for(int i=0;i<nb_elts;++i){
                        ident_vect=rs_export_list_vect_ibfr_monnode
                                (ident_node);
                        CGAL_assertion_msg(rs_export_dim_vect_ibfr
                                                (ident_vect)==1,
                                           "vector dimension must be 1");
                        ident_elt=rs_export_elt_vect_ibfr(ident_vect,0);
                        mpfi_ptr root_pointer=
                                (mpfi_ptr)rs_export_ibfr_mpfi(ident_elt);
                        mp_prec_t root_prec=mpfi_get_prec(root_pointer);
                        // Construct Gmpfr's with pointers to endpoints.
                        Gmpfr left(&(root_pointer->left),root_prec);
                        Gmpfr right(&(root_pointer->right),root_prec);
                        CGAL_assertion(left<=right);
                        // Copy them, to have the data out of RS memory.
                        *x++=Gmpfi(left,right,root_prec+1);
                        ident_node=rs_export_list_vect_ibfr_nextnode
                                (ident_node);
                }
                return x;
        }

}; // struct RS2_calls

} // namespace RS2
} // namespace CGAL

#endif // CGAL_RS_RS2_CALLS_H
