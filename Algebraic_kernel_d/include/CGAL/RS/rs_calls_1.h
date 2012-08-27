// Copyright (c) 2006-2009 Inria Lorraine (France). All rights reserved.
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
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_RS_CALLS_1_H
#define CGAL_RS_RS_CALLS_1_H

#include <CGAL/RS/basic.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>
#include <rs_exports.h>

namespace CGAL{

#define CGALRS_CSTR(S)  ((char*)(S))

#ifdef CGAL_RS_OLD_INCLUDES
#define CGALRS_PTR(a)   long int a
#else
#define CGALRS_PTR(a)   void *a
#endif

// initialize RS solver
inline void init_solver(){
        static bool first=true;
        if(first){
                first=false;
                //rs_version(stderr);
                rs_init_rs();
        }else{
                rs_reset_all();
        }
}

// reset RS memory
inline void reset_solver(){
        rs_reset_all();
}

inline int affiche_sols_eqs(mpfi_ptr *x){
        CGALRS_PTR(ident_sols_eqs);
        CGALRS_PTR(ident_node);
        CGALRS_PTR(ident_vect);
        CGALRS_PTR(ident_elt);
        int nb_elts;
        ident_sols_eqs=rs_get_default_sols_eqs();
        nb_elts=rs_export_list_vect_ibfr_nb(ident_sols_eqs);
        ident_node=rs_export_list_vect_ibfr_firstnode(ident_sols_eqs);
        //x=(mpfi_ptr*)malloc(nb_elts*sizeof(mpfi_ptr));
        for(int i=0;i<nb_elts;++i){
                ident_vect=rs_export_list_vect_ibfr_monnode(ident_node);
                CGAL_assertion_msg(rs_export_dim_vect_ibfr(ident_vect)==1,
                                "the dimension of vector must be 1");
                ident_elt=rs_export_elt_vect_ibfr(ident_vect,0);
                x[i]=(mpfi_ptr)rs_export_ibfr_mpfi(ident_elt);
                ident_node=rs_export_list_vect_ibfr_nextnode(ident_node);
        }
        return nb_elts;
}

inline void affiche_sols_constr(int nr,mpfi_ptr p){
        CGALRS_PTR(ident_sols_eqs);
        CGALRS_PTR(ident_node);
        CGALRS_PTR(ident_vect);
        CGALRS_PTR(ident_elt);
        int nb_elts,nb;
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

// nroot is the number of root of the algebraic number
inline Sign affiche_signs_constr(int nroot){
        CGALRS_PTR(ident_sols_eqs);
        CGALRS_PTR(ident_node);
        CGALRS_PTR(ident_vect);
        CGALRS_PTR(ident_elt);
        int nb_elts,nb;
        mpfi_t tmp;
        mpfi_init(tmp);
        ident_sols_eqs=rs_get_default_sols_ineqs();
        nb_elts=rs_export_list_vect_ibfr_nb(ident_sols_eqs);
        ident_node=rs_export_list_vect_ibfr_firstnode(ident_sols_eqs);
        for(int i=1;i<nb_elts+1;++i){
                ident_vect=rs_export_list_vect_ibfr_monnode(ident_node);
                nb=rs_export_dim_vect_ibfr(ident_vect);
                CGAL_assertion_msg((nb==1),
                                "the vector must contain one element");
                ident_elt=rs_export_elt_vect_ibfr(ident_vect,0);
                if(i==nroot+1){
                        mpfi_set(tmp,(mpfi_ptr)rs_export_ibfr_mpfi(ident_elt));
                        break;
                }
                ident_node=rs_export_list_vect_ibfr_nextnode(ident_node);
        }
        //std::cerr << "\nreturned value: ";
        //mpfi_out_str(stderr,10,0,tmp);
        //std::cerr << std::endl;
        // mpfi_is_zero(tmp) doesn't work. The reason is that MPFR_SIGN in
        // the mpfi code returns 1 when applied to the left and right zeros.
        // This is not surprising because zero is signed in IEEE 754, and MPFR
        // adopts it. Nevertheless, mpfr_sgn returns 0, but mpfi doesn't use
        // it to implement mpfi_is_zero.
        // Here is the difference (from MPFR source code):
        //  define mpfr_sgn(_x)      (mpfr_zero_p(_x) ? 0 : MPFR_SIGN(_x))
        //
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

inline void create_rs_upoly(mpz_t *poly,const int deg,CGALRS_PTR(ident_pol)){
        CGALRS_PTR(ident_mon);
        CGALRS_PTR(ident_coeff);
        rs_import_uppring(CGALRS_CSTR("T"));
        for(int i=0;i<=deg;++i)
                if(mpz_sgn(poly[i])){   // don't add if == 0
                        ident_mon=rs_export_new_mon_upp_bz();
                        ident_coeff=rs_export_new_gmp();
                        rs_import_bz_gmp
                                (ident_coeff,TO_RSPTR_IN(&(poly[i])));
                        rs_dset_mon_upp_bz(ident_mon,ident_coeff,i);
                        rs_dappend_list_mon_upp_bz(ident_pol,ident_mon);
                }
}

inline void create_rs_uconstr(mpz_t **list_constr,
                               const int *list_degs,
                               CGALRS_PTR(ident_list)){
        CGALRS_PTR(ident_poly);
        ident_poly=rs_export_new_list_mon_upp_bz();
        create_rs_upoly(*list_constr,*list_degs,ident_poly);
        rs_dappend_list_sup_bz(ident_list,ident_poly);
}

} // namespace CGAL

#endif  // CGAL_RS_RS_CALLS_1_H
