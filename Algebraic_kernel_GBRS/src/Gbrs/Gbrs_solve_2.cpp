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
// $URL: $
// $Id: $
// 
//
// Author(s)     : Luis Peñaranda <penarand@loria.fr>

#include <gmp.h>
#include <mpfr.h>
#include <CGAL/assertions.h>
#include <CGAL/Gbrs_polynomial_2.h>
#include <CGAL/Gbrs_solve_2.h>
#include <rs_exports.h>

#ifdef USE_FGB
#include "fgb_2_rs.h"
#endif

CGAL_BEGIN_NAMESPACE

void create_rs_bipoly_ring(char *v1,char * v2)
{
  int i=0;
  /* initialise l'anneau de polynômes univariés à coeffs entiers
     de RS : pas obligatoire mais mieux si on veut faire afficher
     des polynômes par RS */
  rs_init_ppring_nb(2,"DRL");
  rs_import_ppring_var(0,v1);
  rs_import_ppring_var(1,v2);
}

int create_rs_bipoly(MP_INT ** poly,int deg1,int deg2,int ident_pol)
{
  int i,j;
  int ident_ring,ident_mon,ident_coeff,ident_pp;
  /* demande à RS l'adresse de son polynome univarie global */
  for (i=0;i<deg1;i++) {
    for (j=0;j<deg2;j++) {
      /* demande a RS de creer un nouveau monôme */
      if (mpz_cmp_ui((&(poly[i][j])),0)) {
	ident_mon=rs_export_new_mon_pp_bz();
	/* demande a RS de creer un entier long */
	ident_coeff=rs_export_new_gmp();
	/* affecte l'entier gmp de RS avec l'entier gmp local */
	rs_import_bz_gmp(ident_coeff,TO_RSPTR_IN(&(poly[i][j])));
	/* affecte le monôme de RS avec l'entier local et le degré qui va bien */
	ident_pp=rs_export_new_pp();
	rs_set_pp_deg(ident_pp,0,i);
	rs_set_pp_deg(ident_pp,1,j);
	rs_dset_mon_pp_bz(ident_mon,ident_coeff,ident_pp);
	/* ajoute le nouveau monome au polynôme univarié par défaut de RS */
	rs_dappend_list_mon_pp_bz(ident_pol,ident_mon);
      }
    }
  }
  return 0;
}

int create_rs_bisys(MP_INT **p1,int deg11, int deg12,MP_INT **p2,int deg21, int deg22)
{
    int  ident_sys=rs_export_new_list_mon_pp_bz();
    int ident_poly=0;
    ident_sys=rs_get_default_grob();
    ident_poly=rs_export_new_list_mon_pp_bz();
    create_rs_bipoly(p1,deg11,deg12,ident_poly);
    rs_dappend_list_smp_bz(ident_sys,ident_poly);
    ident_poly=rs_export_new_list_mon_pp_bz();
    create_rs_bipoly(p2,deg21,deg22,ident_poly);
    rs_dappend_list_smp_bz(ident_sys,ident_poly);
    return 0;
}

/* fonction d'affichage d'un vecteur de RS à entrées des intervalles
   multi-precision (mpfi) */

void affiche_vect_ibfr(const int ident_vect)
{
  int ident_elt;
  int i;
  mpfi_t tmp;
  mpfi_init(tmp);
  /* on recupere la dimension du vecteur */
  int nb=rs_export_dim_vect_ibfr(ident_vect);
  fprintf(stderr,"(%d)[",nb);
  for (i=0;i<nb;i++) {
    /* on recupere le pointeur sur le mpfi de RS */
    ident_elt=rs_export_elt_vect_ibfr(ident_vect,i);
    /* on initialise le mpfi local */
    mpfi_set(tmp,(mpfi_ptr)rs_export_ibfr_mpfi(ident_elt));
    /* on affiche le mpfi local */
    fprintf(stderr,",");
    mpfi_out_str(stderr,10,0,tmp);
  }
  fprintf(stderr,"]\n");  
}

/* fonction de demonstration qui recupere les intervalles 
   d'isolation obtenus par RS.*/

void affiche_sols_eqs()
     /*     cette fonction fonctionne donc pour un polynome en une variable
            autant que pour une système d'équations*/
{
  int ident_sols_eqs,nb_elts,ident_node,ident_vect;
  int i;
  /* les solutions des systèmes d'équations (ou des polynômes en une
     variables calculées par RS sont toujours dans une variable par
     defaut : */ 
  ident_sols_eqs=rs_get_default_sols_eqs();
  /* les solutions sont une liste de vecteurs de longueur le nombre de 
     variables, les entrées de ces vecteurs sont des intervalles
     de type mpfi*/
  /* on attrappe le nombre d'éléments de la liste de RS */
  nb_elts=rs_export_list_vect_ibfr_nb(ident_sols_eqs);
  /* on attrappe le premier element de la liste de RS */
  ident_node=rs_export_list_vect_ibfr_firstnode(ident_sols_eqs);
  for (i=1;i<nb_elts+1;i++) {
    /* on dépacte le i-eme élément de la liste */
    ident_vect=rs_export_list_vect_ibfr_monnode(ident_node);
    affiche_vect_ibfr(ident_vect);
    /* on passe a l'élément suivant */
    ident_node=rs_export_list_vect_ibfr_nextnode(ident_node);
  }
}

void affiche_sols_constr()
{
  int ident_sols_eqs,nb_elts,ident_node,ident_vect;
  int i;
  ident_sols_eqs=rs_get_default_sols_ineqs();
  nb_elts=rs_export_list_vect_ibfr_nb(ident_sols_eqs);
  ident_node=rs_export_list_vect_ibfr_firstnode(ident_sols_eqs);
  for (i=1;i<nb_elts+1;i++) {
    ident_vect=rs_export_list_vect_ibfr_monnode(ident_node);
    fprintf(stderr,"\nContraintes en la solution %d : ",i);
    affiche_vect_ibfr(ident_vect);
    ident_node=rs_export_list_vect_ibfr_nextnode(ident_node);
  }
}

/* Fabrice's original function */
#ifdef kk
int solve_2(Rational_polynomial_2 &f,Rational_polynomial_2 &g)
{
  /* les coefficients d'un polynôme par ordre croissants de degrés 
     pour les besoins de la demo */
  /* le degré du polynôme */
     int deg11=2;
     int deg12=2;
          int poly1[]={1,2,3,4};
     //int poly1[]={0,1,1,0};
  /* le degré du polynôme */
     int deg21=2;
     int deg22=2;
     int poly2[]={1,-1,1,-1};
     //int poly2[]={0,-1,1,0};
     MP_INT ** p1=(MP_INT **)malloc(deg11*sizeof(MP_INT*));
     MP_INT ** p2=(MP_INT **)malloc(deg21*sizeof(MP_INT*));
  /* creation d'un polynôme local à coeffs entiers 
     pour les besoins de la demo */
     int i=0,j=0;
     for (i=0;i<deg11;i++) {
       p1[i]=(MP_INT *)malloc(deg12*sizeof(MP_INT));
       for (j=0;j<deg12;j++) {
	 mpz_init_set_si(&(p1[i][j]),poly1[i*deg11+j]);
	 mpz_out_str(stderr,10,&(p1[i][j])); fprintf(stderr," ");
       }
     }
     fprintf(stderr,"\n 2nd : \n");
     for (i=0;i<deg21;i++) {
       p2[i]=(MP_INT *)malloc(deg22*sizeof(MP_INT));
       for (j=0;j<deg22;j++) {
	 mpz_init_set_si(&(p2[i][j]),poly2[i*deg21+j]);
	 mpz_out_str(stderr,10,&(p2[i][j])); fprintf(stderr," ");
       }
     }

     /*  ca ca initialise les buffers et le gestionnaire 
     de mémoire à faire une seule fois par session */
     rs_init_rs();

     /* ca ca initialise les variables internes */ 
     rs_reset_all();

#ifdef USE_FGB
    set_rs_groebner(fgb_compute_DRL_INT);
#endif

     /* déclare un anneau de polynômes en une varaible */
     create_rs_bipoly_ring("x","y");

     
     /* affecte le polynôme univarié par défaut de RS */ 
     create_rs_bisys(p1,deg11,deg12,p2,deg21,deg22);

     /* affecte la précision d'isolation dans RS */
     set_rs_precisol(23);

     /* affecte le niveau de verbosité de RS */
     set_rs_verbose(2);

     /* demande à RS d'isoler les racines de son polynôme univarié par
	défaut*/
     rs_run_algo("SISOLE");

     /* récupère le résultat dans Rs et l'affiche */

     affiche_sols_eqs();
     affiche_sols_constr();    

     /* pour faire un autre calcul, il suffit de refaire la même chose
	a partir de l'instruction rs_reset_all() */
     return(0);
}
#endif

int solve_2(const Rational_polynomial_2 &f,const Rational_polynomial_2 &g)
{
     mpz_t** f_ptr=f.get_coefs();
     int deg11=f.get_degree_x()+1;
     int deg12=f.get_degree_y()+1;
     int deg21=g.get_degree_x()+1;
     int deg22=g.get_degree_y()+1;
     mpz_t** g_ptr=g.get_coefs();
     MP_INT ** p1=(MP_INT **)malloc(deg11*sizeof(MP_INT*));
     MP_INT ** p2=(MP_INT **)malloc(deg21*sizeof(MP_INT*));
     int i=0,j=0;
     for (i=0;i<deg11;i++) {
       p1[i]=(MP_INT *)malloc(deg12*sizeof(MP_INT));
       for (j=0;j<deg12;j++)
	 mpz_init_set(&(p1[i][j]),f_ptr[i][j]);
     }
     for (i=0;i<deg21;i++) {
       p2[i]=(MP_INT *)malloc(deg22*sizeof(MP_INT));
       for (j=0;j<deg22;j++)
	 mpz_init_set(&(p2[i][j]),g_ptr[i][j]);
     }

     /*  ca ca initialise les buffers et le gestionnaire 
     de mémoire à faire une seule fois par session */
     rs_init_rs();

     /* ca ca initialise les variables internes */ 
     rs_reset_all();

#ifdef USE_FGB
    set_rs_groebner(fgb_compute_DRL_INT);
#endif

     /* déclare un anneau de polynômes en une varaible */
     create_rs_bipoly_ring("x","y");

     
     /* affecte le polynôme univarié par défaut de RS */ 
     create_rs_bisys(p1,deg11,deg12,p2,deg21,deg22);

     /* affecte la précision d'isolation dans RS */
     set_rs_precisol(23);

     /* affecte le niveau de verbosité de RS */
     set_rs_verbose(2);

     /* demande à RS d'isoler les racines de son polynôme univarié par
	défaut*/
     rs_run_algo("SISOLE");

     /* récupère le résultat dans Rs et l'affiche */
	/* TODO: take results back */
     affiche_sols_eqs();
     affiche_sols_constr();    

     /* pour faire un autre calcul, il suffit de refaire la même chose
	a partir de l'instruction rs_reset_all() */
     return(0);
}

CGAL_END_NAMESPACE
