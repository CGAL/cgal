// ======================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.4-I-8 $
// release_date  : $CGAL_Date: 2001/09/11 $
//
// file          : config/testfiles/CGAL_CFG_NET2003_MATCHING_BUG.C
// package       : Configuration (2.12)
// maintainer    : Geert-Jan Giezeman <geert@cs.uu.nl>
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : Radu Ursu
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NET2003_MATCHING_BUG.C
// ---------------------------------------------------------------------
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set, if the compiler does not match member 
//| definition to an existing declaration (eg. cl1310 Beta)

#include <iostream>
#include <vector>

template <class Gt, class Tds>
class Triangulation_3{
  typedef typename Tds::Facet                  Facet;
  typedef typename Tds::Vertex_handle          Vertex_handle;
  typedef typename Tds::Cell_handle            Cell_handle;

};

template <class Gt, class Tds>
class D_Triangulation : public Triangulation_3<Gt, Tds>
{
  typedef Triangulation_3<Gt, Tds> Tr_base;

  typedef typename Tr_base::Cell_handle   Cell_handle;
  typedef typename Tr_base::Vertex_handle Vertex_handle;
  typedef typename Tr_base::Facet Facet;
  
  void make_hole(Vertex_handle, std::vector<Facet>&, 
    std::vector<Cell_handle>&);
  
};

template <class Gt, class Tds>
void
D_Triangulation<Gt, Tds>::make_hole(Vertex_handle v, std::vector<Facet>& f, 
std::vector<Cell_handle>& g){
  std::cout << "test";
}

int main(){
  return 1;
}