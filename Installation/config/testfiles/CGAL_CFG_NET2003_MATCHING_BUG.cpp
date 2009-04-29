// Copyright (c) 1997-2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// Author(s)     : Radu Ursu

//| This flag is set, if the compiler does not match a member 
//| definition to an existing declaration (eg., cl1310 Beta).

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
  return 0;
}
