// ======================================================================
//
// Copyright (c) 2005-2011 GeometryFactory (France).  All Rights Reserved.
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
//
// Author(s): Le-Jeng Shiue <Andy.Shiue@gmail.com>
//
// ======================================================================

#ifndef CGAL_POLYHEDRON_SUBDIVISION_STENCILS_H_01292002
#define CGAL_POLYHEDRON_SUBDIVISION_STENCILS_H_01292002

#include <CGAL/basic.h>
#include <CGAL/Origin.h>

#include <CGAL/circulator.h>

namespace CGAL {

// ======================================================================
/// The stencil of the Primal-Quadrilateral-Quadrisection 
template <class Poly>
class PQQ_stencil_3 {
public:
  typedef Poly                                         Polyhedron;

  typedef typename Polyhedron::Vertex_handle           Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
  typedef typename Polyhedron::Facet_handle            Facet_handle;

  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;

  typedef typename Kernel::FT                          FT;
  typedef typename Kernel::Point_3                     Point;
  typedef typename Kernel::Vector_3                    Vector;

public:
  void facet_node(Facet_handle, Point&) {};
  void edge_node(Halfedge_handle, Point&) {};
  void vertex_node(Vertex_handle, Point&) {};

  void border_node(Halfedge_handle, Point&, Point&) {};
};


// ======================================================================
/// Bi-linear geometry mask for PQQ, PTQ, and Sqrt(3) scheme 
template <class Poly>
class Linear_mask_3 : public PQQ_stencil_3<Poly> {
public:
  typedef Poly                                         Polyhedron;

  typedef typename Polyhedron::Vertex_handle           Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
  typedef typename Polyhedron::Facet_handle            Facet_handle;

  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator 
                                            Halfedge_around_vertex_circulator;

  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;

  typedef typename Kernel::FT                          FT;
  typedef typename Kernel::Point_3                     Point;
  typedef typename Kernel::Vector_3                    Vector;

public:
  //
  void facet_node(Facet_handle facet, Point& pt) {
    Halfedge_around_facet_circulator hcir = facet->facet_begin();
    int n = 0;
    Point p(0,0,0);
    do {
      p = p + (hcir->vertex()->point() - ORIGIN);
      ++n;
    } while (++hcir != facet->facet_begin());
    pt = ORIGIN + (p - ORIGIN)/FT(n);
  }
  //
  void edge_node(Halfedge_handle edge, Point& pt) {
    Point p1 = edge->vertex()->point();
    Point p2 = edge->opposite()->vertex()->point();
    pt = Point((p1[0]+p2[0])/2, (p1[1]+p2[1])/2, (p1[2]+p2[2])/2);
  }
  //
  void vertex_node(Vertex_handle vertex, Point& pt) {
    pt = vertex->point();
  }
  //
  void border_node(Halfedge_handle edge, Point& ept, Point& /*vpt*/){
    edge_node(edge, ept);
  }
};

// ======================================================================
/// The geometry mask of Catmull-Clark subdivision 
template <class Poly>
class CatmullClark_mask_3 : public Linear_mask_3<Poly> {
public:
  typedef Poly                                         Polyhedron;

  typedef typename Polyhedron::Vertex_handle           Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
  typedef typename Polyhedron::Facet_handle            Facet_handle;

  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator 
                                            Halfedge_around_vertex_circulator;

  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;

  typedef typename Kernel::FT                          FT;
  typedef typename Kernel::Point_3                     Point;
  typedef typename Kernel::Vector_3                    Vector;

public:
  //
  void edge_node(Halfedge_handle edge, Point& pt) {
    Point p1 = edge->vertex()->point();
    Point p2 = edge->opposite()->vertex()->point();
    Point f1, f2;
    this->facet_node(edge->facet(), f1);
    this->facet_node(edge->opposite()->facet(), f2);
    pt = Point((p1[0]+p2[0]+f1[0]+f2[0])/4,
	       (p1[1]+p2[1]+f1[1]+f2[1])/4,
	       (p1[2]+p2[2]+f1[2]+f2[2])/4 );
  }
  //
  void vertex_node(Vertex_handle vertex, Point& pt) {
    Halfedge_around_vertex_circulator vcir = vertex->vertex_begin();
    int n = static_cast<int>(circulator_size(vcir));    

    FT Q[] = {0.0, 0.0, 0.0}, R[] = {0.0, 0.0, 0.0};
    Point& S = vertex->point();
    
    Point q;
    for (int i = 0; i < n; i++, ++vcir) {
      Point& p2 = vcir->opposite()->vertex()->point();
      R[0] += (S[0]+p2[0])/2;
      R[1] += (S[1]+p2[1])/2;
      R[2] += (S[2]+p2[2])/2;
      this->facet_node(vcir->facet(), q);
      Q[0] += q[0];      
      Q[1] += q[1];      
      Q[2] += q[2];
    }
    R[0] /= n;    R[1] /= n;    R[2] /= n;
    Q[0] /= n;    Q[1] /= n;    Q[2] /= n;
      
    pt = Point((Q[0] + 2*R[0] + S[0]*(n-3))/n,
	       (Q[1] + 2*R[1] + S[1]*(n-3))/n,
	       (Q[2] + 2*R[2] + S[2]*(n-3))/n );
  }
  //
  void border_node(Halfedge_handle edge, Point& ept, Point& vpt) {
    Point& ep1 = edge->vertex()->point();
    Point& ep2 = edge->opposite()->vertex()->point();
    ept = Point((ep1[0]+ep2[0])/2, (ep1[1]+ep2[1])/2, (ep1[2]+ep2[2])/2);

    Halfedge_around_vertex_circulator vcir = edge->vertex_begin();
    Point& vp1  = vcir->opposite()->vertex()->point();
    Point& vp0  = vcir->vertex()->point();
    Point& vp_1 = (--vcir)->opposite()->vertex()->point();
    vpt = Point((vp_1[0] + 6*vp0[0] + vp1[0])/8,
		(vp_1[1] + 6*vp0[1] + vp1[1])/8,
		(vp_1[2] + 6*vp0[2] + vp1[2])/8 );
  }
};


// ======================================================================
/// The geometry mask of Loop subdivision
template <class Poly>
class Loop_mask_3 : public PQQ_stencil_3<Poly> {
public:
  typedef Poly                                        Polyhedron;

  typedef typename Polyhedron::Vertex_handle           Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
  typedef typename Polyhedron::Facet_handle            Facet_handle;

  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator 
                                            Halfedge_around_vertex_circulator;

  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;

  typedef typename Kernel::FT                          FT;
  typedef typename Kernel::Point_3                     Point;
  typedef typename Kernel::Vector_3                    Vector;

public:
  //
  void edge_node(Halfedge_handle edge, Point& pt) {
    Point& p1 = edge->vertex()->point();
    Point& p2 = edge->opposite()->vertex()->point();
    Point& f1 = edge->next()->vertex()->point();
    Point& f2 = edge->opposite()->next()->vertex()->point();
      
    pt = Point((3*(p1[0]+p2[0])+f1[0]+f2[0])/8,
	       (3*(p1[1]+p2[1])+f1[1]+f2[1])/8,
	       (3*(p1[2]+p2[2])+f1[2]+f2[2])/8 );
  }
  //
  void vertex_node(Vertex_handle vertex, Point& pt) {
    Halfedge_around_vertex_circulator vcir = vertex->vertex_begin();
    size_t n = circulator_size(vcir);    

    FT R[] = {0.0, 0.0, 0.0};
    Point& S = vertex->point();
    
    for (size_t i = 0; i < n; i++, ++vcir) {
      Point& p = vcir->opposite()->vertex()->point();
      R[0] += p[0]; 	R[1] += p[1]; 	R[2] += p[2];
    }
    if (n == 6) {
      pt = Point((10*S[0]+R[0])/16, (10*S[1]+R[1])/16, (10*S[2]+R[2])/16);
    } else {
      FT Cn = (FT) (5.0/8.0 - CGAL::square(3+2*std::cos(2 * CGAL_PI/(double) n))/64.0);
      FT Sw = (double)n*(1-Cn)/Cn;
      FT W = (double)n/Cn;
      pt = Point((Sw*S[0]+R[0])/W, (Sw*S[1]+R[1])/W, (Sw*S[2]+R[2])/W);
    }
  }
  //
  //void facet_node(Facet_handle facet, Point& pt) {};
  //
  void border_node(Halfedge_handle edge, Point& ept, Point& vpt) {
    Point& ep1 = edge->vertex()->point();
    Point& ep2 = edge->opposite()->vertex()->point();
    ept = Point((ep1[0]+ep2[0])/2, (ep1[1]+ep2[1])/2, (ep1[2]+ep2[2])/2);

    Halfedge_around_vertex_circulator vcir = edge->vertex_begin();
    Point& vp1  = vcir->opposite()->vertex()->point();
    Point& vp0  = vcir->vertex()->point();
    Point& vp_1 = (--vcir)->opposite()->vertex()->point();
    vpt = Point((vp_1[0] + 6*vp0[0] + vp1[0])/8,
		(vp_1[1] + 6*vp0[1] + vp1[1])/8,
		(vp_1[2] + 6*vp0[2] + vp1[2])/8 );
  }
};


//==========================================================================
/// The setncil of the Dual-Quadrilateral-Quadrisection 
template <class Poly>
class DQQ_stencil_3 {
public:
  typedef Poly                                        Polyhedron;
  typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;

  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;

  typedef typename Kernel::FT                          FT;
  typedef typename Kernel::Point_3                     Point;
  typedef typename Kernel::Vector_3                    Vector;

public:
  //
  void corner_node(Halfedge_handle /*edge*/, Point& /*pt*/) {};
};


// ======================================================================
/// The geometry mask of Doo-Sabin subdivision
template <class Poly>
class DooSabin_mask_3 : public DQQ_stencil_3<Poly> {
public:
  typedef Poly                                        Polyhedron;
  
  typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;

  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;

  typedef typename Kernel::FT                          FT;
  typedef typename Kernel::Point_3                     Point;
  typedef typename Kernel::Vector_3                    Vector;

public:
  //
  void corner_node(Halfedge_handle he, Point& pt) {
    size_t n =  CGAL::circulator_size(he->facet()->facet_begin()); 

    Vector cv(0,0,0);
    if (n == 4) {
      cv = cv + (he->vertex()->point()-CGAL::ORIGIN)*9;
      cv = cv + (he->next()->vertex()->point()-CGAL::ORIGIN)*3;
      cv = cv + (he->next()->next()->vertex()->point()-CGAL::ORIGIN);
      cv = cv + (he->prev()->vertex()->point()-CGAL::ORIGIN)*3;
      cv = cv/16;
		} else {
			FT a;
			for (size_t k = 0; k < n; ++k, he = he->next()) {
				if (k == 0) a = (FT) ((5.0/n) + 1);
				else a = (FT) (3+2*std::cos(2*k*CGAL_PI/n))/n;
				cv = cv + (he->vertex()->point()-CGAL::ORIGIN)*a;
			}
			cv = cv/4;
		}
    pt = CGAL::ORIGIN + cv;
  }
};

// ======================================================================
// The geometry mask of Sqrt(3) subdivision
template <class Poly>
class Sqrt3_mask_3 : public Linear_mask_3<Poly> {
public:
  typedef Poly                                        Polyhedron;
  
  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;
 

  typedef typename Polyhedron::Vertex_handle           Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
  typedef typename Polyhedron::Facet_handle            Facet_handle;

  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator 
                                            Halfedge_around_vertex_circulator;

  typedef typename Kernel::FT                          FT;
  typedef typename Kernel::Point_3                     Point;
  typedef typename Kernel::Vector_3                    Vector;

public:
  //
  //void edge_node(Halfedge_handle edge, Point& pt) {}
  //
  void vertex_node(Vertex_handle vertex, Point& pt) {
    Halfedge_around_vertex_circulator vcir = vertex->vertex_begin();
    size_t n = circulator_size(vcir);

    FT a = (FT) ((4.0-2.0*std::cos(2.0*CGAL_PI/(double)n))/9.0);

    Vector cv = ((FT)(1.0-a)) * (vertex->point() - CGAL::ORIGIN);
    for (size_t i = 1; i <= n; ++i, --vcir) {
      cv = cv + (a/FT(n))*(vcir->opposite()->vertex()->point()-CGAL::ORIGIN);
    }

    pt = CGAL::ORIGIN + cv;    
  }
  //
  // TODO
  //void border_node(Halfedge_handle edge, Point& ept, Point& vpt) {}
};

} //namespace CGAL

#endif //CGAL_POLYHEDRON_SUBDIVISION_STENCILS_H_01292002
