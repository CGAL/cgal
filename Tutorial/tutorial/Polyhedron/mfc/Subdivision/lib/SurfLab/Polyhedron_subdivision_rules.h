// ======================================================================
//
// Copyright (c) 2002 SurfLab of CISE of University of Florida
//
// File          : libs/src/cgalExt/Polyhedron_subdivision_rules.h
// Description   : Provides the subdivision rules used to template the 
//                 subdivision functions.
// Creation_date : 29 Jan 2002
// Author(s)     : Le-Jeng Shiue <sle-jeng@cise.ufl.edu>
//
// ======================================================================

// $Id$

#ifndef _POLYHEDRON_SUBDIVISION_RULES_H_01292002
#define _POLYHEDRON_SUBDIVISION_RULES_H_01292002

#include <CGAL/circulator.h>
#include <CGAL/Vector_3.h>

// ======================================================================
template <class _Poly>
class quadralize_rule {
public:
  typedef _Poly                                         Polyhedron;

  typedef typename Polyhedron::Vertex_handle           Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
  typedef typename Polyhedron::Facet_handle            Facet_handle;
  typedef typename Polyhedron::Point_3                 Point;

/**@name Class Methods */
//@{
public:
  void face_point_rule(Facet_handle, Point&) {};
  void edge_point_rule(Halfedge_handle, Point&) {};
  void vertex_point_rule(Vertex_handle, Point&) {};

  void border_point_rule(Halfedge_handle, Point&, Point&) {};
//@}
};


// ======================================================================
///
template <class _Poly>
class average_rule : public quadralize_rule<_Poly> {
public:
  typedef _Poly                                         Polyhedron;

  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;

  typedef typename Polyhedron::Vertex_handle           Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
  typedef typename Polyhedron::Facet_handle            Facet_handle;

  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator 
                                            Halfedge_around_vertex_circulator;

  typedef typename Polyhedron::Point_3                 Point;
  typedef typename Kernel::FT                          FT;

/**@name Class Methods */
//@{
public:
  void face_point_rule(Facet_handle facet, Point& pt) {
#ifdef USE_CGAL_POINT
    Halfedge_around_facet_circulator hcir = facet->facet_begin();
    Vector vec = hcir->vertex()->point() - CGAL::ORIGIN;
    ++hcir;
    do {
       vec = vec + hcir->vertex()->point();
    } while (++hcir != facet->facet_begin());
    pt = CGAL::ORIGIN + vec/circulator_size(hcir);
#else
    Halfedge_around_facet_circulator hcir = facet->facet_begin();
    int n = 0;
    FT p[] = {0,0,0};
    do {
      Point t = hcir->vertex()->point();
      p[0] += t[0], p[1] += t[1], p[2] += t[2]; 
      ++n;
    } while (++hcir != facet->facet_begin());
    pt = Point(p[0]/n, p[1]/n, p[2]/n);
#endif
  }

  void edge_point_rule(Halfedge_handle edge, Point& pt) {
    Point p1 = edge->vertex()->point();
    Point p2 = edge->opposite()->vertex()->point();
    pt = Point((p1[0]+p2[0])/2, (p1[1]+p2[1])/2, (p1[2]+p2[2])/2);
  }

  void vertex_point_rule(Vertex_handle vertex, Point& pt) {
    pt = vertex->point();
  }

  void border_point_rule(Halfedge_handle edge, Point& ept, Point& vpt) {
    edge_point_rule(edge, ept);
  }
//@}
};

// ======================================================================
///
template <class _Poly>
class CatmullClark_rule : public average_rule<_Poly> {
public:
  typedef _Poly                                        Polyhedron;

  typedef typename Polyhedron::Vertex_handle           Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
  typedef typename Polyhedron::Facet_handle            Facet_handle;

  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator 
                                            Halfedge_around_vertex_circulator;

  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;

  typedef typename Polyhedron::Point_3                 Point;

/**@name Class Methods */
//@{
public:
  void edge_point_rule(Halfedge_handle edge, Point& pt) {
    Point p1 = edge->vertex()->point();
    Point p2 = edge->opposite()->vertex()->point();
    Point f1, f2;
    face_point_rule(edge->facet(), f1);
    face_point_rule(edge->opposite()->facet(), f2);
    pt = Point((p1[0]+p2[0]+f1[0]+f2[0])/4,
	       (p1[1]+p2[1]+f1[1]+f2[1])/4,
	       (p1[2]+p2[2]+f1[2]+f2[2])/4 );
  }

  void vertex_point_rule(Vertex_handle vertex, Point& pt) {
    Halfedge_around_vertex_circulator vcir = vertex->vertex_begin();
    int n = circulator_size(vcir);    

    float Q[] = {0.0, 0.0, 0.0}, R[] = {0.0, 0.0, 0.0};
    Point& S = vertex->point();
    
    Point q;
    for (int i = 0; i < n; i++, ++vcir) {
      Point& p2 = vcir->opposite()->vertex()->point();
      R[0] += (S[0]+p2[0])/2;
      R[1] += (S[1]+p2[1])/2;
      R[2] += (S[2]+p2[2])/2;
      face_point_rule(vcir->facet(), q);
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

  void border_point_rule(Halfedge_handle edge, Point& ept, Point& vpt) {
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
//@}
};


// ======================================================================
///
template <class _Poly>
class Loop_rule : public quadralize_rule<_Poly> {
public:
  typedef _Poly                                        Polyhedron;

  typedef typename Polyhedron::Vertex_handle           Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
  typedef typename Polyhedron::Facet_handle            Facet_handle;

  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator 
                                            Halfedge_around_vertex_circulator;

  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;

  typedef typename Polyhedron::Point_3                 Point;

/**@name Class Methods */
//@{
public:
  void edge_point_rule(Halfedge_handle edge, Point& pt) {
    Point& p1 = edge->vertex()->point();
    Point& p2 = edge->opposite()->vertex()->point();
    Point& f1 = edge->next()->vertex()->point();
    Point& f2 = edge->opposite()->next()->vertex()->point();
      
    pt = Point((3*(p1[0]+p2[0])+f1[0]+f2[0])/8,
	       (3*(p1[1]+p2[1])+f1[1]+f2[1])/8,
	       (3*(p1[2]+p2[2])+f1[2]+f2[2])/8 );
  }

  void vertex_point_rule(Vertex_handle vertex, Point& pt) {
    Halfedge_around_vertex_circulator vcir = vertex->vertex_begin();
    int n = circulator_size(vcir);    

    float Q[] = {0.0, 0.0, 0.0}, R[] = {0.0, 0.0, 0.0};
    Point& S = vertex->point();
    
    for (int i = 0; i < n; i++, ++vcir) {
      Point& p = vcir->opposite()->vertex()->point();
      R[0] += p[0]; 	R[1] += p[1]; 	R[2] += p[2];
    }
    if (n == 6) {
      pt = Point((10*S[0]+R[0])/16, (10*S[1]+R[1])/16, (10*S[2]+R[2])/16);
    } else {
      double Cn = 5.0/8.0 - std::sqrt(3+2*std::cos(6.283/n))/64.0;
      double Sw = n*(1-Cn)/Cn;
      double W = n/Cn;
      pt = Point((Sw*S[0]+R[0])/W, (Sw*S[1]+R[1])/W, (Sw*S[2]+R[2])/W);
    }
  }

  void border_point_rule(Halfedge_handle edge, Point& ept, Point& vpt) {
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
//@}
};


//==========================================================================
template <class _Poly>
class dualize_rule {
public:
  typedef _Poly                                        Polyhedron;
  typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;

  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;
  typedef typename Polyhedron::Point_3                 Point;

public:
  void point_rule(Halfedge_handle edge, Point& pt) {};
};


// ======================================================================
///
template <class _Poly>
class DooSabin_rule : public dualize_rule<_Poly> {
public:
  typedef _Poly                                        Polyhedron;
  typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;

  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;
  typedef typename Polyhedron::Point_3                 Point;
  typedef CGAL::Vector_3<Kernel>                       Vector;

public:
  void point_rule(Halfedge_handle he, Point& pt) {
    int n =  CGAL::circulator_size(he->facet()->facet_begin()); 

    Vector cv(0,0,0);
    if (n == 4) {
      cv = cv + (he->vertex()->point()-CGAL::ORIGIN)*9;
      cv = cv + (he->next()->vertex()->point()-CGAL::ORIGIN)*3;
      cv = cv + (he->next()->next()->vertex()->point()-CGAL::ORIGIN);
      cv = cv + (he->prev()->vertex()->point()-CGAL::ORIGIN)*3;
      cv = cv/16;
    } else {
      double a;
      for (int k = 0; k < n; ++k, he = he->next()) {
	if (k == 0) a = ((double)5/n) + 1;
	else a = (3+2*std::cos(2*k*3.141593/n))/n;
	cv = cv + (he->vertex()->point()-CGAL::ORIGIN)*a;
      }
      cv = cv/4;
    }
    pt = CGAL::ORIGIN + cv;
  }
};

#endif //_POLYHEDRON_SUBDIVISION_RULES_H_01292002
