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
/**
 * All rule (for Polyhedron_subdivision::quadralize_polyhedron<RULE>(...))
 * should be inheritated from quadralize_rule and implement these three rule 
 * functions.
 */
template <class _Poly>
class quadralize_rule {
public:
  typedef _Poly                                        Polyhedron;

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
  typedef _Poly                                        Polyhedron;

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
//     Halfedge_around_facet_circulator hcir = facet->facet_begin();
//     Vector vec = hcir->vertex()->point() - CGAL::ORIGIN;
//     ++hcir;
//     do {
//       vec = vec + hcir->vertex()->point();
//     } while (++hcir != facet->facet_begin());
//     pt = CGAL::ORIGIN + vec/circulator_size(hcir);

    Halfedge_around_facet_circulator hcir = facet->facet_begin();
    int n = 0;
    FT p[] = {0,0,0};
    do {
      Point t = hcir->vertex()->point();
      p[0] += t[0], p[1] += t[1], p[2] += t[2]; 
      ++n;
    } while (++hcir != facet->facet_begin());
    pt = Point(p[0]/n, p[1]/n, p[2]/n);
  }

  void edge_point_rule(Halfedge_handle edge, Point& pt) {
    Point p1 = edge->vertex()->point();
    Point p2 = edge->opposite()->vertex()->point();
    pt = Point((p1[0]+p2[0])/2, (p1[1]+p2[1])/2, (p1[2]+p2[2])/2);
  }

  void vertex_point_rule(Vertex_handle vertex, Point& pt) {
    pt = vertex->point();
  }

  void border_point_rule(Halfedge_handle edge, Point& ept, Point& vpt){
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
  //typedef typename Kernel::FT                          FT;

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

  //void face_point_rule(Facet_handle facet, Point& pt) {};

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
class QT43_rule : public average_rule<_Poly> {
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
  typedef CGAL::Vector_3<Kernel>                       Vector;
  typedef typename Kernel::FT                          FT;

/**@name Class Methods */
//@{
public:
  void edge_point_rule(Halfedge_handle edge, Point& pt) {
    int f1deg = CGAL::circulator_size(edge->facet()->facet_begin());
    int f2deg = CGAL::circulator_size(edge->opposite()->facet()->facet_begin());
    
    if (f1deg > 3 && f2deg > 3) qe_rule(edge, pt);
    else if (f1deg == 3 && f2deg == 3) te_rule(edge, pt);
    else if (f1deg == 4 && f2deg == 3) qte_rule(edge, pt);
    else if (f1deg == 3 && f2deg == 4) qte_rule(edge->opposite(), pt);
  }

  void vertex_point_rule(Vertex_handle vertex, Point& pt) {
    Halfedge_around_vertex_circulator vcir = vertex->vertex_begin();
    int n = circulator_size(vcir);    
    
    int nt = 0, nq = 0;
    int* fdeg = new int[n];
    for (int i = 0; i < n; ++i, --vcir) {
      fdeg[i] = circulator_size(vcir->facet()->facet_begin());
      if (fdeg[i] == 3) ++nt;
      else              ++nq;
    }
      
    bool reg = false;;

    if (n <= 6 && n > 3) { // possible regular node
      if (n == 6 && nt == n) { tv_rule(vertex, pt); reg = true; }
      else if (n == 4 && nq == n) { qv_rule(vertex, pt); reg = true; }
      else if (n == 5 && nq == 2 && nt == 3) {
	Halfedge_around_vertex_circulator c = vertex->vertex_begin();
	for (int i = 0; i < 5; ++i, --c) 
	  if (fdeg[i] == 4 && fdeg[(i+1) % 5] == 4) {
	    qtv_rule(c, pt);
	    reg = true;  break;
	  }
      }
    }

    if (!reg) eov_rule(n, vcir, pt);

    delete[] fdeg;
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

private:
  ///
  void qe_rule(Halfedge_handle edge, Point& pt) {
    // p3 --- p1 --- p6 
    //        |
    //        |
    // p4 --- p2 --- p5
    Point& p1 = edge->vertex()->point();
    Point& p2 = edge->opposite()->vertex()->point();
    Point& p3 = edge->next()->vertex()->point();
    Point& p4 = edge->prev()->prev()->vertex()->point();
    Point& p5 = edge->opposite()->next()->vertex()->point();
    Point& p6 = edge->opposite()->prev()->prev()->vertex()->point();
    
    pt = CGAL::ORIGIN + ((p1-CGAL::ORIGIN + p2)*6.0 + p3+p4+p5+p6)/16.0;
  }

  ///
  void te_rule(Halfedge_handle edge, Point& pt) {
    //    p3
    //  /    \
    // p1 -- p2
    //  \    /
    //    p4
    Point& p1 = edge->vertex()->point();
    Point& p2 = edge->opposite()->vertex()->point();
    Point& p3 = edge->opposite()->next()->vertex()->point();
    Point& p4 = edge->next()->vertex()->point();
    
    pt = CGAL::ORIGIN + ((p1-CGAL::ORIGIN + p2)*3.0 + p3+p4)/8.0;
  }
  void qte_rule(Halfedge_handle edge, Point& pt) {
    //      p5
    //    /    \
    //   p1 -- p2
    //   |      |
    //   p3     p4
    Point& p1 = edge->vertex()->point();
    Point& p2 = edge->opposite()->vertex()->point();
    Point& p3 = edge->next()->vertex()->point();
    Point& p4 = edge->prev()->opposite()->vertex()->point();
    Point& p5 = edge->opposite()->next()->vertex()->point();
    
    pt = CGAL::ORIGIN + ((p1-CGAL::ORIGIN + p2)*6.0 + p3+p4+(p5-CGAL::ORIGIN)*2)/16.0;
  }

  ///
  void qv_rule(Vertex_handle vertex, Point& pt) {
    //        p2
    //        |
    // p3 --- p1 --- p5
    //        |
    //        p4
    Point& p1 = vertex->point();

    Halfedge_around_vertex_circulator vcir = vertex->vertex_begin();
    Point& p2 = vcir->opposite()->vertex()->point(); --vcir;
    Point& p3 = vcir->opposite()->vertex()->point(); --vcir;
    Point& p4 = vcir->opposite()->vertex()->point(); --vcir;
    Point& p5 = vcir->opposite()->vertex()->point(); 
    
    pt = CGAL::ORIGIN + ((p1-CGAL::ORIGIN)*4.0 + p2+p3+p4+p5)/8.0;
  }

  ///
  void tv_rule(Vertex_handle vertex, Point& pt) {
    Point& p1 = vertex->point();
    
    Halfedge_around_vertex_circulator vcir = vertex->vertex_begin();
    Point& p2 = vcir->opposite()->vertex()->point(); --vcir;
    Point& p3 = vcir->opposite()->vertex()->point(); --vcir;
    Point& p4 = vcir->opposite()->vertex()->point(); --vcir;
    Point& p5 = vcir->opposite()->vertex()->point(); --vcir;
    Point& p6 = vcir->opposite()->vertex()->point(); --vcir;
    Point& p7 = vcir->opposite()->vertex()->point();

    pt = CGAL::ORIGIN + ((p1-CGAL::ORIGIN)*10.0 + p2+p3+p4+p5+p6+p7)/16.0;
  }

  ///
  void qtv_rule(Halfedge_around_vertex_circulator vcir, Point& pt) {
    //        p6
    //        |   p5
    //   p2 - p1
    //        |   p4
    //        p3
    Point& p1 = vcir->vertex()->point(); ++vcir;
    Point& p2 = vcir->opposite()->vertex()->point(); --vcir;
    Point& p3 = vcir->opposite()->vertex()->point(); --vcir;
    Point& p4 = vcir->opposite()->vertex()->point(); --vcir;
    Point& p5 = vcir->opposite()->vertex()->point(); --vcir;
    Point& p6 = vcir->opposite()->vertex()->point(); 

    pt = CGAL::ORIGIN + ((p1-CGAL::ORIGIN)*18.0 + (p2-CGAL::ORIGIN)*4.0 + 
			 (p3-CGAL::ORIGIN + p6)*3.0 + (p4-CGAL::ORIGIN + p5)*2.0)/32.0;
  }

  ///
  void eov_rule(int n, Halfedge_around_vertex_circulator vcir, Point& pt) {
    FT a = 3.0/5.0;
    if (n < 7) a = (1.0+cos(2.0*3.1415926/n)*0.5)*0.5;
    if (n == 3) a = 0.25;

    FT outlet_w = (1.0-a)/n, center_w = a;
    
    Vector cv = (vcir->vertex()->point() - CGAL::ORIGIN) * center_w;
    for (int i = 1; i <= n; ++i, --vcir) {
      cv = cv + (vcir->opposite()->vertex()->point()-CGAL::ORIGIN)*outlet_w;
    }

    pt = CGAL::ORIGIN + cv;
    
    // NOTE: this is not a complete implementation of 4-3 subdivsion!
    //       To simplify the code, disk projection and blending are not implemented.
  }
};




//==========================================================================
/**
 * All rule (for Polyhedron_subdivision::dualize_polyhedron<RULE>(...))
 * should be inheritated from dualize_rule and implement point_rule(...) 
 * functions.
 */
template <class _Poly>
class dualize_rule {
public:
  typedef _Poly                                         Polyhedron;
  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;
  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Point_3                 Point;

public:
  void point_rule(Halfedge_around_facet_circulator cir, Point& pt) {};
};


// ======================================================================
///
template <class _Poly>
class DooSabin_rule : public dualize_rule<_Poly> {
public:
  typedef _Poly                                         Polyhedron;
  
  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;

  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Point_3                 Point;
  typedef CGAL::Vector_3<Kernel>                        Vector;
  //typedef typename Kernel::FT                          FT;

public:
  void point_rule(Halfedge_around_facet_circulator cir, Point& pt) {
    int n =  CGAL::circulator_size(cir); 

    Vector cv(0,0,0), t;
    if (n == 4) {
      cv = cv + (cir->vertex()->point()-CGAL::ORIGIN)*9;
      cv = cv + ((++cir)->vertex()->point()-CGAL::ORIGIN)*3;
      cv = cv + ((++cir)->vertex()->point()-CGAL::ORIGIN);
      cv = cv + ((++cir)->vertex()->point()-CGAL::ORIGIN)*3;
      cv = cv/16;
    } else {
      double a;
      for (int k = 0; k < n; ++k, ++cir) {
	if (k == 0) a = ((double)5/n) + 1;
	else a = (3+2*std::cos(2*k*3.141593/n))/n;
	cv = cv + (cir->vertex()->point()-CGAL::ORIGIN)*a;
      }
      cv = cv/4;
    }
    pt = CGAL::ORIGIN + cv;
  }
};

// ======================================================================
///
template <class _Poly>
class Sqrt3_rule : public average_rule<_Poly> {
public:
  typedef _Poly                                        Polyhedron;
  
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
  typedef CGAL::Vector_3<Kernel>                       Vector;
  typedef typename Kernel::FT                          FT;

/**@name Class Methods */
//@{
public:
  //void edge_point_rule(Halfedge_handle edge, Point& pt) {
  //}

  void vertex_point_rule(Vertex_handle vertex, Point& pt) {
    Halfedge_around_vertex_circulator vcir = vertex->vertex_begin();
    int n = circulator_size(vcir);
    const double pi = 3.1415926;

    FT a = (4.0-2.0*cos(2.0*pi/n))/9.0;

    Vector cv = (1.0-a) * (vertex->point() - CGAL::ORIGIN);
    for (int i = 1; i <= n; ++i, --vcir) {
      cv = cv + (a/n)*(vcir->opposite()->vertex()->point()-CGAL::ORIGIN);
    }

    pt = CGAL::ORIGIN + cv;    
  }

  //void border_point_rule(Halfedge_handle edge, Point& ept, Point& vpt) {
  //}
//@}
};

#endif //_POLYHEDRON_SUBDIVISION_RULES_H_01292002
