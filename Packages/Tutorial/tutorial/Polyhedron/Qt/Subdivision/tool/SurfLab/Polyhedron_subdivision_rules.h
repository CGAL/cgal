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

SURFLAB_BEGIN_NAMESPACE

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

  typedef typename Polyhedron::Vertex_iterator         Vertex_iterator;
  typedef typename Polyhedron::Halfedge_iterator       Halfedge_iterator;
  typedef typename Polyhedron::Facet_iterator          Facet_iterator;

  typedef typename Polyhedron::Point_3                 Point;

/**@name Class Methods */
//@{
public:
  virtual void face_point_rule(Facet_iterator fitr, Point& pt) = 0;
  virtual void edge_point_rule(Halfedge_iterator eitr, Point& pt) = 0;
  virtual void vertex_point_rule(Vertex_iterator vitr, Point& pt) = 0;

  virtual void border_point_rule(Halfedge_iterator eitr, 
				 Point& ept, Point& vpt) {};
//@}
};


// ======================================================================
///
template <class _Poly>
class average_rule : public quadralize_rule<_Poly> {
  typedef _Poly                                        Polyhedron;

  typedef typename Polyhedron::Vertex_iterator         Vertex_iterator;
  typedef typename Polyhedron::Halfedge_iterator       Halfedge_iterator;
  typedef typename Polyhedron::Facet_iterator          Facet_iterator;

  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator 
                                            Halfedge_around_vertex_circulator;

  typedef typename Polyhedron::Point_3                 Point;
  typedef typename Polyhedron::FT                           FT;

/**@name Class Methods */
//@{
public:
  virtual void face_point_rule(Facet_iterator fitr, Point& pt) {
    Halfedge_around_facet_circulator hcir = fitr->facet_begin();
    int n = CGAL::circulator_size(hcir); 
    FT xyz[3] = {0, 0, 0};
    for (int v = 0; v < n; v++, ++hcir) {
      Point& p = hcir->vertex()->point();
      xyz[0] += p[0];    xyz[1] += p[1];     xyz[2] += p[2];
    }
    pt = Point(xyz[0]/n, xyz[1]/n, xyz[2]/n);
  }

  virtual void edge_point_rule(Halfedge_iterator eitr, Point& pt) {
    Point& p1 = eitr->vertex()->point();
    Point& p2 = eitr->opposite()->vertex()->point();
    pt = Point((p1[0]+p2[0])/2, (p1[1]+p2[1])/2, (p1[2]+p2[2])/2);
  }

  virtual void vertex_point_rule(Vertex_iterator vitr, Point& pt) {
    pt = vitr->point();
  }

  virtual void border_point_rule(Halfedge_iterator eitr, 
				 Point& ept, Point& vpt){
    edge_point_rule(eitr, ept);
  }
//@}
};

// ======================================================================
///
template <class _Poly>
class CatmullClark_rule : public average_rule<_Poly> {

public:
  typedef _Poly                                           Polyhedron;
private:
  typedef typename Polyhedron::Vertex_iterator         Vertex_iterator;
  typedef typename Polyhedron::Halfedge_iterator       Halfedge_iterator;
  typedef typename Polyhedron::Facet_iterator          Facet_iterator;

  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator 
                                            Halfedge_around_vertex_circulator;

  typedef typename Polyhedron::Point_3                 Point;
  typedef typename Polyhedron::FT                           FT;

/**@name Class Methods */
//@{
public:
  void edge_point_rule(Halfedge_iterator eitr, Point& pt) {
    Point& p1 = eitr->vertex()->point();
    Point& p2 = eitr->opposite()->vertex()->point();
    Point f1, f2;
    face_point_rule(Facet_iterator(eitr->facet()), f1);
    face_point_rule(Facet_iterator(eitr->opposite()->facet()), f2);
    pt = Point((p1[0]+p2[0]+f1[0]+f2[0])/4,
	       (p1[1]+p2[1]+f1[1]+f2[1])/4,
	       (p1[2]+p2[2]+f1[2]+f2[2])/4 );
  }

  void vertex_point_rule(Vertex_iterator vitr, Point& pt) {
    Halfedge_around_vertex_circulator vcir = vitr->vertex_begin();
    int n = CGAL::circulator_size(vcir);    

    float Q[] = {0.0, 0.0, 0.0}, R[] = {0.0, 0.0, 0.0};
    Point& S = vitr->point();
    
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

  void border_point_rule(Halfedge_iterator eitr, Point& ept, Point& vpt) {
    Point& ep1 = eitr->vertex()->point();
    Point& ep2 = eitr->opposite()->vertex()->point();
    ept = Point((ep1[0]+ep2[0])/2, (ep1[1]+ep2[1])/2, (ep1[2]+ep2[2])/2);

    Halfedge_around_vertex_circulator vcir = eitr->vertex_begin();
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
  typedef _Poly                                           Polyhedron;

private:

  typedef typename Polyhedron::Vertex_iterator         Vertex_iterator;
  typedef typename Polyhedron::Halfedge_iterator       Halfedge_iterator;
  typedef typename Polyhedron::Facet_iterator          Facet_iterator;

  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator 
                                            Halfedge_around_vertex_circulator;

  typedef typename Polyhedron::Point_3                 Point;
  typedef typename Polyhedron::FT                           FT;

/**@name Class Methods */
//@{
public:
  void edge_point_rule(Halfedge_iterator eitr, Point& pt) {
    Point& p1 = eitr->vertex()->point();
    Point& p2 = eitr->opposite()->vertex()->point();
    Point& f1 = eitr->next()->vertex()->point();
    Point& f2 = eitr->opposite()->next()->vertex()->point();
      
    pt = Point((3*(p1[0]+p2[0])+f1[0]+f2[0])/8,
	       (3*(p1[1]+p2[1])+f1[1]+f2[1])/8,
	       (3*(p1[2]+p2[2])+f1[2]+f2[2])/8 );
  }

  void vertex_point_rule(Vertex_iterator vitr, Point& pt) {
    Halfedge_around_vertex_circulator vcir = vitr->vertex_begin();
    int n = CGAL::circulator_size(vcir);    

    float R[] = {0.0, 0.0, 0.0};
    Point& S = vitr->point();
    
    for (int i = 0; i < n; i++, ++vcir) {
      Point& p = vcir->opposite()->vertex()->point();
      R[0] += p[0]; 	R[1] += p[1]; 	R[2] += p[2];
    }
    if (n == 6) {
      pt = Point((10*S[0]+R[0])/16, (10*S[1]+R[1])/16, (10*S[2]+R[2])/16);
    } else {
      FT Cn = (FT)(5.0/8.0 - std::sqrt(3+2*std::cos(6.283/n))/64.0);
      FT Sw = n*(1-Cn)/Cn;
      FT W = n/Cn;
      pt = Point((Sw*S[0]+R[0])/W, (Sw*S[1]+R[1])/W, (Sw*S[2]+R[2])/W);
    }
  }

  void face_point_rule(Facet_iterator fitr, Point& pt) {};

  void border_point_rule(Halfedge_iterator eitr, Point& ept, Point& vpt) {
    Point& ep1 = eitr->vertex()->point();
    Point& ep2 = eitr->opposite()->vertex()->point();
    ept = Point((ep1[0]+ep2[0])/2, (ep1[1]+ep2[1])/2, (ep1[2]+ep2[2])/2);

    Halfedge_around_vertex_circulator vcir = eitr->vertex_begin();
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
/**
 * All rule (for Polyhedron_subdivision::dualize_polyhedron<RULE>(...))
 * should be inheritated from dualize_rule and implement point_rule(...) 
 * functions.
 */
template <class _Poly>
class dualize_rule {
public:
  typedef _Poly                                           Polyhedron;
  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Point_3                 Point;
  typedef typename Polyhedron::FT                          FT;

public:
  virtual void point_rule(Halfedge_around_facet_circulator cir, FT* xyz) = 0;
};


// ======================================================================
///
template <class _Poly>
class DooSabin_rule : public dualize_rule<_Poly> {
public:
  typedef _Poly                                           Polyhedron;
private:
  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Point_3                 Point;
  typedef typename Polyhedron::FT                           FT;

public:
  void point_rule(Halfedge_around_facet_circulator cir, FT* xyz) {
    int n =  CGAL::circulator_size(cir); 

    FT x = 0, y = 0, z = 0;
    if (n == 4) {
      Point fp = cir->vertex()->point(); ++cir;
      x = 9*fp.x(); y = 9*fp.y(); z = 9*fp.z();
      fp = cir->vertex()->point(); ++cir;
      x += 3*fp.x(); y += 3*fp.y(); z += 3*fp.z();
      fp = cir->vertex()->point(); ++cir;
      x += fp.x(); y += fp.y(); z += fp.z();
      fp = cir->vertex()->point();
      x += 3*fp.x(); y += 3*fp.y(); z += 3*fp.z();
      x /= 16; y /= 16; z /= 16;
    } else {
      FT a;
      for (int k = 0; k < n; ++k) {
	if (k == 0) a = ((FT)5/n) + 1;
	else a = (FT)((3+2*std::cos(2*k*3.141593/n))/n);
	Point& fp = cir->vertex()->point(); ++cir;
	x += a*fp.x(); y += a*fp.y(); z += a*fp.z();
      }
      x /= 4; y /= 4; z /= 4;
    }
    xyz[0] = x;    xyz[1] = y;    xyz[2] = z;
  }
};


SURFLAB_END_NAMESPACE

#endif //_POLYHEDRON_SUBDIVISION_RULES_H_01292002
