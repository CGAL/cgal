// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : include/CGAL/Nef_3/Infimaximal_box.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Peter Hachenberger    <hachenberger@mpi-sb.mpg.de>
// maintainer    : Peter Hachenberger    <hachenberger@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// Infimaximal_box.h    answers queries related to the infimaximal box no matter
//                      if it exists (extended kernel) or not (otherwise)
// ============================================================================
#ifndef CGAL_INFIMAXIMAL_BOX_H
#define CGAL_INFIMAXIMAL_BOX_H

#undef _DEBUG
#define _DEBUG 191
#include <CGAL/Nef_3/debug.h>

#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Extended_homogeneous_3.h>

#include <CGAL/Nef_S2/Sphere_point.h>

CGAL_BEGIN_NAMESPACE

template <class T, class Kernel>
class Infimaximal_box {

  typedef typename Kernel::RT               NT;
  typedef Kernel                            Standard_kernel;
  typedef typename Kernel::Point_3          Point_3;
  typedef Point_3                           Standard_point;

 public:
  static bool is_standard(Point_3& p) {
    return true;
  }

  static Point_3 simplify(Point_3& p) {
    return p;
  }

  static Point_3 box_point(Point_3& p, NT d=10000) {
    return p;
  }

  static Standard_point standard_point(Point_3 p, NT d=1) {
    return p;
  }

  static int degree(const typename Kernel::RT& n) {
    return 0;
  }

  template <typename SNC_decorator_>
    static Point_3 target_for_ray_shot(SNC_decorator_& deco, Point_3& p) {
    return deco.vertices_begin()->point();
  }

};

template <class Kernel>
class Infimaximal_box<Tag_true, Kernel > {

  typedef typename Kernel::RT               RT;
  typedef typename Kernel::RT::NT           NT;
  typedef typename Kernel::Standard_kernel  Standard_kernel;
  typedef typename Standard_kernel::Point_3 Standard_point;
  typedef typename Kernel::Point_3          Point_3;
  typedef typename Kernel::Plane_3          Plane_3;
  typedef typename Kernel::Vector_3         Vector_3;

  enum Boundary { EXCLUDED=0, INCLUDED=1 };
  
 public:
  static bool is_standard(Point_3& p) {
    return Kernel::is_standard(p);
  }

  static Point_3 simplify(Point_3& p) {
    CGAL_assertion(p.hw().degree() == 0);
    int deg = p.hx().degree() > p.hy().degree() 
      ? p.hx().degree() 
      : p.hy().degree();
    deg = p.hz().degree() > deg 
      ? p.hz().degree() 
      : deg;
    return Point_3(p.hx()(deg),p.hy()(deg),p.hz()(deg),p.hw()[0]);
  }

  static int degree(const RT& n) {
    return n.degree();
  }

  template <typename SNC_decorator_>
  static Point_3 target_for_ray_shot(SNC_decorator_& deco, Point_3 p) {
    return Kernel::epoint(0, p.hx()[0], 0, p.hy()[0], p.hw()[0], 0, p.hw()[0]);
  }

  static Standard_point standard_point(Point_3 p, NT d=1) {
    return Standard_point(p.hx().eval_at(d),
			  p.hy().eval_at(d),
			  p.hz().eval_at(d),
			  p.hw().eval_at(1));
  }


  static Point_3 box_point(Point_3 p, NT d=10000) {
    return Point_3(p.hx().eval_at(d),
		   p.hy().eval_at(d),
		   p.hz().eval_at(d),
		   p.hw().eval_at(1));
  }

  static Point_3 create_extended_point(NT x, NT y, NT z) {
    return Kernel::epoint(x,0,y,0,z,0,1);
  }


  /*
  template <typename SNC_structure_>
  static void create_vertices_of_box_with_plane(const Plane_3& h, Boundary b, SNC_structure_& snc) {

    typedef typename SNC_structure_::Sphere_point  Sphere_point;

    Point_3 loc(-h.d(),0,0,h.a());
    Vector_3 orth = h.orthogonal_vector();
    
    NT orth_coords[3];
    orth_coords[0] = orth.hx()[0];
    orth_coords[1] = orth.hy()[0];
    orth_coords[2] = orth.hz()[0];

    int add_corners = 0;
    while(orth_coords[add_corners] == 0) add_corners++;
    CGAL_assertion(add_corners < 3);

    std::list<Point_3> points;
    for(int dir=0; dir<3;++dir) {

      NT cnst[3];
      for(int i=0; i<3;i++)
	cnst[i] = (i==dir? -h.d()[0] : 0);

      NT cross[4][4];
      cross[0][dir] = -orth_coords[(dir+1)%3]-orth_coords[(dir+2)%3];
      cross[1][dir] =  orth_coords[(dir+1)%3]-orth_coords[(dir+2)%3];
      cross[2][dir] =  orth_coords[(dir+1)%3]+orth_coords[(dir+2)%3];  
      cross[3][dir] = -orth_coords[(dir+1)%3]+orth_coords[(dir+2)%3];
  
      for(int i=0;i<4;++i)
	cross[i][3] = orth_coords[dir];

      cross[0][(dir+1)%3] = cross[3][(dir+1)%3] =  orth_coords[dir];
      cross[1][(dir+1)%3] = cross[2][(dir+1)%3] = -orth_coords[dir];
      
      cross[0][(dir+2)%3] = cross[1][(dir+2)%3] =  orth_coords[dir];
      cross[2][(dir+2)%3] = cross[3][(dir+2)%3] = -orth_coords[dir];

      for(int i=0; i<4; ++i)
	if(CGAL_NTS abs(RT(cnst[dir],cross[i][dir])) < CGAL_NTS abs(RT(0,orth_coords[dir])) ||
	   (CGAL_NTS abs(RT(cnst[dir],cross[i][dir])) == CGAL_NTS abs(RT(0,orth_coords[dir])) && dir == add_corners))
	  points.push_back(Kernel::epoint(cross[i][0],cnst[0],cross[i][1],cnst[1],cross[i][2],cnst[2],cross[i][3]));
      
    }

    for(int i=0;i<2;i++)
      orth_coords[i] = CGAL_NTS abs(orth_coords[i]);

    int max = 0;
    if(orth_coords[1] > orth_coords[0])
      max = 1;
    if(orth_coords[2] > orth_coords[max])
      max = 2;   

    int min = 0;
    if(orth_coords[1] < orth_coords[0])
      min = 1;
    if(orth_coords[2] < orth_coords[min])
      min = 2;   

    SNC_constructor<SNC_structure_> C(snc);
    points.sort(circle_lt<Point_3>(max));

    typename std::list<Point_3>::const_iterator p,prev,next;
    for(p=points.begin();p!=points.end();p++)
      TRACEN(*p);

    for(p=points.begin();p!=points.end();p++){

      if(p==points.begin()) prev = --points.end();
      else { prev = p; prev--;}
      if(p==--points.end()) next=points.begin();
      else {next = p; ++next;}
      TRACEN("points " << *prev << "           " << *p << "      " << *next);

      Vector_3 v= *prev - *p;
      Sphere_point sp1(v);
      sp1 = sp1.normalized();
      CGAL_assertion(sp1.hx().degree() == 0);
      CGAL_assertion(sp1.hy().degree() == 0);
      CGAL_assertion(sp1.hz().degree() == 0);
      CGAL_assertion(sp1.hw().degree() == 0);

      v= *next - *p;
      Sphere_point sp2(v);
      sp2 = sp2.normalized();
      CGAL_assertion(sp2.hx().degree() == 0);
      CGAL_assertion(sp2.hy().degree() == 0);
      CGAL_assertion(sp2.hz().degree() == 0);
      CGAL_assertion(sp2.hw().degree() == 0);

      TRACEN("sps " << sp1 << "     " << sp2);
      TRACEN(orth_coords[min] << "|" << orth_coords[(min+1)%3] << "|" << orth_coords[(min+2)%3]);

      if(orth_coords[min]==0 && orth_coords[(min+1)%3] == orth_coords[(min+2)%3] && h.d() == 0) 
	C.create_degenerate_corner_frame_point(*p,sp1,sp2,min, max, (b==INCLUDED));
      else if(CGAL_NTS abs(p->hx()) == CGAL_NTS abs(p->hy()) && CGAL_NTS abs(p->hz()) == CGAL_NTS abs(p->hy()))
	C.create_corner_frame_point(*p,sp1,sp2,max,(b==INCLUDED));
      else
	C.create_frame_point(*p,sp1,sp2,h,(b==INCLUDED));
    }

    RT sum= h.a()+h.b()+h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      TRACEN(sum); 
      C.create_extended_box_corner( 1, 1, 1, (sum<0 || (sum == 0 && h.d()<0)));
    }
    sum=-h.a()+h.b()+h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      TRACEN(sum); 
      C.create_extended_box_corner(-1, 1, 1, (sum<0 || (sum == 0 && h.d()<0)));
    }
    sum= h.a()-h.b()+h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      TRACEN(sum); 
      C.create_extended_box_corner( 1,-1, 1, (sum<0 || (sum == 0 && h.d()<0)));
    }
    sum=-h.a()-h.b()+h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      TRACEN(sum); 
      C.create_extended_box_corner(-1,-1, 1, (sum<0 || (sum == 0 && h.d()<0)));
    }
    sum= h.a()+h.b()-h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      TRACEN(sum); 
      C.create_extended_box_corner( 1, 1,-1, (sum<0 || (sum == 0 && h.d()<0)));
    }
    sum=-h.a()+h.b()-h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      TRACEN(sum); 
      C.create_extended_box_corner(-1, 1,-1, (sum<0 || (sum == 0 && h.d()<0)));
    }
    sum= h.a()-h.b()-h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      TRACEN(sum); 
      C.create_extended_box_corner( 1,-1,-1, (sum<0 || (sum == 0 && h.d()<0)));
    }
    sum=-h.a()-h.b()-h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      TRACEN(sum); 
      C.create_extended_box_corner(-1,-1,-1, (sum<0 || (sum == 0 && h.d()<0)));
    }
  }
  */

};

CGAL_END_NAMESPACE
#endif // CGAL_INFIMAXIMAL_BOX_H
