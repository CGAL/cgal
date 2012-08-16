// Copyright (c) 1997-2000  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Thomas Herrmann, Lutz Kettner

#ifndef CGAL_WIDTH_3_H
#define CGAL_WIDTH_3_H

#include <CGAL/basic.h>
#include <cstdlib>
#include <iostream>
#include <CGAL/convex_hull_3.h>

#include <CGAL/Polyhedron_3.h> 
#include <CGAL/HalfedgeDS_list.h>

#include <CGAL/assertions.h>
#include <CGAL/Width_polyhedron_3.h>
#include <CGAL/width_assertions.h>

namespace CGAL {

template<class Traits_>
class Width_3 {
  // +----------------------------------------------------------------------+
  // | Typedef Area                                                         |
  // +----------------------------------------------------------------------+
 private:
  typedef Traits_                   Traits;
  typedef typename Traits::Point_3  Point_3;
  typedef typename Traits::Vector_3 Vector_3;
  typedef typename Traits::Plane_3  Plane_3;
  typedef typename Traits::RT       RT;

  // +----------------------------------------------------------------------+
  // | Variable Declaration                                                 |
  // +----------------------------------------------------------------------+
 private:
  //the current best plane coefficients: e1:Ax+By+Cz+1=0
  //                                     e2:Ax+By+Cz+D=0     
  RT A,B,C,D,K;

  // Width itself
  RT WNum,WDenom;

  // Planes and directions are derived from these variables 

  // A list with all quadruples A/K, B/K, C/K, D/K
  std::vector< std::vector<RT> > allsolutions;

  // A list with all quadruples that give an optimal solution
  std::vector< std::vector<RT> > alloptimal;

  //The traits class object
  Traits tco;
  
  //The new origin to know how to translate back
  Point_3 neworigin;


  // +----------------------------------------------------------------------+
  // | Access to Private Variables                                          |
  // +----------------------------------------------------------------------+
 public:
  void get_width_coefficients(RT& a, RT& b, RT& c, RT& d, RT& k) {
    d=-A*neworigin.hx()-B*neworigin.hy()-C*neworigin.hz()+D*neworigin.hw();
    k=-A*neworigin.hx()-B*neworigin.hy()-C*neworigin.hz()+K*neworigin.hw();
    a=A*neworigin.hw();
    b=B*neworigin.hw();
    c=C*neworigin.hw();
#ifdef GCD_COMPUTATION
    simplify_solution(a,b,c,d,k);
#endif
  }

  void get_squared_width(RT& num, RT& denom) {
    num=WNum;
    denom=WDenom;
  }
  
  Vector_3 get_build_direction() {
    return tco.make_vector(A,B,C);
  }
  
  void get_width_planes(Plane_3& e1, Plane_3& e2) {
    RT a,b,c,d,k;
    get_width_coefficients(a,b,c,d,k);
    e1=tco.make_plane(a,b,c,d);
    e2=tco.make_plane(a,b,c,k);
  }
    
  void get_all_build_directions(std::vector<Vector_3>& alldir) {
    typename std::vector< std::vector<RT> >::iterator it=alloptimal.begin();
    RT a,b,c;
    while(it!=alloptimal.end()) {
      a=(*it)[0];
      b=(*it)[1];
      c=(*it)[2];
#ifdef GCD_COMPUTATION
      RT dummy1=0;
      RT dummy2=0;
      simplify_solution(a,b,c,dummy1,dummy2);
#endif
      Vector_3 dir=tco.make_vector(a,b,c);
      alldir.push_back(dir);
      ++it;
    }
  }

  int get_number_of_optimal_solutions() {
    return int(alloptimal.size());
  }

  int get_number_of_possible_solutions() {
    return int(allsolutions.size());
  }

  void get_all_possible_solutions(std::vector< std::vector<RT> >& allsol) {
    allsol.clear();
    typename std::vector< std::vector<RT> >::iterator it=allsolutions.begin();
    while(it!=allsolutions.end()) {
      RT d=-((*it)[0])*neworigin.hx()-((*it)[1])*neworigin.hy()
	-((*it)[2])*neworigin.hz()+((*it)[3])*neworigin.hw();
      RT k=-((*it)[0])*neworigin.hx()-((*it)[1])*neworigin.hy()
	-((*it)[2])*neworigin.hz()+((*it)[4])*neworigin.hw();
      RT a=((*it)[0])*neworigin.hw();
      RT b=((*it)[1])*neworigin.hw();
      RT c=((*it)[2])*neworigin.hw();
#ifdef GCD_COMPUTATION
      simplify_solution(a,b,c,d,k);
#endif
      std::vector<RT> sol;
      sol.push_back(a);
      sol.push_back(b);
      sol.push_back(c);
      sol.push_back(d);
      sol.push_back(k);
      allsol.push_back(sol);
      ++it;
    }
  }


  // +----------------------------------------------------------------------+
  // | The Con- and Destructors                                             |
  // +----------------------------------------------------------------------+
 public:
  Width_3(): A(0), B(0), C(0), D(2), K(1), WNum(0), WDenom(1) {}

  template<class InputIterator>
  Width_3( InputIterator begin, InputIterator beyond):
    A(0), B(0), C(0), D(2), K(1), WNum(0), WDenom(1) {
    INFOMSG(INFO,"Waiting for new HDS to build class Width_Polyhedron!"
	    <<std::endl<<"Working with extern additional data structures.");
    typedef typename Traits::ChullTraits CHT;
    typedef Width_polyhedron_items_3                        Items;
    typedef Polyhedron_3< Traits, Items, HalfedgeDS_list>   LocalPolyhedron;
    LocalPolyhedron P;
    convex_hull_3( begin, beyond, P, CHT());
    width_3_convex(P);
  }
  
  template <class InputPolyhedron>
  Width_3(InputPolyhedron& Poly):
    A(0), B(0), C(0), D(2), K(1), WNum(0), WDenom(1) {
    // Compute convex hull with new width_polyhedron structure
    INFOMSG(INFO,"Working with extern additional data structures.");
    width_3_convex(Poly);
  }

  ~Width_3() {
    allsolutions.clear();
    alloptimal.clear();
  }
  
  // +----------------------------------------------------------------------+
  // | Begin of private function area                                       |
  // +----------------------------------------------------------------------+
 private:
  // Just to remember:
  // E1: -axh - byh - czh - kwh <= 0          axh + byh + czh + kwh >= 0
  // E2:  axh + byh + czh + dwh <= 0
  // VF-pair: 3xE1 + 1xE2 
  // EE-pair: 2xE1 + 2xE2
  // plane equation in facets: Ax + By + Cz + 1 = 0
  //                           ax + by + cz + k = 0 (A=a/k,...)

  //-----------------------------
  //---Combinatorial functions---
  //-----------------------------
  
  // *** PREPARATION_CHECK ***
  //---------------------------
  //This function determines the next facet if the halfedge e is a
  //valable halfedge over which we can rotate. If so fnext is returned.
  //PRECONDITION: e is the LAST edge in the go_on or impassable list!
  template <class InputDA, class Halfedge_handle_, class Facet_handle_>
    bool preparation_check(InputDA& dao,
			   Halfedge_handle_& e, 
			   Facet_handle_& fnext,
			   std::vector<Halfedge_handle_>& go_on,
			   std::vector<Halfedge_handle_>& imp)
    {
    //If the halfedge flag impassable is set then we can pop e from the stack 
    //of the possibale go_on edges
    DEBUGMSG(PREPARATION_CHECK,"\nBegin PREPARATION_CHECK");
    DEBUGENDL(PREPARATION_CHECK,"Edge e: "<<e->opposite()->vertex()->point()
	      <<" --> ",e->vertex()->point());
    CGAL_precondition(go_on.back()==e);
    DEBUGMSG(ASSERTION_OUTPUT,"e is last element on stack go_on. "
	     <<"ASSERTION OK.");
    if ( dao.is_impassable(e) ) {
      DEBUGMSG(PREPARATION_CHECK," is impassable. Erase from go_on.");
      go_on.pop_back();
      DEBUGMSG(PREPARATION_CHECK,"End PREPARATION_CHECK");
      return false;
    } else {
      //If the opposite halfedge of e is already visited, then we insert
      //e in the impassable list an pop e from the stack of the go_on edges
      typename InputDA::Halfedge_handle h=e->opposite();
      if(dao.is_visited(h)) {
	DEBUGMSG(PREPARATION_CHECK," has a visited opposite edge. Set "
		 <<"impassable flag, push on impassable stack and erase"
		 <<" from go_on");
	imp.push_back(e);
	dao.set_impassable_flag(h,true);
	go_on.pop_back();
	DEBUGMSG(PREPARATION_CHECK,"End PREPARATION_CHECK");
	return false;
      } else {
	DEBUGMSG(PREPARATION_CHECK," is a valable candidate. Compute next "
		 <<"facet and erase from go_on. Set visited flag of all to "
		 <<"fnext incident edges."); 
	//e is a valable candidate. Thus set fnext to the next facet we visit
	fnext=h->facet();
	//Delete e from go_on and insert the edges of fnext (except opposite
	//of e) in the go_on list
	go_on.pop_back();
	typename InputDA::Halfedge_handle h0=h;
	h=h->next();
	while ( h!=h0) {
	  DEBUGENDL(PREPARATION_CHECK,"Adding edge to go_on stack: ",
		    h->opposite()->vertex()->point()<<"  -->  "
		    <<h->vertex()->point());
	  go_on.push_back(h);
	  dao.set_visited_flag(h,true);
	  h=h->next();
	}
	DEBUGMSG(PREPARATION_CHECK,"End PREPARATION_CHECK");
	return true;
      }
    }
  }

  // *** NEIGHBORS_OF ***
  //----------------------
  //To compute the neighbors of a vertex. The vertex is implicitely given
  //as the vertex the halfedge points to.
  template <class Halfedge_handle_, class Vertex_handle_>
    void neighbors_of(const Halfedge_handle_& h, 
		      std::vector<Vertex_handle_>& V) {
    DEBUGMSG(NEIGHBORS_OF,"\nBegin NEIGHBORS_OF"); 
    DEBUGENDL(NEIGHBORS_OF,"Determining the neighbors of: ",
	      h->vertex()->point());
    Halfedge_handle_ e=h;
    Halfedge_handle_ e0=e->opposite();
    V.clear();
    V.push_back(e0->vertex());
    e=e->next();
    //Now go around the vertex and store the neighbor vertices
    while ( e!=e0 ) {
      V.push_back(e->vertex());
      e=e->opposite()->next();
    }
#if NEIGHBORS_OF
    typename std::vector<Vertex_handle_>::iterator vtxit=V.begin();
    while(vtxit!=V.end()) {
      DEBUGENDL(NEIGHBORS_OF,"Neighbor: ",(*vtxit)->point());
      ++vtxit;
    }
#endif
    DEBUGMSG(NEIGHBORS_OF,"End NEIGHBORS_OF");

  }

  //During the algorithm we have to build union and minus set 
  //of two sets and check wheater two sets are cutting each othe ror not

  // *** SETMINUS ***
  //------------------
  //Builds the set A\B where the set A is changed 
  template <class Vertex_handle_>
    void setminus(std::vector<Vertex_handle_>& res, 
		  const std::vector<Vertex_handle_>& 
		  without) {
    DEBUGMSG(SETMINUS,"\nBegin SETMINUS");
    typename std::vector<Vertex_handle_>::iterator resit;
    typename std::vector<Vertex_handle_>::const_iterator 
      withoutit=without.begin();
    //Scan through all elements of without and check if they are also in res.
    //If so delete the element from res.
    while(withoutit!=without.end()) {
      resit=std::find(res.begin(),res.end(),*withoutit);
      if ( resit!=res.end() ) {
	CGAL_assertion((*resit)->point()==(*withoutit)->point());
	DEBUGMSG(ASSERTION_OUTPUT,"Found an element to erase. ASSERTION OK.");
	DEBUGENDL(SETMINUS,"Erase point: ",(*resit)->point());
	res.erase(resit);
      }
      ++withoutit;
    }
    DEBUGMSG(SETMINUS,"End SETMINUS");
  }

  // *** SETUNION ***
  //------------------
  //Builds the union of two sets A and B. The result is stored in A
  //POSTCONDITION: Every element in A is stored once.
  template <class Vertex_handle_> 
    void setunion(std::vector<Vertex_handle_>&res, 
		  std::vector<Vertex_handle_>& uni) {
    DEBUGMSG(SETUNION,"\nBegin SETUNION");
    typename std::vector<Vertex_handle_>::iterator 
      uniit=uni.begin();
    typename std::vector<Vertex_handle_>::iterator resit;
    //Scan the uni set and add every new element in res
    while(uniit!=uni.end()) {
      resit=std::find(res.begin(),res.end(),*uniit);
      if ( resit==res.end() ) {
	DEBUGENDL(SETUNION,"Insert new point: ",(*uniit)->point());
	res.push_back(*uniit);
      }
      ++uniit;
    }
    DEBUGMSG(SETUNION,"End SETUNION");
  }

  // *** SETCUT ***
  //----------------
  //Checks if two sets are cutting each other or not (the common elements are 
  //not determined
  template <class Vertex_handle_>
    bool setcut(std::vector<Vertex_handle_>& AA, 
		std::vector<Vertex_handle_>& BB) {
    DEBUGMSG(SETCUT,"\nBegin SETCUT");
    typename std::vector<Vertex_handle_>::iterator 
      Ait=AA.begin();
    typename std::vector<Vertex_handle_>::iterator Bfindit;
    while(Ait!=AA.end()) {
      Bfindit=std::find(BB.begin(),BB.end(),*Ait);
      if (Bfindit!=BB.end()) {
	DEBUGMSG(SETCUT,"The sets are cutting each other. Return true.");
	DEBUGMSG(SETCUT,"End SETCUT");
	return true;
      }
      ++Ait;
    }
    DEBUGMSG(SETCUT,"No common element detected. Return false");
    DEBUGMSG(SETCUT,"End SETCUT");
    return false;
  }


  // ---Numerical functions---
  // *************************

  // *** COMPUTE_PLANE_EQUATION ***
  //--------------------------------
  //We don't take the standard plane equation computation from CGAL,
  //because in the context of the width the coefficients have to
  //satisfy a system of linear inequations.
  //PRECONDITION: (0,0,0) is strictly inside the convex hull
  //POSTCONDITION:(0,0,0) is on the positive side of the plane <==> the normal
  //              vector of the plane points to the side where (0,0,0) lies
  template<class InputDA, class Facet_handle_>
    void compute_plane_equation(InputDA,
				const Facet_handle_& f) {
    DEBUGMSG(COMPUTE_PLANE_EQUATION,"\nBegin COMPUTE_PLANE_EQUATION");
    DEBUGENDL(COMPUTE_PLANE_EQUATION,"Compute plane equations of facet f: ("
	      <<f->halfedge()->opposite()->vertex()->point()<<"), (",
	      f->halfedge()->vertex()->point()<<"), ("
	      <<f->halfedge()->next()->vertex()->point()<<")");
    typename InputDA::Halfedge_handle e = f->halfedge();
    typename InputDA::PolyPoint p,q,r;
    q = e -> opposite() -> vertex() -> point();
    p = e -> vertex() -> point();
    r = e -> next() -> vertex() -> point();
    CGAL_assertion(r!=p && r!=q && p!=q);
    DEBUGMSG(ASSERTION_OUTPUT,"There are 3 different points. ASSERTION OK.");
    RT a,b,c,k;
    solve_3x3(InputDA(),p,q,r,a,b,c,k);
    f->plane()=tco.make_plane(a,b,c,k);
    DEBUGENDL(COMPUTE_PLANE_EQUATION,"Plane Coefficients: ",f->plane());
    DEBUGMSG(COMPUTE_PLANE_EQUATION,"End COMPUTE_PLANE_EQUATION");
  }
  
  // *** SOLVE_3X3 ***
  //-------------------
  //To solve a special 3x3 system. The rows of the coefficient matrix
  //are the (homogeneous) x,y,z-coordinates of points and the right
  //hand side is the homogeneous part of the point times the provided 
  //coefficient. The system is solved with Cramer's Rule. The sign of
  //the coefficients is chosen in a way that (0,0,0) lies on the positive 
  //side of the plane.
  template<class InputDA, class PolyPoint_>
    void solve_3x3(InputDA,
		   const PolyPoint_& p, 
		   const PolyPoint_& q,
		   const PolyPoint_& r, 
		   RT& a, RT& b, RT& c, RT& k) {
    DEBUGMSG(SOLVE_3X3,"\nBegin SOLVE_3X3");
    RT px,py,pz,ph;
    tco.get_point_coordinates(p,px,py,pz,ph);
    RT qx,qy,qz,qh;
    tco.get_point_coordinates(q,qx,qy,qz,qh);
    RT rx,ry,rz,rh;
    tco.get_point_coordinates(r,rx,ry,rz,rh);
    CGAL_assertion(ph>0 && qh>0 && rh>0);
    DEBUGMSG(ASSERTION_OUTPUT,"All homogeneous parts >0. ASSERTION OK.");
    DEBUGMSG(SOLVE_3X3,"Matrix:");
    DEBUGENDL(SOLVE_3X3,"",px<<" "<<py<<" "<<pz<<" : "<<-ph);
    DEBUGENDL(SOLVE_3X3,"",qx<<" "<<qy<<" "<<qz<<" : "<<-qh);
    DEBUGENDL(SOLVE_3X3,"",rx<<" "<<ry<<" "<<rz<<" : "<<-rh);
    k=px*(qy*rz-ry*qz)-qx*(py*rz-ry*pz)+rx*(py*qz-qy*pz);
    RT sig(1);
    if (k<=0) {
      if(k<0) {
	sig=-1;
	k=-k;
      }
      else 
	CGAL_assertion_msg(k!=0,"Couldn't solve plane equation system");
    }
    a=sig*(-ph*(qy*rz-ry*qz)+qh*(py*rz-ry*pz)-rh*(py*qz-qy*pz));
    b=sig*(px*(rh*qz-qh*rz)-qx*(rh*pz-ph*rz)+rx*(qh*pz-ph*qz));
    c=sig*(px*(ry*qh-qy*rh)-qx*(ry*ph-py*rh)+rx*(qy*ph-py*qh));
#ifdef GCD_COMPUTATION
    RT dummy=0;
    DEBUGENDL(SOLVE_3X3,"Solution of 3x3 (before GCD computation):\n",a
	      <<std::endl
	      <<b<<std::endl<<c<<std::endl<<k<<std::endl); 
    simplify_solution(a,b,c,k,dummy);
#endif
    DEBUGENDL(SOLVE_3X3,"Solution of 3x3:\n",a<<std::endl
	      <<b<<std::endl<<c<<std::endl<<k<<std::endl); 
    DEBUGMSG(SOLVE_3X3,"End SOLVE_3X3");
  }
  
  // *** SOLVE_4X4 ***
  //-------------------
  //To enumerate EE-pairs we need to solve a 4x4 linear equation system
  //The rows of the coefficient matrix
  //are the (homogeneous) x,y,z-coordinates of points and the right
  //hand side is the homogeneous part of the point times the provided 
  //coefficient. The system is solved with Cramer's Rule.
  template<class InputDA, class PolyPoint_>
    bool solve_4x4(InputDA,
		   const PolyPoint_& p, 
		   const PolyPoint_& q,
		   const PolyPoint_& r,
		   const PolyPoint_& v,
		   RT& a, RT& b, RT& c, RT& d, RT& k) {
    DEBUGMSG(SOLVE_4X4,"\nBegin SOLVE_4X4");
    RT px,py,pz,ph;
    tco.get_point_coordinates(p,px,py,pz,ph);
    RT qx,qy,qz,qh;
    tco.get_point_coordinates(q,qx,qy,qz,qh);
    RT rx,ry,rz,rh;
    tco.get_point_coordinates(r,rx,ry,rz,rh);
    RT vx,vy,vz,vh;
    tco.get_point_coordinates(v,vx,vy,vz,vh);
    CGAL_assertion(ph>0 && qh>0 && vh>0 && rh>0);
    DEBUGMSG(ASSERTION_OUTPUT,"All homogeneous parts >0. ASSERTION OK.");
    DEBUGMSG(SOLVE_4X4,"Matrix: ");
    DEBUGENDL(SOLVE_4X4,"",px<<" "<<py<<" "<<pz<<" 0 : "<<-ph);
    DEBUGENDL(SOLVE_4X4,"",qx<<" "<<qy<<" "<<qz<<" 0 : "<<-qh);
    DEBUGENDL(SOLVE_4X4,"",rx<<" "<<ry<<" "<<rz<<" "<<rh<<" : 0");
    DEBUGENDL(SOLVE_4X4,"",vx<<" "<<vy<<" "<<vz<<" "<<vh<<" : 0");
    
    k=-rh*(px*(qy*vz-vy*qz)-qx*(py*vz-vy*pz)+vx*(py*qz-qy*pz))
      +vh*(px*(qy*rz-ry*qz)-qx*(py*rz-ry*pz)+rx*(py*qz-qy*pz));
    RT sig(1);
    if (k<=0) {
      if (k<0) {
	sig=-1;
	k=-k;
	DEBUGMSG(SOLVE_4X4,"Sign of k (and of all other coefficients) "
		 <<"changed.");
      } else {
	DEBUGMSG(SOLVE_4X4,"No proper solution.");
	return false;
      }	
    }
    
    a=sig*(-ph*(qy*(rz*vh-vz*rh)-ry*qz*vh+vy*qz*rh)
      +qh*(py*(rz*vh-vz*rh)-ry*pz*vh+vy*pz*rh));
    b=sig*(ph*(qx*(rz*vh-vz*rh)-rx*qz*vh+vx*qz*rh)
      -qh*(px*(rz*vh-vz*rh)-rx*pz*vh+vx*pz*rh));
    c=sig*(-ph*(qx*(ry*vh-vy*rh)-rx*qy*vh+vx*qy*rh)
      +qh*(px*(ry*vh-vy*rh)-rx*py*vh+vx*py*rh));
    d=sig*(ph*(qx*(ry*vz-vy*rz)-rx*(qy*vz-vy*qz)+vx*(qy*rz-ry*qz))
      -qh*(px*(ry*vz-vy*rz)-rx*(py*vz-vy*pz)+vx*(py*rz-ry*pz)));
    if (d>k) {
      DEBUGMSG(SOLVE_4X4,"d>k: Interchange d and k");
      RT tmp=d;
      d=k;
      k=tmp;
      CGAL_assertion(a*px+b*py+c*pz+d*ph==0);
      CGAL_assertion(a*qx+b*qy+c*qz+d*qh==0);
      CGAL_assertion(a*rx+b*ry+c*rz+k*rh==0);
      CGAL_assertion(a*vx+b*vy+c*vz+k*vh==0);
      DEBUGMSG(ASSERTION_OUTPUT,"Interchanged k and d. All Assertions ok.");
    }

    if (a==0 && b==0 && c==0) {
      DEBUGENDL(SOLVE_4X4,"Solution of 4x4:\n ",a<<std::endl<<b<<std::endl<<c
		<<std::endl
		<<d<<std::endl<<k); 
      CGAL_assertion(a!=0);
      CGAL_error();
    } else {
#ifdef GCD_COMPUTATION
      DEBUGENDL(SOLVE_4X4,"Unique Solution of 4x4 (before GCD computation):\n",
		a<<std::endl<<b<<std::endl<<c<<std::endl<<d<<std::endl<<k); 
      simplify_solution(a,b,c,d,k);
#endif
      DEBUGENDL(SOLVE_4X4,"Unique Solution of 4x4:\n",
		a<<std::endl<<b<<std::endl<<c<<std::endl<<d<<std::endl<<k);
      DEBUGMSG(SOLVE_4X4,"End SOLVE_4X4");
    }
    return true;
  }

  // *** CHECK_FEASIBILITY ***
  //---------------------------
  //This function checks the feasibility of a provided quadruple A/K, B/K,
  //C/K and D/K. Because we do not want to check the feasibility for all 
  //the points the list of points is also expected.
  template<class InputDA, class Vertex_handle_>
    bool check_feasibility(InputDA,
			   const RT& a, const RT& b, const RT& c, 
			   const RT& d, const RT& k,
			   const std::vector<Vertex_handle_>& V) {
    DEBUGMSG(CHECK_FEASIBILITY,"\nBegin CHECK_FEASIBILITY");
    if (d==k) {
      DEBUGMSG(CHECK_FEASIBILITY,"The planes e1 and e2 are the same. "
	       <<"Not a feasible solution.");
      DEBUGMSG(CHECK_FEASIBILITY,"End CHECK_FEASIBILITY");
      return false;
    }

    typename std::vector<Vertex_handle_>::const_iterator 
      it=V.begin();
    RT tmp;
    while ( it!=V.end() ) {
      RT px,py,pz,ph;
      tco.get_point_coordinates((*it)->point(),px,py,pz,ph);
      tmp = a*px+b*py+c*pz;
      //Check if the restrictions according to p are satisfied
      if (tmp+k*ph < 0 || tmp + d*ph > 0) {
#if CHECK_FEASIBILITY
	DEBUGENDL(CHECK_FEASIBILITY,"Restriction to point ",
		  (*it)->point()<<" failed.");
	if (tmp+k*ph < 0){
	  DEBUGENDL(CHECK_FEASIBILITY,"E1 not satisfied: ",tmp+k*ph);
	} else { 
	  DEBUGENDL(CHECK_FEASIBILITY,"E2 not satisfied: ",tmp+d*ph);
	}
#endif
	DEBUGMSG(CHECK_FEASIBILITY,"Feasibility Check failed.");
	DEBUGMSG(CHECK_FEASIBILITY,"End CHECK_FEASIBILITY");
	return false;
      }
      ++it;
    }

    //All restrictions are satisfied, thus the check returns true
    DEBUGMSG(CHECK_FEASIBILITY,"Feasibility Check was successful.");
    DEBUGMSG(CHECK_FEASIBILITY,"End CHECK_FEASIBILITY");
    return true;
  }

#if GCD_COMPUTATION
  // *** GCD ***
  //-------------
  //To compute the gcd of 2 integer numbers
  //PRECONDITION: abs(IntNum) must be defined!
  //              %-operator must be defined!
  template<class IntNum>
    IntNum gcd(const IntNum& a, const IntNum& b) {
    DEBUGMSG(GCD_OUTPUT,"\nBegin GCD");
    DEBUGENDL(GCD_OUTPUT,"Compute gcd of ",a<<" and "<<b);
    IntNum r,s,t;
    if (abs(a)<abs(b)) {
      r=abs(b);
      s=abs(a);
    } else {
      r=abs(a);
      s=abs(b);
    }
    if (s==0) {
      DEBUGMSG(GCD_OUTPUT,"End GCD");
      return r;
    }
    t=r%s;
    while(t!=0){
      r=s;
      s=t;
      DEBUGENDL(GCD_OUTPUT,"New r: ",r<<"  and  new s: "<<s);
      t=r%s;
    }
    DEBUGENDL(GCD_OUTPUT,"Return gcd: ",s);
    DEBUGMSG(GCD_OUTPUT,"End GCD");
    return s;
  }

  // *** SIMPLIFY_SOLUTION ***
  //---------------------------
  //To simplify the solutions
  template<class IntNum>
    void simplify_solution(IntNum& a, IntNum& b, IntNum& c, IntNum& d,
			   IntNum& k) {
    DEBUGMSG(SIMPLIFY_SOLUTION,"\nBegin SIMPLIFY_SOLUTION");
    IntNum r=gcd(a,b);
    IntNum s=gcd(c,d);
    IntNum t=gcd(r,s);
    IntNum g=gcd(t,k);
    CGAL_assertion(g*(a/g)==a);
    a=a/g;
    CGAL_assertion(g*(b/g)==b);
    b=b/g;
    CGAL_assertion(g*(c/g)==c);
    c=c/g;
    CGAL_assertion(g*(d/g)==d);
    d=d/g;
    CGAL_assertion(g*(k/g)==k);
    k=k/g;
    DEBUGENDL(SIMPLIFY_SOLUTION,"Simplified solutions: ",a<<" "<<b<<" "<<c
	      <<" "<<d<<" "<<k);
    DEBUGMSG(SIMPLIFY_SOLUTION,"End SIMPLIFY_SOLUTION");
  }
#endif

  // ---Width functions---
  // *********************
  
  // *** INITIAL_VF_PAIR ***
  //-------------------------
  //After the first initialization phase we have to compute an initial 
  //Vertex-Facet pair to start with the enumeration 
  //PRECONDITION: Normal of initial plane points to the interior of the
  //              convex hull
  template<class InputDA, class Facet_handle_, class Polyhedron_,
           class Halfedge_handle_>
    void initial_VF_pair(InputDA& dao,
			 Facet_handle_& f, 
			 Polyhedron_& P,
			 std::vector<Halfedge_handle_>& go_on)
    {
    DEBUGMSG(INITIAL_VF_PAIR,"\nBegin INITIAL_VF_PAIR");
    DEBUGENDL(INITIAL_VF_PAIR,"Compute initial VF-pair with facet f: ("
	      <<f->halfedge()->opposite()->vertex()->point()<<"), (",
	      f->halfedge()->vertex()->point()<<"), ("
	      <<f->halfedge()->next()->vertex()->point()<<")");
    typedef typename InputDA::Vertex_handle Vertex_handle;
    //Compute the facet. ==> e2 is fixed
    tco.get_plane_coefficients(f->plane(),A,B,C,K);
    CGAL_assertion(K>0);
    DEBUGMSG(ASSERTION_OUTPUT,"K greater (strictly) than 0. ASSERTION OK.");

    //Start with an impossible configuration for the still unknown 
    //coefficient D=K, ie plane E1 == plane E2
    D=K;
    DEBUGENDL(INITIAL_VF_PAIR,"Starting with values:\nA:",A<<std::endl
	      <<"B: "<<B<<std::endl<<"C: "<<C<<std::endl<<"D: "
	      <<D<<std::endl<<"K: "<<K);
    
    std::vector<Vertex_handle> apv;
#if !(defined(CGAL_KERNEL_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_KERNEL_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE))\
  || defined(NDEBUG))
      typename InputDA::Vertex_iterator vtxitass = P.vertices_begin();
      while(vtxitass!=P.vertices_end()) {
	RT px,py,pz,ph;
	tco.get_point_coordinates((*vtxitass).point(),px,py,pz,ph);
	CGAL_expensive_assertion(ph>0);
	CGAL_expensive_assertion(A*px+B*py+C*pz+K*ph>=0);
	++vtxitass;
      }
      DEBUGMSG(ASSERTION_OUTPUT,"All points satisfy restriction "
	       <<"type E1. ASSERTION OK>");
#endif

    typename InputDA::Vertex_iterator vtxit=P.vertices_begin();
    RT maxdist=0;
    RT hompart=1;
    //Try every point to be an/the antipodal vertex of the facet f. Take the
    //one with the bigest distance from E1
    DEBUGENDL(INITIAL_VF_PAIR,"Plane E1:",f->plane());
    while (vtxit != P.vertices_end() ) {
      RT pix, piy, piz, pih;
      tco.get_point_coordinates((*vtxit).point(),pix,piy,piz,pih);
      DEBUGENDL(INITIAL_VF_PAIR,"Try Point: ",(*vtxit).point());

      //Compute the sign of the distance from pi to the current plane e2
      RT distpie1=A*pix + B*piy + C*piz;
      DEBUGENDL(INITIAL_VF_PAIR,"Distance from p to current plane e1: ",
		distpie1*hompart);
      //If pi is not between e1 and e2, compute a new plane e2 through pi
      //If pi is also ON the current plane e2, then insert pi in the list 
      //of current antipodal vertices of the facet f 
      if (hompart*distpie1 >= pih*maxdist) {
	DEBUGMSG(INITIAL_VF_PAIR,"Distance of this point is greater (or equal)"
		 <<" than all the distances before."
		 <<"Change plane antipodal vertices.");
	if (hompart*distpie1 > pih*maxdist) {
	  DEBUGMSG(INITIAL_VF_PAIR,"Compute new plane e2!");
	  apv.clear();
	  hompart=pih;
	  maxdist=distpie1;
	}
	apv.push_back(vtxit);
      }
      ++vtxit;
    }
    A=A*hompart;
    B=B*hompart;
    C=C*hompart;
    D=-maxdist;
    K=K*hompart;
#ifdef GCD_COMPUTATION
    simplify_solution(A,B,C,D,K);
#endif
    DEBUGENDL(INITIAL_VF_PAIR,"Initial Plane E1: ",A<<" "<<B<<" "<<C<<" "<<K);
    DEBUGENDL(INITIAL_VF_PAIR,"Initial Plane E2: ",A<<" "<<B<<" "<<C<<" "<<D);

#if !(defined(CGAL_KERNEL_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_KERNEL_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE))\
  || defined(NDEBUG))
      CGAL_expensive_assertion(D!=K && D!=0);
      DEBUGMSG(ASSERTION_OUTPUT,"A real plane E2 has been computed. "
	       <<"ASSERTION OK.");
      vtxit=P.vertices_begin();
      while (vtxit != P.vertices_end() ) {
	RT px, py, pz, ph;
	tco.get_point_coordinates((*vtxit).point(),px,py,pz,ph);
	CGAL_expensive_assertion(A*px+B*py+C*pz+K*ph>=0);
	CGAL_expensive_assertion(A*px+B*py+C*pz+D*ph<=0);
	DEBUGENDL(ASSERTION_OUTPUT,"Restriction values: E1:",
		  A*px+B*py+C*pz+K*ph<<"  and E2: "<<A*px+B*py+C*pz+D*ph);
	DEBUGENDL(ASSERTION_OUTPUT,"Restrictions E1 and E2 according to "
		  <<"point "<<(*vtxit).point()
		  <<" are both satisfied.","ASSERTION OK.");
	++vtxit;
      }
      DEBUGMSG(ASSERTION_OUTPUT,"All restrictions satisfied. "
	       <<"ASSERTION OK.");
#endif
    DEBUGENDL(INITIAL_VF_PAIR,"Initial plane E1:",A<<" "<<B<<" "<<C<<" "<<K);
    DEBUGENDL(INITIAL_VF_PAIR,"Initial plane E2:",A<<" "<<B<<" "<<C<<" "<<D);
    //set the list of antipodal vertices of f definitly
    dao.set_antipodal_vertices(f,apv);

    //All solutions
    std::vector <RT> sol;
    sol.push_back(A);
    sol.push_back(B);
    sol.push_back(C);
    sol.push_back(D);
    sol.push_back(K);
    allsolutions.push_back(sol);
    alloptimal.push_back(sol);

    //Compute the squared width with the determined coefficients
    WNum=(K-D)*(K-D);
    WDenom=A*A+B*B+C*C;
    DEBUGENDL(INITIAL_VF_PAIR,"Initial squared width: ",
	      WNum<<"/"<<WDenom);

    //Set all halfedges of f to be possible edges for a rotation
    //The set of these edges is used in the third phase of the algorithm
    typename InputDA::Halfedge_handle e = f->halfedge();
    go_on.push_back(e);
    dao.set_visited_flag(e,true);
    typename InputDA::Halfedge_handle e0 = e;
    e = e->next();
    while ( e != e0 ) {
      go_on.push_back(e);
      dao.set_visited_flag(e,true);
      e=e->next();
    }
    DEBUGMSG(INITIAL_VF_PAIR,"End INITIAL_VF_PAIR");
  }
  
  // *** CHECK_ABOUT_VF_PAIRS ***
  //------------------------------
  //This function checks if a facet and a subset of a given set of vertices
  //build a vertex facet pair
  template<class InputDA, class Facet_handle_, class Vertex_handle_>
    bool 
    check_about_VF_pairs(InputDA& dao,
			 Facet_handle_& f, 
			 const std::vector<Vertex_handle_>& V)
    {
    DEBUGMSG(CHECK_ABOUT_VF_PAIRS,"\nBegin CHECK_ABOUT_VF_PAIRS");
    DEBUGMSG(CHECK_ABOUT_VF_PAIRS,"Check, if f has antipodal vertices in "
	     <<"a set.");
    RT a,b,c,d,k;
    typename std::vector<typename InputDA::Vertex_handle>
      ::const_iterator  vtxit=V.begin();
    std::vector<typename InputDA::Vertex_handle> W;
    typename std::vector<typename InputDA::Vertex_handle>
      ::iterator neighborit;
    bool feasible=false;
    std::vector<typename InputDA::Vertex_handle> apv;
    std::vector<typename InputDA::Vertex_handle> visited_points;
    while (vtxit!=V.end()) {
      RT vx,vy,vz,vh;
      tco.get_point_coordinates((*vtxit)->point(),vx,vy,vz,vh);
      tco.get_plane_coefficients(f->plane(),a,b,c,k);
      //assume plane e2 parallel to e1 through v: e2:axh+byh+czh+d=0
      d = -a*vx - b*vy - c*vz;
      CGAL_assertion(vh > 0);
      DEBUGMSG(ASSERTION_OUTPUT,"vh is greater than zero (strictly). "
	       <<"ASSERTION OK.");
      a=tco.get_a(f->plane())*vh;
      b=tco.get_b(f->plane())*vh;
      c=tco.get_c(f->plane())*vh;
      k=tco.get_d(f->plane())*vh;
      CGAL_assertion(a*vx+b*vy+c*vz+k*vh>=0);
      CGAL_assertion(a*vx+b*vy+c*vz+d*vh==0);
      DEBUGMSG(ASSERTION_OUTPUT,"Checked: Point on the right side of e1, and "
	       <<"on e2. ASSERTION OK>");

      //If v lies on plane e1 then we can continue (v is not antipodal)
      if (d == k) {
	CGAL_assertion(a*vx+b*vy+c*vz+k*vh==0);
	DEBUGENDL(CHECK_ABOUT_VF_PAIRS,"Point "<<(*vtxit)->point()
		  <<" lies on plane ",f->plane()<<". Continue.");
	++vtxit;
	continue;
      }
      CGAL_assertion(a*vx+b*vy+c*vz+k*vh>0);
      DEBUGMSG(ASSERTION_OUTPUT,"v not on e1. ASSERTION OK.");
      //Else we look if we can find a witness in the neighborhood of v that 
      //shows that v is not an antipodal vertex of f
      neighbors_of((*vtxit)->halfedge(),W);
      
      //Assume there is no witness for infeasibility
      feasible = true;
      
      //Scan all possible witnesses
      while(!W.empty()) {
	neighborit=W.begin();
	visited_points.push_back(*neighborit);
	RT nx,ny,nz,nh;
	tco.get_point_coordinates((*neighborit)->point(),nx,ny,nz,nh);
	
	//Check if n (neighbor of v) satisfies restriction type e2
	//with the presumed plane e2, ie anx+bny+cnz-Dv<=0
	CGAL_assertion(a*nx+b*ny+c*nz+k*nh>=0);
	DEBUGMSG(ASSERTION_OUTPUT,"Restrictions E1 is satisfied. "
		 <<"ASSERTION OK.");
	if ( a*nx+b*ny+c*nz+d*nh >= 0 ) {
	  //Could be a violation. Now check if v and n lie on the
	  //same plane. If so no violation, othervise we can break
	  if (a*nx+b*ny+c*nz+d*nh == 0 ) {
	    DEBUGMSG(CHECK_ABOUT_VF_PAIRS,"Additional Antipodal Vertex "
		     <<"found. Expanding "<<"set of witnesses.");
	    //v and n are both (so far) antipodal vertices ==>EF-pair
	    apv.push_back(*neighborit);
	    //There could now be more witnesses that give violating
	    //restrictions. Therefore compute the new neighbors of n
	    std::vector<typename InputDA::Vertex_handle> Wnew;
	    neighbors_of((*neighborit)->halfedge(),Wnew);
	    
	    //Erase v from the vertices to be considered as new witnesses
	    typename
	      std::vector<typename InputDA::Vertex_handle>
	      ::iterator  res= std::find(Wnew.begin(),Wnew.end(),*vtxit);
	    if ( res!=Wnew.end() ) 
	      Wnew.erase(res);
	    //Erase all the elements we already considered from the new
	    //set of witnesses
	    setminus(Wnew,visited_points);
	    //Erase the neighbor vertex itself from the set of witnesses
	    W.erase(neighborit);
	    //Compute the new whole set of witnesses, that is add the 
	    //remaining new ones to the old set of witnesses
	    setunion(W,Wnew);
	  } else {
	    DEBUGENDL(CHECK_ABOUT_VF_PAIRS,"Violation found. Not a feasible "
		     <<"solution. Violated Point:",(*neighborit)->point());
	    //there is a violation, so do a break
	    feasible = false;
	    break;
	  }
	} else {
	  DEBUGENDL(CHECK_ABOUT_VF_PAIRS,"Restriction E2 also satisfied by "
		   <<"point:",(*neighborit)->point());
	  //There is no violating restriction according to the vertex n
	  //Erase it from the set of witnesses
	  W.erase(neighborit);
	}
      } //end while(!W.empty())
      
      //Now we can determine if we have a feasible solution or not. Because
      //the feasible flag can only be set to false during the while-loop
      //we can be sure of the feasibility of our solution (no witness found)
      //if feasible is false then we have found a witness that v is not an
      //antipodal vertex. So we go on in the list of possible antipodal 
      //vertices otherwise.
      if (feasible == true) {
	DEBUGMSG(CHECK_ABOUT_VF_PAIRS,"All witnesses checked. "
		 <<"Update width and antipodal vertices. Return true");
	apv.push_back(*vtxit);
#ifdef GCD_COMPUTATION
	simplify_solution(a,b,c,d,k);
#endif
	update_width(a,b,c,d,k);
	dao.set_antipodal_vertices(f,apv);
#ifdef VF_PAIR_OUTPUT
	DEBUGENDL(VF_PAIR_OUTPUT,"Antipodal vertices of plane: ",
		  f->plane());
	typename std::vector<Vertex_handle_>::iterator cavfpit=apv.begin();
	while(cavfpit!=apv.end()) {
	  DEBUGENDL(VF_PAIR_OUTPUT,"Antipodal Vertex: ",
		    (*cavfpit)->point());
	  ++cavfpit;
	}
#endif
	DEBUGMSG(CHECK_ABOUT_VF_PAIRS,"End CHECK_ABOUT_VF_PAIRS");
	return true;
      }
      ++vtxit;
    }//end while(vtxit!=V.end())
    //If we could not return with antipodal vertices we return false
    DEBUGMSG(CHECK_ABOUT_VF_PAIRS,"No new VF-pair found. Return false.");
    DEBUGMSG(CHECK_ABOUT_VF_PAIRS,"End CHECK_ABOUT_VF_PAIRS");
    return false;
  }

  // *** UPDATE_WIDTH ***
  //----------------------
  //This function we use to update the current best width. The old width is
  //compared with a new provided one and the better solution will we taken 
  //as the new width. This function also saves all the possible quadruples 
  //to be the width of the point set.
  void update_width(RT& a, RT& b, RT& c, RT& d, RT& k) {
    //Update the list of all possible solutions
    DEBUGMSG(UPDATE_WIDTH,"\nBegin UPDATE_WIDTH");
    std::vector<RT> sol;
    sol.push_back(a);
    sol.push_back(b);
    sol.push_back(c);
    sol.push_back(d);
    sol.push_back(k);
    allsolutions.push_back(sol);
    
    //Compute the squared width provided by the new solution
    RT tocompareNum=(k-d)*(k-d);
    RT tocompareDenom=(a*a+b*b+c*c);
    DEBUGENDL(UPDATE_WIDTH,"New possible width: ",tocompareNum
	      <<" / "<<tocompareDenom);
    //Compare with old width
    if (WNum*tocompareDenom >= tocompareNum*WDenom) {
      if (WNum*tocompareDenom > tocompareNum*WDenom){
	DEBUGMSG(UPDATE_WIDTH,"Optimal width changes");
	WNum=tocompareNum;
	WDenom=tocompareDenom;
	alloptimal.clear();
	alloptimal.push_back(sol);
	A=a;
	B=b;
	C=c;
	D=d;
	K=k;
      } else {
	//now we have an additional optimal solution
	alloptimal.push_back(sol);
      }
    }//end if equal or better width
    DEBUGMSG(UPDATE_WIDTH,"End UPDATE_WIDTH");
  }

  // *** EE_COMPUTATION ***
  //------------------------
  //During the 3rd phase of the width-algorithm we have to rotate planes to
  //enumerate all possible edge-edge pairs. This rotating (in primal context)
  //resp. tracking edges (in the dual context) is made by the following 
  //function. The edge we rotate about is called e. To ensure only to 
  //enumerate a pair once (...only going forward) we need a set of 
  //already visited vertices (Visited) and a set of vertices from that we know
  //they are antipodal to the first facet (V). In this function we don't
  //know the antipodal vertices of the second facet.
  template <class InputDA, class Halfedge_handle_, class Vertex_handle_>
    void EE_computation(InputDA,
			Halfedge_handle_& e,
			std::vector<Vertex_handle_>& V, 
			std::vector<Vertex_handle_>& Visited,
			std::vector<Vertex_handle_>& Nnew) {
    DEBUGMSG(EE_COMPUTATION,"\nBegin EE_COMPUTATION");
    //Compute end points of e and two witnesses: Each in one of the two
    //facets participating 
    Point_3 p,q;
    p=e->opposite()->vertex()->point();
    q=e->vertex()->point();
    typename InputDA::Vertex_handle w1=e->next()->vertex();
    typename InputDA::Vertex_handle w2=e->opposite()->next()->vertex();

    //prepare for the rotating procedure
    Nnew.clear();
    typename std::vector<typename InputDA::Vertex_handle>
      ::iterator  vtxit=V.begin();

    //Consider all the vertices in V. EE-pairs consist of p,q and the 
    //vertex v in V and another neighbor vertex of v
    while(vtxit != V.end() ) {
      std::vector<typename InputDA::Vertex_handle> R;
      neighbors_of((*vtxit)->halfedge(),R);
      std::vector<typename InputDA::Vertex_handle> Witnesses;
      //The set of witnesses are all neighbor vertices of v (=R) and the
      //two vertices "on the other side" that ensure not rotating too far
      Witnesses.push_back(w1);
      Witnesses.push_back(w2);
      setunion(Witnesses,R);
      //The neighbor vertices of v that are also in the basic set V are of no
      //interest, so we exclude them
      setminus(R,V);
      //The set of all vertices we have already visited is also of no interest
      setminus(R,Visited);

      typename std::vector<typename InputDA::Vertex_handle>
	::iterator rit=R.begin();
      //Now look at the modified set of neighbor vertices. For each neighbor r
      //we assume (p,q) and (v,r) to be an EE-pair and want then to find
      //witnesses that against this quadruple. If no such witness exist
      //(p,q) and (v,r) are a legal EE-pair. In that case we break the 
      //quest and update all the sets Visited Nnew and V. If we have considered
      //all vertices v in V and all the respective r in R and if we have 
      //not found a legal EE-pair, then an error occurs
      while (rit!=R.end()) {
	RT a,b,c,d,k;
	//It could be that the system is not uniquely solvable we only want to 
	//enumerate proper solutions no degenerate ones
	if(solve_4x4(InputDA(),p,q,(*rit)->point(),(*vtxit)->point(),
			      a,b,c,d,k)){
	  DEBUGMSG(EE_COMPUTATION,"Now we check if the provided "
		   <<"solution is a feasible one.");
#if !(defined(CGAL_KERNEL_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_KERNEL_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE))\
  || defined(NDEBUG))
	    RT px,py,pz,ph;
	    tco.get_point_coordinates(p,px,py,pz,ph);
	    CGAL_expensive_assertion(a*px+b*py+c*pz+k*ph>=0);
	    CGAL_expensive_assertion(a*px+b*py+c*pz+d*ph<=0);
	    tco.get_point_coordinates(q,px,py,pz,ph);
	    CGAL_expensive_assertion(a*px+b*py+c*pz+k*ph>=0);
	    CGAL_expensive_assertion(a*px+b*py+c*pz+d*ph<=0);
	    tco.get_point_coordinates((*rit)->point(),px,py,pz,ph);
	    CGAL_expensive_assertion(a*px+b*py+c*pz+k*ph>=0);
	    CGAL_expensive_assertion(a*px+b*py+c*pz+d*ph<=0);
	    tco.get_point_coordinates((*vtxit)->point(),px,py,pz,ph);
	    CGAL_expensive_assertion(a*px+b*py+c*pz+k*ph>=0);
	    CGAL_expensive_assertion(a*px+b*py+c*pz+d*ph<=0);
	    DEBUGMSG(ASSERTION_OUTPUT,"All restrictions to the 4 points "
		     <<"are satisfied. ASSERTION OK.");
#endif
	  if (check_feasibility(InputDA(),a,b,c,d,k,Witnesses)) {
	    DEBUGMSG(EE_COMPUTATION,"Update Width and compute all "
		     <<"active restrictions");
	    //Therefore we update the width
	    update_width(a,b,c,d,k);
	    //The next region we consider (because we only go forward) 
	    //contains (at least) r
	    Nnew.push_back(*rit);

	    //Now we look if we are in a special case, that is we look if
	    //other restrictions according to  neighboring vertices of r 
	    //are also active. If so Nnew is expanded we them.
	    std::vector<typename InputDA::Vertex_handle> S;
	    neighbors_of((*rit)->halfedge(),S);
	    //Because we only go forward we exclude v from the neighbor set S
	    std::vector<typename InputDA::Vertex_handle> vtemp;
	    vtemp.push_back(*vtxit);
	    setminus(S,vtemp);
	    //The check of more than 4 active restrictions begins
	    typename
	      std::vector<typename InputDA::Vertex_handle>
	      ::iterator sit; 
	    while(!S.empty()) {
	      sit=S.begin();
	      RT sx,sy,sz,sh;
	      tco.get_point_coordinates((*sit)->point(),sx,sy,sz,sh);
	      if (a*sx+b*sy+c*sz+d*sh==0) {
		//This special case occurs now. Thus we extend Nnew
		Nnew.push_back(*sit);
		//In the neighborhood of this new active vertex could
		//also be other new active vertices but we are only interested
		//in new ones
		std::vector<typename InputDA::Vertex_handle> T;
		neighbors_of((*sit)->halfedge(),T);
		S.erase(sit);
		setunion(S,T);
		T.clear();
		T.push_back(*vtxit);
		setminus(S,T);
		setminus(S,Nnew);
	      } else {
		//s is not active and we can erase it
		S.erase(sit);
	      }
	    } //end while (!S.empty())
	    //Since we have now enumerated all EE-pairs with the active 
	    //restrictions according to p,q and v we can now leave
	    //Nnew contains now all new active restrictions
	    DEBUGMSG(EE_COMPUTATION,"End EE_COMPUTATION");
	    return;
	  } //end if (feasible)
	}//end if(proper)
	//Try next r
	++rit;
      }//end while(!R.empty())
      //Try new v
      ++vtxit;
    }
    //There must be a new EE-pair. If not, an error occurs
    std::cerr<<"No new EE-pair found!"<<std::endl;
    CGAL_error();
  }  
  
  
  // *** EE_PAIRS ***
  //------------------------
  //This function is similar to EE_computation. The difference is that now
  //we know the antipodal vertices of BOTH participating facets 
  template <class InputDA, class Halfedge_handle_>
    void EE_pairs(InputDA& dao,
		  Halfedge_handle_& e,
		  std::vector<Halfedge_handle_>& impassable) {
    DEBUGMSG(EE_PAIRS,"\nBegin EE_PAIRS");
    Point_3 p,q;
    p=e->opposite()->vertex()->point();
    q=e->vertex()->point();
    typename InputDA::Vertex_handle w1=e->next()->vertex();
    typename InputDA::Vertex_handle w2=e->opposite()->next()->vertex();
    typename InputDA::Facet_handle f1=e->facet();
    typename InputDA::Facet_handle f2=e->opposite()->facet();
    std::vector<typename InputDA::Vertex_handle> V1;
    std::vector<typename InputDA::Vertex_handle> V2;
    dao.get_antipodal_vertices(f1,V1);
    dao.get_antipodal_vertices(f2,V2);
    std::vector<typename InputDA::Vertex_handle> N,V,Visited;
    V=V1;
    bool do_break = false;

    while (!setcut(V,V2)) {
      do_break = false;
      typename std::vector<typename InputDA::Vertex_handle>
	::iterator  vtxit=V.begin();
      std::vector<typename InputDA::Vertex_handle> R;
      while (vtxit!=V.end()) {
	neighbors_of((*vtxit)->halfedge(),R);
	std::vector<typename InputDA::Vertex_handle> Witnesses;
	Witnesses.push_back(w1);
	Witnesses.push_back(w2);
	setunion(Witnesses,R);
	setminus(R,V);
	setminus(R,Visited);
	typename std::vector<typename InputDA::Vertex_handle>
	  ::iterator  rit=R.begin();
	while (rit!=R.end()) {
	  RT a,b,c,d,k;
	  //It could be that the system is not uniquely solvable we only want 
	  //to enumerate proper solutions no degenerate ones
	  if(solve_4x4(InputDA(),p,q,(*rit)->point(),(*vtxit)->point(),
				a,b,c,d,k)){
	    if (check_feasibility(InputDA(),a,b,c,d,k,Witnesses)) {
	      update_width(a,b,c,d,k);
	      N.push_back(*rit);
	      std::vector<typename InputDA::Vertex_handle> S;
	      neighbors_of((*rit)->halfedge(),S);
	      setminus(S,V);
	      typename 
		std::vector<typename InputDA::Vertex_handle>
		::iterator sit; 
	      while(!S.empty()) {
		sit=S.begin();
		RT sx,sy,sz,sh;
		tco.get_point_coordinates((*sit)->point(),sx,sy,sz,sh);
		if (a*sx+b*sy+c*sz+d*sh== 0) {
		  N.push_back(*sit);
		  std::vector<typename InputDA::Vertex_handle>
		    T;
		  neighbors_of((*sit)->halfedge(),T);
		  S.erase(sit);
		  setunion(S,T);
		  T.clear();
		  T.push_back(*vtxit);
		  setminus(S,T);
		  setminus(S,N);
		} else {
		  S.erase(sit);
		}
	      }//end while (!S.empty())
	      do_break=true;
	      break;
	    }//if (feasible)
	  }//if(proper)
	  ++rit;
	}//end while(!R.empty())
	if (do_break == true) 
	  break;
	++vtxit;
      }//end while(!V.end())
      setunion(Visited,V);
      V=N;
    }
    impassable.pop_back();
    //Go on with next edge
    DEBUGMSG(EE_PAIRS,"End EE_PAIRS");
  }  


  // *** ORIGIN_INSIDE_CH ***
  //-------------------------
  // To ensure that zero lies completly inside the convex hull of a point set.
  // Returns true if the point set is not coplanar, false otherwise
  // PRECONDITION: Iterator range has at least 3 points 
  template<class InputDA, class Vertex_iterator_>
    bool origin_inside_CH(Vertex_iterator_& start,
			  Vertex_iterator_& beyond, 
			  InputDA){
    DEBUGMSG(ORIGIN_INSIDE_CH,"\nBegin ORIGIN_INSIDE_CH");
    typename InputDA::Vertex_iterator first=start;
    //Take 4 points that build a tetrahedron. This tetrahedron is also 
    //contained in the convex hull of the points. Thus every point 
    //in/on this tetrahedron is a valable point for a new origin
    typename InputDA::PolyPoint p,q,r,s;
    p=(*first).point();
    ++first;
    q=(*first).point();
    ++first;
    r=(*first).point();
    ++first;
    RT px,py,pz,ph,qx,qy,qz,qh,rx,ry,rz,rh;
    tco.get_point_coordinates(p,px,py,pz,ph);
    tco.get_point_coordinates(q,qx,qy,qz,qh);
    tco.get_point_coordinates(r,rx,ry,rz,rh);
    CGAL_assertion(ph>0 && qh>0 && rh>0);
    RT tmpa,tmpb,tmpc,tmpk;
    tmpk=px*(qy*rz-ry*qz)-qx*(py*rz-ry*pz)+rx*(py*qz-qy*pz);
    tmpa=-ph*(qy*rz-ry*qz)+qh*(py*rz-ry*pz)-rh*(py*qz-qy*pz);
    tmpb=px*(rh*qz-qh*rz)-qx*(rh*pz-ph*rz)+rx*(qh*pz-ph*qz);
    tmpc=px*(ry*qh-qy*rh)-qx*(ry*ph-py*rh)+rx*(qy*ph-py*qh);
#ifdef GCD_COMPUTATION
    RT dummy=0;
    DEBUGENDL(ORIGIN_INSIDE_CH,"Solution of 3x3 (before GCD "
	      <<"computation):\n",tmpa<<std::endl
	      <<tmpb<<std::endl<<tmpc<<std::endl<<tmpk<<std::endl); 
    simplify_solution(tmpa,tmpb,tmpc,tmpk,dummy);
#endif
    if (first==beyond) {
      DEBUGMSG(ORIGIN_INSIDE_CH,"3 coplanar Points. Computed plane through "
	       <<"these points. Width=0.");
      WNum=0;
      WDenom=1;
      A=tmpa;
      B=tmpb;
      C=tmpc;
      K=tmpk;
      DEBUGENDL(ORIGIN_INSIDE_CH,"Solution of 3x3:\n",A<<std::endl
		<<B<<std::endl<<C<<std::endl<<K<<std::endl); 
      D=K;
      std::vector <RT> sol;
      sol.push_back(A);
      sol.push_back(B);
      sol.push_back(C);
      sol.push_back(D);
      sol.push_back(K);
      allsolutions.push_back(sol);
      alloptimal.push_back(sol);
      DEBUGMSG(ORIGIN_INSIDE_CH,"End ORIGIN_INSIDE_CH");
      return false;
    } else {
      s=(*first).point();
      RT sx,sy,sz,sh;
      tco.get_point_coordinates(s,sx,sy,sz,sh);
      //Ensure that the 4 points are not coplanar. If so take another 4th point
      while (tmpa*sx+tmpb*sy+tmpc*sz+tmpk*sh==0 && first!=beyond) {
	s=(*first).point();
	tco.get_point_coordinates(s,sx,sy,sz,sh);
	++first;
      }
      //If we could not find a valable 4th point, then the set of the points
      //is coplanar. Therefore the width is zero and we can terminate the 
      //algorithm
      if (tmpa*sx+tmpb*sy+tmpc*sz+tmpk*sh==0) {
	DEBUGMSG(ORIGIN_INSIDE_CH,"n coplanar Points. Compute plane through "
		 <<"these points. Width=0.");
	WNum=0;
	WDenom=1;
	A=tmpa;
	B=tmpb;
	C=tmpc;
	K=tmpk;
	DEBUGENDL(ORIGIN_INSIDE_CH,"Solution of 3x3:\n",A<<std::endl
		  <<B<<std::endl<<C<<std::endl<<K<<std::endl); 
	D=K;
	std::vector <RT> sol;
	sol.push_back(A);
	sol.push_back(B);
	sol.push_back(C);
	sol.push_back(D);
	sol.push_back(K);
	allsolutions.push_back(sol);
	alloptimal.push_back(sol);
	DEBUGMSG(ORIGIN_INSIDE_CH,"End ORIGIN_INSIDE_CH");
	return false;
      } else {
	//Take center of tetrahedron pqrs
	RT ux,uy,uz,uh,vx,vy,vz,vh,nox,noy,noz,noh;
	ux=px*qh+ph*qx;
	vx=rx*sh+rh*sx;
	uy=py*qh+ph*qy;
	vy=ry*sh+rh*sy;
	uz=pz*qh+ph*qz;
	vz=rz*sh+rh*sz;
	uh=RT(2)*ph*qh;
	vh=RT(2)*rh*sh;
	nox=ux*vh+uh*vx;
	noy=uy*vh+uh*vy;
	noz=uz*vh+uh*vz;
	noh=RT(2)*uh*vh;
	neworigin=tco.make_point(nox,noy,noz,noh);
	CGAL_assertion(noh!=0);
	DEBUGENDL(ORIGIN_INSIDE_CH,"New Origin: ",neworigin);
	//Translate all the points
	first=start;
	while(first!=beyond) {
	  typename InputDA::PolyPoint tmp=(*first).point();
	  RT tmpx,tmpy,tmpz,tmph;
	  tco.get_point_coordinates(tmp,tmpx,tmpy,tmpz,tmph);
	  RT newx,newy,newz,newh;
	  newx=tmpx*noh-tmph*nox;
	  newy=tmpy*noh-tmph*noy;
	  newz=tmpz*noh-tmph*noz;
	  newh=tmph*noh;
	  DEBUGENDL(ORIGIN_INSIDE_CH,"Old Point: ",(*first).point());
#ifdef GCD_COMPUTATION
	  RT dummy=0;
	  simplify_solution(newx,newy,newz,newh,dummy);
#endif	
	  (*first).point()=tco.make_point(newx,newy,newz,newh);
	  DEBUGENDL(ORIGIN_INSIDE_CH,"New Point: ",(*first).point());
	  ++first;
	}
	DEBUGMSG(ORIGIN_INSIDE_CH,"Zero now inside polyhedron.");
	DEBUGMSG(ORIGIN_INSIDE_CH,"End ORIGIN_INSIDE_CH");
	return true;
      }
    }
  }


  /* ****************************************************** */
  /* *** --- *** The main enumeration functions *** --- *** */
  /* ****************************************************** */
  template<class InputPolyhedron>
  void width_3_convex(InputPolyhedron &P) {
    DEBUGMSG(WIDTH_3_CONVEX,"\nBegin WIDTH_3_CONVEX");
    typedef CGAL::Width_3_internal::Data_access<InputPolyhedron,Traits> DA;
    typedef typename DA::Facet_handle Facet_handle;
    typedef typename DA::Vertex_handle Vertex_handle;
    typedef typename DA::Halfedge_handle Halfedge_handle;
    typedef typename DA::Vertex_iterator Vertex_iterator;
    //Ensure that Polyhedron has at least one vertex
    CGAL_assertion_msg(P.size_of_vertices()>2,
         "Can not compute width of a 0, 1 or 2-vertex polyhedron");

    Vertex_iterator first=P.vertices_begin();
    Vertex_iterator beyond=P.vertices_end();

    //Begin with Phase 2
    if (origin_inside_CH(first,beyond,DA())) {
      DEBUGMSG(WIDTH_3_CONVEX,"Origin is now Inside the Polyhedron. "
	       <<std::endl
	       <<"And polyhedron has at least 4 not coplanar vertices");

      DA dao;
      std::vector<Halfedge_handle> go_on;
      std::vector<Halfedge_handle> impassable;

      //Ensure that the plane equations are determined because of the 
      //compare operator in DA
      Facet_handle feq=P.facets_begin();
      while(feq!=P.facets_end()) {
	compute_plane_equation(DA(),feq);
	++feq;
      }
      DEBUGMSG(WIDTH_3_CONVEX,"All plane equations of all facets computed.");

      //ensure all flags are false
#if !(defined(CGAL_KERNEL_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
      || defined(NDEBUG))
	int halfedgecount=0;
#endif
	Halfedge_handle esf=P.halfedges_begin();
	while(esf!=P.halfedges_end()) {
#if !(defined(CGAL_KERNEL_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
      || defined(NDEBUG))
          ++halfedgecount;
#endif
	  DEBUGENDL(EDGE_INITIALIZING,"Edge e: "
		    <<esf->opposite()->vertex()->point()
		    <<" --> ",esf->vertex()->point());
	  dao.set_visited_flag(esf,false);
	  dao.set_impassable_flag(esf,false);
	  ++esf;
	}
      
#if !(defined(CGAL_KERNEL_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
      || defined(NDEBUG))
      CGAL_assertion(int(P.size_of_halfedges())==halfedgecount);
      DEBUGENDL(WIDTH_3_CONVEX,"Visited all ",halfedgecount
		<<" halfedges. ASSERTION OK.");
      CGAL_assertion(dao.size_of_visited()==halfedgecount);
      CGAL_assertion(dao.size_of_impassable()==halfedgecount);
      DEBUGMSG(WIDTH_3_CONVEX,"Map sizes of visited and impassable "
	       <<"halfedges are correct. ASSERTION OK.");
#endif
      DEBUGMSG(WIDTH_3_CONVEX,"All flags set to false.");

      //Now begin with the main enumeration 
      Facet_handle f = P.facets_begin();
      initial_VF_pair(dao,f,P,go_on);

#if !(defined(CGAL_KERNEL_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_KERNEL_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE))\
  || defined(NDEBUG))
	Vertex_iterator vtxass=P.vertices_begin();
	while(vtxass!=P.vertices_end()) {
	  RT px,py,pz,ph;
	  tco.get_point_coordinates(vtxass->point(),px,py,pz,ph);
	  CGAL_expensive_assertion(A*px+B*py+C*pz+K*ph>=0);
	  CGAL_expensive_assertion(A*px+B*py+C*pz+D*ph<=0);
	  ++vtxass;
	}
	//Assert that the initial facet has antipodal vertices
	//and that all incident edges are visited (flag=true) but
	//that the impassable flag is not set yet.
	std::vector<Vertex_handle> avass;
	dao.get_antipodal_vertices(f,avass);
	DEBUGENDL(ASSERTION_OUTPUT,"Size of avass: ",avass.size());
	CGAL_expensive_assertion(avass.size()!=0);
	Halfedge_handle eass=f->halfedge();
	Halfedge_handle eass0=eass;
	CGAL_expensive_assertion(dao.is_visited(eass));
	CGAL_expensive_assertion(!dao.is_impassable(eass));
	eass=eass->next();
	while (eass != eass0) {
	  CGAL_expensive_assertion(dao.is_visited(eass));
	  CGAL_expensive_assertion(!dao.is_impassable(eass));
	  eass=eass->next();
	}
	DEBUGMSG(ASSERTION_OUTPUT,"All edges of the first initial facet "
		 <<"has a visited flag.");
#endif
      // Begin Phase 3
      Facet_handle fnext;
      Halfedge_handle e;
      std::vector<Vertex_handle> Visited;
      std::vector<Vertex_handle> N;
      std::vector<Vertex_handle> Nnew;
      
      //While there still exist an edge we can rotate an incident facet with 
      //known antipodal vertices in the other facet with unknown antipodal 
      //vertices then do this rotation
      while ( !go_on.empty()) {
	DEBUGENDL(WIDTH_3_CONVEX,"Size of go_on: ",go_on.size());
#ifdef GO_ON_OUTPUT
	DEBUGMSG(GO_ON_OUTPUT,"Edges on stack go_on:");
	typename std::vector<Halfedge_handle>::iterator 
	  goonit=go_on.begin();
	while(goonit!=go_on.end()) {
	  DEBUGENDL(GO_ON_OUTPUT,"Edge: ",
		    (*goonit)->opposite()->vertex()->point()<<" --> "
		    <<(*goonit)->vertex()->point());
	  ++goonit;
	}
#endif
	//Take last edge on stack go_on
	e=go_on.back();
	//Check if e is a proper edge or not. If so determine fnext
	if (preparation_check(dao,e,fnext,go_on,impassable)) {
	  DEBUGMSG(WIDTH_3_CONVEX,"Preparation Check successful");
	  //f is the facet of which we know the antipodal vertices
	  f=e->facet();
	  Visited.clear();
	  dao.get_antipodal_vertices(f,N);
	  CGAL_assertion (!N.empty());
	  DEBUGMSG(ASSERTION_OUTPUT,"f has some antipodal vertices. Assertion "
		   <<"successful.");
	  while(!check_about_VF_pairs(dao,fnext,N)) {
	    DEBUGMSG(WIDTH_3_CONVEX,"No new VF-pair. Continue (Begin) "
		     <<"rotation of the planes.");
	    EE_computation(DA(),e,N,Visited,Nnew);
	    DEBUGMSG(WIDTH_3_CONVEX,"Planes have been rotated. Check now "
		     <<"for a new VF-pair");
	    setunion(Visited,N);
	    N=Nnew;
	  }
	}
      }
#if !(defined(CGAL_KERNEL_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_KERNEL_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE))\
  || defined(NDEBUG))
      Facet_handle fass=P.facets_begin();
      while(fass!=P.facets_end()) {
	std::vector<Vertex_handle> apvass;
	dao.get_antipodal_vertices(fass,apvass);
	DEBUGENDL(ASSERTION,"Current checking facet: ",fass->plane());
	CGAL_assertion(!apvass.empty());
	++fass;
      }
      DEBUGMSG(ASSERTION,"All facets have antipodal vertices. "
	       <<"ASSERTION OK.");
      Facet_handle fec = P.facets_begin();
      std::vector<Halfedge_handle> fakego_on;
      DA daoec;
      while(fec!=P.facets_end()) {
	std::vector<Vertex_handle> avec;
	std::vector<Vertex_handle> avivf;
	initial_VF_pair(daoec,fec,P,fakego_on);
	daoec.get_antipodal_vertices(fec,avec);
	dao.get_antipodal_vertices(fec,avivf);
	CGAL_assertion(int(avivf.size())==int(avec.size()));
	CGAL_assertion(int(avec.size())>0);
	DEBUGENDL(EXPENSIVE_CHECKS_OUTPUT,"Antipodal vertices of facet: ("
		  <<fec->halfedge()->opposite()->vertex()->point()
		  <<"), (",fec->halfedge()->vertex()->point()<<"), ("
		  <<fec->halfedge()->next()->vertex()->point()<<")");
	std::vector<Vertex_handle>::iterator vtxit=avec.begin();
	while(vtxit!=avec.end()) {
	  std::vector<Vertex_handle>::iterator it;
	  it=std::find(avivf.begin(),avivf.end(),*vtxit);
	  CGAL_assertion(it!=avivf.end());
	  DEBUGENDL(EXPENSIVE_CHECKS_OUTPUT,"Antipodal vertex: ",
		    (*vtxit)->point());
	  ++vtxit;
	}
	++fec;
      }
      DEBUGMSG(EXPENSIVE_CHECKS_OUTPUT,"All VF-pairs verified. "
	       <<"Expensive Check successful.");
//
//    This assertion should not currently be true since the convex hull
//    polyhedron is triangulated;  no postprocessing is done to merge coplanar
//    neighboring facets.
//
//    CGAL_assertion(dao.size_of_antipodal_vertices()
//                    ==int(P.size_of_facets()));
#endif
      //Begin with phase 4. As long as the set of impassable edges is not empty
      //rotate one of the planes into the other sharing the impassable edge
      while(!impassable.empty()) {
	//Take top edge on stack impassable
	e=impassable.back();
	EE_pairs(dao,e,impassable);
	//In EE_pairs the top element will be removed
      }
    }
    DEBUGMSG(WIDTH_3_CONVEX,"Width computed.");
  }
};  

} //namespace CGAL

#endif
