// ======================================================================
//
// Copyright (c) 1999 The GALIA Consortium
//
// This software and related documentation is part of the
// Computational Geometry Algorithms Library (CGAL).
//
// Every use of CGAL requires a license. Licenses come in three kinds:
//
// - For academic research and teaching purposes, permission to use and
//   copy the software and its documentation is hereby granted free of  
//   charge, provided that
//   (1) it is not a component of a commercial product, and
//   (2) this notice appears in all copies of the software and
//       related documentation.
// - Development licenses grant access to the source code of the library 
//   to develop programs. These programs may be sold to other parties as 
//   executable code. To obtain a development license, please contact
//   the GALIA Consortium (at cgal@cs.uu.nl).
// - Commercialization licenses grant access to the source code and the
//   right to sell development licenses. To obtain a commercialization 
//   license, please contact the GALIA Consortium (at cgal@cs.uu.nl).
//
// This software and documentation is provided "as-is" and without
// warranty of any kind. In no event shall the CGAL Consortium be
// liable for any damage of any kind.
//
// The GALIA Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Free University of Berlin (Germany),
// INRIA Sophia-Antipolis (France), Trier University
// (Germany), Max-Planck-Institute Saarbrucken (Germany),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// file          : demo/GeoWin/CGAL_min_sphere_3.C
//
// ======================================================================

#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) || (__LEDA__ < 400)
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA 4.0 or higher installed!\n";
 std::cout << "A LEDA version >= 4.0 is required to run GeoWin!\n";
 return 0;
}
#else 

#include<CGAL/Cartesian.h>
#include<CGAL/Cartesian_d.h>
#include <CGAL/geowin_support.h>
#include<iostream>
#include<cstdlib>
#include<CGAL/Optimisation_d_traits_d.h>
#include<CGAL/Min_sphere_d.h>
#include<CGAL/leda_rational.h>

#include <LEDA/geowin.h>
#include <LEDA/graph.h>
#include <LEDA/d3_sphere.h>
#include <LEDA/d3_hull.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

using namespace std;

//#if !defined(_MSC_VER)
//#define USE_RAT
//#endif

#if defined(USE_RAT)
typedef leda_rational                         NT;
typedef leda_integer                          NT2;
#else
typedef double                                NT;
typedef double                                NT2;
#endif

typedef CGAL::Cartesian<NT>                   K;
typedef CGAL::Cartesian_d<NT>                 Kd;
typedef CGAL::Optimisation_d_traits_d<Kd>     Traits;
typedef CGAL::Min_sphere_d<Traits>            Min_sphere;
typedef K::Point_3                            Point;
typedef Kd::Point_d                           Pointd;

// dimension
const int dim = 3; 

// --------------------------------------------
// routines for generation of the sphere graph
// --------------------------------------------

void circle_segs(leda_list<leda_segment>& LS, leda_circle C, int n)
{
  leda_list<leda_rat_point> L;
  leda_point p = C.point1();
  leda_point q = C.point2();
  leda_point r = C.point3();
  
  leda_rat_point rp(p), rq(q), rr(r);
  leda_rat_circle R(rp,rq,rr);
  
  double d = (2*LEDA_PI)/n;
  double eps = 0.001;
  double a = 0;
  
  for(int i=0; i < n; i++)
    { 
      leda_rat_point pp = R.point_on_circle(a,eps);
      L.append(pp);
      a += d;
    }
    
  // now generate the segments desribing the circle ...
  list_item lit = L.first();
  
  while(lit && L.succ(lit))
  {
    LS.append(leda_segment(L[lit].to_point(), L[L.succ(lit)].to_point()));
    lit = L.succ(lit);
  }
  LS.append(leda_segment(L.tail().to_point(), L.head().to_point()));  
}

leda_list<leda_segment> zoom_segments(const leda_list<leda_segment>& LS, const leda_point& pm, double zf)
{
 leda_segment siter;
 leda_list<leda_segment> rt;
 forall(siter,LS){
   leda_vector v1 = (siter.source() - pm) * zf;
   leda_vector v2 = (siter.target() - pm) * zf;
   rt.append(leda_segment(leda_point(pm+v1),leda_point(pm+v2)));
 }
 return rt;
}

void generate_sphere_graph(const leda_d3_sphere& Sph, GRAPH<leda_d3_point,int>& G)
{
  leda_list<leda_d3_point> pts;
  
  leda_d3_point ct = Sph.center();
  double r  = Sph.radius();
  double ymin = Sph.center().ycoord() - 0.99*r;
  double ymax = Sph.center().ycoord() + 0.99*r;
  double yakt = ymin, ystep = 0.15*r;
  double xmin;
  leda_point ct2(ct.xcoord(), ct.ycoord());
  leda_circle C(ct2,r);
  leda_circle C2;
  leda_d3_point iter;
  list_item lit;
  leda_list<leda_d3_point> pts1,pts2;
  leda_list<leda_node> NL,NL2,NL_first;
  leda_edge e1,e2;
  leda_list<leda_segment> Lseg;
  circle_segs(Lseg, leda_circle(ct2, 1.0), 30);
  int w;
  
  
  while(yakt <= ymax){
    leda_segment S(leda_point(ct.xcoord()-r-100.0,yakt), leda_point(ct.xcoord()+r+100.0,yakt));
    leda_list<leda_point> res = C.intersection(S);
    leda_point p1 = res.pop(), p2 = res.pop();
    if (p1.xcoord() < p2.xcoord()) { xmin = p1.xcoord(); }
    else { xmin = p2.xcoord(); }
    C2 = leda_circle(ct2,ct.xcoord()-xmin); 
    
    leda_list<leda_segment> LS = zoom_segments(Lseg,ct2,ct.xcoord()-xmin);
    leda_segment siter;
    forall(siter,LS){
      leda_point akt = siter.source();
      pts.append(leda_d3_point(akt.xcoord(),yakt,akt.ycoord()));
    }

    NL2.clear();
    forall(iter,pts) {
      leda_node nn = G.new_node(iter);
      NL2.append(nn);
    }
    
    list_item prev = NL2.first();
    lit = NL2.cyclic_succ(prev);
    
    for(;lit != NL2.first(); prev = lit, lit = NL2.cyclic_succ(lit)){
       // insert edges and set reveral information
       e1 = G.new_edge(NL2[prev],NL2[lit]);
       e2 = G.new_edge(NL2[lit],NL2[prev]);
       G.set_reversal(e1,e2);
    }
    
    leda_node v1 = NL2[NL2.first()];
    leda_node v2 = NL2[NL2.last()];
    leda_edge eh1,eh2;
    
    eh1 = G.first_adj_edge(v1);
    
    e1 = G.new_edge(eh1,v2,w,LEDA::before);
    e2 = G.new_edge(v2,v1);    
    G.set_reversal(e1,e2);
    
    pts.clear();
    
    // link vertices ...
    if (! NL.empty()){
      list_item lit2 = NL2.first();
      forall_items(lit,NL){
        v1 = NL[lit]; 
	v2 = NL2[lit2];
	
	eh1 = G.first_adj_edge(v1);
	eh2 = G.last_adj_edge(v2);
	
        e1 = G.new_edge(eh1, v2);
        e2 = G.new_edge(eh2 ,v1);	
	G.set_reversal(e1,e2);
	lit2 = NL2.cyclic_succ(lit2);
      }
    }
    else { NL_first = NL2; }

    NL = NL2; 
        
    yakt = yakt + ystep;
  }
  leda_node NH1 = G.new_node(leda_d3_point(ct.xcoord(),ct.ycoord()-r,ct.zcoord()));
  leda_node NH2 = G.new_node(leda_d3_point(ct.xcoord(),ct.ycoord()+r,ct.zcoord()));
  leda_node niter;
  
  // low vertex ...
  forall(niter,NL_first) {
       e1 = G.new_edge(niter, NH1);
       e2 = G.new_edge(NH1 ,niter);	
       G.set_reversal(e1,e2);    
  }
  
  // high vertex ...
  forall(niter,NL2) {
       leda_edge eh1 = G.last_adj_edge(NH2);
       leda_edge eh2 = G.first_adj_edge(niter);
       e1 = G.new_edge(eh2, NH2, w, LEDA::after);
       if (eh1==NULL) e2 = G.new_edge(NH2 ,niter);	
       else e2= G.new_edge(eh1, niter, w, LEDA::after);
       G.set_reversal(e1,e2);
  }
}

// --------------------------------------------

#undef list

GRAPH<leda_d3_point,int> SPHGR;

// d3 output function

void show_d3_points(geo_scene sc, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GeoEditScene<std::list<Point> >* esc = (GeoEditScene<std::list<Point> > *) sc;
 std::list<Point>& L= esc->get_objref();
 
 if (L.size()==0) return;
 Point p;
 std::list<Point>::const_iterator it = L.begin();
 
 Pointd*  P = new Pointd[L.size()];    
 NT2      coord[dim];
 int      i=0;
 
 for(; it != L.end();++it,i++) { 
   p=*it; 
   H.new_node(convert_to_leda(p));
#if defined(USE_RAT)
   coord[0] = NT2(p.x().to_double()); coord[1]= NT2(p.y().to_double()); 
   coord[2] = NT2(p.z().to_double()); 
#else
   coord[0] = p.x(); coord[1]= p.y(); coord[2] = p.z();    
#endif
   P[i] = Pointd(dim, coord, coord+dim);
 }  

 // construct min sphere ...
 Min_sphere  ms (P, P+L.size());
 
 delete [] P;
 
 // build output graph describing the sphere ...
#if defined(USE_RAT)
 double rad = CGAL_CLIB_STD::sqrt(CGAL::to_double(ms.squared_radius().normalize()));
 double xc = CGAL::to_double((ms.center())[0].normalize());
 double yc = CGAL::to_double((ms.center())[1].normalize());
 double zc = CGAL::to_double((ms.center())[2].normalize());
#else
 double rad = CGAL_CLIB_STD::sqrt(ms.squared_radius());
 double xc = (ms.center())[0];
 double yc = (ms.center())[1];
 double zc = (ms.center())[2];
#endif
 
 leda_d3_point p1(xc-rad,yc,zc), p2(xc+rad,yc,zc), p3(xc,yc+rad,zc), p4(xc,yc,zc+rad);
 leda_d3_point pcenter(xc,yc,zc);
 leda_d3_sphere sph(p1,p2,p3,p4);
 GRAPH<leda_d3_point,int> G;
 leda_vector v;
 leda_d3_point orig(0,0,0);
 double scale;
 leda_node piter;
 v = pcenter - orig; scale = rad/173.2; G=SPHGR;
    
 forall_nodes(piter,G)  {
      leda_d3_point hp(G[piter].xcoord()*scale, G[piter].ycoord()*scale, G[piter].zcoord()*scale);
      G[piter]= hp + v; }
 
 H.join(G);

 leda_node_array<leda_vector> pos(H);
 leda_node vi;
 forall_nodes(vi,H) pos[vi] = H[vi].to_vector();
 W.init(pos); 
}


int main ()
{
  generate_sphere_graph(leda_d3_sphere(leda_d3_point(100,100,100),leda_d3_point(-100,100,100), \
                                  leda_d3_point(100,-100,100),leda_d3_point(-100,-100,-100)),SPHGR);

  geowin_init_default_type((std::list<Point>*)0, leda_string("CGALPoint_3_List"));
 
  std::list<Point> L;
  GeoWin GW("Minimum enclosing sphere in 3d");
  GW.message("To show the enclosing sphere use 'Show d3 output' in Window menu.");

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_d3_fcn(my_scene, show_d3_points);

  GW.edit(my_scene);
  
  return 0;  
}

#endif
