// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Extended_homogeneous.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// source        : nef_2d/Simple_extended_kernel.lw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Extended homogeneous kernel
// ============================================================================
#ifndef CGAL_EXTENDED_HOMOGENEOUS_3_H
#define CGAL_EXTENDED_HOMOGENEOUS_3_H

#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h> 
#include <CGAL/Point_3.h> 
// #include <CGAL/Line_3_Line_3_intersection.h> 
#include <CGAL/squared_distance_3.h> 
#include <CGAL/number_utils.h>
#if (defined( _MSC_VER) && (_MSC_VER <= 1200))
#include <CGAL/Nef_2/Polynomial_MSC.h>
#define Polynomial Polynomial_MSC
#else
#include <CGAL/Nef_2/Polynomial.h>
#endif
#undef _DEBUG
#define _DEBUG 5
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

template <class T> class Extended_homogeneous_3;

template <class RT_>
class Extended_homogeneous_3 : public 
  CGAL::Homogeneous< CGAL::Polynomial<RT_> > { 

public:

  typedef CGAL::Homogeneous< CGAL::Polynomial<RT_> >  Base;
  typedef Extended_homogeneous_3<RT_> Self;

  typedef CGAL::Homogeneous<RT_>                Standard_kernel;

  typedef RT_ Standard_RT;

  typedef typename Standard_kernel::FT Standard_FT;

  typedef typename Standard_kernel::Point_3     Standard_point_3;

  typedef typename Standard_kernel::Segment_3   Standard_segment_3;

  typedef typename Standard_kernel::Line_3      Standard_line_3;

  typedef typename Standard_kernel::Direction_3 Standard_direction_3;

  typedef typename Standard_kernel::Ray_3       Standard_ray_3;

  typedef typename Standard_kernel::Aff_transformation_3 
    Standard_aff_transformation_3;

  typedef typename Base::RT RT;

  typedef typename Base::Point_3      Point_3;

  typedef typename Base::Segment_3    Segment_3;

  typedef typename Base::Direction_3  Direction_3;

  typedef typename Base::Line_3       Line_3;

  enum Point_type { ERROR=0,
                    STANDARD=1,
		    SWFCORNER, NWFCORNER, NEFCORNER, SEFCORNER,
		    SWBCORNER, NWBCORNER, NEBCORNER, SEBCORNER,
                    BOTTOMFFRAME, TOPFFRAME, RIGHTFFRAME, LEFTFFRAME, 
                    BOTTOMBFRAME, TOPBFRAME, RIGHTBFRAME, LEFTBFRAME,
                    SWFRAME, NWFRAME, NEFRAME, SEFRAME,
                    BOTTOMPLANE, TOPPLANE, RIGHTPLANE, LEFTPLANE, FRONTPLANE, BACKPLANE};


  public:
  static Point_3 epoint(const Standard_RT& m1, const Standard_RT& n1, 
			const Standard_RT& m2, const Standard_RT& n2, 
			const Standard_RT& m3, const Standard_RT& n3, 
			const Standard_RT& n4) 
    { return Point_3(RT(n1,m1),RT(n2,m2),RT(n3,m3),RT(n4));}
  


  /* Vorsicht> vereinfacht, da Ortsvektor nicht mit einbezogen */
  Point_3 construct_point(const Standard_line_3& l, Point_type& t) const
  { 
    //   Standard_direction_3 d = l.direction();
    Point_3 res = epoint(l.dx(),0,l.dy(),0,l.dz(),0);
    t = (Point_type) determine_type(res);
    return res;
  }


  template <class Forward_iterator>
  void determine_frame_radius(Forward_iterator start, Forward_iterator end,
                              Standard_RT& R0) const
  { 
    Standard_RT R, mx, nx, my, ny, mz, nz;
    while ( start != end ) {
      Point_3 p = *start++;
      if ( is_standard(p) ) {
        R = CGAL_NTS max(CGAL_NTS abs(p.hx()[0])/p.hw()[0], 
                         CGAL_NTS max(CGAL_NTS abs(p.hy()[0])/p.hw()[0],
				      CGAL_NTS abs(p.hz()[0])/p.hw()[0]));
      } 
      else 
      {
        RT rx = CGAL_NTS abs(p.hx()), ry = CGAL_NTS abs(p.hy()), rz = CGAL_NTS abs(p.hz());
        mx = ( rx.degree()>0 ? rx[1] : Standard_RT(0) ); nx = rx[0];
        my = ( ry.degree()>0 ? ry[1] : Standard_RT(0) ); ny = ry[0];
        mz = ( rz.degree()>0 ? rz[1] : Standard_RT(0) ); nz = rz[0];

	if ( mx > my )      R = CGAL_NTS abs((ny-nx)/(mx-my));
        else if ( mx < my ) R = CGAL_NTS abs((nx-ny)/(my-mx));
        else /* mx == my */ R = CGAL_NTS abs(nx-ny)/(2*p.hw()[0]);
	if ( mx > mz )      R = CGAL_NTS max(R, CGAL_NTS abs((nz-nx)/(mx-mz)));
	else if ( mx < mz)  R = CGAL_NTS max(R, CGAL_NTS abs((nx-nz)/(mz-mx)));
	else                R = CGAL_NTS max(R, CGAL_NTS abs(nx-ny)/(2*p.hw()[0]));
      }
      R0 = CGAL_NTS max(R+1,R0);
    }
  }

  
  Point_3 construct_point(const Standard_point_3& p) const
    { return Point_3(p.hx(), p.hy(), p.hz(), p.hw()); }

  Point_3 construct_point(const Standard_point_3& p1, 
                          const Standard_point_3& p2, 
                          const Standard_point_3& p3,
                          Point_type& t) const
  { return construct_point(Standard_line_3(p1,p2,p3),t); }

  Point_3 construct_point(const Standard_line_3& l) const
  { Point_type dummy; return construct_point(l,dummy); }

  Point_3 construct_point(const Standard_point_3& p1, 
                          const Standard_point_3& p2,
                          const Standard_point_3& p3) const
  { return construct_point(Standard_line_3(p1,p2,p3)); }

  Point_3 construct_point(const Standard_point_3& p, 
                          const Standard_direction_3& d) const
  { return construct_point(Standard_line_3(p,d)); }

  Point_3 construct_opposite_point(const Standard_line_3& l) const
  { Point_type dummy; return construct_point(l.opposite(),dummy); }


  static Point_type type(const Point_3& p)
  {
    CGAL_assertion(p.hx().degree()>=0 && p.hy().degree()>=0 && p.hz().degree()>=0 );
    CGAL_assertion(p.hw().degree()==0);
 
   if (p.hx().degree() == 0 && p.hy().degree() == 0 && p.hz().degree() == 0) 
      return STANDARD;
    
    RT r[3];
    r[0] = p.hx();
    r[1] = p.hy();
    r[2] = p.hz();

    int s[3];
    for(int i=0;i<3;i++) {
      s[i] = CGAL_NTS sign(r[i]);
      if(s[i] < 0) r[i] = -r[i];
    }
    
    int on_frame = 1;
    int ref = 0;
    for(int i=1,bit = 2;i<3;++i) {
      if(r[i]>r[ref]){
	on_frame = bit;
	ref = i;
      }
      else if(r[i]==r[ref])
	on_frame += bit;
      bit*=2;
    }

    switch(on_frame) 
    {
      case 7:
	CGAL_assertion(r[0]==r[1]);
	CGAL_assertion(r[1]==r[2]);

	if(s[0]>0)
	  if(s[1]>0)
	    if(s[2]>0) return NEFCORNER;
	    else       return SEFCORNER;
	  else
	    if(s[2]>0) return NEBCORNER;
	    else       return SEBCORNER;
	else
	  if(s[1]>0)
	    if(s[2]>0) return NWFCORNER;
	    else       return SWFCORNER;
	  else
	    if(s[2]>0) return NWBCORNER;
	    else       return SWBCORNER;
	break;

      case 6:
	CGAL_assertion(r[2]==r[1]);
	CGAL_assertion(r[2]>r[0]);

	if(s[1]>0)
	  if(s[2]>0) return TOPFFRAME;
	  else       return BOTTOMFFRAME;
	else
	  if(s[2]>0) return TOPBFRAME;
	  else       return BOTTOMBFRAME;
	break;
 
     case 5:
	CGAL_assertion(r[2]==r[0]);
	CGAL_assertion(r[2]>r[1]);	

	if(s[0]>0)
	  if(s[2]>0) return NEFRAME;
	  else       return SEFRAME;
	else
	  if(s[2]>0) return NWFRAME;
	  else       return SWFRAME;
	break;

     case 4:
	CGAL_assertion(r[2]>r[0]);
	CGAL_assertion(r[2]>r[1]);	

	if(s[2]>0) return TOPPLANE;
	else       return BOTTOMPLANE;
	break;

     case 3:
	CGAL_assertion(r[1]==r[0]);
	CGAL_assertion(r[1]>r[2]);	

	if(s[0]>0)
	  if(s[1]>0) return RIGHTFFRAME;
	  else       return LEFTFFRAME;
	else
	  if(s[1]>0) return RIGHTBFRAME;
	  else       return LEFTBFRAME;
	break;

     case 2:
	CGAL_assertion(r[1]>r[0]);
	CGAL_assertion(r[1]>r[2]);	

	if(s[1]>0) return FRONTPLANE;
	else       return BACKPLANE;

	break;

     case 1:
	CGAL_assertion(r[0]>r[1]);
	CGAL_assertion(r[0]>r[2]);	

	if(s[0]>0) return RIGHTPLANE;
	else       return LEFTPLANE;

	break;

      default: 
	CGAL_assertion_msg(0,"EPoint type not correct!");
	return ERROR;
    }
  }


  static bool is_standard(const Point_3& p)
  { return (type(p)==STANDARD);  }

  Standard_point_3 standard_point(const Point_3& p) const
  { CGAL_assertion(type(p)==STANDARD);
    CGAL_assertion(p.hw() > RT(0));
    return Standard_point_3(p.hx()[0],p.hy()[0],p.hz()[0],p.hw()[0]);
  }

  Standard_line_3 standard_line(const Point_3& p) const
  { CGAL_assertion(type(p)!=STANDARD);
    RT hx = p.hx(), hy = p.hy(), hz = p.hz(), hw = p.hw();
    Standard_RT dx,dy,dz;
    if (hx.degree()>0) dx=hx[1]; else dx=0;
    if (hy.degree()>0) dy=hy[1]; else dy=0;
    if (hz.degree()>0) dz=hz[1]; else dz=0;
    Standard_point_3 p0(hx[0],hy[0],hz[0],hw[0]);
    Standard_point_3 p1(hx[0]+dx,hy[0]+dy,hz[0]+dz,hw[0]);
    return Standard_line_3(p0,p1);
  }

  Standard_ray_3 standard_ray(const Point_3& p) const
  { CGAL_assertion(type(p)!=STANDARD);
    Standard_line_3 l = standard_line(p);
    Standard_direction_3 d = l.direction();
    Standard_point_3 q = l.point(0);
    return Standard_ray_3(q,d);
  }

  Point_3 NEF() const { return construct_point(Standard_line_3( 1, 1, 1,0)); }
  Point_3 SEF() const { return construct_point(Standard_line_3( 1, 1,-1,0)); }
  Point_3 NWF() const { return construct_point(Standard_line_3(-1, 1, 1,0)); }
  Point_3 SWF() const { return construct_point(Standard_line_3(-1, 1,-1,0)); }
  Point_3 NEB() const { return construct_point(Standard_line_3( 1,-1, 1,0)); }
  Point_3 SEB() const { return construct_point(Standard_line_3( 1,-1,-1,0)); }
  Point_3 NWB() const { return construct_point(Standard_line_3(-1,-1, 1,0)); }
  Point_3 SWB() const { return construct_point(Standard_line_3(-1,-1,-1,0)); }


  //  Line_3 upper() const { return construct_line(NW(),NE()); }
  //  Line_3 lower() const { return construct_line(SW(),SE()); }
  //  Line_3 left()  const { return construct_line(SW(),NW()); }
  //  Line_3 right() const { return construct_line(SE(),NE()); }


  Point_3 source(const Segment_3& s) const
  { typename Base::Construct_vertex_3 _source = construct_vertex_3_object();
    return _source(s,0); }

  Point_3 target(const Segment_3& s) const
  { typename Base::Construct_vertex_3 _target = construct_vertex_3_object();
    return _target(s,1); }

  Segment_3 construct_segment(const Point_3& p, const Point_3& q) const
  { typename Base::Construct_segment_3 _segment =
      construct_segment_3_object();
    return _segment(p,q); }

  void simplify(Point_3& p) const
  { TRACEN("simplify("<<p<<")");
    RT x=p.hx(), y=p.hy(), z=p.hz(), w=p.hw();
    RT common = x.is_zero() ? z : gcd(x,z);
    common = y.is_zero() ? common : gcd(common,y);
    common = gcd(common,w);
    p = Point_3(x/common,y/common,z/common,w/common);
    TRACEN("canceled="<<p);
  }

  Line_3 construct_line(const Standard_line_3& l)  const
    { return Line_3(l.point(1),l.point(2)); }

  Line_3 construct_line(const Point_3& p1, const Point_3& p2) const
    { return Line_3(p1,p2); }

  int orientation(const Segment_3& s, const Point_3& p) const
  { typename Base::Orientation_3 _orientation =
      orientation_3_object();
    return static_cast<int> ( _orientation(source(s),target(s),p) ); 
  }

  int orientation(const Point_3& p1, const Point_3& p2, const Point_3& p3) const
  { typename Base::Orientation_3 _orientation =
      orientation_3_object();
    return static_cast<int> ( _orientation(p1,p2,p3) ); 
  }

  bool left_turn(const Point_3& p1, const Point_3& p2, const Point_3& p3) 
  const
  { return orientation(p1,p2,p3) > 0; }
   
  bool is_degenerate(const Segment_3& s) const
  { typename Base::Is_degenerate_3 _is_degenerate =
      is_degenerate_3_object();
    return _is_degenerate(s); 
  }

  int compare_xyz(const Point_3& p1, const Point_3& p2) const 
  { 
    typename Base::Compare_xyz_3 _compare_xyz = compare_xyz_3_object();
    return static_cast<int>( _compare_xyz(p1,p2) );
  }
  
  int compare_x(const Point_3& p1, const Point_3& p2) const
  { typename Base::Compare_x_3 _compare_x =
      compare_x_3_object();
    return static_cast<int>( _compare_x(p1,p2) );
  }

  int compare_y(const Point_3& p1, const Point_3& p2) const
  { typename Base::Compare_y_3 _compare_y =
      compare_y_3_object();
    return static_cast<int>( _compare_y(p1,p2) );
  }

  int compare_z(const Point_3& p1, const Point_3& p2) const
  { typename Base::Compare_z_3 _compare_z =
      compare_z_3_object();
    return static_cast<int>( _compare_z(p1,p2) );
  }

  Point_3 intersection(
    const Segment_3& s1, const Segment_3& s2) const
  { typename Base::Intersect_3 _intersect =
      intersect_3_object();
    typename Base::Construct_line_3 _line =
      construct_line_3_object();
    Point_3 p; 
    CGAL::Object result =
      _intersect(_line(s1),_line(s2));
    if ( !CGAL::assign(p, result) )
    CGAL_assertion_msg(false,"intersection: no intersection.");
    simplify(p);
    return p;
  }

  Direction_3 construct_direction(
    const Point_3& p1, const Point_3& p2) const
  { typename Base::Construct_direction_3 _direction =
      construct_direction_3_object();
    return _direction(construct_line(p1,p2)); 
  }

  /*
  bool strictly_ordered_ccw(const Direction_3& d1, 
    const Direction_3& d2, const Direction_3& d3) const
  { 
??    if ( d1 < d2 )  return ( d2 < d3 )||( d3 <= d1 );
??    if ( d1 > d2 )  return ( d2 < d3 )&&( d3 <= d1 );
    return false;
  }
  */

  bool strictly_ordered_along_line(
    const Point_3& p1, const Point_3& p2, const Point_3& p3) const
  { typename Base::Are_strictly_ordered_along_line_3 _ordered =
      are_strictly_ordered_along_line_3_object();
    return _ordered(p1,p2,p3);
  }

  bool contains(const Segment_3& s, const Point_3& p) const
  { typename Base::Has_on_3 _contains = has_on_3_object();
    return _contains(s,p);
  }

  bool first_pair_closer_than_second(
    const Point_3& p1, const Point_3& p2, 
    const Point_3& p3, const Point_3& p4) const
  { return ( squared_distance(p1,p2) < squared_distance(p3,p4) ); }

  /*
  void scale_first_by_second(RT& r1, RT& r2, RT& r3, RT& w) const
  { CGAL_assertion(w.degree()==0&&w!=RT(0)&& r2[1]!=Standard_RT(0));
    Standard_RT w_res = w[0]*r2[1];
    int sm2 = CGAL_NTS sign(r2[1]);
    RT r2_res = RT(Standard_RT(0),sm2 * w_res); 
    RT r1_res = RT(r2[1]*r1[0]-r1[1]*r2[0], w[0]*r1[1]*sm2);
    r1 = r1_res; r2 = r2_res; w = w_res;
  }
 

  Point_3 transform(const Point_3& p, 
                    const Standard_aff_transformation_3& t) const
  {
    RT tpx = t.homogeneous(0,0)*p.hx() + t.homogeneous(0,1)*p.hy() +
      t.homogeneous(0,2)*p.hz()+t.homogeneous(0,3)*p.hw();
    RT tpy = t.homogeneous(1,0)*p.hx() + t.homogeneous(1,1)*p.hy() +
      t.homogeneous(1,2)*p.hz()+t.homogeneous(1,3)*p.hw();
    RT tpw = t.homogeneous(2,2)*p.hw();
    if ( is_standard(p) ) {
      Point_3 res(tpx,tpy,tpz,tpw); simplify(res);
      return res;
    } 
    RT tpxa = CGAL_NTS abs(tpx);
    RT tpya = CGAL_NTS abs(tpy);
    RT tpza = CGAL_NTS abs(tpz);
??    if ( tpxa > tpya ) {
??      scale_first_by_second(tpy,tpx,tpw);
??    } else { // tpxa <= tpya
??      scale_first_by_second(tpx,tpy,tpw);
    }
    Point_3 res(tpx,tpy,tpz, tpw); simplify(res);
    return res;
  }
  */

  const char* output_identifier() const { return "Extended_homogeneous"; }



};


#undef Polynomial
CGAL_END_NAMESPACE
#endif // CGAL_EXTENDED_HOMOGENEOUS_3_H
