// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
#ifndef CGAL_EXTENDED_CARTESIAN_3_H
#define CGAL_EXTENDED_CARTESIAN_3_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h> 
#include <CGAL/number_utils.h>
#include <CGAL/Nef_2/Nef_polynomial.h>
#undef _DEBUG
#define _DEBUG 5
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

template <class RT_>
class Extended_cartesian_3 : public 
  CGAL::Cartesian< CGAL::Nef_polynomial<RT_> > { 

public:
  typedef CGAL::Cartesian< CGAL::Nef_polynomial<RT_> >  Base;
  typedef Extended_cartesian_3<RT_>                     Self;

  typedef CGAL::Cartesian<RT_>                Standard_kernel;
  typedef RT_                                   Standard_RT;
  typedef typename Standard_kernel::FT          Standard_FT;
  typedef typename Standard_kernel::Point_3     Standard_point_3;
  typedef typename Standard_kernel::Segment_3   Standard_segment_3;
  typedef typename Standard_kernel::Line_3      Standard_line_3;
  typedef typename Standard_kernel::Ray_3       Standard_ray_3;
  typedef typename Standard_kernel::Aff_transformation_3 
                                                Standard_aff_transformation_3;

  typedef typename Base::RT           RT;
  typedef typename Base::Point_3      Point_3;
  typedef typename Base::Segment_3    Segment_3;
  typedef typename Base::Line_3       Line_3;

  public:
  static Point_3 epoint(const Standard_RT& m1, const Standard_RT& n1, 
			const Standard_RT& m2, const Standard_RT& n2, 
			const Standard_RT& m3, const Standard_RT& n3, 
			const Standard_RT& n4) 
    { return Point_3(RT(n1,m1),RT(n2,m2),RT(n3,m3),RT(n4));}
  
  
  static bool is_standard(const Point_3& p) { 
    CGAL_assertion(p.hx().degree()>=0 && p.hy().degree()>=0 && p.hz().degree()>=0 );
    CGAL_assertion(p.hw().degree()==0);
    
    if (p.hx().degree() == 0 && p.hy().degree() == 0 && p.hz().degree() == 0) 
      return true;
    return false;
  }

  Standard_point_3 standard_point(const Point_3& p) const
  { CGAL_assertion(!is_standard(p));
    CGAL_assertion(p.hw() > RT(0));
    return Standard_point_3(p.hx()[0],p.hy()[0],p.hz()[0],p.hw()[0]);
  }

  Standard_line_3 standard_line(const Point_3& p) const
  { CGAL_assertion(!is_standard(p));
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
  { CGAL_assertion(!is_standard(p));
    Standard_line_3 l = standard_line(p);
    Standard_point_3 q = l.point(0);
    return Standard_ray_3(q,l);
  }

  Segment_3 construct_segment(const Point_3& p, const Point_3& q) const
  { typename Base::Construct_segment_3 _segment =
      construct_segment_3_object();
    return _segment(p,q); 
  }
};

CGAL_END_NAMESPACE
#endif // CGAL_EXTENDED_CARTESIAN_3_H
