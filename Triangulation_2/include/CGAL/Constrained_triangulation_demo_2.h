// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mariette Yvinec

#ifndef CGAL_CONSTRAINED_TRIANGULATION_DEMO_2_H
#define CGAL_CONSTRAINED_TRIANGULATION_DEMO_2_H

#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Constrained_triangulation_sweep_demo_2.h>
#include <CGAL/IO/Window_stream.h>

CGAL_BEGIN_NAMESPACE

template < class Gt,class Tds>
class Constrained_triangulation_demo_2
  : public Constrained_triangulation_2<Gt,Tds>
{
public:
  typedef Constrained_triangulation_2<Gt,Tds> Constrained_triangulation;
  typedef Constrained_triangulation_sweep_demo_2<Gt,Tds>  Sweep_demo;
  typedef typename Gt::Segment Segment;
  typedef Window_stream Window_stream;

Constrained_triangulation_demo_2() :  
  Constrained_triangulation_2<Gt,Tds>() {} 
  
Constrained_triangulation_demo_2(const Gt& gt=Gt()) 
  : Constrained_triangulation_2<Gt,Tds>(gt) {}
  

Constrained_triangulation_demo_2(Window_stream& W,
				     list<Constraint>& lc, const Gt& gt=Gt()) 
  : Constrained_triangulation_2<Gt,Tds>(gt)
{
  Sweep_demo  sweep(W,lc, gt);
  Constrained_triangulation_2<Gt,Tds> Tr( sweep.vertex(), gt);
  swap(Tr);
  CGAL_triangulation_postcondition( is_valid() );
}

};

CGAL_END_NAMESPACE

#endif //CGAL_CONSTRAINED_TRIANGULATION_DEMO_2_H
