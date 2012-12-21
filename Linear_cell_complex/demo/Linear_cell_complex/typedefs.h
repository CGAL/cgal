// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
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
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/Color.h>
#include <CGAL/Timer.h>
#include <CGAL/Random.h>

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>

// Global random
extern CGAL::Random myrandom;

// Use to define properties on volumes.
#define LCC_DEMO_VISIBLE 1 // if not visible => hidden
#define LCC_DEMO_FILLED  2 // if not filled, wireframe

class Volume_info
{
public:
  Volume_info() : m_color(CGAL::Color(myrandom.get_int(0,256),
                                      myrandom.get_int(0,256),
                                      myrandom.get_int(0,256))),
    m_status( LCC_DEMO_VISIBLE | LCC_DEMO_FILLED )
  {}

  CGAL::Color& color()
  { return m_color; }
  const CGAL::Color& color() const
  { return m_color; }

  std::string color_name() const
  {
    std::ostringstream ss;
    ss<<std::setfill('0');
    ss<<"#"<<std::hex<<std::setw(2)<<(int)m_color.red()
     <<std::setw(2)<<(int)m_color.green()<<std::setw(2)<<(int)m_color.blue();
    return ss.str();
  }

  bool is_visible() const
  { return (m_status & LCC_DEMO_VISIBLE)!=0; }
  bool is_filled() const
  { return (m_status & LCC_DEMO_FILLED)!=0; }
  bool is_filled_and_visible() const
  { return is_filled() && is_visible(); }

  void set_visible(bool val=true)
  {
    if ( is_visible()==val ) return;
    if ( val ) m_status = m_status | LCC_DEMO_VISIBLE;
    else       m_status = m_status ^ LCC_DEMO_VISIBLE;
  }
  void set_filled(bool val=true)
  {
    if ( is_filled()==val ) return;
    if ( val ) m_status = m_status | LCC_DEMO_FILLED;
    else       m_status = m_status ^ LCC_DEMO_FILLED;
  }

  void negate_visible()
  { set_visible(!is_visible()); }
  void negate_filled()
  { set_filled(!is_filled()); }

  private:
  CGAL::Color m_color;
  char        m_status;
};

class Myitems
{
public:
  template < class Refs >
  struct Dart_wrapper 
  {
    typedef CGAL::Dart<3, Refs > Dart;
    
    typedef CGAL::Cell_attribute_with_point< Refs > Vertex_attrib;
    typedef CGAL::Cell_attribute< Refs, Volume_info> Volume_attrib;
    
    typedef CGAL::cpp11::tuple<Vertex_attrib,void,void,
                               Volume_attrib> Attributes;
  };
};

typedef CGAL::Linear_cell_complex_traits
<3,CGAL::Exact_predicates_inexact_constructions_kernel> Mytraits;

typedef CGAL::Linear_cell_complex<3,3,Mytraits,Myitems> LCC;
typedef LCC::Dart_handle      Dart_handle;
typedef LCC::Vertex_attribute Vertex;

typedef LCC::Point    Point_3;
typedef LCC::Vector   Vector_3;

typedef CGAL::Timer Timer;

struct Scene {
  LCC* lcc;
};

#endif
