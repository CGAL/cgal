// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Combinatorial_map_save_load.h>

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
  friend void CGAL::read_cmap_attribute_node<Volume_info>
  (const boost::property_tree::ptree::value_type &v,Volume_info &val);

  friend void CGAL::write_cmap_attribute_node<Volume_info>(boost::property_tree::ptree & node,
                                                           const Volume_info& arg);
public:
  Volume_info() : m_color(CGAL::IO::Color(myrandom.get_int(0,256),
                                      myrandom.get_int(0,256),
                                      myrandom.get_int(0,256))),
    m_status( LCC_DEMO_VISIBLE | LCC_DEMO_FILLED )
  {}

  friend bool operator==(const Volume_info& v1, const Volume_info& v2)
  { return v1.m_color==v2.m_color && v1.m_status==v2.m_status; }

  CGAL::IO::Color& color()
  { return m_color; }
  const CGAL::IO::Color& color() const
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
  CGAL::IO::Color m_color;
  char        m_status;
};

typedef CGAL::Linear_cell_complex_traits
<3,CGAL::Exact_predicates_inexact_constructions_kernel> Mytraits;

namespace CGAL
{

template<>
inline void read_cmap_attribute_node<Volume_info>
(const boost::property_tree::ptree::value_type &v,Volume_info &val)
{
  try
  {
    val.m_status = v.second.get<int>("status");
  }
  catch(const std::exception &  )
  {}

  try
  {
    char r = v.second.get<int>("color-r");
    char g = v.second.get<int>("color-g");
    char b = v.second.get<int>("color-b");
    val.m_color = CGAL::IO::Color(r,g,b);
  }
  catch(const std::exception &  )
  {}
}

// Definition of function allowing to save custom information.
template<>
inline void write_cmap_attribute_node<Volume_info>(boost::property_tree::ptree & node,
                                                   const Volume_info& arg)
{
  boost::property_tree::ptree & nValue = node.add("v","");
  nValue.add("status",(int)arg.m_status);
  nValue.add("color-r",(int)arg.m_color.r());
  nValue.add("color-g",(int)arg.m_color.g());
  nValue.add("color-b",(int)arg.m_color.b());
}

}

struct Myitems
{
public:
#ifdef CMAP_WITH_INDEX
  using Use_index=CGAL::Tag_true;
#endif // CMAP_WITH_INDEX

  template < class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute_with_point< Refs > Vertex_attrib;
    typedef CGAL::Cell_attribute< Refs, Volume_info> Volume_attrib;

    typedef std::tuple<Vertex_attrib,void,void,
                       Volume_attrib> Attributes;
  };
};

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3,3,Mytraits,Myitems> LCC;
typedef LCC::Dart_descriptor       Dart_descriptor;
typedef LCC::Dart_const_descriptor Dart_const_descriptor;
typedef LCC::Vertex_attribute      Vertex;

typedef LCC::Point    Point_3;
typedef LCC::Vector   Vector_3;

typedef CGAL::Timer Timer;

struct Vertex_info
{
  Dart_descriptor dh;
  Vector_3 v;
};

struct Face_info {
  bool exist_edge[3];
  bool is_external;
  bool is_process;
};

typedef CGAL::Projection_traits_3<CGAL::Exact_predicates_inexact_constructions_kernel> P_traits;
typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, P_traits> Vb;

typedef CGAL::Triangulation_face_base_with_info_2<Face_info,P_traits> Fb1;

typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>    Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                   TDS;
typedef CGAL::Exact_predicates_tag                                    Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS,
                                                   Itag>              CDT;

struct Scene {
  LCC* lcc;
};

#endif
