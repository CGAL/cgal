// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)         : Oren Nechushtan <theoren@math.tau.ac.il>
#ifndef CGAL_TRAPEZOIDAL_DECOMPOSITION_2_IOSTREAM_H
#define CGAL_TRAPEZOIDAL_DECOMPOSITION_2_IOSTREAM_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#ifndef CGAL_TRAPEZOIDAL_DECOMPOSITION_2_H
#include <CGAL/Arr_point_location/Trapezoidal_decomposition_2.h>
#endif

namespace CGAL {

template < class Traits>
std::ostream& operator<<(
                         std::ostream &out,
                         const typename Trapezoidal_decomposition_2<Traits>::In_face_iterator& it)
{
  out << "In_face_iterator(";
  if (it.operator->()) out << *it; else out <<"end";
  return out << ")" << std::flush;
}

template < class Traits>
std::ostream& write(
                    std::ostream &out,const Trapezoidal_decomposition_2<Traits>& td)
{
  return out << td;
}

template < class Traits>
std::ostream& operator<<(
                         std::ostream &out,const Trapezoidal_decomposition_2<Traits>& td)
{
  return write(out,td.dag_root(),td.get_traits());
}

template < class Traits>
std::ostream& write(std::ostream &out,const Td_X_trapezoid<Traits>& t,
                    const Traits& traits,bool validate=true)
{
  typedef Trapezoidal_decomposition_2<Traits> TD;
  typedef Td_X_trapezoid<Traits> X_trapezoid;
  typedef X_trapezoid* pointer;
  bool pad=false;

  //MICHAL: this may fail for inactive trapezoids

  out << "(";
  if (!t.is_on_left_boundary()) out << t.left(); else out << "-oo"; //MICHAL: this may fail for removed point
  out << ",";
  if (!t.is_on_right_boundary()) out << t.right(); else out << "+oo";
  out << "," << std::flush;
  if (!t.is_on_bottom_boundary()) out << t.bottom(); else out << "-oo";
  out << ",";
  if (!t.is_on_top_boundary()) out << t.top(); else out << "+oo";
  out << ",neighbours(" << std::flush;

  // debug neighbours equivalence relation
  int max_size=4+1;
  int null_size=2,size=null_size,i,j;

  X_trapezoid* value[] =
    {
      0,
      (X_trapezoid*)CGAL_TD_DELETE_SIGNATURE,
      t.lb(),
      t.lt(),
      t.rb(),
      t.rt()
    };
  typedef char debug_string[256];
  debug_string name[]=
    {
      "none",
      "deleted",
      "lb",
      "lt",
      "rb",
      "rt"
    };
  for (j=null_size;j<=max_size;j++)
    {
      for (i=0;i<size;i++)
        if (value[j]==value[i])
          {
            std::strcat(name[i],"=");
            std::strcat(name[i],name[j]);
            break;
          }
      if (i==size)
        {
          value[size]=value[j];
          std::strcpy(name[size++],name[j]);
        }
    }
  if (size==null_size) {out << "none";}
  for(j=null_size;j<size;j++)
    {
      if (pad) out << " ";
      else pad=true;
      out << name[j];
      // identify neighbours
      if (traits.is_td_vertex(t) && value[j])
        out << "=" << value[j]->top();
    }
  out << ")" << std::flush;

  if (t.is_active())
    {
      if (traits.is_td_trapezoid(t))
        {
          if (t.is_unbounded())
            {
              out << ",U=";
              if (t.is_on_left_boundary())   out << ")";
              if (t.is_on_bottom_boundary()) {
                if (t.is_on_top_boundary())
                  out << "/\\/";
                else
                  out << "/\\";
              }
              else if (t.is_on_top_boundary())    out << "\\/";
              if (t.is_on_right_boundary())  out << "(";
            }
          else
            out << ",T";
        }
      else if (traits.is_td_edge(t))
        out << ",C";
      else // if (t.is_degenerate_point())
        out << ",P";
    }
  else
    out << ",D";

  /*        Calling t.is_valid requires the traits to be initialized with the proper
    bounding box */
  if (validate)
    {
      if (t.is_valid(&traits))
        out << ",+";
      else
        out << ",-";
    }

  out << ")" << std::flush;

  return out;
}

template < class Traits>
std::ostream& operator<<(std::ostream &out,const Td_X_trapezoid<Traits>& t)
{
  typedef Trapezoidal_decomposition_2<Traits> TD;
  typedef Td_X_trapezoid<Traits> X_trapezoid;
  typedef X_trapezoid* pointer;
  Traits traits;

  //MICHAL: this may fail for inactive trapezoids

  out << "(";
  if (!t.is_on_left_boundary()) out << t.left(); else out << "-oo";
  out << ",";
  if (!t.is_on_right_boundary()) out << t.right(); else out << "+oo";
  out << "," << std::flush;
  if (!t.is_on_bottom_boundary()) out << t.bottom(); else out << "-oo";
  out << ",";
  if (!t.is_on_top_boundary()) out << t.top(); else out << "+oo";
  out << ",neighbours(" << std::flush;

  // debug neighbours equivalence relation
  int max_size=4+1;
  int null_size=2,size=null_size,i,j;

  X_trapezoid* value[] =
    {
      0,
      (X_trapezoid*)CGAL_TD_DELETE_SIGNATURE,
      t.lb(),
      t.lt(),
      t.rb(),
      t.rt()
    };
  typedef char debug_string[256];
  debug_string name[]=
    {
      "none",
      "deleted",
      "lb",
      "lt",
      "rb",
      "rt"
    };
  for (j=null_size;j<=max_size;j++)
    {
      for (i=0;i<size;i++)
        if (value[j]==value[i])
          {
            std::strcat(name[i],"=");
            std::strcat(name[i],name[j]);
            break;
          }
      if (i==size)
        {
          value[size]=value[j];
          std::strcpy(name[size++],name[j]);
        }
    }
  if (size==null_size) {out << "none";}
  for(j=null_size;j<size;j++)
    {
      out << name[j];
      // identify neighbours
      if (traits.is_td_vertex(t) && value[j])
        out << "=" << value[j]->top();
      out << " ";
    }
  out << ")" << std::flush;

  if (t.is_active())
    {
      if (traits.is_td_trapezoid(t))
        {
          if (t.is_unbounded())
            {
              out << ",U";
              if (t.is_on_left_boundary())
                out << (char) 174;
              if (t.is_on_bottom_boundary())
                out << (char) 25;
              if (t.is_on_top_boundary())
                out << (char) 24;
              if (t.is_on_right_boundary())
                out << (char) 175;
            }
          else
            out << ",T";
        }
      else if (traits.is_td_edge(t))
        out << ",C";
      else // if (t.is_degenerate_point())
        out << ",P";
    }
  else
    out << ",D";

  /*        Calling t.is_valid requires the traits to be initialized with the proper
    bounding box
    if (traits)
    {
    if (t.is_valid(traits))
    out << ",+";
    else
    out << ",-";
    }
  */
  out << ")" << std::flush;

  return out;
}

} //namespace CGAL

#endif //CGAL_TRAPEZOIDAL_DECOMPOSITION_2_IOSTREAM_H
