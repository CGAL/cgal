// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//                 Mael Rouxel-Labbé

#ifndef CGAL_IO_ALPHA_SHAPE_GEOMVIEW_OSTREAM_3_H
#define CGAL_IO_ALPHA_SHAPE_GEOMVIEW_OSTREAM_3_H

#include <CGAL/license/Alpha_shapes_3.h>

#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/Alpha_shape_3.h>

#include <map>
#include <utility>

// TODO :
// - Check the correctness when dimension < 3.
// - Use the stream color instead of built-in constant/random.
// - If interfaces were more similar, we could think of sharing 2d and 3d ?

//-------------------------------------------------------------------
namespace CGAL {
//-------------------------------------------------------------------

template<class Dt,class EACT>
void
Alpha_shape_3<Dt,EACT>::
show_triangulation_edges(Geomview_stream &gv) const
{
  // To become unordered when CGAL's Point_3/Weighted_point_3 has a hash function
  typedef typename std::map<Point, int> PMap;

  // Used to keep the insertion order in memory
  typedef typename std::vector<typename PMap::iterator> PMapIterVector;

  PMap P;
  PMapIterVector PIV;
  int number_of_points = 0;

  for(Finite_edges_iterator eit=finite_edges_begin(); eit!=finite_edges_end(); ++eit)
  {
    Point p = point(eit->first, eit->second);
    Point q = point(eit->first, eit->third);

    std::pair<typename PMap::iterator, bool> is_p_insert_successful =
      P.insert(std::make_pair(p, number_of_points));

    if(is_p_insert_successful.second)
    {
      ++number_of_points;
      PIV.push_back(is_p_insert_successful.first);
    }

    std::pair<typename PMap::iterator, bool> is_q_insert_successful =
      P.insert(std::make_pair(q, number_of_points));

    if(is_q_insert_successful.second)
    {
      ++number_of_points;
      PIV.push_back(is_q_insert_successful.first);
    }
  }

  // Header
  gv.set_ascii_mode();
  gv << "(geometry " << gv.get_new_id("triangulation_edges")
     << " {appearance {}{ SKEL \n"
     << number_of_points << this->number_of_finite_edges() << "\n";

  // Vertices
  typename PMapIterVector::iterator Pit, Pbegin = PIV.begin(), Pend = PIV.end();
  for(Pit = Pbegin; Pit != Pend; ++Pit)
    gv << this->geom_traits().construct_point_3_object()((*Pit)->first);

  // Finite edges indices
  for(Finite_edges_iterator eit=finite_edges_begin(); eit!=finite_edges_end(); ++eit) {
    gv << 2
       << P[point(eit->first, eit->second)]
       << P[point(eit->first, eit->third)]
       << "\n"; // without color
    // << 4 << drand48() << drand48() << drand48() << 1.0; // random color
  }
}

//-------------------------------------------------------------------
// This function outputs the facets.
template<class Dt,class EACT>
void
Alpha_shape_3<Dt,EACT>::
show_alpha_shape_faces(Geomview_stream &gv) const
{
  // To become unordered when CGAL's Point_3/Weighted_point_3 has a hash function
  typedef typename std::map<Point, int> PMap;

  // Used to keep the insertion order in memory
  typedef typename std::vector<typename PMap::iterator> PMapIterVector;

  typename Alpha_shape_3<Dt,EACT>::Alpha_shape_facets_iterator Flist_it,
    Flist_begin = Alpha_shape_facets_begin(),
    Flist_end = Alpha_shape_facets_end();

  PMap P;
  PMapIterVector PIV;
  int number_of_points = 0;
  int number_of_facets = std::distance(Flist_begin, Flist_end);

  for(Flist_it = Flist_begin; Flist_it != Flist_end; Flist_it++)
  {
    // Collect the points
    // Note that we cannot simply loop over the vertices using e.g.
    // Alpha_shape_vertices_begin() because v->point() might not always
    // be the correct point position (see periodic triangulations for example);
    // we must instead use tr.point(c,i)
    for(int i=0; i<4; i++)
    {
      if(i != Flist_it->second)
      {
        Point p = point(Flist_it->first, i);
        std::pair<typename PMap::iterator, bool> is_insert_successful =
            P.insert(std::make_pair(p, number_of_points));

        if(is_insert_successful.second)
        {
          ++number_of_points;
          PIV.push_back(is_insert_successful.first);
        }
      }
    }
  }

  // Header
  gv.set_binary_mode();
  gv << "(geometry " << gv.get_new_id("alpha_shape")
     << " {appearance {}{ OFF BINARY\n"
     << number_of_points << number_of_facets << 0;

  std::cout << number_of_points << " " << number_of_facets << std::endl;

  // Vertices
  typename PMapIterVector::iterator Pit, Pbegin = PIV.begin(), Pend = PIV.end();
  for(Pit = Pbegin; Pit != Pend; ++Pit)
    gv << this->geom_traits().construct_point_3_object()((*Pit)->first);

  // Finite facets indices
  for(Flist_it = Flist_begin; Flist_it != Flist_end; Flist_it++)
  {
    gv << 3;
    for(int i=0; i<4; i++)
      if(i != Flist_it->second)
        gv << P[point(Flist_it->first, i)];
    gv << 0; // without color
//    gv << 4 << drand48() << drand48() << drand48() << 1.0; // random color
  }
}

//-------------------------------------------------------------------

template < class Dt,class EACT >
Geomview_stream&
operator<<(Geomview_stream &gv, Alpha_shape_3<Dt,EACT>& A)
{
  bool ascii_bak = gv.get_ascii_mode();
  bool raw_bak = gv.set_raw(true);

  if(gv.get_wired())
    A.show_triangulation_edges(gv);
  else
    A.show_alpha_shape_faces(gv);

  // Footer
  gv << "}})";

  gv.set_raw(raw_bak);
  gv.set_ascii_mode(ascii_bak);
  return gv;
}

//-------------------------------------------------------------------
} //namespace CGAL
//-------------------------------------------------------------------

#endif // CGAL_IO_ALPHA_SHAPE_GEOMVIEW_OSTREAM_3_H
