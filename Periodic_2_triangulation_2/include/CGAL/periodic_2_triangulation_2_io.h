// DO NOT REDISTRIBUTE
//
// This file is not yet part of CGAL (www.cgal.org), but will be.
//
// $URL$
// $Id$
// SPDX-License-Identifier:
//
// Author(s)     : GPT-Authors

#ifndef CGAL_PERIODIC_2_TRIANGULATION_2_IO_H
#define CGAL_PERIODIC_2_TRIANGULATION_2_IO_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/array.h>

#include <iostream>
#include <utility>

namespace CGAL {

template <class Stream, class T2>
bool write_DT2_to_OFF(Stream& os, const T2& t2)
{
  typedef typename T2::Point    Point;

  os << "OFF\n";

  std::map<Point, int> ids;
  std::vector<std::array<int, 3> > faces;
  std::vector<Point> points;

  int pid = 0;
  for(auto fit=t2.finite_faces_begin(); fit!=t2.finite_faces_end(); ++fit)
  {
    std::array<int, 3> face;
    for(int i=0; i<3; ++i)
    {
      std::pair<typename std::map<Point, int>::iterator, bool> itb =
          ids.insert(std::make_pair(t2.point(fit, i), pid));
      if(itb.second)
      {
        points.push_back(t2.point(fit, i));
        ++pid;
      }

      face[i] = itb.first->second;
    }
    faces.push_back(face);
  }

  CGAL_assertion(faces.size() == t2.number_of_faces());

  os << ids.size() << " " << t2.number_of_faces() << " 0\n";

  for(const Point& pt : points)
    os << pt << " 0" << std::endl;

  for(const auto& face : faces)
    os << "3 " << face[0] << " " << face[1] << " " << face[2] << std::endl;

  os << std::endl;

  return true;
}

template <class Stream, class P2T2>
bool write_PD2T2_to_OFF(Stream &os, P2T2& p2t2)
{
  typedef typename P2T2::Point    Point;

  os << "OFF\n";

  std::map<Point, int> ids;
  std::vector<std::array<int, 3> > faces;
  std::vector<Point> points;

  int pid = 0;
  for(auto fit=p2t2.faces_begin(); fit!=p2t2.faces_end(); ++fit)
  {
    std::array<int, 3> face;
    for(int i=0; i<3; ++i)
    {
      auto itb = ids.insert(std::make_pair(p2t2.point(fit, i), pid));
      if(itb.second)
      {
        points.push_back(p2t2.point(fit, i));
        ++pid;
      }

      face[i] = itb.first->second;
    }
    faces.push_back(face);
  }

  os << ids.size() << " " << faces.size() << " 0\n";

  for(const Point& pt : points)
    os << pt << " 0" << std::endl;

  for(const auto& face : faces)
    os << "3 " << face[0] << " " << face[1] << " " << face[2] << std::endl;

  os << std::endl;

  return true;
}

template <class Stream, class GPT>
bool write_GP2T2_to_OFF(Stream &os, GPT& gp2t2)
{
  if(gp2t2.is_1_cover())
    return write_PD2T2_to_OFF(os, gp2t2.p2dt2);
  else
    return write_DT2_to_OFF(os, gp2t2.dt2);
}

} // namespace CGAL

#endif //CGAL_PERIODIC_2_TRIANGULATION_2_IO_H
