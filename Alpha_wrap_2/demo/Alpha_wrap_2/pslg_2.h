// Copyright (c) 2019-2020 X, The Moonshot Factory (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Author(s)     : Pierre Alliez pierre.alliez@inria.fr
//               : Michael Hemmer mhsaar@gmail.com
//               : Cedric Portaneri cportaneri@gmail.com
//
#ifndef CGAL_ALPHA_WRAP_2_DEMO_PSLG_H
#define CGAL_ALPHA_WRAP_2_DEMO_PSLG_H

#include <list>

namespace CGAL {
namespace Alpha_wraps_2 {
namespace internal {

template <class GeomTraits>
class Component
  : public std::vector<typename GeomTraits::Point_2>
{
  typedef GeomTraits Geom_traits;
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_2 Point;
  typedef typename Geom_traits::Vector_2 Vector;

  typedef typename std::vector<Point> Base;
  typedef typename Base::iterator iterator;

public:
  bool is_closed = false;

public:
  Component() {}
  Component(const Component& c) : Base(c) { }

  virtual ~Component() {}

 public:
  void set_is_closed(bool b) { is_closed = b; }

  void resample(const FT d) {
    const FT sqd = CGAL::square(d);
    while(resample_once(sqd)) { }
  }

  bool resample_once(const FT sqd)
  {
    for(iterator it = this->begin() ; it != (this->end()-1); it++) {
      const Point& curr = *it;
      const Point& next = next_loop(it);
      if(sqdistance(curr,next) > sqd) {
        this->insert(it+1, CGAL::midpoint(next,curr));
        return true;
      }
    }
    return false;
  }

  FT sqdistance(const Point& a, const Point& b) {
    return CGAL::squared_distance(a, b);
  }

  FT distance(const Point& a, const Point& b) {
    return CGAL::approximate_sqrt(sqdistance(a, b));
  }

  const Point& next_loop(iterator it) {
    it++;
    if(it == this->end())
      return *this->begin();
    else
      return *it;
  }

  const Point& prev_loop(iterator it) {
    if(it == this->begin())
      it = this->end();
    return *(--it);
  }

  void smooth(unsigned int nb_iter) {
    for(unsigned int i=0;i<nb_iter;i++)
      smooth();
  }

  void smooth()
  {
    if(this->size() < 3)
      return;

    Component tmp;
    for(iterator it = this->begin(); it != this->end(); it++)
    {
      Vector v1 = prev_loop(it) - CGAL::ORIGIN;
      Vector v2 = *it - CGAL::ORIGIN;
      Vector v3 = next_loop(it) - CGAL::ORIGIN;
      tmp.push_back(CGAL::ORIGIN + (v1+2*v2+v3)/4);
    }

    if(is_closed) {
      tmp.set_is_closed(true);
      tmp.push_back(tmp[0]);
    }

    *this = tmp;
  }

  void range(FT* xrange, FT* yrange)
  {
    for(iterator it = this->begin(); it != this->end(); it++)
    {
      const Point& p = *it;
      xrange[0] = std::min(xrange[0],p.x());
      xrange[1] = std::max(xrange[1],p.x());
      yrange[0] = std::min(yrange[0],p.y());
      yrange[1] = std::max(yrange[1],p.y());
    }
  }

  void normalize(const FT xmin, const FT ymin, const FT range)
  {
    for(iterator it = this->begin(); it != this->end(); it++)
    {
      Point& p = *it;
      p = Point((p.x()-xmin)/range,
                (p.y()-ymin)/range);
    }
  }

  void translate(const FT dx, const FT dy)
  {
    for(iterator it = this->begin(); it != this->end(); it++)
    {
      Point& p = *it;
      p = p + Vector(dx,dy);
    }
  }

  FT length() {
    FT total_length = 0.;
    for(iterator it = this->begin() ; it != this->end()-1; it++) {
      const Point& curr = *it;
      const Point& next = next_loop(it);
      total_length += distance(curr,next);
    }
    return total_length;
  }
};

template <typename GeomTraits>
class Pslg
  : public std::vector<Component<GeomTraits> >
{
public:
  typedef GeomTraits Geom_traits;
  typedef typename GeomTraits::FT FT;

  typedef internal::Component<GeomTraits> Component;
  typedef typename std::vector<Component> Base;
  typedef typename Base::iterator iterator;
  typedef typename Base::const_iterator citerator;

public:
  Pslg() {}
  virtual ~Pslg() {}

public:
  void normalize() {
    FT xrange[3] = {1e38,-1e38,0};
    FT yrange[3] = {1e38,-1e38,0};

    for(iterator c = this->begin(); c != this->end(); ++c) {
      Component& component = *c;
      component.range(xrange,yrange);
    }
    xrange[2] = xrange[1]-xrange[0];
    yrange[2] = yrange[1]-yrange[0];
    FT range = std::max(xrange[2],yrange[2]);

    for(iterator c = this->begin();
        c != this->end();
        c++) {
      Component& component = *c;
      component.normalize(xrange[0],yrange[0],range);
    }
  }

  void compute_range(FT& xmin, FT& xmax,
                     FT& ymin, FT& ymax,
                     const FT stretch)
  {
    FT xrange[3] = {1e38,-1e38,0};
    FT yrange[3] = {1e38,-1e38,0};

    for(iterator c = this->begin(); c != this->end(); ++c) {
      Component& component = *c;
      component.range(xrange,yrange);
    }

    xrange[2] = xrange[1]-xrange[0];
    yrange[2] = yrange[1]-yrange[0];
    FT range = stretch * std::max(xrange[2],yrange[2]);
    FT xmid = 0.5*(xrange[0]+xrange[1]);
    FT ymid = 0.5*(yrange[0]+yrange[1]);
    xmin = xmid-0.5*range;
    xmax = xmid+0.5*range;
    ymin = ymid-0.5*range;
    ymax = ymid+0.5*range;
  }

  void resample(const FT d)
  {
    for(iterator c = this->begin();
        c != this->end();
        c++) {
      Component& component = *c;
      component.resample(d);
    }
  }

  void smooth(const unsigned int iter = 1)
  {
    for(iterator c = this->begin();
        c != this->end();
        c++) {
      Component& component = *c;
      component.smooth(iter);
    }
  }

  size_t number_of_points() const
  {
    size_t num_points = 0;
    for(citerator c = this->cbegin();
        c != this->cend();
        c++) {
      const Component& component = *c;
      num_points += component.size();
    }
    return num_points;
  }

  CGAL::Bbox_2 bbox_2() const
  {
    CGAL::Bbox_2 bbox;
    for(citerator c = this->cbegin(); c != this->cend(); ++c) {
      const Component& component = *c;
      bbox += CGAL::bbox_2(component.cbegin(), component.cend());
    }
    return bbox;
  }
};

template <class Kernel, class Triangulation>
Pslg<Kernel> extract_pslg_2_soup(const Triangulation& tr)
{
  using Component = Component<Kernel>;
  using Pslg2 = Pslg<Kernel>;

  using Face_handle = typename Triangulation::Face_handle;

  Pslg2 pslg_soup;
  for (const auto &edge : tr.all_edges())
  {
    if (tr.is_infinite(edge)) continue;
    Face_handle fh2 = edge.first;
    int id = edge.second;
    const Face_handle& neighbor = fh2->neighbor(id);
    if (fh2->is_outside() == neighbor->is_outside()) continue;

    Component new_cmp;
    new_cmp.push_back(fh2->vertex((id+1)%3)->point());
    new_cmp.push_back(fh2->vertex((id+2)%3)->point());
    pslg_soup.push_back(new_cmp);
  }

  return pslg_soup;
}

} // namespace internal
} // namespace Alpha_wraps_2
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_2_DEMO_PSLG_H
