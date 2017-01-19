// Copyright (c) 2015  Università della Svizzera italiana.
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
// $URL$
// $Id$
// 
//
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

#ifndef CGAL_POLYCHAIN_2_H
#define CGAL_POLYCHAIN_2_H

#include <CGAL/license/Segment_Delaunay_graph_Linf_2.h>


#include <CGAL/basic.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/basic.h>


// Polychainsegment_2

namespace CGAL {

namespace Qt {
  template <typename K> class PainterOstream;
}

template <class Traits_P, class Container_P
        = std::vector<typename Traits_P::Point_2> >
class Polychainsegment_2 : public Polygon_2<Traits_P, Container_P> {

  public:
    typedef Traits_P Traits;
    typedef Polygon_2<Traits_P, Container_P> Base;

    // constructors

    Polychainsegment_2() : Base() {}

    Polychainsegment_2(
        const Polychainsegment_2<Traits_P,Container_P>& pc)
      : Base((Base) pc) {}

    template <class InputIterator>
    Polychainsegment_2(InputIterator first, InputIterator last,
                       Traits p_traits = Traits())
        : Base(first, last, p_traits)  {}


    template< class K >
    void draw(CGAL::Qt::PainterOstream<K>& stream) const {
      typedef typename K::Segment_2 K_Segment_2;
      typedef typename Polychainsegment_2<
                        Traits_P,Container_P>::Vertex_const_iterator
              VI;

      if (this->size() > 1) {
        VI source = this->vertices_begin();
        VI target = source+1;
        for( ; target!=this->vertices_end(); ++source, ++target)
        {
          stream << K_Segment_2(*source, *target);
        }
      }
    }

    template< class Stream >
    void draw(Stream & str) const
    {
      typedef typename Traits_P::Segment_2 Segment_2;
      typedef typename Polychainsegment_2<
                        Traits_P,Container_P>::Vertex_const_iterator
              VI;
      if (this->size() > 1) {
        VI source = this->vertices_begin();
        VI target = source+1;
        for( ; target!=this->vertices_end(); ++source, ++target)
        {
          str << Segment_2(*source, *target);
        }
      }
    }

} ;




template <class Traits_P, class Container_P>
std::ostream&
operator<<(std::ostream &os,
           const Polychainsegment_2<Traits_P,Container_P>& p)
{
  typename Polychainsegment_2<Traits_P,Container_P>
              ::Vertex_const_iterator
           i;

  switch(get_mode(os)) {
    case IO::ASCII :
      os << p.size() << ' ';
      for (i = p.vertices_begin(); i != p.vertices_end(); ++i) {
        os << *i << ' ';
      }
      return os;

    case IO::BINARY :
      os << p.size();
      for (i = p.vertices_begin(); i != p.vertices_end(); ++i) {
        os << *i;
      }
      return os;

    default:
      os << "Polychainsegment_2(" << std::endl;
      for (i = p.vertices_begin(); i != p.vertices_end(); ++i) {
        os << "  " << *i << std::endl;
      }
      os << ")" << std::endl;
      return os;
  }
}





/* the following not allowed in some compilers: */
/*
template< class Stream,
          class Traits_P,
          class Container_P = std::vector<typename Traits_P::Point_2>
        >
inline
Stream& operator<<(Stream &s,
                   const Polychainsegment_2 <Traits_P, Container_P> &P)
{
  P.draw(s);
  return s;
}
*/


/* instead, we use overloading */

/*
template< class Stream,
          class Traits_P
        >
inline
Stream& operator<<(Stream &s,
                   const Polychainsegment_2 <Traits_P> &P)
{
  P.draw(s);
  return s;
}


template< class Stream,
          class Traits_P,
          class Container_P
        >
inline
Stream& operator<<(Stream &s,
                   const Polychainsegment_2 <Traits_P, Container_P> &P)
{
  P.draw(s);
  return s;
}
*/



// Polychainray_2

template <class Traits_P, class Container_P
        = std::vector<typename Traits_P::Point_2> >
class Polychainray_2 :
  public Polychainsegment_2<Traits_P, Container_P> {

public:
    typedef Traits_P Traits;

private:
    typedef typename Traits_P::Direction_2 OutgoingDirection;
    typedef Polychainsegment_2<Traits_P, Container_P> Base;

    OutgoingDirection outgoing;

public:

    // constructors

    Polychainray_2(): Base(), outgoing() {}

    Polychainray_2(const Polychainray_2<Traits_P,Container_P>& pcr)
      : Base((Base) pcr), outgoing(pcr.outgoing) {}

    template <class InputIterator>
    Polychainray_2(InputIterator first, InputIterator last,
	           OutgoingDirection d,
              Traits p_traits = Traits())
        : Base(first, last, p_traits), outgoing(d)
    {
    }


    // get_outgoing

    OutgoingDirection get_outgoing() const {
	return outgoing;
    }


    // drawing

    template< class K >
    void draw(CGAL::Qt::PainterOstream<K>& stream) const {
      typedef typename K::Segment_2 K_Segment_2;
      typedef typename K::Ray_2 K_Ray_2;
      typedef typename
	  Polychainray_2<Traits_P,Container_P>::Vertex_const_iterator VI;

      CGAL_assertion( this->size() > 0 );

      VI source = this->vertices_begin();
      if (this->size() > 1) {
        VI target = source+1;
        for( ; target!=this->vertices_end(); ++source, ++target)
        {
          stream << K_Segment_2(*source, *target);
        }
      }
      // now source contains the last point;
      // draw outgoing ray from this point
      stream << K_Ray_2(*source, this->get_outgoing());
    }

    template< class Stream >
    void draw(Stream & stream) const {
      typedef typename Traits_P::Segment_2 Segment_2;
      typedef typename Traits_P::Ray_2 Ray_2;
      typedef typename
	  Polychainray_2<Traits_P,Container_P>::Vertex_const_iterator VI;

      CGAL_assertion( this->size() > 0 );

      VI source = this->vertices_begin();
      if (this->size() > 1) {
        VI target = source+1;
        for( ; target!=this->vertices_end(); ++source, ++target)
        {
          stream << Segment_2(*source, *target);
        }
      }
      // now source contains the last point;
      // draw outgoing ray from this point
      stream << Ray_2(*source, this->get_outgoing());
    }


} ;


template <class Traits_P, class Container_P>
std::ostream&
operator<<(std::ostream &os,
           const Polychainray_2<Traits_P,Container_P>& p)
{
  typename Polychainray_2<Traits_P,Container_P>::Vertex_const_iterator i;

  switch(get_mode(os)) {
    case IO::ASCII :
      os << p.size() << ' ';
      for (i = p.vertices_begin(); i != p.vertices_end(); ++i) {
        os << *i << ' ';
      }
      os << ", dout=" << p.get_outgoing() << ' ';
      return os;

    case IO::BINARY :
      os << p.size();
      for (i = p.vertices_begin(); i != p.vertices_end(); ++i) {
        os << *i;
      }
      os << p.get_outgoing();
      return os;

    default:
      os << "Polychainray_2(" << std::endl;
      for (i = p.vertices_begin(); i != p.vertices_end(); ++i) {
        os << "  " << *i << std::endl;
      }
      os << "dout=" << p.get_outgoing() << std::endl;
      os << ")" << std::endl;
      return os;
  }
}





/*
template< class Stream,
          class Traits_P
        >
inline
Stream& operator<<(Stream &s,
                   const Polychainray_2 <Traits_P> &P)
{
  P.draw(s);
  return s;
}


template< class Stream,
          class Traits_P,
          class Container_P
        >
inline
Stream& operator<<(Stream &s,
                   const Polychainray_2 <Traits_P, Container_P> &P)
{
  P.draw(s);
  return s;
}
*/



// Polychainline_2

template <class Traits_P, class Container_P
        = std::vector<typename Traits_P::Point_2> >
class Polychainline_2 : public Polychainray_2<Traits_P, Container_P> {

public:

private:
    typedef Traits_P Traits;
    typedef typename Traits_P::Direction_2 IncomingDirection;
    typedef typename Traits_P::Direction_2 OutgoingDirection;
    typedef typename Traits_P::Line_2      Line_2;
    typedef typename Traits_P::Ray_2       Ray_2;
    typedef typename Traits_P::Segment_2   Segment_2;
    typedef typename Traits_P::Point_2     Point_2;
    typedef Polychainray_2<Traits_P, Container_P>  Base;
    typedef Polychainline_2<Traits_P, Container_P> Self;

    IncomingDirection incoming;
    bool is_line_optimization;

public:

    // constructors

    Polychainline_2() : Base(), incoming(),
      is_line_optimization(false)
    {}

    Polychainline_2(const Self& pcl)
      : Base((Base) pcl), incoming(pcl.incoming),
        is_line_optimization(pcl.is_line_optimization)
    {}

    template <class InputIterator>
    Polychainline_2(IncomingDirection dinc,
	            InputIterator first, InputIterator last,
	            OutgoingDirection dout,
                    Traits p_traits = Traits())
        : Base(first, last, dout, p_traits), incoming(dinc),
          is_line_optimization(false)
    {
    }

    void set_line_optimization() {
      CGAL_assertion(this->size() == 1);
      CGAL_assertion(this->incoming == - this->get_outgoing());
      is_line_optimization = true;
    }

    // get_incoming

    IncomingDirection get_incoming() const {
	return incoming;
    }


    // unary minus operator
    // that reverses the orientation of the polychainline

    Self
    operator-() const {
      // first reverse list of points of polychainline

      std::vector<typename Traits_P::Point_2> reverse;

      unsigned int npts = this->size();

      CGAL_SDG_DEBUG(std::cout << "pcl_reverse npts="
          << npts << std::endl;);

      reverse.resize(npts);

      //std::reverse_copy(this->vertices_begin(),
      //                  this->vertices_end(),
      //                  reverse);

      typedef typename Polychainline_2<
                        Traits_P,Container_P>::Vertex_const_iterator
              VI;

      unsigned int i = npts - 1;

      for (VI vi = this->vertices_begin();
              vi != this->vertices_end();
              ++vi)
      {
        //CGAL_SDG_DEBUG(std::cout << "debug pcl_reverse setting at " << i
        //          << " value " << *vi << std::endl;);
        reverse[i] = *vi;
        --i;
      }


      // then also swap incoming and outgoing directions

      //return Self(this->get_outgoing(),
      //            reverse.begin(),
      //            reverse.end(),
      //            this->get_incoming());

      Self pclreverse(this->get_outgoing(),
                      reverse.begin(),
                      reverse.end(),
                      this->get_incoming());

      if (is_line_optimization) {
        pclreverse.set_line_optimization();
      }

      return pclreverse;
    }

    // line_first_intersection_point_with
    typename Traits_P::Point_2
    line_first_intersection_point_with(
        const Polychainline_2<Traits_P,Container_P>& pcl)
    {
      CGAL_assertion(is_line_optimization);
      Line_2 line((*this)[0], this->get_outgoing());

      typedef typename
	  Polychainline_2<Traits_P,Container_P>::
                     Vertex_const_iterator
          VI;

      CGAL::Object result;

      VI sourcepcl  = pcl.vertices_begin();
      Ray_2 rayincpcl (*sourcepcl,  pcl.get_incoming());
      result = CGAL::intersection(line, rayincpcl);
      if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
        return *ipoint;
      }

      if (pcl.size() > 1) {
        VI targetpcl = sourcepcl+1;
        for( ;
            targetpcl != pcl.vertices_end();
            ++sourcepcl, ++targetpcl)
        {
          Segment_2 testseg(*sourcepcl, *targetpcl);
          result = CGAL::intersection(line, testseg);
          if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
            return *ipoint;
          }
        }
      }
      Ray_2 rayoutpcl(*sourcepcl, pcl.get_outgoing());
      result = CGAL::intersection(line, rayoutpcl);
      if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
        return *ipoint;
      }

      CGAL_SDG_DEBUG(std::cout
          << "debug error: no intersection found for "
          << "this=" << *this << " pcl=" << pcl << std::endl;);

      CGAL_assertion(false);

      return Point_2();

    } // end of line_first_intersection_point_with

    // first_intersection_point_with
    typename Traits_P::Point_2
    first_intersection_point_with(
        const Polychainline_2<Traits_P,Container_P>& pcl)
    {
      // for every piece of object,
      // try intersecting with every piece of pcl

      CGAL_SDG_DEBUG(std::cout
          << "debug first_intersection entering this="
          << *this << " pcl=" << pcl << std::endl;);

      typedef typename
	  Polychainline_2<Traits_P,Container_P>::
                     Vertex_const_iterator
          VI;

      typedef typename std::vector<typename Traits_P::Segment_2>::
               const_iterator SI;

      CGAL_assertion( this->size() > 0 );
      CGAL_assertion( pcl.size() > 0 );

      if (this->is_line_optimization) {
        return line_first_intersection_point_with(pcl);
      }

#if 0
      CGAL_SDG_DEBUG(std::cout << "debug first_intersection "
          << "creating empty vectors" << std::endl;);
#endif

      // create two empty vectors for storing the segments
      std::vector<Segment_2> segmentsthis;
      std::vector<Segment_2> segmentspcl;

      VI sourcethis = this->vertices_begin();
      Ray_2 rayincthis(*sourcethis, this->get_incoming());
      if (this->size() > 1) {
        VI targetthis = sourcethis+1;
        for( ;
            targetthis != this->vertices_end();
            ++sourcethis, ++targetthis)
        {
          segmentsthis.push_back(Segment_2(*sourcethis, *targetthis));
        }
      }
      Ray_2 rayoutthis(*sourcethis, this->get_outgoing());

      CGAL_assertion(this->size() == (segmentsthis.size() + 1));

#if 0
      CGAL_SDG_DEBUG(std::cout << "debug first_intersection "
                     << "segmentsthis computed" << std::endl;);
#endif


      VI sourcepcl  = pcl.vertices_begin();
      Ray_2 rayincpcl (*sourcepcl,  pcl.get_incoming());
      if (pcl.size() > 1) {
        VI targetpcl = sourcepcl+1;
        for( ;
            targetpcl != pcl.vertices_end();
            ++sourcepcl, ++targetpcl)
        {
          segmentspcl.push_back(Segment_2(*sourcepcl, *targetpcl));
        }
      }
      Ray_2 rayoutpcl(*sourcepcl, pcl.get_outgoing());

      CGAL_assertion(pcl.size() == (segmentspcl.size() + 1));

#if 0
      CGAL_SDG_DEBUG(std::cout << "debug first_intersection "
                     << "segmentspcl computed" << std::endl;);
#endif


      CGAL::Object result;

#if 0
      CGAL_SDG_DEBUG(std::cout << "debug first_intersection "
                    << "trying rayincthis with pcl" << std::endl;);
#endif

#if 0
      CGAL_SDG_DEBUG(std::cout << "debug first_intersection "
                << "trying ray " << rayincthis
                     << " with ray " << rayincpcl << std::endl;);
#endif

      result = CGAL::intersection(rayincthis, rayincpcl);
      if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
        return *ipoint;
      }

      for (SI sipcl = segmentspcl.begin();
              sipcl != segmentspcl.end();
              ++sipcl) {
#if 0
        CGAL_SDG_DEBUG(std::cout << "debug first_intersection "
                << "trying ray " << rayincthis
                       << " with segment " << *sipcl << std::endl;);
#endif

        result = CGAL::intersection(rayincthis, *sipcl);
        if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
          return *ipoint;
        }
      }

#if 0
      CGAL_SDG_DEBUG(std::cout << "debug first_intersection "
                << "trying ray " << rayincthis
                     << " with ray " << rayoutpcl << std::endl;);
#endif

      result = CGAL::intersection(rayincthis, rayoutpcl);
      if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
        return *ipoint;
      }

#if 0
      CGAL_SDG_DEBUG(std::cout << "debug first_intersection "
                     << "trying segmentsthis with pcl" << std::endl;);
#endif

      for (SI sithis = segmentsthis.begin();
              sithis != segmentsthis.end();
              ++sithis) {

#if 0
        CGAL_SDG_DEBUG(std::cout << "debug first_intersection "
                << "trying segment " << *sithis
                       << " with ray " << rayincpcl << std::endl;);
#endif

        result = CGAL::intersection(*sithis, rayincpcl);
        if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
          return *ipoint;
        }

        for (SI sipcl = segmentspcl.begin();
            sipcl != segmentspcl.end();
            ++sipcl) {
#if 0
          CGAL_SDG_DEBUG(std::cout << "debug first_intersection "
                << "trying segment " << *sithis
                         << " with segment " << *sipcl << std::endl;);
#endif

          result = CGAL::intersection(*sithis, *sipcl);
          if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
          {
            return *ipoint;
          }
        }

#if 0
        CGAL_SDG_DEBUG(std::cout << "debug first_intersection "
            << "trying segment " << *sithis
            << " with ray " << rayoutpcl << std::endl;);
#endif

        result = CGAL::intersection(*sithis, rayoutpcl);
        if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
          return *ipoint;
        }

      }

#if 0
      CGAL_SDG_DEBUG(std::cout << "debug first_intersection "
                     << "trying rayoutthis with pcl" << std::endl;);
#endif

#if 0
      CGAL_SDG_DEBUG(std::cout << "debug first_intersection "
                << "trying ray " << rayoutthis
                     << " with ray " << rayincpcl << std::endl;);
#endif

      result = CGAL::intersection(rayoutthis, rayincpcl);
      if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
        return *ipoint;
      }

      for (SI sipcl = segmentspcl.begin();
              sipcl != segmentspcl.end();
              ++sipcl) {
#if 0
        CGAL_SDG_DEBUG(std::cout << "debug first_intersection "
                << "trying ray " << rayoutthis
                       << " with segment " << *sipcl << std::endl;);
#endif

        result = CGAL::intersection(rayoutthis, *sipcl);
        if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
          return *ipoint;
        }
      }

#if 0
      CGAL_SDG_DEBUG(std::cout << "debug first_intersection "
                << "trying ray " << rayoutthis
                     << " with ray " << rayoutpcl << std::endl;);
#endif

      result = CGAL::intersection(rayoutthis, rayoutpcl);
      if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
        return *ipoint;
      }

      CGAL_SDG_DEBUG(std::cout
          << "debug error: no intersection found for "
          << "this=" << *this << " pcl=" << pcl << std::endl;);

      CGAL_assertion(false);

      // philaris: added to avoid complaining from some compilers

      return Point_2();

    }


    // drawing

    template< class K >
    void draw(CGAL::Qt::PainterOstream<K>& stream) const {
      typedef typename K::Segment_2 K_Segment_2;
      typedef typename K::Ray_2 K_Ray_2;
      typedef typename
	  Polychainline_2<Traits_P,Container_P>::Vertex_const_iterator VI;

      CGAL_assertion( this->size() > 0 );

      VI source = this->vertices_begin();

      // draw outgoing ray from first point
      stream << K_Ray_2(*source, this->get_incoming());

      if (this->size() > 1) {
        VI target = source+1;
        for( ; target!=this->vertices_end(); ++source, ++target)
        {
          stream << K_Segment_2(*source, *target);
        }
      }

      // now source contains the last point;
      // draw outgoing ray from this point
      stream << K_Ray_2(*source, this->get_outgoing());
    }

    template< class Stream >
    void draw(Stream & stream) const {
      typedef typename
	  Polychainline_2<Traits_P,Container_P>::Vertex_const_iterator VI;

      CGAL_assertion( this->size() > 0 );

      VI source = this->vertices_begin();

      // draw outgoing ray from first point
      stream << Ray_2(*source, this->get_incoming());

      if (this->size() > 1) {
        VI target = source+1;
        for( ; target!=this->vertices_end(); ++source, ++target)
        {
          stream << Segment_2(*source, *target);
        }
      }

      // now source contains the last point;
      // draw outgoing ray from this point
      stream << Ray_2(*source, this->get_outgoing());
    }

} ;



template <class Traits_P, class Container_P>
std::ostream&
operator<<(std::ostream &os,
           const Polychainline_2<Traits_P,Container_P>& p)
{
  typename Polychainline_2<Traits_P,Container_P>::Vertex_const_iterator i;

  switch(get_mode(os)) {
    case IO::ASCII :
      os << p.size() << ' ';
      os << ", dinc=" << p.get_incoming() << ", ";
      for (i = p.vertices_begin(); i != p.vertices_end(); ++i) {
        os << *i << ' ';
      }
      os << ", dout=" << p.get_outgoing() << ' ';
      return os;

    case IO::BINARY :
      os << p.size();
      os << p.get_incoming();
      for (i = p.vertices_begin(); i != p.vertices_end(); ++i) {
        os << *i;
      }
      os << p.get_outgoing();
      return os;

    default:
      os << "Polychainline_2(" << std::endl;
      os << "dinc=" << p.get_incoming() << std::endl;
      for (i = p.vertices_begin(); i != p.vertices_end(); ++i) {
        os << "  " << *i << std::endl;
      }
      os << "dout=" << p.get_outgoing() << std::endl;
      os << ")" << std::endl;
      return os;
  }
}


/*
template< class Stream,
          class Traits_P
        >
inline
Stream& operator<<(Stream &s,
                   const Polychainline_2 <Traits_P> &P)
{
  P.draw(s);
  return s;
}


template< class Stream,
          class Traits_P,
          class Container_P
        >
inline
Stream& operator<<(Stream &s,
                   const Polychainline_2 <Traits_P, Container_P> &P)
{
  P.draw(s);
  return s;
}
*/



} //namespace CGAL

#endif // CGAL_POLYCHAIN_2_H
