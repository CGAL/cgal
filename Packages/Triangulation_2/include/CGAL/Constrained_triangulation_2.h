// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Constrained_triangulation_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$

// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  < Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================


#ifndef CGAL_CONSTRAINED_TRIANGULATION_2_H
#define CGAL_CONSTRAINED_TRIANGULATION_2_H

#include <utility>
#include <list>
#include <vector>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Constrained_triangulation_sweep_2.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds>
class Constrained_triangulation_2
  : public Triangulation_2<Gt,Tds>
{

public:
  typedef Triangulation_2<Gt,Tds> Triangulation;
  typedef Constrained_triangulation_2<Gt,Tds> Constrained_triangulation;
  typedef std::pair<Point,Point> Constraint;

  typedef Constrained_triangulation_sweep_2<Gt,Tds>  Sweep;

  Constrained_triangulation_2() : Triangulation() { }

  Constrained_triangulation_2(const Gt& gt) : Triangulation(gt) { }

  Constrained_triangulation_2(const Constrained_triangulation_2& ct)
    : Triangulation(ct) {}

  Constrained_triangulation_2(const Vertex_handle&  v, const Gt& gt) 
    : Triangulation(v,gt) {}

  Constrained_triangulation_2(std::list<Constraint>& lc, const Gt& gt=Gt())
      : Triangulation_2<Gt,Tds>(gt)
  {
    Sweep sweep(lc,gt);
    Triangulation_2<Gt,Tds> Tr ( sweep.vertex(), gt);
    swap(Tr);
  }

 #ifdef CGAL_CFG_NO_MEMBER_TEMPLATES
  
  #if defined(LIST_H) || defined(__SGI_STL_LIST_H)
  Constrained_triangulation_2(std::list<Constraint>::const_iterator first,
                                   std::list<Constraint>::const_iterator last,
                                   const Gt& gt=Gt() )
     : Triangulation_2<Gt,Tds>(gt)
  {
      std::list<Constraint> lc;
      while(first != last){
          lc.push_back(*first); ++first;
      }
      Sweep sweep(lc,gt);
      Triangulation_2<Gt,Tds> Tr ( sweep.vertex(), gt);
      swap(Tr);
      //init(sweep.vertex());
      CGAL_triangulation_postcondition( is_valid() );
  }
  #endif // LIST_H
  
  #if defined(VECTOR_H) || defined(__SGI_STL_VECTOR_H)
  Constrained_triangulation_2(std::vector<Constraint>::const_iterator first,
                              std::vector<Constraint>::const_iterator last,
			      const Gt& gt=Gt() )
     : Triangulation_2<Gt,Tds>(gt)
  {
      std::list<Constraint> lc;
      while(first != last){
          lc.push_back(*first); ++first;
      }
      Sweep sweep(lc,gt);
      Triangulation_2<Gt,Tds> Tr ( sweep.vertex(), gt);
      swap(Tr);
    //init(sweep.vertex());
      CGAL_triangulation_postcondition( is_valid() );
  }
  #endif // VECTOR_H
  
  #ifdef ITERATOR_H
  Constrained_triangulation_2(
			   std::istream_iterator<Constraint,ptrdiff_t> first,
			   std::istream_iterator<Constraint,ptrdiff_t> last,
			   const Gt& gt=Gt() )
     : Triangulation_2<Gt,Tds>(gt)
  {
      std::list<Constraint> lc;
      while(first != last){
          lc.push_back(*first); ++first;
      }
      Sweep sweep(lc,gt);
      Triangulation_2<Gt,Tds> Tr ( sweep.vertex(), gt);
      swap(Tr);
      //init(sweep.vertex());
      CGAL_triangulation_postcondition( is_valid() );
  }
  #endif // ITERATOR_H
  
  Constrained_triangulation_2(Constraint* first,
                                   Constraint* last,
                                    const Gt& gt=Gt() )
     : Triangulation_2<Gt,Tds>(gt)
  {
      std::list<Constraint> lc;
      while(first != last){
          lc.push_back(*first); ++first;
      }
      Sweep sweep(lc,gt);
      Triangulation_2<Gt,Tds> Tr ( sweep.vertex(), gt);
      swap(Tr);
      //init(sweep.vertex());
      CGAL_triangulation_postcondition( is_valid() );
  }
  
  #else
  
  template<class InputIterator>
  Constrained_triangulation_2(InputIterator first,
                                   InputIterator last,
                                   const Gt& gt=Gt() )
     : Triangulation_2<Gt,Tds>(gt)
  {
      std::list<Constraint> lc;
      while(first != last){
          lc.push_back(*first++);
      }
      Sweep sweep(lc,gt);
      Triangulation_2<Gt,Tds> Tr ( sweep.vertex(), gt);
      swap(Tr);
          //init(sweep.vertex());
      CGAL_triangulation_postcondition( is_valid() );
  }
  
  #endif // CGAL_CFG_NO_MEMBER_TEMPLATES
  
  // private:
  // private part of class Constrained_triangulation_2
};

template < class Gt, class Tds >
std::ostream &
operator<<(std::ostream& os, const Constrained_triangulation_2<Gt,Tds> &Ct)
{
  return os << (const Triangulation_2<Gt,Tds>&)Ct;
}

CGAL_END_NAMESPACE

#endif //CGAL_CONSTRAINED_TRIANGULATION_2_H
