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
// file          : Triangulation/include/CGAL/Constrained_triangulation_2.h
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

#include <pair.h>
#include <list.h>
#include <vector.h>
#include <map.h> 

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Constrained_triangulation_sweep_2.h>


template < class Gt, class Tds>
class CGAL_Constrained_triangulation_2
  : public CGAL_Triangulation_2<Gt,Tds>
{

public:
  typedef CGAL_Triangulation_2<Gt,Tds> Triangulation;
  typedef CGAL_Constrained_triangulation_2<Gt,Tds> Constrained_triangulation;
  typedef pair<Point,Point> Constraint;

  typedef CGAL_Constrained_triangulation_sweep_2<Gt,Tds>  Sweep;

  CGAL_Constrained_triangulation_2() : Triangulation() { }

  CGAL_Constrained_triangulation_2(const Gt& gt) : Triangulation(gt) { }

  CGAL_Constrained_triangulation_2(const CGAL_Constrained_triangulation_2& ct)
    : Triangulation(ct) {}

  CGAL_Constrained_triangulation_2(const Vertex_handle&  v, const Gt& gt) 
    : Triangulation(v,gt) {}

  CGAL_Constrained_triangulation_2(list<Constraint>& lc, const Gt& gt=Gt())
      : CGAL_Triangulation_2<Gt,Tds>(gt)
  {
    Sweep sweep(lc,gt);
    CGAL_Triangulation_2<Gt,Tds> Tr ( sweep.vertex(), gt);
    swap(Tr);
  }

 #ifdef CGAL_CFG_NO_MEMBER_TEMPLATES
  
  #if defined(LIST_H) || defined(__SGI_STL_LIST_H)
  CGAL_Constrained_triangulation_2(list<Constraint>::const_iterator first,
                                   list<Constraint>::const_iterator last,
                                   const Gt& gt=Gt() )
     : CGAL_Triangulation_2<Gt,Tds>(gt)
  {
      list<Constraint> lc;
      while(first != last){
          lc.push_back(*first); ++first;
      }
      Sweep sweep(lc,gt);
      CGAL_Triangulation_2<Gt,Tds> Tr ( sweep.vertex(), gt);
      swap(Tr);
      //init(sweep.vertex());
      CGAL_triangulation_postcondition( is_valid() );
  }
  #endif // LIST_H
  
  #if defined(VECTOR_H) || defined(__SGI_STL_VECTOR_H)
  CGAL_Constrained_triangulation_2(vector<Constraint>::const_iterator first,
                                   vector<Constraint>::const_iterator last,
                                    const Gt& gt=Gt() )
     : CGAL_Triangulation_2<Gt,Tds>(gt)
  {
      list<Constraint> lc;
      while(first != last){
          lc.push_back(*first); ++first;
      }
      Sweep sweep(lc,gt);
      CGAL_Triangulation_2<Gt,Tds> Tr ( sweep.vertex(), gt);
      swap(Tr);
    //init(sweep.vertex());
      CGAL_triangulation_postcondition( is_valid() );
  }
  #endif // VECTOR_H
  
  #ifdef ITERATOR_H
  CGAL_Constrained_triangulation_2(istream_iterator<Constraint,ptrdiff_t> first,
                                   istream_iterator<Constraint,ptrdiff_t> last,
                                    const Gt& gt=Gt() )
     : CGAL_Triangulation_2<Gt,Tds>(gt)
  {
      list<Constraint> lc;
      while(first != last){
          lc.push_back(*first); ++first;
      }
      Sweep sweep(lc,gt);
      CGAL_Triangulation_2<Gt,Tds> Tr ( sweep.vertex(), gt);
      swap(Tr);
      //init(sweep.vertex());
      CGAL_triangulation_postcondition( is_valid() );
  }
  #endif // ITERATOR_H
  
  CGAL_Constrained_triangulation_2(Constraint* first,
                                   Constraint* last,
                                    const Gt& gt=Gt() )
     : CGAL_Triangulation_2<Gt,Tds>(gt)
  {
      list<Constraint> lc;
      while(first != last){
          lc.push_back(*first); ++first;
      }
      Sweep sweep(lc,gt);
      CGAL_Triangulation_2<Gt,Tds> Tr ( sweep.vertex(), gt);
      swap(Tr);
      //init(sweep.vertex());
      CGAL_triangulation_postcondition( is_valid() );
  }
  
  #else
  
  template<class InputIterator>
  CGAL_Constrained_triangulation_2(InputIterator first,
                                   InputIterator last,
                                   const Gt& gt=Gt() )
     : CGAL_Triangulation_2<Gt,Tds>(gt)
  {
      list<Constraint> lc;
      while(first != last){
          lc.push_back(*first++);
      }
      Sweep sweep(lc,gt);
      CGAL_Triangulation_2<Gt,Tds> Tr ( sweep.vertex(), gt);
      swap(Tr);
          //init(sweep.vertex());
      CGAL_triangulation_postcondition( is_valid() );
  }
  
  #endif // CGAL_CFG_NO_MEMBER_TEMPLATES
  
  // private:
  // private part of class CGAL_Constrained_triangulation_2
};

template < class Gt, class Tds >
ostream &
operator<<(ostream& os, const CGAL_Constrained_triangulation_2<Gt,Tds> &Ct)
{
  return os << (const CGAL_Triangulation_2<Gt,Tds>&)Ct;
}



#endif CGAL_CONSTRAINED_TRIANGULATION_2_H
