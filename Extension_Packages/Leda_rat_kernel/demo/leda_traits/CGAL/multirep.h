// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/multirep.h
// package       : 
// maintainer    : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : 
// coordinator   : 
//
// ======================================================================


#ifndef CGAL_MULTIREP_H
#define CGAL_MULTIREP_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Number_type_traits.h>
#include <iostream>


CGAL_BEGIN_NAMESPACE


// templated by first (FR) and second (SR) representation
// by a Converter (Converter) ...
// ... and by a tag class used for tagging lokal NT type definitions

template<class T = Tag_false>
struct MultiRepNTSupport {

  template<class ENT>
  struct Inner {
    typedef Tag_false   Has_gcd;
    typedef Tag_false   Has_division;
    typedef Tag_false   Has_sqrt;    
  };
   
};


template<>
struct MultiRepNTSupport<Tag_true> {
  
  template<class ENT>
  struct Inner {
    typedef Number_type_traits<ENT>     NTT;
    
    typedef typename NTT::Has_gcd       Has_gcd;
    typedef typename NTT::Has_division  Has_division;
    typedef typename NTT::Has_sqrt      Has_sqrt;         
  };  
  
};



template<class FR, class SR, class Converter,
         class NTAG = CGAL::Tag_false>
struct MultiRep : public FR {
protected:

  // second representation ...
  SR second_rep;
  
public:
  typedef typename MultiRepNTSupport<NTAG>::Inner<FR>  HELPER;
  
  typedef typename HELPER::Has_gcd                     Has_gcd;
  typedef typename HELPER::Has_division                Has_division;
  typedef typename HELPER::Has_sqrt                    Has_sqrt;   


  typedef  FR          first_representation;
  typedef  SR          second_representation;
  typedef  Converter   converter;

  // we need a special constructor for two precomputed
  // representations (maybe we have to make this a member function because of VC ?)
  
  MultiRep(const FR& rep1, const SR& rep2) : FR(rep1), second_rep(rep2)
  { }

  // now the templated constructors ...

  MultiRep() : FR()
  { second_rep = Converter()(*((FR*)this)); }

  template<class T1>
  MultiRep(T1 arg1) : FR(arg1)
  { second_rep = Converter()(*((FR*)this)); }
  
  template<class T1,class T2>
  MultiRep(T1 arg1, T2 arg2) : FR(arg1,arg2)
  { second_rep = Converter()(*((FR*)this)); }
  
  template<class T1,class T2,class T3>
  MultiRep(T1 arg1, T2 arg2, T3 arg3) : FR(arg1,arg2,arg3)
  { second_rep = Converter()(*((FR*)this)); }
  
  template<class T1,class T2,class T3,class T4>
  MultiRep(T1 arg1, T2 arg2, T3 arg3, T4 arg4) : FR(arg1,arg2,arg3,arg4)
  { second_rep = Converter()(*((FR*)this)); }
 
  template<class T1,class T2,class T3,class T4,class T5>
  MultiRep(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5) : FR(arg1,arg2,arg3,arg4,arg5)
  { second_rep = Converter()(*((FR*)this)); }
  
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  MultiRep(T1 arg1, T2 arg2, T3 arg3, T4 arg4, 
           T5 arg5, T6 arg6) : FR(arg1,arg2,arg3,arg4,arg5,arg6)
  { second_rep = Converter()(*((FR*)this)); }
  
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  MultiRep(T1 arg1, T2 arg2, T3 arg3, T4 arg4, 
           T5 arg5, T6 arg6, T7 arg7) : FR(arg1,arg2,arg3,arg4,arg5,arg6,arg7)
  { second_rep = Converter()(*((FR*)this)); }
  
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
  MultiRep(T1 arg1, T2 arg2, T3 arg3, T4 arg4, 
           T5 arg5, T6 arg6, T7 arg7, T8 arg8) : FR(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
  { second_rep = Converter()(*((FR*)this)); }       
 
  
  // assignment ...
  
  MultiRep& operator=(const FR& val)
  {
    (*((FR*)this)) = val;
    second_rep = Converter()(*((FR*)this)); 
    return *this;
  }  
  
  const FR&   get_first_rep() const { return *this; }
  const SR&   get_second_rep() const { return second_rep; }
  
  void output(std::ostream& o) const
  {
    o << get_first_rep() << "  " << get_second_rep() << " ";
  }
};

template<class T,class CONV,class TAG>
inline Interval_base to_interval (const MultiRep<T,CGAL::Interval_base,CONV,TAG>& z)
{ 
  return z.get_second_rep(); 
}

CGAL_END_NAMESPACE

#endif
