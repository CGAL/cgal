#ifndef __TREE_TRAITS__
#define __TREE_TRAITS__
// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of the
// Computational Geometry Algorithms Library (CGAL).
//
// Every use of CGAL requires a license. Licenses come in three kinds:
//
// - For academic research and teaching purposes, permission to use and
//   copy the software and its documentation is hereby granted free of  
//   charge, provided that
//   (1) it is not a component of a commercial product, and
//   (2) this notice appears in all copies of the software and
//       related documentation.
// - Development licenses grant access to the source code of the library 
//   to develop programs. These programs may be sold to other parties as 
//   executable code. To obtain a development license, please contact
//   the CGAL Consortium (at cgal@cs.uu.nl).
// - Commercialization licenses grant access to the source code and the
//   right to sell development licenses. To obtain a commercialization 
//   license, please contact the CGAL Consortium (at cgal@cs.uu.nl).
//
// This software and documentation is provided "as-is" and without
// warranty of any kind. In no event shall the CGAL Consortium be
// liable for any damage of any kind.
//
// The CGAL Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Free University of Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany) Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// release       : CGAL-1.2
// release_date  : 1999, January 18
//
// file          : src/test/RangeSegmentTrees/include/Tree_Traits.h
// source        : src/test/RangeSegmentTrees/include/Tree_Traits.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Gabriele Neyer
//
// coordinator   : Peter Widmayer, ETH Zurich
//
//
// email         : cgal@cs.uu.nl
//
// ======================================================================

#include <functional>
#include <utility>

CGAL_BEGIN_NAMESPACE

class Tree_traits_1{
 public:
  typedef double Key;
  typedef double Key_1;
  typedef std::pair<Key,Key> Interval;

  class C_Low_1{
  public:
    Key_1 operator()(const Interval& i)
    { return i.first;}
  };

  class C_High_1{
  public:
    Key_1 operator()(const Interval& i)
    { return i.second;}
  };

  class C_Key_1{
  public:
    Key_1 operator()(const Key& k)
    { return k;}
  };

  class C_Compare_1{
  public:
    bool operator()(Key_1 k1, Key_1 k2)
    {
      return std::less<double>()(k1,k2);
    }
  };

  typedef C_Compare_1 compare_1;
  typedef C_Low_1 low_1;
  typedef C_High_1 high_1;
  typedef C_Key_1 key_1;
};



class Tree_traits_2{
 public:
  typedef std::pair<int, double> Key;
  typedef int Key_1;
  typedef double Key_2;
  typedef std::pair<Key,Key> Interval;

  class C_Low_1{
  public:
    Key_1 operator()(const Interval& i)
    { return i.first.first;}
  };

  class C_High_1{
  public:
    Key_1 operator()(const Interval& i)
    { return i.second.first;}
  };

  class C_Low_2{
  public:
    Key_2 operator()(const Interval& i)
    { return i.first.second;}
  };

  class C_High_2{
  public:
    Key_2 operator()(const Interval& i)
    { return i.second.second;}
  };

  class C_Key_1{
  public:
    Key_1 operator()(const Key& k)
    { return k.first;}
  };

  class C_Key_2{
  public:
    Key_2 operator()(const Key& k)
    { return k.second;}
  };
  
  class C_Compare_1{
  public:
    bool operator()(Key_1 k1, Key_1 k2)
    {
      return std::less<int>()(k1,k2);
    }
  };

  class C_Compare_2{
  public:
    
    bool operator()(Key_2 k1, Key_2 k2)
    {
      return std::less<double>()(k1,k2);
    }
  };

  typedef C_Compare_1 compare_1;
  typedef C_Compare_2 compare_2;
  typedef C_Low_1 low_1;
  typedef C_High_1 high_1;
  typedef C_Key_1 key_1;
  typedef C_Low_2 low_2;
  typedef C_High_2 high_2;
  typedef C_Key_2 key_2;

};

class tuple_3{
 public:
  int key_1;
  double key_2;
  long key_3;
  tuple_3(){}
  tuple_3(int a, double b, long c) : key_1(a), key_2(b), key_3(c) {}
};

class Tree_traits_3{
 public:
  typedef tuple_3 Key;
  typedef int Key_1;
  typedef double Key_2;
  typedef long Key_3;
  typedef std::pair<Key,Key> Interval;

  class C_Low_1{
  public:
    Key_1 operator()(const Interval& i)
    { return i.first.key_1;}
  };

  class C_High_1{
  public:
    Key_1 operator()(const Interval& i)
    { return i.second.key_1;}
  };

  class C_Low_2{
  public:
    Key_2 operator()(const Interval& i)
    { return i.first.key_2;}
  };

  class C_High_2{
  public:
    Key_2 operator()(const Interval& i)
    { return i.second.key_2;}
  };

  class C_Low_3{
  public:
    Key_3 operator()(const Interval& i)
    { return i.first.key_3;}
  };

  class C_High_3{
  public:
    Key_3 operator()(const Interval& i)
    { return i.second.key_3;}
  };

  class C_Key_1{
  public:
    Key_1 operator()(const Key& k)
    { return k.key_1;}
  };

  class C_Key_2{
  public:
    Key_2 operator()(const Key& k)
    { return k.key_2;}
  };

  class C_Key_3{
  public:
    Key_3 operator()(const Key& k)
    { return k.key_3;}
  };
  
  class C_Compare_1{
  public:
    
    bool operator()(Key_1 k1, Key_1 k2)
    {
      return std::less<int>()(k1,k2);
    }
  };

  class C_Compare_2{
  public:
    
    bool operator()(Key_2 k1, Key_2 k2)
    {
      return std::less<double>()(k1,k2);
    }
  };

  class C_Compare_3{
  public:
    
    bool operator()(Key_3 k1, Key_3 k2)
    {
      return std::less<long>()(k1,k2);
    }
  };

  typedef C_Compare_1 compare_1;
  typedef C_Compare_2 compare_2;
  typedef C_Compare_3 compare_3;
  typedef C_Low_1 low_1;
  typedef C_Low_2 low_2;
  typedef C_Low_3 low_3;
  typedef C_High_1 high_1;
  typedef C_High_2 high_2;
  typedef C_High_3 high_3;
  typedef C_Key_1 key_1;
  typedef C_Key_2 key_2;
  typedef C_Key_3 key_3;
};



class tuple_4{
 public:
  int key_1;
  double key_2;
  long key_3;
  double key_4;
  tuple_4(){}
  tuple_4(int a, double b, long c, double d) : key_1(a), key_2(b), key_3(c), key_4(d) {}
};

class Tree_traits_4{
 public:
  typedef tuple_4 Key;
  typedef int Key_1;
  typedef double Key_2;
  typedef long Key_3;
  typedef double Key_4;
  typedef std::pair<Key,Key> Interval;

  class C_Low_1{
  public:
    Key_1 operator()(const Interval& i)
    { return i.first.key_1;}
  };

  class C_High_1{
  public:
    Key_1 operator()(const Interval& i)
    { return i.second.key_1;}
  };

  class C_Low_2{
  public:
    Key_2 operator()(const Interval& i)
    { return i.first.key_2;}
  };

  class C_High_2{
  public:
    Key_2 operator()(const Interval& i)
    { return i.second.key_2;}
  };

  class C_Low_3{
  public:
    Key_3 operator()(const Interval& i)
    { return i.first.key_3;}
  };

  class C_High_3{
  public:
    Key_3 operator()(const Interval& i)
    { return i.second.key_3;}
  };

  class C_Low_4{
  public:
    Key_4 operator()(const Interval& i)
    { return i.first.key_4;}
  };

  class C_High_4{
  public:
    Key_4 operator()(const Interval& i)
    { return i.second.key_4;}
  };

  class C_Key_1{
  public:
    Key_1 operator()(const Key& k)
    { return k.key_1;}
  };

  class C_Key_2{
  public:
    Key_2 operator()(const Key& k)
    { return k.key_2;}
  };

  class C_Key_3{
  public:
    Key_3 operator()(const Key& k)
    { return k.key_3;}
  };

  class C_Key_4{
  public:
    Key_4 operator()(const Key& k)
    { return k.key_4;}
  };
  
  class C_Compare_1{
  public:
    
    bool operator()(Key_1 k1, Key_1 k2)
    {
      return std::less<int>()(k1,k2);
    }
  };

  class C_Compare_2{
  public:
    
    bool operator()(Key_2 k1, Key_2 k2)
    {
      return std::less<double>()(k1,k2);
    }
  };

  class C_Compare_3{
  public:
    
    bool operator()(Key_3 k1, Key_3 k2)
    {
      return std::less<long>()(k1,k2);
    }
  };

  class C_Compare_4{
  public:
    
    bool operator()(Key_4 k1, Key_4 k2)
    {
      return std::less<double>()(k1,k2);
    }
  };

  typedef C_Compare_1 compare_1;
  typedef C_Compare_2 compare_2;
  typedef C_Compare_3 compare_3;
  typedef C_Compare_4 compare_4;
  typedef C_Low_1 low_1;
  typedef C_Low_2 low_2;
  typedef C_Low_3 low_3;
  typedef C_Low_4 low_4;
  typedef C_High_1 high_1;
  typedef C_High_2 high_2;
  typedef C_High_3 high_3;
  typedef C_High_4 high_4;
  typedef C_Key_1 key_1;
  typedef C_Key_2 key_2;
  typedef C_Key_3 key_3;
  typedef C_Key_4 key_4;
};

CGAL_END_NAMESPACE

#endif

