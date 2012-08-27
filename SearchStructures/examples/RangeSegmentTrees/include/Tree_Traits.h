// Copyright (c) 1997  ETH Zurich (Switzerland).
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
// Author(s)     : Gabriele Neyer

#ifndef __TREE_TRAITS__
#define __TREE_TRAITS__

#include <functional>
#include <utility>

namespace CGAL {

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

} //namespace CGAL

#endif

