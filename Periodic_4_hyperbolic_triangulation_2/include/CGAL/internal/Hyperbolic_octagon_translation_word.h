// Copyright (c) 1999-2004,2006-2009,2014-2016   INRIA Nancy - Grand Est (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Iordan Iordanov <iordan.iordanov@loria.fr>

#ifndef CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_WORD
#define CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_WORD

#include <CGAL/license/Periodic_4_hyperbolic_triangulation_2.h>

#include <CGAL/internal/Hyperbolic_octagon_translation_matrix.h>
#include <CGAL/internal/Dehn_hyperbolic_octagon_translation_word.h>

#include <iostream>
#include <vector>

//-------------------------------------------------------
// Global variables -- used only for profiling!
#if defined PROFILING_MODE

extern long calls_apply_identity;
extern long calls_apply_non_identity;
extern long calls_append_identity;
extern long calls_append_non_identity;

#endif
//-------------------------------------------------------

namespace CGAL {

template <class Int>
class Hyperbolic_octagon_translation_word
{
  typedef Dehn_hyperbolic_octagon_translation_word<Int>   Dehn_reductor;
  typedef Hyperbolic_octagon_translation_word<Int>        Self;

private:
  Int   w0 : 3;
  Int   w1 : 3;
  Int   w2 : 3;
  Int   w3 : 3;
  bool  b0 : 1;
  bool  b1 : 1;
  bool  b2 : 1;
  bool  b3 : 1;

  static const Int RELATION_LENGTH = 8; // Length of the group relation
  static const Int INVERSES_DISTANCE = 4; // How many elements between inverses, i.e., inv(a) = (a + INVERSE_DISTANCE) % RELATION_LENGTH

  static std::map<std::string, Int> wmap;

  static std::map<std::string, Int> init_wmap() {
    std::map<std::string, Int> m;

    m["_"] = -1;

    m["0527"] = 0;
    m["052"] = 1;
    m["05"] = 2;
    m["0"] = 3;
    m["03"] = 4;
    m["036"] = 5;

    m["1630"] = 6;
    m["163"] = 7;
    m["16"] = 8;
    m["1"] = 9;
    m["14"] = 10;
    m["147"] = 11;

    m["2741"] = 12;
    m["274"] = 13;
    m["27"] = 14;
    m["2"] = 15;
    m["25"] = 16;
    m["250"] = 17;

    m["3052"] = 18;
    m["305"] = 19;
    m["30"] = 20;
    m["3"] = 21;
    m["36"] = 22;
    m["361"] = 23;

    m["4163"] = 24;
    m["416"] = 25;
    m["41"] = 26;
    m["4"] = 27;
    m["47"] = 28;
    m["472"] = 29;

    m["5274"] = 30;
    m["527"] = 31;
    m["52"] = 32;
    m["5"] = 33;
    m["50"] = 34;
    m["503"] = 35;

    m["6305"] = 36;
    m["630"] = 37;
    m["63"] = 38;
    m["6"] = 39;
    m["61"] = 40;
    m["614"] = 41;

    m["7416"] = 42;
    m["741"] = 43;
    m["74"] = 44;
    m["7"] = 45;
    m["72"] = 46;
    m["725"] = 47;

    return m;
  }

public:
  Hyperbolic_octagon_translation_word()
    : w0(0), w1(0), w2(0), w3(0), b0(false), b1(false), b2(false), b3(false)
  { }

  Hyperbolic_octagon_translation_word(Int x)
    : w0(x), w1(0), w2(0), w3(0), b0(true), b1(false), b2(false), b3(false)
  { }

  Hyperbolic_octagon_translation_word(Int x, Int y)
    : w0(x), w1(y), w2(0), w3(0), b0(true), b1(true), b2(false), b3(false)
  { }

  Hyperbolic_octagon_translation_word(Int x, Int y, Int z)
    : w0(x), w1(y), w2(z), w3(0), b0(true), b1(true), b2(true), b3(false)
  { }

  Hyperbolic_octagon_translation_word(Int x, Int y, Int z, Int t)
    : w0(x), w1(y), w2(z), w3(t), b0(true), b1(true), b2(true), b3(true)
  { }

  Hyperbolic_octagon_translation_word(const Hyperbolic_octagon_translation_word& other)
    : w0(other.w0), w1(other.w1), w2(other.w2), w3(other.w3),
      b0(other.b0), b1(other.b1), b2(other.b2), b3(other.b3)
  { }

  Hyperbolic_octagon_translation_word(const std::vector<Int>& v)
  {
    switch(v.size())
    {
      case 0:
        w0 = 0; b0 = false;
        w1 = 0; b1 = false;
        w2 = 0; b2 = false;
        w3 = 0; b3 = false;
        break;

      case 1:
        w0 = v[0]; b0 = true;
        w1 = 0; b1 = false;
        w2 = 0; b2 = false;
        w3 = 0; b3 = false;
        break;

      case 2:
        w0 = v[0]; b0 = true;
        w1 = v[1]; b1 = true;
        w2 = 0; b2 = false;
        w3 = 0; b3 = false;
        break;

      case 3:
        w0 = v[0]; b0 = true;
        w1 = v[1]; b1 = true;
        w2 = v[2]; b2 = true;
        w3 = 0; b3 = false;
        break;

      default:
        w0 = v[0]; b0 = true;
        w1 = v[1]; b1 = true;
        w2 = v[2]; b2 = true;
        w3 = v[3]; b3 = true;
        break;
    }
  }

  void operator()(Int idx, Int val)
  {
    switch (idx) {
      case 0:
        w0 = val;
        b0 = true;
        break;
      case 1:
        w1 = val;
        b1 = true;
        break;
      case 2:
        w2 = val;
        b2 = true;
        break;
      default:
        b3 = true;
        w3 = val;
    }
  }

  Int operator()(Int i) const
  {
    switch (i) {
      case 0:
        return w0;
      case 1:
        return w1;
      case 2:
        return w2;
      default:
        return w3;
    }
  }

  bool b(Int i) const
  {
    switch (i) {
      case 0:
        return b0;
      case 1:
        return b1;
      case 2:
        return b2;
      default:
        return b3;
    }

  }

  Int length() const
  {
    if(b0 == false) {
      return 0;
    } else {
      if(b1 == false) {
        return 1;
      } else {
        if(b2 == false) {
          return 2;
        } else {
          if(b3 == false) {
            return 3;
          } else {
            return 4;
          }
        }
      }
    }
  }

private:
  void copy_from(Self o)
  {
    this->w0 = o.w0; this->w1 = o.w1; this->w2 = o.w2; this->w3 = o.w3;
    this->b0 = o.b0; this->b1 = o.b1; this->b2 = o.b2; this->b3 = o.b3;
  }

  Int inv(Int v) const
  {
    return ((v + RELATION_LENGTH/2) % RELATION_LENGTH);
  }

  void reduce()
  {
    std::vector<Int> old = this->get_vector();
    std::vector<Int> red;
    Dehn_reductor dehn;
    dehn(red, old);
    Self newone(red);
    this->copy_from(newone);
  }

public:
  int index_in_order() const { return wmap[this->to_string()]; }

public:
  Self inverse() const
  {
    if(b3) {
      Self r(inv(w3), inv(w2), inv(w1), inv(w0));
      r.reduce();
      return r;
    } else {
      if(b2) {
        return Self(inv(w2), inv(w1), inv(w0));
      } else {
        if(b1) {
          return Self(inv(w1), inv(w0));
        } else {
          if(b0) {
            return Self(inv(w0));
          }
        }
      }
    }

    return Self();
  }

  Self append(Self val) const
  {
    if(this->is_identity())
    {
      if(val.is_identity()) {
        return Self();
      } else {
        return Self(val);
      }
    }

    if(val.is_identity())
    {
      if(this->is_identity()) {
        return Self();
      } else {
        return Self(*this);
      }
    }

    std::vector<Int> o1 = this->get_vector();
    std::vector<Int> o2 = val.get_vector();
    for(int i=0; i<o2.size(); ++i)
      o1.push_back(o2[i]);

    std::vector<Int> red;
    Dehn_reductor dehn;
    dehn(red, o1);

    CGAL_assertion(red.size() < 5);

    Self r(red);

    #if defined PROFILING_MODE
      if(r.is_identity()) {
        calls_append_identity++;
      } else {
        calls_append_non_identity++;
      }
    #endif

    return r;
  }

  std::vector<Int> get_vector() const
  {
    std::vector<Int> v;
    if(b0) {
      v.push_back(w0);
    }
    if(b1) {
      v.push_back(w1);
    }
    if(b2) {
      v.push_back(w2);
    }
    if(b3) {
      v.push_back(w3);
    }

    return v;
  }

  std::vector<Int> get_relation_vector() const
  {
    std::vector<Int> v;
    if(b0)
      v.push_back(ridx(w0));

    if(b1)
      v.push_back(ridx(w1));

    if(b2)
      v.push_back(ridx(w2));

    if(b3)
      v.push_back(ridx(w3));

    return v;
  }

  bool is_identity() const { return (b0 == false); }

  std::string to_string() const
  {
    std::string s = "";

    if(b0)
      s += w0 + '0';
    else
      return "_";

    if(b1)
      s += w1 + '0';

    if(b2)
      s += w2 + '0';

    if(b3)
      s += w3 + '0';

    return s;
  }

  Self operator*(const Self& rh) const
  {
    std::vector<Int> old = this->get_vector();
    std::vector<Int> oth = rh.get_vector();
    for(std::size_t i=0; i<oth.size(); ++i)
      old.push_back(oth[i]);

    std::vector<Int> red;
    Dehn_reductor dehn;
    dehn(red, old);
    CGAL_assertion(red.size() < 5);
    Self newone(red);

      return newone;
    }

    Self operator-(const Self& other) const
    {
      Self res = ((*this) * other.inverse());
      res.reduce();
      return res;
    }

  Self& operator=(const Self& other)
  {
    b0 = other.b0; b1 = other.b1; b2 = other.b2; b3 = other.b3;
    w0 = other.w0; w1 = other.w1; w2 = other.w2; w3 = other.w3;
    return *this;
  }

  // Equality
  bool operator==(const Self& other) const
  {
    int N = this->length();
    if(N == other.length())
    {
      for(int i=0; i<N; ++i)
      {
        if(this->operator()(i) != other(i))
          return false;
      }
      return true;
    }
    else
    {
      return false;
    }
  }

  // Inequality
  bool operator!=(const Self& other) const
  {
     return !operator==(other);
  }

  // just to give an ordering
  bool operator<(  const Self& other) const
  {
    return this->index_in_order() < other.index_in_order();
  }
};

template <class Int>
std::map<std::string, Int> Hyperbolic_octagon_translation_word<Int>::wmap = init_wmap();

template <class Int>
std::ostream& operator<<(std::ostream& s, const Hyperbolic_octagon_translation_word<Int>& o) {

  if(o.is_identity())
  {
    s << "_";
    return s;
  }

  for(int i=0; i<4; ++i)
  {
    if(o.b(i)) {
      s << o(i);
    }
  }
  return s;
}

} // namespace CGAL

#endif // CGAL_HYPERBOLIC_OCTAGON_WORD_4
