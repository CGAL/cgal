// Copyright (c) 1999-2004,2006-2009,2014-2018   INRIA Sophia Antipolis, INRIA Nancy - Grand Est (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Iordan Iordanov <iordan.iordanov@loria.fr>

#ifndef CGAL_HYPERBOLIC_OCTAGON_TRANSLATION
#define CGAL_HYPERBOLIC_OCTAGON_TRANSLATION

#include <CGAL/license/Periodic_4_hyperbolic_triangulation_2.h>

#include <CGAL/internal/Hyperbolic_octagon_translation_matrix.h>
#include <CGAL/internal/Hyperbolic_octagon_translation_word.h>
#include <CGAL/internal/Exact_complex.h>
#include <CGAL/Exact_algebraic.h>
#include <CGAL/tss.h>

#include <iostream>
#include <vector>

namespace CGAL {

template <typename FT = Exact_algebraic>
class Hyperbolic_octagon_translation
{
public:
  typedef unsigned short int                                    Word_letter;

  enum Generator {
    A = 0,
    B_BAR = 1,
    C = 2,
    D_BAR = 3,
    A_BAR = 4,
    B = 5,
    C_BAR = 6,
    D = 7
  };

private:
  typedef Exact_complex<FT>                                     ECplx;
  typedef Hyperbolic_octagon_translation<FT>                    Self;
  typedef Hyperbolic_octagon_translation_word<Word_letter>      Word;
  typedef Hyperbolic_octagon_translation_matrix<ECplx>          Matrix;

  Word _wrd;

  static auto initialize_gmap() {
    std::map<std::string, Matrix> m;
    std::vector<Matrix> g;
    Matrix::generators(g);

    m["_"] = Matrix();

    m["0527"] = g[A]*g[B]*g[C]*g[D];
    m["052"] = g[A]*g[B]*g[C];
    m["05"] = g[A]*g[B];
    m["0"] = g[A];
    m["03"] = g[A]*g[D_BAR];
    m["036"] = g[A]*g[D_BAR]*g[C_BAR];

    m["1630"] = g[B_BAR]*g[C_BAR]*g[D_BAR]*g[A];
    m["163"] = g[B_BAR]*g[C_BAR]*g[D_BAR];
    m["16"] = g[B_BAR]*g[C_BAR];
    m["1"] = g[B_BAR];
    m["14"] = g[B_BAR]*g[A_BAR];
    m["147"] = g[B_BAR]*g[A_BAR]*g[D];

    m["2741"] = g[C]*g[D]*g[A_BAR]*g[B_BAR];
    m["274"] = g[C]*g[D]*g[A_BAR];
    m["27"] = g[C]*g[D];
    m["2"] = g[C];
    m["25"] = g[C]*g[B];
    m["250"] = g[C]*g[B]*g[A];

    m["3052"] = g[D_BAR]*g[A]*g[B]*g[C];
    m["305"] = g[D_BAR]*g[A]*g[B];
    m["30"] = g[D_BAR]*g[A];
    m["3"] = g[D_BAR];
    m["36"] = g[D_BAR]*g[C_BAR];
    m["361"] = g[D_BAR]*g[C_BAR]*g[B_BAR];

    m["4163"] = g[A_BAR]*g[B_BAR]*g[C_BAR]*g[D_BAR];
    m["416"] = g[A_BAR]*g[B_BAR]*g[C_BAR];
    m["41"] = g[A_BAR]*g[B_BAR];
    m["4"] = g[A_BAR];
    m["47"] = g[A_BAR]*g[D];
    m["472"] = g[A_BAR]*g[D]*g[C];

    m["5274"] = g[B]*g[C]*g[D]*g[A_BAR];
    m["527"] = g[B]*g[C]*g[D];
    m["52"] = g[B]*g[C];
    m["5"] = g[B];
    m["50"] = g[B]*g[A];
    m["503"] = g[B]*g[A]*g[D_BAR];

    m["6305"] = g[C_BAR]*g[D_BAR]*g[A]*g[B];
    m["630"] = g[C_BAR]*g[D_BAR]*g[A];
    m["63"] = g[C_BAR]*g[D_BAR];
    m["6"] = g[C_BAR];
    m["61"] = g[C_BAR]*g[B_BAR];
    m["614"] = g[C_BAR]*g[B_BAR]*g[A_BAR];

    m["7416"] = g[D]*g[A_BAR]*g[B_BAR]*g[C_BAR];
    m["741"] = g[D]*g[A_BAR]*g[B_BAR];
    m["74"] = g[D]*g[A_BAR];
    m["7"] = g[D];
    m["72"] = g[D]*g[C];
    m["725"] = g[D]*g[C]*g[B];

    { // This block abuses `operator<<` of numbers, to a null stream.
      // That ensures that the following memory pool are correctly
      // initialized:
      //   - `CORE::MemoryPool<CORE::Realbase_for<long, 1024>`
      //   - `CORE::MemoryPool<CORE::Realbase_for<double, 1024>`
      //   - `CORE::MemoryPool<CORE::BigFloatRep, 1024>`
      //   - `CORE::MemoryPool<CORE::BigIntRep, 1024>`
      // otherwise, there is an assertion during the destruction of
      // static (or `thread_local`) objects
      struct NullBuffer : public std::streambuf {
        int overflow(int c) { return c; }
      };
      NullBuffer null_buffer;
      std::ostream null_stream(&null_buffer);
      for(auto& pair: m) null_stream << pair.second;
    }
    return m;
  }

  static const Matrix& gmap(const std::string& s)
  {
    typedef std::map<std::string, Matrix>  M;
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(M, m, initialize_gmap());
    return m[s];
  }

public:
  Hyperbolic_octagon_translation() : _wrd() {}
  Hyperbolic_octagon_translation(Word w) : _wrd(w) {}
  Hyperbolic_octagon_translation(Word_letter w1) : _wrd(w1) {}
  Hyperbolic_octagon_translation(Word_letter w1, Word_letter w2) : _wrd(w1,w2) {}
  Hyperbolic_octagon_translation(Word_letter w1, Word_letter w2, Word_letter w3) : _wrd(w1,w2,w3) {}
  Hyperbolic_octagon_translation(Word_letter w1, Word_letter w2, Word_letter w3, Word_letter w4) : _wrd(w1,w2,w3,w4) {}

  std::pair<FT,FT> alpha() const
  {
    const Matrix& _m = gmap(_wrd.to_string());
    ECplx _a = _m.alpha();
    FT ax = _a.real();
    FT ay = _a.imag();

    return std::pair<FT,FT>(ax,ay);
  }

  std::pair<FT,FT> beta() const
  {
    const Matrix& _m = gmap(_wrd.to_string());
    ECplx _b = _m.beta();
    FT bx = _b.real();
    FT by = _b.imag();
    return std::pair<FT,FT>(bx,by);
  }

  bool is_identity() const { return _wrd.is_identity(); }
  Self inverse() const { return Self(_wrd.inverse()); }

  Self operator*(const Self& rh) const { return Self(this->_wrd * rh._wrd); }
  Self operator-(const Self& other) const { return Self(this->_wrd - other._wrd); }

  bool operator==(const Hyperbolic_octagon_translation<FT>& other) const {
    return this->_wrd == other._wrd;
  }
  bool operator!=(const Hyperbolic_octagon_translation<FT>& other) const {
    return this->_wrd != other._wrd;
  }
  bool operator<(const Hyperbolic_octagon_translation<FT>& other) const {
    return this->_wrd < other._wrd;
  }

  std::string to_string() const { return _wrd.to_string(); }

  static Self generator(const Word_letter wl) { return Self(wl); }

  static void generators(std::vector<Self>& gens)
  {
    gens.push_back(Self(A));
    gens.push_back(Self(B_BAR));
    gens.push_back(Self(C));
    gens.push_back(Self(D_BAR));
    gens.push_back(Self(A_BAR));
    gens.push_back(Self(B));
    gens.push_back(Self(C_BAR));
    gens.push_back(Self(D));
  }
};


template <typename FT>
std::ostream& operator<<(std::ostream& s, const Hyperbolic_octagon_translation<FT>& tr) {
  s << tr.to_string();
  return s;
}

} // namespace CGAL

#endif  // CGAL_HYPERBOLIC_OCTAGON_TRANSLATION
