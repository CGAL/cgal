// Copyright (c) 2010-2018  INRIA Sophia Antipolis, INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Iordan Iordanov
//

#ifndef CGAL_DEHN_HYPERBOLIC_OCTAGON_TRANSLATION_WORD
#define CGAL_DEHN_HYPERBOLIC_OCTAGON_TRANSLATION_WORD

#include <CGAL/license/Periodic_4_hyperbolic_triangulation_2.h>

#include <CGAL/number_utils.h>

#include <iostream>
#include <vector>

namespace CGAL {

template < typename Idx_type >
class Dehn_hyperbolic_octagon_translation_word
{
private:
  typedef Idx_type Word_idx_type;

  // Check whether two Matrix_elements of the group are inverse of each other
  bool are_inverse(Word_idx_type x, Word_idx_type y)
  {
    Word_idx_type idx = x % 4;
    Word_idx_type idy = y % 4;
    bool r = ((idx == idy) && (x != y));
    return r;
  }

  // Recursively eliminate neighboring inverse Matrix_elements present in the word
  void simplify_adjacent_inverses(std::vector<Word_idx_type>& w)
  {
    std::vector<Word_idx_type> t;
    bool reduced = false;
    Word_idx_type N = static_cast<Word_idx_type>(w.size());
    if(N > 1) {
      for(Word_idx_type i = 0; i < N-1; ++i)
      {
        if(!are_inverse(w[i], w[i+1]))
        {
          t.push_back(w[i]);
        }
        else
        {
          reduced = true;
          i++;
        }
      }

      if(!are_inverse(w[N-2], w[N-1]))
        t.push_back(w[N-1]);
      else
        reduced = true;

      if(reduced)
        simplify_adjacent_inverses(t);

      w = t;
    }
  }

  // Computes and returns the next index in the identity Matrix_element chain
  Word_idx_type next_relation_index(Word_idx_type idx)
  {
    return ((idx + 5) % 8);
  }

  // Checks whether y is the next Matrix_element from x
  bool is_next_relation_index(Word_idx_type x, Word_idx_type y)
  {
    return (next_relation_index(y) == x);
  }

  // Given a word, find the largest subsequence of consecutive Matrix_elements it contains.
  // The sequence itself is placed in 'seq', while the index at which it starts is
  // the return argument of the function.
  Word_idx_type longest_relation_subsequence(std::vector<Word_idx_type>& seq,
                                             const std::vector<Word_idx_type>& w)
  {
    Word_idx_type start = 0;
    Word_idx_type mstart = 0;
    Word_idx_type end = 1;
    Word_idx_type max = 1;
    Word_idx_type len = 1;
    std::vector<Word_idx_type> tmp, mvec;
    tmp.push_back(w[0]);
    for(Word_idx_type i=1; i<w.size(); ++i)
    {
      if(is_next_relation_index(w[i], w[i-1]))
      {
        end++;
        len++;
        tmp.push_back(w[i]);
        if(len > max) {
          max = len;
          mvec = tmp;
          mstart = start;
        }
      }
      else
      {
        tmp.clear();
        tmp.push_back(w[i]);
        start = i;
        end = i;
        len = 0;
      }
    }
    seq = mvec;
    return mstart;
  }


  Word_idx_type relation_index_of_inverse(Word_idx_type x)
  {
    return ((x + 4) % 8);
  }

  // Given a word, construct its inverse
  void invert_word(std::vector<Word_idx_type>& w,
                   const std::vector<Word_idx_type>& original)
  {
    w.clear();
    for(int i=static_cast<int>(original.size())-1; i>=0; --i)
      w.push_back(relation_index_of_inverse(original[i]));
  }

  void invert_4_word(std::vector<Word_idx_type>& w,
                     const std::vector<Word_idx_type>& original)
  {
    w.clear();
    for(int i = static_cast<int>(original.size())-1; i>=0; --i)
      w.push_back(original[i]);
  }

  // Given a sequence of consecutive indices, return the complementary set of consecutive indices in mod 8.
  // For instance, if start = 5 and end = 1, the output is the sequence 2, 3, 4
  void complementary_relation_indices(std::vector<Word_idx_type>& v,
                                      Word_idx_type begin, Word_idx_type end)
  {
    std::vector<Word_idx_type> tmp;
    for(Word_idx_type i=next_relation_index(end); i!=begin; i=next_relation_index(i))
      tmp.push_back(i);

    for(int i=static_cast<int>(tmp.size())-1; i>=0; --i)
      v.push_back(tmp[i]);
  }


  void complementary_relation_indices(std::vector<Word_idx_type>& v,
                                      const std::vector<Word_idx_type>& original)
  {
    std::vector<Word_idx_type> tmp;
    complementary_relation_indices(tmp, original[0], original[original.size() - 1]);
    for(int i = static_cast<int>(tmp.size())-1; i>=0; --i)
      v.push_back(tmp[i]);
  }


  bool is_principal(Word_idx_type w)
  {
    return (w == 0 || w == 2 || w == 5 || w == 7);
  }

  // Given a word, identifies the longest subword consisting of consecutive Matrix_elements and substitutes
  // it with its equivalent shorter word. The search is executed in both the original word and its
  // inverse, and the substitution is made considering the longest subword from both cases.
  bool replace_relation_subword(std::vector<Word_idx_type>& w,
                                const std::vector<Word_idx_type>& original)
  {
    bool replaced = false;

    // Handle empty string case
    if(original.size() == 0)
      return replaced; // false

    // Look for longest subword forward
    std::vector<Word_idx_type> lfwd;
    Word_idx_type idxf = longest_relation_subsequence(lfwd, original);
    Word_idx_type Nf = static_cast<Word_idx_type>(lfwd.size());

    // Get inverse of the original word
    std::vector<Word_idx_type> inv;
    invert_word(inv, original);

    // Look for longest subword backwards
    std::vector<Word_idx_type> lbwd;
    Word_idx_type idxb = longest_relation_subsequence(lbwd, inv);
    Word_idx_type Nb = static_cast<Word_idx_type>(lbwd.size());

    // Assign parameters based on results to homogenise the logic
    std::vector<Word_idx_type> word, sub;
    bool is_inverse;
    Word_idx_type N, idx;
    //cout << "Nb = " << Nb << ", Nf = " << Nf << std::endl;
    if(Nb > Nf)
    {
      word = inv;
      sub = lbwd;
      idx = idxb;
      is_inverse = true;
      N = Nb;
    }
    else
    {
      word = original;
      sub = lfwd;
      idx = idxf;
      is_inverse = false;
      N = Nf;
    }

    // Take care of sequences with length greater or equal to 8 -- each chain of length 8
    // is by default equal to the identity Matrix_element and can be directly eliminated.
    while(N >= 8)
    {
      replaced = true;
      std::vector<Word_idx_type> ttt = word;
      word.clear();
      for(Word_idx_type i = 0; i < idx; ++i)
        word.push_back(ttt[i]);

      for(Word_idx_type i = idx + 8; i < ttt.size(); ++i)
        word.push_back(ttt[i]);

      w = word;

      ttt = sub;
      sub.clear();
      for(Word_idx_type i = 8; i < ttt.size(); ++i)
        sub.push_back(ttt[i]);

      N -= 8;
    }

    // Dehn's algorithm substitutes only chains longer that half-circle.
    // Considering equality may lead to infinite loop -- the equivalent
    // of a chain with length 4 is also a chain of length 4, so the
    // substitution becomes cyclic.
    if(N > 4)
    {
      replaced = true;
      std::vector<Word_idx_type> cmpl;
      complementary_relation_indices(cmpl, sub);

      if(is_inverse)
      {
        std::vector<Word_idx_type> tmp;
        invert_word(tmp, cmpl);
        cmpl = tmp;
      }

      for(Word_idx_type i = 0; i < idx; ++i)
        w.push_back(word[i]);

      for(Word_idx_type i = 0; i < cmpl.size(); ++i)
        w.push_back(cmpl[i]);

      for(Word_idx_type i = N + idx; i < word.size(); ++i)
        w.push_back(word[i]);
    }
    else if(N == 4)
    {
      if(is_inverse)
      {
        std::vector<Word_idx_type> tmp;
        invert_4_word(tmp, word);
        w = tmp;
        replaced = true;
      }
    }

    // If we have been working with the inverse of the original, the result has to be inverted.
    //if(is_inverse) {
      std::vector<Word_idx_type> tmp;
      invert_word(tmp, w);
      w = tmp;
    //}

    return replaced;
  }

public:
  Dehn_hyperbolic_octagon_translation_word() {}

  // Applies Dehn's algorithm to a given word. The result is the equivalent irreducible word
  // of the original. The boolean return argument of the function indicates whether the
  // resulting equivalent irreducible word is the identity Matrix_element or not.
  bool operator()(std::vector<Word_idx_type>& w, std::vector<Word_idx_type> const original)
  {
    std::vector<Word_idx_type> tmp;
    tmp = original;
    while(tmp.size() > 0) {
      simplify_adjacent_inverses(tmp);
      std::vector<Word_idx_type> sub;
      bool replaced = replace_relation_subword(sub, tmp);
      if(!replaced) {
        w = tmp;
        return false;
      }
      tmp = sub;
    }

    w = tmp;

    return true;
  }
}; // class Dehn_hyperbolic_octagon_translation_word

} // namespace CGAL

#endif  // CGAL_DEHN_HYPERBOLIC_OCTAGON_TRANSLATION_WORD
