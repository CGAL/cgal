// Copyright (c) 2010   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
// 
//
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_MATRIX_H
#define CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_MATRIX_H

#include <complex>
#include <iostream>
#include <vector>
#include <iterator>
#include <map>
#include <set>
#include <assert.h>
#include <string>
#include <fstream>

#include <CGAL/Square_root_2_field.h>
#include <CGAL/number_utils.h>

typedef unsigned short int                                                  Word_idx_type;


template <class Word_idx_type>
std::ostream& operator<<(std::ostream& s, const std::vector<Word_idx_type>& o) {
  for (int i = 0; i < o.size(); i++) {
    s << o[i];
  }
  return s;
} 


template<class GT>
class Hyperbolic_octagon_translation_matrix
{

  typedef typename GT::FT                                   FT;

public:
  

  typedef Square_root_2_field<int>                          Field_number;
  typedef Hyperbolic_octagon_translation_matrix<GT>         Self;  
  typedef std::complex<Field_number>                        Matrix_element;

  Field_number factor;    // The multiplicative factor present in the general formula: sqrt(2) - 1

  Matrix_element  A;
  Matrix_element  B;
  std::string     label;  

  Hyperbolic_octagon_translation_matrix(const Matrix_element& A_, const Matrix_element& B_, 
    const std::string& label_ = std::string("") ) :
  A( A_ ), B( B_ ), label( label_ ), factor(-1, 1) {}

  Hyperbolic_octagon_translation_matrix() :
  A( Matrix_element(1,0) ), B( Matrix_element(0,0) ), label( std::string("") ), factor(-1, 1) {}

  std::string merge_labels(const Self& rh) const
  {
    return label + "*" + rh.label;
  }

  Self operator*(const Self& rh) const
  {
    return Self(  A*rh.A + factor*B*conj(rh.B), 
      A*rh.B + B*conj(rh.A), merge_labels(rh) );
  }
  
  Self inverse() const
  {
    Self inv = Self(conj(A), -B);
    std::string inv_label;
    for(long i = label.size() - 1; i >= 0; i--) {
      if(label[i] >= 'a' && label[i] <= 'd') {
        inv_label.push_back('\\');
        inv_label.push_back('b');
        inv_label.push_back('a');
        inv_label.push_back('r');
        inv_label.push_back('{');
        inv_label.push_back(label[i]);
        inv_label.push_back('}');
      }
      if(label[i] == '*') {
        inv_label.push_back('*');
      }
      if(label[i] == '1') {
        assert(i - 3 >= 0);
        assert(label[i - 3] >= 'a' && label[i - 3] <= 'd');
        inv_label.push_back(label[i - 3]);
        i = i - 3;
      }
    }
    inv.label = inv_label;
    return inv;
  }

  // rotation \pi/4
  Hyperbolic_octagon_translation_matrix rotate() const
  {
    Field_number B1 = real(B);
    Field_number B2 = imag(B);
    
    // sqrt(2)
    Field_number k = Field_number(0, 1);
    
    Field_number BB1 = (B1 - B2)*k;
    Field_number BB2 = (B1 + B2)*k;

    assert(BB2.l % 2 == 0 && BB2.r % 2 == 0);

    BB1.l = BB1.l/2;
    BB1.r = BB1.r/2;
    BB2.l = BB2.l/2;
    BB2.r = BB2.r/2;

    return Hyperbolic_octagon_translation_matrix(A, Matrix_element(BB1, BB2)); 
  }

  Field_number trace() const 
  {
    return Field_number(2, 0)*real(A);
  }

  double length() const
  {
    typedef long double ld;

    ld l = real(A).l;
    ld r = real(A).r;
    ld tr = l + sqrt(2.)*r;
    if (tr < 0) {
      tr = -tr;
    }

    return 2.*acosh(tr);
  }

  // determinant == 1
  Matrix_element det() const
  {
    return norm(A) - factor * norm(B);
  }
  
  static std::complex<double> toComplexDouble(Matrix_element M) //const
  {
    Field_number rl  = real(M);
    Field_number img = imag(M);
    
    return std::complex<double>(rl.l + sqrt(2.)*rl.r, img.l + sqrt(2.)*img.r);
  }

  std::pair<double, double> apply(double x, double y)
  {
    typedef std::complex<double> Cmpl;
    Cmpl Aa = toComplexDouble(A);
    Cmpl Bb = toComplexDouble(B);

    double ax = sqrt(factor.l + sqrt(2.)*factor.r) ;
    
    Cmpl z(x, y);
    Cmpl res = (Aa*z + ax*Bb)/(ax*(conj(Bb)*z) + conj(Aa));
    return std::pair<double, double>(real(res), imag(res)); 
  }


  template<class Point>
  Point apply(Point p) {
    pair<double, double> arg = apply(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
    return Point(arg.first, arg.second);
  }


};



//template< class GT, template <class> class Hyperbolic_octagon_translation_matrix>
//typename Hyperbolic_octagon_translation_matrix<GT>::Field_number Hyperbolic_octagon_translation_matrix<GT>::factor = Hyperbolic_octagon_translation_matrix<GT>:: Field_number(-1, 1); 


// just to give an order(ing)
template< class GT, template <class> class Hyperbolic_octagon_translation_matrix >
bool operator < (const std::complex< typename Hyperbolic_octagon_translation_matrix<GT>:: template Field_number<GT::FT> >& lh,
  const std::complex< typename Hyperbolic_octagon_translation_matrix<GT>:: template Field_number<GT::FT> >& rh)
{
  if (real(lh) < real(rh)) {
    return true;
  }

  if (real(lh) == real(rh)) {
    if (imag(lh) < imag(rh)) {
      return true;
    }
  }

  return false;   
}

// just to order octagon_matrices 
template<class Field_number>
bool operator < (const Hyperbolic_octagon_translation_matrix<Field_number>& lh, 
  const Hyperbolic_octagon_translation_matrix<Field_number>& rh)
{
  if (lh.A < rh.A) {
    return true;
  }
  
  if (lh.A == rh.A ) {
    if (lh.B < rh.B) {
      return true;
    }
  }

  return false;
}

template<class Field_number>
bool operator == (const Hyperbolic_octagon_translation_matrix<Field_number>& lh, 
  const Hyperbolic_octagon_translation_matrix<Field_number>& rh)
{
  return (lh.A == rh.A && lh.B == rh.B);
}

template<class Field_number>
std::ostream& operator<<(std::ostream& os, const Hyperbolic_octagon_translation_matrix<Field_number>& m)
{
  os << m.A << " " << m.B;
  return os;
}


enum Direction {
  DIRECTION_A = 0,  // 0
  DIRECTION_B_BAR,  // 1
  DIRECTION_C,      // 2
  DIRECTION_D_BAR,  // 3
  DIRECTION_A_BAR,  // 4
  DIRECTION_B,      // 5
  DIRECTION_C_BAR,  // 6
  DIRECTION_D       // 7
};



template < class GT, template <class> class Hyperbolic_octagon_translation_matrix >
void get_generators(std::vector< Hyperbolic_octagon_translation_matrix<GT> >& gens)
{
  typedef Hyperbolic_octagon_translation_matrix<GT>           Matrix;
  typedef typename Matrix::Matrix_element                     Matrix_element;
  typedef typename Matrix::Field_number                       Field_number;
  // This is a in the matrix, equal to sqrt(2) + 1
  Matrix_element A = Matrix_element(Field_number(1, 1), Field_number(0, 0));

  // This vector holds all other Matrix_elements, results of the exponentials for various k
  std::vector<Matrix_element> B(8, Matrix_element(Field_number(0, 0), Field_number(0, 0)));

  // Corrected ordering (initial ordering is present in backup file)
  B[DIRECTION_A]     = A * Matrix_element(Field_number( 0,  1), Field_number( 0,  0));
  B[DIRECTION_B]     = A * Matrix_element(Field_number(-1,  0), Field_number(-1,  0));
  B[DIRECTION_C]     = A * Matrix_element(Field_number( 0,  0), Field_number( 0,  1));
  B[DIRECTION_D]     = A * Matrix_element(Field_number( 1,  0), Field_number(-1,  0));
  B[DIRECTION_A_BAR] = A * Matrix_element(Field_number( 0, -1), Field_number( 0,  0));
  B[DIRECTION_B_BAR] = A * Matrix_element(Field_number( 1,  0), Field_number( 1,  0));
  B[DIRECTION_C_BAR] = A * Matrix_element(Field_number( 0,  0), Field_number( 0, -1));
  B[DIRECTION_D_BAR] = A * Matrix_element(Field_number(-1,  0), Field_number( 1,  0));
  
  std::string labels[8] = {string("a"), string("\\bar{b}"), string("c"), string("\\bar{d}"),
  string("\\bar{a}"), string("b"), string("\\bar{c}"), string("d")};
  for(int i = 0; i < 8; i++) {
    gens.push_back(Matrix(A, B[i], labels[i]));
  }
}





/***************************************************************/



  // Check whether two Matrix_elements of the group are inverse of each other
  bool Dehn_are_inverse(Word_idx_type x, Word_idx_type y) {
    Word_idx_type idx = x % 4;
    Word_idx_type idy = y % 4;
    bool r = ((idx == idy) && (x != y));
    return r;
  }



  // Recursively eliminate neighboring inverse Matrix_elements present in the word
  void Dehn_simplify_adjacent_inverses(vector<Word_idx_type>& w) {
    vector<Word_idx_type> t;
    bool reduced = false;
    Word_idx_type N = w.size();
    if (N > 1) {
      for (Word_idx_type i = 0; i < N-1; i++) {
        if (!Dehn_are_inverse(w[i], w[i+1])) {
          t.push_back(w[i]);
        } else {
          reduced = true;
          i++;
        }
      }
      if (!Dehn_are_inverse(w[N-2], w[N-1])) {
        t.push_back(w[N-1]);
      } else {
        reduced = true;
      }
    
      if (reduced) {
        Dehn_simplify_adjacent_inverses(t);
      }

      w = t;
    }
  }


  // Computes and returns the next index in the identity Matrix_element chain
  Word_idx_type Dehn_next_relation_index(Word_idx_type idx) {
    return ( (idx + 5) % 8 );
  }


  // Checks whether y is the next Matrix_element from x
  bool Dehn_is_next_relation_index(Word_idx_type x, Word_idx_type y) {
    return ( Dehn_next_relation_index(y) == x);
  }


  // Given a word, find the largest subsequence of consecutive Matrix_elements it contains.
  // The sequence itself is placed in 'seq', while the index at which it starts is
  // the return argument of the function.
  Word_idx_type Dehn_longest_relation_subsequence(vector<Word_idx_type>& seq, vector<Word_idx_type> const w) {
    Word_idx_type start = 0; 
    Word_idx_type mstart = 0;
    Word_idx_type end = 1;
    Word_idx_type max = 1;
    Word_idx_type len = 1;
    vector<Word_idx_type> tmp, mvec;
    tmp.push_back(w[0]);
    for (Word_idx_type i = 1; i < w.size(); i++) {
      if ( Dehn_is_next_relation_index( w[i], w[i-1] ) ) {
        end++;
        len++;
        tmp.push_back(w[i]);
        if (len > max) {
          max = len;
          mvec = tmp;
          mstart = start;
        }
      } else {
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


  Word_idx_type Dehn_relation_index_of_inverse(Word_idx_type x) {
    return ((x + 4) % 8);
  }


  // Given a word, construct its inverse
  void Dehn_invert_word(vector<Word_idx_type>& w, vector<Word_idx_type> const original) {
    w.clear();
    for (int i = original.size() - 1; i >= 0; i--) {
      w.push_back( Dehn_relation_index_of_inverse(original[i]) );
    }
  }



  // Given a sequence of consecutive indices, return the complementary set of consecutive indices in mod 8.
  // For instance, if start = 5 and end = 1, the output is the sequence 2, 3, 4
  void Dehn_complementary_relation_indices(vector<Word_idx_type>& v, Word_idx_type begin, Word_idx_type end) {
    vector<Word_idx_type> tmp;
    for (Word_idx_type i = Dehn_next_relation_index(end); i != begin; i = Dehn_next_relation_index(i)) {
      tmp.push_back(i);
    }

    for (Word_idx_type i = tmp.size() - 1; i >= 0; i--) {
      v.push_back(tmp[i]);
    }
  }


  void Dehn_complementary_relation_indices(vector<Word_idx_type>& v, vector<Word_idx_type> const original) {
    std::vector<Word_idx_type> tmp;
    Dehn_complementary_relation_indices(tmp, original[0], original[original.size() - 1]);
    for (int i = tmp.size() - 1; i >= 0; i--) {
      v.push_back(tmp[i]);
    }
  }


  // Given a word, identifies the longest subword consisting of consecutive Matrix_elements and substitutes 
  // it with its equivalent shorter word. The search is executed in both the original word and its
  // inverse, and the substitution is made considering the longest subword from both cases.
  bool Dehn_replace_relation_subword(vector<Word_idx_type>& w, const vector<Word_idx_type> original) {

    bool replaced = false;

    // Handle empty string case
    if (original.size() == 0) {
      return replaced; // false
    }

    // Look for longest subword forward
    vector<Word_idx_type> lfwd;
    Word_idx_type idxf = Dehn_longest_relation_subsequence(lfwd, original);
    Word_idx_type Nf = lfwd.size();

    // Get inverse of the original word
    vector<Word_idx_type> inv;
    Dehn_invert_word(inv, original);

    // Look for longest subword backwards
    vector<Word_idx_type> lbwd;
    Word_idx_type idxb = Dehn_longest_relation_subsequence(lbwd, inv);
    Word_idx_type Nb = lbwd.size();

    // Assign parameters based on results to homogenise the logic
    vector<Word_idx_type> word, sub;
    bool is_inverse;
    Word_idx_type N, idx;
    //cout << "Nb = " << Nb << ", Nf = " << Nf << endl;
    if (Nb > Nf) {
      word = inv;
      sub = lbwd;
      idx = idxb;
      is_inverse = true;
      N = Nb;
    } else {
      word = original;
      sub = lfwd;
      idx = idxf;
      is_inverse = false;
      N = Nf;
    }

    // Take care of sequences with length greater or equal to 8 -- each chain of length 8
    // is by default equal to the identity Matrix_element and can be directly eliminated.
    while (N >= 8) {
      replaced = true;
      vector<Word_idx_type> ttt = word;
      word.clear();
      for (Word_idx_type i = 0; i < idx; i++) {
        word.push_back(ttt[i]);
      }
      for (Word_idx_type i = idx + 8; i < ttt.size(); i++) {
        word.push_back(ttt[i]);
      }
      w = word;

      ttt = sub;
      sub.clear();
      for (Word_idx_type i = 8; i < ttt.size(); i++) {
        sub.push_back(ttt[i]);
      }
      N -= 8;
    }

    // Dehn's algorithm substitutes only chains longer that half-circle. 
    // Considering equality may lead to infinite loop -- the equivalent 
    // of a chain with length 4 is also a chain of length 4, so the 
    // substitution becomes cyclic.
    if (N > 4) {
      replaced = true;
      vector<Word_idx_type> cmpl;
      Dehn_complementary_relation_indices(cmpl, sub);

      if(is_inverse) {
        vector<Word_idx_type> tmp;
        Dehn_invert_word(tmp, cmpl);
        cmpl = tmp;
      }

      for (Word_idx_type i = 0; i < idx; i++) {
        w.push_back(word[i]);
      }
      for (Word_idx_type i = 0; i < cmpl.size(); i++) {
        w.push_back(cmpl[i]);
      }
      for (Word_idx_type i = N + idx; i < word.size(); i++) {
        w.push_back(word[i]);
      }
    } else if (N == 4) {
      if (is_inverse) {
        vector<Word_idx_type> tmp;
        Dehn_invert_word(tmp, word);
        w = tmp;
        replaced = true;
      }
    }

    // If we have been working with the inverse of the original, the result has to be inverted.
    if (is_inverse) {
      vector<Word_idx_type> tmp;
      Dehn_invert_word(tmp, w);
      w = tmp;
    }

    return replaced;
  }


  // Applies Dehn's algorithm to a given word. The result is the equivalent ireducible word
  // of the original. The boolean return argument of the function indicates whether the
  // resulting equivalent irreducible word is the identity Matrix_element or not.
  bool Dehn_reduce_word(vector<Word_idx_type>& w, vector<Word_idx_type> const original) {
    bool is_identity = false;
    vector<Word_idx_type> tmp;
    tmp = original;
    while (tmp.size() > 0) {
      Dehn_simplify_adjacent_inverses(tmp);
      vector<Word_idx_type> sub;
      bool replaced = Dehn_replace_relation_subword(sub, tmp);
      if (!replaced) {
        w = tmp;
        return false;
      }
      tmp = sub;
    }

    w = tmp;

    return true;
  }


  // Applies Dehn's algorithm to a given word. The result is the equivalent ireducible word
  // of the original. The boolean return argument of the function indicates whether the
  // resulting equivalent irreducible word is the identity Matrix_element or not.
  // Overload for generic type.
  template <class Word>
  bool Dehn_reduce_word(Word& w, Word const original) {

    std::vector<Word_idx_type> in = original.get_vector();
    std::vector<Word_idx_type> out;
    bool result = Dehn_reduce_word(out, in);
    w = Word(out);
    return result;

  }

/**************************************************************/


template < class GT, template <class> class Hyperbolic_octagon_translation_matrix >
void generate_words(const vector< Hyperbolic_octagon_translation_matrix<GT> > gens,
                          set<    Hyperbolic_octagon_translation_matrix<GT> >& words, 
                          vector< Hyperbolic_octagon_translation_matrix<GT> >& prev, 
                          int depth, double threshold  )
{
  if (depth == 1) {
    for(int i = 0; i < 8; i++) {
      words.insert( gens[i] );
      prev.push_back( gens[i] );
    }
    return;
  }

  typedef Hyperbolic_octagon_translation_matrix<GT>           Matrix;
  typedef typename Matrix::Matrix_element                     Matrix_element;
  typedef typename Matrix::Field_number                       Field_number;

  vector<Matrix> els;
  generate_words( words, els, depth - 1, threshold);
  
  Matrix temp = Matrix(Matrix_element(), Matrix_element());
  int size = els.size();
  bool is_new = false;
  for(int k = 0; k < size; k++) {
    for(int i = 0; i < 8; i++) {
      temp = els[k]*gens[i];

      if(temp.length() > threshold /*15.*/) {
      continue;
    }     

    is_new = words.insert(temp).second;
    if(is_new == true) {
      prev.push_back(temp);
    }
  }
}
}




// does the axis of a given matrix go through the fundamental octagon 
template < class GT, template <class> class Hyperbolic_octagon_translation_matrix >
bool IsCanonical(const Hyperbolic_octagon_translation_matrix<GT> & m)
{

  typedef Hyperbolic_octagon_translation_matrix<GT>           Matrix;
  typedef typename Matrix::Matrix_element                     Matrix_element;
  typedef typename Matrix::Field_number                       Field_number;

  Matrix temp = m;
  
  // rotate while |B1| < |B2|
  Field_number B1, B2;
  Field_number C = Field_number(-1, -1);
  for(int i = 0; i < 8 && C != C.abs(); i++) {
    B1 = real(temp.B).abs();
    B2 = imag(temp.B).abs();
    C = B1 - B2;
    
    temp = temp.rotate();
  }
  assert(C == C.abs());

  // (2 - sqrt(2))(|B1| + (sqrt(2) -  1)|B2|)
  Field_number right = Field_number(2, -1)*(B1 + Field_number(-1, 1)*B2);
  
  // |A2|
  Field_number left = imag(temp.A).abs();

  // left <= right -> true
  C = right - left;
  return C == C.abs();
}


template < class GT, template <class> class Hyperbolic_octagon_translation_matrix >
void dfs( const vector< Hyperbolic_octagon_translation_matrix<GT> > gens,
          const Hyperbolic_octagon_translation_matrix<GT>  & m, 
          set< Hyperbolic_octagon_translation_matrix<GT> >& visited)
{

  typedef Hyperbolic_octagon_translation_matrix<GT>           Matrix;
  typedef typename Matrix::Matrix_element                     Matrix_element;
  typedef typename Matrix::Field_number                       Field_number;

  assert(IsCanonical(m));
  visited.insert(m);

  Matrix candidate = m;
  for(int i = 0; i < 8; i++) {
    candidate = gens[i]*m*gens[(i + 4) % 8];
    if(IsCanonical(candidate) == true && visited.find(candidate) == visited.end()) {
      dfs(candidate, visited);
    } 
  }
}

// map<Hyperbolic_octagon_translation_matrix m, Hyperbolic_octagon_translation_matrix Aux>, m = Aux * origin * Aux^{-1} 
template < class GT, template <class> class Hyperbolic_octagon_translation_matrix >
void dfs_with_info(const vector< Hyperbolic_octagon_translation_matrix<GT> > gens,
                   const pair<Hyperbolic_octagon_translation_matrix<GT>, Hyperbolic_octagon_translation_matrix<GT> >& new_pair, 
                         map <Hyperbolic_octagon_translation_matrix<GT>, Hyperbolic_octagon_translation_matrix<GT> >& visited)
{

  typedef Hyperbolic_octagon_translation_matrix<GT>           Matrix;
  typedef typename Matrix::Matrix_element                     Matrix_element;
  typedef typename Matrix::Field_number                       Field_number;

  assert(IsCanonical(new_pair.first));  
  visited.insert(new_pair);
  
  const Matrix& current = new_pair.first;
  const Matrix& current_factor = new_pair.second;
  Matrix candidate = current, candidate_factor = current_factor;
  for(int i = 0; i < 8; i++) {
    candidate = gens[i]*current*gens[(i + 4) % 8];
    if(IsCanonical(candidate) == true && visited.find(candidate) == visited.end()) {
      candidate_factor = gens[i]*current_factor;
      dfs_with_info(pair<Matrix, Matrix>(candidate, candidate_factor), visited);
    } 
  }
}


template < class GT, template <class> class Hyperbolic_octagon_translation_matrix >
void dfs_with_info(const vector< Hyperbolic_octagon_translation_matrix<GT> > gens,
                   const         Hyperbolic_octagon_translation_matrix<GT>& origin,
                            map< Hyperbolic_octagon_translation_matrix<GT>, Hyperbolic_octagon_translation_matrix<GT> >& visited)
{

  typedef Hyperbolic_octagon_translation_matrix<GT>           Matrix;
  typedef typename Matrix::Matrix_element                     Matrix_element;
  typedef typename Matrix::Field_number                       Field_number;

  assert(IsCanonical(origin));
  Matrix id = Matrix(Matrix_element(Field_number(1, 0), Field_number(0, 0)), 
                     Matrix_element(Field_number(0, 0), Field_number(0, 0)));
  pair<Matrix, Matrix> new_pair(origin, id);
  
  dfs_with_info(new_pair, visited);
}


template < class GT, template <class> class Hyperbolic_octagon_translation_matrix >
class IntersectionNumber
{

  typedef Hyperbolic_octagon_translation_matrix<GT>           Matrix;
  typedef typename Matrix::Matrix_element                     Matrix_element;
  typedef typename Matrix::Field_number                       Field_number;

public:

  struct Point
  {
    Point(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}

    double x, y;
  };
  
  Matrix m;

  IntersectionNumber(const Matrix& m_) : m(m_)
  {}
  
  long operator() () const
  {
    set<Matrix> visited;
    typename set<Matrix>::iterator it, it2;
    map<long, long> nb_map;
    map<long, long>::iterator mit;

    dfs(m, visited);

    set<pair< Matrix, Matrix> > common;
    for(it = visited.begin(); it != visited.end(); ++it) {
      for(it2 = it; it2 != visited.end(); ++it2) {
        if(*it == *it2) {
          continue;
        }
        if(haveIntersection(*it, *it2) == true) {
          common.clear();
          count_nb(*it, *it2, common);

          mit = nb_map.find(common.size());
          if(mit != nb_map.end()) {
            mit->second += 1;
          } else {
            nb_map.insert(pair<long, long>(common.size(), 1));
          }
        }
      }
    }

    long nb = 0;
    for(mit = nb_map.begin(); mit != nb_map.end(); mit++) {
      assert( mit->second % mit->first == 0 ); 
      nb += mit->second/mit->first;
    }
    return nb;
  }
  
  
  void count_nb(const vector<Matrix> gens, const Matrix& m1, const Matrix& m2, set< pair< Matrix, Matrix> >& visited) const
  {
    typedef pair<Matrix, Matrix> matrix_pair;
    visited.insert(matrix_pair(m1, m2));
    Matrix c1 = m1, c2 = m2;
    for(int i = 0; i < 8; i++) {
      c1 = gens[i]*m1*gens[(i + 4) % 8];
      c2 = gens[i]*m2*gens[(i + 4) % 8];
      if(IsCanonical(c1) == true && IsCanonical(c2) == true && visited.find(matrix_pair(c1, c2)) == visited.end()) {
        count_nb(c1, c2, visited);
      } 
    }
  }

//private:

// check whether two axis have intersection
  bool haveIntersection(const Matrix& m1, const Matrix& m2) const
  {
    Point p1, p2;
    intersectWithInfinity(m1, p1, p2);

    Point p3, p4;
    intersectWithInfinity(m2, p3, p4);

  // orientation test
    double sign1 = (p1.x - p3.x)*(p2.y - p3.y) - (p1.y - p3.y)*(p2.x - p3.x);
    double sign2 = (p1.x - p4.x)*(p2.y - p4.y) - (p1.y - p4.y)*(p2.x - p4.x);

    assert( sign1 * sign2 != 0);
    return (sign1 * sign2 < 0);
  }
  
  void intersectWithInfinity(const Matrix& m, Point& p1, Point& p2) const
  {
    Matrix_element a = m.A, b = m.B, factor = m.factor;

    Matrix_element four = Matrix_element(Field_number(4, 0), Field_number(0, 0));
    Matrix_element two = Matrix_element(Field_number(2, 0), Field_number(0, 0));

    Matrix_element D = (a - conj(a))*(a - conj(a));
    D += four*b*conj(b)*factor;
    Matrix_element T1 = conj(a) - a;
    Matrix_element T2 = two*conj(b);

    complex<double> d = m.toComplexDouble(D);
    complex<double> t1 = m.toComplexDouble(T1);
    complex<double> t2 = m.toComplexDouble(T2);
    complex<double> au = complex<double>(m.factor.l + sqrt(2.)*m.factor.r, 0);

    complex<double> z1 = (t1 + sqrt(d))/(t2*sqrt(au));
    complex<double> z2 = (t1 - sqrt(d))/(t2*sqrt(au));

    p1 = Point(real(z1), imag(z1));
    p2 = Point(real(z2), imag(z2));

    assert(p1.x*p1.x + p1.y*p1.y > 0.99 && p1.x*p1.x + p1.y*p1.y < 1.01);
    assert(p2.x*p2.x + p2.y*p2.y > 0.99 && p2.x*p2.x + p2.y*p2.y < 1.01);
  }

  
};


template < class GT, template <class> class Hyperbolic_octagon_translation_matrix >
void Delete(const set< Hyperbolic_octagon_translation_matrix<GT> >& canonical_set, 
               vector< Hyperbolic_octagon_translation_matrix<GT> >& output)
{

  typedef Hyperbolic_octagon_translation_matrix<GT>           Matrix;
  typedef typename Matrix::Matrix_element                     Matrix_element;
  typedef typename Matrix::Field_number                       Field_number;

  set<Matrix> redundant;

  typename set<Matrix>::iterator it;
  for(it = canonical_set.begin(); it != canonical_set.end(); ++it) {
    if(redundant.find(*it) != redundant.end()) { 
      continue;
    }

    set<Matrix> visited;
    dfs(*it, visited);
    visited.erase(*it);
    
    redundant.insert(visited.begin(), visited.end());    
    output.push_back(*it);
  }
}


template < class GT, template <class> class Hyperbolic_octagon_translation_matrix >
void generate_unique_words( const vector< Hyperbolic_octagon_translation_matrix<GT> >  gens,
                                  vector< Hyperbolic_octagon_translation_matrix<GT> >& output, 
                                  double threshold = 10, int word_length = 13)
{

  typedef Hyperbolic_octagon_translation_matrix<GT>           Matrix;
  typedef typename Matrix::Matrix_element                     Matrix_element;
  typedef typename Matrix::Field_number                       Field_number;

  get_generators(gens);
  
  set<Matrix> unique_words;
  vector<Matrix> temp;
  generate_words(unique_words, temp, word_length, threshold);
  
  double l = 0;
  typename set<Matrix>::iterator uit;
  for(uit = unique_words.begin(); uit != unique_words.end(); ++uit) {
    l = uit->length();
    if(0. < l && l < threshold) {
      output.push_back( *uit );
    }
  }
  
  cout << "nb of unique words " << output.size() << endl;
}

// words that correspond to union of 1-cycles
template < class GT, template <class> class Hyperbolic_octagon_translation_matrix >
void generate_words_union_1_cycles( const vector< Hyperbolic_octagon_translation_matrix<GT> >  gens,
                                          vector< Hyperbolic_octagon_translation_matrix<GT> >& out)
{

  typedef Hyperbolic_octagon_translation_matrix<GT>           Matrix;
  typedef typename Matrix::Matrix_element                     Matrix_element;
  typedef typename Matrix::Field_number                       Field_number;

  if(gens.size() == 0) {
    get_generators(gens);
  }
  Matrix f[4] = {gens[0], gens[5], gens[2], gens[7]};
  
  Matrix F[8] = {
    f[0]*f[0].inverse(),
    f[0],
    f[0]*f[1],
    f[0]*f[1]*f[2],
    f[0]*f[1]*f[2]*f[3],
    f[3]*f[2]*f[1],
    f[3]*f[2],
    f[3]
  };
  
  long counter = 0;
  for(int i = 1; i < 5; i++) {
    for(int j = 0; j < 8; j++) {
      for(int k = 0; k < 8; k++) {
        if (j == k) {
          continue;
        }
        // intersection
        if ((0 < j && j < i) && (i < k)) {
          continue;
        }
        if ((0 < k && k < i) && (i < j)) {
          continue;
        }
        counter++;
        
        Matrix Tr = F[i]*F[k].inverse()*F[j];
        // check that Tr != Id
        if(Tr.length() > 2.) {
          out.push_back(Tr);
        }
      }
    }
  }
  cout << counter << endl;
  
}



#endif  // CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_MATRIX_H



