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

#ifndef CGAL_HYPERBOLIC_OCTAGON_GROUP_H
#define CGAL_HYPERBOLIC_OCTAGON_GROUP_H


#include <complex>
#include <iostream>
#include <vector>
#include <iterator>
#include <map>
#include <set>
#include <assert.h>
#include <string>
#include <fstream>


template<class Square_root_2_field>
class Hyperbolic_octagon_group
{

public:
  static Square_root_2_field factor;    // The multiplicative factor present in the general formula: sqrt(2) - 1

  typedef Hyperbolic_octagon_group<Square_root_2_field>   Self;  
  typedef complex<Square_root_2_field>                    Matrix_element;

  Matrix_element  A;
  Matrix_element  B;
  string          label;  

  Hyperbolic_octagon_group(const Matrix_element& A_, const Matrix_element& B_, 
    const string& label_ = string("") ) :
  A(A_), B(B_), label(label_) {}

  string merge_labels(const Self& rh) const
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
    string inv_label;
    for(long i = label.size() - 1; i >= 0; i--) {
      if(label[i] >= 'a' && label[i] <= 'd') {
        inv_label.push_back(label[i]);
        inv_label.push_back('^');
        inv_label.push_back('-');
        inv_label.push_back('1');
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
  Hyperbolic_octagon_group rotate() const
  {
    Square_root_2_field B1 = real(B);
    Square_root_2_field B2 = imag(B);
    
    // sqrt(2)
    Square_root_2_field k = Square_root_2_field(0, 1);
    
    Square_root_2_field BB1 = (B1 - B2)*k;
    Square_root_2_field BB2 = (B1 + B2)*k;

    assert(BB2.l % 2 == 0 && BB2.r % 2 == 0);

    BB1.l = BB1.l/2;
    BB1.r = BB1.r/2;
    BB2.l = BB2.l/2;
    BB2.r = BB2.r/2;

    return Hyperbolic_octagon_group(A, Matrix_element(BB1, BB2)); 
  }

  Square_root_2_field trace() const 
  {
    return Square_root_2_field(2, 0)*real(A);
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
  
  static complex<double> toComplexDouble(Matrix_element M) //const
  {
    Square_root_2_field rl  = real(M);
    Square_root_2_field img = imag(M);
    
    return complex<double>(rl.l + sqrt(2.)*rl.r, img.l + sqrt(2.)*img.r);
  }

  pair<double, double> apply(double x, double y)
  {
    typedef complex<double> Cmpl;
    Cmpl Aa = toComplexDouble(A);
    Cmpl Bb = toComplexDouble(B);

    double ax = sqrt(factor.l + sqrt(2.)*factor.r);
    
    Cmpl z(x, y);
    Cmpl res = (Aa*z + ax*Bb)/(ax*(conj(Bb)*z) + conj(Aa));
    return pair<double, double>(real(res), imag(res)); 
  }

};



template<class Square_root_2_field>
Square_root_2_field Hyperbolic_octagon_group<Square_root_2_field>::factor = Square_root_2_field(-1, 1); /*Square_root_2_field(-1, 1);*/

// just to give an order(ing)
template<class Int>
bool operator < (const complex<Square_root_2_field<Int> >& lh,
  const complex<Square_root_2_field<Int> >& rh)
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
template<class Square_root_2_field>
bool operator < (const Hyperbolic_octagon_group<Square_root_2_field>& lh, 
  const Hyperbolic_octagon_group<Square_root_2_field>& rh)
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

template<class Square_root_2_field>
bool operator == (const Hyperbolic_octagon_group<Square_root_2_field>& lh, 
  const Hyperbolic_octagon_group<Square_root_2_field>& rh)
{
  return (lh.A == rh.A && lh.B == rh.B);
}

template<class Square_root_2_field>
ostream& operator<<(ostream& os, const Hyperbolic_octagon_group<Square_root_2_field>& m)
{
  os << m.A << " " << m.B;
  return os;
}

typedef long long                                   ll;
typedef Square_root_2_field<ll>                     Sqrt_field;
typedef Hyperbolic_octagon_group<Sqrt_field>        Octagon_group;
typedef pair<Octagon_group, int>         Octagon_group_with_index;
typedef Octagon_group::Matrix_element    Element;


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




void get_generators(vector<Octagon_group>& gens)
{
  // This is a in the matrix, equal to sqrt(2) + 1
  Element A = Element(Sqrt_field(1, 1), Sqrt_field(0, 0));

  // This vector holds all other elements, results of the exponentials for various k
  vector<Element> B(8, Element(Sqrt_field(0, 0), Sqrt_field(0, 0)));

  // Corrected ordering (initial ordering is present in backup file)
  
  /* This here produces the correct ordering, and identity is given by g[G_A]*g[G_B]*g[G_C]*g[G_D]*g[G_InvA]*g[G_InvB]*g[G_InvC]*g[G_InvD] */
  B[DIRECTION_A]     = A * Element(Sqrt_field( 0,  1), Sqrt_field( 0,  0));
  B[DIRECTION_B]     = A * Element(Sqrt_field(-1,  0), Sqrt_field(-1,  0));
  B[DIRECTION_C]     = A * Element(Sqrt_field( 0,  0), Sqrt_field( 0,  1));
  B[DIRECTION_D]     = A * Element(Sqrt_field( 1,  0), Sqrt_field(-1,  0));
  B[DIRECTION_A_BAR] = A * Element(Sqrt_field( 0, -1), Sqrt_field( 0,  0));
  B[DIRECTION_B_BAR] = A * Element(Sqrt_field( 1,  0), Sqrt_field( 1,  0));
  B[DIRECTION_C_BAR] = A * Element(Sqrt_field( 0,  0), Sqrt_field( 0, -1));
  B[DIRECTION_D_BAR] = A * Element(Sqrt_field(-1,  0), Sqrt_field( 1,  0));
  
  string labels[8] = {string("a"), string("\\bar{b}"), string("c"), string("\\bar{d}"),
  string("\\bar{a}"), string("b"), string("\\bar{c}"), string("d")};
  for(int i = 0; i < 8; i++) {
    gens.push_back(Octagon_group(A, B[i], labels[i]));
  }
}





/***************************************************************/



  // Check whether two elements of the group are inverse of each other
  bool Are_inverse(Octagon_group_with_index x, Octagon_group_with_index y) {
    int idx = x.second % 4;
    int idy = y.second % 4;
    bool r = ((idx == idy) && (x.second != y.second));
    return r;
  }



  // Recursively eliminate neighboring inverse elements present in the word
  void Simplify_word(vector<Octagon_group_with_index>& w) {
    vector<Octagon_group_with_index> t;
    bool reduced = false;
    int N = w.size();
    if (N > 0) {
      for (int i = 0; i < N-1; i++) {
        if (!Are_inverse(w[i], w[i+1])) {
          t.push_back(w[i]);
        } else {
          reduced = true;
          i++;
        }
      }
      if (!Are_inverse(w[N-2], w[N-1])) {
        t.push_back(w[N-1]);
      } else {
        reduced = true;
      }
    
      if (reduced) {
        Simplify_word(t);
      }

      w = t;
    }
  }


  // Checks whether y is the next element from x
  bool Is_next(Octagon_group_with_index x, Octagon_group_with_index y) {
    return (((y.second + 5) % 8) == x.second);
  }


  // Given a word, find the largest subsequence of consecutive elements it contains.
  // The sequence itself is placed in 'seq', while the index at which it starts is
  // the return argument of the function.
  int Longest_sequence(vector<Octagon_group_with_index>& seq, vector<Octagon_group_with_index> const w) {
    int start = 0; 
    int mstart = 0;
    int end = 1;
    int max = 1;
    int len = 1;
    vector<Octagon_group_with_index> tmp, mvec;
    tmp.push_back(w[0]);
    for (int i = 1; i < w.size(); i++) {
      if ( Is_next( w[i], w[i-1] ) ) {
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



  Octagon_group_with_index Invert(Octagon_group_with_index x) {
    Octagon_group y = x.first.inverse();
    int idx = (x.second + 4) % 8;
    return Octagon_group_with_index(y, idx);
  }


  // Given a word, construct its inverse
  void Invert(vector<Octagon_group_with_index>& w, vector<Octagon_group_with_index> const original) {
    w.clear();
    for (int i = original.size() - 1; i >= 0; i--) {
      w.push_back( Invert(original[i]) );
    }
  }


  // Computes and returns the next index in the identity element chain
  int Next_index(int idx) {
    return ( (idx + 5) % 8 );
  }


  // Given a sequence of consecutive indices, return the complementary set of consecutive indices in mod 8.
  // For instance, if start = 5 and end = 1, the output is the sequence 2, 3, 4
  void Get_complement_indices(vector<int>& v, int begin, int end) {
    vector<int> tmp;
    for (int i = Next_index(end); i != begin; i = Next_index(i)) {
      tmp.push_back(i);
    }

    for (int i = tmp.size() - 1; i >= 0; i--) {
      v.push_back(tmp[i]);
    }
  }


  int Inverse(int idx) {
    return ((idx + 4) % 8);
  }


  // Given the start and end indices of a sequence of consecutive elements, construct the complementary sequence of elements
  void Get_complement_word(vector<Octagon_group_with_index>& c, int start, int end) {
    vector<int> idx;
    vector<Octagon_group> gens;
    get_generators(gens);

    Get_complement_indices(idx, start, end);    

    for (int i = idx.size()-1; i >= 0; i--) {
      int ii = idx[i]; 
      c.push_back(Octagon_group_with_index(gens[ii], ii));
    }
  }


  // Given a word consisting of consecutive elements, construct its complementary word
  void Get_complement_word(vector<Octagon_group_with_index>& c, vector<Octagon_group_with_index> w) {
    if (w.size() > 0) {
      Get_complement_word(c, w[0].second, w[w.size()-1].second);
    } else {
      c = w;
    }
  }




  // Given a word, identifies the longest subword consisting of consecutive elements and substitutes 
  // it with its equivalent shorter word. The search is executed in both the original word and its
  // inverse, and the substitution is made considering the longest subword from both cases.
  bool Replace_subword(vector<Octagon_group_with_index>& w, vector<Octagon_group_with_index> original) {

    bool replaced = false;

    // Handle empty string case
    if (original.size() == 0) {
      return replaced; // false
    }

    // Look for longest subword forward
    vector<Octagon_group_with_index> lfwd;
    int idxf = Longest_sequence(lfwd, original);
    int Nf = lfwd.size();

    // Get inverse of the original word
    vector<Octagon_group_with_index> inv;
    Invert(inv, original);

    // Look for longest subword backwards
    vector<Octagon_group_with_index> lbwd;
    int idxb = Longest_sequence(lbwd, inv);
    int Nb = lbwd.size();

    // Assign parameters based on results to homogenise the logic
    vector<Octagon_group_with_index> word, sub;
    bool is_inverse;
    int N, idx;
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
    // is by default equal to the identity element and can be directly eliminated.
    while (N >= 8) {
      replaced = true;
      vector<Octagon_group_with_index> ttt = word;
      word.clear();
      for (int i = 0; i < idx; i++) {
        word.push_back(ttt[i]);
      }
      for (int i = idx + 8; i < ttt.size(); i++) {
        word.push_back(ttt[i]);
      }
      w = word;

      ttt = sub;
      sub.clear();
      for (int i = 8; i < ttt.size(); i++) {
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
      vector<Octagon_group_with_index> cmpl;
      Get_complement_word(cmpl, sub);

      vector<Octagon_group_with_index> tmp;
      Invert(tmp, cmpl);
      cmpl = tmp;

      for (int i = 0; i < idx; i++) {
        w.push_back(word[i]);
      }
      for (int i = 0; i < cmpl.size(); i++) {
        w.push_back(cmpl[i]);
      }
      for (int i = N + idx; i < word.size(); i++) {
        w.push_back(word[i]);
      }
    }

    // If we have been working with the inverse of the original, the result has to be inverted.
    if (is_inverse) {
      vector<Octagon_group_with_index> tmp;
      Invert(tmp, w);
      w = tmp;
    }

    return replaced;
  }



  // Applies Dehn's algorithm to a given word. The result is the equivalent ireducible word
  // of the original. The boolean return argument of the function indicates whether the
  // resulting equivalent irreducible word is the identity element or not.
  bool Apply_Dehn(vector<Octagon_group_with_index>& w, vector<Octagon_group_with_index> const original) {
    bool is_identity = false;
    vector<Octagon_group_with_index> tmp;
    tmp = original;

    while (tmp.size() > 0) {
      Simplify_word(tmp);
      vector<Octagon_group_with_index> sub;
      bool replaced = Replace_subword(sub, tmp);
      if (!replaced) {
        w = tmp;
        return false;
      }
      tmp = sub;
    }

    w = tmp;

    return true;
  }



/**************************************************************/



// a, \bar{b}, c, \bar{d}, \bar{a}, b, \bar{c}, d
vector<Octagon_group> gens;

bool IsCanonical(const Octagon_group& m);

void generate_words( set<Octagon_group>& words, vector<Octagon_group>& prev, int depth, double threshold  )
{
  if (depth == 1) {
    for(int i = 0; i < 8; i++) {
      words.insert( gens[i] );
      prev.push_back( gens[i] );
    }
    return;
  }

  vector<Octagon_group> els;
  generate_words( words, els, depth - 1, threshold);
  
  Octagon_group temp = Octagon_group(Element(), Element());
  ll size = els.size();
  bool is_new = false;
  for(ll k = 0; k < size; k++) {
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
bool IsCanonical(const Octagon_group& m)
{
  Octagon_group temp = m;
  
  // rotate while |B1| < |B2|
  Sqrt_field B1, B2;
  Sqrt_field C = Sqrt_field(-1, -1);
  for(int i = 0; i < 8 && C != C.abs(); i++) {
    B1 = real(temp.B).abs();
    B2 = imag(temp.B).abs();
    C = B1 - B2;
    
    temp = temp.rotate();
  }
  assert(C == C.abs());

  // (2 - sqrt(2))(|B1| + (sqrt(2) -  1)|B2|)
  Sqrt_field right = Sqrt_field(2, -1)*(B1 + Sqrt_field(-1, 1)*B2);
  
  // |A2|
  Sqrt_field left = imag(temp.A).abs();

  // left <= right -> true
  C = right - left;
  return C == C.abs();
}

void dfs(const Octagon_group& m, set<Octagon_group>& visited)
{
  assert(IsCanonical(m));
  visited.insert(m);

  Octagon_group candidate = m;
  for(int i = 0; i < 8; i++) {
    candidate = gens[i]*m*gens[(i + 4) % 8];
    if(IsCanonical(candidate) == true && visited.find(candidate) == visited.end()) {
      dfs(candidate, visited);
    } 
  }
}

// map<Octagon_group m, Octagon_group Aux>, m = Aux * origin * Aux^{-1} 
void dfs_with_info(const pair<Octagon_group, Octagon_group>& new_pair, 
 map <Octagon_group, Octagon_group>& visited)
{
  assert(IsCanonical(new_pair.first));  
  visited.insert(new_pair);
  
  const Octagon_group& current = new_pair.first;
  const Octagon_group& current_factor = new_pair.second;
  Octagon_group candidate = current, candidate_factor = current_factor;
  for(int i = 0; i < 8; i++) {
    candidate = gens[i]*current*gens[(i + 4) % 8];
    if(IsCanonical(candidate) == true && visited.find(candidate) == visited.end()) {
      candidate_factor = gens[i]*current_factor;
      dfs_with_info(pair<Octagon_group, Octagon_group>(candidate, candidate_factor), visited);
    } 
  }
}


void dfs_with_info(const               Octagon_group& origin,
 map<Octagon_group, Octagon_group>& visited)
{
  assert(IsCanonical(origin));
  Octagon_group id = Octagon_group(Element(Sqrt_field(1, 0), Sqrt_field(0, 0)), 
   Element(Sqrt_field(0, 0), Sqrt_field(0, 0)));
  pair<Octagon_group, Octagon_group> new_pair(origin, id);
  
  dfs_with_info(new_pair, visited);
}


class IntersectionNumber
{
public:

  struct Point
  {
    Point(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}

    double x, y;
  };
  
  Octagon_group m;

  IntersectionNumber(const Octagon_group& m_) : m(m_)
  {}
  
  long operator() () const
  {
    set<Octagon_group> visited;
    set<Octagon_group>::iterator it, it2;
    map<long, long> nb_map;
    map<long, long>::iterator mit;

    dfs(m, visited);

    set<pair< Octagon_group, Octagon_group> > common;
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
  
  
  void count_nb(const Octagon_group& m1, const Octagon_group& m2, set<pair< Octagon_group, Octagon_group> >& visited) const
  {
    typedef pair<Octagon_group, Octagon_group> matrix_pair;
    visited.insert(matrix_pair(m1, m2));

    Octagon_group c1 = m1, c2 = m2;
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
  bool haveIntersection(const Octagon_group& m1, const Octagon_group& m2) const
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
  
  void intersectWithInfinity(const Octagon_group& m, Point& p1, Point& p2) const
  {
    Element a = m.A, b = m.B, factor = m.factor;

    Element four = Element(Sqrt_field(4, 0), Sqrt_field(0, 0));
    Element two = Element(Sqrt_field(2, 0), Sqrt_field(0, 0));

    Element D = (a - conj(a))*(a - conj(a));
    D += four*b*conj(b)*factor;
    Element T1 = conj(a) - a;
    Element T2 = two*conj(b);

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

void Delete(const set<Octagon_group>& canonical_set, vector<Octagon_group>& output)
{
  set<Octagon_group> redundant;

  set<Octagon_group>::iterator it;
  for(it = canonical_set.begin(); it != canonical_set.end(); ++it) {
    if(redundant.find(*it) != redundant.end()) { 
      continue;
    }

    set<Octagon_group> visited;
    dfs(*it, visited);
    visited.erase(*it);
    
    redundant.insert(visited.begin(), visited.end());    
    output.push_back(*it);
  }
}

void generate_unique_words(vector<Octagon_group>& output, double threshold = 10, int word_length = 13)
{
  get_generators(gens);
  
  set<Octagon_group> unique_words;
  vector<Octagon_group> temp;
  generate_words(unique_words, temp, word_length, threshold);
  
  double l = 0;
  set<Octagon_group>::iterator uit;
  for(uit = unique_words.begin(); uit != unique_words.end(); ++uit) {
    l = uit->length();
    if(0. < l && l < threshold) {
      output.push_back( *uit );
    }
  }
  
  cout << "nb of unique words " << output.size() << endl;
}

// words that correspond to union of 1-cycles
void generate_words_union_1_cycles(vector<Octagon_group>& out)
{
  if(gens.size() == 0) {
    get_generators(gens);
  }
  Octagon_group f[4] = {gens[0], gens[5], gens[2], gens[7]};
  
  Octagon_group F[8] = {
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
        
        Octagon_group Tr = F[i]*F[k].inverse()*F[j];
        // check that Tr != Id
        if(Tr.length() > 2.) {
          out.push_back(Tr);
        }
      }
    }
  }
  cout << counter << endl;
  
}



#endif  // CGAL_HYPERBOLIC_OCTAGON_GROUP_H



