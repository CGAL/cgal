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
Square_root_2_field Hyperbolic_octagon_group<Square_root_2_field>::factor = Square_root_2_field(-1, 1);

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

typedef long long                                 ll;
typedef Square_root_2_field<ll>                   SqrtField;
typedef Hyperbolic_octagon_group<SqrtField>       HyperbolicOctagonGroup;
typedef HyperbolicOctagonGroup::Matrix_element    Element;


enum Direction {
  G_A = 0,  // 0
  G_InvB,   // 1
  G_C,      // 2
  G_InvD,   // 3
  G_InvA,   // 4
  G_B,      // 5
  G_InvC,   // 6
  G_D       // 7
};


void get_generators(vector<HyperbolicOctagonGroup>& gens)
{
  // This is a in the matrix, equal to sqrt(2) + 1
  Element A = Element(SqrtField(1, 1), SqrtField(0, 0));

  // This vector holds all other elements, results of the exponentials for various k
  vector<Element> B(8, Element(SqrtField(0, 0), SqrtField(0, 0)));

  // Corrected ordering (initial ordering is present in backup file)
  B[G_A]    = A * Element(SqrtField( 0,  0), SqrtField( 0, -1));
  B[G_B] = A * Element(SqrtField( 1,  0), SqrtField(-1,  0));
  B[G_C]    = A * Element(SqrtField( 0,  1), SqrtField( 0,  0));
  B[G_D] = A * Element(SqrtField( 1,  0), SqrtField( 1,  0));
  B[G_InvA] = A * Element(SqrtField( 0,  0), SqrtField( 0,  1));
  B[G_InvB]    = A * Element(SqrtField(-1,  0), SqrtField( 1,  0));
  B[G_InvC] = A * Element(SqrtField( 0, -1), SqrtField( 0,  0));
  B[G_InvD]    = A * Element(SqrtField(-1,  0), SqrtField(-1,  0));
  
  

  string labels[8] = {string("a"), string("\\bar{b}"), string("c"), string("\\bar{d}"),
  string("\\bar{a}"), string("b"), string("\\bar{c}"), string("d")};
  for(int i = 0; i < 8; i++) {
    gens.push_back(HyperbolicOctagonGroup(A, B[i], labels[i]));
  }
}

// a, \bar{b}, c, \bar{d}, \bar{a}, b, \bar{c}, d
vector<HyperbolicOctagonGroup> gens;

bool IsCanonical(const HyperbolicOctagonGroup& m);

void generate_words( set<HyperbolicOctagonGroup>& words, vector<HyperbolicOctagonGroup>& prev, int depth, double threshold  )
{
  if (depth == 1) {
    for(int i = 0; i < 8; i++) {
      words.insert( gens[i] );
      prev.push_back( gens[i] );
    }
    return;
  }

  vector<HyperbolicOctagonGroup> els;
  generate_words( words, els, depth - 1, threshold);
  
  HyperbolicOctagonGroup temp = HyperbolicOctagonGroup(Element(), Element());
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
bool IsCanonical(const HyperbolicOctagonGroup& m)
{
  HyperbolicOctagonGroup temp = m;
  
  // rotate while |B1| < |B2|
  SqrtField B1, B2;
  SqrtField C = SqrtField(-1, -1);
  for(int i = 0; i < 8 && C != C.abs(); i++) {
    B1 = real(temp.B).abs();
    B2 = imag(temp.B).abs();
    C = B1 - B2;
    
    temp = temp.rotate();
  }
  assert(C == C.abs());

  // (2 - sqrt(2))(|B1| + (sqrt(2) -  1)|B2|)
  SqrtField right = SqrtField(2, -1)*(B1 + SqrtField(-1, 1)*B2);
  
  // |A2|
  SqrtField left = imag(temp.A).abs();

  // left <= right -> true
  C = right - left;
  return C == C.abs();
}

void dfs(const HyperbolicOctagonGroup& m, set<HyperbolicOctagonGroup>& visited)
{
  assert(IsCanonical(m));
  visited.insert(m);

  HyperbolicOctagonGroup candidate = m;
  for(int i = 0; i < 8; i++) {
    candidate = gens[i]*m*gens[(i + 4) % 8];
    if(IsCanonical(candidate) == true && visited.find(candidate) == visited.end()) {
      dfs(candidate, visited);
    } 
  }
}

// map<HyperbolicOctagonGroup m, HyperbolicOctagonGroup Aux>, m = Aux * origin * Aux^{-1} 
void dfs_with_info(const pair<HyperbolicOctagonGroup, HyperbolicOctagonGroup>& new_pair, 
 map <HyperbolicOctagonGroup, HyperbolicOctagonGroup>& visited)
{
  assert(IsCanonical(new_pair.first));  
  visited.insert(new_pair);
  
  const HyperbolicOctagonGroup& current = new_pair.first;
  const HyperbolicOctagonGroup& current_factor = new_pair.second;
  HyperbolicOctagonGroup candidate = current, candidate_factor = current_factor;
  for(int i = 0; i < 8; i++) {
    candidate = gens[i]*current*gens[(i + 4) % 8];
    if(IsCanonical(candidate) == true && visited.find(candidate) == visited.end()) {
      candidate_factor = gens[i]*current_factor;
      dfs_with_info(pair<HyperbolicOctagonGroup, HyperbolicOctagonGroup>(candidate, candidate_factor), visited);
    } 
  }
}


void dfs_with_info(const               HyperbolicOctagonGroup& origin,
 map<HyperbolicOctagonGroup, HyperbolicOctagonGroup>& visited)
{
  assert(IsCanonical(origin));
  HyperbolicOctagonGroup id = HyperbolicOctagonGroup(Element(SqrtField(1, 0), SqrtField(0, 0)), 
   Element(SqrtField(0, 0), SqrtField(0, 0)));
  pair<HyperbolicOctagonGroup, HyperbolicOctagonGroup> new_pair(origin, id);
  
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
  
  HyperbolicOctagonGroup m;

  IntersectionNumber(const HyperbolicOctagonGroup& m_) : m(m_)
  {}
  
  long operator() () const
  {
    set<HyperbolicOctagonGroup> visited;
    set<HyperbolicOctagonGroup>::iterator it, it2;
    map<long, long> nb_map;
    map<long, long>::iterator mit;

    dfs(m, visited);

    set<pair< HyperbolicOctagonGroup, HyperbolicOctagonGroup> > common;
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
  
  
  void count_nb(const HyperbolicOctagonGroup& m1, const HyperbolicOctagonGroup& m2, set<pair< HyperbolicOctagonGroup, HyperbolicOctagonGroup> >& visited) const
  {
    typedef pair<HyperbolicOctagonGroup, HyperbolicOctagonGroup> matrix_pair;
    visited.insert(matrix_pair(m1, m2));

    HyperbolicOctagonGroup c1 = m1, c2 = m2;
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
  bool haveIntersection(const HyperbolicOctagonGroup& m1, const HyperbolicOctagonGroup& m2) const
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
  
  void intersectWithInfinity(const HyperbolicOctagonGroup& m, Point& p1, Point& p2) const
  {
    Element a = m.A, b = m.B, factor = m.factor;

    Element four = Element(SqrtField(4, 0), SqrtField(0, 0));
    Element two = Element(SqrtField(2, 0), SqrtField(0, 0));

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

void Delete(const set<HyperbolicOctagonGroup>& canonical_set, vector<HyperbolicOctagonGroup>& output)
{
  set<HyperbolicOctagonGroup> redundant;

  set<HyperbolicOctagonGroup>::iterator it;
  for(it = canonical_set.begin(); it != canonical_set.end(); ++it) {
    if(redundant.find(*it) != redundant.end()) { 
      continue;
    }

    set<HyperbolicOctagonGroup> visited;
    dfs(*it, visited);
    visited.erase(*it);
    
    redundant.insert(visited.begin(), visited.end());    
    output.push_back(*it);
  }
}

void generate_unique_words(vector<HyperbolicOctagonGroup>& output, double threshold = 10, int word_length = 13)
{
  get_generators(gens);
  
  set<HyperbolicOctagonGroup> unique_words;
  vector<HyperbolicOctagonGroup> temp;
  generate_words(unique_words, temp, word_length, threshold);
  
  double l = 0;
  set<HyperbolicOctagonGroup>::iterator uit;
  for(uit = unique_words.begin(); uit != unique_words.end(); ++uit) {
    l = uit->length();
    if(0. < l && l < threshold) {
      output.push_back( *uit );
    }
  }
  
  cout << "nb of unique words " << output.size() << endl;
}

// words that correspond to union of 1-cycles
void generate_words_union_1_cycles(vector<HyperbolicOctagonGroup>& out)
{
  if(gens.size() == 0) {
    get_generators(gens);
  }
  HyperbolicOctagonGroup f[4] = {gens[0], gens[5], gens[2], gens[7]};
  
  HyperbolicOctagonGroup F[8] = {
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
        
        HyperbolicOctagonGroup Tr = F[i]*F[k].inverse()*F[j];
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



