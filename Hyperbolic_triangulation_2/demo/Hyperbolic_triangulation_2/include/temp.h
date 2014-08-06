// to compile with boost headers
// g++ main.cpp -I /opt/local/include -o main

#include <complex>
#include <iostream>
#include <vector>
#include <iterator>
#include <map>
#include <set>
#include <assert.h>
#include <unordered_set>
#include <string>
#include <fstream>

//#include <boost/optional.hpp>

using namespace std;

template<class Int>
class Sqrt_field 
{
public:
  typedef Int Integer;

  Int l;
  Int r;

public:
  Sqrt_field(Int l_ = 0, Int r_ = 0) : l(l_), r(r_) {}
  
  Sqrt_field operator + (const Sqrt_field& rh) const
  {
    Sqrt_field temp = *this;
    temp += rh;

    return temp;
  } 

  Sqrt_field& operator += (const Sqrt_field& rh)
  {
    l += rh.l;
    r += rh.r;
    return *this;
  }
  
  Sqrt_field& operator -= (const Sqrt_field& rh)
  {
    l -= rh.l;
    r -= rh.r;
    return *this;
  }

  Sqrt_field operator - (const Sqrt_field& rh) const
  {
    return Sqrt_field( l - rh.l, r - rh.r );
  }

  Sqrt_field operator - () const
  {
    return Sqrt_field(-l, -r);
  }

  Sqrt_field operator * (const Sqrt_field& rh) const
  {
    Sqrt_field temp = *this;
    temp *= rh;

    return temp;
  }

  Sqrt_field& operator *= (const Sqrt_field& rh)
  {
    Int l1 = l * rh.l + 2*r*rh.r;
    Int r1 = r*rh.l + l*rh.r;
    l = l1;
    r = r1;

    return *this; 
  }

  Sqrt_field abs() const
  {
    if (l >= 0 && r >= 0) {
      return *this;
    }

    if (l <= 0 && r <= 0) {
      return Sqrt_field(-l, -r); 
    }
  
    if(l*l*l - 2*r*r*l >= 0) {
      return *this;
    }
    return Sqrt_field(-l, -r);
  }
};

template<class Int>
bool operator < (const Sqrt_field<Int>& lh, const Sqrt_field<Int>& rh)
{
  if( (lh.l - rh.l)*(lh.l - rh.l) > 2*(rh.r - lh.r)*(rh.r - lh.r) ) {
    if ((lh.l - rh.l) < 0) {
      return true;
    }
  }
  
  if( (lh.l - rh.l)*(lh.l - rh.l) < 2*(rh.r - lh.r)*(rh.r - lh.r) ) {
    if ((rh.r - lh.r) >= 0) {
      return true;
    }
  }
  return false;
  
  /*
  if (lh.l < rh.l) {
    return true;
  }
  
  if (lh.l == rh.l) {
    if (lh.r < rh.r) {
      return true;
    }
  }*/
  
  return false;
}
/*Seems correct!
template<class Int>
bool operator > (const Sqrt_field<Int>& lh, const Sqrt_field<Int>& rh)
{
  if( (lh.l - rh.l)*(lh.l - rh.l) > 2*(rh.r - lh.r)*(rh.r - lh.r) ) {
    if ((lh.l - rh.l) >= 0 && (rh.r - lh.r) >= 0) {
      return true;
    }
  }
  
  if( (lh.l - rh.l)*(lh.l - rh.l) < 2*(rh.r - lh.r)*(rh.r - lh.r) ) {
    if ((lh.l - rh.l) <= 0 && (rh.r - lh.r) <= 0) {
      return true;
    }
  }
  return false;
}*/

template<class Int>
bool operator >= (const Sqrt_field<Int>& lh, const Sqrt_field<Int>& rh)
{
  if( (lh.l - rh.l)*(lh.l - rh.l) >= 2*(rh.r - lh.r)*(rh.r - lh.r) ) {
    if ((lh.l - rh.l) >= 0) {
      return true;
    }
  }
  
  if( (lh.l - rh.l)*(lh.l - rh.l) < 2*(rh.r - lh.r)*(rh.r - lh.r) ) {
    if ((rh.r - lh.r) < 0) {
      return true;
    }
  }
  return false;
}

template<class Int>
bool operator == (const Sqrt_field<Int>& lh, const Sqrt_field<Int>& rh)
{
  return (lh.l == rh.l && lh.r == rh.r);
}

template<class Int>
bool operator != (const Sqrt_field<Int>& lh, const Sqrt_field<Int>& rh)
{
  return !(lh == rh);
}

template<class Int>
Sqrt_field<Int> operator * (const Int& val, const Sqrt_field<Int>& rh)
{
  Sqrt_field<Int> temp = rh;
  temp.l *= val;
  temp.r *= val;

  return temp;
} 

template<class Int>
ostream& operator<<(ostream& os, const Sqrt_field<Int>& nb)
{
  os << nb.l << " " << nb.r;
  return os;
}

template<class Sqrt_field>
class Octagon_matrix
{
public:
  typedef Octagon_matrix<Sqrt_field> Self;  
  typedef complex<Sqrt_field> Extended_field;

  Extended_field M11;
  Extended_field M12;
  string label;  

  Octagon_matrix(const Extended_field& M11_, const Extended_field& M12_, 
    const string& label_ = string("") ) :
    M11(M11_), M12(M12_), label(label_) {}

  string merge_labels(const Self& rh) const
  {
    return label + "*" + rh.label;
  }

  Self operator*(const Self& rh) const
  {
    return Self( M11*rh.M11 + aux*M12*conj(rh.M12), 
      M11*rh.M12 + M12*conj(rh.M11), merge_labels(rh) );
  }
  
  Self inverse() const
  {
    Self inv = Self(conj(M11), -M12);
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
  Octagon_matrix rotate() const
  {
    Sqrt_field B1 = real(M12);
    Sqrt_field B2 = imag(M12);
    
    // sqrt(2)
    Sqrt_field k = Sqrt_field(0, 1);
    
    Sqrt_field BB1 = (B1 - B2)*k;
    Sqrt_field BB2 = (B1 + B2)*k;

    assert(BB2.l % 2 == 0 && BB2.r % 2 == 0);
    BB1.l = BB1.l/2;
    BB1.r = BB1.r/2;
    BB2.l = BB2.l/2;
    BB2.r = BB2.r/2;

    return Octagon_matrix(M11, Extended_field(BB1, BB2)); 
  }

  Sqrt_field trace() const 
  {
    return Sqrt_field(2, 0)*real(M11);
  }

  double length() const
  {
    typedef long double ld;

    ld l = real(M11).l;
    ld r = real(M11).r;
    ld tr = l + sqrt(2.)*r;
    if (tr < 0) {
      tr = -tr;
    }

    return 2.*acosh(tr);
  }

  // determinant == 1
  Extended_field det() const
  {
    return norm(M11) - aux * norm(M12);
  }
  
  static complex<double> toComplexDouble(Extended_field M) //const
  {
    Sqrt_field rl = real(M);
    Sqrt_field img = imag(M);
    
    return complex<double>(rl.l + sqrt(2.)*rl.r, img.l + sqrt(2.)*img.r);
  }

  pair<double, double> apply(double x, double y)
  {
    typedef complex<double> Cmpl;
    Cmpl m11 = toComplexDouble(M11);
    Cmpl m12 = toComplexDouble(M12);

    double ax = sqrt(aux.l + sqrt(2.)*aux.r);
    
    Cmpl z(x, y);
    Cmpl res = (m11*z + ax*m12)/(ax*(conj(m12)*z) + conj(m11));
    return pair<double, double>(real(res), imag(res)); 
  }

//private:
  static Sqrt_field aux;
};


template<class Sqrt_field>
Sqrt_field Octagon_matrix<Sqrt_field>::aux = Sqrt_field(-1, 1);

// just to give an order
template<class Int>
bool operator < (const complex<Sqrt_field<Int> >& lh,
  const complex<Sqrt_field<Int> >& rh)
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
template<class Sqrt_field>
bool operator < (const Octagon_matrix<Sqrt_field>& lh, 
  const Octagon_matrix<Sqrt_field>& rh)
{
  if (lh.M11 < rh.M11) {
    return true;
  }
  
  if (lh.M11 == rh.M11 ) {
    if (lh.M12 < rh.M12) {
      return true;
    }
  }

  return false;
}

template<class Sqrt_field>
bool operator == (const Octagon_matrix<Sqrt_field>& lh, 
  const Octagon_matrix<Sqrt_field>& rh)
{
  return (lh.M11 == rh.M11 && lh.M12 == rh.M12);
}

template<class Sqrt_field>
ostream& operator<<(ostream& os, const Octagon_matrix<Sqrt_field>& m)
{
  os << m.M11 << " " << m.M12;
  return os;
}

typedef long long ll;
typedef Sqrt_field<ll> SqrtField;
typedef Octagon_matrix<SqrtField> OctagonMatrix;
typedef OctagonMatrix::Extended_field Entry;

void get_generators(vector<OctagonMatrix>& gens)
{
  Entry M11 = Entry(SqrtField(1, 1), SqrtField(0, 0));

  vector<Entry> M12(8, Entry(SqrtField(0, 0), SqrtField(0, 0)));
  M12[0] = M11 * Entry(SqrtField(0, 1), SqrtField(0, 0));
  M12[1] = M11 * Entry(SqrtField(1, 0), SqrtField(1, 0));
  M12[2] = M11 * Entry(SqrtField(0, 0), SqrtField(0, 1));
  M12[3] = M11 * Entry(SqrtField(-1, 0), SqrtField(1, 0));
  M12[4] = M11 * Entry(SqrtField(0, -1), SqrtField(0, 0));
  M12[5] = M11 * Entry(SqrtField(-1, 0), SqrtField(-1, 0));
  M12[6] = M11 * Entry(SqrtField(0, 0), SqrtField(0, -1));
  M12[7] = M11 * Entry(SqrtField(1, 0), SqrtField(-1, 0));

  string labels[8] = {string("a"), string("b^-1"), string("c"), string("d^-1"),
    string("a^-1"), string("b"), string("c^-1"), string("d")};
  for(int i = 0; i < 8; i++) {
    gens.push_back(OctagonMatrix(M11, M12[i], labels[i]));
  }
}

// a, b, c, d, a^-1, b^-1, c^-1, d^-1
vector<OctagonMatrix> gens;

bool IsCanonical(const OctagonMatrix& m);

void generate_words( set<OctagonMatrix>& words, vector<OctagonMatrix>& prev, int depth  )
{
  if (depth == 1) {
    for(int i = 0; i < 8; i++) {
      words.insert( gens[i] );
      prev.push_back( gens[i] );
    }
    return;
  }

  vector<OctagonMatrix> els;
  generate_words( words, els, depth - 1);
  
  OctagonMatrix temp = OctagonMatrix(Entry(), Entry());
  ll size = els.size();
  bool is_new = false;
  for(ll k = 0; k < size; k++) {
    for(int i = 0; i < 8; i++) {
      temp = els[k]*gens[i];

      if(temp.length() > 15.) {
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
bool IsCanonical(const OctagonMatrix& m)
{
  OctagonMatrix temp = m;
  
  // rotate while |B1| < |B2|
  SqrtField B1, B2;
  SqrtField C = SqrtField(-1, -1);
  for(int i = 0; i < 8 && C != C.abs(); i++) {
    B1 = real(temp.M12).abs();
    B2 = imag(temp.M12).abs();
    C = B1 - B2;
    
    temp = temp.rotate();
  }
  assert(C == C.abs());

  // (2 - sqrt(2))(|B1| + (sqrt(2) -  1)|B2|)
  SqrtField right = SqrtField(2, -1)*(B1 + SqrtField(-1, 1)*B2);
  
  // |A2|
  SqrtField left = imag(temp.M11).abs();
 
  // left <= right -> true
  C = right - left;
  return C == C.abs();
}

void dfs(const OctagonMatrix& m, set<OctagonMatrix>& visited)
{
  assert(IsCanonical(m));
  visited.insert(m);

  OctagonMatrix candidate = m;
  for(int i = 0; i < 8; i++) {
    candidate = gens[i]*m*gens[(i + 4) % 8];
    if(IsCanonical(candidate) == true && visited.find(candidate) == visited.end()) {
      dfs(candidate, visited);
    } 
  }
}

// map<OctagonMatrix m, OctagonMatrix Aux>, m = Aux * origin * Aux^{-1} 
void dfs_with_info(const pair<OctagonMatrix, OctagonMatrix>& new_pair, 
  map<OctagonMatrix, OctagonMatrix>& visited)
{
  assert(IsCanonical(new_pair.first));  
  visited.insert(new_pair);
  
  const OctagonMatrix& current = new_pair.first;
  const OctagonMatrix& current_aux = new_pair.second;
  OctagonMatrix candidate = current, candidate_aux = current_aux;
  for(int i = 0; i < 8; i++) {
    candidate = gens[i]*current*gens[(i + 4) % 8];
    if(IsCanonical(candidate) == true && visited.find(candidate) == visited.end()) {
      candidate_aux = gens[i]*current_aux;
      dfs_with_info(pair<OctagonMatrix, OctagonMatrix>(candidate, candidate_aux), visited);
    } 
  }
}


void dfs_with_info(const OctagonMatrix& origin, map<OctagonMatrix, OctagonMatrix>& visited)
{
  assert(IsCanonical(origin));
  OctagonMatrix id = OctagonMatrix(Entry(SqrtField(1, 0), SqrtField(0, 0)), 
                                   Entry(SqrtField(0, 0), SqrtField(0, 0)));
  pair<OctagonMatrix, OctagonMatrix> new_pair(origin, id);
  
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
  
OctagonMatrix m;

IntersectionNumber(const OctagonMatrix& m_) : m(m_)
{}
  
long operator() () const
{
  set<OctagonMatrix> visited;
  set<OctagonMatrix>::iterator it, it2;
  map<long, long> nb_map;
  map<long, long>::iterator mit;
  
  dfs(m, visited);
  
  set<pair< OctagonMatrix, OctagonMatrix> > common;
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
  
  
void count_nb(const OctagonMatrix& m1, const OctagonMatrix& m2, set<pair< OctagonMatrix, OctagonMatrix> >& visited) const
{
  typedef pair<OctagonMatrix, OctagonMatrix> matrix_pair;
  visited.insert(matrix_pair(m1, m2));
  
  OctagonMatrix c1 = m1, c2 = m2;
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
bool haveIntersection(const OctagonMatrix& m1, const OctagonMatrix& m2) const
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
  
void intersectWithInfinity(const OctagonMatrix& m, Point& p1, Point& p2) const
{
  Entry a = m.M11, b = m.M12, aux = m.aux;
  
  Entry four = Entry(SqrtField(4, 0), SqrtField(0, 0));
  Entry two = Entry(SqrtField(2, 0), SqrtField(0, 0));
  
  Entry D = (a - conj(a))*(a - conj(a));
  D += four*b*conj(b)*aux;
  Entry T1 = conj(a) - a;
  Entry T2 = two*conj(b);
  
  complex<double> d = m.toComplexDouble(D);
  complex<double> t1 = m.toComplexDouble(T1);
  complex<double> t2 = m.toComplexDouble(T2);
  complex<double> au = complex<double>(m.aux.l + sqrt(2.)*m.aux.r, 0);
  
  complex<double> z1 = (t1 + sqrt(d))/(t2*sqrt(au));
  complex<double> z2 = (t1 - sqrt(d))/(t2*sqrt(au));
  
  p1 = Point(real(z1), imag(z1));
  p2 = Point(real(z2), imag(z2));
  
  assert(p1.x*p1.x + p1.y*p1.y > 0.99 && p1.x*p1.x + p1.y*p1.y < 1.01);
  assert(p2.x*p2.x + p2.y*p2.y > 0.99 && p2.x*p2.x + p2.y*p2.y < 1.01);
}

  
};

void Delete(const set<OctagonMatrix>& canonical_set, vector<OctagonMatrix>& output)
{
  set<OctagonMatrix> redundant;

  set<OctagonMatrix>::iterator it;
  for(it = canonical_set.begin(); it != canonical_set.end(); ++it) {
    if(redundant.find(*it) != redundant.end()) { 
      continue;
    }
  
    set<OctagonMatrix> visited;
    dfs(*it, visited);
    visited.erase(*it);
    
    redundant.insert(visited.begin(), visited.end());    
    output.push_back(*it);
  }
}

void generate_unique_words(vector<OctagonMatrix>& output, double threshold = 10, int word_length = 13)
{
  get_generators(gens);
  
  set<OctagonMatrix> unique_words;
  vector<OctagonMatrix> temp;
  generate_words(unique_words, temp, word_length);
  
  double l = 0;
  set<OctagonMatrix>::iterator uit;
  for(uit = unique_words.begin(); uit != unique_words.end(); ++uit) {
    l = uit->length();
    if(0. < l && l < threshold) {
      output.push_back( *uit );
    }
  }
  
  cout << "nb of unique words " << output.size() << endl;
}

// words that correspond to union of 1-cycles
void generate_words_union_1_cycles(vector<OctagonMatrix>& out)
{
  if(gens.size() == 0) {
    get_generators(gens);
  }
  OctagonMatrix f[4] = {gens[0], gens[5], gens[2], gens[7]};
  
  OctagonMatrix F[8] = {
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
        
        OctagonMatrix Tr = F[i]*F[k].inverse()*F[j];
        // check that Tr != Id
        if(Tr.length() > 2.) {
          out.push_back(Tr);
        }
      }
    }
  }
  cout << counter << endl;
  
}

