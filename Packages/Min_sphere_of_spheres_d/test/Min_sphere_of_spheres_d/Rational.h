// A wrapper around the GMP mpq_t rational number-type.
// You need to have GNU GMP installed in order to use this.

#ifndef RATIONAL_WRAPPER
#define RATIONAL_WRAPPER

#include <utility>
#include <iostream>
#include <string>
#include <gmp.h>

class Rational {
private:
  mpq_t r;
  
private:
  void create() {
    mpq_init(r);
  }
  
public:
  // Construct zero.
  Rational() {
    create();
  }
  
  Rational(const Rational& x) {
    create();
    *this = x;
  }
  
  Rational& operator=(const Rational& x) {
    mpq_set(r,x.r);
    return *this;
  }
  
  ~Rational() {
    mpq_clear(r);
  }
  
public:
  // Construct the rational u.
  Rational(int u) {
    create();
    *this = std::pair<int,int>(u,1);
  }
  
  // Construct the rational u/v.
  Rational(int u,int v) {
    create();
    *this = std::pair<int,int>(u,v);
    mpq_canonicalize(r);
  }
  
  // Construct the rational from a double.
  Rational(double d) {
    create();
    mpq_set_d(r,d);
  }
  
public: // conversion:
  double toDouble() const {
    return mpq_get_d(r);
  }
  
  std::string toString() const {
    const int Base = 10;
    const int size = mpz_sizeinbase (mpq_numref(r),Base) + 
      mpz_sizeinbase (mpq_denref(r),Base) + 3;
    char *str = new char[size];
    mpq_get_str(str,Base,r);
    std::string s(str);
    delete str;
    return s;
  }
  
public: // comparision:
  bool operator==(const Rational& x) const {
    return mpq_cmp(r,x.r) == 0;
  }
  
  bool operator!=(const Rational& x) const {
    return mpq_cmp(r,x.r) != 0;
  }
  
  bool operator<(const Rational& x) const {
    return mpq_cmp(r,x.r) < 0;
  }
  
  bool operator>(const Rational& x) const {
    return mpq_cmp(r,x.r) > 0;
  }
  
  bool operator<=(const Rational& x) const {
    return mpq_cmp(r,x.r) <= 0;
  }
  
  bool operator>=(const Rational& x) const {
    return mpq_cmp(r,x.r) >= 0;
  }
  
public: // arithmetic:
  friend Rational operator+(const Rational& x,const Rational& y);
  friend Rational operator-(const Rational& x,const Rational& y);
  friend Rational operator*(const Rational& x,const Rational& y);
  friend Rational operator/(const Rational& x,const Rational& y);
  
  Rational& operator+=(const Rational& x) {
    return *this = *this + x;
  }
  
  Rational& operator-=(const Rational& x) {
    return *this = *this - x;
  }
  
  Rational& operator*=(const Rational& x) {
    return *this = *this * x;
  }
  
  Rational& operator/=(const Rational& x) {
    assert(x != 0);
    return *this = *this / x;
  }
  
  Rational operator-() const {
    Rational z;
    mpq_neg(z.r,r);
    return z;
  }
  
public:
  Rational& operator=(std::pair<int,int> frac) {
    if (frac.second < 0) {
      frac.second = -frac.second;
      frac.first  = -frac.first;
    }
    mpq_set_si(r,frac.first,frac.second);
    return *this;
  }
};

Rational operator+(const Rational& x,const Rational& y) {
  Rational z;
  mpq_add(z.r,x.r,y.r);
  return z;
}

Rational operator-(const Rational& x,const Rational& y) {
  Rational z;
  mpq_sub(z.r,x.r,y.r);
  return z;
}

Rational operator*(const Rational& x,const Rational& y) {
  Rational z;
  mpq_mul(z.r,x.r,y.r);
  return z;
}

Rational operator/(const Rational& x,const Rational& y) {
  assert(y != 0);
  Rational z;
  mpq_div(z.r,x.r,y.r);
  return z;
}

std::ostream& operator<<(std::ostream& out,const Rational& r) {
  return out << r.toString();
}

double to_double(const Rational& x) {
  return x.toDouble();
}

namespace CGAL {
  namespace NTS {
    double to_double(const Rational& x) {
      return ::to_double(x);
    }
  }
}

#endif // RATIONAL_WRAPPER
