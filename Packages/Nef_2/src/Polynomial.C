// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : src/Polynomial.C
// package       : Nef_2
// author(s)     : Michael Seel
//                 Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Susan Heart
//
// ============================================================================



#include <CGAL/Nef_2/Polynomial.h>

namespace CGAL{

#if defined(_MSC_VER)
Polynomial<int> operator - (const Polynomial<int>& p)
{
  CGAL_assertion(p.degree()>=0);
  Polynomial<int> res(p.coeffs().begin(),p.coeffs().end());
  Polynomial<int>::iterator it, ite=res.coeffs().end();
  for(it=res.coeffs().begin(); it!=ite; ++it) *it = -*it;
  return res;
}

Polynomial<double> operator - (const Polynomial<double>& p)
{
  CGAL_assertion(p.degree()>=0);
  Polynomial<double> res(p.coeffs().begin(),p.coeffs().end());
  Polynomial<double>::iterator it, ite=res.coeffs().end();
  for(it=res.coeffs().begin(); it!=ite; ++it) *it = -*it;
  return res;
}
#endif

// SPECIALIZE_IMPLEMENTATION(NT,int double) START
// SPECIALIZING to :

void Polynomial<int>::euclidean_div(
  const Polynomial<int>& f, const Polynomial<int>& g,
  Polynomial<int>& q, Polynomial<int>& r)
{
  r = f; r.copy_on_write();
  int rd=r.degree(), gd=g.degree(), qd;
  if ( rd < gd ) { q = Polynomial<int>(int(0)); }
  else { qd = rd-gd+1; q = Polynomial<int>(std::size_t(qd)); }
  while ( rd >= gd ) {
    int S = r[rd] / g[gd];
    qd = rd-gd;
    q.coeff(qd) += S;
    r.minus_offsetmult(g,S,qd);
    rd = r.degree();
  }
  CGAL_postcondition( f==q*g+r );
}



void Polynomial<int>::pseudo_div(
  const Polynomial<int>& f, const Polynomial<int>& g, 
  Polynomial<int>& q, Polynomial<int>& r, int& D)
{
  TRACEN("pseudo_div "<<f<<" , "<< g);
  int fd=f.degree(), gd=g.degree();
  if ( fd<gd ) 
  { q = Polynomial<int>(0); r = f; D = 1; 
    CGAL_postcondition(Polynomial<int>(D)*f==q*g+r); return; 
  }
  // now we know fd >= gd and f>=g
  int qd=fd-gd, delta=qd+1, rd=fd;
  { q = Polynomial<int>( std::size_t(delta) ); }; // workaround for SUNPRO
  int G = g[gd]; // highest order coeff of g
  D = G; while (--delta) D*=G; // D = G^delta
  Polynomial<int> res = Polynomial<int>(D)*f;
  TRACEN("  pseudo_div start "<<res<<" "<<qd<<" "<<q.degree());
  while (qd >= 0) {
    int F = res[rd]; // highest order coeff of res
    int t = F/G;     // ensured to be integer by multiplication of D
    q.coeff(qd) = t;    // store q coeff
    res.minus_offsetmult(g,t,qd); 
    if (res.is_zero()) break;
    rd = res.degree();
    qd = rd - gd;
  }
  r = res;
  CGAL_postcondition(Polynomial<int>(D)*f==q*g+r);
  TRACEN("  returning "<<q<<", "<<r<<", "<< D);
}



Polynomial<int> Polynomial<int>::gcd(
  const Polynomial<int>& p1, const Polynomial<int>& p2)
{ TRACEN("gcd("<<p1<<" , "<<p2<<")");
  if ( p1.is_zero() )
    if ( p2.is_zero() ) return Polynomial<int>(int(1));
    else return p2.abs();
  if ( p2.is_zero() )
    return p1.abs();

  Polynomial<int> f1 = p1.abs();
  Polynomial<int> f2 = p2.abs();
  int f1c = f1.content(), f2c = f2.content();
  f1 /= f1c; f2 /= f2c;
  int F = CGAL_NTS gcd(f1c,f2c);
  Polynomial<int> q,r; int M=1,D;
  bool first = true;
  while ( ! f2.is_zero() ) { 
    Polynomial<int>::pseudo_div(f1,f2,q,r,D);
    if (!first) M*=D;
    TRACEV(f1);TRACEV(f2);TRACEV(q);TRACEV(r);TRACEV(M);
    r /= r.content();
    f1=f2; f2=r;
    first=false;
  }
  TRACEV(f1.content());
  return Polynomial<int>(F)*f1.abs();
}


// SPECIALIZING to :

void Polynomial<double>::euclidean_div(
  const Polynomial<double>& f, const Polynomial<double>& g,
  Polynomial<double>& q, Polynomial<double>& r)
{
  r = f; r.copy_on_write();
  int rd=r.degree(), gd=g.degree(), qd;
  if ( rd < gd ) { q = Polynomial<double>(double(0)); }
  else { qd = rd-gd+1; q = Polynomial<double>(std::size_t(qd)); }
  while ( rd >= gd ) {
    double S = r[rd] / g[gd];
    qd = rd-gd;
    q.coeff(qd) += S;
    r.minus_offsetmult(g,S,qd);
    rd = r.degree();
  }
  CGAL_postcondition( f==q*g+r );
}



void Polynomial<double>::pseudo_div(
  const Polynomial<double>& f, const Polynomial<double>& g, 
  Polynomial<double>& q, Polynomial<double>& r, double& D)
{
  TRACEN("pseudo_div "<<f<<" , "<< g);
  int fd=f.degree(), gd=g.degree();
  if ( fd<gd ) 
  { q = Polynomial<double>(0); r = f; D = 1; 
    CGAL_postcondition(Polynomial<double>(D)*f==q*g+r); return; 
  }
  // now we know fd >= gd and f>=g
  int qd=fd-gd, delta=qd+1, rd=fd;
  q = Polynomial<double>( std::size_t(delta) );
  double G = g[gd]; // highest order coeff of g
  D = G; while (--delta) D*=G; // D = G^delta
  Polynomial<double> res = Polynomial<double>(D)*f;
  TRACEN("  pseudo_div start "<<res<<" "<<qd<<" "<<q.degree());
  while (qd >= 0) {
    double F = res[rd]; // highest order coeff of res
    double t = F/G;     // ensured to be integer by multiplication of D
    q.coeff(qd) = t;    // store q coeff
    res.minus_offsetmult(g,t,qd); 
    if (res.is_zero()) break;
    rd = res.degree();
    qd = rd - gd;
  }
  r = res;
  CGAL_postcondition(Polynomial<double>(D)*f==q*g+r);
  TRACEN("  returning "<<q<<", "<<r<<", "<< D);
}

Polynomial<double> Polynomial<double>::gcd(
  const Polynomial<double>& p1, const Polynomial<double>& p2)
{ TRACEN("gcd("<<p1<<" , "<<p2<<")");
  if ( p1.is_zero() )
    if ( p2.is_zero() ) return Polynomial<double>(double(1));
    else return p2.abs();
  if ( p2.is_zero() )
    return p1.abs();

  Polynomial<double> f1 = p1.abs();
  Polynomial<double> f2 = p2.abs();
  double f1c = f1.content(), f2c = f2.content();
  f1 /= f1c; f2 /= f2c;
  Polynomial<double> q,r; double M=1,D;
  bool first = true;
  while ( ! f2.is_zero() ) { 
    Polynomial<double>::pseudo_div(f1,f2,q,r,D);
    if (!first) M*=D;
    TRACEV(f1);TRACEV(f2);TRACEV(q);TRACEV(r);TRACEV(M);
    r /= r.content();
    f1=f2; f2=r;
    first=false;
  }
  TRACEV(f1.content());
  return Polynomial<double>(1)*f1.abs();
}

int    Polynomial<int>::R_;
double Polynomial<double>::R_;

}//end namespace CGAL
