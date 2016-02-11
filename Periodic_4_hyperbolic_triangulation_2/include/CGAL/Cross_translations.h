#ifndef CGAL_HYPERBOLIC_CROSS_TRANSLATIONS_2_H
#define CGAL_HYPERBOLIC_CROSS_TRANSLATIONS_2_H

#include <CGAL/Hyperbolic_isometry_2.h>


template<typename Gt>
class Cross_translations
{
public:

  enum Direction {
    A = 0,  // 0
    B,      // 1
    InvA,   // 2
    InvB,   // 3
    C,      // 4
    D,      // 5
    InvC,   // 6
    InvD    // 7
  };

  typedef typename Gt::FT FT;
  typedef typename std::complex<FT> complex;
  typedef CGAL::Hyperbolic_isometry_2<Gt> Hyperbolic_isometry;
  typedef std::vector<Hyperbolic_isometry> Vector;
  typedef typename Vector::iterator Vector_iterator;
  

  Cross_translations() {
    computed = false;
    compute();
  }
  
  Hyperbolic_isometry& a()
  {
    compute();
    return g[A];
  }
    
  Hyperbolic_isometry& b()
  {
    compute();
    return g[B];
  } 
    
  Hyperbolic_isometry& c()
  {
    compute();
    return g[C];
  }
    
  Hyperbolic_isometry& d()
  {
    compute();
    return g[D];
  }
  
  const Vector& get_vector_of_translations()
  {
    compute();
    return g; 
  }
  
  Vector_iterator vector_begin()
  {
    compute();
    return g.begin();
  }
  
  Vector_iterator vector_end()
  {
    compute();
    return g.end();
  }


private:
  
  void compute_g()
  {

    const FT k1 = (FT(2) + CGAL::sqrt(2.))/FT(2);
    const FT k2 = CGAL::sqrt(CGAL::sqrt(2.));
    const FT k3 = (CGAL::sqrt(2.)*k2)/FT(2);
    
    std::complex<FT> m(k1, k1);
    std::complex<FT> n(k2*k1, -k3);
    
    g.resize(8);
    
    // a
    g[A] = Hyperbolic_isometry(conj(m), conj(n));
    g[InvA] = g[A].inverse();    
    
    // b
    g[B] = Hyperbolic_isometry(m, -n);
    g[InvB] = g[B].inverse();

    // c
    g[C] = Hyperbolic_isometry(conj(m), -conj(n));
    g[InvC] = g[C].inverse();

    // d
    g[D] = Hyperbolic_isometry(m, n);
    g[InvD] = g[D].inverse();

  }
  

  void compute()
  {
    if(!computed) {
      compute_g();
      computed = true;
    }
  }
  
  bool    computed;
  Vector  g;

};

 
#endif // CGAL_HYPERBOLIC_CROSS_TRANSLATIONS_2_H
