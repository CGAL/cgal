#ifndef CGAL_HYPERBOLIC_DIAMETRIC_TRANSLATIONS_2_H
#define CGAL_HYPERBOLIC_DIAMETRIC_TRANSLATIONS_2_H

#include <CGAL/Hyperbolic_isometry_2.h>  

template<typename Gt>
class Diametric_translations
{
public:

  enum Direction {
    A = 0,  // 0
    InvB,   // 1
    C,      // 2
    InvD,   // 3
    InvA,   // 4
    B,      // 5
    InvC,   // 6
    D       // 7
  };

  typedef typename Gt::FT FT;
  typedef typename std::complex<FT> complex;
  typedef CGAL::Hyperbolic_isometry_2<Gt> Hyperbolic_isometry;
  typedef std::vector<Hyperbolic_isometry> Vector;
  typedef typename Vector::iterator Vector_iterator;
  

  Diametric_translations() {
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

    const FT k1 = FT(1) + CGAL::sqrt(2.);
    const FT k2 = FT(0);
        
    std::complex<FT> m(k1, k2);
        
    // This is the multiplicative factor for all b's
    const FT k3 = CGAL::sqrt(2.) * CGAL::sqrt( k1 );
        
    std::vector< std::complex<FT> > n;
        
    // Euler's identity: exp(i k \theta) = cos(k \theta) + i sin(k \theta)
    for (int kk = 0; kk < 8; kk++) {
        n.push_back( std::complex<FT>( k3 * cos( FT(kk) * CGAL_PI / FT(4)), k3 * sin( FT(kk) * CGAL_PI / FT(4) ) ) );
    }
    
    g.resize(8);
    
    // a
    g[A] = Hyperbolic_isometry(m, n[A]);
    g[InvA] = g[A].inverse();    
    
    // b
    g[B] = Hyperbolic_isometry(m, n[B]);
    g[InvB] = g[B].inverse();

    // c
    g[C] = Hyperbolic_isometry(m, n[C]);
    g[InvC] = g[C].inverse();

    // d
    g[D] = Hyperbolic_isometry(m, n[D]);
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

 
#endif // CGAL_HYPERBOLIC_DIAMETRIC_TRANSLATIONS_2_H
