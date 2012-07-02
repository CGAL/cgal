#ifndef CGAL_HYPERBOLIC_TRANSLATIONS_2_H
#define CGAL_HYPERBOLIC_TRANSLATIONS_2_H

#include <CGAL/Hyperbolic_isometry_2.h>

template<typename Translation>
struct Element
{
  typedef typename std::list<Element> List;
  typedef CGAL::Circulator_from_container<List> Circulator;
  
  Translation g;
  
  // circulator iterator to an inverse translation in the list
  Circulator inverse;
};

template<typename Gt>
class Translations
{
public:
  typedef typename Gt::FT FT;
  typedef typename std::complex<FT> complex;
  typedef CGAL::Hyperbolic_isometry_2<Gt> Hyperbolic_isometry;
  
  typedef Element<Hyperbolic_isometry> Element;
  typedef std::list<Element> List;
  typedef typename List::iterator List_iterator;
  typedef CGAL::Circulator_from_container<List> Circulator;
  
  typedef std::pair<Hyperbolic_isometry, int> Node;
  typedef std::vector<Node> Vector;
  typedef typename Vector::iterator Vector_iterator;
  
  
  static Hyperbolic_isometry& a()
  {
    compute();
    return g[0].first;
  }
  
  static Hyperbolic_isometry& b()
  {
    compute();
    return g[1].first;
  }
  
  static Hyperbolic_isometry& c()
  {
    compute();
    return g[2].first;
  }
  
  static Hyperbolic_isometry& d()
  {
    compute();
    return g[3].first;
  }
  
  static const Vector& get_vector_of_translations()
  {
    compute();
    return g; 
  }
  
  static Vector_iterator vector_begin()
  {
    compute();
    return g.begin();
  }
  
  static Vector_iterator vector_end()
  {
    compute();
    return g.end();
  }
  
  static List_iterator list_begin()
  {
    compute();
    return l.begin();
  }
  
  static List_iterator list_end()
  {
    compute();
    return l.end();
  }
  
  static List& list()
  {
    compute();
    return l;
  }
    
private:
  
  static void compute_g()
  {
    const FT k1 = (FT(2) + CGAL::sqrt(2.))/FT(2);
    const FT k2 = CGAL::sqrt(CGAL::sqrt(2.));
    const FT k3 = (CGAL::sqrt(2.)*k2)/FT(2);
    
    std::complex<FT> m(k1, k1);
    std::complex<FT> n(k2*k1, -k3);
    
    g.resize(8);
    
    // a
    g[0].first = Hyperbolic_isometry(conj(m), conj(n));
    g[0].second = 2;
    
    // b
    g[3].first = Hyperbolic_isometry(m, -n);
    g[3].second = 1;
    // c
    g[4].first = Hyperbolic_isometry(conj(m), -conj(n));
    g[4].second = 6;
    // d
    g[7].first = Hyperbolic_isometry(m, n);
    g[7].second = 5;
    
    int index = g[0].second;
    g[index].first = g[0].first.inverse();
    g[index].second = 0;
    
    index = g[3].second;
    g[index].first = g[3].first.inverse();
    g[index].second = 3;
    
    index = g[4].second;
    g[index].first = g[4].first.inverse();
    g[index].second = 4;
    
    index = g[7].second;
    g[index].first = g[7].first.inverse();
    g[index].second = 7;
  }
  
  
  static void compute_l()
  {
    l.resize(g.size());
    
    std::vector<Circulator> aux_list;
    aux_list.reserve(8);
    
    for(List_iterator li = l.begin(); li != l.end(); li++) {
      aux_list.push_back( Circulator(&l, li) );
    }
    
    for(typename List::size_type i = 0; i < aux_list.size(); i++) {
      aux_list[i]->g = g[i].first;
      aux_list[i]->inverse = aux_list[g[i].second];
    }
  }
  
  static void compute()
  {
    static bool computed = false;
    if(!computed) {
      compute_g();
      compute_l();
      computed = true;
    }
  }
  
  static Vector g;
  
  static List l;
};

// default initialization

template<typename Gt>
typename Translations<Gt>::Vector
Translations<Gt>::g;

template<typename Gt>
typename Translations<Gt>::List
Translations<Gt>::l;
 
#endif // CGAL_HYPERBOLIC_TRANSLATIONS_2_H
