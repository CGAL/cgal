#ifndef CGAL_HYPERBOLIC_BASE_TRANSLATIONS_2_H
#define CGAL_HYPERBOLIC_BASE_TRANSLATIONS_2_H

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
class Base_hyperbolic_translations
{
public:
  typedef typename Gt::FT FT;
  typedef typename std::complex<FT> complex;
  typedef CGAL::Hyperbolic_isometry_2<Gt> Hyperbolic_isometry;
  
  //typedef Element<Hyperbolic_isometry> Element_t;
  //typedef std::list<Element_t> List;
  //typedef typename List::iterator List_iterator;
  //typedef CGAL::Circulator_from_container<List> Circulator;
  
  typedef std::pair<Hyperbolic_isometry, int> Node;
  typedef std::vector<Node> Vector;
  typedef typename Vector::iterator Vector_iterator;
  
  
  virtual Hyperbolic_isometry& a() = 0;
  virtual Hyperbolic_isometry& b() = 0;
  virtual Hyperbolic_isometry& c() = 0;
  virtual Hyperbolic_isometry& d() = 0;
  
  const Vector& get_vector_of_translations();
  Vector_iterator vector_begin();
  Vector_iterator vector_end();
  //List_iterator list_begin();
  //List_iterator list_end();
  //List& list();
    
private:
  
  void compute_g();
  //void compute_l();
  void compute();
  
  Vector g;
  //List l;
};
 
#endif // CGAL_HYPERBOLIC_BASE_TRANSLATIONS_2_H
