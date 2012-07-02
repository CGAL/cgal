#ifndef CGAL_GROUP_OF_INDEX_2
#define CGAL_GROUP_OF_INDEX_2

#include <Translations.h>

template<class Translation>
class Group_of_index_2
{
public:
  typedef Element<Translation> Element;
  typedef typename Translation::FT FT;
  
  typedef std::list<Element> List;
  typedef CGAL::Circulator_from_container<List> Circulator;
  typedef typename List::iterator List_iterator;
  typedef typename List::size_type Size_type;
  
  typedef std::pair<Translation, int> Node;
  typedef std::vector<Node> Vector;
  
  Group_of_index_2() : is_computed(false), l(&Translations<typename Translation::Geom_traits>::list())  
  {
  }
  
  void add_group(List* l_)
  {
    l = l_;
    
    generate();
  }
  
  List& list()
  {
    return l2;
  }
  
  const List& list() const
  {
    return l2;
  }
  
  /*
  template<class OutputIterator>
  OutputIterator group_of_index_2(OutputIterator oit) const
  {
    return std::copy(g2.begin(), g2.end(), oit);
  }*/
  
  void generate()
  {
    if(!is_computed) {
      compute();
      refine();
      
      is_computed = true;
    }
  }
  
  List_iterator begin()
  {
    return l2.begin();
  }
  
  List_iterator end() 
  {
    return l2.end();
  }
  
  Size_type number_of_elements() const
  {
    return l2.size();
  }
  
private:
  // the lists can not be copied.
  Group_of_index_2(const Group_of_index_2&);
  Group_of_index_2& operator= (const Group_of_index_2&) const;
  
  bool is_computed;
  void compute();
  void compute_via_vector();
  
  void refine(FT e = 0.000001);
  
  List* l;
  List l2;
  Vector g;
  Vector g2;
};


template<class Translation>
void Group_of_index_2<Translation>::refine(FT e)
{
  typedef typename Translation::Geom_traits Gt;
  typedef typename Translation::complex complex;
  
  Circulator begin(&l2);
  Circulator next = boost::next(begin);
  
  if(begin == 0 || begin == next) {
    return;
  }
  
  int steps_nb = l2.size();
  
  for(int i = 0; i < steps_nb; i++, begin = next, next = boost::next(begin)) {
    complex dif_m = begin->g.m() - next->g.m();
    complex dif_n = begin->g.n() - next->g.n();
    
    if(std::norm(dif_m) < e && std::norm(dif_n) < e) {
      l2.erase(begin->inverse.current_iterator());
      l2.erase(begin.current_iterator());
      
      steps_nb--;
    }
  }
}

template<class Translation>
void Group_of_index_2<Translation>::compute_via_vector()
{
  int nb = g.size();
  g2.resize(2*nb);
  
  for(int i = 0; i < nb; i++) {
    int inv = g[i].second;
    int next_inv = (inv + 1) & nb-1;
    int prev_inv = (inv - 1) & nb-1;
    
    g2[2*i].first = g[i].first*g[next_inv].first;
    // compute position of the inverse element of g2[2*i]
    int inv_pos = g[next_inv].second;
    g2[2*i].second = 2*inv_pos + 1;
    
    g2[2*i + 1].first = g[i].first*g[prev_inv].first;
    // compute position of the inverse element of g2[2*i + 1]
    inv_pos = g[prev_inv].second;
    g2[2*i + 1].second = 2*inv_pos;
  }
}


template<class Translation>
void Group_of_index_2<Translation>::compute()
{
  typedef std::map<Element*, std::pair<Circulator, Circulator> > MapLtoL2;
  MapLtoL2 l_to_l2;
  typedef typename MapLtoL2::iterator Map_it;
  
  l2.resize(2*l->size());
  
  Circulator li(l);
  Circulator Next, Prev;
  Circulator l2i(&l2);

  if(li == 0) {
    return;
  }
  
  // add g
  
  do {
    Circulator item1 = l2i, item2 = boost::next(l2i);
    
    Next = li->inverse;
    Next = boost::next(Next);
    item1->g = li->g * Next->g;  
    
    // aux value
    item1->inverse = Next->inverse;
    
    Prev = li->inverse;
    Prev = boost::prior(Prev);
    item2->g = li->g * Prev->g;
    
    // aux value
    item2->inverse = Prev->inverse;
    
    l_to_l2.insert(std::make_pair(&*li, std::make_pair(item1, item2)));
    
    li = boost::next(li);
    l2i = boost::next(l2i, 2);
  } while( li != Circulator( l ) );
  
  // add inverse
  
  int i = 0;
  l2i = Circulator(&l2);
  do {  
    Map_it mit = l_to_l2.find(&*l2i->inverse);
    
    assert(mit != l_to_l2.end());
    std::pair<Circulator, Circulator> val = mit->second;
    
    if(i % 2 == 0) {
      l2i->inverse = val.second;
    } else {
      l2i->inverse = val.first;
    }
    
    l2i++;
    i++;
  } while( l2i != Circulator( &l2 ) );
}

#endif // CGAL_GROUP_OF_INDEX_2

