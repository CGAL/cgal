#ifndef CGAL_FILTRED_CONTAINER_H
#define CGAL_FILTRED_CONTAINER_H

CGAL_BEGIN_NAMESPACE

template <class Cont, class Pred>
class Filtred_container
{
  Cont cont;
  Pred test;
public:
  Filtred_container(Pred p=Pred()) : cont(), test(p) {};
  Filtred_container(Cont& c, Pred p=Pred()) : cont(c), test(p) {};

  typedef typename Cont::reference reference;
  typedef typename Cont::iterator iterator;
  typedef typename Cont::value_type value_type;

  inline
  reference front()
    {
      iterator r=cont.begin();
      while(!test(*r))
	{
	  cont.erase(r);
	  r=cont.begin();
	}
      return *r;
    }

  inline
  bool empty()
    {
      if(cont.empty())
	return true;
      else
	{
	  while(!cont.empty())
	    {
	      iterator r=cont.begin();
	    if(!test(*r))
	      pop_front();
	    else
	      return false;
	    }
	  return true;
	}
    }

  inline
  void pop_front() { cont.erase(cont.begin()); }

  inline
  void push_back(const value_type& e) { cont.push_back(e); }

  inline
  void insert(const value_type& e) { cont.insert(e); }

  inline
  void clear() { cont.clear(); }
};

CGAL_END_NAMESPACE

#endif
