#ifndef OBJECT_INDEX_H
#define OBJECT_INDEX_H

#include <CGAL/basic.h>
#include <CGAL/Hash_map.h>
#include <string>
#include <strstream>

CGAL_BEGIN_NAMESPACE

template <typename I>
class Object_index {
  char _prefix;
  CGAL::Hash_map<I,int> _index;
public:
  Object_index() : _prefix('\0'), _index(-1) {}
  Object_index(I first, I beyond, char c=' ') : _prefix(c), _index(-1)
  { for(int i=0 ; first!=beyond; ++i,++first) _index[first]=i; }
  int operator[](const I& it) const { return _index[it]; } 
  int& operator[](const I& it) { return _index[it]; } 

  void index(I first, I beyond, char c=' ')
  { _prefix=c;
    for(int i=0 ; first!=beyond; ++i,++first) _index[first]=i;
  }
  std::string operator()(const I& it, bool verbose=true) const
  { if (verbose && _index[it]==-1) return "nil";
    if (verbose && _index[it]==-2) return "end";
    std::ostrstream os; 
    if (verbose) os << _prefix;
    os << _index[it] << '\0';    
    std::string res(os.str()); os.freeze(0); return res; }
 
};

CGAL_END_NAMESPACE

#endif //OBJECT_INDEX_H




