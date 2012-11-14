
#include <vector>

namespace CGAL {
  template <typename T, int N=10>
class Small_stack {

private:
  std::vector<T> V;

  int TOP;
  
  T ts[N];

public:
  
  Small_stack()
    : TOP(-1)
  {}

  void push_back(const T& t)
  {
    if(TOP < (N-1)){
      ++TOP;
      ts[TOP] = t;
    } else {
      V.push_back(t);
    }
  }

  const T& back() const
  {
    if(TOP < N){
      return ts[TOP];
    } else {
      return V.back();
    }    
  }

  void pop_back()
  {
    if(TOP < N){
      --TOP;
    } else {
      V.pop_back();
    }    
  }

  bool empty() const
  {
    return TOP == -1;
  }
  
};

} // namespace
