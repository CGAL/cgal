#include <CGAL/assertions.h>
#include <CGAL/algorithm.h>
#include <boost/heap/fibonacci_heap.hpp>
#include <iostream>
#include <functional> 

//first is the index (id) and second the value compared, 
typedef std::pair<std::size_t,int> Type;

struct More{
  bool operator()(const Type* t1,const Type* t2) const {
    return t1->second > t2->second;
  }
};

int main()
{
  #ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_USE_RELAXED_HEAP
  //testing min-heap
  typedef boost::heap::fibonacci_heap<Type*,boost::heap::compare<More> > Queue;
  typedef Queue::handle_type Handle;
  typedef std::vector<Queue>::iterator Iterator;
  Queue q;
  assert( q.size() == 0 );
  assert( q.empty() );
  
  std::vector<Type> data;
  std::map<Type,Handle> h;
  
  data.push_back(Type(0,10));
  data.push_back(Type(1,20));
  data.push_back(Type(2,30));
  data.push_back(Type(3,40));
  data.push_back(Type(4,1));
  data.push_back(Type(5,2));
  
  h[data[0]] = (q.push(&data[0]));
  h[data[1]] = (q.push(&data[0]+1));
  h[data[2]] = (q.push(&data[0]+2));
  h[data[3]] = (q.push(&data[0]+3));
  
  assert( q.top()->first == 0 );
  assert( q.size() == 4 );
  
  q.pop();
  assert( q.top()->first == 1 );
  assert( q.size() == 3 );
  
  h[data[4]]=q.push(&data[0]+4);
  assert( q.top()->first == 4 );
  assert( q.size() == 4 );
  
  q.erase(h[data[4]]);
  assert( q.top()->first == 1 );
  assert( q.size() == 3 ); 
  
  h[data[5]]=q.push(&data[0]+5);
  assert( q.top()->first == 5 );
  assert( q.size() == 4 );
  Handle ex_h5 = h[data[5]];
  data[5].second=43;
  q.update(ex_h5);
  h[data[5]] = ex_h5;
  assert( q.top()->first == 1 );
  assert( q.size() == 4 );  
  
  q.pop();
  assert( q.top()->first == 2 );
  assert( q.size() == 3 );  

  q.pop();
  assert( q.top()->first == 3 );
  assert( q.size() == 2 );  

  q.pop();
  assert( q.top()->first == 5 );
  assert( q.size() == 1 );  
  
  q.pop();
  assert( q.size() == 0 );
  assert( q.empty() );
  
  h[data[0]] = q.push(&data[0]);
  h[data[1]] = q.push(&data[0]+1);
  h[data[2]] = q.push(&data[0]+2);
  h[data[3]] = q.push(&data[0]+3);
  
  assert( q.top()->first == 0 );
  assert( q.size() == 4 );
  
  q.erase(h[data[1]]);
  assert( q.top()->first == 0 );
  assert( q.size() == 3 );

  q.erase(h[data[2]]);
  assert( q.top()->first == 0 );
  assert( q.size() == 2 );

  q.erase(h[data[0]]);
  assert( q.top()->first == 3 );
  assert( q.size() == 1 );
  
  q.erase(h[data[3]]);
  assert( q.size() == 0 );
  assert( q.empty() );  

//testing correctness of the order
  int array[10] = {0,1,2,3,4,5,6,7,8,9};
  CGAL::cpp98::random_shuffle(array,array+10);
  data.clear();
  data.reserve(10);
  for (int i=0;i<10;++i){
    data.push_back(Type(i,array[i]));
    q.push(&data[0]+i);
  }
  
  for (int i=0;i<10;++i){
    assert(q.top()->second==i);
    q.pop();
  }
  assert( q.empty() );
  
//testing update (increase key)
  h.clear();
  data.clear();
  data.reserve(10);
  for (int i=0;i<10;++i){
    data.push_back(Type(i,10+i));
    h[data[i]] = q.push(&data[0]+i);
  }

  for (unsigned int i=0;i<10;++i){
    Handle ex_hi = h[data[i]];
    data[i].second=9-i;
    q.update(ex_hi);
    h[data[i]] = ex_hi;
    assert(q.top()->first==i);
  }

//testing contains
  for (int i=0;i<10;++i){
    q.erase(h[data[i]]);
    assert(q.size()==9-i);
  }
  
//testing update (decrease key of top)
  h.clear();
  data.clear();
  data.reserve(10);
  for (int i=0;i<10;++i){
    data.push_back(Type(i,i));
    h[data[i]] = q.push(&data[0]+i);
  }

  for (unsigned int i=0;i<9;++i){
    Handle ex_hi = h[data[i]];
    data[i].second=10+i;
    q.update(ex_hi);
    h[data[i]] = ex_hi;
    assert(q.top()->first==i+1);
  }
  
//revert order
  for (unsigned int i=0;i<10;++i){
    Handle ex_hi = h[data[9-i]];
    data[9-i].second=i;
    q.update(ex_hi);
    h[data[9-i]] = ex_hi;
    assert(q.top()->first==9);
  }  
//testing remove (emulate pop)  
  for (std::size_t i=0;i<10;++i){
    assert(q.top()->first==9-i);
    q.erase(h[data[9-i]]);
  }
  assert( q.empty() );
  
  //testing remove+contains
  h.clear();
  data.clear();
  data.reserve(10);
  for (int i=0;i<10;++i){
    data.push_back(Type(i,i));
    h[data[i]] = q.push(&data[0]+i);
  }
  
  for (std::size_t i=0;i<10;++i){
    assert(q.top()->first==0);
    q.erase(h[data[-i+9]]);
  }
  assert( q.empty() );
  
  std::cout << "OK" << std::endl;
 
  return 0;
  #else
  std::cerr << "ERROR: Nothing is tested" << std::endl;
  return 1;
  #endif
}
