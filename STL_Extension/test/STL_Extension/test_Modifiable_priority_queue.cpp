#include <CGAL/assertions.h>
#include <CGAL/algorithm.h>
#include <CGAL/Modifiable_priority_queue.h>
#include <boost/unordered_map.hpp>
#include <iostream>
#include <functional> 

//first is the index (id) and second the value compared, 
typedef std::pair<std::size_t,int> Type;

//property map
struct First_of_pair{
  //classical typedefs
  typedef Type* key_type;
  typedef std::size_t value_type;
  typedef std::size_t reference;
  typedef boost::readable_property_map_tag category;
};
//get function for property map
First_of_pair::value_type
get(const First_of_pair&, const First_of_pair::key_type& k) {
  return k->first;
}

struct Less{
  bool operator()(const Type* t1,const Type* t2) const {
    return t1->second < t2->second;
  };
};

template <class Map>
int queue_size(Map& m, int n){
  return m.size();
}

int main()
{
  #ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_USE_RELAXED_HEAP
  //testing min-heap
  typedef CGAL::Modifiable_priority_queue<Type*,Less,First_of_pair> Queue;
  typedef boost::unordered_map<Type*, Queue::handle> Heap_map;
  Queue q(45,Queue::Compare(),Queue::ID());
  Heap_map h;
  assert( queue_size(h,45) == 0 );
  assert( q.empty() );
  
  std::vector<Type> data;
  
  data.push_back(Type(0,10));
  data.push_back(Type(1,20));
  data.push_back(Type(2,30));
  data.push_back(Type(3,40));
  data.push_back(Type(4,1));
  data.push_back(Type(5,2));
  
  h[&data[0]] = q.push(&data[0]);
  h[&data[0]+1] = q.push(&data[0]+1);
  h[&data[0]+2] = q.push(&data[0]+2);
  h[&data[0]+3] = q.push(&data[0]+3);
  
  assert( q.top()->first == 0 );
  assert( queue_size(h,45) == 4 );
  
  Type *top = q.top();
  q.pop();
  h.erase(top);
  assert( q.top()->first == 1 );
  assert( queue_size(h,45) == 3 );
  
  h[&data[0]+4] = q.push(&data[0]+4);
  
  assert( q.top()->first == 4 );
  assert( queue_size(h,45) == 4 );

  q.erase(&data[0]+4, h[&data[0]+4]);
  h.erase(&data[0]+4);
  
  assert( q.top()->first == 1 );
  assert( queue_size(h,45) == 3 ); 
  
  h[&data[0]+5] = q.push(&data[0]+5);
  
  assert( q.top()->first == 5 );
  assert( queue_size(h,45) == 4 );
  
  Queue::handle ex_h = h[&data[0]+5];
  data[5].second=43;
  q.update(&data[0]+5,ex_h);
  h[&data[0]+5] = ex_h;
  
  assert( q.top()->first == 1 );
  assert( queue_size(h,45) == 4 );  
  
  top = q.top();
  q.pop();
  h.erase(top);
  
  assert( q.top()->first == 2 );
  assert( queue_size(h,45) == 3 );  

  top=q.top();
  q.pop();
  h.erase(top);
  assert( q.top()->first == 3 );
  assert( queue_size(h,45) == 2 );  

  top = q.top();
  q.pop();
  h.erase(top);
  assert( q.top()->first == 5 );
  assert( queue_size(h,45) == 1 );  
  
  top=q.top();
  q.pop();
  h.erase(top);
  assert( queue_size(h,45) == 0 );
  assert( q.empty() );
  
  h[&data[0]] = q.push(&data[0]);
  h[&data[0]+1] = q.push(&data[0]+1);
  h[&data[0]+2] = q.push(&data[0]+2);
  h[&data[0]+3] = q.push(&data[0]+3);
  
  assert( q.top()->first == 0 );
  assert( queue_size(h,45) == 4 );
  
  q.erase(&data[0]+1,h[&data[0]+1]);
  h.erase(&data[0]+1);
  assert( q.top()->first == 0 );
  assert( queue_size(h,45) == 3 );

  q.erase(&data[0]+2,h[&data[0]+2]);
  h.erase(&data[0]+2);
  assert( q.top()->first == 0 );
  assert( queue_size(h,45) == 2 );

  q.erase(&data[0],h[&data[0]]);
  h.erase(&data[0]);
  assert( q.top()->first == 3 );
  assert( queue_size(h,45) == 1 );
  
  q.erase(&data[0]+3,h[&data[0]+3]);
  h.erase(&data[0]+3);
  assert( queue_size(h,45) == 0 );
  assert( q.empty() );  

//testing correctness of the order
  int array[10] = {0,1,2,3,4,5,6,7,8,9};
  CGAL::cpp98::random_shuffle(array,array+10);
  data.clear();
  data.reserve(10);
  for (int i=0;i<10;++i){
    data.push_back(Type(i,array[i]));
    h[&data[0]+i] = q.push(&data[0]+i);
  }
  
  for (int i=0;i<10;++i){
    top = q.top();
    assert(top->second==i);
    q.pop();
    h.erase(top);
  }
  assert( q.empty() );
  
//testing update (increase key)
  data.clear();
  data.reserve(10);
  for (int i=0;i<10;++i){
    data.push_back(Type(i,10+i));
    h[&data[0]+i] = q.push(&data[0]+i);
  }
  
  for (unsigned int i=0;i<10;++i){
    ex_h = h[&data[0]+i];
    data[i].second=9-i;
    q.update(&data[0]+i,ex_h);
    h[&data[0]+i] = ex_h;
    assert(q.top()->first==i);
  }

//testing contains
  for (int i=0;i<10;++i){
    q.erase(&data[0]+i,h[&data[0]+i]);
    h.erase(&data[0]+i);
    assert(queue_size(h,45)==9-i);
  }
  
//testing update (decrease key of top)
  data.clear();
  data.reserve(10);
  for (int i=0;i<10;++i){
    data.push_back(Type(i,i));
    h[&data[0]+i] = q.push(&data[0]+i);
  }

  for (unsigned int i=0;i<9;++i){
    ex_h = h[&data[0]+i];
    data[i].second=10+i;
    q.update(&data[0]+i,ex_h);
    h[&data[0]+i] = ex_h;
    assert(q.top()->first==i+1);
  }
  
//revert order
  for (unsigned int i=0;i<10;++i){
    ex_h = h[&data[0]+9-i];
    data[9-i].second=i;
    q.update(&data[0]+9-i,ex_h);
    h[&data[0]+9-i] = ex_h;
    assert(q.top()->first==9);
  }  
//testing remove (emulate pop)  
  for (std::size_t i=0;i<10;++i){
    assert(q.top()->first==9-i);
    q.erase(&data[0]-i+9,h[&data[0]+9-i]);
    h.erase(&data[0]+9-i);
  }
  assert( q.empty() );
  
  //testing remove+contains
  data.clear();
  data.reserve(10);
  for (int i=0;i<10;++i){
    data.push_back(Type(i,i));
    h[&data[0]+i] = q.push(&data[0]+i);
  }
  
  for (std::size_t i=0;i<10;++i){
    assert(q.top()->first==0);
    q.erase(&data[0]-i+9,h[&data[0]-i+9]);
    h.erase(&data[0]-i+9);
    for (std::size_t k=0;k<9-i;++k)
      assert(h.find(&data[0]+k) != h.end());
    for (std::size_t k=0;k<i+1;++k)
      assert(h.find(&data[0]+9-k) == h.end());
  }
  assert( q.empty() );
  
  std::cout << "OK" << std::endl;
 
  return 0;
  #else
  std::cerr << "ERROR: Nothing is tested" << std::endl;
  return 1;
  #endif
}
