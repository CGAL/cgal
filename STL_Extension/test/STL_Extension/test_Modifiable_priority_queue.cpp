#include <CGAL/assertions.h>
#include <CGAL/Modifiable_priority_queue.h>
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

template <class Queue>
int queue_size(Queue& q,int n){
  int k=0;
  Type* t=new Type();
  for (int i=0;i<n;++i){
    *t = Type(i,0);
    if ( q.contains(t) ) ++k;
  }
  delete t;
  return k;
}

int main()
{
  #ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_USE_RELAXED_HEAP
  //testing min-heap
  typedef CGAL::Modifiable_priority_queue<Type*,Less,First_of_pair> Queue;
  Queue q(45,Queue::Compare(),Queue::ID());
  assert( queue_size(q,45) == 0 );
  assert( q.empty() );
  
  std::vector<Type> data;
  
  data.push_back(Type(0,10));
  data.push_back(Type(1,20));
  data.push_back(Type(2,30));
  data.push_back(Type(3,40));
  data.push_back(Type(4,1));
  data.push_back(Type(5,2));
  
  q.push(&data[0]);
  q.push(&data[0]+1);
  q.push(&data[0]+2);
  q.push(&data[0]+3);
  
  assert( q.top()->first == 0 );
  assert( queue_size(q,45) == 4 );
  
  q.pop();
  assert( q.top()->first == 1 );
  assert( queue_size(q,45) == 3 );
  
  q.push(&data[0]+4);
  assert( q.top()->first == 4 );
  assert( queue_size(q,45) == 4 );
  
  q.erase(&data[0]+4,false);
  assert( q.top()->first == 1 );
  assert( queue_size(q,45) == 3 ); 
  
  q.push(&data[0]+5);
  assert( q.top()->first == 5 );
  assert( queue_size(q,45) == 4 );
  
  data[5].second=43;
  q.update(&data[0]+5,true);
  assert( q.top()->first == 1 );
  assert( queue_size(q,45) == 4 );  
  
  q.pop();
  assert( q.top()->first == 2 );
  assert( queue_size(q,45) == 3 );  

  q.pop();
  assert( q.top()->first == 3 );
  assert( queue_size(q,45) == 2 );  

  q.pop();
  assert( q.top()->first == 5 );
  assert( queue_size(q,45) == 1 );  
  
  q.pop();
  assert( queue_size(q,45) == 0 );
  assert( q.empty() );
  
  q.push(&data[0]);
  q.push(&data[0]+1);
  q.push(&data[0]+2);
  q.push(&data[0]+3);
  
  assert( q.top()->first == 0 );
  assert( queue_size(q,45) == 4 );
  
  q.erase(&data[0]+1,true);
  assert( q.top()->first == 0 );
  assert( queue_size(q,45) == 3 );

  q.erase(&data[0]+2,true);
  assert( q.top()->first == 0 );
  assert( queue_size(q,45) == 2 );

  q.erase(&data[0],true);
  assert( q.top()->first == 3 );
  assert( queue_size(q,45) == 1 );
  
  q.erase(&data[0]+3,true);
  assert( queue_size(q,45) == 0 );
  assert( q.empty() );  

//testing correctness of the order
  int array[10] = {0,1,2,3,4,5,6,7,8,9};
  std::random_shuffle(array,array+10);
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
  data.clear();
  data.reserve(10);
  for (int i=0;i<10;++i){
    data.push_back(Type(i,10+i));
    q.push(&data[0]+i);
  }

  for (unsigned int i=0;i<10;++i){
    data[i].second=9-i;
    q.update(&data[0]+i,true);
    assert(q.top()->first==i);
  }

//testing contains
  for (int i=0;i<10;++i){
    q.erase(&data[0]+i,true);
    assert(queue_size(q,45)==9-i);
  }
  
//testing update (decrease key of top)
  data.clear();
  data.reserve(10);
  for (int i=0;i<10;++i){
    data.push_back(Type(i,i));
    q.push(&data[0]+i);
  }

  for (unsigned int i=0;i<9;++i){
    data[i].second=10+i;
    q.update(&data[0]+i,true);
    assert(q.top()->first==i+1);
  }
  
//revert order
  for (unsigned int i=0;i<10;++i){
    data[9-i].second=i;
    q.update(&data[0]+9-i,true);
    assert(q.top()->first==9);
  }  
//testing remove (emulate pop)  
  for (std::size_t i=0;i<10;++i){
    assert(q.top()->first==9-i);
    q.erase(&data[0]-i+9,true);
  }
  assert( q.empty() );
  
  //testing remove+contains
  data.clear();
  data.reserve(10);
  for (int i=0;i<10;++i){
    data.push_back(Type(i,i));
    q.push(&data[0]+i);
  }
  
  for (std::size_t i=0;i<10;++i){
    assert(q.top()->first==0);
    q.erase(&data[0]-i+9,true);
    for (std::size_t k=0;k<9-i;++k)
      assert(q.contains(&data[0]+k)==true);
    for (std::size_t k=0;k<i+1;++k)
      assert(q.contains(&data[0]+9-k)==false);
  }
  assert( q.empty() );
  
  std::cout << "OK" << std::endl;
 
  return 0;
  #else
  std::cerr << "ERROR: Nothing is tested" << std::endl;
  return 1;
  #endif
}
