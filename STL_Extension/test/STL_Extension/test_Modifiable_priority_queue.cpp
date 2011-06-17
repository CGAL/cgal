#include <CGAL/Cartesian.h>
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
  CGAL_assertion( queue_size(q,45) == 0 );
  CGAL_assertion( q.empty() );
  
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
  
  CGAL_assertion( q.top()->first == 0 );
  CGAL_assertion( queue_size(q,45) == 4 );
  
  q.pop();
  CGAL_assertion( q.top()->first == 1 );
  CGAL_assertion( queue_size(q,45) == 3 );
  
  q.push(&data[0]+4);
  CGAL_assertion( q.top()->first == 4 );
  CGAL_assertion( queue_size(q,45) == 4 );
  
  q.erase(&data[0]+4,false);
  CGAL_assertion( q.top()->first == 1 );
  CGAL_assertion( queue_size(q,45) == 3 ); 
  
  q.push(&data[0]+5);
  CGAL_assertion( q.top()->first == 5 );
  CGAL_assertion( queue_size(q,45) == 4 );
  
  data[5].second=43;
  q.update(&data[0]+5,true);
  CGAL_assertion( q.top()->first == 1 );
  CGAL_assertion( queue_size(q,45) == 4 );  
  
  q.pop();
  CGAL_assertion( q.top()->first == 2 );
  CGAL_assertion( queue_size(q,45) == 3 );  

  q.pop();
  CGAL_assertion( q.top()->first == 3 );
  CGAL_assertion( queue_size(q,45) == 2 );  

  q.pop();
  CGAL_assertion( q.top()->first == 5 );
  CGAL_assertion( queue_size(q,45) == 1 );  
  
  q.pop();
  CGAL_assertion( queue_size(q,45) == 0 );
  CGAL_assertion( q.empty() );
  
  q.push(&data[0]);
  q.push(&data[0]+1);
  q.push(&data[0]+2);
  q.push(&data[0]+3);
  
  CGAL_assertion( q.top()->first == 0 );
  CGAL_assertion( queue_size(q,45) == 4 );
  
  q.erase(&data[0]+1,true);
  CGAL_assertion( q.top()->first == 0 );
  CGAL_assertion( queue_size(q,45) == 3 );

  q.erase(&data[0]+2,true);
  CGAL_assertion( q.top()->first == 0 );
  CGAL_assertion( queue_size(q,45) == 2 );

  q.erase(&data[0],true);
  CGAL_assertion( q.top()->first == 3 );
  CGAL_assertion( queue_size(q,45) == 1 );
  
  q.erase(&data[0]+3,true);
  CGAL_assertion( queue_size(q,45) == 0 );
  CGAL_assertion( q.empty() );  

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
    CGAL_assertion(q.top()->second==i);
    q.pop();
  }
  CGAL_assertion( q.empty() );
  
//testing update
  data.clear();
  data.reserve(10);
  for (int i=0;i<10;++i){
    data.push_back(Type(i,10+i));
    q.push(&data[0]+i);
  }

  for (std::size_t i=0;i<10;++i){
    data[i].second=9-i;
    q.update(&data[0]+i,true);
    CGAL_assertion(q.top()->first==i);
  }
  
//testing remove (emulate pop)
  for (std::size_t i=0;i<10;++i){
    CGAL_assertion(q.top()->first==9-i);
    q.erase(&data[0]-i+9,true);
  }
  CGAL_assertion( q.empty() );
  
  //testing remove+contains
  data.clear();
  data.reserve(10);
  for (int i=0;i<10;++i){
    data.push_back(Type(i,i));
    q.push(&data[0]+i);
  }
  
  for (std::size_t i=0;i<10;++i){
    CGAL_assertion(q.top()->first==0);
    q.erase(&data[0]-i+9,true);
    for (std::size_t k=0;k<9-i;++k)
      CGAL_assertion(q.contains(&data[0]+k)==true);
    for (std::size_t k=0;k<i+1;++k)
      CGAL_assertion(q.contains(&data[0]+9-k)==false);
  }
  CGAL_assertion( q.empty() );
  
  std::cout << "OK" << std::endl;
 
  return 0;
  #else
  std::cerr << "ERROR: Nothing is tested" << std::endl;
  return 1;
  #endif
}