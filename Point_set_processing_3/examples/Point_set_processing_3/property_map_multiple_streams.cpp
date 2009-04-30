
#include <CGAL/Simple_cartesian.h>
#include <CGAL/property_map.h>
#include <algorithm>
#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;

typedef std::pair<Point_3, Vector_3> PointVectorPair; 

typedef boost::tuple<int, Point_3, bool, Vector_3> IndexedPointWithOrientableNormalTuple;

namespace CGAL
{
  // stores the information needed to access the elements of one kind
  template<typename value_type,int _id=-1>
  struct pointer_property_map
  {
    enum { id = _id };
    char* data;       // pointer to first element
    const int stride; // number of bytes between two elements
    pointer_property_map(value_type* d, int s=sizeof(value_type))
      : data((char*)d), stride(s)
    {}
  };
  
  struct NullType {};

  // forward declaration
  template<typename pmap_cluster_type> class pointer_property_map_iterator;
  template<typename pmap_cluster_type> class pointer_pmap_object;
  template<typename pmap_cluster_type> class pointer_pmap_reference;
  template<typename pmap_cluster_type> class pointer_pmap_const_reference;
  
  // represents a cluster of up to 4 property streams
  template<typename T0, typename T1=NullType, typename T2=NullType, typename T3=NullType>
  class pointer_property_map_cluster
  {
    int m_size;
  public:
    typedef pointer_pmap_reference<pointer_property_map_cluster> pointer_pmap_reference_type;
    
    pointer_property_map_cluster(int size) : m_size(size) {}
    
    typedef pointer_property_map_iterator<pointer_property_map_cluster> iterator;
    iterator begin() { return iterator(this,0); }
    iterator end() { return iterator(this,m_size); }
    #define pointer_property_map_cluster_code(I) \
      public: \
        typedef T##I value_type##I; \
        pointer_property_map<T##I,I> pmap##I() { return pointer_property_map<T##I,I>(m_data##I,m_stride##I); } \
        void attach##I(T##I* data, int stride=sizeof(T##I)) { m_data##I=data; m_stride##I=stride; } \
        const T##I& get##I(int i) const { return *(const T##I*)((const char*)m_data##I+i*m_stride##I); } \
        T##I& get##I(int i) { return *(T##I*)((char*)m_data##I+i*m_stride##I); } \
      protected: \
        int m_stride##I; \
        T##I* m_data##I; \

    pointer_property_map_cluster_code(0)
    pointer_property_map_cluster_code(1)
    pointer_property_map_cluster_code(2)
    pointer_property_map_cluster_code(3)
    
    void swap(pointer_pmap_reference_type& a, pointer_pmap_reference_type& b)
    {
      if(!boost::is_same<T0,NullType>::value) std::swap(a.cluster->get0(a.index),b.cluster->get0(b.index));
      if(!boost::is_same<T1,NullType>::value) std::swap(a.cluster->get1(a.index),b.cluster->get1(b.index));
      if(!boost::is_same<T2,NullType>::value) std::swap(a.cluster->get2(a.index),b.cluster->get2(b.index));
      if(!boost::is_same<T3,NullType>::value) std::swap(a.cluster->get3(a.index),b.cluster->get3(b.index));
    }
  };
  
  // a reference to one object
  template<typename pmap_cluster_type>
  class pointer_pmap_reference
  {
  public:
    pointer_pmap_reference(pmap_cluster_type* c, int i)
      : cluster(c), index(i)
    {}
    explicit pointer_pmap_reference(const pointer_pmap_object<pmap_cluster_type>& other)
    {
      *this = other;
    }
    pointer_pmap_reference& operator=(const pointer_pmap_object<pmap_cluster_type>& other);
    pointer_pmap_reference& operator=(const pointer_pmap_reference<pmap_cluster_type>& other);
    pointer_pmap_reference& operator=(const pointer_pmap_const_reference<pmap_cluster_type>& other);
    
    pmap_cluster_type* cluster;
    int index;
  }; 
  
  // a reference to one object (read only)
  template<typename pmap_cluster_type>
  class pointer_pmap_const_reference
  {
  public:
    pointer_pmap_const_reference(pmap_cluster_type* c, int i)
      : cluster(c), index(i)
    {}
    const pmap_cluster_type* cluster;
    int index;
  };
  
  // a temporary object with storage
  template<typename pmap_cluster_type>
  class pointer_pmap_object
  {
  public:
    pointer_pmap_object(const pointer_pmap_reference<pmap_cluster_type>& other)
    { *this = other; }
    pointer_pmap_object(const pointer_pmap_const_reference<pmap_cluster_type>& other)
    { *this = other; }
    
    pointer_pmap_object& operator=(const pointer_pmap_reference<pmap_cluster_type>& other)
    {
      #define pointer_pmap_object_copy(I) \
        if(!boost::is_same<typename pmap_cluster_type::value_type0,NullType>::value) \
          m_value##I = other.cluster->get##I(other.index);
      pointer_pmap_object_copy(0);
      pointer_pmap_object_copy(1);
      pointer_pmap_object_copy(2);
      pointer_pmap_object_copy(3);
      #undef pointer_pmap_object_copy
      return *this;
    }
    
    pointer_pmap_object& operator=(const pointer_pmap_const_reference<pmap_cluster_type>& other)
    {
      #define pointer_pmap_object_copy(I) \
        if(!boost::is_same<typename pmap_cluster_type::value_type0,NullType>::value) \
          m_value##I = other.cluster->get##I(other.index);
      pointer_pmap_object_copy(0);
      pointer_pmap_object_copy(1);
      pointer_pmap_object_copy(2);
      pointer_pmap_object_copy(3);
      #undef pointer_pmap_object_copy
      return *this;
    }
    
    #define pointer_pmap_object_value(I) \
      protected: \
        typename pmap_cluster_type::value_type##I m_value##I; \
      public: \
        const typename pmap_cluster_type::value_type##I&  \
        value##I() const { return m_value##I; } \
        typename pmap_cluster_type::value_type##I&  \
        value##I() { return m_value##I; }
    pointer_pmap_object_value(0);
    pointer_pmap_object_value(1);
    pointer_pmap_object_value(2);
    pointer_pmap_object_value(3);
    #undef pointer_pmap_object_value
  };
  
  template<typename pmap_cluster_type>
  pointer_pmap_reference<pmap_cluster_type>&
  pointer_pmap_reference<pmap_cluster_type>::operator=(const pointer_pmap_object<pmap_cluster_type>& other)
  {
    #define pointer_pmap_reference_copy(I) \
      if(!boost::is_same<typename pmap_cluster_type::value_type0,NullType>::value) \
        cluster->get##I(index) = other.value##I();
    pointer_pmap_reference_copy(0);
    pointer_pmap_reference_copy(1);
    pointer_pmap_reference_copy(2);
    pointer_pmap_reference_copy(3);
    #undef pointer_pmap_reference_copy
    return *this;
  }
  
  template<typename pmap_cluster_type>
  pointer_pmap_reference<pmap_cluster_type>&
  pointer_pmap_reference<pmap_cluster_type>::operator=(const pointer_pmap_reference<pmap_cluster_type>& other)
  {
    #define pointer_pmap_reference_copy(I) \
      if(!boost::is_same<typename pmap_cluster_type::value_type0,NullType>::value) \
        cluster->get##I(index) = other.cluster->get##I(other.index);
    pointer_pmap_reference_copy(0);
    pointer_pmap_reference_copy(1);
    pointer_pmap_reference_copy(2);
    pointer_pmap_reference_copy(3);
    #undef pointer_pmap_reference_copy
    return *this;
  }
  
  template<typename pmap_cluster_type>
  pointer_pmap_reference<pmap_cluster_type>&
  pointer_pmap_reference<pmap_cluster_type>::operator=(const pointer_pmap_const_reference<pmap_cluster_type>& other)
  {
    #define pointer_pmap_reference_copy(I) \
      if(!boost::is_same<typename pmap_cluster_type::value_type0,NullType>::value) \
        cluster->get##I(index) = other.cluster->get##I(other.index);
    pointer_pmap_reference_copy(0);
    pointer_pmap_reference_copy(1);
    pointer_pmap_reference_copy(2);
    pointer_pmap_reference_copy(3);
    #undef pointer_pmap_reference_copy
    return *this;
  }
    
  // STL compatible iterator
  template<typename pmap_cluster_type>
  class pointer_property_map_iterator
  {
  public:
  
    typedef std::bidirectional_iterator_tag iterator_category;
    typedef pointer_pmap_object<pmap_cluster_type> value_type;
    typedef size_t difference_type ;
    typedef value_type* pointer;
    typedef pointer_pmap_reference<pmap_cluster_type> reference;
    
    pointer_property_map_iterator(pmap_cluster_type* c, int i=0) : cluster(c), index(i) {}
    
    pointer_property_map_iterator& operator++()
    { ++index; return *this; }
    pointer_property_map_iterator operator++(int)
    { return pointer_property_map_iterator(cluster,index++); }
    pointer_property_map_iterator operator+(difference_type i) const
    { return pointer_property_map_iterator(cluster,index+i); }
    pointer_property_map_iterator& operator+=(difference_type i)
    { index+=i; return *this; }
    
    pointer_property_map_iterator& operator--()
    { --index; return *this; }
    pointer_property_map_iterator operator--(int)
    { return pointer_property_map_iterator(cluster,index--); }
    pointer_property_map_iterator operator-(difference_type i) const
    { return pointer_property_map_iterator(cluster,index-i); }
    pointer_property_map_iterator& operator-=(difference_type i)
    { index+=i; return *this; }
    difference_type operator-(const pointer_property_map_iterator& other) const
    { return (index-other.index); }
    
    bool operator!=(const pointer_property_map_iterator& other) const
    { assert(cluster==other.cluster); return index!=other.index; }
    bool operator==(const pointer_property_map_iterator& other) const
    { assert(cluster==other.cluster); return index==other.index; }
    bool operator<(const pointer_property_map_iterator& other) const
    { assert(cluster==other.cluster); return index<other.index; }
    bool operator>(const pointer_property_map_iterator& other) const
    { assert(cluster==other.cluster); return index>other.index; }
    
    reference operator*() { return pointer_pmap_reference<pmap_cluster_type>(cluster,index); }
  protected:
    pmap_cluster_type* cluster;
    size_t index;
  };
}

namespace std
{
  // fast std::swap
  template<typename pmap_cluster_type>
  void swap(CGAL::pointer_pmap_reference<pmap_cluster_type>& a, CGAL::pointer_pmap_reference<pmap_cluster_type>& b)
  {
    pmap_cluster_type::swap(a,b);
  }
}

namespace boost
{
  
  template <typename T,int I,typename C>
  const T& get(CGAL::pointer_property_map<T,I> pm, const CGAL::pointer_pmap_const_reference<C>& w)
  {
    return *(const T*)(pm.data+w.index*pm.stride);
  }
  
  template <typename T,int I,typename C>
  T& get(CGAL::pointer_property_map<T,I> pm, const CGAL::pointer_pmap_reference<C>& w)
  {
    return *(T*)(pm.data+w.index*pm.stride);
  }
  
  #define pointer_pmap_object_get(I) \
    template <typename T,typename C> \
    const T& get(CGAL::pointer_property_map<T,I> pm, const CGAL::pointer_pmap_object<C>& w) \
    { return w.value##I(); }
  pointer_pmap_object_get(0);
  pointer_pmap_object_get(1);
  pointer_pmap_object_get(2);
  pointer_pmap_object_get(3);
  #undef pointer_pmap_object_get

  template <typename T,int I,typename C>
  void put(CGAL::pointer_property_map<T,I> pm, const CGAL::pointer_pmap_reference<C>& w, const T& v)
  {
    *(T*)(pm.data+w.index*pm.stride) = v;
  }
}


template <typename PointPmap>
struct MyLess {

  PointPmap pm;

  MyLess(const PointPmap& p)
    : pm(p)
  {}
  
  template<typename T0, typename T1>
  bool operator()(T0 t0, T1 t1) const
  {
    return boost::get(pm,t0) < boost::get(pm,t1);
  }

};


// In this example we have a function that only operates on the point part
// It sorts them lexicographically
template <typename Iterator, typename PointPmap >
void process_point_set(Iterator beg, Iterator end, PointPmap pm)
{
  MyLess<PointPmap> less(pm);
  std::sort(beg,end,less);
}


// We can call it just with points. Then interally we use a property map
// that maps points on points
// The identity_property_map is buggy so we use the one of CGAL
template <typename Iterator>
void process_point_set(Iterator beg, Iterator end)
{
  typedef CGAL::identity_property_map<typename std::iterator_traits<Iterator>::value_type > Identity;
  Identity identity;
  process_point_set(beg,end,identity);
}


// Here comes a function that changes the orientation and the normal
template <typename Iterator, typename PointPmap, typename OrientationPmap, typename NormalPmap >
void orient_normals(Iterator beg, Iterator end, PointPmap ppm, OrientationPmap opm, NormalPmap npm)
{
  for(;beg!= end;++beg){
    Vector_3& v = boost::get(npm,*beg);
    
    boost::put(opm, *beg, (v == CGAL::NULL_VECTOR));

    if(v.x() < 0){
      v = -v;
      boost::put(npm,*beg, v);
    }
  }
}


int main()
{
  CGAL::set_pretty_mode(std::cout);

  // an example with the data stored into different std::vectors
  {
    std::vector<Point_3> points;
    std::vector<std::pair<Vector_3,bool> > orientableNormals;
    for(int i = 0; i < 10; i++){
      double x = (i%2)?i:-i;
      points.push_back(Point_3(9-i,0,0));
      orientableNormals.push_back(std::make_pair(Vector_3(x,0,0),false));
    }
    
    CGAL::pointer_property_map_cluster<Point_3, Vector_3, bool> cluster(points.size());
    cluster.attach0(&points[0]);
    cluster.attach1(&orientableNormals[0].first, sizeof(std::pair<Vector_3,bool>));
    cluster.attach2(&orientableNormals[0].second, sizeof(std::pair<Vector_3,bool>));
    
    process_point_set(cluster.begin(),
                      cluster.end(),
                      cluster.pmap0());
    
    std::cout << "\n\nBefore orient_normals\n";
    for(int i = 0; i < 10; i++){
      std::cout << points[i] << " " << orientableNormals[i].first << " " << orientableNormals[i].second << std::endl;
    }
    
    orient_normals(cluster.begin(),
                   cluster.end(),
                   cluster.pmap0(),
                   cluster.pmap2(),
                   cluster.pmap1());
    
    std::cout << "\nAfter orient_normals\n";
    for(int i = 0; i < 10; i++){
      std::cout << points[i] << " " << orientableNormals[i].first << " " << orientableNormals[i].second << std::endl;
    } 
  }

  return 0;
}
