#ifndef GEOMETRY_CONTAINER_H
#define GEOMETRY_CONTAINER_H
namespace CGAL{

template <typename PointRange, typename TAG>
struct Geometry_container{
  typedef std::pair<PointRange, TAG> type;
  typedef typename PointRange::difference_type difference_type; 

  PointRange range;
  Geometry_container() {}
  Geometry_container(PointRange& range)
    :range(range){}
  
  typename PointRange::iterator begin()
  { return range.begin(); }
  
  typename PointRange::iterator end()
  { return range.end(); }
  
  typename PointRange::const_iterator begin()const
  { return range.begin(); }
  
  typename PointRange::const_iterator end()const
  { return range.end(); }
  
  void clear(){ range.clear(); }
  
  template<typename size_t>
  void resize(size_t n){ range.resize(n); }
  
  template<typename T>
  void push_back(const T& t){ range.push_back(t); }
  
};//end Geometry_container
}//end CGAL

namespace boost{

template< class T, typename TAG >
struct range_iterator<CGAL::Geometry_container<T, TAG> >
{ typedef typename T::iterator type; };

template< class T, typename TAG >
struct range_iterator<const CGAL::Geometry_container<T, TAG> >
{ typedef typename T::const_iterator type; };

template< class T, typename TAG >
struct range_mutable_iterator<CGAL::Geometry_container<T, TAG> >
{ typedef typename range_mutable_iterator<T>::type type; };

}//end boost
#endif // GEOMETRY_CONTAINER_H
