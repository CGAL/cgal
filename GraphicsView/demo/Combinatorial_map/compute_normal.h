#ifndef _COMPUTE_NORMAL_
#define _COMPUTE_NORMAL_

template <class Map>
typename Map::Vector compute_facet_normal(Map& m, typename Map::Dart_handle d)
{
  typedef typename Map::Point Point;
  typedef typename Map::Vector Vector;
  Vector normal = CGAL::NULL_VECTOR;
  for (typename Map::Dart_iterator_of_beta1 it(m,d); it.cont(); ++it)
    {
      if (!it->is_free(0)&&!it->is_free(1))
	{
	  const Point& prev = it->beta(0)->vertex()->point();
	  const Point& curr = it->vertex()->point();
	  const Point& next = it->beta(1)->vertex()->point();
	  Vector n = CGAL::cross_product(next-curr,prev-curr);
	  normal = normal + (n / std::sqrt(n*n));
	}
    }
  return normal / std::sqrt(normal * normal);
}

template <class Map>
typename Map::Vector compute_vertex_normal(Map& m, typename Map::Dart_handle d)
{
  typedef typename Map::Point Point;
  typedef typename Map::Vector Vector;
  Vector normal = CGAL::NULL_VECTOR;

   for (typename Map::Dart_iterator_of_vertex it(m,d); it.cont(); ++it)
     {
       // ?? if(!he->is_border())
       {
	 Vector n = compute_facet_normal(m,*it);
	 normal = normal + (n / std::sqrt(n*n));
       }
     }
  return normal / std::sqrt(normal * normal);
}

#endif // _COMPUTE_NORMAL_
