#ifndef SVD_INSERT_H
#define SVD_INSERT_H


template<class SVD, class Point>
typename SVD::Vertex_handle
insert_point(SVD& svd, const Point& p)
{
  return svd.insert(p);
}

template<class SVD, class Point>
typename SVD::Vertex_handle
insert_segment(SVD& svd, const Point& p1, const Point& p2)
{
  return svd.insert(p1, p2);
}

template<class SVD, class Point>
typename SVD::Vertex_handle
insert_segment(SVD& svd, const Point& p1, const Point& p2,
	       typename SVD::Vertex_handle v)
{
  return svd.insert(p1, p2, v);
}


template<class SVD, class Polygon>
typename SVD::Vertex_handle
insert_polygon(SVD& svd, const Polygon& pgn)
{
  typename SVD::Vertex_handle v;
  int psize = pgn.size();
  for (int i = 0; i < psize; i++ ) {
    v = svd.insert( pgn[i], pgn[(i+1)%psize] );
  }
  return v;
}



#endif // SVD_INSERT_H
