#ifndef CGAL_BGL_DUAL_H
#define CGAL_BGL_DUAL_H

#include <CGAL/boost/graph/properties.h>
#include <boost/range/distance.hpp>

namespace CGAL {

template <typename Primal_>
class Dual
{
  
  const Primal_& primal_;

public:
  typedef Primal_ Primal;

  Dual(const Primal& primal)
  : primal_(primal)
  {}

  const Primal& primal() const
  {
    return primal_;
  }
};

} // namespace CGAL

namespace boost {
  
template <typename Primal>
class graph_traits<CGAL::Dual<Primal> >
{
  typedef boost::graph_traits<Primal> GTP;
  struct Dual_traversal_category : public virtual boost::bidirectional_graph_tag,
                                   public virtual boost::vertex_list_graph_tag,
                                   public virtual boost::edge_list_graph_tag
  {};

public:
  typedef typename GTP::face_descriptor     vertex_descriptor;
  typedef typename GTP::vertex_descriptor   face_descriptor;
  typedef typename GTP::halfedge_descriptor halfedge_descriptor;
  typedef typename GTP::edge_descriptor     edge_descriptor;
  typedef typename GTP::directed_category   directed_category;
  typedef boost::allow_parallel_edge_tag    edge_parallel_category; 
  typedef Dual_traversal_category           traversal_category;

  typedef typename GTP::faces_size_type          vertices_size_type;
  typedef typename GTP::vertices_size_type       faces_size_type;
  typedef typename GTP::edges_size_type          edges_size_type;
  typedef typename GTP::halfedges_size_type      halfedges_size_type;
  typedef typename GTP::degree_size_type         degree_size_type;

  typedef typename GTP::face_iterator     vertex_iterator;
  typedef typename GTP::vertex_iterator   face_iterator;
  typedef typename GTP::halfedge_iterator halfedge_iterator;
  typedef typename GTP::edge_iterator     edge_iterator;

  typedef CGAL::Halfedge_around_face_iterator<Primal> in_edge_iterator;
  typedef CGAL::Opposite_edge_around_face_iterator<Primal> out_edge_iterator;

  static vertex_descriptor   null_vertex()   { return vertex_descriptor(); }
  static face_descriptor     null_face()     { return face_descriptor(); }
  static halfedge_descriptor null_halfedge() { return halfedge_descriptor(); }
};
 
template<typename P>
struct graph_traits< const CGAL::Dual<P> >  
  : public graph_traits< CGAL::Dual<P> >
{}; 

template <typename P>
struct property_map<CGAL::Dual<P>, boost::vertex_index_t>
{
  typedef typename property_map<P, CGAL::face_index_t>::type type; 
  typedef typename property_map<P, CGAL::face_index_t>::const_type const_type; 
};

} // namespace boost


namespace CGAL {

template <typename P>
typename boost::property_map<P, boost::face_index_t>::type
get(boost::vertex_index_t, const Dual<P>& dual)
{
  return get(CGAL::face_index, dual.primal());
}


template <typename P>
typename boost::graph_traits<CGAL::Dual<P> >::vertices_size_type
num_vertices(const CGAL::Dual<P>& dual)
{
  return num_faces(dual.primal());
}
     
template <typename P>
Iterator_range<typename boost::graph_traits<Dual<P> >::vertex_iterator>
vertices(const CGAL::Dual<P>& dual)
{
  return faces(dual.primal()); 
}
    
template <typename P>
Iterator_range<typename boost::graph_traits<Dual<P> >::halfedge_iterator>
halfedges(const CGAL::Dual<P>& dual)
{
  return halfedges(dual.primal()); 
}
  
template <typename P>
Iterator_range<typename boost::graph_traits<Dual<P> >::edge_iterator>
edges(const CGAL::Dual<P>& dual)
{
  return edges(dual.primal()); 
}

template <typename P>
typename boost::graph_traits<Dual<P> >::vertex_descriptor
source(typename boost::graph_traits<Dual<P> >::halfedge_descriptor h,
       const Dual<P>& dual)
{
  const Dual<P>::Primal& primal = dual.primal();
  return face(opposite(h,primal),primal);
}
 
template <typename P>
typename boost::graph_traits<Dual<P> >::vertex_descriptor
target(typename boost::graph_traits<Dual<P> >::halfedge_descriptor h,
       const Dual<P>& dual)
{
  const Dual<P>::Primal& primal = dual.primal();
  return face(h,primal);
}
 

template <typename P>
typename boost::graph_traits<Dual<P> >::vertex_descriptor
source(typename boost::graph_traits<Dual<P> >::edge_descriptor h,
       const Dual<P>& dual)
{
  const Dual<P>::Primal& primal = dual.primal();
  return face(opposite(halfedge(h,primal),primal),primal);
}
 
template <typename P>
typename boost::graph_traits<Dual<P> >::vertex_descriptor
target(typename boost::graph_traits<Dual<P> >::edge_descriptor h,
       const Dual<P>& dual)
{
  const Dual<P>::Primal& primal = dual.primal();
  return face(halfedge(h,primal),primal);
}

template <typename P>
typename boost::graph_traits<Dual<P> >::halfedge_descriptor
halfedge(typename boost::graph_traits<Dual<P> >::vertex_descriptor v,
         const Dual<P>& dual)
{
  return halfedge(v, dual.primal());
}

template <typename P>
typename boost::graph_traits<Dual<P> >::halfedge_descriptor
halfedge(typename boost::graph_traits<Dual<P> >::face_descriptor f,
         const Dual<P>& dual)
{
  return halfedge(f, dual.primal());
}

template <typename P>
typename boost::graph_traits<Dual<P> >::halfedge_descriptor
halfedge(typename boost::graph_traits<Dual<P> >::edge_descriptor e,
         const Dual<P>& dual)
{
  return halfedge(e, dual.primal());
}
template <typename P>
typename boost::graph_traits<Dual<P> >::face_descriptor
face(typename boost::graph_traits<Dual<P> >::halfedge_descriptor h,
       const Dual<P>& dual)
{
  const Dual<P>::Primal& primal = dual.primal();
  return target(h,primal);
}

template <typename P>
typename boost::graph_traits<Dual<P> >::halfedge_descriptor
opposite(typename boost::graph_traits<Dual<P> >::halfedge_descriptor h,
         const Dual<P>& dual)
{
  return opposite(h, dual.primal());
}

template <typename P>
typename boost::graph_traits<Dual<P> >::halfedge_descriptor
next(typename boost::graph_traits<Dual<P> >::halfedge_descriptor h,
         const Dual<P>& dual)
{
  const Dual<P>::Primal& primal = dual.primal();
  return opposite(prev(h,primal),primal);
}  

template <typename P>
typename boost::graph_traits<Dual<P> >::halfedge_descriptor
prev(typename boost::graph_traits<Dual<P> >::halfedge_descriptor h,
         const Dual<P>& dual)
{
  const Dual<P>::Primal& primal = dual.primal();
  return next(opposite(h,primal),primal);
}  

template <typename P>
Iterator_range<typename boost::graph_traits<Dual<P> >::out_edge_iterator>
out_edges(typename boost::graph_traits<Dual<P> >::vertex_descriptor v,
          const Dual<P>& dual)
{
  typename const Dual<P>::Primal& primal = dual.primal();
  return opposite_edges_around_face(halfedge(v,primal),primal);
}

template <typename P>
Iterator_range<typename boost::graph_traits<Dual<P> >::out_edge_iterator>
in_edges(typename boost::graph_traits<Dual<P> >::vertex_descriptor v,
         const Dual<P>& dual)
{
  typename const Dual<P>::Primal& primal = dual.primal();
  return halfedges_around_face(halfedge(v,primal),primal);
}
       
template <typename P>
typename boost::graph_traits<Dual<P> >::degree_size_type
out_degree(typename boost::graph_traits<Dual<P> >::vertex_descriptor v,
           const Dual<P>& dual)
{
  typename const Dual<P>::Primal& primal = dual.primal();
  return boost::distance(halfedges_around_face(halfedge(v,primal),primal));
}

 template <typename P>
typename boost::graph_traits<Dual<P> >::degree_size_type
in_degree(typename boost::graph_traits<Dual<P> >::vertex_descriptor v,
           const Dual<P>& dual)
{
  return out_degree(v,dual);
}
         
        
} // namespace CGAL

#endif
