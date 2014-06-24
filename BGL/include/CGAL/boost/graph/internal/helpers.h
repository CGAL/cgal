#ifndef CGAL_BOOST_GRAPH_INTERNAL_HELPERS_H
#define CGAL_BOOST_GRAPH_INTERNAL_HELPERS_H

#include <boost/graph/graph_traits.hpp>
#include <boost/optional.hpp>
#include <CGAL/boost/graph/iterator.h>

namespace CGAL {

// breaks a dependency loop between <CGAL/boost/graph/helpers.h>
// and <CGAL/boost/graph/iterator.h>
template <typename Graph> class Halfedge_around_target_iterator;

namespace internal {

template <typename Graph>
void
set_border(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
          , Graph& g)
{
  set_face(h, boost::graph_traits<Graph>::null_face(), g);
}

template <typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
copy(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
                    , Graph& g)
{
  typename boost::graph_traits<Graph>::edge_descriptor e = add_edge(g);
  typename boost::graph_traits<Graph>::halfedge_descriptor res = halfedge(e,g);
  typename boost::graph_traits<Graph>::halfedge_descriptor ropp = opposite(res, g);
  typename boost::graph_traits<Graph>::halfedge_descriptor hopp = opposite(h, g);
  set_target(res, target(h, g), g);
  set_target(hopp, target(hopp, g), g);
  set_face(res, face(h, g), g);
  set_face(ropp, face(hopp, g), g);
  // note that we cannot call set_next as it then would call set_prev on the  original
  return res;
 }


template <typename Graph>
void
set_vertex_halfedge(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
                    , Graph& g)
{ set_halfedge(target(h, g), h, g); }


template <typename Graph>
void
close_tip(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
          , Graph& g)
{
  // makes `opposite(h,g)' the successor of h.
  set_next( h, opposite(h, g), g);
}


template <typename Graph>
void
close_tip(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
          , typename boost::graph_traits<Graph>::vertex_descriptor const& v
          , Graph& g)
{
  // makes `h->opposite()' the successor of h and sets the incident
  // vertex of h to v.
  set_next(h, opposite(h, g), g);
  set_target(h, v, g);
  set_halfedge(v, h, g);
}

template <typename Graph>
void
insert_tip(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
           , typename boost::graph_traits<Graph>::halfedge_descriptor const& h2
           , Graph& g)
{
  set_next(h, next(h2,g), g);
  set_next(h2, opposite(h, g), g);
  set_target(h, target(h2, g), g);
}


template <typename Graph>
void
remove_tip(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
           , Graph& g)
{ 
  set_next(h, next(opposite(next(h, g), g), g), g);
}


template <typename Graph>
void 
set_face_in_face_loop(typename boost::graph_traits<Graph>::halfedge_descriptor h, 
                      typename boost::graph_traits<Graph>::face_descriptor f, 
                      Graph& g) 
{
  typename boost::graph_traits<Graph>::halfedge_descriptor end = h;
  do {
    set_face(h, f, g);
    h = next(h, g);
  } while ( h != end);
}
    

template <typename Graph>
void insert_halfedge(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
                     , typename boost::graph_traits<Graph>::halfedge_descriptor const& f
                     , Graph& g)
{
  set_next(h, next(f, g), g);
  set_next(f, h, g);
  set_face(h, face(f, g), g);
}

template <typename Graph>
std::size_t
exact_num_vertices(const Graph& g)
{ 
  typename boost::graph_traits<Graph>::vertex_iterator beg, end;
  boost::tie(beg,end) = vertices(g);
  return std::distance(beg,end);
 }

template <typename Graph>
std::size_t
exact_num_halfedges(const Graph& g)
{ 
  typename boost::graph_traits<Graph>::halfedge_iterator beg, end;
  boost::tie(beg,end) = halfedges(g);
  return std::distance(beg,end);
 }

template <typename Graph>
std::size_t
exact_num_edges(const Graph& g)
{ 
  typename boost::graph_traits<Graph>::edge_iterator beg, end;
  boost::tie(beg,end) = edges(g);
  return std::distance(beg,end);
 }

template <typename Graph>
std::size_t
exact_num_faces(const Graph& g)
{ 
  typename boost::graph_traits<Graph>::face_iterator beg, end;
  boost::tie(beg,end) = faces(g);
  return std::distance(beg,end);
 }

  /*
template <typename Graph, typename BorderMap>
boost::optional<typename boost::graph_traits<Graph>::halfedge_descriptor>
is_border_vertex(const Graph& g
                 , typename boost::graph_traits<Graph>::vertex_descriptor v
                 , BorderMap bm)
{
  CGAL::Halfedge_around_target_iterator<Graph> havib, havie;
  for(boost::tie(havib, havie) = halfedges_around_target(halfedge(v, g), g); havib != havie; ++havib) {
    if(get(bm, *havib)) {
      typename boost::graph_traits<Graph>::halfedge_descriptor h = *havib;
      return h;
    }
  }
  // empty
  return boost::optional<typename boost::graph_traits<Graph>::halfedge_descriptor>();
}
  */


template <typename Graph>
bool is_valid(const Graph& g)
{
  return g.is_valid();
#if 0
  typedef boost::graph_traits<Graph> traits;

  if(num_halfedges(g) % 2 !=0) {
    return false;
  }

  bool valid = true;
  typename traits::halfedge_iterator it, end;

  for(boost::tie(it, end) = halfedges(g); it != end; ++it) {
    // should we introduce null_halfedge?
    // valid = valid && (next(*it) != traits::null_halfedge());
    // valid = valid && (opposite(*it) != traits::null_halfedge());
    // if(!valid) {
    //   std::cerr << "Integrity of halfedge " << *it << " corrupted."  << std::endl;
    //   break;
    // }
    valid = valid && (opposite(*it, g) != *it);
    valid = valid && (opposite(opposite(*it, g), g) == *it);
    if(!valid) {
      std::cerr << "Integrity of opposite halfedge of " << *it << " corrupted."  << std::endl;
      break;
    }

    valid = valid && (next(prev(*it, g), g) == *it);
    if(!valid) {
      std::cerr << "Integrity of previous halfedge of " << *it << " corrupted."  << std::endl;
      break;
    }

    valid = valid && (prev(next(*it, g), g) == *it);
    if(!valid) {
      std::cerr << "Integrity of next halfedge of " << *it << " corrupted."  << std::endl;
      break;
    }

    valid = valid && (target(*it, g) != traits::null_vertex());
    if(!valid) {
      std::cerr << "Integrity of vertex of halfedge " << *it << " corrupted."  << std::endl;
      break;
    }

    valid = valid && (target(*it, g) == target(opposite(next(*it, g), g), g));
    if(!valid) {
      std::cerr << "Halfedge vertex of next opposite is not the same for " << *it << "."  << std::endl;
      break;
    }
    
    valid = valid && (face(*it, g) == face(next(*it, g), g));
    if(!valid) {
      std::cerr << "Face of " << *it << " and face of next(*it, g) are not identical."  << std::endl;
      break;
    }
  }

  typename traits::face_iterator fit, fend;
  for(boost::tie(fit, fend) = faces(g); fit != fend; ++fit) {
    valid = valid && (halfedge(*fit, g) != traits::null_halfedge());
    if(!valid) {
      std::cerr << "Face " << *fit << " is associated to the null_halfedge." << std::endl;
      break;
    }

    valid = valid && (face(halfedge(*fit, g), g) == *fit);
    if(!valid) {
      std::cerr << "Face <-> Halfedge connectivity of " << *fit << " corrupted." << std::endl;
      break;
    }

    CGAL::Halfedge_around_face_iterator<Graph> hafit, hafiend;
    typename traits::halfedges_size_type n = 0;
    for(boost::tie(hafit, hafiend) = halfedges_around_face(halfedge(*fit, g), g); 
        hafit != hafiend; ++hafit, ++n) {
      if(n >= num_halfedges(g)) {
        valid = false;
        break;
      }
    }
    if(!valid) {
      std::cerr << "Too many halfedges around face " << *fit << std::endl;
      break;
    }
  }
  return valid;
#endif
}


} // internal
} // CGAL


#endif // CGAL_BOOST_GRAPH_INTERNAL_HELPERS_H
