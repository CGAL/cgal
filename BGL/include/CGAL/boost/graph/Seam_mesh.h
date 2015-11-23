#ifndef CGAL_SEAM_MESH_H
#define CGAL_SEAM_MESH_H

#include <set>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/graph_traits_Seam_mesh.h>

namespace CGAL {

template <class TM>
  class Seam_mesh {

  typedef Seam_mesh<TM> Self;
  typedef typename boost::graph_traits<TM>::halfedge_descriptor TM_halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::edge_descriptor TM_edge_descriptor;
  typedef typename boost::graph_traits<TM>::vertex_descriptor TM_vertex_descriptor;
  
public:
typedef typename boost::graph_traits<TM>::vertex_descriptor vertex_descriptor;

  const TM& mesh()const
  {
    return tm;
  }

  struct halfedge_descriptor {
    TM_halfedge_descriptor tmhd;
    bool seam;

    halfedge_descriptor()
      : tmhd(), seam(false)
    {}

    halfedge_descriptor(const halfedge_descriptor& other)
      : tmhd(other.tmhd), seam(other.seam)
    {}


    halfedge_descriptor(TM_halfedge_descriptor tmhd, bool seam=false)
      : tmhd(tmhd),seam(seam)
    {}


    bool operator ==(const halfedge_descriptor& other) const
    {
      return (tmhd == other.tmhd) && (seam == other.seam); 
    }


    bool operator !=(const halfedge_descriptor& other) const
    {
      return (tmhd != other.tmhd) || (seam != other.seam); 
    }


    bool operator<(const halfedge_descriptor& other) const
    {
      return tmhd < other.tmhd;
    }


    operator TM_halfedge_descriptor() const
    {
      return tmhd;
    }

    friend
    std::ostream& operator<<(std::ostream& os, const halfedge_descriptor& hd)
  {
    os << hd.tmhd  << ((hd.seam)?" on seam":"");
    return os;
  }
  };

  const TM& tm;
  std::set<TM_edge_descriptor> seam_edges;

  int index;
  
public:
  template <typename EdgeRange, typename HalfedgeAsVertexIndexMap>
  Seam_mesh(const TM& tm, EdgeRange er, typename boost::graph_traits<TM>::halfedge_descriptor smhd, HalfedgeAsVertexIndexMap hvipm)
    : tm(tm), seam_edges(er.begin(), er.end()), index(0)
  {
    Self& mesh=*this;
     // Initialize all indices with -1
  BOOST_FOREACH(TM_halfedge_descriptor hd, halfedges(tm)){
    put(hvipm,hd,-1);
  }

  halfedge_descriptor bhd(smhd);
  bhd = opposite(bhd,mesh);
  // Walk along the seam which may contain real border edges
  CGAL::Halfedge_around_face_circulator<Mesh> hafc(bhd,mesh), prev(hafc), done;
  ++hafc;
  done = hafc;
  do {
    halfedge_descriptor ohd = opposite(*hafc,mesh);
    assert(! ohd.seam);
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(ohd,mesh)){
      if(! hd.seam){
        TM_halfedge_descriptor shd(hd);
        if(get(hvipm,shd) == -1){
          put(hvipm,shd,index);
        }
      }
      if(hd == *prev){
        break;
      }
    }
    ++index;
    prev = hafc;
    ++hafc;
  }while(hafc != done);

 
  // now as all halfedges incident to seam vertices are handled
  // we look at the not yet marked halfedges
  
  BOOST_FOREACH(TM_halfedge_descriptor hd, halfedges(tm)){
    if(get(hvipm,hd) == -1){
      BOOST_FOREACH(halfedge_descriptor hav, halfedges_around_target(hd,tm)){
        put(hvipm,TM_halfedge_descriptor(hav),index);
      }
      ++index;
    }
  } 
  }

  // this is the number of different halfedge indices
  int m_num_vertices() const
  {
    return index;
  }


  bool is_on_seam(const halfedge_descriptor hd) const
  {
    return seam_edges.find(edge(hd, tm)) != seam_edges.end();
  }


  halfedge_descriptor m_next(const halfedge_descriptor& hd) const
  {
    if((! hd.seam)&& (! is_border(hd.tmhd,tm))){
      return halfedge_descriptor(next(hd.tmhd, tm));
    }
    Halfedge_around_target_circulator<TM> hatc(hd.tmhd,tm);
    do {
      --hatc;
    }while((! is_on_seam(*hatc))&&(! is_border(opposite(*hatc,tm),tm)));
    return halfedge_descriptor(opposite(*hatc,tm), ! is_border(opposite(*hatc,tm),tm));
  }
  

  halfedge_descriptor m_prev(const halfedge_descriptor& hd) const
  {
    if((! hd.seam)&& (! is_border(hd.tmhd,tm))){
      return halfedge_descriptor(prev(hd.tmhd, tm));
    }
    Halfedge_around_source_circulator<TM> hatc(hd.tmhd,tm);
    do {
      ++hatc;
    }while((! is_on_seam(*hatc))&&(! is_border(opposite(*hatc,tm),tm)));
    return halfedge_descriptor(opposite(*hatc,tm), ! is_border(opposite(*hatc,tm),tm));
  }
  

 halfedge_descriptor m_opposite(const halfedge_descriptor& hd) const
  {
    if(! hd.seam){
      return halfedge_descriptor(opposite(hd.tmhd,tm), is_on_seam(hd));
    }
    
    return halfedge_descriptor(opposite(hd.tmhd,tm));
  }


  vertex_descriptor m_target(const halfedge_descriptor& hd) const
  {
    return target(hd.tmhd, tm);
  }


  vertex_descriptor m_source(const halfedge_descriptor& hd) const
  {
    return source(hd.tmhd, tm);
  }
};



} // namespace

#endif CGAL_SEAM_MESH_H
