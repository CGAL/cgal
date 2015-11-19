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
      if(tmhd < other.tmhd) return true;
      if(tmhd > other.tmhd) return false;
      if( (! seam) && other.seam) return true;
      return false;
    }


    operator TM_halfedge_descriptor() const
    {
      return tmhd;
    }
  };

  const TM& tm;
  std::set<TM_edge_descriptor> seam_edges;
public:
  template <typename EdgeRange>
  Seam_mesh(const TM& tm, EdgeRange er)
    : tm(tm), seam_edges(er.begin(), er.end())
  {
   
  }
 

  bool is_on_seam(const halfedge_descriptor hd) const
  {
    return seam_edges.find(edge(hd, tm)) != seam_edges.end();
  }


  halfedge_descriptor next(const halfedge_descriptor& hd) const
  {
    if((! hd.seam)&& (! is_border(hd.tmhd,tm))){
      return halfedge_descriptor(CGAL::next(hd.tmhd, tm));
    }
    Halfedge_around_target_circulator<TM> hatc(hd.tmhd,tm);
    do {
      --hatc;
    }while((! is_on_seam(*hatc))&&(! is_border(CGAL::opposite(*hatc,tm),tm)));
    return halfedge_descriptor(CGAL::opposite(*hatc,tm), ! is_border(CGAL::opposite(*hatc,tm),tm));
  }
  

  halfedge_descriptor prev(const halfedge_descriptor& hd) const
  {
    if((! hd.seam)&& (! is_border(hd.tmhd,tm))){
      return halfedge_descriptor(CGAL::prev(hd.tmhd, tm));
    }
    Halfedge_around_source_circulator<TM> hatc(hd.tmhd,tm);
    do {
      ++hatc;
    }while((! is_on_seam(*hatc))&&(! is_border(CGAL::opposite(*hatc,tm),tm)));
    return halfedge_descriptor(CGAL::opposite(*hatc,tm), ! is_border(CGAL::opposite(*hatc,tm),tm));
  }
  

 halfedge_descriptor opposite(const halfedge_descriptor& hd) const
  {
    if(! hd.seam){
      return halfedge_descriptor(CGAL::opposite(hd.tmhd,tm), is_on_seam(hd));
    }
    
    return halfedge_descriptor(CGAL::opposite(hd.tmhd,tm));
  }


  vertex_descriptor target(const halfedge_descriptor& hd) const
  {
    return CGAL::target(hd.tmhd, tm);
  }


  vertex_descriptor source(const halfedge_descriptor& hd) const
  {
    return CGAL::source(hd.tmhd, tm);
  }
};

} // namespace

#endif CGAL_SEAM_MESH_H
