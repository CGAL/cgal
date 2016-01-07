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
  typedef typename boost::graph_traits<TM>::halfedge_iterator TM_halfedge_iterator;
  typedef typename boost::graph_traits<TM>::edge_descriptor TM_edge_descriptor;
  typedef typename boost::graph_traits<TM>::vertex_descriptor TM_vertex_descriptor;
  
public:

typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;

  const TM& mesh()const
  {
    return tm;
  }

  struct halfedge_descriptor;

  /**  A vertex 
   *
   **/
  struct vertex_descriptor{

    vertex_descriptor()
    {}

    vertex_descriptor(const halfedge_descriptor& h)
      :hd(h)
    {}

    vertex_descriptor(const vertex_descriptor& other)
      : hd(other.hd)
    {}


    bool operator ==(const vertex_descriptor& other) const
    {
      return (hd == other.hd);
    }


    bool operator !=(const vertex_descriptor& other) const
    {
      return (hd != other.hd);
    }


    bool operator<(const vertex_descriptor& other) const
    {
      return hd < other.hd;
    }


    operator TM_halfedge_descriptor() const
    {
      return hd;
    }

    friend std::ostream& operator<<(std::ostream& os, const vertex_descriptor vd)
    {
      os << "seam mesh vertex: " <<  vd.hd;
      return os;
    }

    TM_halfedge_descriptor hd;
  };


    class vertex_iterator
      : public boost::iterator_facade< vertex_iterator,
                                       vertex_descriptor,
                                       std::forward_iterator_tag
                                       >
    {
        typedef boost::iterator_facade< vertex_iterator,
                                        vertex_descriptor,
                                        std::forward_iterator_tag
                                        > Facade;

    public:
      vertex_iterator() : hd(), end(), mesh_(NULL) {}

      vertex_iterator(const Iterator_range<TM_halfedge_iterator>& ir, const Self* m)
        : hd(ir.first), end(ir.second), mesh_(m) {
        }

      // constructor for the past the end iterator 
      vertex_iterator(const TM_halfedge_iterator& hd, const Self* m)
        : hd(hd), end(hd), mesh_(m) {
        }

      vertex_iterator(const vertex_iterator& other)
        : hd(other.hd), end(other.end), mesh_(other.mesh_)
      {}

    private:
        friend class boost::iterator_core_access;

        void increment()
        {
          do{
            ++hd;
            if(hd == end) return;
            TM_vertex_descriptor tvd = target(*hd,mesh_->mesh());
            if( (! mesh_->has_on_seam(tvd))&& (halfedge(tvd,mesh_->mesh()) == *hd)) return;
            if(mesh_->has_on_seam(edge(*hd,mesh_->mesh()))) return;
          }while(true);
        }
  
        bool equal(const vertex_iterator& other) const
        {
          return (this->hd == other.hd) && (this->mesh_ == other.mesh_);
        }

        vertex_descriptor dereference() const { return vertex_descriptor(*hd); }


      TM_halfedge_iterator hd, end;
      const Self* mesh_;
    };




  bool has_on_seam(TM_vertex_descriptor vd) const
  {
    return seam_vertices.find(vd) != seam_vertices.end();
  }


  bool has_on_seam(TM_edge_descriptor ed) const
  {
    return seam_edges.find(ed) != seam_edges.end();
  }


  Iterator_range<vertex_iterator> m_vertices() const
  {
    Iterator_range<TM_halfedge_iterator> ir = halfedges(tm);
    vertex_iterator beg(ir,this);
    vertex_iterator end(ir.second,this);
  return make_range(beg,end);
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

  
  struct edge_descriptor {
    halfedge_descriptor hd;

    edge_descriptor(const halfedge_descriptor& hd)
      : hd(hd)
    {}

  };
  
  const TM& tm;
  std::set<TM_edge_descriptor> seam_edges;
  std::set<TM_vertex_descriptor> seam_vertices;

  int index;
  
public:
  template <typename EdgeRange>
  Seam_mesh(const TM& tm, EdgeRange er) 
    : tm(tm), seam_edges(er.begin(), er.end()), index(0)
    {}

  template <typename VertexIndexMap>
  void initialize_vertex_index_map(halfedge_descriptor bhd, VertexIndexMap& vipm)
  {
    Self& mesh=*this;
    // Initialize all indices with -1
    BOOST_FOREACH(TM_halfedge_descriptor hd, halfedges(tm)){
      put(vipm,halfedge_descriptor(hd),-1);
    }
    
    // Walk along the seam which may contain real border edges
    CGAL::Halfedge_around_face_circulator<Self> hafc(bhd,mesh), prev(hafc), done;
    ++hafc;
    done = hafc;
    do {
      halfedge_descriptor ohd = opposite(*hafc,mesh);
      TM_halfedge_descriptor tohd(ohd);
      seam_vertices.insert(target(tohd,tm));
      assert(! ohd.seam);
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(ohd,mesh)){
        if(! hd.seam){
          vertex_descriptor vd(hd);
          if(get(vipm,vd) == -1){
            put(vipm,vd,index);
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
      if(get(vipm,vertex_descriptor(halfedge_descriptor(hd))) == -1){
        BOOST_FOREACH(halfedge_descriptor hav, halfedges_around_target(hd,tm)){
          put(vipm,vertex_descriptor(halfedge_descriptor(hav)), index);
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
    return vertex_descriptor(hd.tmhd);
  }


  vertex_descriptor m_source(const halfedge_descriptor& hd) const
  {
    return vertex_descriptor(opposite(hd.tmhd, tm));
  }
};



} // namespace

#endif CGAL_SEAM_MESH_H
