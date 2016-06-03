#ifndef CGAL_SEAM_MESH_H
#define CGAL_SEAM_MESH_H

#include <boost/unordered_set.hpp>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/graph_traits_Seam_mesh.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

namespace CGAL {

/*!
\ingroup PkgBGLHelper

The class template `Seam_mesh` is an adaptor for a triangle mesh and some marked edges
which look like boundary edges, when exploring the seam mesh with the generic BGL style
functions, and `boost::graph_traits<Seam_mesh<TM> >`. 
*/
template <class TM>
  class Seam_mesh {

private:
  typedef Seam_mesh<TM> Self;
  typedef typename boost::graph_traits<TM>::halfedge_descriptor TM_halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::halfedge_iterator TM_halfedge_iterator;
  typedef typename boost::graph_traits<TM>::edge_descriptor TM_edge_descriptor;
  typedef typename boost::graph_traits<TM>::vertex_descriptor TM_vertex_descriptor;
  
public:
/// @cond CGAL_DOCUMENT_INTERNALS
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

   friend std::size_t hash_value(const vertex_descriptor&  vd)
  {
    return hash_value(vd.hd.tmhd);
  }

    halfedge_descriptor hd;
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
        : hd(ir.first), end(ir.second), mesh_(m)
      {
        //std::cerr << "vertex_iterator(..)\n";
        //std::cerr << *hd << std::endl;
        if(hd == end) return;
        TM_vertex_descriptor tvd = target(*hd,mesh_->mesh());
        if( (! mesh_->has_on_seam(tvd))&& (halfedge(tvd,mesh_->mesh()) == *hd)) return;
        if(mesh_->has_on_seam(edge(*hd,mesh_->mesh()))) return;
        if(mesh_->has_on_seam(tvd) && is_border(opposite(*hd,mesh_->mesh()),mesh_->mesh())) return;
        increment();
        //std::cerr << *hd << "  after increment" << std::endl;
        //std::cerr << "leave vertex_iterator(..)\n";
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
          //std::cerr << "increment\n";
        if(hd == end) return;
          do{
            ++hd;
            //std::cerr << *hd << "  ++" << std::endl;
            if(hd == end) return;
            TM_vertex_descriptor tvd = target(*hd,mesh_->mesh());
            //std::cerr << "tvd = " << tvd << std::endl;
            if( (! mesh_->has_on_seam(tvd))&& (halfedge(tvd,mesh_->mesh()) == *hd)){
              //std::cerr <<"return as not on seam and reverse incidence\n";
              return;
            }
            if(mesh_->has_on_seam(edge(*hd,mesh_->mesh()))){
              //std::cerr <<"return as edge on seam\n";
              return;
            }
            if(mesh_->has_on_seam(tvd) && is_border(opposite(*hd,mesh_->mesh()),mesh_->mesh())){
              //std::cerr <<"return as edge on border and target on seam\n";
              return;
            }
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
  boost::unordered_set<TM_edge_descriptor> seam_edges;
  boost::unordered_set<TM_vertex_descriptor> seam_vertices;

  int index;
  
  /// @endcond

public:


  /// Constructs a seam mesh for a triangle mesh and a range of edges of the triangle mesh.
  template <typename EdgeRange>
  Seam_mesh(const TM& tm, EdgeRange er) 
    : tm(tm), seam_edges(er.begin(), er.end()), index(0)
    {
      BOOST_FOREACH(TM_edge_descriptor ed, seam_edges){
        seam_vertices.insert(source(ed,tm));
        seam_vertices.insert(target(ed,tm));
      }
    }

  /// Sets indices to 0,1,2,... for vertices in the connected component with the boundary on which lies `bhd`.
  /// The values are written into a property map with keytype `vertex_dscriptor` and
  /// value type `boost::graph_traits<TM>::vertices_size_type`.
  /// 
  /// The highest index is cached and returned when calling `num_vertices(sm)`.
  template <typename VertexIndexMap>
  void initialize_vertex_index_map(halfedge_descriptor bhd, VertexIndexMap& vipm)
  {
    Self& mesh=*this;
   
    index = 0;
    std::vector<face_descriptor> faces;
    boost::graph_traits<Seam_mesh>::halfedge_descriptor shd(opposite(bhd,*this));
    CGAL::Polygon_mesh_processing::connected_component(face(shd,*this),
                                                       *this,
                                                       std::back_inserter(faces));

    BOOST_FOREACH(face_descriptor fd, faces){
      BOOST_FOREACH(TM_halfedge_descriptor tmhd , halfedges_around_face(halfedge(fd,tm),tm)){
        halfedge_descriptor hd(tmhd);
        vertex_descriptor vd = target(hd,mesh);
          put(vipm,vd,-1);
      }
    }
    BOOST_FOREACH(face_descriptor fd, faces){
      BOOST_FOREACH(TM_halfedge_descriptor tmhd , halfedges_around_face(halfedge(fd,tm),tm)){
        halfedge_descriptor hd(tmhd);
        vertex_descriptor vd = target(hd,mesh);
        if(get(vipm,vd) == -1){
          put(vipm,vd,index);
          ++index;
        }
      }
    }
  }

/// @cond CGAL_DOCUMENT_INTERNALS

  // this is the number of different halfedge indices
  int m_num_vertices() const
  {
    return index;
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


  vertex_descriptor m_target(halfedge_descriptor hd) const
  {
    TM_halfedge_descriptor tmhd(hd);
    if(! has_on_seam(target(tmhd,tm))){
      tmhd = halfedge(target(tmhd,tm),tm);
      return vertex_descriptor(halfedge_descriptor(tmhd));
    }

    if(hd.seam){
      return m_target(halfedge_descriptor(prev(opposite(tmhd,tm),tm)));
    }

    while((! has_on_seam(tmhd)) && (! is_border(opposite(tmhd,tm),tm))){
      tmhd = prev(opposite(tmhd,tm),tm);
    }
    
  return vertex_descriptor(halfedge_descriptor(tmhd));
}

  /*
  vertex_descriptor m_target(const halfedge_descriptor& hd) const
  {
    return m_target(hd.tmhd);
  }
  */

  vertex_descriptor m_source(const halfedge_descriptor& hd) const
  {
    return m_target(opposite(hd.tmhd, tm));
  }

  /// @endcond

  /// returns `true` if the halfedge is on the seam.
  bool is_on_seam(const halfedge_descriptor& hd) const
  {
    return seam_edges.find(edge(hd, tm)) != seam_edges.end();
  }

};



} // namespace

#endif CGAL_SEAM_MESH_H
