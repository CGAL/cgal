#ifndef CGAL_SM_ITEMS_H
#define CGAL_SM_ITEMS_H

#include <CGAL/basic.h>
#include <CGAL/In_place_list.h>
#include <CGAL/Object.h>
#include <string>
#include <strstream>

CGAL_BEGIN_NAMESPACE

template <typename K, typename M> class SM_items;
template <typename K> class Sphere_map;
template <typename SM, typename K> class SM_const_decorator;
template <typename SM, typename K> class SM_decorator;
template <typename EH> struct move_edge_around_vertex;
template <typename EH> struct move_edge_around_face;

template <typename Kernel_, typename Mark_>
struct SM_items {
public:
  typedef Kernel_                         Kernel;
  typedef typename Kernel::Sphere_point   Sphere_point;
  typedef typename Kernel::Sphere_segment Sphere_segment;
  typedef typename Kernel::Sphere_circle  Sphere_circle;
  typedef Mark_                           Mark;
  typedef SM_items<Kernel_,Mark_>         Self;
  typedef void*                           GenPtr;

  template <typename Refs>
  class Vertex : public CGAL::In_place_list_base< Vertex<Refs> >
  { 
    typedef typename Refs::Items Items;
    friend class Sphere_map<Kernel_>;
    friend class SM_const_decorator<Refs,Kernel_>;
    friend class SM_decorator<Refs,Kernel_>;

    typedef typename Refs::Vertex_handle   Vertex_handle;
    typedef typename Refs::Halfedge_handle Halfedge_handle;
    typedef typename Refs::Face_handle     Face_handle;

    Sphere_point    point_; 
    Mark            mark_;
    Halfedge_handle edge_;
    Face_handle     face_;
    GenPtr          info_;
    // temporary information:

  public:
    Vertex() : 
      point_(), mark_(), edge_(), face_(), info_() {}
    Vertex(const Mark& m) : 
      point_(), mark_(m), edge_(), face_(), info_() {}
    Vertex(const Sphere_point& p) : 
      point_(p), mark_(), edge_(), face_(), info_() {}

    ~Vertex() {}

    Vertex(const Vertex<Refs>& v)
    { point_ = v.point_;
      mark_ = v.mark_;
      edge_ = v.edge_;
      face_ = v.face_;
      info_ = 0;
    }

    Vertex<Refs>& operator=(const Vertex<Refs>& v)
    { point_ = v.point_;
      mark_ = v.mark_;
      edge_ = v.edge_;
      face_ = v.face_;
      info_ = 0;
      return *this;
    }
                          
    public:
    std::string debug() const
    { std::ostrstream os; set_pretty_mode(os);
      os<<"V"<<point_<<' '<<info_<<'\0';
      std::string res(os.str()); os.freeze(0); return res;
    }
 
  }; // Vertex    


  template <typename Refs>
  class Halfedge : public CGAL::In_place_list_base< Halfedge<Refs> >
  { 
    typedef typename Refs::Items Items;
    typedef typename Refs::Vertex_handle   Vertex_handle;
    typedef typename Refs::Halfedge_handle Halfedge_handle;
    typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Refs::Face_handle     Face_handle;

    friend class Sphere_map<Kernel_>;
    friend class SM_const_decorator<Refs,Kernel_>;
    friend class SM_decorator<Refs,Kernel_>;
    friend class move_edge_around_vertex<Halfedge_handle>;
    friend class move_edge_around_face<Halfedge_handle>;
    friend class move_edge_around_vertex<Halfedge_const_handle>;
    friend class move_edge_around_face<Halfedge_const_handle>;

    // Role within local graph:
    Sphere_circle     circle_;
    Mark              mark_;
    Halfedge_handle   twin_, prev_, next_;
    Vertex_handle     source_;
    Face_handle       face_;
    GenPtr            info_;

  public:
    Halfedge() : circle_(), mark_(), twin_(), prev_(), next_(),
		 source_(), face_(), info_() {}

    ~Halfedge() {}

    Halfedge(const Halfedge<Refs>& e)
    {
      circle_ = e.circle_;
      mark_ = e.mark_;
      twin_ = e.twin_;
      prev_ = e.prev_;
      next_ = e.next_;
      source_ = e.source_;
      face_ = e.face_;
      info_ = 0;
    }
    Halfedge<Refs>& operator=(const Halfedge<Refs>& e)
    {
      circle_ = e.circle_;
      mark_ = e.mark_;
      twin_ = e.twin_;
      prev_ = e.prev_;
      next_ = e.next_;
      source_ = e.source_;
      face_ = e.face_;
      info_ = 0;
      return *this;
    }

    std::string debug() const
    { std::ostrstream os; set_pretty_mode(os); 
      os <<"e["<<source_->debug()<<", "
         <<twin_->source_->debug()<<" "<<info_<<"]"<<'\0';
      std::string res(os.str()); os.freeze(0); return res;
    }

  }; // Halfedge

  template <typename Refs> 
  class Halfloop : public CGAL::In_place_list_base< Halfloop<Refs> >
  {
    typedef typename Refs::Halfloop_handle Halfloop_handle;
    typedef typename Refs::Face_handle Face_handle;
    typedef typename Refs::Items Items;
    friend class Sphere_map<Kernel_>;
    friend class SM_const_decorator<Refs,Kernel_>;
    friend class SM_decorator<Refs,Kernel_>;
    friend class Self::Vertex<Refs>;

    Sphere_circle   circle_;
    Mark            mark_;
    Halfloop_handle twin_;
    Face_handle     face_;
    GenPtr          info_;
    // temporary needed:

  public:
    Halfloop() : circle_(), mark_(), twin_(), face_(), info_() {}
    ~Halfloop() {}
    Halfloop(const Halfloop<Refs>& l)
    { 
      circle_ = l.circle_;
      mark_ = l.mark_;
      twin_ = l.twin_;
      face_ = l.face_;
      info_ = 0;
    }
    Halfloop<Refs>& operator=(const Halfloop<Refs>& l)
    { 
      circle_ = l.circle_;
      mark_ = l.mark_;
      twin_ = l.twin_;
      face_ = l.face_;
      info_ = 0;
      return *this;
    }

    std::string debug() const
    { std::ostrstream os; set_pretty_mode(os); 
      os<<"l"<<circle_<<' '<<info_<<'\0';
      std::string res(os.str()); os.freeze(0); return res;
    }

  }; // Halfloop

  template <typename Refs>
  class Face : public CGAL::In_place_list_base< Face<Refs> >
  { 
    typedef typename Refs::Items Items;
    friend class Sphere_map<Kernel_>;
    friend class SM_const_decorator<Refs,Kernel_>;
    friend class SM_decorator<Refs,Kernel_>;

    typedef typename Refs::Object_handle Object_handle;
    typedef typename Refs::Object_list Object_list;
    typedef typename Refs::Face_cycle_iterator Face_cycle_iterator;
    typedef typename Refs::Face_cycle_const_iterator 
      Face_cycle_const_iterator;

    Mark             mark_;
    Object_list      boundary_; // Halfedges, Halfloops, Vertices
    GenPtr           info_;
    // temporary needed:

    public:
    Face() : mark_(), info_() {}
    ~Face() {}

    Face(const Face<Refs>& f)
    { mark_ = f.mark_;
      boundary_ = f.boundary_;
      info_ = 0;
    }
    Face<Refs>& operator=(const Face<Refs>& f)
    { if (this == &f) return *this;
      mark_ = f.mark_;
      boundary_ = f.boundary_;
      info_ = 0;
      return *this;
    }

    Face_cycle_iterator face_cycles_begin() 
    { return boundary_.begin(); }
    Face_cycle_iterator face_cycles_end()
    { return boundary_.end(); }

    Face_cycle_const_iterator face_cycles_begin() const
    { return boundary_.begin(); }
    Face_cycle_const_iterator face_cycles_end() const
    { return boundary_.end(); }

  }; // Face

}; // SM_items


CGAL_END_NAMESPACE
#endif // CGAL_SM_ITEMS_H


