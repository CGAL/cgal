#ifndef CGAL_SM_ITEMS_H
#define CGAL_SM_ITEMS_H

#include <CGAL/basic.h>
#include <CGAL/In_place_list.h>
#include <CGAL/Object.h>
#include <string>
#include <strstream>

CGAL_BEGIN_NAMESPACE

template <typename K, typename M> class SM_items;
template <typename K, typename I> class Sphere_map;
template <typename SM> class SM_const_decorator;
template <typename SM> class SM_decorator;
template <typename EH> struct move_edge_around_svertex;
template <typename EH> struct move_edge_around_sface;

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
  class SVertex : public CGAL::In_place_list_base< SVertex<Refs> >
  { 
    typedef typename Refs::Items Items;
    friend class Sphere_map<Kernel_, Items>;
    friend class SM_const_decorator<Refs>;
    friend class SM_decorator<Refs>;

    typedef typename Refs::SVertex_handle   SVertex_handle;
    typedef typename Refs::SHalfedge_handle SHalfedge_handle;
    typedef typename Refs::SFace_handle     SFace_handle;

    Sphere_point    point_; 
    Mark            mark_;
    SHalfedge_handle out_sedge_;
    SFace_handle     incident_sface_;
    GenPtr          info_;
    // temporary information:

  public:
    SVertex() : 
      point_(), mark_(), out_sedge_(), incident_sface_(), info_() {}
    SVertex(const Mark& m) : 
      point_(), mark_(m), out_sedge_(), incident_sface_(), info_() {}
    SVertex(const Sphere_point& p) : 
      point_(p), mark_(), out_sedge_(), incident_sface_(), info_() {}

    ~SVertex() {}

    SVertex(const SVertex<Refs>& v)
    { point_ = v.point_;
      mark_ = v.mark_;
      out_sedge_ = v.out_sedge_;
      incident_sface_ = v.incident_sface_;
      info_ = 0;
    }

    SVertex<Refs>& operator=(const SVertex<Refs>& v)
    { point_ = v.point_;
      mark_ = v.mark_;
      out_sedge_ = v.out_sedge_;
      incident_sface_ = v.incident_sface_;
      info_ = 0;
      return *this;
    }
                          
    public:
    std::string debug() const
    { std::ostrstream os; set_pretty_mode(os);
      os<<"V"<<point_<<' '<<info_<<'\0';
      std::string res(os.str()); os.freeze(0); return res;
    }
 
  }; // SVertex    


  template <typename Refs>
  class SHalfedge : public CGAL::In_place_list_base< SHalfedge<Refs> >
  { 
    typedef typename Refs::Items Items;
    typedef typename Refs::SVertex_handle   SVertex_handle;
    typedef typename Refs::SHalfedge_handle SHalfedge_handle;
    typedef typename Refs::SHalfedge_const_handle SHalfedge_const_handle;
    typedef typename Refs::SFace_handle     SFace_handle;

    friend class Sphere_map<Kernel_, Items>;
    friend class SM_const_decorator<Refs>;
    friend class SM_decorator<Refs>;
    friend class move_edge_around_svertex<SHalfedge_handle>;
    friend class move_edge_around_sface<SHalfedge_handle>;
    friend class move_edge_around_svertex<SHalfedge_const_handle>;
    friend class move_edge_around_sface<SHalfedge_const_handle>;

    // Role within local graph:
    Sphere_circle     circle_;
    Mark              mark_;
    SHalfedge_handle   twin_, sprev_, snext_;
    SVertex_handle     source_;
    SFace_handle       incident_sface_;
    GenPtr            info_;

  public:
    SHalfedge() : circle_(), mark_(), twin_(), sprev_(), snext_(),
		 source_(), incident_sface_(), info_() {}

    ~SHalfedge() {}

    SHalfedge(const SHalfedge<Refs>& e)
    {
      circle_ = e.circle_;
      mark_ = e.mark_;
      twin_ = e.twin_;
      sprev_ = e.sprev_;
      snext_ = e.snext_;
      source_ = e.source_;
      incident_sface_ = e.incident_sface_;
      info_ = 0;
    }
    SHalfedge<Refs>& operator=(const SHalfedge<Refs>& e)
    {
      circle_ = e.circle_;
      mark_ = e.mark_;
      twin_ = e.twin_;
      sprev_ = e.sprev_;
      snext_ = e.snext_;
      source_ = e.source_;
      incident_sface_ = e.incident_sface_;
      info_ = 0;
      return *this;
    }

    bool is_twin() const { return (&*twin_ < this); }

    std::string debug() const
    { std::ostrstream os; set_pretty_mode(os); 
      os <<"e["<<source_->debug()<<", "
         <<twin_->source_->debug()<<" "<<info_<<"]"<<'\0';
      std::string res(os.str()); os.freeze(0); return res;
    }

  }; // SHalfedge

  template <typename Refs> 
  class SHalfloop : public CGAL::In_place_list_base< SHalfloop<Refs> >
  {
    typedef typename Refs::SHalfloop_handle SHalfloop_handle;
    typedef typename Refs::SFace_handle SFace_handle;
    typedef typename Refs::Items Items;
    friend class Sphere_map<Kernel_, Items>;
    friend class SM_const_decorator<Refs>;
    friend class SM_decorator<Refs>;
    friend class Self::SVertex<Refs>;

    Sphere_circle   circle_;
    Mark            mark_;
    SHalfloop_handle twin_;
    SFace_handle     incident_sface_;
    GenPtr          info_;
    // temporary needed:

  public:
    SHalfloop() : circle_(), mark_(), twin_(), incident_sface_(), info_() {}
    ~SHalfloop() {}
    SHalfloop(const SHalfloop<Refs>& l)
    { 
      circle_ = l.circle_;
      mark_ = l.mark_;
      twin_ = l.twin_;
      incident_sface_ = l.incident_sface_;
      info_ = 0;
    }
    SHalfloop<Refs>& operator=(const SHalfloop<Refs>& l)
    { 
      circle_ = l.circle_;
      mark_ = l.mark_;
      twin_ = l.twin_;
      incident_sface_ = l.incident_sface_;
      info_ = 0;
      return *this;
    }

    bool is_twin() const { return (&*twin_ < this); }

    std::string debug() const
    { std::ostrstream os; set_pretty_mode(os); 
      os<<"l"<<circle_<<' '<<info_<<'\0';
      std::string res(os.str()); os.freeze(0); return res;
    }

  }; // SHalfloop

  template <typename Refs>
  class SFace : public CGAL::In_place_list_base< SFace<Refs> >
  { 
    typedef typename Refs::Items Items;
    friend class Sphere_map<Kernel_, Items>;
    friend class SM_const_decorator<Refs>;
    friend class SM_decorator<Refs>;

    typedef typename Refs::Object_handle Object_handle;
    typedef typename Refs::Object_list Object_list;
    typedef typename Refs::SFace_cycle_iterator SFace_cycle_iterator;
    typedef typename Refs::SFace_cycle_const_iterator 
      SFace_cycle_const_iterator;

    Mark             mark_;
    Object_list      boundary_entry_objects_; 
    // SHalfedges, SHalfloops, Vertices
    GenPtr           info_;
    // temporary needed:

    public:
    SFace() : mark_(), info_() {}
    ~SFace() {}

    SFace(const SFace<Refs>& f)
    { mark_ = f.mark_;
      boundary_entry_objects_ = f.boundary_entry_objects_;
      info_ = 0;
    }
    SFace<Refs>& operator=(const SFace<Refs>& f)
    { if (this == &f) return *this;
      mark_ = f.mark_;
      boundary_entry_objects_ = f.boundary_entry_objects_;
      info_ = 0;
      return *this;
    }

    SFace_cycle_iterator sface_cycles_begin() 
    { return boundary_entry_objects_.begin(); }
    SFace_cycle_iterator sface_cycles_end()
    { return boundary_entry_objects_.end(); }

    SFace_cycle_const_iterator sface_cycles_begin() const
    { return boundary_entry_objects_.begin(); }
    SFace_cycle_const_iterator sface_cycles_end() const
    { return boundary_entry_objects_.end(); }

  }; // SFace

}; // SM_items


CGAL_END_NAMESPACE
#endif // CGAL_SM_ITEMS_H


