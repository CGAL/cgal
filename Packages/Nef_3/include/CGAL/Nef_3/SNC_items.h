// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : include/CGAL/Nef_3/SNC_items.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// maintainer    : Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// SNC_items.h             items classes (Vertex, Halfedge, Facet, Volume objects)
// ============================================================================
#ifndef CGAL_SNC_ITEMS_H
#define CGAL_SNC_ITEMS_H

#include <CGAL/basic.h>
#include <CGAL/In_place_list.h>
#include <CGAL/Nef_2/Object_handle.h>
#include <string>

CGAL_BEGIN_NAMESPACE

template <typename K, typename M> class SNC_items;
template <typename I> class SNC_structure;
template <typename R> class SNC_decorator;
template <typename R> class SNC_io_parser;
template <typename R> class SNC_SM_decorator;
template <typename R> class SNC_SM_const_decorator;
template <typename R> class SNC_FM_decorator;
template <typename EH>
  struct move_shalfedge_around_svertex;
template <typename EH>
  struct move_shalfedge_around_sface;
template <typename HE>
  struct move_shalfedge_around_facet; 

template <typename Kernel_, typename Mark_>
class SNC_items {
public:
  typedef SNC_items<Kernel_,Mark_>  Self;
  typedef Kernel_                   Kernel;
  typedef typename Kernel::Point_3  Point_3;
  typedef typename Kernel::Plane_3  Plane_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef CGAL::Sphere_geometry<Kernel> Sphere_kernel;
  typedef typename Sphere_kernel::Sphere_point     Sphere_point;
  typedef typename Sphere_kernel::Sphere_segment   Sphere_segment;
  typedef typename Sphere_kernel::Sphere_circle    Sphere_circle;
  typedef typename Sphere_kernel::Sphere_direction Sphere_direction;
  typedef typename Sphere_kernel::Sphere_triangle  Sphere_triangle;
  typedef Mark_                                    Mark;
  typedef void*                                    GenPtr;

//-----------------------------------------------------------------------------

  template <typename Refs>
  class Vertex : public CGAL::In_place_list_base< Vertex<Refs> >
  {
    typedef typename Refs::Items Items;
    friend class SNC_structure<Items>;
    friend class SNC_decorator<Refs>;
    friend class SNC_io_parser<Refs>;
    friend class SNC_SM_decorator<Refs>;
    friend class SNC_SM_const_decorator<Refs>;
    typedef typename Refs::SVertex_iterator SVertex_iterator;
    typedef typename Refs::SHalfedge_iterator SHalfedge_iterator;
    typedef typename Refs::SHalfloop_iterator SHalfloop_iterator;
    typedef typename Refs::SFace_iterator SFace_iterator;

    typedef typename Refs::Vertex_iterator Vertex_iterator;
    typedef typename Refs::SHalfloop_handle SHalfloop_handle;
    typedef typename Refs::SHalfloop SHalfloop;
    typedef typename Refs::SVertex_const_iterator SVertex_const_iterator;
    typedef typename Refs::SHalfedge_const_iterator SHalfedge_const_iterator;
    typedef typename Refs::SFace_const_iterator SFace_const_iterator;
    typedef typename Refs::SHalfloop_const_handle SHalfloop_const_handle;
    typedef typename Refs::Vertex_const_iterator Vertex_const_iterator;

    Point_3            point_at_center_;
    Mark               mark_;
    // local view (surface graph):
    Refs*              sncp_;
    SVertex_iterator   svertices_begin_, svertices_last_;
    SHalfedge_iterator shalfedges_begin_, shalfedges_last_;
    SFace_iterator     sfaces_begin_, sfaces_last_;
    SHalfloop_iterator shalfloop_;
    Mark               m_neg_, m_pos_;
    GenPtr             info_;

  public:

    Vertex() : point_at_center_(), mark_(), sncp_(), 
      svertices_begin_(), svertices_last_(),
      shalfedges_begin_(), shalfedges_last_(),
      sfaces_begin_(), sfaces_last_(), shalfloop_(),
      m_neg_(), m_pos_(), info_() {}

    Vertex(const Point_3& p, Mark m) : 
      point_at_center_(p), mark_(m), sncp_(), 
      svertices_begin_(), svertices_last_(),
      shalfedges_begin_(), shalfedges_last_(),
      sfaces_begin_(), sfaces_last_(), shalfloop_(),
      m_neg_(), m_pos_() , info_() {}

    Vertex(const Vertex<Refs>& v) 
    { 
      point_at_center_ = v.point_at_center_;
      mark_ = v.mark_;
      sncp_ = v.sncp_;
      svertices_begin_ = v.svertices_begin_;
      svertices_last_ = v.svertices_last_;
      shalfedges_begin_ = v.shalfedges_begin_;
      shalfedges_last_ = v.shalfedges_last_;
      sfaces_begin_ = v.sfaces_begin_;
      sfaces_last_ = v.sfaces_last_;
      shalfloop_ = v.shalfloop_;
      m_neg_ = v.m_neg_; m_pos_ = v.m_pos_;
      info_ = 0;
    }
      
    Vertex<Refs>& operator=(const Vertex<Refs>& v)
    { if (this == &v) return *this;
      point_at_center_ = v.point_at_center_;
      mark_ = v.mark_;
      sncp_ = v.sncp_;
      svertices_begin_ = v.svertices_begin_;
      svertices_last_ = v.svertices_last_;
      shalfedges_begin_ = v.shalfedges_begin_;
      shalfedges_last_ = v.shalfedges_last_;
      sfaces_begin_ = v.sfaces_begin_;
      sfaces_last_ = v.sfaces_last_;
      shalfloop_ = v.shalfloop_;
      m_neg_ = v.m_neg_; m_pos_ = v.m_pos_; 
      return *this;
    }

    Refs* sncp() const { return sncp_; }

    /* all sobjects of the local graph are stored in a global list
       where each vertex has a continous range in each list for its
       sobjects. All objects of the range [sxxx_begin_,sxxx_last_]
       belong to a vertex. This range is empty iff 
       sxxx_begin_ == sxxx_last_ == sncp()->sxxx_end()
       ( the latter being the standard end iterator of the 
         corresponding list )
       for the past the end iterator we have to increment sxxx_last_
       once iff the range is non-empty. */

    void init_range(SVertex_iterator it)
    { svertices_begin_ = svertices_last_ = it; }
    void init_range(SHalfedge_iterator it)
    { shalfedges_begin_ = shalfedges_last_ = it; }
    void init_range(SFace_iterator it)
    { sfaces_begin_ = sfaces_last_ = it; }

    SVertex_iterator svertices_begin() 
    { return svertices_begin_; }
    SVertex_iterator svertices_last() 
    { return svertices_last_; }
    SVertex_iterator svertices_end() 
    { if ( svertices_last_ == sncp()->svertices_end() ) 
        return svertices_last_; 
      else 
	return ++SVertex_iterator(svertices_last_); }

    SHalfedge_iterator shalfedges_begin() 
    { return shalfedges_begin_; }
    SHalfedge_iterator shalfedges_last() 
    { return shalfedges_last_; }
    SHalfedge_iterator shalfedges_end() 
    { if ( shalfedges_last_ == sncp()->shalfedges_end() ) 
        return shalfedges_last_;
      else 
	return ++SHalfedge_iterator(shalfedges_last_); }

    SFace_iterator sfaces_begin() 
    { return sfaces_begin_; }
    SFace_iterator sfaces_last()
    { return sfaces_last_; }
    SFace_iterator sfaces_end()
    { if ( sfaces_last_ == sncp()->sfaces_end() ) 
        return sfaces_last_;
      else 
	return ++SFace_iterator(sfaces_last_); }

    SHalfloop_handle shalfloop() 
    { return shalfloop_; }

    SVertex_const_iterator svertices_begin() const
    { return svertices_begin_; }
    SVertex_const_iterator svertices_last() const
    { return svertices_last_; }
    SVertex_const_iterator svertices_end() const
    { if ( svertices_last_ == sncp()->svertices_end() ) 
        return svertices_last_; 
      else 
	return ++SVertex_const_iterator(svertices_last_); }

    SHalfedge_const_iterator shalfedges_begin() const
    { return shalfedges_begin_; }
    SHalfedge_const_iterator shalfedges_last() const
    { return shalfedges_last_; }
    SHalfedge_const_iterator shalfedges_end() const
    { if ( shalfedges_last_ == sncp()->shalfedges_end() ) 
        return shalfedges_last_;
      else 
	return ++SHalfedge_const_iterator(shalfedges_last_); }

    SFace_const_iterator sfaces_begin() const
    { return sfaces_begin_; }
    SFace_const_iterator sfaces_last() const
    { return sfaces_last_; }
    SFace_const_iterator sfaces_end() const
    { if ( sfaces_last_ == sncp()->sfaces_end() ) 
        return sfaces_last_;
      else 
	return ++SFace_const_iterator(sfaces_last_); }

    SHalfloop_const_handle shalfloop() const
    { return shalfloop_; }


    /*{\Xtext Vertices provide access to their local graphs via
    the iterator ranges:
    \begin{Mverb}
    SVertex_iterator     svertices_begin()/svertices_end()
    SHalfedge_iterator   shalfedges_begin()/shalfedges_end()
    SFace_iterator       sfaces_begin()/sfaces_end() 
    SHalfloop_handle     shalfloop()
    \end{Mverb}
    }*/

    void clear_local_graph() 
    /*{\Xop clears the local graph.}*/
    { if ( shalfloop() != sncp()->shalfloops_end() ) {
        sncp()->delete_shalfloop_only(shalfloop_->twin_);
        sncp()->delete_shalfloop_only(shalfloop_);
        shalfloop_ = sncp()->shalfloops_end();
      }
      SVertex_iterator vit = svertices_begin(),
                       vend = svertices_end();
      while (vit != vend) {
        SVertex_iterator vdel = vit++;
        sncp()->delete_svertex_only(vdel);
      }
      svertices_begin_ = svertices_last_ = sncp()->halfedges_end();

      SHalfedge_iterator eit = shalfedges_begin(),
                         eend = shalfedges_end();
      while (eit != eend) {
        SHalfedge_iterator edel = eit++;
        sncp()->delete_shalfedge_only(edel);
      }
      shalfedges_begin_ = shalfedges_last_ = sncp()->shalfedges_end();

      SFace_iterator fit = sfaces_begin(),
                     fend = sfaces_end();
      while (fit != fend) {
        SFace_iterator fdel = fit++;
        sncp()->delete_sface_only(fdel);
      }
      sfaces_begin_ = sfaces_last_ = sncp()->sfaces_end();
    }

    Point_3& point() { return point_at_center_; }
    Mark& mark() { return mark_; }
    GenPtr& info() { return info_; }

  }; // Vertex


//-----------------------------------------------------------------------------

  template <typename Refs>
  class Halfedge : public CGAL::In_place_list_base< Halfedge<Refs> >
  { // == SVertex
    typedef typename Refs::Items Items;
    friend class SNC_structure<Items>;
    friend class SNC_decorator<Refs>;
    friend class SNC_io_parser<Refs>;
    friend class SNC_SM_decorator<Refs>;
    friend class SNC_SM_const_decorator<Refs>;
    typedef typename Refs::Vertex_handle    Vertex_handle;
    typedef typename Refs::SVertex_handle   SVertex_handle;
    typedef typename Refs::SHalfedge_handle SHalfedge_handle;
    typedef typename Refs::SFace_handle     SFace_handle;

    Vertex_handle      center_vertex_;
    Mark               mark_;
    SVertex_handle     twin_;
    SHalfedge_handle   out_sedge_; 
    SFace_handle       incident_sface_;
    GenPtr             info_;
    // temporary information:
    Sphere_point       point_on_surface_;

  public:

    Halfedge() : center_vertex_(), mark_(), twin_(),
      out_sedge_(), incident_sface_(),
      info_(), point_on_surface_() {}

    Halfedge(Mark m) :  center_vertex_(), mark_(m), twin_(),
      out_sedge_(), incident_sface_(),
      info_(), point_on_surface_() {}

    ~Halfedge() {}

    Halfedge(const Halfedge<Refs>& e) 
    { center_vertex_ = e.center_vertex_;
      point_on_surface_ = e.point_on_surface_;
      mark_ = e.mark_;
      twin_ = e.twin_;
      out_sedge_ = e.out_sedge_;
      incident_sface_ = e.incident_sface_;
      info_ = 0;
    }

    Halfedge<Refs>& operator=(const Halfedge<Refs>& e) 
    { center_vertex_ = e.center_vertex_;
      point_on_surface_ = e.point_on_surface_;
      mark_ = e.mark_;
      twin_ = e.twin_;
      out_sedge_ = e.out_sedge_;
      incident_sface_ = e.incident_sface_;
      info_ = 0;
      return *this;
    }


    Mark& mark() 
    { return mark_; }

    Sphere_point& tmp_point()
    { return point_on_surface_; }
    const Sphere_point& tmp_point() const
    { return point_on_surface_; }

  public:
    std::string debug() const
    { std::ostrstream os; set_pretty_mode(os);
      os<<"v"<<tmp_point()<<info_<<'\0';
      std::string res(os.str()); os.freeze(0); return res;
    }

    bool is_twin() const { return (&*twin_ < this); }

  }; // Halfedge


//-----------------------------------------------------------------------------

  template <typename Refs>
  class Halffacet : public CGAL::In_place_list_base< Halffacet<Refs> >
  {
    typedef typename Refs::Items Items;
    friend class SNC_structure<Items>;
    friend class SNC_decorator<Refs>;
    friend class SNC_io_parser<Refs>;
    friend class SNC_SM_decorator<Refs>;
    friend class SNC_SM_const_decorator<Refs>;
    friend class SNC_FM_decorator<Refs>;
    typedef typename Refs::Halffacet_handle   Halffacet_handle;
    typedef typename Refs::Volume_handle  Volume_handle;
    typedef typename Refs::Object_list    Object_list;
    typedef typename Refs::Halffacet_cycle_iterator
                                          Halffacet_cycle_iterator;
    typedef typename Refs::Halffacet_cycle_const_iterator
                                          Halffacet_cycle_const_iterator;

    Plane_3              supporting_plane_;
    Mark                 mark_;
    Halffacet_handle         twin_;
    Volume_handle        volume_;
    Object_list          boundary_entry_objects_; // SEdges, SLoops

  public:

    Halffacet() : supporting_plane_(), mark_() {}

    Halffacet(const Plane_3& h, Mark m) :
      supporting_plane_(h), mark_(m) {}

    ~Halffacet() {}

    Halffacet(const Halffacet<Refs>& f)
    { supporting_plane_ = f.supporting_plane_;
      mark_ = f.mark_;
      twin_ = f.twin_;
      volume_ = f.volume_;
      boundary_entry_objects_ = f.boundary_entry_objects_;
    }

    Halffacet<Refs>& operator=(const Halffacet<Refs>& f)
    { if (this == &f) return *this;
      supporting_plane_ = f.supporting_plane_;
      mark_ = f.mark_;
      twin_ = f.twin_;
      volume_ = f.volume_;
      boundary_entry_objects_ = f.boundary_entry_objects_;
      return *this;
    }

    Mark& mark() 
    { if ( this < &*twin_ ) return mark_; 
      else return twin_->mark_; }

    Plane_3& plane()
    { return supporting_plane_; }

    Halffacet_cycle_iterator facet_cycles_begin()
    { return boundary_entry_objects_.begin(); }
    Halffacet_cycle_iterator facet_cycles_end()
    { return boundary_entry_objects_.end(); }
    Halffacet_cycle_const_iterator facet_cycles_begin() const
    { return boundary_entry_objects_.begin(); }
    Halffacet_cycle_const_iterator facet_cycles_end() const
    { return boundary_entry_objects_.end(); }

    bool is_twin() const { return (&*twin_ < this); }


  }; // Halffacet


//-----------------------------------------------------------------------------

  template <typename Refs>
  class Volume : public CGAL::In_place_list_base< Volume<Refs> >
  {
    typedef typename Refs::Items Items;
    friend class SNC_structure<Items>;
    friend class SNC_decorator<Refs>;
    friend class SNC_io_parser<Refs>;
    friend class SNC_SM_decorator<Refs>;
    friend class SNC_SM_const_decorator<Refs>;
    typedef typename Refs::Object_handle  Object_handle;
    typedef typename Refs::Volume_handle  Volume_handle;
    typedef typename Refs::SObject_list   SObject_list;
    typedef typename Refs::Shell_entry_iterator
                                          Shell_entry_iterator;
    typedef typename Refs::Shell_entry_const_iterator
                                          Shell_entry_const_iterator;

    Mark         mark_;
    SObject_list shell_entry_objects_; // SFaces

  public:

    Volume() {}

    Volume(Mark m) : mark_(m) {}

    ~Volume() {}

    Volume(const Volume<Refs>& v)
    { mark_ = v.mark_;
      shell_entry_objects_ = v.shell_entry_objects_;
    }

    Volume<Refs>& operator=(const Volume<Refs>& v)
    { if (this == &v) return *this;
      mark_ = v.mark_;
      shell_entry_objects_ = v.shell_entry_objects_;
      return *this;
    }

    Mark& mark() 
    { return mark_; }

    Shell_entry_iterator shells_begin()
    { return shell_entry_objects_.begin(); }
    Shell_entry_iterator shells_end()
    { return shell_entry_objects_.end(); }
    Shell_entry_const_iterator shells_begin() const
    { return shell_entry_objects_.begin(); }
    Shell_entry_const_iterator shells_end() const
    { return shell_entry_objects_.end(); }



  }; // Volume



//-----------------------------------------------------------------------------

  template <typename Refs>
  class SHalfedge : public CGAL::In_place_list_base< SHalfedge<Refs> >
  { 
    typedef typename Refs::Items Items;
    typedef typename Refs::SVertex_handle SVertex_handle;
    typedef typename Refs::SHalfedge_handle SHalfedge_handle;
    typedef typename Refs::SHalfedge_const_handle SHalfedge_const_handle;
    typedef typename Refs::SFace_handle SFace_handle;
    typedef typename Refs::Halffacet_handle Halffacet_handle;
    friend class SNC_structure<Items>;
    friend class SNC_decorator<Refs>;
    friend class SNC_io_parser<Refs>;
    friend class SNC_SM_decorator<Refs>;
    friend class SNC_SM_const_decorator<Refs>;
    friend class move_shalfedge_around_svertex<SHalfedge_handle>;
    friend class move_shalfedge_around_sface<SHalfedge_handle>;
    friend class move_shalfedge_around_facet<SHalfedge_handle>;
    friend class move_shalfedge_around_svertex<SHalfedge_const_handle>;
    friend class move_shalfedge_around_sface<SHalfedge_const_handle>;
    friend class move_shalfedge_around_facet<SHalfedge_const_handle>;

    // Role within local graph:
    SVertex_handle     source_;
    SHalfedge_handle   sprev_, snext_;
    SFace_handle       incident_sface_;
    SHalfedge_handle   twin_;
    // Topology within global Nef structure:
    SHalfedge_handle   prev_, next_; 
    Halffacet_handle       incident_facet_;
    GenPtr             info_;
    // temporary needed:
    Mark               mark_;
    Sphere_circle      circle_;

  public:

    SHalfedge() : source_(), sprev_(), snext_(),
     incident_sface_(), twin_(), 
     prev_(), next_(), incident_facet_(), 
     info_(), mark_(), circle_() {}

    ~SHalfedge() {}

    SHalfedge(const SHalfedge<Refs>& e)
    {
      source_ = e.source_;
      sprev_ = e.sprev_;
      snext_ = e.snext_;
      incident_sface_ = e.incident_sface_;
      twin_ = e.twin_;
      prev_ = e.prev_;
      next_ = e.next_;
      incident_facet_ = e.incident_facet_;
      info_ = 0;
      mark_ = e.mark_;
      circle_ = e.circle_;
    }

    SHalfedge<Refs>& operator=(const SHalfedge<Refs>& e)
    {
      source_ = e.source_;
      sprev_ = e.sprev_;
      snext_ = e.snext_;
      incident_sface_ = e.incident_sface_;
      twin_ = e.twin_;
      prev_ = e.prev_;
      next_ = e.next_;
      incident_facet_ = e.incident_facet_;
      info_ = 0;
      mark_ = e.mark_;
      circle_ = e.circle_;
      return *this;
    }

    Mark& mark() 
    { return incident_facet_->mark(); }

    // Plane_3& plane() { return incident_facet_->plane(); }

    Mark& tmp_mark() 
    { if ( this < &*twin_ ) return mark_;
      else return twin_->mark_; }

    Sphere_circle& tmp_circle() { return circle_; }
    const Sphere_circle& tmp_circle() const { return circle_; }

  public:
    std::string debug() const
    { std::ostrstream os; set_pretty_mode(os); 
      os <<"e["<<source_->debug()<<","
         <<twin_->source_->debug()<<" "<<info_<<"]"<<'\0';
      std::string res(os.str()); os.freeze(0); return res;
    }

  }; // SHalfedge


//-----------------------------------------------------------------------------


  template <typename Refs> 
  class SHalfloop : public CGAL::In_place_list_base< SHalfloop<Refs> >
  {
    typedef typename Refs::Vertex_handle Vertex_handle;
    typedef typename Refs::SHalfloop_handle SHalfloop_handle;
    typedef typename Refs::SFace_handle SFace_handle;
    typedef typename Refs::Halffacet_handle Halffacet_handle;
    typedef typename Refs::Items Items;
    friend class SNC_structure<Items>;
    friend class SNC_decorator<Refs>;
    friend class SNC_io_parser<Refs>;
    friend class SNC_SM_decorator<Refs>;
    friend class SNC_SM_const_decorator<Refs>;
    friend class Self::Vertex<Refs>;

    SHalfloop_handle   twin_;
    SFace_handle       incident_sface_;
    Halffacet_handle       incident_facet_;
    GenPtr             info_;
    // temporary needed:
    Mark               mark_;
    Sphere_circle      circle_;

  public:

    SHalfloop() : twin_(), incident_sface_(), incident_facet_(), 
      info_(), mark_(), circle_() {}

    ~SHalfloop() {}

    SHalfloop(const SHalfloop<Refs>& l)
    { twin_ = l.twin_;
      incident_sface_ = l.incident_sface_;
      incident_facet_ = l.incident_facet_;
      info_ = 0;
      mark_ = l.mark_;
      circle_ = l.circle_;
    }

    SHalfloop<Refs>& operator=(const SHalfloop<Refs>& l)
    { twin_ = l.twin_;
      incident_sface_ = l.incident_sface_;
      incident_facet_ = l.incident_facet_;
      info_ = 0;
      mark_ = l.mark_;
      circle_ = l.circle_;
      return *this;
    }

    Mark& mark() 
    { return incident_facet_->mark(); }

    Mark& tmp_mark() 
    { if ( this < &*twin_ ) return mark_;
      else return twin_->mark_; }

    Sphere_circle& tmp_circle() { return circle_; }
    const Sphere_circle& tmp_circle() const { return circle_; }

  public:
    std::string debug() const
    { std::ostrstream os; set_pretty_mode(os); 
      os<<"l"<<tmp_circle()<<'\0';
      std::string res(os.str()); os.freeze(0); return res;
    }

  }; // SHalfloop


//-----------------------------------------------------------------------------

  template <typename Refs>
  class SFace : public CGAL::In_place_list_base< SFace<Refs> >
  { 
    typedef typename Refs::Items Items;
    friend class SNC_structure<Items>;
    friend class SNC_decorator<Refs>;
    friend class SNC_io_parser<Refs>;
    friend class SNC_SM_decorator<Refs>;
    friend class SNC_SM_const_decorator<Refs>;
    typedef typename Refs::Vertex_handle  Vertex_handle;
    typedef typename Refs::SFace_handle   SFace_handle;
    typedef typename Refs::SObject_handle SObject_handle;
    typedef typename Refs::Volume_handle  Volume_handle;
    typedef typename Refs::SObject_list   SObject_list;
    typedef typename Refs::SFace_cycle_iterator 
                                          SFace_cycle_iterator;
    typedef typename Refs::SFace_cycle_const_iterator 
                                          SFace_cycle_const_iterator;
    Vertex_handle  center_vertex_;
    Volume_handle  incident_volume_;
    SObject_list   boundary_entry_objects_; // SEdges, SLoops, SVertices
    GenPtr         info_;
    // temporary needed:
    Mark           mark_;

  public:

    SFace() : center_vertex_(), incident_volume_(), info_(), mark_() {}

    ~SFace() {}

    SFace(const SFace<Refs>& f)
    { center_vertex_ = f.center_vertex_;
      incident_volume_ = f.incident_volume_;
      boundary_entry_objects_ = f.boundary_entry_objects_;
      info_ = 0;
      mark_ = f.mark_;
    }

    SFace<Refs>& operator=(const SFace<Refs>& f)
    { if (this == &f) return *this;
      center_vertex_ = f.center_vertex_;
      incident_volume_ = f.incident_volume_;
      boundary_entry_objects_ = f.boundary_entry_objects_;
      info_ = 0;
      mark_ = f.mark_;
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

    Mark& mark() 
    { return incident_volume_->mark(); }

    Mark& tmp_mark() 
    { return mark_; }

  }; // SFace

}; // SNC_items



CGAL_END_NAMESPACE
#endif //CGAL_SNC_ITEMS_H
