// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Michael Seel        <seel@mpi-sb.mpg.de>
//                 Miguel Granados     <granados@mpi-sb.mpg.de>
//                 Susan Hert          <hert@mpi-sb.mpg.de>
//                 Lutz Kettner        <kettner@mpi-sb.mpg.de>
//                 Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_SNC_ITEMS_H
#define CGAL_SNC_ITEMS_H

#include <CGAL/basic.h>
#include <string>
#include <sstream>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/Nef_3/SNC_iteration.h>
#include <CGAL/Nef_3/SNC_list.h>

#undef _DEBUG
#define _DEBUG 83
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

template <typename K, typename I> class SNC_sphere_map;
template <typename R> class SM_decorator;

class SNC_items {
 public:

//-----------------------------------------------------------------------------

  template <typename Refs>
  class Vertex
  {
    typedef void* GenPtr;
    typedef typename Refs::Mark  Mark;
    typedef typename Refs::Point_3 Point_3;
    typedef typename Refs::Sphere_map Sphere_map;
    typedef typename Refs::SVertex_iterator SVertex_iterator;
    typedef typename Refs::SHalfedge_iterator SHalfedge_iterator;
    typedef typename Refs::SHalfloop_iterator SHalfloop_iterator;
    typedef typename Refs::SFace_iterator SFace_iterator;

    typedef typename Refs::Vertex_handle Vertex_handle;
    typedef typename Refs::Vertex_iterator Vertex_iterator;
    typedef typename Refs::SHalfloop_handle SHalfloop_handle;
    typedef typename Refs::SHalfloop SHalfloop;
    typedef typename Refs::SVertex_const_iterator SVertex_const_iterator;
    typedef typename Refs::SHalfedge_const_iterator SHalfedge_const_iterator;
    typedef typename Refs::SFace_const_iterator SFace_const_iterator;
    typedef typename Refs::SHalfloop_const_handle SHalfloop_const_handle;
    typedef typename Refs::Vertex_const_iterator Vertex_const_iterator;

    typedef typename Refs::Size_type  Size_type;

    Point_3            point_at_center_;
    Mark               mark_;
    // local view (surface graph):
  public:
    Refs*              sncp_;
    SVertex_iterator   svertices_begin_, svertices_last_;
    SHalfedge_iterator shalfedges_begin_, shalfedges_last_;
    SFace_iterator     sfaces_begin_, sfaces_last_;
    SHalfloop_iterator shalfloop_;
    GenPtr             info_;

  public:

    Vertex() : point_at_center_(), mark_(), sncp_(), 
      svertices_begin_(), svertices_last_(),
      shalfedges_begin_(), shalfedges_last_(),
      sfaces_begin_(), sfaces_last_(), shalfloop_(),
      info_() 
      // , sm_(Vertex_handle((SNC_in_place_list_vertex<Vertex>*) this)) 
      {}

    Vertex(const Point_3& p, Mark m) : 
      point_at_center_(p), mark_(m), sncp_(), 
      svertices_begin_(), svertices_last_(),
      shalfedges_begin_(), shalfedges_last_(),
      sfaces_begin_(), sfaces_last_(), shalfloop_(),
      info_() 
      //    , sm_(Vertex_handle((SNC_in_place_list_vertex<Vertex>*) this)) 
      {}

    Vertex(const Vertex<Refs>& v) 
      //      : sm_(Vertex_handle((SNC_in_place_list_vertex<Vertex>*)this)) 
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
      return *this;
    }

    Refs* sncp() const { return sncp_; }
    Refs*& sncp() { return sncp_; }

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

    SVertex_iterator& svertices_begin() 
    { return svertices_begin_; }
    SVertex_iterator& svertices_last() 
    { return svertices_last_; }
    SVertex_iterator svertices_end() 
    { if ( svertices_last_ == sncp()->svertices_end() ) 
        return svertices_last_; 
      else 
	return ++SVertex_iterator(svertices_last_); }

    SHalfedge_iterator& shalfedges_begin() 
    { return shalfedges_begin_; }
    SHalfedge_iterator& shalfedges_last() 
    { return shalfedges_last_; }
    SHalfedge_iterator shalfedges_end() 
    { if ( shalfedges_last_ == sncp()->shalfedges_end() ) 
        return shalfedges_last_;
      else 
	return ++SHalfedge_iterator(shalfedges_last_); }

    SFace_iterator& sfaces_begin() 
    { return sfaces_begin_; }
    SFace_iterator& sfaces_last()
    { return sfaces_last_; }
    SFace_iterator sfaces_end()
    { if ( sfaces_last_ == sncp()->sfaces_end() ) 
        return sfaces_last_;
      else 
	return ++SFace_iterator(sfaces_last_); }

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

    SHalfloop_handle& shalfloop() { return shalfloop_; }
    SHalfloop_handle shalfloop() const { return shalfloop_; }

    bool has_shalfloop() const {
      return shalfloop_ != sncp()->shalfloops_end();
    }

    Size_type number_of_svertices() const 
      /*{\Mop returns the number of vertices.}*/
      { Size_type n(0);
      SVertex_const_iterator vit;
      CGAL_forall_svertices(vit, *this) ++n;
      return n; }
    
    Size_type number_of_shalfedges() const 
      /*{\Mop returns the number of halfedges.}*/
      { Size_type n(0);
      SHalfedge_const_iterator eit;
      CGAL_forall_shalfedges(eit, *this) ++n;
      return n;}
    
    Size_type number_of_sedges() const 
      /*{\Mop returns the number of edges.}*/
      { return number_of_shalfedges()/2; }
    
    Size_type number_of_shalfloops() const 
      /*{\Mop returns the number of halfloops.}*/
      { return ( has_shalfloop() ? 2 : 0); }
    
    Size_type number_of_sloops() const 
      /*{\Mop returns the number of loops.}*/
      { return number_of_shalfloops()/2; }
    
    Size_type number_of_sfaces() const    
      /*{\Mop returns the number of faces.}*/
      { Size_type n(0);
      SFace_const_iterator fit;
      CGAL_forall_sfaces(fit, *this) ++n;
      return n; }


    /*{\Xtext Vertices provide access to their local graphs via
    the iterator ranges:
    \begin{Mverb}
    SVertex_iterator     svertices_begin()/svertices_end()
    SHalfedge_iterator   shalfedges_begin()/shalfedges_end()
    SFace_iterator       sfaces_begin()/sfaces_end() 
    SHalfloop_handle     shalfloop()
    \end{Mverb}
    }*/

    void clear() 
    /*{\Xop clears the local graph.}*/ { 
      SFace_iterator fit = sfaces_begin(),
                     fend = sfaces_end();
      while (fit != fend) {
        SFace_iterator fdel = fit++;
	/* TO VERIFY: next statement needs access to a private attribute */
	sncp()->reset_sm_object_list(fdel->boundary_entry_objects_);
        sncp()->delete_sface_only(fdel);
      }
      sfaces_begin_ = sfaces_last_ = sncp()->sfaces_end();

      if ( shalfloop() != sncp()->shalfloops_end() ) {
        sncp()->delete_shalfloop_only(shalfloop_->twin());
        sncp()->delete_shalfloop_only(shalfloop_);
        shalfloop_ = sncp()->shalfloops_end();
      }

      SHalfedge_iterator eit = shalfedges_begin(),
                         eend = shalfedges_end();
      while (eit != eend) {
        SHalfedge_iterator edel = eit++;
        sncp()->delete_shalfedge_only(edel);
      }
      shalfedges_begin_ = shalfedges_last_ = sncp()->shalfedges_end();

      SVertex_iterator vit = svertices_begin(),
                       vend = svertices_end();
      while (vit != vend) {
        SVertex_iterator vdel = vit++;
        sncp()->delete_halfedge_only(vdel);
      }
      svertices_begin_ = svertices_last_ = sncp()->halfedges_end();
    }

    Point_3& point() { return point_at_center_; }
    const Point_3& point() const { return point_at_center_; }
    Mark& mark() { return mark_; }
    const Mark& mark() const { return mark_;}
    GenPtr& info() { return info_; }
    const GenPtr& info() const { return info_; }

    ~Vertex() {
      TRACEN("  destroying Vertex item "<<&*this);
    }

  public:
    std::string debug() const
    { std::stringstream os; 
      set_pretty_mode(os); 
      os<<"{ addr, point, mark, snc, svb, sve, seb, see, sfb, sfe, sl,"
	<<" info }"<<std::endl;
      os<<"{ "<<this<<", "<<point_at_center_<<", "<<mark_<<", "<<&*sncp_<<", "
	<<&*svertices_begin_ <<", "<<&*svertices_last_ <<", "
	<<&*shalfedges_begin_<<", "<<&*shalfedges_last_<<", "
	<<&*sfaces_begin_    <<", "<<&*sfaces_last_    <<", "
	<<&*shalfloop_       <<", "
	<<info_<<" }";
      return os.str();
    }

    bool check_basic_functions() {
      /*
      if(svertices_begin_ == sncp()->svertices_end())
        CGAL_assertion(svertices_end() == sncp()->svertices_end());
      if(shalfedges_begin_ == sncp()->shalfedges_end())
        CGAL_assertion(shalfedges_end() == sncp()->shalfedges_end());
      if(sfaces_begin_ == sncp()->sfaces_end())
        CGAL_assertion(sfaces_end() == sncp()->sfaces_end());
      */
    }

    bool is_valid( bool verb = false, int level = 0) const {
      
      Verbose_ostream verr(verb);
      verr << "begin CGAL::SNC_items<...>::Vertex::is_valid( verb=true, "
	"level = " << level << "):" << std::endl;
      
      bool valid = (sncp_ != NULL);
      valid = valid && (svertices_begin_ != NULL && svertices_begin_ != SVertex_iterator());
      valid = valid && (svertices_last_  != NULL && svertices_last_  != SVertex_iterator());
      valid = valid && (shalfedges_begin_ != NULL && shalfedges_begin_ != SHalfedge_iterator());
      valid = valid && (shalfedges_last_  != NULL && shalfedges_last_  != SHalfedge_iterator());
      valid = valid && (sfaces_begin_ != NULL && sfaces_begin_ != SFace_iterator());
      valid = valid && (sfaces_last_  != NULL && sfaces_last_  != SFace_iterator());
      valid = valid && (shalfloop_ != NULL && shalfloop_ != SHalfloop_iterator());

      if(shalfedges_begin_ == sncp()->shalfedges_end()) {         // point in volume or on plane, which is either isolated or has one outgoing edge
	if(shalfloop_ != sncp()->shalfloops_end())
	  valid = valid && (++SFace_const_iterator(sfaces_begin_) == sfaces_last_);
	else 
	  valid = valid && (sfaces_begin_ == sfaces_last_);
      }

      valid = valid && (sfaces_begin_ != sncp()->sfaces_end());
      if(sfaces_begin_ == sfaces_last_) {
	valid = valid && (shalfloop_ == sncp()->shalfloops_end());
      }
      else
	valid = valid && (sfaces_begin_->sface_cycles_begin() != 
			  sfaces_begin_->sface_cycles_end());
     
      verr << "end of CGAL::SNC_items<...>::Vertex::is_valid(): structure is "
	   << ( valid ? "valid." : "NOT VALID.") << std::endl;
      return valid;
    }
    
  }; // Vertex


//-----------------------------------------------------------------------------

  template <typename Refs>
  class SVertex
  { // == Halfedge
    typedef void* GenPtr;
    typedef typename Refs::Mark  Mark;
    typedef typename Refs::Sphere_point  Sphere_point;
    typedef typename Refs::Sphere_map Sphere_map;
    friend class SM_decorator<Refs>;
    typedef typename Refs::Vertex_handle    Vertex_handle;
    typedef typename Refs::Halfedge_handle  Halfedge_handle;
    typedef typename Refs::SVertex_handle   SVertex_handle;
    typedef typename Refs::SHalfedge_handle SHalfedge_handle;
    typedef typename Refs::SFace_handle     SFace_handle;

    Vertex_handle      center_vertex_;
    Mark               mark_;
    SVertex_handle     twin_;              
    SHalfedge_handle   out_sedge_;           
    SFace_handle       incident_sface_;    
    GenPtr             info_;                 
    Sphere_point       point_;   

  public:

    SVertex() : center_vertex_(), mark_(), twin_(),
      out_sedge_(), incident_sface_(),
      info_(), point_() {}

    SVertex(Mark m) :  center_vertex_(), mark_(m), twin_(),
      out_sedge_(), incident_sface_(),
      info_(), point_() {}

    ~SVertex() {
      TRACEN("  destroying Halfedge item "<<&*this);
    }

    SVertex(const SVertex<Refs>& e) 
    { center_vertex_ = e.center_vertex_;
      point_ = e.point_;
      mark_ = e.mark_;
      twin_ = e.twin_;
      out_sedge_ = e.out_sedge_;
      incident_sface_ = e.incident_sface_;
      info_ = 0;
    }

    SVertex<Refs>& operator=(const SVertex<Refs>& e) 
    { center_vertex_ = e.center_vertex_;
      point_ = e.point_;
      mark_ = e.mark_;
      twin_ = e.twin_;
      out_sedge_ = e.out_sedge_;
      incident_sface_ = e.incident_sface_;
      info_ = 0;
      return *this;
    }

    Vertex_handle& center_vertex() { return center_vertex_; }
    Vertex_handle center_vertex() const { return center_vertex_; }

    Mark& mark() { return mark_; }
    const Mark& mark() const { return mark_; }

    Sphere_point& vector(){ return point_; }
    const Sphere_point& vector() const { return point_; }
    Sphere_point& point(){ return point_; }
    const Sphere_point& point() const { return point_; }

    SVertex_handle& twin() { return twin_; }
    SVertex_handle twin()  const { return twin_; }

    SHalfedge_handle& out_sedge() { return out_sedge_; }
    SHalfedge_handle out_sedge() const { return out_sedge_; }

    SFace_handle& incident_sface() { return incident_sface_; } 
    SFace_handle incident_sface() const { return incident_sface_; } 

    GenPtr& info() { return info_; }
    const GenPtr& info() const { return info_; }

  public:
    std::string debug() const
    { std::stringstream os; 
      set_pretty_mode(os);
      os<<"sv [ "<<point_<<info_<<" ] ";
      return os.str();
    }

    bool is_twin() const { return (&*twin_ < this); }

    bool is_valid( bool verb = false, int level = 0) const {
      
      Verbose_ostream verr(verb);
      verr << "begin CGAL::SNC_items<...>::SVertex::is_valid( verb=true, "
	"level = " << level << "):" << std::endl;

      bool valid = (center_vertex_ != NULL && center_vertex_ != Vertex_handle());
      valid = valid && (twin_ != NULL && twin_ != SVertex_handle() &&
			                 twin_ != SVertex_handle());
      //      valid = valid && (out_sedge_ != NULL);
      //      valid = valid && (incident_sface_ != SFace_handle());
      
      valid = valid &&((out_sedge_ != NULL && incident_sface_ == NULL) ||
		       (out_sedge_ == NULL && incident_sface_ != NULL));
      
      verr << "end of CGAL::SNC_items<...>::SVertex::is_valid(): structure is "
	   << ( valid ? "valid." : "NOT VALID.") << std::endl;

      return valid;
    }

  }; // SVertex


//-----------------------------------------------------------------------------

  template <typename Refs>
  class Halffacet  {

    typedef void* GenPtr;
    typedef typename Refs::Mark  Mark;
    typedef typename Refs::Plane_3   Plane_3;
    typedef typename Refs::Halffacet_handle   Halffacet_handle;
    typedef typename Refs::Volume_handle  Volume_handle;
    typedef typename Refs::SHalfedge_handle SHalfedge_handle;
    typedef typename Refs::SHalfloop_handle SHalfloop_handle;
    typedef typename Refs::Object_list    Object_list;
    typedef typename Refs::Object_const_iterator Object_iterator;
    typedef typename Refs::Halffacet_cycle_iterator
                                          Halffacet_cycle_iterator;
    typedef typename Refs::Halffacet_cycle_const_iterator
                                          Halffacet_cycle_const_iterator;
  public: // TODO: TO REMOVE!!!


    Plane_3              supporting_plane_;
    Mark                 mark_;
    Halffacet_handle     twin_;
    Volume_handle        volume_;
    Object_list          boundary_entry_objects_; // SEdges, SLoops

  public:

    Halffacet() : supporting_plane_(), mark_() {}

    Halffacet(const Plane_3& h, Mark m) :
      supporting_plane_(h), mark_(m) {}

    ~Halffacet() {
      TRACEN("  destroying Halffacet item "<<&*this);
    }

    Halffacet(const Halffacet<Refs>& f)
    { supporting_plane_ = f.supporting_plane_;
      mark_ = f.mark_;
      twin_ = f.twin_;
      TRACEN("VOLUME const");
      volume_ = f.volume_;
      boundary_entry_objects_ = f.boundary_entry_objects_;
    }

    Halffacet<Refs>& operator=(const Halffacet<Refs>& f)
    { if (this == &f) return *this;
      supporting_plane_ = f.supporting_plane_;
      mark_ = f.mark_;
      twin_ = f.twin_;
      TRACEN("VOLUME op=");
      volume_ = f.volume_;
      boundary_entry_objects_ = f.boundary_entry_objects_;
      return *this;
    }

    Mark& mark() { return mark_; }
    const Mark& mark() const { return mark_; }

    Halffacet_handle& twin() { return twin_; }
    Halffacet_handle twin() const { return twin_; }

    Plane_3& plane() { return supporting_plane_; }
    const Plane_3& plane() const { return supporting_plane_; }

    Volume_handle& volume() { return volume_; }
    Volume_handle volume() const { return volume_; }

    Object_list& boundary_entry_objects() { return boundary_entry_objects_; }
    const Object_list& boundary_entry_objects() const { return boundary_entry_objects_; }

    GenPtr& info() { return info_; }
    const GenPtr& info() const { return info_; }

    Halffacet_cycle_iterator facet_cycles_begin()
    { return boundary_entry_objects_.begin(); }
    Halffacet_cycle_iterator facet_cycles_end()
    { return boundary_entry_objects_.end(); }
    Halffacet_cycle_const_iterator facet_cycles_begin() const
    { return boundary_entry_objects_.begin(); }
    Halffacet_cycle_const_iterator facet_cycles_end() const
    { return boundary_entry_objects_.end(); }

    bool is_twin() const { return (&*twin_ < this); }

    bool is_valid( bool verb = false, int level = 0) const {
      
      Verbose_ostream verr(verb);
      verr << "begin CGAL::SNC_items<...>::Halffacet::is_valid( verb=true, "
	"level = " << level << "):" << std::endl;

      bool valid = (twin_ != NULL && twin_ != Halffacet_handle());
      valid = valid && (volume_ != NULL && volume_ != Volume_handle());
      
      valid = valid && (supporting_plane_.a() != 0 || 
			supporting_plane_.b() != 0 ||
			supporting_plane_.c() != 0);
      
      valid = valid && (!boundary_entry_objects_.empty());

      verr << "end of CGAL::SNC_items<...>::Halffacet::is_valid(): structure is "
	   << ( valid ? "valid." : "NOT VALID.") << std::endl;

      return valid;
    }

  }; // Halffacet


//-----------------------------------------------------------------------------

  template <typename Refs>
  class Volume  {
    typedef void* GenPtr;
    typedef typename Refs::Mark  Mark;
    typedef typename Refs::Object_handle  Object_handle;
    typedef typename Refs::Volume_handle  Volume_handle;
    typedef typename Refs::Object_list   Object_list;
    typedef typename Refs::Shell_entry_iterator
                                          Shell_entry_iterator;
    typedef typename Refs::Shell_entry_const_iterator
                                          Shell_entry_const_iterator;

    Mark         mark_;
    Object_list shell_entry_objects_; // SFaces

  public:

    Volume() {}

    Volume(Mark m) : mark_(m) {}

    ~Volume() {
      TRACEN("  destroying Volume item "<<&*this);
    }

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

    Mark& mark() { return mark_; }
    const Mark& mark() const { return mark_; }

    GenPtr& info() { return info_; }
    const GenPtr& info() const { return info_; }

    Object_list& shell_entry_objects() { return shell_entry_objects_; }
    const Object_list& shell_entry_objects() const { 
      return shell_entry_objects_; 
    }

    Shell_entry_iterator shells_begin()
    { return shell_entry_objects_.begin(); }
    Shell_entry_iterator shells_end()
    { return shell_entry_objects_.end(); }
    Shell_entry_const_iterator shells_begin() const
    { return shell_entry_objects_.begin(); }
    Shell_entry_const_iterator shells_end() const
    { return shell_entry_objects_.end(); }

    bool is_valid( bool verb = false, int level = 0) const {
      
      Verbose_ostream verr(verb);
      verr << "begin CGAL::SNC_items<...>::Volume::is_valid( verb=true, "
	"level = " << level << "):" << std::endl;

      bool valid = (!shell_entry_objects_.empty());

      verr << "end of CGAL::SNC_items<...>::Volume::is_valid(): structure is "
	   << ( valid ? "valid." : "NOT VALID.") << std::endl;

      return valid;
    }

  }; // Volume



//-----------------------------------------------------------------------------

  template <typename Refs>
  class SHalfedge  { 
    typedef void* GenPtr;
    typedef typename Refs::Mark  Mark;
    typedef typename Refs::Sphere_circle  Sphere_circle;
    typedef typename Refs::Sphere_map Sphere_map;
    typedef typename Refs::Halfedge_handle Halfedge_handle;
    typedef typename Refs::SVertex_handle SVertex_handle;
    typedef typename Refs::SHalfedge_handle SHalfedge_handle;
    typedef typename Refs::SHalfedge_const_handle SHalfedge_const_handle;
    typedef typename Refs::SFace_handle SFace_handle;
    typedef typename Refs::Halffacet_handle Halffacet_handle;
    friend class SM_decorator<Refs>;

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

    ~SHalfedge() {
      TRACEN("  destroying SHalfedge item "<<&*this);
    }

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

    Mark& mark() { return mark_; }
    const Mark& mark() const { return mark_; }

    SHalfedge_handle& twin() { return twin_; }
    SHalfedge_handle twin() const { return twin_; }

    SVertex_handle& source() { return source_; }
    SVertex_handle source() const { return source_; }

    SHalfedge_handle& prev() { return prev_; }
    SHalfedge_handle prev() const { return prev_; }

    SHalfedge_handle& next() { return next_; }
    SHalfedge_handle next() const { return next_; }

    SHalfedge_handle& sprev() { return sprev_; }
    SHalfedge_handle sprev() const { return sprev_; }

    SHalfedge_handle& snext() { return snext_; }
    SHalfedge_handle snext() const { return snext_; }

    Sphere_circle& circle() { return circle_; }
    const Sphere_circle& circle() const { return circle_; }
    
    SFace_handle& incident_sface() { return incident_sface_; }
    SFace_handle incident_sface() const { return incident_sface_; }

    Halffacet_handle& incident_facet() { return incident_facet_; }
    Halffacet_handle incident_facet() const { return incident_facet_; }

    GenPtr& info() { return info_; }
    const GenPtr& info() const { return info_; }

  public:
    std::string debug() const
    { std::stringstream os; 
      set_pretty_mode(os); 
      os <<"e[ "<<source_->debug()<<", "
         <<twin_->source_->debug()<<" "<<info_<<" ] ";
      return os.str();
    }

    bool is_twin() const { return (&*twin_ < this); }

    bool is_valid( bool verb = false, int level = 0) const {
      
      Verbose_ostream verr(verb);
      verr << "begin CGAL::SNC_items<...>::SHalfedge::is_valid( verb=true, "
	"level = " << level << "):" << std::endl;

      bool valid = (source_ != SVertex_handle() &&
		    source_ != NULL &&
		    source_ != Halfedge_handle());
      valid = valid && (twin_  != SHalfedge_handle() && twin_  != NULL);
      valid = valid && (sprev_ != SHalfedge_handle() && sprev_ != NULL);
      valid = valid && (snext_ != SHalfedge_handle() && snext_ != NULL);
      valid = valid && (prev_  != SHalfedge_handle() && prev_  != NULL);
      valid = valid && (next_  != SHalfedge_handle() && next_  != NULL);
      
      valid = valid && (incident_sface_ != SFace_handle() && 
			incident_sface_ != NULL);
      valid = valid && (incident_facet_ != Halffacet_handle() &&
			incident_facet_ != NULL);
      valid = valid && (circle_.d() == 0);
      valid = valid && (circle_.a() != 0 || circle_.b() != 0 || circle_.c() !=0);
      
      verr << "end of CGAL::SNC_items<...>::SHalfedge::is_valid(): structure is "
	   << ( valid ? "valid." : "NOT VALID.") << std::endl;

      return valid;
    }

  }; // SHalfedge


//-----------------------------------------------------------------------------


  template <typename Refs> 
  class SHalfloop {
    typedef void* GenPtr;
    typedef typename Refs::Mark  Mark;
    typedef typename Refs::Sphere_circle  Sphere_circle;
    typedef typename Refs::Vertex_handle Vertex_handle;
    typedef typename Refs::SHalfloop_handle SHalfloop_handle;
    typedef typename Refs::SFace_handle SFace_handle;
    typedef typename Refs::Halffacet_handle Halffacet_handle;
    typedef typename Refs::Sphere_map Sphere_map;
    friend class SM_decorator<Refs>;

    SHalfloop_handle   twin_;
    SFace_handle       incident_sface_;
    Halffacet_handle   incident_facet_;
    GenPtr             info_;
    // temporary needed:
    Mark               mark_;
    Sphere_circle      circle_;

  public:

    SHalfloop() : twin_(), incident_sface_(), incident_facet_(), 
      info_(), mark_(), circle_() {}

    ~SHalfloop() {
      TRACEN("  destroying SHalfloop item "<<&*this);
    }
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

    Mark& mark() { return mark_;}
    const Mark& mark() const { return mark_; }

    SHalfloop_handle& twin() { return twin_; }
    SHalfloop_handle twin() const { return twin_; }

    Sphere_circle& circle() { return circle_; }
    const Sphere_circle& circle() const { return circle_; }

    SFace_handle& incident_sface() { return incident_sface_; }
    SFace_handle incident_sface() const { return incident_sface_; }

    Halffacet_handle& incident_facet() { return incident_facet_; }
    Halffacet_handle incident_facet() const { return incident_facet_; }

    GenPtr& info() { return info_; }
    const GenPtr& info() const { return info_; }

  public:
    std::string debug() const
    { std::stringstream os; 
      set_pretty_mode(os); 
      os<<"sl [ "<<circle_<<" ] ";
      return os.str();
    }

    bool is_twin() const { return (&*twin_ < this); }

    bool is_valid( bool verb = false, int level = 0) const {
      
      Verbose_ostream verr(verb);
      verr << "begin CGAL::SNC_items<...>::SHalfloop::is_valid( verb=true, "
	"level = " << level << "):" << std::endl;

      bool valid = (twin_  != SHalfloop_handle() && twin_  != NULL);
      valid = valid && (incident_sface_ != SFace_handle() && 
			incident_sface_ != NULL);
      valid = valid && (incident_facet_ != Halffacet_handle() &&
			incident_facet_ != NULL);
      valid = valid && (circle_.d() == 0);
      valid = valid && (circle_.a() != 0 || circle_.b() != 0 || circle_.c() !=0);
      
      verr << "end of CGAL::SNC_items<...>::SHalfloop::is_valid(): structure is "
	   << ( valid ? "valid." : "NOT VALID.") << std::endl;

      return valid;
    }

  }; // SHalfloop


//-----------------------------------------------------------------------------

  template <typename Refs>
  class SFace { 
    typedef void* GenPtr;
    typedef typename Refs::Mark  Mark;
    typedef typename Refs::Sphere_map Sphere_map;
    friend class SM_decorator<Refs>;
    typedef typename Refs::Vertex_handle  Vertex_handle;
    typedef typename Refs::SFace_handle   SFace_handle;
    typedef typename Refs::Object_handle Object_handle;
    typedef typename Refs::Volume_handle  Volume_handle;
    typedef typename Refs::Object_list   Object_list;
    typedef typename Refs::SFace_cycle_iterator 
                                          SFace_cycle_iterator;
    typedef typename Refs::SFace_cycle_const_iterator 
                                          SFace_cycle_const_iterator;
    Vertex_handle  center_vertex_;
    Volume_handle  incident_volume_;
    //    Object_list   boundary_entry_objects_; // SEdges, SLoops, SVertices
    GenPtr         info_;
    // temporary needed:
    Mark           mark_;

  public:
    Object_list   boundary_entry_objects_; // SEdges, SLoops, SVertices

    SFace() : center_vertex_(), incident_volume_(), info_(), mark_() {}

    ~SFace() {
      TRACEN("  destroying SFace item "<<&*this);
    }

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

    Mark& mark() { return mark_; }
    const Mark& mark() const { return mark_; }

    Vertex_handle& center_vertex() { return center_vertex_; }
    Vertex_handle center_vertex() const { return center_vertex_; }

    Volume_handle& incident_volume() { return incident_volume_; }
    Volume_handle incident_volume() const { return incident_volume_; }

    Object_list& boundary_entry_objects() { return boundary_entry_objects_; }
    const Object_list& boundary_entry_objects() const { return boundary_entry_objects_; }

    GenPtr& info() { return info_; }
    const GenPtr& info() const { return info_; }

    bool is_valid( bool verb = false, int level = 0) const {
      
      Verbose_ostream verr(verb);
      verr << "begin CGAL::SNC_items<...>::SFace::is_valid( verb=true, "
	"level = " << level << "):" << std::endl;

      bool valid =(center_vertex_ != Vertex_handle() && center_vertex_ != NULL);
      valid = valid && (incident_volume_ != Volume_handle() &&
			incident_volume_ != NULL);

      if(boundary_entry_objects_.empty()) {
	valid = valid && 
	  (center_vertex_->shalfedges_begin() == center_vertex_->shalfedges_end());
      }
      verr << "end of CGAL::SNC_items<...>::SFace::is_valid(): structure is "
	   << ( valid ? "valid." : "NOT VALID.") << std::endl;

      return valid;
    }

  }; // SFace

}; // SNC_items



CGAL_END_NAMESPACE
#endif //CGAL_SNC_ITEMS_H
