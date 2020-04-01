// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel        <seel@mpi-sb.mpg.de>
//                 Miguel Granados     <granados@mpi-sb.mpg.de>
//                 Susan Hert          <hert@mpi-sb.mpg.de>
//                 Lutz Kettner        <kettner@mpi-sb.mpg.de>
//                 Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_NEF_VERTEX_H
#define CGAL_NEF_VERTEX_H

#include <CGAL/license/Nef_3.h>


#include <string>
#include <sstream>
#include <CGAL/IO/io.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/Nef_3/SNC_iteration.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 83
#include <CGAL/Nef_2/debug.h>

#ifndef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <boost/any.hpp>
#endif

namespace CGAL {

template<class Refs>
class Vertex_base {

  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  typedef void* GenPtr;
  #else
  typedef boost::any GenPtr;
  #endif
  typedef typename Refs::Mark  Mark;
  typedef typename Refs::Point_3 Point_3;

  typedef typename Refs::Vertex_handle Vertex_handle;
  typedef typename Refs::SHalfloop_handle SHalfloop_handle;

  typedef typename Refs::Vertex_iterator Vertex_iterator;
  typedef typename Refs::SVertex_iterator SVertex_iterator;
  typedef typename Refs::SHalfedge_iterator SHalfedge_iterator;
  typedef typename Refs::SHalfloop_iterator SHalfloop_iterator;
  typedef typename Refs::SFace_iterator SFace_iterator;

  typedef typename Refs::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Refs::SVertex_const_iterator SVertex_const_iterator;
  typedef typename Refs::SHalfedge_const_iterator SHalfedge_const_iterator;
  typedef typename Refs::SFace_const_iterator SFace_const_iterator;
  typedef typename Refs::SHalfloop_const_handle SHalfloop_const_handle;

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

  Vertex_base() : point_at_center_(), mark_(), sncp_(),
    svertices_begin_(), svertices_last_(),
    shalfedges_begin_(), shalfedges_last_(),
    sfaces_begin_(), sfaces_last_(), shalfloop_(),
    info_()
    // , sm_(Vertex_handle((SNC_in_place_list_vertex<Vertex_base>*) this))
      {}

    Vertex_base(const Point_3& p, Mark m) :
      point_at_center_(p), mark_(m), sncp_(),
      svertices_begin_(), svertices_last_(),
      shalfedges_begin_(), shalfedges_last_(),
      sfaces_begin_(), sfaces_last_(), shalfloop_(),
      info_()
      //    , sm_(Vertex_handle((SNC_in_place_list_vertex<Vertex_base>*) this))
        {}

      Vertex_base(const Vertex_base<Refs>& v)
      //      : sm_(Vertex_handle((SNC_in_place_list_vertex<Vertex_base>*)this))
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

      Vertex_base<Refs>& operator=(const Vertex_base<Refs>& v)
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

      ~Vertex_base() {
        CGAL_NEF_TRACEN("  destroying Vertex item "<<&*this);
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
        return true;
      }

      bool is_valid( bool verb = false, int level = 0) const {

        Verbose_ostream verr(verb);
        verr << "begin CGAL::SNC_items<...>::Vertex_base::is_valid( verb=true, "
          "level = " << level << "):" << std::endl;

        bool valid = (sncp_ != nullptr);
        valid = valid && (svertices_begin_ != nullptr && svertices_begin_ != SVertex_iterator());
        valid = valid && (svertices_last_  != nullptr && svertices_last_  != SVertex_iterator());
        valid = valid && (shalfedges_begin_ != nullptr && shalfedges_begin_ != SHalfedge_iterator());
        valid = valid && (shalfedges_last_  != nullptr && shalfedges_last_  != SHalfedge_iterator());
        valid = valid && (sfaces_begin_ != nullptr && sfaces_begin_ != SFace_iterator());
        valid = valid && (sfaces_last_  != nullptr && sfaces_last_  != SFace_iterator());
        valid = valid && (shalfloop_ != nullptr && shalfloop_ != SHalfloop_iterator());

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

        verr << "end of CGAL::SNC_items<...>::Vertex_base::is_valid(): structure is "
             << ( valid ? "valid." : "NOT VALID.") << std::endl;
        return valid;
      }

}; // Vertex_base


} //namespace CGAL
#endif // CGAL_NEF_VERTEX_H
