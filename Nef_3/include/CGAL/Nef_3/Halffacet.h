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
#ifndef CGAL_NEF_HALFFACET_H
#define CGAL_NEF_HALFFACET_H

#include <CGAL/license/Nef_3.h>


#include <string>
#include <sstream>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/Nef_3/SNC_iteration.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 83
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template <typename Refs>
class Halffacet_base  {

  typedef typename Refs::Mark  Mark;
  typedef typename Refs::Plane_3   Plane_3;
  typedef typename Refs::Halffacet_handle         Halffacet_handle;
  typedef typename Refs::Halffacet_const_handle   Halffacet_const_handle;
  typedef typename Refs::Volume_handle            Volume_handle;
  typedef typename Refs::Volume_const_handle      Volume_const_handle;
  typedef typename Refs::Object_list    Object_list;
  typedef typename Refs::Halffacet_cycle_iterator
    Halffacet_cycle_iterator;
  typedef typename Refs::Halffacet_cycle_const_iterator
    Halffacet_cycle_const_iterator;

  Plane_3              supporting_plane_;
  Mark                 mark_;
  Halffacet_handle     twin_;
  Volume_handle        volume_;
  Object_list          boundary_entry_objects_; // SEdges, SLoops

 public:

  Halffacet_base() : supporting_plane_(), mark_() {}

    Halffacet_base(const Plane_3& h, Mark m) :
      supporting_plane_(h), mark_(m) {}

      ~Halffacet_base() {
        CGAL_NEF_TRACEN("  destroying Halffacet_base item "<<&*this);
      }

      Halffacet_base(const Halffacet_base<Refs>& f)
        { supporting_plane_ = f.supporting_plane_;
          mark_ = f.mark_;
          twin_ = f.twin_;
          CGAL_NEF_TRACEN("VOLUME const");
          volume_ = f.volume_;
          boundary_entry_objects_ = f.boundary_entry_objects_;
        }

      Halffacet_base<Refs>& operator=(const Halffacet_base<Refs>& f)
        { if (this == &f) return *this;
          supporting_plane_ = f.supporting_plane_;
          mark_ = f.mark_;
          twin_ = f.twin_;
          CGAL_NEF_TRACEN("VOLUME op=");
          volume_ = f.volume_;
          boundary_entry_objects_ = f.boundary_entry_objects_;
          return *this;
        }

      Mark& mark() { return mark_; }
      const Mark& mark() const { return mark_; }

      Halffacet_handle& twin() { return twin_; }
      Halffacet_const_handle twin() const { return twin_; }

      Plane_3& plane() { return supporting_plane_; }
      const Plane_3& plane() const { return supporting_plane_; }

      Volume_handle& incident_volume() { return volume_; }
      Volume_const_handle incident_volume() const { return volume_; }

      Object_list& boundary_entry_objects() { return boundary_entry_objects_; }
      const Object_list& boundary_entry_objects() const { return boundary_entry_objects_; }

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
        verr << "begin CGAL::SNC_items<...>::Halffacet_base::is_valid( verb=true, "
          "level = " << level << "):" << std::endl;

        bool valid = (twin_ != nullptr && twin_ != Halffacet_handle());
        valid = valid && (volume_ != nullptr && volume_ != Volume_handle());

        valid = valid && (supporting_plane_.a() != 0 ||
                          supporting_plane_.b() != 0 ||
                          supporting_plane_.c() != 0);

        valid = valid && (!boundary_entry_objects_.empty());

        verr << "end of CGAL::SNC_items<...>::Halffacet_base::is_valid(): structure is "
             << ( valid ? "valid." : "NOT VALID.") << std::endl;

        return valid;
      }

}; // Halffacet_base

} //namespace CGAL
#endif //CGAL_NEF_HALFFACET_H
