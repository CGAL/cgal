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
#ifndef CGAL_NEF_SHALFLOOP_H
#define CGAL_NEF_SHALFLOOP_H

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
class SHalfloop_base {

  typedef typename Refs::Mark  Mark;
  typedef typename Refs::Sphere_circle  Sphere_circle;
  typedef typename Refs::SHalfloop_handle SHalfloop_handle;
  typedef typename Refs::SHalfloop_const_handle SHalfloop_const_handle;
  typedef typename Refs::SFace_handle SFace_handle;
  typedef typename Refs::SFace_const_handle SFace_const_handle;
  typedef typename Refs::Halffacet_handle Halffacet_handle;
  typedef typename Refs::Halffacet_const_handle Halffacet_const_handle;

  SHalfloop_handle   twin_;
  SFace_handle       incident_sface_;
  Halffacet_handle   facet_;

  // temporary needed:
  Mark               mark_;
  Sphere_circle      circle_;

 public:

  SHalfloop_base() : twin_(), incident_sface_(), facet_(),
    mark_(), circle_() {}

    ~SHalfloop_base() {
      CGAL_NEF_TRACEN("  destroying SHalfloop_base item "<<&*this);
    }
    SHalfloop_base(const SHalfloop_base<Refs>& l)
      { twin_ = l.twin_;
        incident_sface_ = l.incident_sface_;
        facet_ = l.facet_;
        mark_ = l.mark_;
        circle_ = l.circle_;
      }

    SHalfloop_base<Refs>& operator=(const SHalfloop_base<Refs>& l)
      { twin_ = l.twin_;
        incident_sface_ = l.incident_sface_;
        facet_ = l.facet_;
        mark_ = l.mark_;
        circle_ = l.circle_;
        return *this;
      }

    SHalfloop_base<Refs>& operator=(SHalfloop_base<Refs>&& l) noexcept
      { twin_ = std::move(l.twin_);
        incident_sface_ = std::move(l.incident_sface_);
        facet_ = std::move(l.facet_);
        mark_ = std::move(l.mark_);
        circle_ = std::move(l.circle_);
        return *this;
      }

    Mark& mark() { return mark_;}
    const Mark& mark() const { return mark_; }

    SHalfloop_handle& twin() { return twin_; }
    SHalfloop_const_handle twin() const { return twin_; }

    Sphere_circle& circle() { return circle_; }
    const Sphere_circle& circle() const { return circle_; }

    SFace_handle& incident_sface() { return incident_sface_; }
    SFace_const_handle incident_sface() const { return incident_sface_; }

    Halffacet_handle& facet() { return facet_; }
    Halffacet_const_handle facet() const { return facet_; }

 public:
    std::string debug() const
      { std::stringstream os;
        CGAL::IO::set_pretty_mode(os);
        os<<"sl [ "<<circle_<<" ] ";
        return os.str();
      }

    bool is_twin() const { return (&*twin_ < this); }

    bool is_valid( bool verb = false, int level = 0) const {

      Verbose_ostream verr(verb);
      verr << "begin CGAL::SNC_items<...>::SHalfloop_base::is_valid( verb=true, "
        "level = " << level << "):" << std::endl;

      bool valid = (twin_  != SHalfloop_handle() && twin_  != nullptr);
      valid = valid && (incident_sface_ != SFace_handle() &&
                        incident_sface_ != nullptr);
      valid = valid && (facet_ != Halffacet_handle() &&
                        facet_ != nullptr);
      valid = valid && (circle_.d() == 0);
      valid = valid && (circle_.a() != 0 || circle_.b() != 0 || circle_.c() !=0);

      verr << "end of CGAL::SNC_items<...>::SHalfloop_base::is_valid(): structure is "
           << ( valid ? "valid." : "NOT VALID.") << std::endl;

      return valid;
    }

}; // SHalfloop_base

} //namespace CGAL
#endif //CGAL_NEF_SHALFLOOP_H
