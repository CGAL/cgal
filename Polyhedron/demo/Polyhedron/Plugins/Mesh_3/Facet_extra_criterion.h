// Copyright (c) 2015 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau
#ifndef CGAL_FACET_EXTRA_CRITERION_H
#define CGAL_FACET_EXTRA_CRITERION_H

#include <CGAL/Mesh_3/mesh_standard_criteria.h>


template <typename Tr, typename Domain, typename Visitor_>
class Facet_extra_criterion :
  public CGAL::Mesh_3::Abstract_criterion<Tr, Visitor_>
{
private:
  typedef typename Tr::Facet Facet;

  typedef CGAL::Mesh_3::Abstract_criterion<Tr,Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Badness Badness;

  typedef Facet_extra_criterion<Tr,Domain,Visitor_> Self;

  const Domain& domain;

public:
  /// Constructor
  Facet_extra_criterion(const Domain& domain) : domain(domain) {}
  /// Destructor
  ~Facet_extra_criterion() {}

protected:
  virtual void do_accept(Visitor_& v) const
  {
    v.visit(*this);
  }

  virtual Self* do_clone() const
  {
    // Call copy ctor on this
    return new Self(*this);
  }


  virtual Badness do_is_bad (const Tr& /*tr*/, const Facet& f) const
  {
    typedef typename Tr::Vertex_handle  Vertex_handle;
    typedef typename Tr::Cell_handle    Cell_handle;
    typedef typename Domain::Surface_patch_index  Surface_patch_index;

    const Cell_handle& ch = f.first;
    const int& i = f.second;

    const Vertex_handle& v1 = ch->vertex((i+1)&3);
    const Vertex_handle& v2 = ch->vertex((i+2)&3);
    const Vertex_handle& v3 = ch->vertex((i+3)&3);

    Surface_patch_index index = Surface_patch_index();
    bool is_index_initialized = false;

    std::cerr << typeid(index).name() << std::endl;
    if ( v1->in_dimension() == 2 )
    {
      index = domain.surface_patch_index(v1->index());
      if(index != 0) {    //   (index.first != 0 && index.second != 0)
        //index = Surface_patch_index(std::make_pair(1,1));
      }
      is_index_initialized = true;
    }

    if ( v2->in_dimension() == 2 )
    {
      Surface_patch_index index2 = domain.surface_patch_index(v2->index());
      if(index2) { // (index2.first != 0 && index2.second != 0){
        //index2 = Surface_patch_index(1,1);
      }
      if ( is_index_initialized )
      {
        if ( !(index2 == index) )
        {
          return Badness(Quality(1));
        }
      }
      else
      {
        index = index2;
        is_index_initialized = true;
      }
    }

    if ( v3->in_dimension() == 2 )
    {
      Surface_patch_index index3 = domain.surface_patch_index(v3->index());
      if(index3 != 0) { // (index3.first != 0 && index3.second != 0)
          // index3 = Surface_patch_index(1,1);
      }
      if ( is_index_initialized && !(index3 == index) )
      {
        return Badness(Quality(1));
      }
    }

    return  Badness();
  }

}; // end class Facet_extra_criterion


#endif // CGAL_FACET_EXTRA_CRITERION_H
