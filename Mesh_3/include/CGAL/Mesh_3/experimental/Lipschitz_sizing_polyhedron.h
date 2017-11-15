// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Jane Tournois
//

#ifndef CGAL_LIPSCHITZ_SIZING_POLYHEDRON_H
#define CGAL_LIPSCHITZ_SIZING_POLYHEDRON_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/Mesh_3/experimental/Lipschitz_sizing_parameters.h>

#include <CGAL/Default.h>
#include <CGAL/array.h>
#include <CGAL/Bbox_3.h>

#include <boost/shared_ptr.hpp>

#include <list>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>

namespace CGAL
{
namespace Mesh_3
{

template <class Kernel, class MeshDomain
        , typename AABBTreeTemplate = CGAL::Default
>
class Lipschitz_sizing
{
public:
  typedef Kernel                      K;
  typedef typename Kernel::FT         FT;
  typedef typename Kernel::Triangle_3 Triangle;
  typedef typename Kernel::Point_3    Point_3;

  typedef typename std::list<Triangle>::iterator        Tr_iterator;
  typedef CGAL::AABB_triangle_primitive<K, Tr_iterator> Primitive;
  typedef CGAL::AABB_traits<K, Primitive>               AABB_tr_traits;
  typedef CGAL::AABB_tree<AABB_tr_traits>               AABB_tree;

  typedef typename CGAL::Default::Get<AABBTreeTemplate, AABB_tree>::type Tree;

  typedef typename MeshDomain::Index                    Index;
  typedef typename MeshDomain::Subdomain_index          Subdomain_index;
  typedef typename MeshDomain::Surface_patch_index      Surface_patch_index;

  typedef CGAL::Lipschitz_sizing_parameters<MeshDomain, FT> Parameters;


private:
  const Tree* m_ptree;
  boost::shared_ptr<Tree> m_own_ptree;

  const MeshDomain& m_domain;
  Parameters m_params;

public:
  Lipschitz_sizing(const MeshDomain& domain)
    : m_ptree(NULL)
    , m_own_ptree()
    , m_domain(domain)
    , m_params(domain)
  {
  }

  Lipschitz_sizing(const MeshDomain& domain
    , const Tree* ptree
    )
    : m_ptree(ptree)
    , m_own_ptree()
    , m_domain(domain)
    , m_params(domain)
  {
  }

  FT operator()(const Point_3& p, const int dim, const Index& index) const
  {
    CGAL_assertion(!m_params.empty());
    if (dim == 3)
    {
      return size_in_subdomain(p, m_domain.subdomain_index(index));
    }
    else if (dim == 2)
    {
      CGAL_assertion(false);//should not be used for dimension 2
      return 0.;
    }
    else if (dim == 1)
    {
      CGAL_assertion(false);//should not be used for dimension 1
      return 0.;
    }
    else if (dim == 0)
    {
      CGAL_assertion(false);//should not be used for dimension 0
      return 0;
    }
    CGAL_assertion(false);
    return 0.;
  }

public:
  template <typename C3T3>
  void init_aabb_tree_from_c3t3(const C3T3* p_c3t3)
  {
    static std::list<Triangle> triangles;
    for (typename C3T3::Facets_in_complex_iterator
      fit = p_c3t3->facets_in_complex_begin();
      fit != p_c3t3->facets_in_complex_end();
      ++fit)
    {
        triangles.push_back(p_c3t3->triangulation().triangle(*fit));
    }

    m_own_ptree.reset(new Tree(triangles.begin(), triangles.end()));
    m_own_ptree->build();
    m_own_ptree->accelerate_distance_queries();
  }

private:
  FT size_in_subdomain(const Point_3& p, const Subdomain_index& index) const
  {
    FT k, size_min, size_max;
    m_params.get_parameters(index, k, size_min, size_max);

    FT sqdist = 0.;
    if(m_ptree == NULL)
    {
      sqdist = m_own_ptree->squared_distance(p);
    }
    else
    {
      Point_3 closest = compute_closest_point(p);
      sqdist = CGAL::squared_distance(p, closest);
    }

    FT size = k * CGAL::sqrt(sqdist) + size_min;
    return (std::min)(size, size_max);
  }

private:
  Point_3 compute_closest_point(const Point_3& p) const
  {
    return m_ptree->closest_point(p);
  }

public:
  void add_parameters_for_subdomain(const Subdomain_index& id
                                  , const FT& k
                                  , const FT& size_min
                                  , const FT& size_max)
  {
    m_params.add_subdomain(id, k, size_min, size_max);
  }

};

}//namespace Mesh_3
}//namespace CGAL

#endif // CGAL_LIPSCHITZ_SIZING_POLYHEDRON_H
