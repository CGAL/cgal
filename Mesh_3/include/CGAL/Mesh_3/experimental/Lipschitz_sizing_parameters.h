// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois
//

#ifndef CGAL_LIP_SIZING_PARAMETERS_H
#define CGAL_LIP_SIZING_PARAMETERS_H

#include <CGAL/license/Mesh_3.h>

#include <map>
#include <limits>

namespace CGAL
{
template <typename MeshDomain, typename FT>
class Lipschitz_sizing_parameters
{
  struct SubdomainParam
  {
    FT m_k;
    FT m_size_min;//min size in subdomain
    FT m_size_max;//max size in subdomain

  public:
    SubdomainParam()
      : m_k(0.)
      , m_size_min(0.)
      , m_size_max(0.)
    {}
    SubdomainParam(const FT& k
                 , const FT& size_min
                 , const FT& size_max)
      : m_k(k)
      , m_size_min(size_min)
      , m_size_max(size_max)
    {}

    const FT& k() const                   { return m_k; }
    const FT& size_min() const       { return m_size_min; }
    const FT& size_max() const            { return m_size_max; }
  };

private:
  typedef typename MeshDomain::Subdomain_index     Subdomain_index;
  typedef typename MeshDomain::Surface_patch_index Surface_patch_index;

  typedef std::map<Subdomain_index, SubdomainParam> Parameters_map;
  Parameters_map m_parameters;
  SubdomainParam m_default_params;
  const MeshDomain* p_domain_;

public:
  Lipschitz_sizing_parameters(const MeshDomain& domain)
    : p_domain_(&domain)
  {}

  bool empty()       const {  return m_parameters.empty(); }
  std::size_t size() const { return m_parameters.size(); }

  void add_subdomain(const Subdomain_index& index
                   , const FT& k
                   , const FT& size_min
                   , const FT& size_max)
  {
    typename Parameters_map::iterator it = m_parameters.find(index);
    if (it != m_parameters.end())
      std::cout << "Warning : sizing parameters for subdomain " << index
                << "will be overwritten." << std::endl;

    if (index == -1)
      m_default_params
        = SubdomainParam(k, size_min, size_max);
    else
      m_parameters[index]
        = SubdomainParam(k, size_min, size_max);
  }

  void get_parameters(const Subdomain_index& index
                    , FT& k
                    , FT& size_min
                    , FT& size_max) const
  {
    typename Parameters_map::const_iterator it
      = m_parameters.find(index);

    SubdomainParam p = (it != m_parameters.end())
                     ? (*it).second
                     : m_default_params;
    k               = p.k();
    size_min        = p.size_min();
    size_max        = p.size_max();
  }

#ifdef CGAL_MESH_3_LIPSCHITZ_SIZING_EXPERIMENTAL
  void get_parameters(const Surface_patch_index& spi
                    , FT& size_max) const
  {
    const std::pair<Subdomain_index, Subdomain_index>& index
      = incident_subdomains(spi);

    typename Parameters_map::const_iterator it1
      = m_parameters.find(index.first);
    typename Parameters_map::const_iterator it2
      = m_parameters.find(index.second);

    SubdomainParam p1
      = (it1 != m_parameters.end()) ? (*it1).second : m_default_params;
    SubdomainParam p2
      = (it2 != m_parameters.end()) ? (*it2).second : m_default_params;

    bool boundary1 = (index.first == 0);   //boundary of the domain, inside the cube
    bool boundary2 = (index.second == 0);  //boundary of the domain, inside the cube

#ifdef CGAL_MESH_3_IMAGE_EXAMPLE
    boundary1 = boundary1 || (index.first == INT_MIN); //boundary of the cube
    boundary2 = boundary2 || (index.second == INT_MIN);//boundary of the cube
#endif
    CGAL_assertion(!boundary1 || !boundary2);

    if (!boundary1)
    {
      if (!boundary2)
        size_max = (std::min)(p1.size_min(), p2.size_min());
      else
        size_max = p1.size_min();
    }
    else
      size_max = p2.size_min();
  }

  const std::pair<Subdomain_index, Subdomain_index> &
    incident_subdomains(const Surface_patch_index& index) const
  {
//#ifdef CGAL_MESH_3_IMAGE_EXAMPLE
//    return index;
//#else //POLYHEDRAL_EXAMPLE
    return p_domain_->incident_subdomains_indices(index);
//#endif
  }
#endif

};
}//namespace CGAL

#endif //CGAL_LIP_SIZING_PARAMETERS_H
