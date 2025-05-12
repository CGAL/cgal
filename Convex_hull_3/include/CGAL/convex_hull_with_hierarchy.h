// Copyright (c) 2025 Geometry Factory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Léo Valque
//

#ifndef CGAL_CONVEX_HULL_WITH_HIERARCHY_3_H
#define CGAL_CONVEX_HULL_WITH_HIERARCHY_3_H

#include <CGAL/license/Convex_hull_3.h>

#include <CGAL/Surface_mesh.h>

#include <vector>

namespace CGAL{

std::array<size_t, 8> nb_visited_in_hierarchy;

template <class P>
struct Convex_hull_with_hierarchy{
  // parameterization of the hierarchy
  constexpr static size_t RATIO = 32;
  constexpr static size_t MINSIZE_FOR_NEXT_HIERARCHY = RATIO*4;
  constexpr static size_t MAXSIZE_FOR_NAIVE_SEARCH = RATIO;

  typedef Surface_mesh< P > Mesh;
  typedef typename Mesh::Vertex_index Vertex_index;

//   const Base &sm;
  std::vector< Mesh > hierarchy_sm;

  std::size_t maxlevel() const{
    return hierarchy_sm.size()-1;
  }

  const Mesh& hierarchy_mesh(std::size_t level) const{
    return hierarchy_sm[level];
  }

  void init_hierarchy(){
    size_t size=hierarchy_sm[0].vertices().size();
    size_t level=0;

    if(size<=MINSIZE_FOR_NEXT_HIERARCHY)
      return;

    Random r;
    while(size>MINSIZE_FOR_NEXT_HIERARCHY){
      std::vector<P> select_points;
      std::map<P, Vertex_index> select_vertices;
      select_points.reserve(2*size/RATIO);
      Mesh& above_sm=hierarchy_sm[level];
      ++level;

      for(Vertex_index v: above_sm.vertices())
        if(r.get_int(0,RATIO-1)==0){
          select_points.push_back(above_sm.point(v));
          select_vertices[above_sm.point(v)]=v;
        }

      Mesh new_sm;
      convex_hull_3(select_points.begin(), select_points.end(), new_sm);
      typename Mesh::Property_map<Vertex_index, Vertex_index> pm=new_sm.template add_property_map<Vertex_index>("v:next_in_hierarchy", *(above_sm.vertices().begin())).first;
      for(Vertex_index v : new_sm.vertices())
        pm[v]=select_vertices[new_sm.point(v)];
      size=new_sm.vertices().size();

      hierarchy_sm.push_back(std::move(new_sm));
    }
  }

  Convex_hull_with_hierarchy(Mesh &sm_){
    hierarchy_sm.push_back(sm_);
    init_hierarchy();
  };

  template <class IK, class Converter, class Vector_3>
  const typename IK::Point_3 extreme_point(const Vector_3 dir, const Converter& converter) const {
    using Point_3 = typename Kernel_traits<Vector_3>::Kernel::Point_3;
    using FT= typename Kernel_traits<Vector_3>::Kernel::FT;

    size_t level=maxlevel();

    const Mesh &sm = hierarchy_sm[maxlevel()];

    Vertex_index argmax=*sm.vertices().begin();
    FT tmax=Vector_3(ORIGIN, converter(sm.point(argmax)))*dir;

    if(sm.vertices().size() <= MAXSIZE_FOR_NAIVE_SEARCH){
      //If maxlevel is small, we simply go through all its vertices
      for(auto vh=++(sm.vertices().begin()); vh!=sm.vertices().end(); ++vh){
        Vertex_index v=*vh;
        ++(nb_visited_in_hierarchy[level]);
        FT p=Vector_3(ORIGIN, converter(sm.point(v)))*dir;
        if(compare(tmax, p)==SMALLER){
          tmax=p;
          argmax=v;
        }
      }

      // Go to under level
      if(level>0){
        argmax=((sm.template property_map<Vertex_index, Vertex_index>("v:next_in_hierarchy")).first)[argmax];
        --level;
      }
    }

    for(; level>=0; --level){
      // Starting from the vertex of the previous level, we walk on the graph
      // along neighbor that increase the "score"
      const Mesh &csm = hierarchy_mesh(level);
      while(true){
        loop_walk_on_the_graph:;
        for(Vertex_index v: vertices_around_target(argmax ,csm)){
          ++(nb_visited_in_hierarchy[level]);
          FT p=Vector_3(ORIGIN, converter(csm.point(v)))*dir;
          if(compare(tmax, p)==SMALLER){
            tmax=p;
            argmax=v;
            goto loop_walk_on_the_graph; // repeat with the new vertex
          }
        }
        break;
      }
      if(level>0)
        argmax=((csm.template property_map<Vertex_index, Vertex_index>("v:next_in_hierarchy")).first)[argmax];
      else
        return csm.point(argmax);
    }
  }
};

template <class IK, class Converter, class Vector_3>
const typename IK::Point_3 extreme_point(const Convex_hull_with_hierarchy<typename IK::Point_3> &C, const Vector_3 dir, const Converter &c){
  return C.template extreme_point<IK, Converter>(dir, c);
}

}

#endif