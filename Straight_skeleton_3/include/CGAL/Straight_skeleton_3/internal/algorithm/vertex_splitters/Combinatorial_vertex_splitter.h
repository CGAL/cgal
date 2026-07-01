// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * file   algo/3d/CombiVertexSplitter.h
 * author Gernot Walzl
 * date   2012-07-26
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_COMBINATORIAL_VERTEX_SPLITTER_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_COMBINATORIAL_VERTEX_SPLITTER_H

#include <CGAL/license/Straight_skeleton_3.h>

#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/Straight_skeleton_3.h>

#include <boost/shared_array.hpp>

#include <list>
#include <memory>
#include <vector>
#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename GeomTraits>
class Combi_vertex_splitter
  : public Abstract_vertex_splitter<GeomTraits>
{
  using Base = Abstract_vertex_splitter<GeomTraits>;
  using Combi_vertex_splitter_sptr = std::shared_ptr<Combi_vertex_splitter>;

private:
  using Polyhedron = HDS::Polyhedron<GeomTraits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using Vertex = typename Polyhedron::Vertex;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using Edge = typename Polyhedron::Edge;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using Facet = typename Polyhedron::Facet;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

  using Skeleton_vertex_data = typename Polyhedron::Skeleton_vertex_data;
  using SkelVertexDataSPtr = typename Polyhedron::SkelVertexDataSPtr;
  using Skeleton_facet_data = typename Polyhedron::Skeleton_facet_data;
  using SkelFacetDataSPtr = typename Polyhedron::SkelFacetDataSPtr;

private:
  using Straight_skeleton_3 = CGAL::Straight_skeleton_3<GeomTraits>;

  using NodeWPtr = typename Straight_skeleton_3::NodeWPtr;
  using NodeSPtr = typename Straight_skeleton_3::NodeSPtr;

private:
  using vec2i = boost::shared_array<int>;
  using combi = std::vector<vec2i>;

public:
  Combi_vertex_splitter()
  {
    this->type_ = Base::COMBI_VERTEX_SPLITTER;
    ConfigurationSPtr config = Configuration::get_instance();
    selected_combi_ = config->get_int("Algorithm", "selected_combinatorial_split");
  }

  virtual ~Combi_vertex_splitter() { /*intentionally does nothing*/ }

  static Combi_vertex_splitter_sptr create()
  {
    return std::make_shared<Combi_vertex_splitter>();
  }

  static vec2i create_split(int begin, int end)
  {
    vec2i result(new int[2]);
    result[0] = begin;
    result[1] = end;
    return result;
  }

  static int compare_splits(const vec2i& split1, const vec2i& split2)
  {
    int result = 0;
    if (split1[0] < split2[0] || (split1[0] == split2[0] && split1[1] < split2[1])) {
      result = 1;
    } else if (split1[0] > split2[0] || (split1[0] == split2[0] && split1[1] > split2[1])) {
      result = -1;
    }
    return result;
  }

  static std::vector<int> init_labels(unsigned int degree)
  {
    std::vector<int> result;
    for (unsigned int i = 0; i < degree; ++i) {
      result.push_back(i);
    }
    return result;
  }

  static std::vector<int> split_labels(std::vector<int>& labels,
                                       const vec2i& split)
  {
    std::vector<int> result;
    int begin = split[0];
    int end = split[1];
    bool inside = false;
    std::vector<int>::const_iterator it = labels.begin();
    while (it != labels.end()) {
      std::vector<int>::const_iterator it_current = it;
      int label = *it++;
      if (label == begin) {
        inside = true;
      }
      if (inside) {
        result.push_back(label);
        if (label != begin && label != end) {
          it = labels.erase(it_current);
        }
      }
      if (label == end) {
        inside = false;
      }
      if (inside) {
        if (it == labels.end()) {
          it = labels.begin();
        }
      }
    }
    return result;
  }

  static std::list<vec2i> create_single_split_combiations(const std::vector<int>& labels)
  {
    std::list<vec2i> result;
    unsigned int degree = labels.size();
    for (unsigned int i = 0; i < degree-1; ++i) {
      for (unsigned int j = i+2; j < degree; ++j) {
        if (i == 0 && j == degree-1) {
          continue;
        }
        vec2i split = create_split(labels[i], labels[j]);
        result.push_back(split);
      }
    }
    return result;
  }

  static std::list<combi> append_split_combinations(const combi& history,
                                                    const std::list<vec2i>& splits)
  {
    std::list<combi> result;
    if (history.size() == 0) {
      std::list<vec2i>::const_iterator it_splits = splits.begin();
      while (it_splits != splits.end()) {
        vec2i split = *it_splits++;
        combi combination;
        combination.push_back(split);
        result.push_back(combination);
      }
    } else {
      vec2i last_split = history.back();
      std::list<vec2i>::const_iterator it_splits = splits.begin();
      while (it_splits != splits.end()) {
        vec2i split = *it_splits++;
        if (compare_splits(last_split, split) > 0) {
          combi combination(history);
          combination.push_back(split);
          result.push_back(combination);
        }
      }
    }
    return result;
  }

  static std::list<combi> merge_combinations(const combi& history,
                                             const std::list<combi>& combis1,
                                             const std::list<combi>& combis2)
  {
    std::list<combi> result;
    unsigned int history_size = history.size();
    std::list<combi>::const_iterator it_combis1 = combis1.begin();
    while (it_combis1 != combis1.end()) {
      combi combi1 = *it_combis1++;
      std::list<combi>::const_iterator it_combis2 = combis2.begin();
      while (it_combis2 != combis2.end()) {
        combi combi2 = *it_combis2++;
        combi combi_merged(history);
        std::vector<vec2i>::iterator it_combi1 = combi1.begin();
        std::vector<vec2i>::iterator it_combi2 = combi2.begin();
        for (unsigned int i = 0; i < history_size; ++i) {
          it_combi1++;
          it_combi2++;
        }
        while (it_combi1 != combi1.end() || it_combi2 != combi2.end()) {
          if (it_combi1 == combi1.end()) {
            combi_merged.push_back(*it_combi2++);
            continue;
          }
          if (it_combi2 == combi2.end()) {
            combi_merged.push_back(*it_combi1++);
            continue;
          }
          vec2i split1 = *it_combi1;
          vec2i split2 = *it_combi2;
          if (compare_splits(split1, split2) > 0) {
            combi_merged.push_back(split1);
            it_combi1++;
          } else {
            combi_merged.push_back(split2);
            it_combi2++;
          }
        }
        result.push_back(combi_merged);
      }
    }
    return result;
  }

  static std::list<combi> generate_combinations_rec(const combi& history,
                                                    const std::vector<int>& labels)
  {
    std::list<combi> result;
    std::list<vec2i> splits = create_single_split_combiations(labels);
    std::list<combi> combis = append_split_combinations(history, splits);
    if (labels.size() <= 4) {
      result = combis;
    } else {
      std::list<combi>::iterator it_combis = combis.begin();
      while (it_combis != combis.end()) {
        combi split_combi = *it_combis++;
        std::list<combi> combis_next;
        vec2i split = split_combi.back();
        std::vector<int> labels1(labels);
        std::vector<int> labels2 = split_labels(labels1, split);
        if (labels1.size() > 3 && labels2.size() > 3) {
          std::list<combi> combis1 = generate_combinations_rec(split_combi, labels1);
          std::list<combi> combis2 = generate_combinations_rec(split_combi, labels2);
          combis_next = merge_combinations(split_combi, combis1, combis2);
        } else if (labels1.size() > 3) {
          combis_next = generate_combinations_rec(split_combi, labels1);
        } else if (labels2.size() > 3) {
          combis_next = generate_combinations_rec(split_combi, labels2);
        }
        result.insert(result.end(), combis_next.begin(), combis_next.end());
      }
    }
    return result;
  }

  static std::list<combi> generate_all_combinations(unsigned int degree)
  {
    combi history;
    std::vector<int> labels = init_labels(degree);
    std::list<combi> result = generate_combinations_rec(history, labels);
    // CGAL_SS3_SPLITTER_TRACE("degree = " << degree);
    // CGAL_SS3_SPLITTER_TRACE_CODE(for (const combi& combination : result))
    // CGAL_SS3_SPLITTER_TRACE("  " << combi_to_string(combination));
    return result;
  }

  static PolyhedronSPtr copy_vertex(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    PolyhedronSPtr result = Polyhedron::create();
    VertexSPtr vertex_c = vertex->clone();
    result->add_vertex(vertex_c);
    FacetSPtr facet = vertex->first_facet();
    EdgeSPtr edge_first = vertex->find_edge(facet);
    EdgeSPtr edge;
    EdgeSPtr edge_prev;
    EdgeSPtr edge_prev_c;
    EdgeSPtr edge_first_c;
    FacetSPtr facet_next = facet;
    FacetSPtr facet_c;
    SkelFacetDataSPtr facet_c_data;
    while (edge != edge_first) {
      if (!edge) {
        edge = edge_first;
      }
      EdgeSPtr edge_c;
      VertexSPtr vertex_tgt_c;
      if (edge->source() == vertex) {
        vertex_tgt_c = edge->target()->clone();
      } else if (edge->target() == vertex) {
        vertex_tgt_c = edge->source()->clone();
      }
      result->add_vertex(vertex_tgt_c);
      edge_c = Edge::create(vertex_c, vertex_tgt_c);
      if (edge == edge_first) {
        edge_first_c = edge_c;
      }
      result->add_edge(edge_c);
      if (edge_prev_c) {
        facet_c = Facet::create();
        facet_c->set_plane(facet->get_plane());
        edge_prev_c->set_facet_L(facet_c);
        facet_c->add_edge(edge_prev_c);
        edge_c->set_facet_R(facet_c);
        facet_c->add_edge(edge_c);
        facet_c_data = Skeleton_facet_data::create(facet_c);
        facet_c_data->set_facet_origin(facet);
        if (facet->has_data()) {
          facet_c_data->set_speed(std::dynamic_pointer_cast<Skeleton_facet_data>(facet->get_data())->get_speed());
        }
        result->add_facet(facet_c);
      }
      edge_prev_c = edge_c;
      edge_prev = edge;
      edge = edge->next(vertex);
      facet = facet_next;
      facet_next = edge->other(facet);
    }
    facet_c = Facet::create();
    facet_c->set_plane(facet->get_plane());
    edge_prev_c->set_facet_L(facet_c);
    facet_c->add_edge(edge_prev_c);
    edge_first_c->set_facet_R(facet_c);
    facet_c->add_edge(edge_first_c);
    facet_c_data = Skeleton_facet_data::create(facet_c);
    facet_c_data->set_facet_origin(facet);
    if (facet->has_data()) {
        facet_c_data->set_speed(std::dynamic_pointer_cast<Skeleton_facet_data>(facet->get_data())->get_speed());
    }
    result->add_facet(facet_c);
    result->initialize_all_IDs();
    return result;
  }

  static PolyhedronSPtr split_vertex(const VertexSPtr& vertex,
                                     const combi& combination)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_precondition(vertex->degree() - 3 == combination.size());
    PolyhedronSPtr polyhedron = vertex->get_polyhedron();
    std::list<VertexSPtr> vertices_tosplit;
    vertices_tosplit.push_back(vertex);

    std::vector<FacetSPtr> facets;
    FacetSPtr facet = vertex->first_facet();
    EdgeSPtr edge_first = vertex->find_edge(facet);
    EdgeSPtr edge;
    while (edge != edge_first) {
      if (!edge) {
        edge = edge_first;
      }
      facets.push_back(facet);
      edge = edge->next(vertex);
      facet = edge->other(facet);
    }

    std::list<EdgeSPtr> edges_toremove;
    for (unsigned int i = 0; i < combination.size(); ++i) {
      vec2i split = combination[i];
      FacetSPtr facet_right = facets[split[0]];
      FacetSPtr facet_left = facets[split[1]];
      typename std::list<VertexSPtr>::iterator it_v = vertices_tosplit.begin();
      while (it_v != vertices_tosplit.end()) {
        typename std::list<VertexSPtr>::iterator it_current = it_v;
        VertexSPtr vertex = *it_v++;
        if (vertex->has_incident_facet(facet_right) &&
            vertex->has_incident_facet(facet_left)) {
          vertices_tosplit.erase(it_current);
          VertexSPtr vertex2 = vertex->split(facet_left, facet_right);
          if (facet_left->get_plane() == facet_right->get_plane()) {
            EdgeSPtr edge = vertex->find_edge(vertex2);
            edges_toremove.push_back(edge);
          }
          if (vertex->has_data()) {
            SkelVertexDataSPtr data = std::dynamic_pointer_cast<Skeleton_vertex_data>(vertex->get_data());
            SkelVertexDataSPtr data2 = Skeleton_vertex_data::create(vertex2);
            data2->set_wnode(data->get_wnode());
          }
          if (vertex->degree() > 3) {
            vertices_tosplit.push_back(vertex);
          }
          if (vertex2->degree() > 3) {
            vertices_tosplit.push_back(vertex2);
          }
          break;
        }
      }
    }

    typename std::list<EdgeSPtr>::iterator it_e = edges_toremove.begin();
    while (it_e != edges_toremove.end()) {
      EdgeSPtr edge = *it_e++;
      VertexSPtr vertex_src = edge->source();
      VertexSPtr vertex_tgt = edge->target();
      FacetSPtr facet_l = edge->get_facet_L();
      FacetSPtr facet_r = edge->get_facet_R();
      facet_l->remove_edge(edge);
      facet_r->remove_edge(edge);
      polyhedron->remove_edge(edge);
      if (facet_l != facet_r) {
        facet_l->merge(facet_r);
        polyhedron->remove_facet(facet_r);
      }

      // remove vertices of degree 2
      for (unsigned int i = 0; i < 2; ++i) {
        VertexSPtr vertex;
        if (i == 0) {
          vertex = vertex_src;
        } else if (i == 1) {
          vertex = vertex_tgt;
        }
        VertexSPtr vertex_merged_src = vertex->prev(facet_l);
        VertexSPtr vertex_merged_tgt = vertex->next(facet_l);
        FacetSPtr facet_r = facet_l->next(vertex);
        typename std::list<EdgeWPtr>::iterator it_ew = vertex->edges().begin();
        while (it_ew != vertex->edges().end()) {
          EdgeWPtr edge_wptr = *it_ew++;
          if (EdgeSPtr edge_toremove = edge_wptr.lock()) {
            facet_l->remove_edge(edge_toremove);
            facet_r->remove_edge(edge_toremove);
            polyhedron->remove_edge(edge_toremove);
          }
        }
        typename std::list<FacetWPtr>::iterator it_fw = vertex->facets().begin();
        while (it_fw != vertex->facets().end()) {
          FacetWPtr facet_wptr = *it_fw++;
          if (FacetSPtr facet = facet_wptr.lock()) {
            facet->remove_vertex(vertex);
          }
        }
        polyhedron->remove_vertex(vertex);
        EdgeSPtr edge_merged = Edge::create(vertex_merged_src, vertex_merged_tgt);
        edge_merged->set_facet_L(facet_l);
        edge_merged->set_facet_R(facet_r);
        facet_l->add_edge(edge_merged);
        facet_r->add_edge(edge_merged);
        polyhedron->add_edge(edge_merged);
      }
    }
    return polyhedron;
  }

  static PolyhedronSPtr apply(PolyhedronSPtr poly_split,
                              const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(poly_split);
    CGAL_SS3_DEBUG_SPTR(vertex);
    PolyhedronSPtr polyhedron = vertex->get_polyhedron();
    NodeWPtr node;
    if (vertex->has_data()) {
      SkelVertexDataSPtr data = std::dynamic_pointer_cast<Skeleton_vertex_data>(vertex->get_data());
      node = data->get_wnode();
    }
    std::map<VertexSPtr, VertexSPtr> vertices;
    typename std::list<VertexSPtr>::iterator it_v = poly_split->vertices().begin();
    while (it_v != poly_split->vertices().end()) {
      VertexSPtr vertex_ps = *it_v++;
      if (vertex_ps->degree() > 1) {
        VertexSPtr vertex_vs = Vertex::create(vertex_ps->point());
        vertices[vertex_ps] = vertex_vs;
        polyhedron->add_vertex(vertex_vs);
        SkelVertexDataSPtr data = Skeleton_vertex_data::create(vertex_vs);
        data->set_wnode(node);
      }
    }

    typename std::list<EdgeSPtr>::iterator it_e = poly_split->edges().begin();
    while (it_e != poly_split->edges().end()) {
      EdgeSPtr edge_ps = *it_e++;
      VertexSPtr vertex_ps_src = edge_ps->source();
      VertexSPtr vertex_ps_tgt = edge_ps->target();
      FacetSPtr facet_ps_l = edge_ps->get_facet_L();
      FacetSPtr facet_ps_r = edge_ps->get_facet_R();
      if (vertex_ps_tgt->degree() == 1) {
        EdgeSPtr edge_vs;
        typename std::list<EdgeWPtr>::iterator it_ve = vertex->edges().begin();
        while (it_ve != vertex->edges().end()) {
          EdgeWPtr edge_wptr = *it_ve++;
          if (EdgeSPtr edge = edge_wptr.lock()) {
            if (edge->source()->point() == vertex_ps_tgt->point() ||
                edge->target()->point() == vertex_ps_tgt->point()) {
              edge_vs = edge;
              break;
            }
          }
        }
        VertexSPtr vertex_vs = vertices[vertex_ps_src];
        if (edge_vs->source() == vertex) {
          edge_vs->replace_vertex_src(vertex_vs);
        } else {
          edge_vs->replace_vertex_tgt(vertex_vs);
        }
        if (!edge_vs->get_facet_L()->contains_vertex(vertex_vs)) {
          edge_vs->get_facet_L()->add_vertex(vertex_vs);
        }
        if (!edge_vs->get_facet_R()->contains_vertex(vertex_vs)) {
          edge_vs->get_facet_R()->add_vertex(vertex_vs);
        }
      } else if (vertex_ps_src->degree() == 1) {
        EdgeSPtr edge_vs;
        typename std::list<EdgeWPtr>::iterator it_ve = vertex->edges().begin();
        while (it_ve != vertex->edges().end()) {
          EdgeWPtr edge_wptr = *it_ve++;
          if (EdgeSPtr edge = edge_wptr.lock()) {
            if (edge->source()->point() == vertex_ps_src->point() ||
                edge->target()->point() == vertex_ps_src->point()) {
              edge_vs = edge;
              break;
            }
          }
        }
        VertexSPtr vertex_vs = vertices[vertex_ps_tgt];
        if (edge_vs->source() == vertex) {
          edge_vs->replace_vertex_src(vertex_vs);
        } else {
          edge_vs->replace_vertex_tgt(vertex_vs);
        }
        if (!edge_vs->get_facet_L()->contains_vertex(vertex_vs)) {
          edge_vs->get_facet_L()->add_vertex(vertex_vs);
        }
        if (!edge_vs->get_facet_R()->contains_vertex(vertex_vs)) {
          edge_vs->get_facet_R()->add_vertex(vertex_vs);
        }
      } else {
        EdgeSPtr edge_vs = Edge::create(vertices[vertex_ps_src], vertices[vertex_ps_tgt]);
        SkelFacetDataSPtr data_l = std::dynamic_pointer_cast<Skeleton_facet_data>(facet_ps_l->get_data());
        SkelFacetDataSPtr data_r = std::dynamic_pointer_cast<Skeleton_facet_data>(facet_ps_r->get_data());
        FacetSPtr facet_vs_l = data_l->get_facet_origin();
        FacetSPtr facet_vs_r = data_r->get_facet_origin();
        edge_vs->set_facet_L(facet_vs_l);
        edge_vs->set_facet_R(facet_vs_r);
        facet_vs_l->add_edge(edge_vs);
        facet_vs_r->add_edge(edge_vs);
        polyhedron->add_edge(edge_vs);
      }
    }

    typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
    while (it_f != vertex->facets().end()) {
      FacetWPtr facet_wptr = *it_f++;
      if (FacetSPtr facet = facet_wptr.lock()) {
        facet->remove_vertex(vertex);
      }
    }
    polyhedron->remove_vertex(vertex);
    return polyhedron;
  }

  virtual bool split_vertex(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    PolyhedronSPtr polyhedron = vertex->get_polyhedron();
    if (vertex->degree() <= 3) {
      return true;
    }
    vertex->sort();
    std::list<combi> combinations = generate_all_combinations(vertex->degree());
    std::list<combi> combinations_valid;
    std::list<PolyhedronSPtr> polys_split;
    std::list<combi>::iterator it_combi = combinations.begin();
    while (it_combi != combinations.end()) {
      combi combination = *it_combi++;
      PolyhedronSPtr poly_c = copy_vertex(vertex);
      VertexSPtr vertex_c = poly_c->vertices().front();
      split_vertex(vertex_c, combination);
      if (Base::check_splitted(poly_c)) {
        CGAL_SS3_SPLITTER_TRACE_V(16, "Valid split-combination found: " << combi_to_string(combination));
        combinations_valid.push_back(combination);
        polys_split.push_back(poly_c);
        CGAL_SS3_SPLITTER_TRACE_V(16, "Found valid combination");
        // 'selected_combi_' is usually 0 => stop on the first valid combination
        if (polys_split.size() > selected_combi_) {
          break;
        }
      }
    }

    CGAL_SS3_SPLITTER_TRACE_V(16, "Valid split-combination #: " << combinations_valid.size());

    if (combinations_valid.empty()) {
      CGAL_SS3_SPLITTER_TRACE_V(1, "Error: No valid split-combination found!");
      return false;
    }

    CGAL_assertion(!polys_split.empty());

    unsigned int selected_combinatorial_split = selected_combi_ % combinations_valid.size();
    it_combi = combinations_valid.begin();
    typename std::list<PolyhedronSPtr>::iterator it_p = polys_split.begin();
    for (unsigned int i = 0; i < combinations_valid.size(); ++i) {
      combi combination_valid = *it_combi++;
      PolyhedronSPtr poly_split = *it_p++;
      if (i == selected_combinatorial_split) {
        CGAL_SS3_SPLITTER_TRACE_V(16, "Selected split-combination: " << combi_to_string(combination_valid));
        apply(poly_split, vertex);
        break;
      }
    }
    return true;
  }

  static std::string combi_to_string(const combi& combination)
  {
    std::stringstream sstr;
    sstr << "{ ";
    for (unsigned int i = 0; i < combination.size(); ++i) {
      vec2i split = combination[i];
      if (i > 0) {
        sstr << ", ";
      }
      sstr << "{" << split[0] << ", " << split[1] << "}";
    }
    sstr << " }";
    return sstr.str();
  }

  virtual std::string to_string() const
  {
    std::stringstream sstr;
    sstr << "Combi_vertex_splitter(" << selected_combi_ << ")";
    return sstr.str();
  }

protected:
  unsigned int selected_combi_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_COMBINATORIAL_VERTEX_SPLITTER_H */
