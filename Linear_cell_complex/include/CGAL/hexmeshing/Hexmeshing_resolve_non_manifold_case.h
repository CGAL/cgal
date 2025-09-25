// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Soichiro Yamazaki <soichiro19998@gmail.com>
//
#ifndef HEXMESHING_RESOLVE_NON_MANIFOLD_CASE_H
#define HEXMESHING_RESOLVE_NON_MANIFOLD_CASE_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_two_refinement_utils.h>
#include <CGAL/hexmeshing/Hexmeshing_set_attributes.h>
#include <array>
#include <vector>
#include <queue>
#include <cassert>


namespace CGAL::internal::Hexmeshing {
  // all volume cells need to have attribute
  int __get_signal(LCC& lcc, Dart_handle vertex, size_type inner_mark) {
    std::array<Dart_handle, 8> volumes = volumes_around_node(lcc, vertex);
    int inner_signal = 0;
    for(int i = 0; i < 8; i++) {
      auto it = volumes[i];
      if(it == nullptr or lcc.attribute<3>(it) == nullptr) {
        return -1;
      }
      inner_signal <<= 1;
      inner_signal |= int(lcc.is_marked(it, inner_mark));
    }
    return inner_signal;
  }

  std::array<std::vector<int>, 7> __get_seven_non_manifold_templates() {
    return {{{0b10000010, 0b01000001, 0b00101000, 0b00010100, 0b01111101, 0b11101011, 0b11010111, 0b10111110},
            {0b10000100, 0b01000010, 0b00100001, 0b00011000, 0b10000001, 0b01001000, 0b00100100, 0b00010010, 0b10100000, 0b01010000, 0b00001010, 0b00000101, 0b01111011, 0b01111110, 0b01011111, 0b11100111, 0b10101111, 0b10110111, 0b11011110, 0b11011011, 0b11111010, 0b10111101, 0b11110101, 0b11101101},
            {0b00010110, 0b01000011, 0b00101010, 0b00101100, 0b00010101, 0b01100001, 0b01001001, 0b00111000, 0b01010100, 0b10000011, 0b10100010, 0b10000110, 0b01101000, 0b00011100, 0b01010001, 0b11000001, 0b10101000, 0b10010100, 0b10010010, 0b11000010, 0b10001010, 0b00110100, 0b01000101, 0b00101001, 0b11101001, 0b10111100, 0b11010101, 0b11010011, 0b11101010, 0b10011110, 0b10110110, 0b11000111, 0b10101011, 0b01111100, 0b01011101, 0b01111001, 0b10010111, 0b11100011, 0b10101110, 0b00111110, 0b01010111, 0b01101011, 0b01101101, 0b00111101, 0b01110101, 0b11001011, 0b10111010, 0b11010110},
            {0b01001010, 0b00100101, 0b00011010, 0b10000101, 0b10100100, 0b01010010, 0b10100001, 0b01011000, 0b10110101, 0b01111010, 0b11100101, 0b11011010, 0b10101101, 0b01011011, 0b10100111, 0b01011110},
            {0b01010101, 0b10101010, 0b10010110, 0b01101001, 0b11000011, 0b00111100},
            {0b00011110, 0b01010011, 0b01101010, 0b00101101, 0b00110101, 0b01100101, 0b01001011, 0b00111010, 0b01010110, 0b10000111, 0b10100011, 0b10100110, 0b01111000, 0b01011100, 0b01011001, 0b11100001, 0b10101100, 0b10010101, 0b11010010, 0b11001010, 0b10011010, 0b10110100, 0b11000101, 0b10101001},
            {0b01011010, 0b10100101}}};
  }

  std::array<int, 256> __get_non_manifold_template_list() {
    return {{0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,
            0,0,1,0,1,1,1,0,1,0,1,0,1,0,1,0,
            0,1,0,0,1,1,0,0,1,1,1,0,1,1,0,0,
            0,0,0,0,1,1,0,0,1,0,1,0,1,1,1,0,
            0,1,1,1,0,1,0,0,1,1,1,1,0,0,0,0,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            0,1,0,0,0,1,0,0,1,1,1,1,0,1,0,0,
            0,0,0,0,0,1,0,0,1,1,1,1,1,1,1,0,
            0,1,1,1,1,1,1,1,0,0,1,0,0,0,0,0,
            0,0,1,0,1,1,1,1,0,0,1,0,0,0,1,0,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            0,0,0,0,1,1,1,1,0,0,1,0,1,1,1,0,
            0,1,1,1,0,1,0,1,0,0,1,1,0,0,0,0,
            0,0,1,1,0,1,1,1,0,0,1,1,0,0,1,0,
            0,1,0,1,0,1,0,1,0,1,1,1,0,1,0,0,
            0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0}};
  }

  std::array<std::vector<int>, 256> __get_solution_to_non_manifold_templates_list() {
    std::array<std::vector<int>, 7> seven_templates = __get_seven_non_manifold_templates();
    std::array<int, 256> is_in_templates = __get_non_manifold_template_list();

    std::array<std::vector<int>, 256> resolve_templates = {};
    for(int i = 0; i < 7; i++) {
      const std::vector<int> &templates = seven_templates[i];
      int template_size = templates.size();
      for(auto temp: templates) {
        for(int signal_num = 1; signal_num <= 8; signal_num++) {
          for(int signal = (1<<signal_num)-1; signal < (1<<8);) {
            int after = temp ^ signal;

            if(is_in_templates[after]^1) resolve_templates[temp].emplace_back(signal);

            int x = signal & -signal, y = signal + x;
            signal = (((signal & ~y) / x) >> 1) | y;
          }
        }
      }
    }

    return resolve_templates;
  }

  // Algorithm 1 in Owen et al. (2014)
  // prerequisite: If the vertex is non-manifold, the manifold condition must be satisfiable through a modification of the fraction within 0.37.
  void resolve_non_manifold_case(LCC& lcc, double s, size_type inner_mark) {
    std::array<int, 256> is_in_templates = __get_non_manifold_template_list();

    std::array<std::vector<int>, 256> resolve_templates = __get_solution_to_non_manifold_templates_list();

    auto vertices = lcc.one_dart_per_cell<0>();
    std::queue<std::pair<Dart_handle, int>> que;
    for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {
      int inner_signal = __get_signal(lcc, vertex, inner_mark);
      if(inner_signal == -1) continue;
      if(is_in_templates[inner_signal]) que.emplace(vertex, inner_signal);
    }

    const int NUMBER_OF_EPSILONS = 15;
    std::array<double, NUMBER_OF_EPSILONS> epsilons = {0.05, 0.07, 0.11, 0.13, 0.17, 0.19, 0.23, 0.29, 0.31, 0.37, 0.41, 0.43, 0.47, 0.53, 0.59};
    std::vector<int> counter(set_vertex_ids(lcc));
    while(!que.empty()) {
      auto [vertex, temp] = que.front();
      que.pop();
      int inner_signal = __get_signal(lcc, vertex, inner_mark);
      if(temp != inner_signal) continue;

      int id = lcc.attribute<0>(vertex)->id;
      assert(counter[id] < NUMBER_OF_EPSILONS);
      double eps = epsilons[counter[id]++];
      std::array<Dart_handle, 8> volumes = volumes_around_node(lcc, vertex);

      int able_to_change = 0;
      for(auto volume: volumes) {
        able_to_change <<= 1;
        able_to_change |= int(abs(lcc.attribute<3>(volume)->info().fraction - s) < eps);
      }

      auto& solutions = resolve_templates[temp];
      int solution = -1;
      for(auto candidate: solutions) {
        if((candidate|able_to_change) != able_to_change) continue;
        solution = candidate;
        break;
      }

      if(solution == -1) {
        que.emplace(vertex, temp);
      }
      else {
        int L_sub = solution & inner_signal;
        int L_add = solution ^ L_sub;
        for(int i = 0; i < 8; i++) {
          if(L_sub>>(7-i)&1) {
            lcc.attribute<3>(volumes[i])->info().fraction -= eps;
            lcc.unmark_cell<3>(volumes[i], inner_mark);
          }
          if(L_add>>(7-i)&1) {
            lcc.attribute<3>(volumes[i])->info().fraction += eps;
            lcc.mark_cell<3>(volumes[i], inner_mark);
          }
        }
        for(int i = 0; i < 8; i++) {
          if(solution>>(7-i)&1) {
            auto vertices = lcc.one_dart_per_incident_cell<0, 3>(volumes[i]);
            for(auto it = vertices.begin(); it != vertices.end(); it++) {
              if(lcc.attribute<0>(vertex)->id == lcc.attribute<0>(it)->id) continue;
              int now_signal = __get_signal(lcc, it, inner_mark);
              if(is_in_templates[now_signal]) que.emplace(it, now_signal);
            }
          }
        }
      }
    }
  }
}




#endif