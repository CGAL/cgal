// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Yannick Verdie, Clément Jamin, Pierre Alliez
//

#ifndef CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_OCTREE_H
#define CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_OCTREE_H

#include <CGAL/license/Shape_detection.h>

#include <stack>
#include <limits>

#include <CGAL/Random.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Shape_detection/Efficient_RANSAC/Shape_base.h>

#include <CGAL/Octree.h>

namespace CGAL {
namespace Shape_detection {

namespace internal {

const std::size_t size_t_max = (std::numeric_limits<std::size_t>::max)();

template<class Point_3>
struct Cell {
  std::size_t first, last;
  Cell *child[8];
  Point_3 center;
  std::size_t level;

  Cell(std::size_t first, std::size_t last, Point_3 center, std::size_t level)
          : first(first), last(last), center(center), level(level) {
    memset(child, 0, sizeof(Cell *) * 8);
  }

  bool isLeaf() const {
    for (std::size_t i = 0; i < 8; i++) {
      if (child[i])
        return false;
    }
    return true;
  }

  std::size_t size() const {
    if (first == size_t_max || last == size_t_max)
      return 0;
    else return (last - first + 1);
  }
};

template<class Sdt>
class Direct_octree {

  typedef Sdt Sd_traits;
  typedef typename Sd_traits::Input_range::iterator Input_iterator;
  typedef typename Sd_traits::Point_3 Point_3;
  typedef typename Sd_traits::Vector_3 Vector_3;
  typedef typename Sd_traits::FT FT;
  typedef typename Sd_traits::Point_map Point_map;

public:

  typedef internal::Cell<Point_3> Cell;

public:

  Direct_octree(Sd_traits const &traits)
          : m_traits(traits), m_root(nullptr) {}

  Direct_octree(Sd_traits const &traits,
                const Input_iterator &begin,
                const Input_iterator &beyond,
                Point_map &point_pmap,
                std::size_t offset = 0) :
          m_traits(traits),
          m_root(nullptr),
          m_point_pmap(point_pmap),
          m_first(begin),
          m_offset(offset) {

    m_beyond = (beyond == begin) ? begin : beyond - 1;
  }

  ~Direct_octree() {
    if (!m_root)
      return;

    std::stack<Cell *> stack;
    stack.push(m_root);
    while (!stack.empty()) {
      Cell *cell = stack.top();
      stack.pop();

      for (std::size_t i = 0; i < 8; i++)
        if (cell->child[i])
          stack.push(cell->child[i]);

      delete cell;
    }
  }

  // Sorting data in a way such that points of one cell
  // are always in one range and ordered child-wise:
  // +---+---+
  // | 1.| 0.|
  // +---+---+
  // | 3.| 2.|
  // +---+---+
  // z max before z min, then y max before y min, then x max before x min
  void createTree(double cluster_epsilon_for_max_level_recomputation = -1., std::size_t bucketSize = 2,
                  std::size_t maxLevel = 10) {
    buildBoundingCube();
    std::size_t count = 0;
    m_max_level = 0;

    if (cluster_epsilon_for_max_level_recomputation > 0.) {
      FT bbox_diagonal = (FT) CGAL::sqrt(
              (m_bBox.xmax() - m_bBox.xmin()) * (m_bBox.xmax() - m_bBox.xmin())
              + (m_bBox.ymax() - m_bBox.ymin()) * (m_bBox.ymax() - m_bBox.ymin())
              + (m_bBox.zmax() - m_bBox.zmin()) * (m_bBox.zmax() - m_bBox.zmin()));

      maxLevel = std::size_t(std::log(bbox_diagonal
                                      / cluster_epsilon_for_max_level_recomputation)
                             / std::log(2.0));
    }

    std::stack<Cell *> stack;
    m_root = new Cell(0, this->size() - 1, m_center, 0);
    stack.push(m_root);
    while (!stack.empty()) {
      Cell *cell = stack.top();
      stack.pop();

      m_max_level = std::max<std::size_t>(m_max_level, cell->level);
      if (cell->level == maxLevel)
        continue;

      std::size_t zLowYHighXSplit, zLowYLowXSplit, zLowYSplit;
      std::size_t zHighYSplit, zHighYHighXSplit, zHighYLowXSplit;

      std::size_t zSplit = split(cell->first, cell->last, 2, cell->center.z());

      if (zSplit != size_t_max) {

        zLowYSplit = split(cell->first, zSplit, 1, cell->center.y());
        if (zLowYSplit != size_t_max) {
          zLowYLowXSplit = split(cell->first,
                                 zLowYSplit, 0, cell->center.x());
          zLowYHighXSplit = split(zLowYSplit + 1,
                                  zSplit, 0, cell->center.x());
        } else {
          zLowYLowXSplit = size_t_max;
          zLowYHighXSplit = split(cell->first, zSplit, 0, cell->center.x());
        }

        zHighYSplit = split(zSplit + 1, cell->last, 1, cell->center.y());

        if (zHighYSplit != size_t_max) {
          zHighYHighXSplit = split(zHighYSplit + 1,
                                   cell->last, 0, cell->center.x());
          zHighYLowXSplit = split(zSplit + 1,
                                  zHighYSplit, 0, cell->center.x());
        } else {
          zHighYLowXSplit = size_t_max;
          zHighYHighXSplit = split(zSplit + 1,
                                   cell->last, 0, cell->center.x());
        }
      } else {
        zLowYSplit = size_t_max;
        zLowYLowXSplit = size_t_max;
        zLowYHighXSplit = size_t_max;

        zHighYSplit = split(cell->first,
                            cell->last,
                            1,
                            cell->center.y());

        if (zHighYSplit != size_t_max) {
          zHighYHighXSplit = split(zHighYSplit + 1,
                                   cell->last,
                                   0,
                                   cell->center.x());

          zHighYLowXSplit = split(cell->first,
                                  zHighYSplit,
                                  0,
                                  cell->center.x());
        } else {
          zHighYLowXSplit = size_t_max;
          zHighYHighXSplit = split(cell->first,
                                   cell->last,
                                   0,
                                   cell->center.x());
        }
      }

      FT width = m_width / (1 << (cell->level + 1));

      if (zSplit != size_t_max) {
        if (zLowYSplit != size_t_max) {
          if (zLowYLowXSplit != size_t_max) {

            if (cell->first <= zLowYLowXSplit) {
              //---
              cell->child[7] = new Cell(cell->first,
                                        zLowYLowXSplit,
                                        cell->center + Vector_3(ORIGIN, Point_3(-width, -width, -width)),
                                        cell->level + 1);

              if (cell->child[7]->size() > bucketSize)
                stack.push(cell->child[7]);
            }
          } else zLowYLowXSplit = cell->first - 1;

          if (zLowYLowXSplit < zLowYSplit || zLowYLowXSplit == size_t_max) {
            //+--
            cell->child[6] = new Cell(zLowYLowXSplit + 1,
                                      zLowYSplit,
                                      cell->center + Vector_3(ORIGIN, Point_3(width, -width, -width)),
                                      cell->level + 1);

            if (cell->child[6]->size() > bucketSize)
              stack.push(cell->child[6]);
          }
        } else zLowYSplit = cell->first - 1;

        if (zLowYHighXSplit != size_t_max) {

          if (zLowYSplit < zLowYHighXSplit || zLowYSplit == size_t_max) {
            //-+-
            cell->child[5] = new Cell(zLowYSplit + 1,
                                      zLowYHighXSplit,
                                      cell->center + Vector_3(ORIGIN, Point_3(-width, width, -width)),
                                      cell->level + 1);

            if (cell->child[5]->size() > bucketSize)
              stack.push(cell->child[5]);
          }
        } else zLowYHighXSplit = zLowYSplit;

        if (zLowYHighXSplit < zSplit || zLowYHighXSplit == size_t_max) {
          //++-
          cell->child[4] = new Cell(zLowYHighXSplit + 1,
                                    zSplit,
                                    cell->center + Vector_3(ORIGIN, Point_3(width, width, -width)),
                                    cell->level + 1);

          if (cell->child[4]->size() > bucketSize)
            stack.push(cell->child[4]);
        }
      } else zSplit = cell->first - 1;

      if (zHighYSplit != size_t_max) {
        if (zHighYLowXSplit != size_t_max) {

          if (zSplit < zHighYLowXSplit || zSplit == size_t_max) {
            //--+
            cell->child[3] = new Cell(zSplit + 1,
                                      zHighYLowXSplit,
                                      cell->center + Vector_3(ORIGIN, Point_3(-width, -width, width)),
                                      cell->level + 1);

            if (cell->child[3]->size() > bucketSize)
              stack.push(cell->child[3]);
          }
        } else zHighYLowXSplit = zSplit;

        if (zHighYLowXSplit < zHighYSplit || zHighYLowXSplit == size_t_max) {
          //+-+
          cell->child[2] = new Cell(zHighYLowXSplit + 1,
                                    zHighYSplit,
                                    cell->center + Vector_3(ORIGIN, Point_3(width, -width, width)),
                                    cell->level + 1);

          if (cell->child[2]->size() > bucketSize)
            stack.push(cell->child[2]);
        }

      } else zHighYSplit = zSplit;

      if (zHighYHighXSplit != size_t_max) {
        if (zHighYSplit < zHighYHighXSplit || zHighYSplit == size_t_max) {
          //-++
          cell->child[1] = new Cell(zHighYSplit + 1,
                                    zHighYHighXSplit,
                                    cell->center + Vector_3(ORIGIN, Point_3(-width, width, width)),
                                    cell->level + 1);

          if (cell->child[1]->size() > bucketSize)
            stack.push(cell->child[1]);
        }
      } else zHighYHighXSplit = zHighYSplit;

      if (zHighYHighXSplit <= cell->last || zHighYHighXSplit == size_t_max) {
        if (zHighYHighXSplit < cell->last || zHighYHighXSplit == size_t_max) {
          //+++
          cell->child[0] = new Cell(zHighYHighXSplit + 1,
                                    cell->last,
                                    cell->center + Vector_3(ORIGIN, Point_3(width, width, width)),
                                    cell->level + 1);

          if (cell->child[0]->size() > bucketSize)
            stack.push(cell->child[0]);
        }
      }

      std::size_t sum = 0;
      for (std::size_t i = 0; i < 8; i++)
        sum += (cell->child[i]) ? cell->child[i]->size() : 0;

      count++;
    }
  }

  std::size_t maxLevel() const { return m_max_level; }

  const Bbox_3 &boundingBox() const { return m_bBox; }

  const Cell *root() const { return m_root; }

  FT width() const { return m_width; }

  Input_iterator at(std::size_t i) {
    return m_first + i;
  }

  std::size_t index(std::size_t i) const {
    return i + m_offset;
  }

  std::size_t offset() {
    return m_offset;
  }

  std::size_t size() {
    return m_beyond - m_first + 1;
  }

  void swap(std::size_t a, std::size_t b) {
    typename std::iterator_traits<Input_iterator>::value_type tmp;
    tmp = m_first[a];
    m_first[a] = m_first[b];
    m_first[b] = tmp;
  }

private:

  Sd_traits m_traits;
  Bbox_3 m_bBox;
  Cell *m_root;
  Point_3 m_center;
  FT m_width;
  std::size_t m_max_level;
  Point_map m_point_pmap;

  Input_iterator m_first;
  Input_iterator m_beyond;
  std::size_t m_offset;

private:

  // returns index of last point below threshold
  std::size_t split(std::size_t first, std::size_t last, std::size_t dimension, FT threshold) {
    if (last == size_t_max || first == size_t_max)
      return size_t_max;

    if (first > last)
      return first - 1;

    std::size_t origFirst = first;

    while (first < last) {
      // find first above threshold
      while (get(m_point_pmap, *this->at(first))[dimension] < threshold
             && first < last) {
        first++;
      }

      // check if last has been reached
      if (first == last) {
        return (get(m_point_pmap, *this->at(first))[dimension] < threshold) ?
               first : (first == origFirst) ? size_t_max : first - 1;
      }

      // find last below threshold
      while (get(m_point_pmap, *this->at(last))[dimension] >= threshold
             && last > first) {
        last--;
      }

      // check if first has been reached
      if (last == first) {
        return (get(m_point_pmap, *this->at(first))[dimension] < threshold) ?
               first : (first == origFirst) ? size_t_max : first - 1;
      }

      this->swap(first, last);
      first++;
      last--;
    }

    return (get(m_point_pmap, *this->at(first))[dimension] < threshold) ?
           first : (first == origFirst) ? size_t_max : first - 1;
  }

  const Bbox_3 &buildBoundingCube() {
    FT min[] = {std::numeric_limits<FT>::infinity(),
                std::numeric_limits<FT>::infinity(),
                std::numeric_limits<FT>::infinity()};
    FT max[] = {-std::numeric_limits<FT>::infinity(),
                -std::numeric_limits<FT>::infinity(),
                -std::numeric_limits<FT>::infinity()};

    for (std::size_t i = 0; i < this->size(); i++) {
      Point_3 p = get(m_point_pmap, *this->at(i));
      min[0] = (std::min<FT>)(min[0], p.x());
      min[1] = (std::min<FT>)(min[1], p.y());
      min[2] = (std::min<FT>)(min[2], p.z());
      max[0] = (std::max<FT>)(max[0], p.x());
      max[1] = (std::max<FT>)(max[1], p.y());
      max[2] = (std::max<FT>)(max[2], p.z());
    }

    m_bBox = Bbox_3(min[0], min[1], min[2], max[0], max[1], max[2]);

    m_width = (std::max)(max[0] - min[0],
                         (std::max)(max[1] - min[1], max[2] - min[2])) * (FT) 0.5;
    m_center = Point_3((min[0] + max[0]) * (FT) 0.5,
                       (min[1] + max[1]) * (FT) 0.5,
                       (min[2] + max[2]) * (FT) 0.5);

    return m_bBox;
  }

};

template<class Sdt>
class Indexed_octree {


  typedef Sdt Sd_traits;
  typedef typename Sd_traits::Input_range::iterator Input_iterator;
  typedef typename Sd_traits::Point_3 Point_3;
  typedef typename Sd_traits::Vector_3 Vector_3;
  typedef typename Sd_traits::FT FT;
  typedef typename Sd_traits::Point_map Point_map;

public:

  typedef internal::Cell<Point_3> Cell;

public:

  Indexed_octree(Sd_traits const &traits)
          : m_traits(traits), m_root(nullptr) {}

  Indexed_octree(Sd_traits
                 const &traits,
                 const Input_iterator &begin,
                 const Input_iterator &beyond,
                 Point_map
                 &point_pmap,
                 std::size_t offset = 0) :
          m_traits(traits),
          m_root(nullptr),
          m_point_pmap(point_pmap),
          m_first(begin) {

    m_beyond = (beyond == begin) ? begin : beyond - 1;
    m_indices.resize(size());
    for (std::size_t i = 0; i < size(); i++)
      m_indices[i] = i;
  }

  ~Indexed_octree() {
    if (!m_root)
      return;

    std::stack<Cell *> stack;
    stack.push(m_root);
    while (!stack.empty()) {
      Cell *cell = stack.top();
      stack.pop();

      for (std::size_t i = 0; i < 8; i++)
        if (cell->child[i])
          stack.push(cell->child[i]);

      delete cell;
    }
  }

  // Sorting data in a way such that points of one cell
  // are always in one range and ordered child-wise:
  // +---+---+
  // | 1.| 0.|
  // +---+---+
  // | 3.| 2.|
  // +---+---+
  // z max before z min, then y max before y min, then x max before x min
  void createTree(double cluster_epsilon_for_max_level_recomputation = -1., std::size_t bucketSize = 2,
                  std::size_t maxLevel = 10) {
    buildBoundingCube();
    std::size_t count = 0;
    m_max_level = 0;

    if (cluster_epsilon_for_max_level_recomputation > 0.) {
      FT bbox_diagonal = (FT) CGAL::sqrt(
              (m_bBox.xmax() - m_bBox.xmin()) * (m_bBox.xmax() - m_bBox.xmin())
              + (m_bBox.ymax() - m_bBox.ymin()) * (m_bBox.ymax() - m_bBox.ymin())
              + (m_bBox.zmax() - m_bBox.zmin()) * (m_bBox.zmax() - m_bBox.zmin()));

      maxLevel = std::size_t(std::log(bbox_diagonal
                                      / cluster_epsilon_for_max_level_recomputation)
                             / std::log(2.0));
    }

    std::stack<Cell *> stack;
    m_root = new Cell(0, this->size() - 1, m_center, 0);
    stack.push(m_root);
    while (!stack.empty()) {
      Cell *cell = stack.top();
      stack.pop();

      m_max_level = std::max<std::size_t>(m_max_level, cell->level);
      if (cell->level == maxLevel)
        continue;

      std::size_t zLowYHighXSplit, zLowYLowXSplit, zLowYSplit;
      std::size_t zHighYSplit, zHighYHighXSplit, zHighYLowXSplit;

      std::size_t zSplit = split(cell->first, cell->last, 2, cell->center.z());

      if (zSplit != size_t_max) {

        zLowYSplit = split(cell->first, zSplit, 1, cell->center.y());
        if (zLowYSplit != size_t_max) {
          zLowYLowXSplit = split(cell->first,
                                 zLowYSplit, 0, cell->center.x());
          zLowYHighXSplit = split(zLowYSplit + 1,
                                  zSplit, 0, cell->center.x());
        } else {
          zLowYLowXSplit = size_t_max;
          zLowYHighXSplit = split(cell->first, zSplit, 0, cell->center.x());
        }

        zHighYSplit = split(zSplit + 1, cell->last, 1, cell->center.y());

        if (zHighYSplit != size_t_max) {
          zHighYHighXSplit = split(zHighYSplit + 1,
                                   cell->last, 0, cell->center.x());
          zHighYLowXSplit = split(zSplit + 1,
                                  zHighYSplit, 0, cell->center.x());
        } else {
          zHighYLowXSplit = size_t_max;
          zHighYHighXSplit = split(zSplit + 1,
                                   cell->last, 0, cell->center.x());
        }
      } else {
        zLowYSplit = size_t_max;
        zLowYLowXSplit = size_t_max;
        zLowYHighXSplit = size_t_max;

        zHighYSplit = split(cell->first,
                            cell->last,
                            1,
                            cell->center.y());

        if (zHighYSplit != size_t_max) {
          zHighYHighXSplit = split(zHighYSplit + 1,
                                   cell->last,
                                   0,
                                   cell->center.x());

          zHighYLowXSplit = split(cell->first,
                                  zHighYSplit,
                                  0,
                                  cell->center.x());
        } else {
          zHighYLowXSplit = size_t_max;
          zHighYHighXSplit = split(cell->first,
                                   cell->last,
                                   0,
                                   cell->center.x());
        }
      }

      FT width = m_width / (1 << (cell->level + 1));

      if (zSplit != size_t_max) {
        if (zLowYSplit != size_t_max) {
          if (zLowYLowXSplit != size_t_max) {

            if (cell->first <= zLowYLowXSplit) {
              //---
              cell->child[7] = new Cell(cell->first,
                                        zLowYLowXSplit,
                                        cell->center + Vector_3(ORIGIN, Point_3(-width, -width, -width)),
                                        cell->level + 1);

              if (cell->child[7]->size() > bucketSize)
                stack.push(cell->child[7]);
            }
          } else zLowYLowXSplit = cell->first - 1;

          if (zLowYLowXSplit < zLowYSplit || zLowYLowXSplit == size_t_max) {
            //+--
            cell->child[6] = new Cell(zLowYLowXSplit + 1,
                                      zLowYSplit,
                                      cell->center + Vector_3(ORIGIN, Point_3(width, -width, -width)),
                                      cell->level + 1);

            if (cell->child[6]->size() > bucketSize)
              stack.push(cell->child[6]);
          }
        } else zLowYSplit = cell->first - 1;

        if (zLowYHighXSplit != size_t_max) {

          if (zLowYSplit < zLowYHighXSplit || zLowYSplit == size_t_max) {
            //-+-
            cell->child[5] = new Cell(zLowYSplit + 1,
                                      zLowYHighXSplit,
                                      cell->center + Vector_3(ORIGIN, Point_3(-width, width, -width)),
                                      cell->level + 1);

            if (cell->child[5]->size() > bucketSize)
              stack.push(cell->child[5]);
          }
        } else zLowYHighXSplit = zLowYSplit;

        if (zLowYHighXSplit < zSplit || zLowYHighXSplit == size_t_max) {
          //++-
          cell->child[4] = new Cell(zLowYHighXSplit + 1,
                                    zSplit,
                                    cell->center + Vector_3(ORIGIN, Point_3(width, width, -width)),
                                    cell->level + 1);

          if (cell->child[4]->size() > bucketSize)
            stack.push(cell->child[4]);
        }
      } else zSplit = cell->first - 1;

      if (zHighYSplit != size_t_max) {
        if (zHighYLowXSplit != size_t_max) {

          if (zSplit < zHighYLowXSplit || zSplit == size_t_max) {
            //--+
            cell->child[3] = new Cell(zSplit + 1,
                                      zHighYLowXSplit,
                                      cell->center + Vector_3(ORIGIN, Point_3(-width, -width, width)),
                                      cell->level + 1);

            if (cell->child[3]->size() > bucketSize)
              stack.push(cell->child[3]);
          }
        } else zHighYLowXSplit = zSplit;

        if (zHighYLowXSplit < zHighYSplit || zHighYLowXSplit == size_t_max) {
          //+-+
          cell->child[2] = new Cell(zHighYLowXSplit + 1,
                                    zHighYSplit,
                                    cell->center + Vector_3(ORIGIN, Point_3(width, -width, width)),
                                    cell->level + 1);

          if (cell->child[2]->size() > bucketSize)
            stack.push(cell->child[2]);
        }

      } else zHighYSplit = zSplit;

      if (zHighYHighXSplit != size_t_max) {
        if (zHighYSplit < zHighYHighXSplit || zHighYSplit == size_t_max) {
          //-++
          cell->child[1] = new Cell(zHighYSplit + 1,
                                    zHighYHighXSplit,
                                    cell->center + Vector_3(ORIGIN, Point_3(-width, width, width)),
                                    cell->level + 1);

          if (cell->child[1]->size() > bucketSize)
            stack.push(cell->child[1]);
        }
      } else zHighYHighXSplit = zHighYSplit;

      if (zHighYHighXSplit <= cell->last || zHighYHighXSplit == size_t_max) {
        if (zHighYHighXSplit < cell->last || zHighYHighXSplit == size_t_max) {
          //+++
          cell->child[0] = new Cell(zHighYHighXSplit + 1,
                                    cell->last,
                                    cell->center + Vector_3(ORIGIN, Point_3(width, width, width)),
                                    cell->level + 1);

          if (cell->child[0]->size() > bucketSize)
            stack.push(cell->child[0]);
        }
      }

      std::size_t sum = 0;
      for (std::size_t i = 0; i < 8; i++)
        sum += (cell->child[i]) ? cell->child[i]->size() : 0;

      count++;
    }
  }

  std::size_t maxLevel() const { return m_max_level; }

  const Bbox_3 &boundingBox() const { return m_bBox; }

  const Cell *root() const { return m_root; }

  FT width() const { return m_width; }


  Input_iterator at(std::size_t i) {
    return m_first + m_indices[i];
  }

  std::size_t index(std::size_t i) const {
    return m_indices[i];
  }

  std::size_t offset() {
    return 0;
  }

  std::size_t size() {
    return m_beyond - m_first + 1;
  }

  void swap(std::size_t a, std::size_t b) {
    std::size_t tmp = m_indices[a];
    m_indices[a] = m_indices[b];
    m_indices[b] = tmp;
  }

private:

  Sd_traits m_traits;
  Bbox_3 m_bBox;
  Cell *m_root;
  Point_3 m_center;
  FT m_width;
  std::size_t m_max_level;
  Point_map m_point_pmap;

  std::vector<std::size_t> m_indices;
  Input_iterator m_first;
  Input_iterator m_beyond;

private:

  // returns index of last point below threshold
  std::size_t split(std::size_t first, std::size_t last, std::size_t dimension, FT threshold) {
    if (last == size_t_max || first == size_t_max)
      return size_t_max;

    if (first > last)
      return first - 1;

    std::size_t origFirst = first;

    while (first < last) {
      // find first above threshold
      while (get(m_point_pmap, *this->at(first))[dimension] < threshold
             && first < last) {
        first++;
      }

      // check if last has been reached
      if (first == last) {
        return (get(m_point_pmap, *this->at(first))[dimension] < threshold) ?
               first : (first == origFirst) ? size_t_max : first - 1;
      }

      // find last below threshold
      while (get(m_point_pmap, *this->at(last))[dimension] >= threshold
             && last > first) {
        last--;
      }

      // check if first has been reached
      if (last == first) {
        return (get(m_point_pmap, *this->at(first))[dimension] < threshold) ?
               first : (first == origFirst) ? size_t_max : first - 1;
      }

      this->swap(first, last);
      first++;
      last--;
    }

    return (get(m_point_pmap, *this->at(first))[dimension] < threshold) ?
           first : (first == origFirst) ? size_t_max : first - 1;
  }

  const Bbox_3 &buildBoundingCube() {
    FT min[] = {std::numeric_limits<FT>::infinity(),
                std::numeric_limits<FT>::infinity(),
                std::numeric_limits<FT>::infinity()};
    FT max[] = {-std::numeric_limits<FT>::infinity(),
                -std::numeric_limits<FT>::infinity(),
                -std::numeric_limits<FT>::infinity()};

    for (std::size_t i = 0; i < this->size(); i++) {
      Point_3 p = get(m_point_pmap, *this->at(i));
      min[0] = (std::min<FT>)(min[0], p.x());
      min[1] = (std::min<FT>)(min[1], p.y());
      min[2] = (std::min<FT>)(min[2], p.z());
      max[0] = (std::max<FT>)(max[0], p.x());
      max[1] = (std::max<FT>)(max[1], p.y());
      max[2] = (std::max<FT>)(max[2], p.z());
    }

    m_bBox = Bbox_3(min[0], min[1], min[2], max[0], max[1], max[2]);

    m_width = (std::max)(max[0] - min[0],
                         (std::max)(max[1] - min[1], max[2] - min[2])) * (FT) 0.5;
    m_center = Point_3((min[0] + max[0]) * (FT) 0.5,
                       (min[1] + max[1]) * (FT) 0.5,
                       (min[2] + max[2]) * (FT) 0.5);

    return m_bBox;
  }

};

}
}
}

#endif
