// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
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
//
//
// Author(s)     : Sven Oesau, Yannick Verdie, Clément Jamin, Pierre Alliez
//

#ifndef CGAL_SHAPE_DETECTION_3_OCTREE_H
#define CGAL_SHAPE_DETECTION_3_OCTREE_H

#include <CGAL/license/Point_set_shape_detection_3.h>


#include <limits>
#include <stack>

#include <CGAL/Random.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Shape_detection_3/Shape_base.h>


extern int scoreTime;

namespace CGAL {
  namespace Shape_detection_3 {
  
    template<class Traits> 
    class Efficient_RANSAC;
    
  namespace internal {

    const std::size_t size_t_max = (std::numeric_limits<std::size_t>::max)();
  
    template<class Sdt>
    class DirectPointAccessor {
    public:
      typedef Sdt Sd_traits;
      typedef typename Sd_traits::Input_range::iterator Input_iterator;

      DirectPointAccessor() {}
      DirectPointAccessor(const Input_iterator &begin,
                          const Input_iterator &beyond, 
                          std::size_t offset) : m_first(begin), m_offset(offset) {
        m_beyond = (beyond == begin) ? begin : beyond - 1;
      }

      Input_iterator at(std::size_t i) {
        return m_first + i;
      }

      std::size_t index(std::size_t i) {
        return i + m_offset;
      }

      std::size_t offset() {
        return m_offset;
      }

      std::size_t size() {
        return m_beyond - m_first + 1;
      }

      Input_iterator first() {
        return m_first;
      }

      Input_iterator beyond() {
        return m_beyond;
      }

      void setData(Input_iterator &begin, Input_iterator &beyond) {
        m_beyond = (beyond == begin) ? begin : beyond - 1;
      }

      void swap(std::size_t a, std::size_t b) {
        typename std::iterator_traits<Input_iterator>::value_type tmp;
        tmp = m_first[a];
        m_first[a] = m_first[b];
        m_first[b] = tmp;
      }

    protected:
      Input_iterator m_first;

    private:
      Input_iterator m_beyond;
      std::size_t m_offset;
    };

    template<class Sdt>
    class IndexedPointAccessor {
    public:
      typedef Sdt Sd_traits;
      typedef typename Sd_traits::Input_range::iterator Input_iterator;

      IndexedPointAccessor() {}
      IndexedPointAccessor(const Input_iterator &begin,
                           const Input_iterator &beyond, std::size_t)
                           : m_first(begin) {
        m_beyond = (beyond == begin) ? begin : beyond - 1;
        m_indices.resize(size());
        for (std::size_t i = 0;i<size();i++)
          m_indices[i] = i;
      }

      Input_iterator at(std::size_t i) {
        return m_first + m_indices[i];
      }

      std::size_t index(std::size_t i) {
        return m_indices[i];
      }

      std::size_t offset() {
        return 0;
      }

      Input_iterator first() {
        return m_first;
      }

      Input_iterator beyond() {
        return m_beyond;
      }

      void setData(Input_iterator &begin, Input_iterator &beyond) {
        m_beyond = (beyond == begin) ? begin : beyond - 1;
        m_indices.resize(size());
        for (std::size_t i = 0;i<size();i++)
          m_indices[i] = i;
      }

      std::size_t size() {
        return m_beyond - m_first + 1;
      }

      void swap(std::size_t a, std::size_t b) {
        std::size_t tmp = m_indices[a];
        m_indices[a] = m_indices[b];
        m_indices[b] = tmp;
      }

    protected:
      Input_iterator m_first;

    private:
      std::vector<std::size_t> m_indices;
      Input_iterator m_beyond;
    };

    template<class PointAccessor>
    class Octree : public PointAccessor {

      typedef typename PointAccessor::Sd_traits Sd_traits;
      typedef typename Sd_traits::Input_range::iterator Input_iterator;
      typedef Shape_base<Sd_traits> Shape;
      typedef typename Sd_traits::Point_3 Point_3;
      typedef typename Sd_traits::Vector_3 Vector_3;
      typedef typename Sd_traits::FT FT;
      typedef typename Sd_traits::Point_map Point_map;
      typedef typename Sd_traits::Normal_map Normal_map;

      template<class Sd_traits>
        friend class ::CGAL::Shape_detection_3::Efficient_RANSAC;

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
          for (std::size_t i = 0;i<8;i++) {
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

      // --------------------------------------------------------------------------
      // Utilities
      // --------------------------------------------------------------------------
      FT get_x(const Vector_3& v){ return m_traits.compute_x_3_object()(v); }
      FT get_y(const Vector_3& v){ return m_traits.compute_y_3_object()(v); }
      FT get_z(const Vector_3& v){ return m_traits.compute_z_3_object()(v); }
      FT get_x(const Point_3& p){ return m_traits.compute_x_3_object()(p); }
      FT get_y(const Point_3& p){ return m_traits.compute_y_3_object()(p); }
      FT get_z(const Point_3& p){ return m_traits.compute_z_3_object()(p); }
      FT get_coord(const Point_3& p, unsigned int d)
      {
        CGAL_assertion(d >= 0 && d < 3);
        switch (d)
        {
          case 0: return get_x(p);
          case 1: return get_y(p);
          case 2: return get_z(p);
          default: return FT(0);
        }
      }

      Point_3 constr_pt(FT x, FT y, FT z) const
      { return m_traits.construct_point_3_object()(x, y, z); }
      Vector_3 constr_vec(const Point_3& p, const Point_3& q) const
      { return m_traits.construct_vector_3_object()(p, q); }

      Point_3 transl(const Point_3& p, const Vector_3 &v)
      { return m_traits.construct_translated_point_3_object()(p, v); }
        
    public:
      Octree(Sd_traits const& traits) 
        : m_traits(traits), m_bucket_size(20), m_set_max_level(10), m_root(NULL) {}
      Octree(Sd_traits const& traits,
             const Input_iterator &first,
             const Input_iterator &beyond,
             Point_map& point_pmap,
             Normal_map& normal_pmap,
             std::size_t offset = 0,
             std::size_t bucketSize = 20, 
             std::size_t maxLevel = 10)
             : PointAccessor(first, beyond, offset),
               m_traits(traits),
               m_root(NULL),
               m_bucket_size(bucketSize),
               m_set_max_level(maxLevel),
               m_point_pmap (point_pmap),
               m_normal_pmap (normal_pmap) {}

      ~Octree() {
        if (!m_root)
          return;

        std::stack<Cell *> stack;
        stack.push(m_root);
        while (!stack.empty()) {
          Cell *cell = stack.top();
          stack.pop();

          for (std::size_t i = 0;i<8;i++)
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
      void createTree() {
        buildBoundingCube();
        std::size_t count = 0;
        m_max_level = 0;

        std::stack<Cell *> stack;
        m_root = new Cell(0, this->size() - 1, m_center, 0);
        stack.push(m_root);
        while (!stack.empty()) {
          Cell *cell= stack.top();
          stack.pop();

          m_max_level = std::max<std::size_t>(m_max_level, cell->level);
          if (cell->level == m_set_max_level)
            continue;
          
          std::size_t zLowYHighXSplit, zLowYLowXSplit, zLowYSplit;
          std::size_t zHighYSplit, zHighYHighXSplit, zHighYLowXSplit;

          std::size_t zSplit = split(cell->first, cell->last, 2, get_z(cell->center));

          if (zSplit != size_t_max) {

            zLowYSplit = split(cell->first, zSplit, 1, get_y(cell->center));
            if (zLowYSplit != size_t_max) {
              zLowYLowXSplit = split(cell->first,
                                     zLowYSplit, 0, get_x(cell->center));
              zLowYHighXSplit = split(zLowYSplit + 1,
                                      zSplit, 0, get_x(cell->center)); 
            }
            else {
              zLowYLowXSplit = size_t_max;
              zLowYHighXSplit = split(cell->first, zSplit, 0, get_x(cell->center));
            }

            zHighYSplit = split(zSplit + 1, cell->last, 1, get_y(cell->center));  

            if (zHighYSplit != size_t_max) {
              zHighYHighXSplit = split(zHighYSplit + 1,
                                       cell->last, 0, get_x(cell->center));
              zHighYLowXSplit = split(zSplit + 1,
                                      zHighYSplit, 0, get_x(cell->center));
            }
            else {
              zHighYLowXSplit = size_t_max;
              zHighYHighXSplit = split(zSplit + 1,
                                       cell->last, 0, get_x(cell->center));
            }
          }
          else {
            zLowYSplit = size_t_max;
            zLowYLowXSplit = size_t_max;
            zLowYHighXSplit = size_t_max;

            zHighYSplit = split(cell->first,
                                cell->last,
                                1,
                                get_y(cell->center)); 

            if (zHighYSplit != size_t_max) {
              zHighYHighXSplit = split(zHighYSplit + 1, 
                                       cell->last, 
                                       0, 
                                       get_x(cell->center));

              zHighYLowXSplit = split(cell->first, 
                                      zHighYSplit,
                                      0, 
                                      get_x(cell->center)); 
            }
            else {
              zHighYLowXSplit = size_t_max;
              zHighYHighXSplit = split(cell->first,
                                       cell->last, 
                                       0,
                                       get_x(cell->center));
            }
          }


          FT width = m_width / (1<<(cell->level + 1));

          if (zSplit != size_t_max) {
            if (zLowYSplit != size_t_max) {
              if (zLowYLowXSplit != size_t_max) {

                if (cell->first <= zLowYLowXSplit) {
                  //---
                  cell->child[7] = new Cell(cell->first,
                                            zLowYLowXSplit, 
                                            transl(cell->center, constr_vec(
                                              ORIGIN, constr_pt(-width,-width,-width))),
                                            cell->level + 1);

                  if (cell->child[7]->size() > m_bucket_size)
                    stack.push(cell->child[7]);
                }
              }
              else zLowYLowXSplit = cell->first - 1;

              if (zLowYLowXSplit < zLowYSplit || zLowYLowXSplit == size_t_max) {
                //+--
                cell->child[6] = new Cell(zLowYLowXSplit + 1,
                                          zLowYSplit, 
                                          transl(cell->center, constr_vec(
                                            ORIGIN, constr_pt(width,-width,-width))),
                                          cell->level + 1);

                if (cell->child[6]->size() > m_bucket_size)
                  stack.push(cell->child[6]);
              }
            }
            else zLowYSplit = cell->first - 1;

            if (zLowYHighXSplit != size_t_max) {

              if (zLowYSplit < zLowYHighXSplit || zLowYSplit == size_t_max) {
                //-+-
                cell->child[5] = new Cell(zLowYSplit + 1, 
                                          zLowYHighXSplit, 
                                          transl(cell->center, constr_vec(
                                            ORIGIN, constr_pt(-width, width,-width))),
                                          cell->level + 1);

                if (cell->child[5]->size() > m_bucket_size)
                  stack.push(cell->child[5]);
              }
            }
            else zLowYHighXSplit = zLowYSplit;

            if (zLowYHighXSplit < zSplit || zLowYHighXSplit == size_t_max) {
              //++-
              cell->child[4] = new Cell(zLowYHighXSplit + 1,
                                        zSplit, 
                                        transl(cell->center, constr_vec(
                                          ORIGIN, constr_pt(width, width,-width))),
                                        cell->level + 1);

              if (cell->child[4]->size() > m_bucket_size)
                stack.push(cell->child[4]);
            }
          }
          else zSplit = cell->first - 1;

          if (zHighYSplit != size_t_max) {
            if (zHighYLowXSplit != size_t_max) {

              if (zSplit < zHighYLowXSplit || zSplit == size_t_max) {
                //--+
                cell->child[3] = new Cell(zSplit + 1,
                                          zHighYLowXSplit,
                                          transl(cell->center, constr_vec(
                                            ORIGIN, constr_pt(-width,-width, width))),
                                          cell->level + 1);

                if (cell->child[3]->size() > m_bucket_size)
                  stack.push(cell->child[3]);
              }
            }
            else zHighYLowXSplit = zSplit;

            if (zHighYLowXSplit < zHighYSplit || zHighYLowXSplit == size_t_max) {
              //+-+
              cell->child[2] = new Cell(zHighYLowXSplit + 1,
                                        zHighYSplit, 
                                        transl(cell->center, constr_vec(
                                          ORIGIN, constr_pt(width,-width, width))),
                                        cell->level + 1);

              if (cell->child[2]->size() > m_bucket_size)
                stack.push(cell->child[2]);
            }

          }
          else zHighYSplit = zSplit;

          if (zHighYHighXSplit != size_t_max) {
            if (zHighYSplit < zHighYHighXSplit || zHighYSplit == size_t_max) {
              //-++
              cell->child[1] = new Cell(zHighYSplit + 1,
                                        zHighYHighXSplit, 
                                        transl(cell->center, constr_vec(
                                          ORIGIN, constr_pt(-width, width, width))), 
                                        cell->level + 1);

              if (cell->child[1]->size() > m_bucket_size)
                stack.push(cell->child[1]);
            }
          }
          else zHighYHighXSplit = zHighYSplit;

          if (zHighYHighXSplit <= cell->last || zHighYHighXSplit == size_t_max) {
            if (zHighYHighXSplit < cell->last || zHighYHighXSplit == size_t_max) {
              //+++
              cell->child[0] = new Cell(zHighYHighXSplit + 1, 
                                        cell->last,
                                        transl(cell->center, constr_vec(
                                          ORIGIN, constr_pt(width, width, width))),
                                        cell->level + 1);

              if (cell->child[0]->size() > m_bucket_size)
                stack.push(cell->child[0]);
            }
          }

          std::size_t sum = 0;
          for (std::size_t i = 0;i<8;i++)
            sum += (cell->child[i]) ? cell->child[i]->size() : 0;

          count++;
        }
      }

      bool drawSamplesFromCellContainingPoint(const Point_3 &p,
                                              std::size_t level,
                                              std::set<std::size_t>& indices,
                                              const std::vector<int>& shapeIndex,
                                              std::size_t requiredSamples) {

        bool upperZ, upperY, upperX;
        Cell *cur = m_root;

        while (cur && cur->level < level) {
          upperX = get_x(cur->center) <= get_x(p);
          upperY = get_y(cur->center) <= get_y(p);
          upperZ = get_z(cur->center) <= get_z(p);

          if (upperZ) {
            if (upperY)
              cur = (upperX) ? cur->child[0] : cur->child[1];
            else cur = (upperX) ? cur->child[2] : cur->child[3];
          }
          else {
            if (upperY)
              cur = (upperX) ? cur->child[4] : cur->child[5];
            else cur = (upperX) ? cur->child[6] : cur->child[7];
          }
        }

        if (cur) {
          std::size_t enough = 0;
          for (std::size_t i = cur->first;i<=cur->last;i++) {
            std::size_t j = this->index(i);
            if (shapeIndex[j] == -1) {
              enough++;
              if (enough >= requiredSamples)
                break;
            }
          }
          if (enough >= requiredSamples) {
            do {
              std::size_t p = CGAL::get_default_random().
                uniform_int<std::size_t>(0, cur->size() - 1);
              std::size_t j = this->index(cur->first + p);
              if (shapeIndex[j] == -1)
                indices.insert(j);
            } while (indices.size() < requiredSamples);

            return true;
          }
          else return false;
        }
        else return false;
      }

      std::size_t maxLevel() {
        return m_max_level;
      }

      std::size_t fullScore(Shape *candidate,
                       std::vector<int> &shapeIndex,
                       FT epsilon, 
                       FT normal_threshold) {

        std::vector<std::size_t> indices(m_root->size());
        for (std::size_t i = 0;i<m_root->size();i++) {
          indices[i] = index(m_root->first + i);
        }
        
        candidate->cost_function(this->begin() + m_root->first, 
                                 this->begin() + m_root->last,
                                 shapeIndex,
                                 epsilon, 
                                 normal_threshold, 
                                 indices);

        return candidate->m_indices.size();
      }

      std::size_t score(Shape *candidate, 
                   std::vector<int> &shapeIndex,
                   FT epsilon,
                   FT normal_threshold) {

        std::stack<Cell *> stack;
        stack.push(m_root);

        while(!stack.empty()) {
          Cell *cell = stack.top();
          stack.pop();

          FT width = m_width / (1<<(cell->level));

          FT diag = CGAL::sqrt(FT(3) * width * width) + epsilon;

          FT dist = candidate->squared_distance(cell->center);

          if (dist > (diag * diag))
            continue;

          // differ between full or partial overlap?
          // if full overlap further traversal of this branch is not necessary
          if (cell->isLeaf()) {
            std::vector<std::size_t> indices;
            indices.reserve(cell->size());
            for (std::size_t i = 0;i<cell->size();i++) {
              if (shapeIndex[this->index(cell->first + i)] == -1) {
                indices.push_back(this->index(cell->first + i));
              }
            }

            candidate->cost_function(epsilon, 
                                     normal_threshold, 
                                     indices);
          }
          else {
            for (std::size_t i = 0;i<8;i++)
              if (cell->child[i])
                stack.push(cell->child[i]);
          }

        }

        return candidate->m_indices.size();
      }

      void setBucketSize(std::size_t bucketSize) {
        m_bucket_size = bucketSize;
      }

      const Bbox_3 &boundingBox() {
        return m_bBox;
      }
        
      const Bbox_3 &buildBoundingCube() {
        FT min[] = {(std::numeric_limits<FT>::max)(),
                    (std::numeric_limits<FT>::max)(),
                    (std::numeric_limits<FT>::max)()};
        FT max[] = {(std::numeric_limits<FT>::min)(),
                    (std::numeric_limits<FT>::min)(),
                    (std::numeric_limits<FT>::min)()};

        for (std::size_t i = 0;i<this->size();i++) {
          Point_3 p = get(m_point_pmap, *this->at(i));
          min[0] = (std::min<FT>)(min[0], get_x(p));
          min[1] = (std::min<FT>)(min[1], get_y(p));
          min[2] = (std::min<FT>)(min[2], get_z(p));
          max[0] = (std::max<FT>)(max[0], get_x(p));
          max[1] = (std::max<FT>)(max[1], get_y(p));
          max[2] = (std::max<FT>)(max[2], get_z(p));
        }

        m_bBox = Bbox_3(min[0], min[1], min[2], max[0], max[1], max[2]);

        m_width = (std::max)(max[0] - min[0], 
          (std::max)(max[1] - min[1], max[2] - min[2])) * (FT) 0.5;
        m_center = constr_pt((min[0] + max[0]) * (FT) 0.5,
                             (min[1] + max[1]) * (FT) 0.5,
                             (min[2] + max[2]) * (FT) 0.5);

        return m_bBox;
      }

      // returns index of last point below threshold
      std::size_t split(std::size_t first, std::size_t last, std::size_t dimension, FT threshold) {
        if (last == size_t_max || first == size_t_max) 
          return size_t_max;

        if (first > last)
          return first - 1;

        std::size_t origFirst = first;

        while(first < last) {
          // find first above threshold
          while (get_coord(
                   get(m_point_pmap, *this->at(first)),
                   static_cast<unsigned int>(dimension)) < threshold
                 && first < last) {
            first++;
          }

          // check if last has been reached
          if (first == last) {
            return (get_coord(
                      get(m_point_pmap, *this->at(first)),
                      static_cast<unsigned int>(dimension)) < threshold) ?
                   first : (first == origFirst) ? size_t_max : first - 1;
          }

          // find last below threshold
          while (get_coord(
                   get(m_point_pmap, *this->at(last)),
                   static_cast<unsigned int>(dimension)) >= threshold
                 && last > first) {
            last--;
          }

          // check if first has been reached
          if (last == first) {
            return (get_coord(
                      get(m_point_pmap, *this->at(first)),
                      static_cast<unsigned int>(dimension)) < threshold) ?
                    first : (first == origFirst) ? size_t_max : first - 1;
          }

          this->swap(first, last);
          first++;
          last--;
        }

        return (get_coord(
                  get(m_point_pmap, *this->at(first)),
                  static_cast<unsigned int>(dimension)) < threshold) ?
               first : (first == origFirst) ? size_t_max : first - 1;
      }

      Sd_traits m_traits;
      Bbox_3 m_bBox;
      Cell *m_root;
      Point_3 m_center;
      FT m_width;
      std::size_t m_bucket_size;
      std::size_t m_set_max_level;
      std::size_t m_max_level;
      Point_map m_point_pmap;
      Normal_map m_normal_pmap;
    };
  }
  }
}

#endif
