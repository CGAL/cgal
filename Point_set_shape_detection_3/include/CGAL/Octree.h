#ifndef CGAL_SHAPE_DETECTION_3_OCTREE_H
#define CGAL_SHAPE_DETECTION_3_OCTREE_H

#include <limits>
#include <random>
#include <stack>

#include <CGAL/Bbox_3.h>

#include "Shape_base.h"

extern int scoreTime;

namespace CGAL {
  
    template<class Sd_traits> 
    class Shape_detection_3;
    
  namespace internal {
  
    template<class Sdt>
    class DirectPointAccessor {
      // think if constructor could be protected (and the rest of the methods)
    public:
      typedef Sdt Sd_traits;
      typedef typename Sd_traits::Input_iterator Input_iterator;

      DirectPointAccessor() {}
      DirectPointAccessor(const Input_iterator &begin, const Input_iterator &beyond, unsigned int offset) : m_first(begin), m_offset(offset) {
        m_beyond = (beyond == begin) ? begin : beyond - 1;
      }

      Input_iterator at(unsigned int i) {
        return m_first + i;
      }

      unsigned int index(unsigned int i) {
        return i + m_offset;
      }

      unsigned int offset() {
        return m_offset;
      }

      unsigned int size() {
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

      void swap(unsigned int a, unsigned int b) {
        typename std::iterator_traits<Input_iterator>::value_type tmp;
        tmp = m_first[a];
        m_first[a] = m_first[b];
        m_first[b] = tmp;
      }

    protected:
      Input_iterator m_first;

    private:
      Input_iterator m_beyond;
      unsigned int m_offset;
    };

    template<class Sdt>
    class IndexedPointAccessor {
    public:
      typedef Sdt Sd_traits;
      typedef typename Sd_traits::Input_iterator Input_iterator;

      IndexedPointAccessor() {}
      IndexedPointAccessor(const Input_iterator &begin, const Input_iterator &beyond, unsigned int) : m_first(begin) {
        m_beyond = (beyond == begin) ? begin : beyond - 1;
        m_indices.resize(size());
        for (unsigned int i = 0;i<size();i++)
          m_indices[i] = i;
      }

      Input_iterator at(unsigned int i) {
        return m_first + m_indices[i];
      }

      unsigned int index(unsigned int i) {
        return m_indices[i];
      }

      unsigned int offset() {
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
        for (unsigned int i = 0;i<size();i++)
          m_indices[i] = i;
      }

      unsigned int size() {
        return m_beyond - m_first + 1;
      }

      void swap(unsigned int a, unsigned int b) {
        unsigned int tmp = m_indices[a];
        m_indices[a] = m_indices[b];
        m_indices[b] = tmp;
      }

    protected:
      Input_iterator m_first;

    private:
      std::vector<unsigned int> m_indices;
      Input_iterator m_beyond;
    };

    template<class PointAccessor>
    class Octree : public PointAccessor {

      typedef typename PointAccessor::Sd_traits Sd_traits;
      typedef typename Sd_traits::Input_iterator Input_iterator;
      typedef Shape_base<Sd_traits> Primitive;
      typedef typename Sd_traits::Geom_traits::Point_3 Point;
      typedef typename Sd_traits::Geom_traits::Vector_3 Vector;
      typedef typename Sd_traits::Geom_traits::FT FT;
      typedef typename Sd_traits::Point_pmap Point_pmap;
      typedef typename Sd_traits::Normal_pmap Normal_pmap;

      template<class Sd_traits>
        friend class ::CGAL::Shape_detection_3;

      struct Cell {
        int first, last;
        Cell *child[8];
        Point center;
        unsigned int level;
        Cell(int first, int last, Point center, unsigned int level) : first(first), last(last), center(center), level(level) {memset(child, 0, sizeof(Cell *) * 8);}
        bool isLeaf() const {
          for (unsigned int i = 0;i<8;i++) {
            if (child[i])
              return false;
          }
          return true;
        }

        unsigned int size() const {
          if (last < first) {
            std::cout << "first > last!" << std::endl;
          }
          if (first == -1 || last == -1)
            return -1;
          else return (last - first + 1);
        }
      };
        
    public:
      Octree() : m_bucketSize(20), m_setMaxLevel(10), m_root(NULL) {}
      Octree(const Input_iterator &first, const Input_iterator &beyond, unsigned int offset = 0, unsigned int bucketSize = 20, unsigned int maxLevel = 10) : PointAccessor(first, beyond, offset), m_root(NULL), m_bucketSize(bucketSize), m_setMaxLevel(maxLevel) {}
      ~Octree() {
        if (!m_root)
          return;

        std::stack<Cell *> stack;
        stack.push(m_root);
        while (!stack.empty()) {
          Cell *cell = stack.top();
          stack.pop();

          for (unsigned int i = 0;i<8;i++)
            if (cell->child[i])
              stack.push(cell->child[i]);

          delete cell;
        }
      }

      // Sorting data in a way such that points of one cell are always in one range and ordered child-wise:
      // +---+---+
      // | 1.| 0.|
      // +---+---+
      // | 3.| 2.|
      // +---+---+
      // z max before z min, then y max before y min, then x max before x min
      void createTree() {
        buildBoundingCube();
        int count = 0;
        m_maxLevel = 0;

        std::stack<Cell *> stack;
        m_root = new Cell(0, this->size() - 1, m_center, 0);
        stack.push(m_root);
        while (!stack.empty()) {
          Cell *cell= stack.top();
          stack.pop();

          m_maxLevel = std::max<unsigned int>(m_maxLevel, cell->level);
          if (cell->level == m_setMaxLevel)
            continue;

          //verifyCell(cell);

          int zLowYHighXSplit, zLowYLowXSplit, zLowYSplit, zHighYSplit, zHighYHighXSplit, zHighYLowXSplit;

          int zSplit = split(cell->first, cell->last, 2, cell->center[2]);

          if (zSplit != -1) {
//             verifyRange(cell->first, zSplit, 2, cell->center[2], true);
//             verifyRange(zSplit + 1, cell->last, 2, cell->center[2], false);

            zLowYSplit = split(cell->first, zSplit, 1, cell->center[1]);
            if (zLowYSplit != -1) {
//               verifyRange(cell->first, zLowYSplit, 1, cell->center[1], true);
//               verifyRange(zLowYSplit + 1, zSplit, 1, cell->center[1], false);

              zLowYLowXSplit = split(cell->first, zLowYSplit, 0, cell->center[0]);
//               if (zLowYLowXSplit != -1) {
//                 verifyRange(cell->first, zLowYLowXSplit, 0, cell->center[0], true);
//                 verifyRange(zLowYLowXSplit + 1, zLowYSplit, 0, cell->center[0], false);
//               }

              zLowYHighXSplit = split(zLowYSplit + 1, zSplit, 0, cell->center[0]);              
//               if (zLowYHighXSplit != -1) {
//                 verifyRange(zLowYSplit + 1, zLowYHighXSplit, 0, cell->center[0], true);
//                 verifyRange(zLowYHighXSplit + 1, zSplit, 0, cell->center[0], false);
//               }
            }
            else {
              zLowYLowXSplit = -1;
              zLowYHighXSplit = split(cell->first, zSplit, 0, cell->center[0]);           
//               if (zLowYHighXSplit != -1) {
//                 verifyRange(cell->first, zLowYHighXSplit, 0, cell->center[0], true);
//                 verifyRange(zLowYHighXSplit + 1, zSplit, 0, cell->center[0], false);
//               }
            }

            zHighYSplit = split(zSplit + 1, cell->last, 1, cell->center[1]);         
            if (zHighYSplit != -1) {
//               verifyRange(zSplit + 1, zHighYSplit, 1, cell->center[1], true);
//               verifyRange(zHighYSplit + 1, cell->last, 1, cell->center[1], false);

              zHighYHighXSplit = split(zHighYSplit + 1, cell->last, 0, cell->center[0]);       
//               if (zHighYHighXSplit != -1) {
//                 verifyRange(zHighYSplit + 1, zHighYHighXSplit, 0, cell->center[0], true);
//                 verifyRange(zHighYHighXSplit + 1, cell->last, 0, cell->center[0], false);
//               }

              zHighYLowXSplit = split(zSplit + 1, zHighYSplit, 0, cell->center[0]);       
//               if (zHighYLowXSplit != -1) {
//                 verifyRange(zSplit + 1, zHighYLowXSplit, 0, cell->center[0], true);
//                 verifyRange(zHighYLowXSplit + 1, zHighYSplit, 0, cell->center[0], false);
//               }
            }
            else {
              zHighYLowXSplit = -1;
              zHighYHighXSplit = split(zSplit + 1, cell->last, 0, cell->center[0]);    
//               if (zHighYLowXSplit != -1) {
//                 verifyRange(zSplit + 1, zHighYLowXSplit, 0, cell->center[0], true);
//                 verifyRange(zHighYLowXSplit + 1, cell->last, 0, cell->center[0], false);
//               }
            }
          }
          else {
//             verifyRange(cell->first, cell->last, 2, cell->center[2], false);
            zLowYSplit = -1;
            zLowYLowXSplit = -1;
            zLowYHighXSplit = -1;

            zHighYSplit = split(cell->first, cell->last, 1, cell->center[1]);     
//             if (zHighYSplit != -1) {
//               verifyRange(cell->first, zHighYSplit, 1, cell->center[1], true);
//               verifyRange(zHighYSplit + 1, cell->last, 1, cell->center[1], false);
//             }

            if (zHighYSplit != -1) {
              zHighYHighXSplit = split(zHighYSplit + 1, cell->last, 0, cell->center[0]); 
//               if (zHighYHighXSplit != -1) {
//                 verifyRange(zHighYSplit + 1, zHighYHighXSplit, 0, cell->center[0], true);
//                 verifyRange(zHighYHighXSplit + 1, cell->last, 0, cell->center[0], false);
//               }

              zHighYLowXSplit = split(cell->first, zHighYSplit, 0, cell->center[0]); 
//               if (zHighYLowXSplit != -1) {
//                 verifyRange(cell->first, zHighYLowXSplit, 0, cell->center[0], true);
//                 verifyRange(zHighYLowXSplit + 1, zHighYSplit, 0, cell->center[0], false);
//               }
            }
            else {
              zHighYLowXSplit = -1;
              zHighYHighXSplit = split(cell->first, cell->last, 0, cell->center[0]); 
//               if (zHighYHighXSplit != -1) {
//                 verifyRange(cell->first, zHighYHighXSplit, 0, cell->center[0], true);
//                 verifyRange(zHighYHighXSplit + 1, cell->last, 0, cell->center[0], false);
//               }
            }
          }

          //verifyCell(cell);

          FT width = m_width / (1<<(cell->level + 1));

          if (zSplit != -1) {
            if (zLowYSplit != -1) {
              if (zLowYLowXSplit != -1) {

                if (cell->first <= zLowYLowXSplit) {
                  //---
                  cell->child[7] = new Cell(cell->first, zLowYLowXSplit, cell->center + Vector(-width,-width,-width), cell->level + 1);
                  //verifyCell(cell->child[7]);
                  if (cell->child[7]->size() > m_bucketSize)
                    stack.push(cell->child[7]);
                }
              }
              else zLowYLowXSplit = cell->first - 1;

              if (zLowYLowXSplit < zLowYSplit) {
                //+--
                cell->child[6] = new Cell(zLowYLowXSplit + 1, zLowYSplit, cell->center + Vector(width,-width,-width), cell->level + 1);
                //verifyCell(cell->child[6]);
                if (cell->child[6]->size() > m_bucketSize)
                  stack.push(cell->child[6]);
              }
            }
            else zLowYSplit = cell->first - 1;

            if (zLowYHighXSplit != -1) {

              if (zLowYSplit < zLowYHighXSplit) {
                //-+-
                cell->child[5] = new Cell(zLowYSplit + 1, zLowYHighXSplit, cell->center + Vector(-width, width,-width), cell->level + 1);
                //verifyCell(cell->child[5]);
                if (cell->child[5]->size() > m_bucketSize)
                  stack.push(cell->child[5]);
              }
            }
            else zLowYHighXSplit = zLowYSplit;

            if (zLowYHighXSplit < zSplit) {
              //++-
              cell->child[4] = new Cell(zLowYHighXSplit + 1, zSplit, cell->center + Vector(width, width,-width), cell->level + 1);
              //verifyCell(cell->child[4]);
              if (cell->child[4]->size() > m_bucketSize)
                stack.push(cell->child[4]);
            }
          }
          else zSplit = cell->first - 1;

          if (zHighYSplit != -1) {
            if (zHighYLowXSplit != -1) {

              if (zSplit < zHighYLowXSplit) {
                //--+
                cell->child[3] = new Cell(zSplit + 1, zHighYLowXSplit, cell->center + Vector(-width,-width, width), cell->level + 1);
                //verifyCell(cell->child[3]);
                if (cell->child[3]->size() > m_bucketSize)
                  stack.push(cell->child[3]);
              }
            }
            else zHighYLowXSplit = zSplit;

            if (zHighYLowXSplit < zHighYSplit) {
              //+-+
              cell->child[2] = new Cell(zHighYLowXSplit + 1, zHighYSplit, cell->center + Vector(width,-width, width), cell->level + 1);
              //verifyCell(cell->child[2]);
              if (cell->child[2]->size() > m_bucketSize)
                stack.push(cell->child[2]);
            }

          }
          else zHighYSplit = zSplit;

          if (zHighYHighXSplit != -1) {

            if (zHighYSplit < zHighYHighXSplit) {
              //-++
              cell->child[1] = new Cell(zHighYSplit + 1, zHighYHighXSplit, cell->center + Vector(-width, width, width), cell->level + 1);
              //verifyCell(cell->child[1]);
              if (cell->child[1]->size() > m_bucketSize)
                stack.push(cell->child[1]);
            }
          }
          else zHighYHighXSplit = zHighYSplit;

          if (zHighYHighXSplit <= cell->last) {
            if (zHighYHighXSplit < cell->last) {
              //+++
              cell->child[0] = new Cell(zHighYHighXSplit + 1, cell->last, cell->center + Vector(width, width, width), cell->level + 1);
              //verifyCell(cell->child[0]);
              if (cell->child[0]->size() > m_bucketSize)
                stack.push(cell->child[0]);
            }
          }

          int sum = 0;
          for (unsigned int i = 0;i<8;i++)
            sum += (cell->child[i]) ? cell->child[i]->size() : 0;

          if (sum != cell->size()) {
            std::cout << "sum of size of childs wrong!" << std::endl;
          }

          count++;
        }
      }

      bool drawSamplesFromCellContainingPoint(const Point &p, unsigned int level, std::set<int> &indices, const std::vector<int> &shapeIndex, int requiredSamples) {
        static std::random_device rd;
        static std::mt19937 rng(rd());
        static std::uniform_int_distribution<> dis(0, this->size() - 1);

        bool upperZ, upperY, upperX;
        Cell *cur = m_root;

        while (cur && cur->level < level) {
          upperX = cur->center[0] <= p[0];
          upperY = cur->center[1] <= p[1];
          upperZ = cur->center[2] <= p[2];

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
          int enough = 0;
          for (unsigned int i = cur->first;i<=cur->last;i++) {
            int j = this->index(i);
            if (shapeIndex[j] == -1) {
              enough++;
              if (enough >= requiredSamples)
                break;
            }
          }
          if (enough >= requiredSamples) {
            do {
              int p = dis(rng) % cur->size();
              int j = this->index(cur->first + p);
              if (shapeIndex[j] == -1)
                indices.insert(j);
            } while (indices.size() < requiredSamples);

            return true;
          }
          else return false;
        }
        else return false;
      }

      unsigned int maxLevel() {
        return m_maxLevel;
      }

      unsigned int fullScore(Primitive *candidate, std::vector<int> &shapeIndex, FT epsilon, FT normal_threshold) {
        std::vector<unsigned int> indices(m_root->size());
        for (unsigned int i = 0;i<m_root->size();i++) {
          indices[i] = index(m_root->first + i);
        }
        
        candidate->cost_function(this->begin() + m_root->first, this->begin() + m_root->last, shapeIndex, epsilon, normal_threshold, indices);

        return candidate->m_indices.size();
      }

      unsigned int score(Primitive *candidate, std::vector<int> &shapeIndex, FT epsilon, FT normal_threshold) {
        //clock_t s = clock();
        std::stack<Cell *> stack;
        stack.push(m_root);

        while(!stack.empty()) {
          Cell *cell = stack.top();
          stack.pop();

          FT width = m_width / (1<<(cell->level));

          FT diag = ::sqrt(3 * width * width) + epsilon;

          FT dist = candidate->squared_distance(cell->center);

          if (dist > (diag * diag))
            continue;

          // differ between full or partial overlap?
          // if full overlap further traversal of this branch is not necessary
          if (cell->isLeaf()) {
            std::vector<unsigned int> indices;
            indices.reserve(cell->size());
            for (unsigned int i = 0;i<cell->size();i++) {
              if (shapeIndex[this->index(cell->first + i)] == -1) {
                indices.push_back(this->index(cell->first + i));
              }
            }

            candidate->cost_function(shapeIndex, epsilon, normal_threshold, indices);
          }
          else {
            for (unsigned int i = 0;i<8;i++)
              if (cell->child[i])
                stack.push(cell->child[i]);
          }

        }

        //scoreTime += clock() - s;

        return candidate->m_indices.size();
      }

      void setBucketSize(unsigned int bucketSize) {
        m_bucketSize = bucketSize;
      }

      void verify() {
        std::set<unsigned int> indices;
        for (unsigned int i = m_root->first;i<=m_root->last;i++) {
          indices.insert(this->index(i));
        }
        std::cout << "size: " << m_root->size() << " unique indices: " << indices.size() << std::endl;
        std::stack<Cell *> stack;
        stack.push(m_root);
        int highestLvl = 0;

        while (!stack.empty()) {
          Cell *cell = stack.top();
          stack.pop();

          highestLvl = std::max<int>(highestLvl, cell->level);

          FT width = m_width / (1<<(cell->level));

          FT diag = sqrt(3 * width * width);

          for (unsigned int i = cell->first;i<cell->last;i++) {
            FT d = sqrt(CGAL::squared_distance(cell->center, this->at(i).first));
            if (d > diag) {
              std::cout << "points out of range" << std::endl;
            }
          }
          for (unsigned int i = 0;i<8;i++) {
            if (cell->child[i])
              stack.push(cell->child[i]);
          }
        }

        std::cout << "highestLvl: " << highestLvl << std::endl;
      }

      void verifyCell(Cell *cell) {
        FT width = m_width / (1<<(cell->level));

        FT diag = sqrt(3 * width * width);

        Point c = cell->center;

        for (int i = cell->first;i<cell->last;i++) {
          FT d = sqrt(CGAL::squared_distance(cell->center, get(m_pointPMap, this->at(i).first)));
          if (d > diag) {
            std::cout << "points out of range" << std::endl;
          }
        }
      }

      void verifyRange(int first, int last, int dimension, FT threshold, bool below) {
        Point p;
        FT v;

        if (below)
          threshold *= 1.000001;
        else
          threshold *= 0.999999;

        for (;first <= last;first++) {
          bool fail = false;
          p = get(m_pointPMap, this->at(first));
          v = p[dimension];
          if (below){
            if (get(m_pointPMap, this->at(first))[dimension] > threshold)
              fail = true;
          }
          else
            if (get(m_pointPMap, this->at(first))[dimension] <= threshold)
              fail = true;

          if (fail) {
            std::cout << "range check failed" << std::endl;
          }
        }
      }

        
      void buildBoundingCube() {
        FT min[] = {(std::numeric_limits<FT>::max)(), (std::numeric_limits<FT>::max)(), (std::numeric_limits<FT>::max)()};
        FT max[] = {(std::numeric_limits<FT>::min)(), (std::numeric_limits<FT>::min)(), (std::numeric_limits<FT>::min)()};

        for (unsigned int i = 0;i<this->size();i++) {
            Point p = get(m_pointPMap, *this->at(i));
          min[0] = (std::min<FT>)(min[0], p.x());
          min[1] = (std::min<FT>)(min[1], p.y());
          min[2] = (std::min<FT>)(min[2], p.z());
          max[0] = (std::max<FT>)(max[0], p.x());
          max[1] = (std::max<FT>)(max[1], p.y());
          max[2] = (std::max<FT>)(max[2], p.z());
        }

        m_width = (std::max)(max[0] - min[0], (std::max)(max[1] - min[1], max[2] - min[2])) * 0.5;

        //         m_center[0] = (min[0] + max[0]) * 0.5;
        //         m_center[1] = (min[1] + max[1]) * 0.5;
        //         m_center[2] = (min[2] + max[2]) * 0.5;

        m_center = Point((min[0] + max[0]) * 0.5, (min[1] + max[1]) * 0.5, (min[2] + max[2]) * 0.5);

        //m_bBox = Bbox_3(m_center[0] - m_width, m_center[1] - m_width, m_center[2] - m_width, m_center[0] + m_width, m_center[1] + m_width, m_center[2] + m_width);
      }

      // returns index of last point below threshold
      int split(int first, int last, unsigned int dimension, FT threshold) {
        if (last == -1 || first == -1)
          return -1;

        if (first > last)
          return first - 1;

        int origFirst = first;
        int origLast = last;

        while(first < last) {
          // find first above threshold
          Point p1 = get(m_pointPMap, *this->at(first));
          FT v1 = p1[dimension];
          while (get(m_pointPMap, *this->at(first))[dimension] < threshold && first < last) {
            first++;
          }

          // check if last has been reached
          if (first == last) {
            return (get(m_pointPMap, *this->at(first))[dimension] < threshold) ? first : (first == origFirst) ? -1 : first - 1;
          }

          // find last below threshold
          p1 = get(m_pointPMap, *this->at(last));
          v1 = p1[dimension];
          while (get(m_pointPMap, *this->at(last))[dimension] >= threshold && last > first) {
            last--;
          }

          // check if first has been reached
          if (last == first) {
            return (get(m_pointPMap, *this->at(first))[dimension] < threshold) ? first : (first == origFirst) ? -1 : first - 1;
          }

          // swap entries
          if (first < origFirst || first > origLast || last < origFirst || last > origLast) {
            std::cout << "split: swap out of bounds!" << std::endl;
          }
          this->swap(first, last);
          p1 = get(m_pointPMap, *this->at(first));
          v1 = p1[dimension];
          p1 = get(m_pointPMap, *this->at(last));
          v1 = p1[dimension];
          first++;
          last--;
        }

        return (get(m_pointPMap, *this->at(first))[dimension] < threshold) ? first : (first == origFirst) ? -1 : first - 1;
      }

      //Bbox_3 m_bBox;
      Cell *m_root;
      Point m_center;
      FT m_width;
      unsigned int m_bucketSize;
      unsigned int m_setMaxLevel;
      unsigned int m_maxLevel;
      Point_pmap m_pointPMap;
      Normal_pmap m_normalPMap;
    };
  }
}

#endif
