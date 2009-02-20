// Copyright (c) 2007-2008  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://ggael@scm.gforge.inria.fr/svn/cgal/trunk/Surface_reconstruction_3/include/CGAL/Neighborhood/point_kdtree_3.h $
// $Id: point_kdtree_3.h 47904 2009-02-04 18:46:54Z ggael $
//
//
// Author(s)     : Gael Guennebaud

#ifndef CGAL_POINT_KDETREE_3_H
#define CGAL_POINT_KDETREE_3_H

#include <vector>
#include <limits>

namespace CGAL {

/** Implements a bounded-size max priority queue using a heap */
template <typename Index, typename Weight>
class HeapMaxPriorityQueue
{
  struct Element
  {
      Weight weight;
      Index index;
  };

public:

  HeapMaxPriorityQueue(void)
  {
    mElements = 0;
    mMaxSize = 0;
    mpOffsetedElements = NULL;
  }

  inline void setMaxSize(int maxSize)
  {
    if (mMaxSize!=maxSize)
    {
      mMaxSize = maxSize;
      delete[] mElements;
      mElements = new Element[mMaxSize];
      mpOffsetedElements = (mElements-1);
    }
    init();
  }

  inline void init() { mCount = 0; }

  inline bool isFull() const { return mCount == mMaxSize; }

  /** returns number of elements inserted in queue */
  inline int getNofElements() const { return mCount; }

  inline Weight getWeight(int i) const { return mElements[i].weight; }
  inline Index  getIndex(int i) const { return mElements[i].index; }

  inline Weight getTopWeight() const { return mElements[0].weight; }

  inline void insert(Index index, Weight weight)
  {
    if (mCount==mMaxSize)
    {
      if (weight<mElements[0].weight)
      {
        register int j, k;
        j = 1;
        k = 2;
        while (k <= mMaxSize)
        {
          Element* z = &(mpOffsetedElements[k]);
          if ((k < mMaxSize) && (z->weight < mpOffsetedElements[k+1].weight))
            z = &(mpOffsetedElements[++k]);

          if(weight >= z->weight)
            break;
          mpOffsetedElements[j] = *z;
          j = k;
          k = 2 * j;
        }
        mpOffsetedElements[j].weight = weight;
        mpOffsetedElements[j].index = index;
      }
    }
    else
    {
      int i, j;
      i = ++mCount;
      while (i >= 2)
      {
        j = i >> 1;
        Element& y = mpOffsetedElements[j];
        if(weight <= y.weight)
          break;
        mpOffsetedElements[i] = y;
        i = j;
      }
      mpOffsetedElements[i].index = index;
      mpOffsetedElements[i].weight = weight;
    }
  }

protected:

  int mCount;
  int mMaxSize;
  Element* mElements;
  Element* mpOffsetedElements;
};

template<typename _DataType>
class ConstDataWrapper
{
public:
    typedef _DataType DataType;
    inline ConstDataWrapper()
      : mpData(0), mStride(0), mSize(0)
    {}
    inline ConstDataWrapper(const DataType* pData, int size, int stride = sizeof(DataType))
      : mpData(reinterpret_cast<const unsigned char*>(pData)), mStride(stride), mSize(size)
    {}
    inline const DataType& operator[] (int i) const
    {
      return *reinterpret_cast<const DataType*>(mpData + i*mStride);
    }
    inline size_t size() const { return mSize; }
protected:
    const unsigned char* mpData;
    int mStride;
    size_t mSize;
};

template<typename SearchTraits>
class KnnNeighborhood : private HeapMaxPriorityQueue<int,typename SearchTraits::FT>
{
public:
  typedef HeapMaxPriorityQueue<int,typename SearchTraits::FT> Base;
  typedef typename SearchTraits::FT FT;

  using Base::init;
  using Base::insert;

  inline KnnNeighborhood(int k)
  {
    Base::setMaxSize(k);
  }

  inline void finalize() {}
  inline int size() const { return Base::getNofElements(); }
  inline int index(int i) const { return Base::getIndex(i); }
  inline FT maxValue() const { return Base::getTopWeight(); }
  inline FT value(int i) const { return Base::getWeight(i); }
};

template<typename SearchTraits>
class PointKdtree_3
{
public:

  typedef typename SearchTraits::FT FT;
  typedef typename SearchTraits::Point_d ObjectType;
  typedef typename SearchTraits::PointType PointType;

protected:
  struct Element : public PointType
  {
    Element(){}
    Element(const PointType& pos, int id) : PointType(pos), index(id) {}
    int index;
  };
  typedef std::vector<Element> ElementList;

  struct Node
  {
    union {
      struct {
        FT splitValue;
        unsigned int firstChildId:24;
        unsigned int dim:2;
        unsigned int leaf:1;
      };
      struct {
        unsigned int start;
        unsigned int size;
      };
    };
  };
  typedef std::vector<Node> NodeList;

public:

  PointKdtree_3(const ConstDataWrapper<PointType>& points, const ObjectType* pObjects = 0, unsigned int nofPointsPerCell = 16, unsigned int maxDepth = 64)
    : mObjects(pObjects), mElements(points.size())
  {
    // compute the AABB of the input
    mElements[0] = Element(points[0],0);
    for (int k=0; k<3; ++k)
      mAABB_max[k] = mAABB_min[k] = mElements[0][k];
    for (unsigned int i=1 ; i<mElements.size() ; ++i)
    {
      mElements[i] = Element(points[i],i);
      for (int k=0; k<3; ++k)
      {
        if (mElements[i][k]<mAABB_min[k])
          mAABB_min[k] = mElements[i][k];
        if (mElements[i][k]>mAABB_max[k])
          mAABB_max[k] = mElements[i][k];
      }
    }

    mNodes.reserve(4*mElements.size()/nofPointsPerCell);
    mNodes.resize(1);
    mNodes.back().leaf = 0;
    createTree(0, 0, mElements.size(), 1, nofPointsPerCell, maxDepth);
  }

  ~PointKdtree_3() {}

  const ObjectType& object(int i) { assert(mObjects); return mObjects[i]; }

  /* Performs the kNN query.
    *
    * This algorithm uses the simple distance to the split plane to prune nodes.
    * A more elaborated approach consists to track the closest corner of the cell
    * relatively to the current query point. This strategy allows to save about 5%
    * of the leaves. However, in practice the slight overhead due to this tracking
    * reduces the overall performance.
    *
    * This algorithm also use a simple stack while a priority queue using the squared
    * distances to the cells as a priority values allows to save about 10% of the leaves.
    * But, again, priority queue insertions and deletions are quite involved, and therefore
    * a simple stack is by far much faster.
    */
  template<typename NeighborList>
  void doQueryK(const PointType& queryPoint, NeighborList& list)
  {
    list.init();
    list.insert(0xffffffff, (std::numeric_limits<float>::max)());

    mNodeStack[0].nodeId = 0;
    mNodeStack[0].sq = 0.f;
    unsigned int count = 1;

    while (count)
    {
      QueryNode& qnode = mNodeStack[count-1];
      Node& node = mNodes[qnode.nodeId];

      if (qnode.sq < list.maxValue())
      {
        if (node.leaf)
        {
          --count; // pop
          unsigned int end = node.start + node.size;
          for (unsigned int i=node.start ; i<end ; ++i)
          {
            FT aux = queryPoint[0] - mElements[i][0];
            FT d2 = aux*aux;
            aux = queryPoint[1] - mElements[i][1];
            d2 += aux*aux;
            aux = queryPoint[2] - mElements[i][2];
            d2 += aux*aux;
            list.insert(mElements[i].index, d2);
          }
        }
        else
        {
          // replace the stack top by the farthest and push the closest
          float new_off = queryPoint[node.dim] - node.splitValue;
          if (new_off < 0.)
          {
              mNodeStack[count].nodeId  = node.firstChildId;
              qnode.nodeId = node.firstChildId+1;
          }
          else
          {
              mNodeStack[count].nodeId  = node.firstChildId+1;
              qnode.nodeId = node.firstChildId;
          }
          mNodeStack[count].sq = qnode.sq;
          qnode.sq = new_off*new_off;
          ++count;
        }
      }
      else
      {
        // pop
        --count;
      }
    }
    list.finalize();
  //   std::cout << "doQueryK OK\n";
  }

protected:

  // element of the stack
  struct QueryNode
  {
      QueryNode() {}
      QueryNode(unsigned int id) : nodeId(id) {}
      unsigned int nodeId;  // id of the next node
      FT sq;                // squared distance to the next node
  };

  // used to build the tree: split the subset [start..end[ according to dim and splitValue,
  // and returns the index of the first element of the second subset
  unsigned int split(int start, int end, unsigned int dim, float splitValue);

  void createTree(unsigned int nodeId, unsigned int start, unsigned int end, unsigned int level, unsigned int targetCellsize, unsigned int targetMaxDepth);

protected:

  const ObjectType* mObjects;
  //AxisAlignedBoxType mAABB;
  FT mAABB_min[3];
  FT mAABB_max[3];
  NodeList mNodes;
  ElementList mElements;
  QueryNode mNodeStack[64];
};

template<typename SearchTraits>
unsigned int PointKdtree_3<SearchTraits>::split(int start, int end, unsigned int dim, float splitValue)
{
  int l(start), r(end-1);
  for ( ; l<r ; ++l, --r)
  {
    while (l < end && mElements[l][dim] < splitValue)
      l++;
    while (r >= start && mElements[r][dim] >= splitValue)
      r--;
    if (l > r)
      break;
    std::swap(mElements[l],mElements[r]);
  }
  return (mElements[l][dim] < splitValue ? l+1 : l);
}

/** recursively builds the kdtree
  *
  *  The heuristic is the following:
  *   - if the number of points in the node is lower than targetCellsize then make a leaf
  *   - else compute the AABB of the points of the node and split it at the middle of
  *     the largest AABB dimension.
  *
  *  This strategy might look not optimal because it does not explicitly prune empty space,
  *  unlike more advanced SAH-like techniques used for RT. On the other hand it leads to a shorter tree,
  *  faster to traverse and our experience shown that in the special case of kNN queries,
  *  this strategy is indeed more efficient (and much faster to build). Moreover, for volume data
  *  (e.g., fluid simulation) pruning the empty space is useless.
  *
  *  Actually, storing at each node the exact AABB (we therefore have a binary BVH) allows
  *  to prune only about 10% of the leaves, but the overhead of this pruning (ball/ABBB intersection)
  *  is more expensive than the gain it provides and the memory consumption is x4 higher !
  */
template<typename SearchTraits>
void PointKdtree_3<SearchTraits>::createTree(unsigned int nodeId, unsigned int start, unsigned int end, unsigned int level, unsigned int targetCellSize, unsigned int targetMaxDepth)
{
  Node& node = mNodes[nodeId];
  FT aabb_min[3], aabb_max[3];
  for (int k=0; k<3; ++k)
  {
    aabb_min[k] = aabb_max[k] = mElements[start][k];
  }
  for (unsigned int i=start+1 ; i<end ; ++i)
  {
    for (int k=0; k<3; ++k)
    {
      if (FT(mElements[i][k]) < aabb_min[k])
        aabb_min[k] = mElements[i][k];
      if (FT(mElements[i][k]) > aabb_max[k])
        aabb_max[k] = mElements[i][k];
    }
  }

  PointType diag(aabb_max[0] - aabb_min[0], aabb_max[1] - aabb_min[1], aabb_max[2] - aabb_min[2]);
  unsigned int dim = diag[0]>diag[1] ? (diag[0]>diag[2] ? 0 : 2) : (diag[1]>diag[2] ? 1 : 2);
  node.dim = dim;
  node.splitValue = FT(0.5*(aabb_max[dim] + aabb_min[dim]));

  unsigned int midId = split(start, end, dim, node.splitValue);

  node.firstChildId = mNodes.size();
  mNodes.resize(mNodes.size()+2);

  {
    // left child
    unsigned int childId = mNodes[nodeId].firstChildId;
    Node& child = mNodes[childId];
    if (midId-start <= targetCellSize || level>=targetMaxDepth)
    {
        child.leaf = 1;
        child.start = start;
        child.size = midId-start;
    }
    else
    {
        child.leaf = 0;
        createTree(childId, start, midId, level+1, targetCellSize, targetMaxDepth);
    }
  }

  {
    // right child
    unsigned int childId = mNodes[nodeId].firstChildId+1;
    Node& child = mNodes[childId];
    if (end-midId <= targetCellSize || level>=targetMaxDepth)
    {
      child.leaf = 1;
      child.start = midId;
      child.size = end-midId;
    }
    else
    {
      child.leaf = 0;
      createTree(childId, midId, end, level+1, targetCellSize, targetMaxDepth);
    }
  }
}

}

#endif
