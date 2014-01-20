#ifndef CGAL_EFFICIENT_RANSAC_PRIMITIVE_H
#define CGAL_EFFICIENT_RANSAC_PRIMITIVE_H

#include <vector>
#include <set>
#include <boost/tuple/tuple.hpp>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_3.h>

extern int ccTime;
extern int ccCount;

namespace CGAL {

  namespace Efficient_ransac {

    template <typename Kernel, class inputDataType>
    class Primitive_ab
    {
      template <typename K, class DT>
      friend class Efficient_ransac;
      template <typename K, class PA, class DT>
      friend class Octree;
    public:
      //--------------------------------------------typedef
      typedef typename std::vector<inputDataType>::iterator inputIterator;
      typedef typename std::vector<inputDataType>::const_iterator InputConstIterator;

      typedef typename Kernel::FT FT;

      // basic geometric types
      typedef typename Kernel::Point_3 Point;
      typedef typename Kernel::Vector_3 Vector;
      typedef typename Kernel::Vector_2 Vector2;
      typedef typename Kernel::Plane_3 Plane_3;
      
      //more advance container for geometric types
      typedef typename std::pair<Point, Vector> Pwn;

      typedef Primitive_ab<Kernel, inputDataType> Primitive;

      Primitive_ab() {init();}
      Primitive_ab(FT _ep, FT _t) {m_epsilon=_ep; m_normalThresh =_t; init();}
      static Primitive_ab* create(int id, FT _m_epsilon = 0.9f, FT _normalThresh = 0.9f);

      enum TYPE : int {PLANE = 0, SPHERE = 1, CYLINDER = 2, CONE = 3, TORUS = 4};

    public:
      FT m_epsilon;
      FT m_normalThresh;	 //deviation of normal, used during first check of the 3 normal

    protected:
      TYPE m_type;
      std::string m_type_name;
      bool m_isValid;
      FT m_lowerBound;
      FT m_upperBound;

      unsigned int m_score;

      FT m_sum_ExpectedValue;
      int m_nb_subset_used;		//count the number of subset used so far for the score, and thus indicate the next one to use
      bool m_hasConnectedComponent;

      std::vector<int> m_indices;	//indices of the points fitting to the candidate

    public:

      TYPE type() {return m_type;}
      FT ExpectedValue() const {return (m_lowerBound + m_upperBound) / 2.f;}
      FT inline minBound() const {return  m_lowerBound;}
      FT inline maxBound() const {return  m_upperBound;}
      FT inline score() const {return m_score;} //return last computed score, or -1 if no score yet
      int inline subsets() const {return m_nb_subset_used;}

      operator FT() const {return ExpectedValue();}//so we can sort by expected value

      void updatePoints(const std::vector<int> &shapeIndex);
      bool equals(InputConstIterator first, const Primitive *other) const {
        FT e2 = m_epsilon * m_epsilon;
        int fit = 0;
        int tested = 0;
        int toTest = std::min<unsigned int>(9, m_indices.size());
        for (unsigned int i = 0;i<toTest;i++) {
          int idx = getRandomInt() % m_indices.size();
          if (other->squared_distance(first[m_indices[idx]].first) <= e2)
            if (other->cos_to_normal(first[m_indices[idx]].first, first[m_indices[idx]].second) > m_normalThresh)
              fit++;
        }
        tested = toTest;
        toTest = std::min<unsigned int>(9, other->m_indices.size());

        for (unsigned int i = 0;i<toTest;i++) {
          int idx = getRandomInt() % other->m_indices.size();
          if (squared_distance(first[other->m_indices[idx]].first) <= e2)
            if (cos_to_normal(first[other->m_indices[idx]].first, first[other->m_indices[idx]].second) > m_normalThresh)
              fit++;
        }
        tested += toTest;

        return fit >= 2 * tested / 3;
      }

      bool isValid() const {return m_isValid;}
      std::vector<int> *getPointsIndices() {return &m_indices;}

      void save(const char* _n, const inputIterator _data/*, const bool _withAllPoints = false*/)
      {
        std::ofstream plyFile(_n);

        plyFile << "ply" << std::endl;
        plyFile << "format ascii 1.0" << std::endl;
        plyFile << "element vertex " << m_indices.size() << std::endl;
        plyFile << "property float x" << std::endl;
        plyFile << "property float y" << std::endl;
        plyFile << "property float z" << std::endl;
        plyFile << "property uchar red" << std::endl;
        plyFile << "property uchar green" << std::endl;
        plyFile << "property uchar blue" << std::endl;
        plyFile << "end_header" << std::endl;

        plyFile << std::setprecision(6);
        unsigned char r, g, b;
        r = 64 + getRandomInt()%192;
        g = 64 + getRandomInt()%192;
        b = 64 + getRandomInt()%192;

        for (unsigned int i = 0;i<m_indices.size();i++) {
          Point p = (*(_data+ m_indices[i])).first;
          plyFile << p[0] << " " << p[1] << " " << p[2];
          plyFile << " " << (int)r << " " << (int)g << " " << (int)b;
          plyFile << std::endl;
        }

        plyFile.close();
      }

      virtual FT squared_distance(const Point &_p) const = 0;
      virtual void squared_distance(InputConstIterator first, std::vector<FT> &dists, const std::vector<int> &shapeIndex, const std::vector<unsigned int> &indices) = 0;
      virtual FT cos_to_normal(const Point &_p, const Vector &_n) const = 0;
      virtual void cos_to_normal(InputConstIterator first, std::vector<FT> &angles, const std::vector<int> &shapeIndex, const std::vector<unsigned int> &indices) const = 0;

      virtual void parameterExtend(const Point &center, FT width, FT min[2], FT max[2]) const = 0;
      virtual void parameters(InputConstIterator first, std::vector<std::pair<FT, FT>> &parameterSpace, const std::vector<int> &indices, FT min[2], FT max[2]) const = 0;
      
      unsigned int connectedComponent(InputConstIterator first, FT m_bitmapEpsilon, const Point &center, FT width);

      virtual void compute(std::set<int> &l_list_index_selected, InputConstIterator &m_it_Point_Normal) = 0;
      virtual std::string info() = 0;
      virtual std::string type_str() const = 0;
      inline bool operator<(const Primitive &c) const {
        return ExpectedValue() < c.ExpectedValue();
      }
      
    protected:
      unsigned int cost_function(InputConstIterator &first, const std::vector<int> &shapeIndex, FT epsilon, FT normal_threshold, const std::vector<unsigned int> &indices);

      static bool default_cost_function (const Primitive_ab<Kernel, inputDataType> *_this, 
        const Point &_p, const Vector &_n, const FT epsilon_distance, const FT delta_dev_normal);

      template<typename T> bool isfinite(T arg)
      {
        return arg == arg && 
          arg != std::numeric_limits<T>::infinity() &&
          arg != -std::numeric_limits<T>::infinity();
      }

      void init()
      {
        m_isValid = true;
        m_hasConnectedComponent = false;
        m_score = 0;
        m_nb_subset_used = 0;
        m_sum_ExpectedValue = 0;
        m_lowerBound = std::numeric_limits<FT>::min();
        m_upperBound = std::numeric_limits<FT>::min();
      };

      void computeBound(const int sizeS1, const int sizeP) {
        hypergeometricalDist(-2 - sizeS1, -2 - sizeP, -1 - signed(m_indices.size()), m_lowerBound, m_upperBound);
        m_lowerBound = -1 - m_lowerBound;
        m_upperBound = -1 - m_upperBound;
      }

      void hypergeometricalDist(const int UN, const int x, const FT n, FT &low, FT &high) 
      {
        if (UN == 1 || UN == 0)
          printf("something wrong here, denominator is zero (UN %d)!! \n", UN);
        if (x > UN)
          printf("SizeP1 smaller than sizeP, something wrong (and sqrt may be negative !!!!");
        FT sq = sqrtf(x * n * (UN- x) * (UN - n) / (UN - 1));
        low  = (x * n - sq) / UN;
        high = (x * n + sq)/UN;

        if (!isfinite<FT>( low ) || !isfinite<FT>( high ))
        {
          low = high = 0;
        }
      }

      virtual bool supportsConnectedComponent() = 0;
      virtual bool wrapsU() const = 0;
      virtual bool wrapsV() const = 0;
    };


    template <typename Kernel, class inputDataType>
    Primitive_ab<Kernel, inputDataType>* Primitive_ab<Kernel, inputDataType>::create(int id, FT _epsilon, FT _normalThresh)
    {
      switch (id)
      {
      case CONE:
        return new Cone<Kernel, inputDataType>(_epsilon, _normalThresh);

      case TORUS:
        return new Torus<Kernel, inputDataType>(_epsilon, _normalThresh);

      case PLANE:
        return new Plane<Kernel, inputDataType>(_epsilon, _normalThresh);

      case SPHERE:
        return new Sphere<Kernel, inputDataType>(_epsilon, _normalThresh);

      case CYLINDER:
        return new Cylinder<Kernel, inputDataType>(_epsilon, _normalThresh);

      default:
        return NULL;
      }
    };

    template <typename Kernel, class inputDataType>
    unsigned int Primitive_ab<Kernel, inputDataType>::cost_function(InputConstIterator &first, const std::vector<int> &shapeIndex, FT epsilon, FT normal_threshold, const std::vector<unsigned int> &indices) {
      std::vector<FT> dists, angles;
      dists.resize(indices.size());
      squared_distance(first, dists, shapeIndex, indices);
      angles.resize(indices.size());
      cos_to_normal(first, angles, shapeIndex, indices);

      unsigned int scoreBefore = m_indices.size();

      FT eps = epsilon * epsilon;
      int i = 0;
      for (unsigned int i = 0;i<indices.size();i++) {
        if (shapeIndex[indices[i]] == -1) {
          if (dists[i] <= eps && angles[i] > normal_threshold)
            m_indices.push_back(indices[i]);
        }
      }

      return m_indices.size() - scoreBefore;
    }

    template <typename Kernel, class inputDataType>
    bool Primitive_ab<Kernel, inputDataType>::default_cost_function( const Primitive_ab<Kernel, inputDataType> *_this, 
      const Point &_p, const Vector &_n, const FT epsilon_distance, const FT delta_dev_normal)
    {
      if (_this->squared_distance(_p)  > epsilon_distance * epsilon_distance) return false;
      if (_this->cos_to_normal(_p, _n) >  delta_dev_normal ) return false;
      return true;
    }

    template <typename Kernel, class inputDataType>
    void Primitive_ab<Kernel, inputDataType>::updatePoints(const std::vector<int> &shapeIndex) {
      if (!m_indices.size())
        return;
      int start = 0, end = m_indices.size() - 1;
      while (start < end) {
        while (shapeIndex[m_indices[start]] == -1 && start < end) start++;
        while (shapeIndex[m_indices[end]] != -1 && start < end) end--;
        if (shapeIndex[m_indices[start]] != -1 && shapeIndex[m_indices[end]] == -1 && start < end) {
          unsigned int tmp = m_indices[start];
          m_indices[start] = m_indices[end];
          m_indices[end] = tmp;
        }
      }
      m_indices.resize(end);
      m_score = m_indices.size();
    }
    
    template <typename Kernel, class inputDataType>
    unsigned int Primitive_ab<Kernel, inputDataType>::connectedComponent(InputConstIterator first, FT m_bitmapEpsilon, const Point &center, FT width) {
      if (m_indices.size() == 0)
        return 0;

//       if (m_hasConnectedComponent)
//         return m_score;

      m_hasConnectedComponent = true;
      if (!supportsConnectedComponent())
        return m_indices.size();

      ccCount++;
      clock_t s, e;
      s = clock();

      FT min[2], max[2];
      //parameterExtend(center, width, min, max);
      std::vector<std::pair<FT, FT>> parameterSpace;
      parameterSpace.resize(m_indices.size());

      parameters(first, parameterSpace, m_indices, min, max);
      int iMin[2], iMax[2];
      iMin[0] = min[0] / m_bitmapEpsilon;
      iMin[1] = min[1] / m_bitmapEpsilon;
      iMax[0] = max[0] / m_bitmapEpsilon;
      iMax[1] = max[1] / m_bitmapEpsilon;

      int uExtend = abs(iMax[0] - iMin[0]) + 2;
      int vExtend = abs(iMax[1] - iMin[1]) + 2;

      std::vector<std::vector<int>> bitmap;
      std::vector<bool> visited;
      bitmap.resize(uExtend * vExtend);
      visited.resize(uExtend * vExtend, false);

      int maxIndex = uExtend * vExtend;

      bool wrapU = wrapsU();
      bool wrapV = wrapsV();

      for (unsigned int i = 0;i<parameterSpace.size();i++) {
        int u = (parameterSpace[i].first - min[0]) / m_bitmapEpsilon;
        int v = (parameterSpace[i].second - min[1]) / m_bitmapEpsilon;
        if (u < 0 || u >= uExtend) {
          if (wrapU) {
            while (u < 0) u += uExtend;
            while (u >= uExtend) u-= uExtend;
          }
          else {
            std::cout << "cc: u out of bounds: " << u << std::endl;
            u = (u < 0) ? 0 : (u >= uExtend) ? uExtend - 1 : u;
          }
        }
        if (v < 0 || v >= vExtend) {
          if (wrapV) {
            while (v < 0) v += vExtend;
            while (v >= vExtend) v-= vExtend;
          }
          else {
            std::cout << "cc: v out of bounds: " << u << std::endl;
            v = (v < 0) ? 0 : (v >= vExtend) ? vExtend - 1 : v;
          }
        }
        bitmap[v * uExtend + u].push_back(m_indices[i]);
      }

      std::vector<std::vector<int>> cluster;
      for (unsigned int i = 0;i<(uExtend * vExtend);i++) {
        cluster.push_back(std::vector<int>());
        if (bitmap[i].empty())
          continue;
        if (visited[i])
          continue;

        std::stack<int> fields;
        fields.push(i);
        while (!fields.empty()) {
          int f = fields.top();
          fields.pop();
          if (visited[f])
            continue;
          visited[f] = true;
          if (bitmap[f].empty())
            continue;

          // copy indices
          std::copy(bitmap[f].begin(), bitmap[f].end(), std::back_inserter(cluster.back()));

          // grow 8-neighborhood
          int vIndex = f / uExtend;
          int uIndex = f % uExtend;
          bool upperBorder = vIndex == 0;
          bool lowerBorder = vIndex == (vExtend - 1);
          bool leftBorder = uIndex == 0;
          bool rightBorder = uIndex == (uExtend - 1);

          int n;
          if (!upperBorder) {
            n = f - uExtend;
            if (!visited[n])
              fields.push(n);
          }
          else if (wrapV) {
            n = f + (vExtend - 1) * uExtend;
            if (!visited[n]) fields.push(n);
          }

          if (!leftBorder) {
            n = f - 1;
            if (!visited[n]) fields.push(n);
          }
          else if (wrapU) {
            n = f + uExtend - 1;
            if (!visited[n]) fields.push(n);
          }

          if (!lowerBorder) {
            n = f + uExtend;
            if (!visited[n]) fields.push(n);
          }
          else if (wrapV) {
            n = f - (vExtend - 1) * uExtend;
            if (!visited[n]) fields.push(n);
          }

          if (!rightBorder) {
            n = f + 1;
            if (!visited[n]) fields.push(n);
          }
          else if (wrapU) {
            n = f - uExtend + 1;
            if (!visited[n]) fields.push(n);
          }
        }
      }

      int maxCluster = 0;
      for (unsigned int i = 1;i<cluster.size();i++) {
        if (cluster[i].size() > cluster[maxCluster].size()) {
          maxCluster = i;
        }
      }

      m_indices = cluster[maxCluster];

      e = clock();
      ccTime += e - s;

      return m_score = m_indices.size();
    }
  }
}
#endif
