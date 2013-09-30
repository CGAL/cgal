#ifndef CGAL_EFFICIENT_RANSAC_PRIMITIVE_H
#define CGAL_EFFICIENT_RANSAC_PRIMITIVE_H

#include <vector>
#include <set>
#include <boost/tuple/tuple.hpp>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_3.h>

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
      //typedef typename type_primitive int;
      typedef typename std::vector<inputDataType>::iterator inputIterator;
      typedef typename std::vector<inputDataType>::const_iterator InputConstIterator;
      typedef std::vector<int>::const_iterator IntConstIterator;

      typedef typename Kernel::FT FT;

      // basic geometric types
      typedef typename Kernel::Point_3 Point;
      typedef typename Kernel::Point_2 Point_2d;
      typedef typename Kernel::Vector_3 Vector;
      typedef typename Kernel::Plane_3 Plane_3;
      typedef typename boost::tuple<Point, int> Point_and_int;

      typedef typename Kernel::Iso_cuboid_3 Iso_cuboid_3;
      typedef typename Kernel::Iso_rectangle_2 bbox_2d;

      //more advance container for geometric types
      typedef typename std::vector<Point>::const_iterator Point_const_iterator;
      typedef typename std::pair<Point, Vector> Pwn;
      typedef typename boost::tuple<Pwn, int> Pwn_and_int;
      typedef typename std::vector<Pwn> Pwn_vector;
      typedef typename Pwn_vector::iterator Pwn_iterator;
      typedef typename std::vector<Pwn_and_int>::iterator Pwn_and_int_iterator;

      //for kdTree
      typedef typename CGAL::Search_traits_3<Kernel> Traits_base;
      typedef typename CGAL::Search_traits_adapter<Point_and_int, My_point_property_map<Kernel>, Traits_base>    Traits;
      typedef typename CGAL::Kd_tree<Traits> Tree;
      typedef typename CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

      typedef Primitive_ab<Kernel, inputDataType> Primitive;

      typedef boost::function<Point(Pwn&)> getPointFunc;
      typedef boost::transform_iterator<getPointFunc, Pwn_iterator> point_iterator;

      typedef boost::function<Vector(Pwn&)> getNormalFunc;
      typedef boost::transform_iterator<getNormalFunc, Pwn_iterator> normal_iterator;

      Primitive_ab() {init();};
      Primitive_ab(FT _ep, FT _t) {m_epsilon=_ep; m_normalThresh =_t; init();};
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

      //std::vector<int> m_score_by_subset;	//save the result of the score function (i.e, the number of points that belong to the candidate) by subset
      FT m_sum_ExpectedValue;
      int m_nb_subset_used;		//count the number of subset used so far for the score, and thus indicate the next one to use

      std::vector<int> m_indices;	//indices of the points fitting to the candidate

    public:
      /*
      auto callback = [Primitive* p, ](int x)  -> int
      { 
      if (CGAL::squared_distance ( m_plane, *(m_point+_elem)) > m_epsilon*m_epsilon) return -1;

      normal_iterator it  = 	(m_normal +_elem);
      //if ( abs( *(it) * m_plane.orthogonal_vector () ) <  m_ndev*sqrtf(it->squared_length()*m_plane.orthogonal_vector ().squared_length()))  return -1;
      if ( (*(it) * m_plane.orthogonal_vector ())*( *(it) * m_plane.orthogonal_vector () ) <  m_ndev*m_ndev*(it->squared_length()*m_plane.orthogonal_vector ().squared_length()))  return -1;


      return _elem;
      };
      */

      TYPE type() {return m_type;}
      FT ExpectedValue() const {return (m_lowerBound + m_upperBound) / 2.f;}
      FT inline minBound() const {return  m_lowerBound;}
      FT inline maxBound() const {return  m_upperBound;}
      FT inline score() const {return m_score;} //return last computed score, or -1 if no score yet

      operator FT() const {return ExpectedValue();}//so we can sort by expected value

      void updatePoints(const std::vector<int> &shapeIndex);

      //bool hasCommonPoints(const std::vector<int> *indices_points_best_candidate);
      //void attachPoints(const point_iterator _it_Point, const normal_iterator _it_normal, const std::vector< std::vector< std::vector<voxel> > >& _v, const Iso_cuboid_3& _bb, const FT _size_c, const std::vector<bool> &_p_available) ;

      //basically add more subset to compute score
      /*bool inline improveBound(const std::pair<std::vector<Tree*>*, std::vector<int>> _t, const int _SizeP, const point_iterator _it_Points, const normal_iterator _it_Normal, const std::vector<bool> &available_points, unsigned int max_subset, unsigned int min_points)
      {
        return improveBound(_t, m_query_ball, _SizeP, _it_Points, _it_Normal, available_points, max_subset, min_points);
      }*/

      //bool improveBound(const std::pair<std::vector<Tree*> *, std::vector<int>> _t, const Fuzzy_sphere _query, const int _SizeP, const point_iterator _it_Points, const normal_iterator _it_Normal, const std::vector<bool> &is_available, unsigned int max_subset, unsigned int min_points);

      bool isValid() const {return m_isValid;}
      std::vector<int> *getPointsIndices() {return &m_indices;}

      void save(const char* _n, const inputIterator _data/*, const bool _withAllPoints = false*/)
      {
        int nbPoints = m_indices.size();
        std::ofstream f(_n);
        f << "OFF\n";
        f << nbPoints << " 0 0\n";
        for (int i=0;i<	nbPoints ;i++)
          f << 	(*(_data+ m_indices[i])).first << "\n";
      }

/*
      virtual Point projection() const = 0;
      virtual Point projection(const Point &_p) const = 0;*/
      virtual FT squared_distance(const Point &_p) const = 0;
      virtual void squared_distance(InputConstIterator first, std::vector<FT> &dists, const std::vector<int> &shapeIndex, const std::vector<unsigned int> &indices) = 0;
      virtual FT cos_to_normal(const Point &_p, const Vector &_n) const = 0;
      virtual void cos_to_normal(InputConstIterator first, std::vector<FT> &angles, const std::vector<int> &shapeIndex, const std::vector<unsigned int> &indices) const = 0;

      virtual void parameterExtend(const Point &center, FT width, FT min[2], FT max[2]) const = 0;
      virtual void parameters(InputConstIterator first, std::vector<std::pair<FT, FT>> &parameterSpace, const std::vector<int> &indices) const = 0;
      
      unsigned int connectedComponent(InputConstIterator first, FT m_bitmapEpsilon, const Point &center, FT width);

      virtual void compute(std::set<int> &l_list_index_selected, Pwn_iterator &m_it_Point_Normal) = 0;

      //used only during region growing (at the end to attach the points)
      //virtual int score(const point_iterator _it_Points, const normal_iterator _it_Normals, std::vector<int> &index_point_to_process, const std::vector<bool> &is_available, const FT coeff_ep = 1)	= 0;
      //int score(const point_iterator _it_Points, const normal_iterator _it_Normals, std::vector<int> &index_point_to_process, const std::vector<bool> &is_available, const FT coeff_ep = 1);

      //used by computeBound
      //virtual int score(const std::vector<Tree*> *_t, const Fuzzy_sphere _query, const point_iterator _it_Points, const normal_iterator _it_Normals, const std::vector<bool> &is_available) = 0;
      //int score(const std::vector<Tree*> *_t, const Fuzzy_sphere _query, const point_iterator _it_Points, const normal_iterator _it_Normals, const std::vector<bool> &is_available, const FT coeff_ep = 1);

      //int score(Pwn_iterator begin, Pwn_iterator end, unsigned int offset, const std::vector<int> &shapeIndex);
      
      //virtual int score(FT gam, std::vector<Point_and_int>* l_points, normal_iterator l_it_Normal)  = 0;
      virtual std::string info() = 0;
      virtual std::string type_str() const = 0;
      
    private:
      // int (*pt2Func)(FT, char, char)
      //typedef bool (*type_cost_function)(const Primitive_ab<Kernel, inputDataType> *_this, 
      //const Point &_p, const Vector &_n, const FT epsilon_distance, const FT delta_dev_normal);

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
          printf("UN%d  S1%d    x%d  P%d  n%f   score%f             %f\n", UN, -2 - UN, x, -2 - x, n, -1 - n, x * n * (UN - x) * (UN - n) / (UN - 1));
          system("pause");
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
    /*
    template <typename Kernel, class inputDataType>
    bool Primitive_ab<Kernel, inputDataType>::hasCommonPoints(const std::vector<int> *indices_points_best_candidate)
    {
      bool commonPoint = false;

#pragma omp parallel for
      for (int i = 0;i < indices_points_best_candidate->size();++i)
      {
#pragma omp flush (commonPoint)	    //abort at the first common index met
        if (!commonPoint) 
        {			 
          for (int j = 0;j < m_indices.size();++j)
          {
            if (indices_points_best_candidate->at(i) == m_indices[j])
            {
              commonPoint = true;
              break;
            }
          }
        }
      }
      return commonPoint;
    }
    /*
    template <typename Kernel, class inputDataType>
    void Primitive_ab<Kernel, inputDataType>::attachPoints(const point_iterator _it_Point, const normal_iterator _it_Normal, const std::vector<std::vector<std::vector<voxel>>>& _v, const Iso_cuboid_3& _bb, const FT _size_c, const std::vector<bool> &_p_available)
    {
      if (m_indices.size() != 0) return;//point already attached
      //1 take a point on the surface and propagate
      Point l_start_p = projection();

      int l_w = int(ceil( (_bb.xmax () - _bb.xmin ()) / _size_c)) ; 
      int l_h = int(ceil( (_bb.ymax () - _bb.ymin ()) / _size_c)) ;
      int l_d = int(ceil( (_bb.zmax () - _bb.zmin ()) / _size_c)) ;

      int id_voxel_start_x = int((l_start_p.x() - _bb.xmin ()) / _size_c) ;
      int id_voxel_start_y = int((l_start_p.y() - _bb.ymin ()) / _size_c) ;
      int id_voxel_start_z = int((l_start_p.z() - _bb.zmin ()) / _size_c) ;
      Point voxel_start(id_voxel_start_x, id_voxel_start_y, id_voxel_start_z);

      std::vector <Point> l_result;
      l_result.push_back(voxel_start);
      std::vector<int> l_index_points;

      std::vector<Point> l_voxel_added;
      l_voxel_added.push_back(voxel_start);

      std::vector<std::vector< std::vector<bool> > > is_free(l_w, std::vector< std::vector<bool> >(l_h, std::vector<bool>(l_d, true)));
      is_free[(l_voxel_added[0].x())][(l_voxel_added[0].y())][(l_voxel_added[0].z())]=false;

      bool new_voxel_to_propagate = true;

      FT x, x2, y, y2, z, z2;
      int id_current_voxel;
      do
      {
        std::vector < Point > voxel_to_add;
        new_voxel_to_propagate=false;
        for (int ii = 0;ii < l_voxel_added.size();ii++)
        {
          {
            //look at neighborg
            int i, j, k;
            //todo, run me in parallel !

            //std::vector<bool> add_me(3*3*3, false);
            //#pragma omp parallel for
            for (int jj=0;jj<3*3*3 ;jj++)
            {
              i = int(jj / 9);		//(9x0, 9x1, 9x2)
              j = (jj - i*9) / 3;		//(3x0, 3x1, 3x2, 3x0, ...)
              k = jj % 3;				//(0, 1, 2, 0, ...)
              i = i - 1; j= j - 1; k= k - 1;

              if (i ==0 && j == 0 && k == 0) continue;   //the current one is skipped

              x2 = (l_voxel_added[ii].x()+k)*_size_c  + _bb.xmin (); 
              y2 = (l_voxel_added[ii].y()+j)*_size_c  + _bb.ymin (); 
              z2 = (l_voxel_added[ii].z()+i)*_size_c  + _bb.zmin (); 

              if (x2 >= _bb.xmin () && x2 <= _bb.xmax () && y2 >= _bb.ymin () && y2 <= _bb.ymax () && z2 >= _bb.zmin () && z2 <= _bb.zmax ())
              {
                Point l_center_p(x2, y2, z2);

                //the voxel somehow intersects the primitive
                if (is_free[(l_voxel_added[ii].x()+k)][(l_voxel_added[ii].y()+j)][(l_voxel_added[ii].z()+i)] )							//if (_v[id_current_voxel].size() > 0)
                  if (CGAL::squared_distance ( l_center_p, projection(l_center_p))*2 <= _size_c*_size_c) 	//d < sqrt(2)*c/2
                  {
                    //add_me[jj]=true;
                    l_center_p = Point((l_voxel_added[ii].x()+k), (l_voxel_added[ii].y()+j), (l_voxel_added[ii].z()+i));
                    voxel_to_add.push_back(l_center_p);
                    l_result.push_back(l_center_p);
                    new_voxel_to_propagate=true;
                    is_free[(l_voxel_added[ii].x()+k)][(l_voxel_added[ii].y()+j)][(l_voxel_added[ii].z()+i)] = false;

                    //get the points' index
                    for (int kk=0;kk<_v[(l_voxel_added[ii].x()+k)][(l_voxel_added[ii].y()+j)][(l_voxel_added[ii].z()+i)].size();kk++)
                      l_index_points.push_back(_v[(l_voxel_added[ii].x()+k)][(l_voxel_added[ii].y()+j)][(l_voxel_added[ii].z()+i)][kk]);
                  }
              }
            }


            //for (int jj=0;jj<3*3*3 ;jj++)
            //	if (add_me[jj]) 
            //	{
            //		i = int(jj / 9);		//(9x0, 9x1, 9x2)
            //		j = (jj - i*9) / 3;		//(3x0, 3x1, 3x2, 3x0, ...)
            //		k = jj % 3;				//(0, 1, 2, 0, ...)
            //		i = i - 1; j= j - 1; k= k - 1;

            //		Point l_center_p = Point((l_voxel_added[ii].x()+k), (l_voxel_added[ii].y()+j), (l_voxel_added[ii].z()+i));
            //		voxel_to_add.push_back(l_center_p);
            //		l_result.push_back(l_center_p);
            //		new_voxel_to_propagate=true;
            //		is_free[(l_voxel_added[ii].x()+k)][(l_voxel_added[ii].y()+j)][(l_voxel_added[ii].z()+i)] = false;
            //	}
          }
        }

        l_voxel_added = voxel_to_add;

      } while(new_voxel_to_propagate);


      //printf("end propagate %d\n", l_result.size());

      //now project the points from those voxels on the primitive and save the index of the matching points !
      //REMOVED BECAUSE OF SCOREint result = score(_it_Point, _it_Normal, l_index_points, _p_available, 3);
      int result = 10;
      m_indices = std::vector<int>(result);//score is the number of point != -1   
      std::vector<int>::iterator new_end = std::remove_if(l_index_points.begin(), l_index_points.end(), std::bind2nd(std::equal_to<int>(), -1));
      std::copy(l_index_points.begin(), new_end, m_indices.begin());

    }

    /*
    template <typename Kernel, class inputDataType>
    class cost_RANSAC
    {
    private:
    Primitive_ab<Kernel, inputDataType>	_primitive;
    point_iterator _it_Points;
    normal_iterator _it_Normals;
    std::vector<bool> _is_available;
    }
    */

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
      if (!supportsConnectedComponent())
        return m_indices.size();

      FT min[2], max[2];
      parameterExtend(center, width, min, max);
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

      std::vector<std::pair<FT, FT>> parameterSpace;
      parameterSpace.resize(m_indices.size());

      parameters(first, parameterSpace, m_indices);

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

      return m_indices.size();
    }
    /*
    template <typename Kernel, class inputDataType>
    int Primitive_ab<Kernel, inputDataType>::score(const point_iterator _it_Points, const normal_iterator _it_Normals, std::vector<int> &index_point_to_process, const std::vector<bool> &_is_available, const FT coeff_ep)
    {
      const Primitive_ab<Kernel, inputDataType> *_this = this;
      FT epsilon_distance = coeff_ep*m_epsilon;
      FT delta_dev_normal = m_normalThresh;

      auto callback =  [_this, &_it_Points, &_it_Normals, &_is_available, &epsilon_distance, &delta_dev_normal](const int& _elem)  -> int 							
      {
        if (!_is_available[_elem]) return -1;
        if (!Primitive_ab<Kernel, inputDataType>::default_cost_function(_this, *(_it_Points + _elem), *(_it_Normals + _elem), epsilon_distance, delta_dev_normal)) return -1;
        //if (!Primitive_ab<Kernel, inputDataType>::_cost_function(_this, *(_it_Points+_elem), *(_it_Normals +_elem), epsilon_distance, delta_dev_normal)) return -1;
        return _elem;
      };

      std::transform (index_point_to_process.begin(), index_point_to_process.end(), index_point_to_process.begin(), callback);

      //for now we skip the cc part, and count the nb point not equal to -1
      int l_score =  std:: count_if(index_point_to_process.begin(), 
        index_point_to_process.end(), 
        bind2nd(std::not_equal_to<int>(), -1));
      return l_score;
    }	

    template <typename Kernel, class inputDataType>
    int Primitive_ab<Kernel, inputDataType>::score(const std::vector<Tree*> *_t, const Fuzzy_sphere _query, const point_iterator _it_Points, const normal_iterator _it_Normals, const std::vector<bool> &_is_available, const FT coeff_ep)
    {
      //get the points as Point_and_int
      std::vector<Point_and_int> l_points_s1;
      _t->at(m_nb_subset_used)->search(std::back_inserter(l_points_s1), _query);

      //get the int only
      // See "Binding an overloaded function" in boost::bind documentation
      int const& (Point_and_int::*get_int_from_PaI) () const = &Point_and_int::get<1>;
      auto l_it_int_start = boost::make_transform_iterator(l_points_s1.begin(), boost::bind( get_int_from_PaI, _1));
      auto l_it_int_end = boost::make_transform_iterator(l_points_s1.end(), boost::bind( get_int_from_PaI, _1));


      std::vector<int> index_point_to_process(l_points_s1.size(), 0);
      std::copy(l_it_int_start, l_it_int_end, index_point_to_process.begin());

      return score(_it_Points, _it_Normals, index_point_to_process, _is_available, coeff_ep);
    }

    template <typename Kernel, class inputDataType>
    int Primitive_ab<Kernel, inputDataType>::score(Pwn_iterator begin, Pwn_iterator end, unsigned int offset, const std::vector<int> &shapeIndex) {
      unsigned int sum = 0;

      for (int i = 0;i<=(end-begin);i++) {
        if (shapeIndex[offset + i] != -1)
          continue;

        if (default_cost_function(this, begin[i].first, begin[i].second, m_epsilon, m_normalThresh)) {
          sum++;
          m_indices.push_back(offset + i);
        }
      }

      return sum;
    }
    */
    //basically add more subset to compute score
    //first parameter:
    //first the subset
    //second the number of used point in this subset
    /*
    template <typename Kernel, class inputDataType>
    bool Primitive_ab<Kernel, inputDataType>::improveBound(const std::pair<std::vector<Tree*>*, std::vector<int>> _t, const Fuzzy_sphere _query, const int _SizeP, const point_iterator _it_Points, const normal_iterator _it_Normal, const std::vector<bool> &is_available, unsigned int max_subset, unsigned int min_points)
    {
      //std::cout << "." << std::flush;
      //std::string inf = info();
      //if (_t.first->size() == m_nb_subset_used) {// YANNICK
      if (m_nb_subset_used >= max_subset || _t.first->size() <= m_nb_subset_used) {
        //printd("no more subset available\n");
        //std::cout << ";" << std::endl;
        return false;
      };


      //what it does is add another subset and recompute lower and upper bound
      //the next subset to include is provided by m_nb_subset_used

      //1. sum nb points of previous subset, 	  and score
      int l_sum_score = 0;
      int l_nb_total_points_subsets = 0;
      for (int i=0;i<m_nb_subset_used;i++)
      {
        l_nb_total_points_subsets += _t.second[i];	   //no some points are forbidden now
        l_sum_score += m_score_by_subset[i] ;
      }

      //2. need score of new subset as well as sum of the score of the previous concidered subset
      int l_new_score = 0;

      do
      {
        l_new_score = score(_t.first, _query, _it_Points, _it_Normal, is_available);
        //int tmp = score(octree)
        //printf("subset %d exp value %f\n", m_nb_subset_used, l_new_score);
        m_score_by_subset.push_back(l_new_score);	//add to vector of score

        l_nb_total_points_subsets+= _t.second[m_nb_subset_used]; //add point new subset
        l_sum_score += l_new_score ;

        m_nb_subset_used++;
      }
      while (l_sum_score < min_points && m_nb_subset_used < _t.first->size());
      // while (l_new_score == 0 && m_nb_subset_used < _t.first->size());YANNICK	   

      if (l_new_score == 0) {
        //std::cout << "c" << std::endl;
        //printd(info().c_str());
        return false;
      };

      //printf("->>>>>  %d %d\n" , l_nb_total_points_subsets, _SizeP);
      computeBound(l_nb_total_points_subsets, _SizeP, l_sum_score);//estimate the bound

      //m_nb_subset_used++;
      //std::cout << info() << " d" << std::endl;
      return true;
    }*/
  }
}
#endif
