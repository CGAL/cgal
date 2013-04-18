// Copyright (c) 2013 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Yannick Verdié, Clément Jamin
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************

#ifndef CGAL_EFFICIENT_RANSAC_H
#define CGAL_EFFICIENT_RANSAC_H

//for octree ------------------------------
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_3.h>
//#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/basic.h>
#include <CGAL/Search_traits_adapter.h>
#include <boost/iterator/zip_iterator.hpp>
#include <CGAL/bounding_box.h> //----------

#include <utility> // defines std::pair
#include <set>	    // defines std::set
#include <list>
#include <vector>
#include <fstream>
#include <random>
#define  _USE_MATH_DEFINES
#include <math.h>


//boost --------------
#include <boost/iterator/counting_iterator.hpp>
#include <boost/bind/make_adaptable.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <functional>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
//---------------------


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL {

namespace Efficient_ransac
{
#ifdef DEBUG
    #define printd(...) printf(__VA_ARGS__)
#else
    #define printd(...)
#endif


class Option_RANSAC
{
public:
    Option_RANSAC(): m_probability(0.001f),    m_minNbPoints(10), m_epsilon(0.01f), m_normalThresh(0.90f), m_bitmapEpsilon(0.01f)   {};

    //variable
public:
  float m_probability;//parameter for stopping generating candidate
  int m_minNbPoints;//parameter for minimum nb of points for a candidate
  float m_epsilon;  //distance gamma-band around the primitive
  float m_normalThresh;	  // this is the cos of the maximal normal deviation
  float m_bitmapEpsilon;	    //for cc, not used now
};




  //definition of the property map
template <typename Kernel>
class My_point_property_map
{
  typedef typename  Kernel::Point_3 Point;
  typedef typename  boost::tuple<Point,int> Point_and_int;

  typedef typename Point value_type;
  typedef const value_type& reference;
  typedef const boost::tuple<Point,int>& key_type;
  typedef boost::readable_property_map_tag category;

  //get function for the property map
  friend reference get(My_point_property_map,key_type p){return boost::get<0>(p);}
};

template <typename T,int i>
class accessTupleElement
{ 
public: 
  typedef typename boost::tuples::element<i,T>::type result_type;

  template <typename U> result_type operator()(const U& u) const
  { 
  // Select the ith element of the tuple T and return it 
  return boost::get<i>(u); 
  } 
};

typedef std::vector<int> voxel ;

template <typename Kernel, class inputDataType>
class Primitive_ab
{	 
  public:
//--------------------------------------------typedef
  //typedef typename type_primitive int;
  typedef typename std::vector<inputDataType>::iterator inputIterator;

  typedef typename Kernel::FT FT;
    
  // basic geometric types
  typedef typename  Kernel::Point_3 Point;
  typedef typename  Kernel::Point_2 Point_2d;
  typedef typename  Kernel::Vector_3 Vector;
  typedef typename  Kernel::Plane_3 Plane_3;
  typedef typename  boost::tuple<Point,int> Point_and_int;

  typedef typename  Kernel::Iso_cuboid_3 Iso_cuboid_3;
  typedef typename  Kernel::Iso_rectangle_2 bbox_2d;

  //more advance container for geometric types
  typedef typename std::vector<Point>::const_iterator Point_const_iterator;
  typedef typename std::pair<Point, Vector> Pwn;
  typedef typename  boost::tuple<Pwn,int> Pwn_and_int;
  typedef typename std::vector<Pwn> Pwn_vector;
  typedef typename Pwn_vector::iterator Pwn_iterator;
  typedef typename std::vector<Pwn_and_int>::iterator Pwn_and_int_iterator;

  //for kdTree
  typedef typename CGAL::Search_traits_3<Kernel> Traits_base;
  typedef typename CGAL::Search_traits_adapter<Point_and_int,My_point_property_map<Kernel>,Traits_base>    Traits;
  typedef typename CGAL::Kd_tree<Traits> Tree;
  typedef typename CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

  typedef Primitive_ab<Kernel,inputDataType> Primitive;

  typedef boost::function<Point(Pwn&)> getPointFunc;
  typedef boost::transform_iterator<getPointFunc, Pwn_iterator> point_iterator;


  typedef boost::function<Vector(Pwn&)> getNormalFunc;
  typedef boost::transform_iterator<getNormalFunc, Pwn_iterator> normal_iterator;

  Primitive_ab(){init();};
  Primitive_ab(float _ep, float _t) {m_epsilon=_ep;m_normalThresh =_t;init();};
  static Primitive_ab* create(int id, float _m_epsilon = 0.9f, float _normalThresh = 0.9f);

  enum TYPE : int {PLANE = 0/*, SPHERE = 1*/};//SPHERE = 1 and so on
  public:
    float m_epsilon;
    float m_normalThresh;	 //deviation of normal, used during first check of the 3 normal

  protected:
    TYPE m_type;
    std::string m_type_name;
    bool m_isValid;
    float m_lowerBound;
    float m_upperBound;

    std::vector<int> m_score_by_subset;	//save the result of the score function (i.e, the number of points that belong to the candidate) by subset
    float m_sum_ExpectedValue;
    int m_nb_subset_used;		//count the number of subset used so far for te score, and thus indicate the next one to use

    Fuzzy_sphere m_query_ball; //the query ball use to create this candidate
    std::vector<int> m_Indices;	  //indices of the points attached !
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




    TYPE type(){ return m_type;}; 
    float ExpectedValue() const { return (m_lowerBound + m_upperBound) / 2.f; }
    float inline minBound() const {return  m_lowerBound;};
    float inline maxBound() const {return  m_upperBound;};
    float inline getLastScore() const {if (m_score_by_subset.size() == 0) return -1; return m_score_by_subset.back();};//return last computed score, or -1 if no score yet

    operator float() const{return  ExpectedValue();}//so we can sort by expected value

    bool hasCommonPoints(const std::vector<int> *indices_points_best_candidate);
    void attachPoints(const point_iterator _it_Point, const normal_iterator _it_normal, const std::vector< std::vector< std::vector<voxel> > >& _v, const Iso_cuboid_3& _bb, const float _size_c, const std::vector<bool> &_p_available) ;

    //basically add more subset to compute score
    bool inline improveBound(const std::pair<std::vector<Tree*>* ,std::vector<int> > _t,const int _SizeP,const point_iterator _it_Points, const normal_iterator _it_Normal, const std::vector<bool> &available_points)
    {
      return improveBound(_t,  m_query_ball, _SizeP, _it_Points, _it_Normal, available_points);
    };	

    void save(const char* _n, const inputIterator _data/*, const bool _withAllPoints = false*/)
    {
      int nbPoints = m_Indices.size();
      std::ofstream f(_n);
      f<<"OFF\n";
      f<<nbPoints<<" 0 0\n";
      for (int i=0;i<	nbPoints ;i++)
        f<<	(*(_data+ m_Indices[i])).first<<"\n";
    }

    bool improveBound(const std::pair<std::vector<Tree*>* ,std::vector<int> > _t, const Fuzzy_sphere _query, const int _SizeP,const point_iterator _it_Points, const normal_iterator _it_Normal, const std::vector<bool> &is_available)	;

    bool isValid() const {return m_isValid;};
    std::vector<int> *getPointsIndices(){return &m_Indices;};

    virtual inline Point projection() const = 0;
    virtual inline Point projection(Point _p) const = 0;
    virtual void compute(std::set<int> &l_list_index_selected, const Fuzzy_sphere _query, Pwn_iterator &m_it_Point_Normal) = 0;

    //used only during region growing (at the end to attach the points)
    virtual int score(const point_iterator _it_Points, const normal_iterator _it_Normals, std::vector<int> &index_point_to_process, const std::vector<bool> &is_available, const float coeff_ep = 1)	= 0;
    //used by  computeBound
    virtual int score(const std::vector<Tree*> *_t, const Fuzzy_sphere _query,const point_iterator _it_Points, const normal_iterator _it_Normals, const std::vector<bool> &is_available) = 0;
      
    //virtual int score(float gam,  std::vector<Point_and_int>* l_points, normal_iterator l_it_Normal)  = 0;
    virtual void info() = 0;
    virtual std::string type_str() const = 0;

private:
    
    template<typename T> bool isfinite(T arg)
    {
      return arg == arg && 
            arg != std::numeric_limits<T>::infinity() &&
            arg != -std::numeric_limits<T>::infinity();
    }

      void inline init(){m_isValid = true;m_nb_subset_used = 0;m_sum_ExpectedValue = 0;m_lowerBound = std::numeric_limits<float>::min();m_upperBound = std::numeric_limits<float>::min();};

    void computeBound(const int sizeS1,const int sizeP,const float scoreS1){ hypergeometricalDist(-2-sizeS1, -2-sizeP, -1-scoreS1, m_lowerBound, m_upperBound);	m_lowerBound = -1-m_lowerBound;	m_upperBound = -1-m_upperBound;};
    void hypergeometricalDist(const int UN,const int x,const float n, float &low, float &high) 
    {
      if (UN == 1 || UN == 0) printf("something wrong here, denomminator is zero (UN %d)!! \n",UN);
      if (x > UN)printf("SizeP1 smaller than sizeP, something nwrong (and sqrt may be negative !!!!");
      float sq = sqrtf(x*n*(UN-x)*(UN-n)/(UN-1));
      low  = (x*n-sq)/UN;
      high = (x*n+sq)/UN;

      if (!isfinite<float>( low ) || !isfinite<float>( high ))
      {
        printf("UN%d  S1%d    x%d  P%d  n%f   score%f             %f\n",UN, -2-UN, x, -2-x, n,-1-n, x*n*(UN-x)*(UN-n)/(UN-1));system("pause");
      }
    };
};


template <typename Kernel, class inputDataType>
class valid_points
{
    //--------------------------------------------typedef
  typedef typename Kernel::FT FT;
    
  // basic geometric types
  typedef typename  Kernel::Point_3 Point;
  typedef typename  Kernel::Point_2 Point_2d;
  typedef typename  Kernel::Vector_3 Vector;
  typedef typename  Kernel::Plane_3 Plane_3;
  typedef typename  boost::tuple<Point,int> Point_and_int;

  typedef typename  Kernel::Iso_cuboid_3 Iso_cuboid_3;
  typedef typename  Kernel::Iso_rectangle_2 bbox_2d;

  //more advance container for geometric types
  typedef typename std::vector<Point>::const_iterator Point_const_iterator;
  typedef typename std::pair<Point, Vector> Pwn;
  typedef typename  boost::tuple<Pwn,int> Pwn_and_int;
  typedef typename std::vector<Pwn> Pwn_vector;
  typedef typename Pwn_vector::iterator Pwn_iterator;
  typedef typename std::vector<Pwn_and_int>::iterator Pwn_and_int_iterator;

  //for kdTree
  typedef typename CGAL::Search_traits_3<Kernel> Traits_base;
  typedef typename  CGAL::Search_traits_adapter<Point_and_int,My_point_property_map<Kernel>,Traits_base>    Traits;
  typedef typename CGAL::Kd_tree<Traits> Tree;
  typedef typename  CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

  typedef Primitive_ab<Kernel,inputDataType> Primitive;

  typedef boost::function<Point(Pwn&)> getPointFunc;
  typedef boost::transform_iterator<getPointFunc, Pwn_iterator> point_iterator;

  typedef boost::function<Vector(Pwn&)> getNormalFunc;
  typedef boost::transform_iterator<getNormalFunc, Pwn_iterator> normal_iterator;



  private:
    point_iterator m_point;
    normal_iterator  m_normal;
    const float m_epsilon;
    const float  m_ndev;
    const Plane_3 m_plane;

  public:	 
    //template <typename T> T operator ( ) (T& _elem ){return _elem;} 
      
    // Constructor initializes the value to multiply by
      
    valid_points  ( const Plane_3 _p, const float _gam, const float _al, const point_iterator _point, const normal_iterator _Val) : m_point(_point), m_normal ( _Val ), m_epsilon(_gam), m_plane(_p), m_ndev(_al){}
      

    int operator ( ) (int& _elem ) 
    {
      if (CGAL::squared_distance ( m_plane, *(m_point+_elem)) > m_epsilon*m_epsilon) return -1;

      normal_iterator it  = 	(m_normal +_elem);
      //if ( abs( *(it) * m_plane.orthogonal_vector () ) <  m_ndev*sqrtf(it->squared_length()*m_plane.orthogonal_vector ().squared_length()))  return -1;
      if ( (*(it) * m_plane.orthogonal_vector ())*( *(it) * m_plane.orthogonal_vector () ) <  m_ndev*m_ndev*(it->squared_length()*m_plane.orthogonal_vector ().squared_length()))  return -1;
        
        
      return _elem;
    } ;	  	
};

template <typename Kernel, class inputDataType>
class Plane : public Primitive_ab<Kernel,inputDataType>
{
  private:
    Plane_3	m_plane;
    Point m_point_on_primitive;

  protected:

  public:
    Plane() :  Primitive_ab<Kernel,inputDataType>(0.1f, 0.9f) {m_type = Primitive_ab::TYPE::PLANE; m_type_name ="Plane";}
    Plane(float _a,float _b)  :  Primitive_ab<Kernel,inputDataType>(_a,_b)  {m_type = Primitive_ab::TYPE::PLANE;m_type_name ="Plane";}
    void compute(std::set<int> &l_list_index_selected, const Fuzzy_sphere _query, Pwn_iterator &m_it_Point_Normal) 
    {
        m_isValid = true;
        if ( l_list_index_selected.size() < 3){	  m_isValid = false; return;};
      m_query_ball = _query;

        std::vector<int> output(l_list_index_selected.begin(), l_list_index_selected.end()); 

        m_plane = Plane_3(  (m_it_Point_Normal + output[0]   )->first, 
                  (m_it_Point_Normal + output[1]   )->first,
                  (m_it_Point_Normal + output[2]   )->first
                  );

        m_point_on_primitive =  (m_it_Point_Normal + output[0]   )->first;

        //check deviation of the 3 normal
      Vector l_v;
      //float cosTheta;
      for (int i=0;i<3;i++)
      {
        l_v = (m_it_Point_Normal + output[i]   )->second;
        //cosTheta= abs( l_v * m_plane.orthogonal_vector () )/sqrtf(l_v.squared_length()*m_plane.orthogonal_vector ().squared_length());
        //printf("val %d %f\n",i,cosTheta); 
        //if ( cosTheta < m_normalThresh)  m_isValid = false;
        if ( abs( l_v * m_plane.orthogonal_vector () ) < m_normalThresh*sqrtf(l_v.squared_length()*m_plane.orthogonal_vector ().squared_length()))  m_isValid = false;
      }
    }
    operator Plane_3 () { return m_plane;};
    void info()
    {
      std::cout<<"Type: "<<m_type_name<<std::endl;
      std::cout <<"Equation: "<< m_plane.a()<<"x + "<<m_plane.b()<<"b + "<<m_plane.c()<<"c + "<<m_plane.d()<<"= 0"<<std::endl;
      std::cout <<"Expected Value: "<<ExpectedValue()<<" with "<<m_nb_subset_used<<" subsets"<< std::endl;
      std::cout <<"Number of points attached (so far): "<<  m_Indices.size()	<< std::endl;
      std::cout<<"---------------------"<< std::endl;
    };
    std::string type_str() const{return m_type_name;};

    inline Point  pointOnPrimitive()const{/*return m_plane.point();*/return m_point_on_primitive;};
    inline Point projection() const {return projection(pointOnPrimitive()) ;}
    Point projection(Point _p) const {return m_plane.projection (_p);}
    
    int score(const point_iterator _it_Points, const normal_iterator _it_Normals, std::vector<int> &index_point_to_process, const std::vector<bool> &_is_available, const float coeff_ep)
    {

      // valid_points<Kernel,inputDataType> m_pred(m_plane, m_epsilon, m_normalThresh,_it_Points,_it_Normals); 
         
        auto callback =  [this,&_it_Points, &_it_Normals, &_is_available, &coeff_ep](const int& _elem)  -> int 							
      {
        if (!_is_available[_elem]) return -1;
        if (CGAL::squared_distance ( m_plane, *(_it_Points+_elem)) > coeff_ep*coeff_ep*m_epsilon*m_epsilon) return -1;
        normal_iterator it  = 	(_it_Normals +_elem);
        if ( (*(it) * m_plane.orthogonal_vector ())*( *(it) * m_plane.orthogonal_vector () ) <   m_normalThresh*m_normalThresh*(it->squared_length()*m_plane.orthogonal_vector ().squared_length()))  return -1;
        return _elem;
        };


        std::transform (index_point_to_process.begin(), index_point_to_process.end(), index_point_to_process.begin(), callback);

          
      //for now we skip the cc part, and count the nb point not equal to -1
        int l_score =  std:: count_if(index_point_to_process.begin(),
                          index_point_to_process.end(), 
                        bind2nd(std::not_equal_to<int>(),-1));
        return l_score;
         
        return 0;
    }	
      

    int score(const std::vector<Tree*> *_t, const Fuzzy_sphere _query, const point_iterator _it_Points, const normal_iterator _it_Normals, const std::vector<bool> &is_available)
    {
        //get the points
        std::vector<Point_and_int> l_points_s1;
        _t->at(m_nb_subset_used)->search(std::back_inserter(l_points_s1),_query);

        // See "Binding an overloaded function" in boost::bind documentation
        int const& (Point_and_int::*get_int_from_PaI) () const = &Point_and_int::get<1>;
        auto l_it_int_start = boost::make_transform_iterator(l_points_s1.begin(),boost::bind( get_int_from_PaI, _1));
        auto l_it_int_end = boost::make_transform_iterator(l_points_s1.end(),boost::bind( get_int_from_PaI, _1));


        //valid_points<Kernel,inputDataType> m_pred(m_plane, m_epsilon, m_normalThresh,_it_Points,_it_Normals); 
         
        std::vector<int> index_point_to_process(l_points_s1.size(),0);
        //std::transform (l_it_int_start, l_it_int_end, index_point_to_process.begin(), m_pred);
         
        auto callback =  [this,&_it_Points, &_it_Normals, &is_available](const int& _elem)  -> int 							
              {
                if (!is_available[_elem]) return -1;
                if (CGAL::squared_distance ( m_plane, *(_it_Points+_elem)) > m_epsilon*m_epsilon) return -1;
                normal_iterator it  = 	(_it_Normals +_elem);
                if ( (*(it) * m_plane.orthogonal_vector ())*( *(it) * m_plane.orthogonal_vector () ) <   m_normalThresh*m_normalThresh*(it->squared_length()*m_plane.orthogonal_vector ().squared_length()))  return -1;
                return _elem;
                };

        std::transform (l_it_int_start, l_it_int_end, index_point_to_process.begin(),callback);
        //parallel_transform(l_it_int_start, l_it_int_end, index_point_to_process.begin(),callback);		ppl only on vs !

      //for now we skip the cc part, and count the nb point not equal to -1
        int l_score =  std:: count_if(index_point_to_process.begin(),
                          index_point_to_process.end(), 
                        bind2nd(std::not_equal_to<int>(),-1));
        return l_score;
    }


    /*
    int score(float gam,  std::vector<Point_and_int>* _points, normal_iterator _it_Normal)
    {
        std::vector<int> l_points( _points->size());

        printf("before %d\n",l_points.size());
        //points inside the gamma-band around the primitive and normal of those point within alpha 
        //1 map (put -1 on every wrong point, index otherwise
        is_gamma_band_and_normal_ok pred(m_plane,gam, m_normalThresh,_it_Normal); 
        std::transform (_points->begin(),  _points->end(), l_points.begin(), pred);

        //for now we skip the cc part
        m_score = l_points.size() - (int) count (l_points.begin(), l_points.end(), -1);
        return m_score;


        //2. remove all the -1	  (not needed anymore, we test -1 and skip during 3.)
        //std::vector<int>::iterator new_end = std::remove (l_points.begin(), l_points.end(), -1) ;
        //l_points.erase(new_end, l_points.end());

        //TODO implement me
        //3. return largest cc
        //3.1   project all the point in parameter space (in the plane)
        std::vector<Point_2d> l_2D_points;
        int count = 0;

        for (	std::vector<int>::iterator it = l_points.begin(); it != l_points.end();it++, count++)
        {	  
          int index = *it;
          if (index != -1)
          {
          Point p =  boost::get<0>((*_points)[count]);
          l_2D_points.push_back(m_plane.to_2d(p));
          }
        }
        //3.2 create bitmap
        bbox_2d bb2d = bounding_box (l_2D_points.begin(), l_2D_points.end()) ;
        printf("after %d\n",l_2D_points.size());
        return 0;
    };
    */
};

//m for members
//l for local
//c for constants/readonlys
//p for pointer (and pp for pointer to pointer)
//v for volatile
//s for static
//i for indexes and iterators
//e for events


template <typename Kernel, class inputDataType>
class Efficient_ransac
{

  typedef typename std::vector<inputDataType>::iterator inputIterator;
//--------------------------------------------typedef
  typedef typename Kernel::FT FT;
    
  // basic geometric types
  typedef typename  Kernel::Point_3 Point;
  typedef typename  Kernel::Vector_3 Vector;
  typedef typename  boost::tuple<Point,int> Point_and_int;

  typedef typename  Kernel::Iso_cuboid_3 Iso_cuboid_3;

  //more advance container for geometric types
  typedef typename std::vector<Point>::const_iterator Point_const_iterator;
  typedef typename  boost::tuple<inputDataType,int> Pwn_and_int;

  typedef typename std::vector<Pwn_and_int>::iterator Pwn_and_int_iterator;

  //for kdTree
  typedef typename CGAL::Search_traits_3<Kernel> Traits_base;
  typedef typename  CGAL::Search_traits_adapter<Point_and_int,My_point_property_map<Kernel>,Traits_base>    Traits;
  typedef typename CGAL::Kd_tree<Traits> Tree;
  typedef typename  CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

  typedef Primitive_ab<Kernel, inputDataType> Primitive;

  //boost
  typedef boost::function<Point(inputDataType&)> getPointFunc;
  typedef boost::transform_iterator<getPointFunc, inputIterator> point_iterator;
  point_iterator m_it_Point, m_last_it_Point;

  typedef boost::function<Vector(inputDataType&)> getNormalFunc;
  typedef boost::transform_iterator<getNormalFunc, inputIterator> normal_iterator;
  normal_iterator m_it_Normal;
//--------------------------------------------typedef



 
//--------------------------------------------Functions
public:
  Efficient_ransac(){};
  ~Efficient_ransac(){/*if (m_global_tree) delete m_global_tree;*/};
  Efficient_ransac(inputIterator first, inputIterator beyond);

  void /*std::pair<std::vector<Primitive *>,std::vector<int> >*/ run();		 //the primitive and the index of the point sorted
  void addPrimitives(Primitive* p){m_list_sought_primitives.insert(p->type());delete p;};
  void setOptions(Option_RANSAC _m){m_options = _m;};

private:
  //Kernel:: Iso_cuboid_3 getCell(int l_index_point_p1);
  int getLevelOctree();
  Primitive* getBestCandidate(std::vector<Primitive* >& l_list_candidates,  const std::pair<std::vector<Tree*>* ,std::vector<int> >& _t, const int _SizeP,const point_iterator _it_Point, const normal_iterator _it_Normal, const std::vector<bool> &_points_available);
  inline float getRadiusSphere_from_level(int l_level_octree){return m_max_radiusSphere_Octree/powf(2,l_level_octree);}
  inline float StopProbability(float _sizeC,float _np,float _dC, float _l) const
  {
    //printf("stop without thr %f\n",	 std::pow(1.f - _sizeC/ (_np * _l * 4),_dC));
    return std::min(std::pow(1.f - _sizeC/ (_np * _l * 4),_dC), 1.f);		//4 is (1 << (m_reqSamples - 1))) with m_reqSamples=3 (min number of points to create a candidate)
  }

//--------------------------------------------Functions


//--------------------------------------------Variables
public:
  static const unsigned int scm_max_depth_octree = 4;

private:
  Option_RANSAC m_options;

  //Tree* m_global_tree;
  std::vector<Tree*> m_subsets_tree;
  std::pair<std::vector<Tree*>*,std::vector<int> > m_subsets_tree_int;
  std::vector<int> m_index_subsets;  //give the index of the subset of point i

  std::vector<std::vector<std::vector<voxel> > > voxelize_space;
  //std::vector<voxel> voxelize_space;
  float m_size_voxel; int m_width; int m_height;  int m_depth;
  Iso_cuboid_3 m_bbox;
  std::set<int> m_list_sought_primitives;
  std::vector<int> m_index_available_points;
  inputIterator m_it_Point_Normal, m_last_it_Point_Normal; //copy iterator of input points

  float m_max_radiusSphere_Octree;
  std::vector<float> m_level_weighting;  	//sum must be 1
//--------------------------------------------Variables
};

//int current(0);
//int increment () { return current++; }


//requirement: object inputDataType has x(), y(), and z() implemented

template <typename Kernel, class inputDataType>
Efficient_ransac<Kernel, inputDataType>::Efficient_ransac(inputIterator first,inputIterator beyond) 
{
  srand ( time(NULL) );

  m_it_Point_Normal = first;
  m_last_it_Point_Normal  =  beyond;

  //Index of available points
  //sequence 0, 1 , ..., nbPoints - 1
  int count = 0;
  for(inputIterator it = first; it != beyond;it++)
    m_index_available_points.push_back(count++);

    //std::iota(m_index_available_points.begin(),m_index_available_points.end(), 0);


  //init octree -------------------------
  getPointFunc f = boost::bind( &Pwn::first, _1);
  m_it_Point = boost::make_transform_iterator(first, f);
  m_last_it_Point  = boost::make_transform_iterator(beyond, f);
      
  getNormalFunc f2 = boost::bind( &Pwn::second, _1);
  m_it_Normal = boost::make_transform_iterator(first, f2);
  
  //create subsets ------------------------------------------------
  //how many subsets ?
  int l_nb_subsets = std::max( (int) std::floor(std::logf(m_index_available_points.size())/std::logf(2) )-9,   2);
  printd("number of subtree: %d\n",l_nb_subsets);


  //create a random vector of number in [0, l_nb_subsets-1] thus as this indicates to which subset the point belongs to.
  srand(time(NULL));
  std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, l_nb_subsets-1);
  std::generate_n(back_inserter(m_index_subsets), m_index_available_points.size(),[&]{ return  dis(gen); });

  //create the tree for each subset
  std::vector<std::vector<Point_and_int> > l_vector_subsets(l_nb_subsets);
  count = 0;
  for(point_iterator it =  m_it_Point; it !=  m_last_it_Point;it++, count++)
    l_vector_subsets[m_index_subsets[count]].push_back(boost::make_tuple(*it,count));

   std::vector<int> l_size_subsets;
   for (int i=0;i<l_nb_subsets;i++)
   {
     m_subsets_tree.push_back(new Tree(l_vector_subsets[i].begin(),l_vector_subsets[i].end()) );
     l_size_subsets.push_back(m_subsets_tree[i]->size());
   }

   //group subset with  number of available points in the subset (need for the scoring)
   m_subsets_tree_int	= std::make_pair(&m_subsets_tree, l_size_subsets);
   //---------------------------------------------------------------


  //bounding box and parameters for octree  
  m_bbox = bounding_box (m_it_Point,m_last_it_Point); 
  int l_max_sizeCell = ceil(max (max ( m_bbox.xmax() - m_bbox.xmin() , m_bbox.ymax() - m_bbox.ymin()), m_bbox.zmax() - m_bbox.zmin() ));			
  m_max_radiusSphere_Octree = l_max_sizeCell * sqrtf(3)/2 ;
  m_level_weighting = std::vector<float>(scm_max_depth_octree, 1.f/scm_max_depth_octree); 
  //-------------------------------------
  

  
  //init voxelize for attaching point to primitive at the end.
  m_size_voxel = sqrtf(2)*m_options.m_epsilon;		//round down

   m_width = int(ceil( (m_bbox.xmax () - m_bbox.xmin ())/m_size_voxel)) ; 
   m_height= int(ceil( (m_bbox.ymax () - m_bbox.ymin ())/m_size_voxel)) ;
   m_depth = int(ceil( (m_bbox.zmax () - m_bbox.zmin ())/m_size_voxel)) ;
  
  //voxelize_space =  std::vector< voxel >(m_width*m_height*m_depth);
  voxelize_space = std::vector< std::vector< std::vector< voxel > > >(m_width , 
                  std::vector< std::vector< voxel > >  (m_height ,
                  std::vector< voxel >				   (m_depth )));
     
  //voxelize_space.reserve(m_index_available_points.size()); 
  count = 0;
   //now put indices of input data in the voxel
  for(point_iterator it = m_it_Point; it != m_last_it_Point;it++)
  { 
    //voxelize_space[int((it->x() - m_bbox.xmin ())/m_size_voxel) + int((it->y() - m_bbox.ymin ())/m_size_voxel)*m_width + int((it->z() - m_bbox.zmin ())/m_size_voxel)*m_width*m_height ].push_back(count++); 
    voxelize_space[(it->x() - m_bbox.xmin ())/m_size_voxel][(it->y() - m_bbox.ymin ())/m_size_voxel][(it->z() - m_bbox.zmin ())/m_size_voxel].push_back(count++); 
  }

  printd("init Ransac done\n");
};	 



template <typename Kernel,  class inputDataType>
int Efficient_ransac<Kernel, inputDataType>::getLevelOctree()
{
  int l_level_octree = 0;
  float selected = float(rand())/RAND_MAX;	 //between [0, 1]

  float l_sum_CDF = m_level_weighting[0];
    while(selected > l_sum_CDF){l_sum_CDF += m_level_weighting[l_level_octree++]; if (l_level_octree == m_level_weighting.size() -1)break; };
  return l_level_octree;
};	

template <typename Kernel, class inputDataType>
Primitive_ab<Kernel,inputDataType>* Efficient_ransac<Kernel, inputDataType>::getBestCandidate(std::vector<Primitive* >& l_list_candidates,const std::pair<std::vector<Tree*>* ,std::vector<int> >& _t, const int _SizeP, const point_iterator _it_Point, const normal_iterator _it_Normal, const std::vector<bool> &_points_available)
{
  //printf("in getBestCandidate\n");
  if (l_list_candidates.size() == 1)	return l_list_candidates.back();

  int index_worse_candidate = 0;
  bool improved = true;
  while (index_worse_candidate < l_list_candidates.size()-1 && improved)		//quit if find best candidate or no more improvement possible
  {
    improved = false;
    std::sort(l_list_candidates.begin()+index_worse_candidate,l_list_candidates.end(), [](Primitive* a, Primitive* b)
    {
       return a->maxBound() < b->maxBound();
        //return a->ExpectedValue() < b->ExpectedValue();
    });
    
    /*
    for (int i=0;i<  l_list_candidates.size();i++)
       l_list_candidates.at(i)->info();
     
    system("pause");
      */

    //printf("best min/max = %f %f | min/max second %f %f\n",
    //l_list_candidates.at(l_list_candidates.size()-1)->minBound(),	l_list_candidates.at(l_list_candidates.size()-1)->maxBound(),
    //l_list_candidates.at(l_list_candidates.size()-2)->minBound(),	l_list_candidates.at(l_list_candidates.size()-2)->maxBound() );
      //l_list_candidates.back()->info();  l_list_candidates.at(l_list_candidates.size()-2)->info();

    //refine the best one 
    l_list_candidates.back()->improveBound(_t,_SizeP,_it_Point,_it_Normal, _points_available);

    int position_stop;
    //printf("-- sorted here %d, exp %f %f\n", l_list_candidates.size(), l_list_candidates.back()->minBound(),  l_list_candidates.at( l_list_candidates.size()-1)->ExpectedValue());


    //printf("best min/max = %f %f | min/max second %f %f\n",
    //l_list_candidates.at(l_list_candidates.size()-1)->minBound(),	l_list_candidates.at(l_list_candidates.size()-1)->maxBound(),
    //l_list_candidates.at(l_list_candidates.size()-2)->minBound(),	l_list_candidates.at(l_list_candidates.size()-2)->maxBound() );
      //l_list_candidates.back()->info();  l_list_candidates.at(l_list_candidates.size()-2)->info();
    //system("pause");
     //Take all those intersecting the best one
    for (position_stop= l_list_candidates.size()-2;position_stop>=index_worse_candidate;position_stop--)
    {
       if  (l_list_candidates.back()->minBound() > l_list_candidates.at(position_stop)->maxBound() ) break;//the intervals do not overlaps anymore
 
        if  (l_list_candidates.at(position_stop)->maxBound() <= m_options.m_minNbPoints) break;  //the following candidate doesnt have enough points !

        //printf("%d intersect (over %d)",position_stop, l_list_candidates.size());
       //if we reach this point, there is an overlaps between best one and position_stop
       //so request refining bound on position_stop
       improved |= l_list_candidates.at(position_stop)->improveBound(_t,_SizeP,_it_Point, _it_Normal, _points_available);//this include the next subset for computing bounds, -> the score is then returned by ExpectedValue()	

       //printf("%d improved %d\n",	position_stop,	improved);
       //if (!improved) break; //all the subsets have been used, but no improvement, quit !

        //test again after refined
       if  (l_list_candidates.back()->minBound() > l_list_candidates.at(position_stop)->maxBound() ) break;//the intervals do not overlaps anymore
    }


    //remove the condidates not processed so far
    //l_list_candidates.erase(l_list_candidates.begin(),l_list_candidates.begin()+position_stop+1);
    index_worse_candidate = position_stop;
  }

  if (index_worse_candidate == l_list_candidates.size()-1 && improved) printf("delete everything, should never happen (%d) !!!!!!\n", index_worse_candidate);
  //printf("best with value %f\n", l_list_candidates.at(0)->ExpectedValue());
   //l_list_candidates.back()->info();
  // system("pause");
  //printf("out getBestCandidate\n");



  return l_list_candidates.back();
};


template <typename Kernel, class inputDataType>
void Efficient_ransac<Kernel, inputDataType>::run()
{
   //no primitives added, exit
   if (m_list_sought_primitives.size() == 0) return;

   std::vector<Primitive* > l_result;  // Y <- 0	   (no result)
   std::vector<Primitive* > l_list_candidates;  // C <- 0	   (no candidate)

   int l_index_point_p1;
   int l_level_octree ;		//level 0 is full size, i.e 
   float l_radius_Sphere;

   float bestExp = 0;

   int numInvalid = 0; // number of points that have been assigned to a shape
   //somehow, we also need a vector<bool> to tell if the points is available
   std::vector<bool> l_points_available(m_index_available_points.size(),true);

   int nbNewCandidates = 0;
   int nbFailedCandidates = 0;  bool forceExit = false;

   printf("Starting...\n");
   do	  //main loop
     {
      bestExp = 0;
      //nbNewCandidates = 0;
     //TODO: do this in bach of n candidate, and thus can be parallelized in n threads
     do	  //generate candidate
     {
      //1. pick a point p1 randomly among available points

      l_index_point_p1 = m_index_available_points[numInvalid + rand() % (m_index_available_points.size()-numInvalid)];
      //printf("îck %d\n", 	l_index_point_p1);
  
      //2. pic a cell size randomly among cell that include the point p1
      //in our implementation, we use sphere centered on p1, thus no need to check that the cell includes p1
      l_level_octree = getLevelOctree();		//level 0 is full size, i.e 

      l_radius_Sphere = getRadiusSphere_from_level(l_level_octree);
      Fuzzy_sphere s_query(*(m_it_Point + l_index_point_p1), l_radius_Sphere, 0.1);
      std::vector<Point_and_int> l_points;

      for(int i=0;i< m_subsets_tree.size();i++)
        m_subsets_tree[i]->search(std::back_inserter(l_points), s_query);


      //3. Pick 2 more points in the l_points and give it to each primitives to create candidate
      int l_nb_try = 0;
      std::set<int> l_list_index_selected;l_list_index_selected.insert(l_index_point_p1 );
      while (	l_list_index_selected.size() < 3 &&  l_nb_try++ < l_points.size())
      {
        int id_p2 = boost::get<1>(l_points[rand() % l_points.size()]);	//pick a points in the sphere randomly
        if ( l_points_available[id_p2])												//is it free ?
          l_list_index_selected.insert(id_p2);									//put in a set (avoid duplicate with p1)

        //printf("id picked : %d %d\n",id_p2, int(l_points_available[id_p2]));
      }
      if (l_list_index_selected.size() < 3 &&  l_nb_try < l_points.size()) {printd("error !, impossible to pick 3 random points\n"); return;};

      //add candidate for each type of primitives
      for(std::set<int>::iterator it =  m_list_sought_primitives.begin(); it != m_list_sought_primitives.end();it++)		  //1 type
      {
        Primitive *p = Primitive::create(*it, m_options.m_epsilon,m_options.m_normalThresh);
        p->compute(	l_list_index_selected, s_query, m_it_Point_Normal);	//compute the primitive and says if the candidate is valid
        
        if (p->isValid())
        {
          //printf("in valid\n");  
          //NOTE: octree NEED to discard points already assigned !
          p->improveBound(m_subsets_tree_int, (m_index_available_points.size()-numInvalid), m_it_Point, m_it_Normal, l_points_available);//this include the next subset for computing bounds, -> the score is then returned by ExpectedValue()

           //printf("(%f %f    ----- %f)\n", p->minBound(), p->maxBound(), p->getLastScore());
       
           //update level sampling
          //paper does that way
          //(*sampleLevelScores)[node->Level()].first += cand.ExpectedValue();
          //++(*sampleLevelScores)[node->Level()].second;

          //evaluate the candidate
          if(p->maxBound() >= m_options.m_minNbPoints)
          {
            //printf("add new candidate %f %f\n", p->minBound(),  p->maxBound());
            if (bestExp < p->ExpectedValue()) bestExp = p->ExpectedValue();
            l_list_candidates.push_back(p);
            nbNewCandidates++;
          }else
          {
            nbFailedCandidates++;
              delete p;
          }
        
        }else{
           nbFailedCandidates++;
          delete p;
        }
      }

    if (nbFailedCandidates >= 10000)  forceExit = true;
      //printf("%d %d\n",   nbNewCandidates, nbFailedCandidates);

    // printf("%d %f ->>> %f %f\n",   nbNewCandidates, bestExp,
    //	    StopProbability(bestExp, (m_index_available_points.size()-numInvalid), nbNewCandidates, scm_max_depth_octree),
    //		StopProbability(m_options.m_minNbPoints, (m_index_available_points.size()-numInvalid), nbNewCandidates, scm_max_depth_octree));

     }while( !forceExit
         && StopProbability(bestExp, (m_index_available_points.size()-numInvalid), nbNewCandidates /*l_list_candidates.size()*/, scm_max_depth_octree)                 > m_options.m_probability 
         && StopProbability(m_options.m_minNbPoints, (m_index_available_points.size()-numInvalid), nbNewCandidates /*l_list_candidates.size()*/, scm_max_depth_octree) > m_options.m_probability);
  
    if (forceExit) break;
    // printf("criteria stop:%f\n",
    //StopProbability(bestExp, (m_index_available_points.size()-numInvalid), nbNewCandidates /*l_list_candidates.size()*/, scm_max_depth_octree)                 ,
    //StopProbability(m_options.m_minNbPoints, m_index_available_points.size() - numInvalid,  nbNewCandidates /*l_list_candidates.size()*/, scm_max_depth_octree) );

    // printf("size list candidates %d (best%f)\n", l_list_candidates.size(),bestExp);
    // for (int i=0;i<  l_list_candidates.size();i++)
    //	 l_list_candidates[i]->info();

    //now get the best candidate in the current set of all candidates
     //Note that the function sort the candidates: the best candidate is always the last element of the vector
    Primitive *best_Candidate = getBestCandidate(l_list_candidates, m_subsets_tree_int, (m_index_available_points.size()-numInvalid),m_it_Point, m_it_Normal,  l_points_available);

     //if the bestCandidate is good enough (proba of overloock lower than  m_options.m_probability)
    if (StopProbability(best_Candidate->ExpectedValue(), (m_index_available_points.size()-numInvalid), 
                                                      nbNewCandidates/*l_list_candidates.size()*/, scm_max_depth_octree) <= m_options.m_probability)
    {
      //First, put best Candidate in first in the vector
      //printf("find best candidate\n");
      //1 find the points attached to the primitives (return immediately if points already attached)
      best_Candidate->attachPoints(m_it_Point,m_it_Normal,voxelize_space, m_bbox, m_size_voxel, l_points_available);

      //we keep it
      if (best_Candidate->getPointsIndices()->size() >=  m_options.m_minNbPoints) 
      {
        l_list_candidates.back() = NULL;	//put null like that when we delete the vector, the object is not lost (pointer saved in bestCandidate)
    
        //1. add best candidate to final result.
        l_result.push_back(best_Candidate);

        //2. remove the points
        //2.1 update boolean
        std::vector<int> *indices_points_best_candidate = best_Candidate->getPointsIndices();
        //#pragma omp parallel for
        for (int i=0;i<indices_points_best_candidate->size();i++)
        {
          l_points_available[	indices_points_best_candidate->at(i) ] = false;
          m_subsets_tree_int.second[m_index_subsets[indices_points_best_candidate->at(i)]]--;	//the subset indicated by m_index_subsets has one point less
        }
      
        //2.2 swap index from m_index_available_points
        std::vector<int >::iterator bound = std::partition (m_index_available_points.begin() + numInvalid, m_index_available_points.end(), [&l_points_available](int const &p){return !l_points_available[p];});
        //printf("nbInvalid %d %d\n",numInvalid ,numInvalid +	std::distance(m_index_available_points.begin() + numInvalid, bound));
        numInvalid += std::distance(m_index_available_points.begin() + numInvalid, bound);

        //2.3 block also the points for the subtrees

         

        nbNewCandidates--;
        nbFailedCandidates = 0;
          bestExp = 0;
        //nbNewCandidates = l_result.size();

        //3. refit	(not implemented yet)
        //best_Candidate->LSfit();

        //4. save primitive
        best_Candidate->info(); 
        std::stringstream ss;
        ss<<best_Candidate->type_str()<<"_"<<l_result.size()<<".off";
        best_Candidate->save(ss.str().c_str(), m_it_Point_Normal);
        //system("pause");
      }

      
      //5. delete ALL the candidate and empty the list
       #pragma omp parallel for
       for (int i=0;i< l_list_candidates.size()-1;i++)   //-1 because we keep BestCandidate
       {
        if (  l_list_candidates[i] != NULL)
        {
          delete l_list_candidates[i];
          l_list_candidates[i] = NULL;
        }
       }

      l_list_candidates = std::vector<Primitive*>(0);
      


      /*
      //<<<<<<<<<<<<<<here what I understsood for the paper BUT it is slow, so for now try nomething else	>>>>>>>>>>>>>

      best_Candidate->save("plane.off", m_it_Point_Normal);
      system("pause");

      //compute for the other candidate or delete it
      printf("nb candidate %d\n",	 l_list_candidates.size());
       //#pragma omp parallel for
       for (int i=0;i< l_list_candidates.size()-1;i++)
       {
        l_list_candidates[i]->attachPoints(m_it_Point,m_it_Normal,voxelize_space, m_bbox, m_size_voxel);
       }
       //-------------------------------------------------------------------------------------------


      //2. remove the points
      //instead of removing, we can move the point in front, count the nb of point attached to a candidate (numInvalid), and disregard those points during future iteration
      //2.1 First we update l_points_available and will use it for 2.2
      std::vector<int> *indices_points_best_candidate = best_Candidate->getPointsIndices();
      #pragma omp parallel for
      for (int i=0;i<indices_points_best_candidate->size();i++)
      {
        l_points_available[	indices_points_best_candidate->at(i) ] = false;
      }
      
      //2.2 even better, we are going to swap index from m_index_available_points
      //we process only the part that has not been partitioned before
      //here,  all indices (from m_index_available_points) that are false in l_points_available are on the left partition
      std::vector<int >::iterator bound = std::partition (m_index_available_points.begin() + numInvalid, m_index_available_points.end(), [&l_points_available](int const &p){return !l_points_available[p];});
      numInvalid += std::distance(m_index_available_points.begin() + numInvalid, bound);

      //3. remove the candidates that have at least one point in common with bestCandidate
       #pragma omp parallel for
       for (int i=0;i< l_list_candidates.size()-1;i++)
       {
        bool deleteMe = l_list_candidates[i]->hasCommonPoints(indices_points_best_candidate);
        if (deleteMe)
        {
          delete l_list_candidates[i];
           l_list_candidates[i] = NULL;
        }
       }

       //and resize the vector of candidate afterward
      std::vector<Primitive* >::iterator bound2 = std::partition (l_list_candidates.begin(), l_list_candidates.end(), [](Primitive const *c){return !c==NULL;});
      l_list_candidates.resize(std::distance(l_list_candidates.begin(), bound2));


      //4. refit	(not implemented yet)
      //best_Candidate->LSfit();


      best_Candidate->save("plane.off", m_it_Point_Normal);
      system("pause");
      */

    }

  }
  while( 	 
      StopProbability(m_options.m_minNbPoints, m_index_available_points.size() - numInvalid,  nbNewCandidates /*l_list_candidates.size()*/, scm_max_depth_octree) > m_options.m_probability
      && float(m_index_available_points.size() - numInvalid) >= m_options.m_minNbPoints
    );




  //printf("criteria stop:%f %f\n",
  //	StopProbability(m_options.m_minNbPoints, m_index_available_points.size() - numInvalid,  nbNewCandidates /*l_list_candidates.size()*/, scm_max_depth_octree),
  //	  float(m_index_available_points.size() - numInvalid) );


  std::cout<< "run Ransac done" << std::endl;

  //return std::pair<std::vector<Primitive *>,std::vector<int> > (std::vector<Primitive *>(),std::vector<int>());
};

template <typename Kernel, class inputDataType>
Primitive_ab<Kernel,inputDataType>* Primitive_ab<Kernel,inputDataType>::create(int id, float _epsilon, float _normalThresh)
{
	switch (id)
	{
		case PLANE:
			return new Plane<Kernel,inputDataType>(_epsilon,_normalThresh);

		default:
			return NULL;
	}
};


template <typename Kernel, class inputDataType>
bool Primitive_ab<Kernel,inputDataType>::hasCommonPoints(const std::vector<int> *indices_points_best_candidate)
{
	bool commonPoint = false;
 
	#pragma omp parallel for
	for (int i=0; i<indices_points_best_candidate->size(); ++i )
	{
		#pragma omp flush (commonPoint)	    //abort at the first common indice met
		if (!commonPoint) 
		{			 
			 for (int j=0; j<m_Indices.size(); ++j )
			 {
				if (indices_points_best_candidate->at(i) == m_Indices[j])
				{
					commonPoint = true;
					break;
				}
			 }
		}
	}
	return commonPoint;
};

template <typename Kernel, class inputDataType>
void Primitive_ab<Kernel,inputDataType>::attachPoints(const point_iterator _it_Point, const normal_iterator _it_Normal, const std::vector< std::vector< std::vector<voxel> > >& _v, const Iso_cuboid_3& _bb, const float _size_c, const std::vector<bool> &_p_available)
{
	if (m_Indices.size() != 0) return;//point already attached
	//1 take a point on the surface and propagate
	Point l_start_p = projection();
		
	int l_w = int(ceil( (_bb.xmax () - _bb.xmin ())/_size_c)) ; 
	int l_h = int(ceil( (_bb.ymax () - _bb.ymin ())/_size_c)) ;
	int l_d = int(ceil( (_bb.zmax () - _bb.zmin ())/_size_c)) ;

	int id_voxel_start_x = int((l_start_p.x() - _bb.xmin ())/_size_c) ;
	int id_voxel_start_y = int((l_start_p.y() - _bb.ymin ())/_size_c) ;
	int id_voxel_start_z = int((l_start_p.z() - _bb.zmin ())/_size_c) ;
	Point voxel_start(id_voxel_start_x,id_voxel_start_y,id_voxel_start_z);

	std::vector <Point> l_result;
	l_result.push_back(voxel_start);
	std::vector<int> l_index_points;

	std::vector<Point> l_voxel_added;
	l_voxel_added.push_back(voxel_start);

	std::vector<std::vector< std::vector<bool> > > is_free(l_w,std::vector< std::vector<bool> >(l_h,std::vector<bool>(l_d,true)));
	is_free[(l_voxel_added[0].x())][(l_voxel_added[0].y())][(l_voxel_added[0].z())]=false;

	bool new_voxel_to_propagate = true;

   float x, x2, y, y2, z, z2;
   int id_current_voxel;
	do
	{
		std::vector < Point > voxel_to_add;
		new_voxel_to_propagate=false;
		for (int ii=0;ii<l_voxel_added.size();ii++)
		{
			{
				//look at neighborg
				int i,j,k;
				//todo, run me in parallel !

				//std::vector<bool> add_me(3*3*3,false);
				//#pragma omp parallel for
				for (int jj=0;jj<3*3*3 ;jj++)
				{
					i = int(jj / 9);		//(9x0,9x1,9x2)
					j = (jj - i*9) / 3;		//(3x0,3x1,3x2, 3x0,...)
					k = jj % 3;				//(0,1,2,0,...)
					i = i - 1; j= j - 1; k= k - 1;

					if (i ==0 && j == 0 && k == 0) continue;   //the current one is skipped

					x2 = (l_voxel_added[ii].x()+k)*_size_c  + _bb.xmin (); 
					y2 = (l_voxel_added[ii].y()+j)*_size_c  + _bb.ymin (); 
					z2 = (l_voxel_added[ii].z()+i)*_size_c  + _bb.zmin (); 

					if (x2 >= _bb.xmin () && x2 <= _bb.xmax () && y2 >= _bb.ymin () && y2 <= _bb.ymax () && z2 >= _bb.zmin () && z2 <= _bb.zmax ())
					{
						Point l_center_p(x2,y2,z2);

						//the voxel somehow intersects the primitive
						if (is_free[(l_voxel_added[ii].x()+k)][(l_voxel_added[ii].y()+j)][(l_voxel_added[ii].z()+i)] )							//if (_v[id_current_voxel].size() > 0)
								if (CGAL::squared_distance ( l_center_p, projection(l_center_p))*2 <= _size_c*_size_c) 	//d < sqrt(2)*c/2
								{
									//add_me[jj]=true;
									l_center_p = Point((l_voxel_added[ii].x()+k),(l_voxel_added[ii].y()+j),(l_voxel_added[ii].z()+i));
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
				//		i = int(jj / 9);		//(9x0,9x1,9x2)
				//		j = (jj - i*9) / 3;		//(3x0,3x1,3x2, 3x0,...)
				//		k = jj % 3;				//(0,1,2,0,...)
				//		i = i - 1; j= j - 1; k= k - 1;

				//		Point l_center_p = Point((l_voxel_added[ii].x()+k),(l_voxel_added[ii].y()+j),(l_voxel_added[ii].z()+i));
				//		voxel_to_add.push_back(l_center_p);
				//		l_result.push_back(l_center_p);
				//		new_voxel_to_propagate=true;
				//		is_free[(l_voxel_added[ii].x()+k)][(l_voxel_added[ii].y()+j)][(l_voxel_added[ii].z()+i)] = false;
				//	}
			}
		}

		l_voxel_added = voxel_to_add;

	}while(new_voxel_to_propagate)	;


	 //printf("end propagate %d\n",l_result.size());
	
	 //now project the points from those voxels on the primitive and save the index of the matching points !
	 int result = score(_it_Point, _it_Normal, l_index_points, _p_available, 3);
	 m_Indices = std::vector<int>(result);//score is the number of point != -1   
	 std::vector<int>::iterator new_end = std::remove_if(l_index_points.begin(), l_index_points.end(), std::bind2nd(std::equal_to<int>(), -1));
	 std::copy(l_index_points.begin(),new_end,m_Indices.begin());

};



//basically add more subset to compute score
//first parameter:
//first the subset
//second the number of used point in this subset
template <typename Kernel, class inputDataType>
bool Primitive_ab<Kernel, inputDataType>::improveBound(const std::pair<std::vector<Tree*>* ,std::vector<int> > _t, const Fuzzy_sphere _query, const int _SizeP,const point_iterator _it_Points, const normal_iterator _it_Normal, const std::vector<bool> &is_available)
{
	if (_t.first->size() == m_nb_subset_used) {printd("no more subset available\n"); return false;};


	//what it does is add another subset and recompute lower and upper bound
	//the next subset to include is provided by m_nb_subset_used

	//1. sum nb points of previous subset,	  and score
	float l_sum_score = 0;
	int l_nb_total_points_subsets = 0;
	for (int i=0;i<m_nb_subset_used;i++)
	{
		l_nb_total_points_subsets+= _t.second[i];	   //no some points are forbidden now
		l_sum_score += m_score_by_subset[i] ;
	}

	//2. need score of new subset as well as sum of the score of the previous concidered subset
	float l_new_score = 0;
				
	do
	{
		l_new_score = score(_t.first, _query,_it_Points,_it_Normal, is_available);
		//printf("subset %d exp value %f\n", m_nb_subset_used, l_new_score);
		m_score_by_subset.push_back(l_new_score);	//add to vector of score

		l_nb_total_points_subsets+= _t.second[m_nb_subset_used]; //add point new subset
		l_sum_score += l_new_score ;

		m_nb_subset_used++;
	}
	while (l_new_score == 0 && m_nb_subset_used < _t.first->size());	   

	if (l_new_score == 0) {printd("cannot improve the score\n");  return false;};

	//printf("->>>>>  %d %d\n" , l_nb_total_points_subsets, _SizeP);
	computeBound(l_nb_total_points_subsets, _SizeP, l_sum_score);//estimate the bound

	//m_nb_subset_used++;
	return true;
};

}}

#endif