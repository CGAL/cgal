// Copyright (c) 2011 Tel-Aviv University (Israel), INRIA Sophia-Antipolis (France).
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
// Author(s)     : Oren Salzman <orenzalz@post.tau.ac.il >
//                 Michael Hemmer <Michael.Hemmer@sophia.inria.fr>

#ifndef CGAL_RATIONAL_ARC_CACHE
#define CGAL_RATIONAL_ARC_CACHE

#include <CGAL/Arr_rat_arc/Base_rational_arc_ds_1.h>
#include <CGAL/Arr_rat_arc/Rational_function.h>
#include <CGAL/Arr_rat_arc/Rational_function_pair.h>
#include <CGAL/Arr_rat_arc/Rational_function_canonicalized_pair.h>

namespace CGAL {
namespace Arr_rational_arc {

//-------------------
//Cache 
//-------------------
template <typename AlgebraicKernel_d_1>
class Cache : public Base_rational_arc_ds_1<AlgebraicKernel_d_1>
{
public:
  typedef AlgebraicKernel_d_1                           Algebraic_kernel_d_1;
  typedef Base_rational_arc_ds_1<Algebraic_kernel_d_1>  Base;
  typedef Cache<Algebraic_kernel_d_1>                   Self;
 
  typedef typename Base::Polynomial_1                   Polynomial_1;
  typedef typename Base::Rational                       Rational;
  typedef typename Base::Algebraic_real_1               Algebraic_real_1;
  typedef typename Base::Integer                        Integer ;
  typedef typename Base::FT_rat_1                       FT_rat_1;
  typedef typename Base::Polynomial_traits_1            Polynomial_traits_1;

  typedef CGAL::Arr_rational_arc::Rational_function<Algebraic_kernel_d_1>
    Rational_function;
  typedef CGAL::Arr_rational_arc::Rational_function_canonicalized_pair<Algebraic_kernel_d_1>
    Rational_function_canonicalized_pair;
  typedef CGAL::Arr_rational_arc::Rational_function_pair<Algebraic_kernel_d_1>
    Rational_function_pair;

  typedef std::pair<Polynomial_1,Polynomial_1>          Rational_function_key;
  class Less_compare_rational_function_key
  {
  public:
    Comparison_result compare(const Rational_function_key& k1,
                              const Rational_function_key& k2) const
    {
      Comparison_result cr =
        typename Polynomial_traits_1::Compare()(k1.first, k2.first);
      if (cr != CGAL::EQUAL)
        return (cr);
      cr = typename Polynomial_traits_1::Compare()(k1.second, k2.second);

      return cr;
    }

    bool operator()(const Rational_function_key& k1,
                    const Rational_function_key& k2) const
    {
      Comparison_result cr = compare(k1,k2);
      return ((cr == CGAL::LARGER) ? true : false);
    }
  };

  typedef typename std::map<Rational_function_key,
                            Rational_function,
                            Less_compare_rational_function_key>
                                                      Rational_function_map;

  typedef typename Rational_function::Id_type         Rational_function_id_type;
  typedef std::pair<Rational_function_id_type, Rational_function_id_type>
    Rational_function_canonicalized_pair_key;

  class Less_compare_rational_function_pair_key
  {
  public:
    bool operator()(const Rational_function_canonicalized_pair_key& k1,
        const Rational_function_canonicalized_pair_key& k2) const
    {
      return (k1<k2);
    }
  };

  typedef typename std::map< Rational_function_canonicalized_pair_key,
                             Rational_function_canonicalized_pair,
                             Less_compare_rational_function_pair_key>
    Rational_function_canonicalized_pair_map;

public:
  Cache() :
    _rat_func_map_watermark(128), _rat_pair_map_watermark(128), _ak_ptr(NULL){};

  void initialize(Algebraic_kernel_d_1* ak_ptr)
  {
    _ak_ptr = ak_ptr;
  }
  void initialize(const Self& other, Algebraic_kernel_d_1* ak_ptr)
  {
    //copy kernel pointer
    _ak_ptr = ak_ptr;

    //copy rational function map
    typename Rational_function_map::const_iterator iter1;
    for ( iter1 = other.rat_func_map().begin();
          iter1 != other.rat_func_map().end();
          ++iter1)
    {
      if (iter1->second.is_shared())
      { 
        Rational_function_key key   = iter1->first;
        //construct new instance
        Rational_function f(iter1->second.numer(), iter1->second.denom(),
                            _ak_ptr);
        _rat_func_map.insert(std::make_pair(key,f)); 
      }
    }
    
    //copy rational function pair map
    typename Rational_function_canonicalized_pair_map::const_iterator iter2;
    for ( iter2  = other.rat_pair_map().begin();
          iter2 != other.rat_pair_map().end();
          ++iter2)
    {
      if (iter2->second.is_shared())
      {
        Rational_function_canonicalized_pair_key key  = iter2->first;
        //construct new instance
        Rational_function_canonicalized_pair p(iter2->second.f(),
                                               iter2->second.g(), _ak_ptr);
        _rat_pair_map.insert(std::make_pair(key,p)); 
      }
    }

  }
  const Rational_function_map& rat_func_map() const
  {
    return _rat_func_map;
  }

  const Rational_function_canonicalized_pair_map& rat_pair_map() const
  {
    return _rat_pair_map;
  }
  
  const Rational_function& get_rational_function(const Polynomial_1& numer,
                                                 const Polynomial_1& denom) const
  {
    CGAL_precondition (_ak_ptr != NULL);
    Rational_function_key key  = get_key(numer,denom);

    //look if element exists in cache already
    typename Rational_function_map::iterator it = _rat_func_map.lower_bound(key);

    if(it != _rat_func_map.end() && !(_rat_func_map.key_comp()(key, it->first)))
      {
        return it->second;  //iterator to a type <key,value>
      }
    else    //element does not exist, create it & insert to cache
      {
        //first check if to clean up cache
         if (_rat_func_map.size() > _rat_func_map_watermark)
          rat_func_map_clean_up();

        //then insert the new element
        Rational_function f(numer,denom,_ak_ptr);
        typename Rational_function_map::iterator it2 =
          _rat_func_map.insert(it,std::make_pair(key,f)); 
        return it2->second; 
      } 
  }
  const Rational_function&  get_rational_function( const Rational& rat) const
  {
    Integer  numer,denom;
    typename FT_rat_1::Decompose()(rat,numer,denom);

    Polynomial_1 numer_poly =
      typename Polynomial_traits_1::Construct_polynomial()(numer);
    Polynomial_1 denom_poly =
      typename Polynomial_traits_1::Construct_polynomial()(denom);

    return get_rational_function (numer_poly,denom_poly);
  }

  const Rational_function_pair get_rational_pair(const Rational_function& f, 
                                                 const Rational_function& g) const
  {
    CGAL_precondition (_ak_ptr != NULL);
    CGAL_precondition(!(f==g));
    Rational_function_canonicalized_pair_key key  = get_key(f,g);
    bool is_opposite = (f.id() < g.id()) ? false : true ; 

    //look if element exists in cache already
    typename Rational_function_canonicalized_pair_map::iterator it =
      _rat_pair_map.lower_bound(key);
  
    if(it != _rat_pair_map.end() && !(_rat_pair_map.key_comp()(key, it->first)))
    {
      return (Rational_function_pair(it->second,is_opposite));
    }
    else    //element does not exist, 
    {
      //first check if to clean up cache
      if (_rat_pair_map.size() > _rat_pair_map_watermark)
        rat_pair_map_clean_up();

      //create it & insert to cache
      Rational_function_canonicalized_pair p(f, g, _ak_ptr);
      std::pair<typename Rational_function_canonicalized_pair_map::const_iterator, bool> res = 
        _rat_pair_map.insert(std::make_pair(key,p));
      return (Rational_function_pair(res.first->second,is_opposite));
    }
  }

  void cleanup() const
  {
    rat_pair_map_clean_up();
    rat_func_map_clean_up();
  }
private:
  Rational_function_key get_key(const Polynomial_1& numer,
                                const Polynomial_1& denom) const 
  {
    return Rational_function_key(numer, denom);
  }

  Rational_function_key get_key(const Rational_function& f) const
  {
    return get_key(f.numer(), f.denom());
  }

  Rational_function_canonicalized_pair_key
  get_key(const Rational_function& f, const Rational_function& g) const
  {
    return (f.id() < g.id()) ?
      Rational_function_canonicalized_pair_key( f.id(),g.id() ):
      Rational_function_canonicalized_pair_key( g.id(),f.id() );

  }

  void rat_func_map_clean_up() const 
  {                                             
   
    //find eraseable rational functions
    std::vector<Rational_function_key> eraseable;
    typename Rational_function_map::iterator iter1;
    for ( iter1 = _rat_func_map.begin();
          iter1 != _rat_func_map.end();
          ++iter1)
    {
      if (iter1->second.is_shared() == false)
        eraseable.push_back(iter1->first);
    }

    //erase functions
    typename std::vector<Rational_function_key>::iterator iter2;
    for ( iter2 = eraseable.begin();
          iter2 != eraseable.end();
          ++iter2)
    {
      _rat_func_map.erase(*iter2);
    }

    //re-set watermark
    _rat_func_map_watermark = (std::max)(
        2*_rat_func_map.size(),
        typename Rational_function_map::size_type(128));
    
    return;
  }
  void rat_pair_map_clean_up() const 
  {
    //find eraseable rational functions
    std::vector<Rational_function_canonicalized_pair_key> eraseable;
    typename Rational_function_canonicalized_pair_map::iterator iter1;
    for ( iter1 = _rat_pair_map.begin();
          iter1 != _rat_pair_map.end();
          ++iter1)
    {
      if (iter1->second.is_shared() == false)
        eraseable.push_back(iter1->first);
    }

    //erase functions
    typename std::vector<Rational_function_canonicalized_pair_key>::iterator iter2;
    for ( iter2 = eraseable.begin();
          iter2 != eraseable.end();
          ++iter2)
    {
      _rat_pair_map.erase(*iter2);
    }

    //re-set watermark
    _rat_pair_map_watermark = 
      (std::max)(2*_rat_pair_map.size(),
          typename Rational_function_canonicalized_pair_map::size_type(128));
  }

private:
  mutable Rational_function_map   _rat_func_map;
  mutable unsigned int                    _rat_func_map_watermark;
  mutable Rational_function_canonicalized_pair_map  _rat_pair_map;
  mutable unsigned int                    _rat_pair_map_watermark;
  mutable Algebraic_kernel_d_1*   _ak_ptr;
}; //Cache
 
}   //namespace Arr_rational_arc
}   //namespace CGAL {
#endif //CGAL_RATIONAL_ARC_CACHE
