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
template < class Algebraic_kernel_ >
class Cache : public Base_rational_arc_ds_1<Algebraic_kernel_>
{
public:
  typedef Algebraic_kernel_                        Algebraic_kernel;
  typedef Base_rational_arc_ds_1<Algebraic_kernel> Base;
 
  typedef typename Base::Polynomial_1         Polynomial_1;
  typedef typename Base::Rational             Rational;
  typedef typename Base::Algebraic_real_1     Algebraic_real_1;
  typedef typename Base::Integer              Integer ;
  typedef typename Base::FT_rat_1             FT_rat_1;
  typedef typename Base::Polynomial_traits_1  Polynomial_traits_1;

  typedef CGAL::Arr_rational_arc::Rational_function<Algebraic_kernel>    Rational_function;
  typedef CGAL::Arr_rational_arc::Rational_function_canonicalized_pair<Algebraic_kernel>  Rational_function_canonicalized_pair;
  typedef CGAL::Arr_rational_arc::Rational_function_pair<Algebraic_kernel> Rational_function_pair;

  typedef std::pair<Polynomial_1,Polynomial_1>      Rational_function_key;
  class Less_compare_rational_function_key
  {
  public:
    Comparison_result compare(const Rational_function_key& k1, const Rational_function_key& k2) const
    {
      Comparison_result cr = typename Polynomial_traits_1::Compare()(k1.first,k2.first);
      if (cr != CGAL::EQUAL)
        return (cr);
      cr = typename Polynomial_traits_1::Compare()(k1.second,k2.second);

      return cr;
    }
    bool operator()(const Rational_function_key& k1, const Rational_function_key& k2) const
    {
      Comparison_result cr = compare(k1,k2);
      return ((cr == CGAL::LARGER) ? true : false);
    }
  };

  typedef typename std::map< Rational_function_key,
                             Rational_function,
                             Less_compare_rational_function_key>
                                                      Rational_function_map;

  typedef typename Rational_function::Id_type         Rational_function_id_type;
  typedef std::pair<Rational_function_id_type,Rational_function_id_type>
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
        Less_compare_rational_function_pair_key> Rational_function_canonicalized_pair_map;
public:
  Cache() : _rat_func_map_watermark(128),_rat_pair_map_watermark(128){};
  
  const Rational_function&  get_rational_function(const Polynomial_1& numer,
                                                  const Polynomial_1& denom,
                                                  Algebraic_kernel kernel = Algebraic_kernel()) 
  {
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

        //then inert the new element
        Rational_function f(numer,denom,kernel);
        typename Rational_function_map::iterator it2 = _rat_func_map.insert(it,std::make_pair(key,f)); 
        return it2->second; 
      } 
  }
  const Rational_function&  get_rational_function( const Rational& rat,
      Algebraic_kernel kernel = Algebraic_kernel()) 
  {
    Integer  numer,denom;
    typename FT_rat_1::Decompose()(rat,numer,denom);

    Polynomial_1 numer_poly = typename Polynomial_traits_1::Construct_polynomial()(numer);
    Polynomial_1 denom_poly = typename Polynomial_traits_1::Construct_polynomial()(denom);

    return get_rational_function (numer_poly,denom_poly,kernel);
  }

  const Rational_function_pair get_rational_pair ( const Rational_function& f, 
      const Rational_function& g,
      Algebraic_kernel kernel = Algebraic_kernel()) 
  {
    CGAL_precondition(!(f==g));
    Rational_function_canonicalized_pair_key key  = get_key(f,g);
    bool is_opposite = (f.id() < g.id()) ? false : true ; 

    //look if element exists in cache already
    typename Rational_function_canonicalized_pair_map::iterator it = _rat_pair_map.lower_bound(key);
  
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
        Rational_function_canonicalized_pair p(f,g,kernel);
        typename Rational_function_canonicalized_pair_map::const_iterator it2 = 
          _rat_pair_map.insert(it,std::make_pair(key,p));
        return (Rational_function_pair(it2->second,is_opposite));
      }
  }

  void cleanup() 
  {
    rat_func_map_clean_up();
    rat_pair_map_clean_up();
  }
private:
  Rational_function_key get_key(const Polynomial_1& numer, const Polynomial_1& denom) const 
  {
    return Rational_function_key(numer, denom);
  }
  Rational_function_key get_key(const Rational_function& f) const
  {
    return get_key( f.numer(),f.denom());
  }
  Rational_function_canonicalized_pair_key get_key(const Rational_function& f, const Rational_function& g) const
  {
    return (f.id() < g.id()) ?
      Rational_function_canonicalized_pair_key( f.id(),g.id() ):
      Rational_function_canonicalized_pair_key( g.id(),f.id() );

  }

  void rat_func_map_clean_up()
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
    _rat_func_map_watermark = 2*_rat_func_map.size();
    return;
  }
  void rat_pair_map_clean_up()
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
    _rat_pair_map_watermark = 2*_rat_pair_map.size();
    return;
  }
private:
  mutable Rational_function_map                     _rat_func_map;
  unsigned int _rat_func_map_watermark;
  mutable Rational_function_canonicalized_pair_map  _rat_pair_map;
  unsigned int _rat_pair_map_watermark;
}; //Cache
 
}   //namespace Arr_rational_arc
}   //namespace CGAL {
#endif //CGAL_RATIONAL_ARC_CACHE
