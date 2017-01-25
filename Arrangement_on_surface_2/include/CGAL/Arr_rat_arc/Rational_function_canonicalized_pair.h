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

#ifndef CGAL_RATIONAL_FUNCTION_CANONICALIZED_PAIR_H
#define CGAL_RATIONAL_FUNCTION_CANONICALIZED_PAIR_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <CGAL/Arr_rat_arc/Base_rational_arc_ds_1.h>
#include <CGAL/Arr_rat_arc/Rational_function.h>
#include <CGAL/Handle_with_policy.h>
#include <CGAL/assertions.h>

namespace CGAL {
namespace Arr_rational_arc {

template <typename AlgebraicKernel_d_1>
class Rational_function_canonicalized_pair_rep:
    public Base_rational_arc_ds_1<AlgebraicKernel_d_1>
{
public:
  typedef AlgebraicKernel_d_1                           Algebraic_kernel_d_1;
  typedef Base_rational_arc_ds_1<Algebraic_kernel_d_1>  Base;

  typedef CGAL::Arr_rational_arc::Rational_function<Algebraic_kernel_d_1>
                                                        Rational_function;
  typedef typename Base::Polynomial_1                   Polynomial_1;
  typedef typename Base::Algebraic_real_1               Algebraic_real_1;
  typedef typename Base::Algebraic_vector               Algebraic_vector;
  typedef typename Base::Multiplicity                   Multiplicity;
  typedef typename Base::Multiplicity_vector            Multiplicity_vector;
  typedef typename Base::Root_multiplicity_vector       Root_multiplicity_vector;
  typedef typename Base::Solve_1                        Solve_1;
  typedef typename Base::Bound                          Bound;
  typedef typename Base::Coefficient                    Coefficient;
 
public:
  Rational_function_canonicalized_pair_rep(const Rational_function& f, 
                                           const Rational_function& g,
                                           Algebraic_kernel_d_1* ak_ptr)
    :_f(f),_g(g),_ak_ptr(ak_ptr)
  {
    CGAL_precondition(_ak_ptr != NULL);
    //canonicalized representation
    if ( !(f.id() < g.id()) ) 
      std::swap(_f,_g);
    _resultant  = (_f.numer() * _g.denom() - _f.denom() * _g.numer());

    //f and g are not the same...
    CGAL_precondition(CGAL::is_zero(_resultant) == false);
                
    Solve_1 solve_1(_ak_ptr->solve_1_object());   
    Root_multiplicity_vector rm_vec;
    
#if 1
    solve_1(_resultant,std::back_inserter(rm_vec));
#else
    CGAL_assertion(false);
    // don't use this code yet since g and resultant/g may still have a
    // common root
    if( CGAL::internal::may_have_common_factor(_f.denom(),_g.denom())){
      Polynomial_1 g = CGAL::gcd(_f.denom(), _g.denom());
      solve_1(CGAL::integral_division(_resultant, g),
              std::back_inserter(rm_vec));
      solve_1(g, std::back_inserter(rm_vec));
      std::sort(rm_vec.begin(), rm_vec.end());
    }else{
      solve_1(_resultant, std::back_inserter(rm_vec));
    }
#endif

    _roots.reserve(rm_vec.size());
    _multiplicities.reserve(rm_vec.size());
  
    for (typename Root_multiplicity_vector::iterator it = rm_vec.begin();
         it != rm_vec.end() ; 
         ++it)
    {
      _roots.push_back(it->first);
      _multiplicities.push_back(it->second);
    }

    Algebraic_vector tmp;
    tmp.reserve(_f.poles().size() + _g.poles().size());
    _event_roots.reserve(rm_vec.size() + _f.poles().size() + _g.poles().size());
    //merge the roots of f,g & resultant
    std::merge( _f.poles().begin(), _f.poles().end(),
                _g.poles().begin(), _g.poles().end(),
                std::back_inserter(tmp));
    std::merge(tmp.begin(), tmp.end(),
               _roots.begin(), _roots.end(),
               std::back_inserter(_event_roots));
  
    //remove duplicate entries
    typename Algebraic_vector::iterator new_end,old_end;
    new_end = std::unique(_event_roots.begin(), old_end=_event_roots.end());
    _event_roots.erase(new_end,old_end);

    //compute the values of _is_above
    bool curr_is_above = get_is_above_near_minus_infinity();
    _is_above.push_back(curr_is_above);

    // TBD: is_sorted is provided both by stl and boost. The interface of
    // boost version changed, and as of version ? it accepts a single
    // argument, that is a model of the SinglePassRange concept. The stl
    // version, on the other hand, is on it's way out.
    // (The name is wrong-should be is_ordered).
    //
    // CGAL_precondition(std::is_sorted(_f.poles().begin(),_f.poles().end()));
    // CGAL_precondition(std::is_sorted(_g.poles().begin(),_g.poles().end()));
    // CGAL_precondition(std::is_sorted(_roots.begin(),_roots.end()));

    typename Algebraic_vector ::const_iterator it_f_pole = _f.poles().begin();
    typename Multiplicity_vector::const_iterator it_f_mult =
      _f.pole_multiplicities().begin();
    typename Algebraic_vector ::const_iterator it_g_pole = _g.poles().begin();
    typename Multiplicity_vector::const_iterator it_g_mult =
      _g.pole_multiplicities().begin();
    typename Algebraic_vector ::const_iterator it_r_root = _roots.begin();
    typename Multiplicity_vector::const_iterator it_r_mult =
      _multiplicities.begin() ; 
    while ( (it_f_pole != _f.poles().end()) || 
        (it_g_pole != _g.poles().end()) || 
        (it_r_root != _roots.end()))
    {
      //std::cout << "_f._poles.size() " << _f._poles.size() <<std::endl;
      //std::cout << "_g._poles.size() " << _g._poles.size()<< std::endl;
      //std::cout << "_roots.size()    " << _roots.size() << std::endl;
      //if current event is only a pole of f
      if ((it_f_pole != _f.poles().end())         && 
          ((it_g_pole == _g.poles().end()) || (*it_f_pole < *it_g_pole)) &&
          ((it_r_root == _roots.end()) || (*it_f_pole < *it_r_root)))
      {
        if (*it_f_mult % 2 == 1)
          curr_is_above = (curr_is_above == true ) ? false : true;
        _is_above.push_back(curr_is_above);
        ++it_f_pole;
        ++it_f_mult;
      }
      //if current event is only a pole of g
      else if ((it_g_pole != _g.poles().end())         && 
               ((it_f_pole == _f.poles().end()) || (*it_g_pole < *it_f_pole)) &&
               ((it_r_root == _roots.end()) || (*it_g_pole < *it_r_root)))
      {
        if (*it_g_mult % 2 == 1)
          curr_is_above = (curr_is_above == true ) ? false : true;
        _is_above.push_back(curr_is_above);
        ++it_g_pole;
        ++it_g_mult;
      }
      //if current event is only a pole of r
      else if ((it_r_root != _roots.end())         && 
               ((it_f_pole == _f.poles().end()) || (*it_r_root < *it_f_pole)) &&
               ((it_g_pole == _g.poles().end()) || (*it_r_root < *it_g_pole)))
      {
        if (*it_r_mult % 2 == 1)
          curr_is_above = (curr_is_above == true ) ? false : true;
        _is_above.push_back(curr_is_above);
        ++it_r_root;
        ++it_r_mult;
      }
      //if current event is a pole of f g and r
      else if ( (it_f_pole != _f.poles().end()) && 
                (it_g_pole != _g.poles().end()) && 
                (it_r_root != _roots.end()) &&
                (*it_r_root == *it_f_pole)    && 
                (*it_r_root == *it_g_pole))
      {
        //both functions switch signs
        if (_f.sign_at(*it_r_root, CGAL::NEGATIVE) ==
            _g.sign_at(*it_r_root, CGAL::NEGATIVE))
          curr_is_above = !curr_is_above;
        if (*it_r_mult % 2 == 1)
          curr_is_above = !curr_is_above;
        if (_f.sign_at(*it_r_root, CGAL::POSITIVE) ==
            _g.sign_at(*it_r_root, CGAL::POSITIVE))
          curr_is_above = !curr_is_above;
        _is_above.push_back(curr_is_above);
        ++it_f_pole;
        ++it_f_mult;
        ++it_g_pole;
        ++it_g_mult;
        ++it_r_root;
        ++it_r_mult;
      }
      else
      {
        //should not be reached
        CGAL_postcondition_msg(false,"invalid case in computing _is_above");
      }
    }
    //std::cout << "_is_above.size()            " << _is_above.size() <<std::endl;
    //std::cout << " _event_roots.size() + 1    " <<  _event_roots.size() + 1 << std::endl;
    CGAL_postcondition(_is_above.size() == _event_roots.size() + 1);
    //check for validity using explicit computation
    CGAL_postcondition_code(std::vector<bool> tmp_is_above =
                            compute_is_above_explicitly();
                            );
    CGAL_postcondition(_is_above == tmp_is_above);   
  }

  Comparison_result compare_f_g_at(const Algebraic_real_1& x,
                                   CGAL::Sign epsilon = CGAL::ZERO) const
  {
    //f and g must be different 
    CGAL_precondition(!CGAL::is_zero(_resultant));

    //find interval 
    typename Algebraic_vector::const_iterator iter =
      std::lower_bound(_event_roots.begin(), _event_roots.end(),x);
  
    //case of a value larger than largest root
    if (iter == _event_roots.end())
      return (_is_above.back() ? CGAL:: LARGER : CGAL::SMALLER);

    typename Algebraic_vector::iterator::difference_type dist =
      iter - _event_roots.begin();

    //if x is not a root, ignore epsilons 
    if (*iter != x){
      return (_is_above[dist] ? CGAL:: LARGER : CGAL::SMALLER);
    }

    //x is a root 
    if (epsilon == CGAL::ZERO){
      CGAL_precondition(_f.poles().end() ==
                        std::find(_f.poles().begin(), _f.poles().end(),x));
      CGAL_precondition(_g.poles().end() ==
                        std::find(_g.poles().begin(), _g.poles().end(),x));
      return (CGAL::EQUAL);
    }
    
    if (epsilon == CGAL::NEGATIVE)
      return (_is_above[dist] ? CGAL:: LARGER : CGAL::SMALLER);
    else // CGAL::POSITIVE
      return (_is_above[dist+1] ? CGAL:: LARGER : CGAL::SMALLER);
  }
  Comparison_result compare_f_g_at(Arr_parameter_space boundary) const
  {
    CGAL_precondition((boundary == ARR_LEFT_BOUNDARY) ||
        (boundary == ARR_RIGHT_BOUNDARY) );
  
    //f and g are the same...
    if (CGAL::is_zero(_resultant))
      return CGAL::EQUAL;

    if (boundary == ARR_LEFT_BOUNDARY)
      return _is_above.front() ? CGAL::LARGER : CGAL::SMALLER ;
    else //  boundary = ARR_RIGHT_BOUNDARY
      return _is_above.back()  ? CGAL::LARGER : CGAL::SMALLER ;
  }

  bool is_intersecting_in_range(const Arr_parameter_space left_parameter_space,
                                const Algebraic_real_1 left,
                                const Arr_parameter_space right_parameter_space,
                                const Algebraic_real_1 right) const
  {
    //the two function intersect iff the left index and the right index
    //of the array _is_above are different

    //get left index
    typename Algebraic_vector::const_iterator left_index =
      (left_parameter_space == ARR_LEFT_BOUNDARY) ? 
      _event_roots.begin():
      std::lower_bound(_event_roots.begin(), _event_roots.end(),left);
    //check if intersect at left index
    if (*left_index == left) 
      return true;

    //get right index
    typename Algebraic_vector::const_iterator right_index =
      (right_parameter_space == ARR_RIGHT_BOUNDARY) ? 
      _event_roots.begin():
      std::lower_bound(_event_roots.begin(), _event_roots.end(),right);
    //check if intersect at right index
    if (*right_index == right) 
      return true;
  
    //check if indices are the same
    return (left_index == right_index);
  }

  const Rational_function& f() const 
  {
    return _f;
  }

  const Rational_function& g() const 
  {
    return _g;
  }

  const Algebraic_vector & roots() const
  {
    return _roots;
  }

  const Multiplicity_vector & multiplicities() const
  {
    return _multiplicities;
  }

private:
  bool get_is_above_near_minus_infinity()
  {
    bool r = _get_is_above_near_minus_infinity();
    CGAL_postcondition(r == __get_is_above_near_minus_infinity());
    return r; 
  }

  bool _get_is_above_near_minus_infinity()
  {
    int f_deg = CGAL::degree(_f.numer()) - CGAL::degree(_f.denom());
    int g_deg = CGAL::degree(_g.numer()) - CGAL::degree(_g.denom());

    if (f_deg > g_deg)  //f is stronger than g
      {
        CGAL::Sign sign = _f.sign_near_minus_infinity();
        return  (sign == CGAL::NEGATIVE) ? false  :
          (sign == CGAL::POSITIVE) ? true   :
          (_g.sign_near_minus_infinity() == CGAL::NEGATIVE);  // _f == zero;
      }
    if (f_deg < g_deg)  //g is stronger than f
      {
        CGAL::Sign sign = _g.sign_near_minus_infinity();
        return  (sign == CGAL::NEGATIVE) ? true  :
          (sign == CGAL::POSITIVE) ? false   :
          (_f.sign_near_minus_infinity() == CGAL::POSITIVE);  // _g == zero;
      }
            
    //both have the same degree difference, 
    //check who's leading coeeficient ratio is larger
    Coefficient lead_coeff_ratio =
      CGAL::abs(CGAL::leading_coefficient(_f.numer())) *
      CGAL::leading_coefficient(_g.denom()) -
      CGAL::abs(CGAL::leading_coefficient(_g.numer())) *
      CGAL::leading_coefficient(_f.denom());
    if (lead_coeff_ratio > 0)//f is stronger than g
      return (_f.sign_near_minus_infinity() == CGAL::NEGATIVE ) ? false : true;
    if (lead_coeff_ratio < 0)  //g is stronger than f
      return (_g.sign_near_minus_infinity() == CGAL::NEGATIVE ) ? true : false;

    //ratio is the same, instead of continuing with taylor expansion,
    //compute explicitly

    return __get_is_above_near_minus_infinity();
  }

  bool __get_is_above_near_minus_infinity()
  {
    Bound b;
    if (_event_roots.empty())
      b = Bound(0);
    else
      {
        b = (_ak_ptr->approximate_relative_1_object()(_event_roots.front(), 0)).first - 1;  //lower bound of first root
      }
    return is_above_at(b);
  }

  bool is_above_at(const Bound& b)
  {
    //return true if f is above g at b which means return (sign  == positive) of :
    //
    //  numer_f * denom _g - numer_f * denom _g
    //f-g= -------------------------------------- at b
    //    denom _f  * denom _g 

    //TODO: unnescecary construction of real
    Algebraic_real_1 x(_ak_ptr->construct_algebraic_real_1_object()(b));
  
    CGAL::Sign   numer = _ak_ptr->sign_at_1_object()(_resultant, x);
    CGAL::Sign   denom_f = _ak_ptr->sign_at_1_object()(_f.denom(),x);
    CGAL::Sign   denom_g = _ak_ptr->sign_at_1_object()(_g.denom(),x);
    CGAL::Sign   denom = denom_f * denom_g ;
    CGAL::Sign   s = numer*denom;

    CGAL_precondition(s != CGAL::ZERO);

    return (s == CGAL::POSITIVE);
  }
  std::vector<bool> compute_is_above_explicitly()
  {
    std::vector<bool> tmp_is_above;
    if (_event_roots.size()== 0)
    {
      Bound b = 1; //all bound are legal, choose 1 for simplicity
      tmp_is_above.push_back(is_above_at(b));  
      return tmp_is_above;
    }
  
    tmp_is_above.reserve(_event_roots.size()+1);  
    //left boundary
    Bound b  = (_ak_ptr->approximate_relative_1_object()
      (_event_roots.front(),0)).first - 1;  //lower bound of first root
    tmp_is_above.push_back(is_above_at(b));  
  
    //mid intervals
    typename Algebraic_vector::size_type i;
    for (i = 0; i < _event_roots.size()-1; ++i)
    {
      b = _ak_ptr->bound_between_1_object()(_event_roots[i],_event_roots[i+1]);
      tmp_is_above.push_back(is_above_at(b));      
    }

    //right boundary
    b = (_ak_ptr->approximate_relative_1_object()(_event_roots.back(), 0)).second + 1;  //lower bound of last root
    tmp_is_above.push_back(is_above_at(b));
    return tmp_is_above;
  }
   
private:   
  Rational_function _f,_g;
  Polynomial_1  _resultant;
  Algebraic_vector _roots;
  //roots of resultant merged with roots of f & g's denomenators
  Algebraic_vector _event_roots;
  Multiplicity_vector _multiplicities;
  //is f above g in the interval induced by index i of _roots
  std::vector<bool> _is_above;
  Algebraic_kernel_d_1* _ak_ptr;
}; // Rational_function_canonicalized_pair_rep

template <typename Algebraic_kernel_>
class Rational_function_canonicalized_pair:
    public Handle_with_policy<Rational_function_canonicalized_pair_rep
                              <Algebraic_kernel_> >
{
public:
  typedef Algebraic_kernel_                           Algebraic_kernel_d_1;
  typedef Rational_function_canonicalized_pair_rep<Algebraic_kernel_d_1>
                                                      Rep;
  typedef Handle_with_policy<Rep>                     Base;
  typedef Rational_function_canonicalized_pair<Algebraic_kernel_d_1>
                                                      Self;
  typedef typename Rep::Rational_function             Rational_function;
  typedef typename Rep::Algebraic_real_1              Algebraic_real_1;
  typedef typename Rep::Polynomial_1                  Polynomial_1;
  typedef typename Rep::Algebraic_vector              Algebraic_vector;
  typedef typename Rep::Multiplicity_vector           Multiplicity_vector;
  typedef typename Rep::Root_multiplicity_vector      Root_multiplicity_vector;

  Rational_function_canonicalized_pair(const Rational_function& f, 
                                       const Rational_function& g,
                                       Algebraic_kernel_d_1* ak_ptr) :
    Base(f, g, ak_ptr) {}
  //used to solve VS bug...
  Rational_function_canonicalized_pair(const Self & p) :
    Base(static_cast<const Base &> (p)) {}

  Comparison_result compare_f_g_at(const Algebraic_real_1& x,
                                   CGAL::Sign epsilon = CGAL::ZERO) const
  {
    return this->ptr()->compare_f_g_at(x,epsilon);
  }

  Comparison_result compare_f_g_at(Arr_parameter_space boundary) const
  {
    return this->ptr()->compare_f_g_at(boundary);
  }

  bool is_intersecting_in_range(const Arr_parameter_space left_parameter_space,
                                const Algebraic_real_1 left,
                                const Arr_parameter_space right_parameter_space,
                                const Algebraic_real_1 right) const
  {
    return this->ptr()->is_intersecting_in_range(left_parameter_space, left,
                                                 right_parameter_space, right);
  }
    
  const Rational_function& f() const 
  {
    return this->ptr()->f();
  }
    
  const Rational_function& g() const 
  {
    return this->ptr()->g();
  }

  const Algebraic_vector & roots() const
  {
    return this->ptr()->roots();
  }

  const Multiplicity_vector & multiplicities() const
  {
    return this->ptr()->multiplicities();
  }
};  //Rational_function_canonicalized_pair

}   //namespace Arr_rational_arc
}   //namespace CGAL {   

#endif //CGAL_RATIONAL_FUNCTION_CANONICALIZED_PAIR_H
