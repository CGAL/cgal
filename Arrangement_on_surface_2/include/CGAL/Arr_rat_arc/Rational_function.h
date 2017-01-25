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

#ifndef CGAL_RATIONAL_FUNCTION_H
#define CGAL_RATIONAL_FUNCTION_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <CGAL/Arr_rat_arc/Base_rational_arc_ds_1.h>
#include <CGAL/Handle_with_policy.h>
namespace CGAL {
namespace Arr_rational_arc {

template <class AlgebraicKernel_d_1 >
class Rational_function_rep : public Base_rational_arc_ds_1<AlgebraicKernel_d_1>
{
public:
  typedef AlgebraicKernel_d_1                          Algebraic_kernel_d_1;
  typedef Base_rational_arc_ds_1<Algebraic_kernel_d_1> Base;

  typedef typename Base::Polynomial_1                  Polynomial_1;
  typedef typename Base::Algebraic_real_1              Algebraic_real_1;
  typedef typename Base::Algebraic_vector              Algebraic_vector;
  typedef typename Base::Multiplicity                  Multiplicity;
  typedef typename Base::Multiplicity_vector           Multiplicity_vector;
  typedef typename Base::Root_multiplicity_vector      Root_multiplicity_vector;
  typedef typename Base::Solve_1                       Solve_1;
    
public:
  Rational_function_rep() : _ak_ptr(NULL){}
  Rational_function_rep(const Polynomial_1& numer,
                        const Polynomial_1& denom, 
                        Algebraic_kernel_d_1* ak_ptr):
    _numer(numer), _denom(denom),_ak_ptr(ak_ptr)
  {
    initialize();
  }
 
  CGAL::Sign sign_at (const Algebraic_real_1& x,
                      CGAL::Sign epsilon = CGAL::ZERO) const
  {
    //find interval 
    typename Algebraic_vector::const_iterator iter =
      std::lower_bound(_event_roots.begin(), _event_roots.end(),x);
  
    //case of a value larger than largest root
    if (iter == _event_roots.end())
      return (_sign.back());

    typename Algebraic_vector::iterator::difference_type dist =
      iter - _event_roots.begin();

    //if x is not a root, ignore epsilons 
    if (*iter != x)
      return (_sign[dist]);

    //x is a root 
    if (epsilon == CGAL::ZERO)
      return (CGAL::EQUAL);
    else if (epsilon == CGAL::NEGATIVE)
      return (_sign[dist] );
    else // CGAL::POSITIVE
      return (_sign[dist+1]);
  }

  CGAL::Sign sign_near_minus_infinity() const
  {
    return _sign.front();
  }

  bool operator==(const Rational_function_rep& other) const
  {
    return ((this->_numer == other.numer()) &&
            (this->_denom == other.denom()) );
  }

  const Polynomial_1& numer() const
  {
    return _numer;
  }

  const Polynomial_1& denom() const
  {
    return _denom;
  }

  const Algebraic_vector& poles() const
  {
    return _poles;
  }

  const Multiplicity_vector& pole_multiplicities() const
  {
    return _pole_multiplicities;
  }
    
private:
  void initialize()
  {
    CGAL_precondition(_ak_ptr != NULL);
    CGAL_precondition(CGAL::is_zero(_denom) == false);
    if (CGAL::is_zero(_numer))
    {
      //function does not change sign
      _sign.push_back(CGAL::ZERO);
      return;
    }

    Solve_1 solve_1 (_ak_ptr->solve_1_object());
    Root_multiplicity_vector rm_poles_vec,rm_intersctions_vec;
    solve_1(_denom, std::back_inserter(rm_poles_vec));  //poles
    solve_1(_numer, std::back_inserter(rm_intersctions_vec)); //intersections with zero

    //reserve memory
    typename Root_multiplicity_vector::size_type num_of_poles =
      rm_poles_vec.size();
    typename Root_multiplicity_vector::size_type num_of_intersections =
      rm_intersctions_vec.size();
  
    _poles.reserve(num_of_poles);
    _pole_multiplicities.reserve(num_of_poles);

    _event_roots.reserve(num_of_poles + num_of_intersections);
    _pole_multiplicities.reserve(num_of_poles + num_of_intersections);

    //initialize poles
    for ( typename Root_multiplicity_vector::iterator it = rm_poles_vec.begin(); 
          it != rm_poles_vec.end() ;
          ++it)
    {
      _poles.push_back(it->first);
      _pole_multiplicities.push_back(it->second);
    }

    //initialize events
    Root_multiplicity_vector events_mult_vec;
    std::merge( rm_poles_vec.begin()  , rm_poles_vec.end()  ,
                rm_intersctions_vec.begin() , rm_intersctions_vec.end() ,
                std::back_inserter(events_mult_vec));

    typename Root_multiplicity_vector::iterator it;
    for (it = events_mult_vec.begin(); it != events_mult_vec.end(); ++it)
    {
      _event_roots.push_back(it->first);
      _event_multiplicities.push_back(it->second);
    }

    //initialize left most interval (at -oo)
    //for (ax^n+.../x^m+...) the sign is:  sign(a) * (-1)^(n +m)
    CGAL::Sign curr_sign = CGAL::sign(CGAL::leading_coefficient(_numer));
    if ((CGAL::degree(_numer) + CGAL::degree(_denom)) % 2 == 1)
      curr_sign = curr_sign * CGAL::NEGATIVE;
    _sign.push_back(curr_sign);

    typename Multiplicity_vector::iterator it2;
    for (it2 = _event_multiplicities.begin(); 
         it2 != _event_multiplicities.end(); 
         ++it2)
    {
      if (*it2 % 2 == 1)
        curr_sign = curr_sign * CGAL::NEGATIVE;
      _sign.push_back(curr_sign);
    }
  }

private:
  Polynomial_1        _numer;
  Polynomial_1        _denom;
  Algebraic_vector    _poles;     //roots of the denominator
  Multiplicity_vector _pole_multiplicities; //multiplicities of the poles
  Algebraic_vector    _event_roots;   //poles and intersection points with y=0, function can change signs in these events
  Multiplicity_vector _event_multiplicities; //multiplicities of the events
  std::vector<CGAL::Sign> _sign;    //function's sign in the corresponding interval induced by _event_roots (if no roots then only one value)
  mutable Algebraic_kernel_d_1*   _ak_ptr;

};//Rational_function_rep 

template < class Algebraic_kernel_ >
class Rational_function:
    public Handle_with_policy<Rational_function_rep<Algebraic_kernel_> >
{
public:
  typedef Algebraic_kernel_                             Algebraic_kernel_d_1;
  typedef Handle_with_policy<Rational_function_rep<Algebraic_kernel_> >
                                                        Base;
  typedef Rational_function<Algebraic_kernel_d_1>       Self;
  typedef Rational_function_rep<Algebraic_kernel_>      Rep;
  typedef typename Rep::Algebraic_real_1                Algebraic_real_1;
  typedef typename Rep::Polynomial_1                    Polynomial_1;
  typedef typename Rep::Algebraic_vector                Algebraic_vector;
  typedef typename Rep::Multiplicity_vector             Multiplicity_vector;

  typedef typename Base::Id_type                        Id_type;
private:
  static Self& get_default_instance()
  {
    static Algebraic_kernel_d_1 kernel;
    static Self x = Self(Polynomial_1(0), Polynomial_1(1), &kernel); 
    return x; 
  } 
public:
  Rational_function(const Polynomial_1& numer,
                    const Polynomial_1& denom, 
                    Algebraic_kernel_d_1* ak_ptr) :
    Base(numer,denom,ak_ptr) {}

  //used to solve VS bug...
  Rational_function () :
    Base(static_cast<const Base &> (get_default_instance())) {}

  // explicit copy-constructor, required by VC9
  Rational_function (const Self & r)
    : Base(static_cast<const Base &> (r)) {}

  CGAL::Sign sign_at(const Algebraic_real_1& x,
                     CGAL::Sign epsilon = CGAL::ZERO) const
  {
    return this->ptr()->sign_at(x,epsilon);
  }

  CGAL::Sign sign_near_minus_infinity() const
  {
    return this->ptr()->sign_near_minus_infinity ();
  }
  bool operator== (const Self & other) const
  {
    if (this->is_identical (other))
      return true;
    return (*(this->ptr()) == *(other.ptr()));
  }

  const Polynomial_1& numer () const
  {
    return this->ptr()->numer();
  }

  const Polynomial_1& denom () const
  {
    return this->ptr()->denom();
  }

  const Algebraic_vector& poles () const
  {
    return this->ptr()->poles();
  }

  const Multiplicity_vector& pole_multiplicities () const
  {
    return this->ptr()->pole_multiplicities();
  }
};  //Rational_function

}   //namespace Arr_rational_arc
}   //namespace CGAL {

#endif //CGAL_RATIONAL_FUNCTION_H
