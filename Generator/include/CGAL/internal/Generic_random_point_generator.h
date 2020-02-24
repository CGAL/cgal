// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno
//
#ifndef CGAL_INTERNAL_GENERIC_RANDOM_POINT_GENERATOR_H
#define CGAL_INTERNAL_GENERIC_RANDOM_POINT_GENERATOR_H

#include <CGAL/assertions.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/generators.h>
#include <CGAL/Random.h>
#include <CGAL/property_map.h>


#include <vector>

namespace CGAL {

template <typename Id, class ObjectFromIdMap, class GeneratorOnObject, class P>
class Generic_random_point_generator : public Random_generator_base<P>
{
  typedef Generic_random_point_generator<Id, ObjectFromIdMap, GeneratorOnObject, P> This;
  typedef typename cpp11::result_of<ObjectFromIdMap(Id)>::type                      Geometric_object;

  std::vector<Id> ids;
  std::vector<double> weights;
  ObjectFromIdMap object_from_id_map;
  Random& random;

protected:
  void generate_point();

public:

  template<class InputRange, class ComputeWeight>
  Generic_random_point_generator( const InputRange& input,
                                  const ObjectFromIdMap& object_from_id_map,
                                  const ComputeWeight& compute_weight,
                                  Random& rnd = get_default_random() )
    : Random_generator_base<P>()
    , object_from_id_map(object_from_id_map)
    , random(rnd)
  {
    std::size_t input_size = input.size();
    CGAL_precondition(input_size > 0);

    ids.reserve(input_size);
    weights.reserve(input_size);

    // fill the weights
    double total_weight = 0;
    for(Id id : input)
    {
      //create a geometric object
      Geometric_object object = object_from_id_map(id);
      ids.push_back(id);
      //compute the weight of a face
      total_weight += to_double( compute_weight(object) );
      weights.push_back(total_weight);
    }

    //generate the first point
    generate_point();
  }

  This& operator++()
  {
    generate_point();
    return *this;
  }
  This  operator++(int)
  {
    This tmp = *this;
    ++(*this);
    return tmp;
  }

  double sum_of_weights() const
  {
    if (weights.empty())
      return 0;
    return weights.back();
  }
};

template < typename Id, class ObjectFromIdMap, class GeneratorOnObject, class P >
void Generic_random_point_generator<Id, ObjectFromIdMap,  GeneratorOnObject, P>::generate_point()
{
  //shoot a random value in weights
  std::size_t target = std::distance(
    weights.begin(),
    std::upper_bound(
      weights.begin(),
      weights.end(),
      random.get_double(0, weights.back())
    )
  );

  // generate the points
  GeneratorOnObject pointCreator(object_from_id_map(ids[target]));
  this->d_item = *pointCreator;
}

namespace internal{
template< class Functor >
struct Apply_approx_sqrt
  : public Functor
{
  Apply_approx_sqrt() : Functor() { }
  Apply_approx_sqrt(const Functor& f) : Functor(f) { }

  template <class T>
  typename boost::remove_reference<
             typename cpp11::result_of<Functor(const T&)>::type>::type
  operator()(const T& t) const
  {
    return approximate_sqrt( static_cast<const Functor&>(*this)(t) );
  }
};
} //namespace internal

}//namesape CGAL

#endif // CGAL_INTERNAL_GENERIC_RANDOM_POINT_GENERATOR_H
