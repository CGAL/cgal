// Copyright (c) 2012 GeometryFactory (France).
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
// Author(s)     : Maxime Gimeno
//
#ifndef GENERIC_RANDOM_POINT_GENERATOR_H
#define GENERIC_RANDOM_POINT_GENERATOR_H

#include <CGAL/generators.h>
#include <CGAL/Random.h>
#include <vector>
#include <boost/foreach.hpp>

namespace CGAL{
template <typename T, class ConstructObject, class GeneratorOnObject, class P>
class Generic_random_generator_on_object: public Random_generator_base<P>
{
  typedef Generic_random_generator_on_object<T, ConstructObject, GeneratorOnObject, P> This;
  typedef typename ConstructObject::reference Geometric_object;
  typedef typename GeneratorOnObject::result_type result_type;
  std::vector<T> faces;
  ConstructObject constructObject;
  //vector of weights
  std::vector<double> weights;
  Random& random;
  void generate_point();
public:
  template<class InputRange, class ComputeObjectWeight>
  Generic_random_generator_on_object(InputRange input, ConstructObject constructObject,
                                     ComputeObjectWeight computeObjectWeight, Random& rnd = get_default_random())
   : Random_generator_base<P>(),
     constructObject(constructObject),
     random(rnd)
  {
    // fill the weights
    double last_weight = 0;
    BOOST_FOREACH(T fit, input)
    {
      //create a face
      Geometric_object face = get(constructObject, fit);
      faces.push_back(fit);
      //compute the weight of a face
      double weight = computeObjectWeight(face);
      last_weight+=weight;
      weights.push_back(last_weight);
    }
    //generate the points
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
};

template < typename T, class ConstructObject, class GeneratorOnObject, class P >
void Generic_random_generator_on_object<T, ConstructObject,  GeneratorOnObject, P>::generate_point()
{
  //shoot a random value in weights
  int target = std::distance(
     weights.begin(),
     std::upper_bound(
      weights.begin(),
      weights.end(),
      random.get_double(0, weights.back())
      )
     );

  // generate the points
  GeneratorOnObject pointCreator(get(constructObject,faces[target]));
  this->d_item = *pointCreator;
}
}//namesape CGAL
#endif // GENERIC_RANDOM_POINT_GENERATOR_H


