#ifndef CGAL_EFFICIENT_RANSAC_TYPES_H
#define CGAL_EFFICIENT_RANSAC_TYPES_H

#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/bind/make_adaptable.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <functional>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <random>

namespace CGAL {

  namespace Efficient_ransac {

    unsigned int getRandomInt() {
      static std::mt19937 rng(23);
      return rng();
    }

    class Option_RANSAC
    {
    public:
      Option_RANSAC(): m_probability(0.001f), m_minNbPoints(10), m_epsilon(0.01f), m_normalThresh(0.90f), m_bitmapEpsilon(0.01f)   {};

      //variable
    public:
      float m_probability;//parameter for stopping generating candidate
      unsigned int m_minNbPoints;//parameter for minimum nb of points for a candidate
      float m_epsilon;  //distance gamma-band around the primitive
      float m_normalThresh;	  // this is the cos of the maximal normal deviation
      float m_bitmapEpsilon;	    //for cc, not used now
    };
  }
}

#endif
