// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Weighted_Minkowski_distance.h
// package       : ASPAS
// revision      : 1.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

// Note: this implementation is based on a weight vector,
// another traits class for weighted Minkowski distances based on a 
// weight matrix will be added to ASPAS in the future

// Note: Use p=0 to denote the weighted Linf-distance 
// For 0<p<1 Lp is not a metric

#ifndef CGAL_WEIGHTED_MINKOWSKI_DISTANCE_H
#define CGAL_WEIGHTED_MINKOWSKI_DISTANCE_H
#include <math.h>
#include <CGAL/Box.h>

namespace CGAL {

  template <class Tree_traits>
  class Weighted_Minkowski_distance {

    public:

    typedef typename Tree_traits::Item Item;
    typedef typename Tree_traits::NT NT;
    typedef std::vector<NT> Weight_vector;

    private:

    NT p; // p denotes the power

    unsigned int The_dimension;

    Weight_vector *The_weights;

    public:

    // default constructor
    // default case dim=2, L2.
    
    Weighted_Minkowski_distance() : p(2.0), The_dimension(2) 
    { 
         The_weights = new Weight_vector(2);
         The_weights[0]=1.0;
         The_weights[1]=1.0;
    }

	Weighted_Minkowski_distance (float power, int dim,
		Weight_vector Weights) : p(power), The_dimension(dim)
	{
		assert(p >= 0.0);
		assert(The_dimension==Weights.size());
		for (unsigned int i = 0; i < Weights.size(); ++i)
                assert(Weights[i]>=0.0);
		The_weights = new Weight_vector(Weights.size());
		for (unsigned int i = 0; i < Weights.size(); ++i) (*The_weights)[i]=Weights[i];
	}

	~Weighted_Minkowski_distance() {
		delete The_weights;
	};

	inline NT distance(const Item& p1, const Item& p2) {
		NT distance = 0.0;
		if (p == 0.0) {
				for (unsigned int i = 0; i < The_dimension; ++i)
			        if ((*The_weights)[i] * fabs(p1[i] - p2[i]) > distance)
				        distance = (*The_weights)[i] * fabs(p1[i]-p2[i]);
		}
		else
			for (unsigned int i = 0; i < The_dimension; ++i)
				distance += (*The_weights)[i] * pow(fabs(p1[i]-p2[i]),p);
        return distance;
	}


	inline NT lower_bound_distance_to_box(const Item& Point,
					      const Box<NT>& b) {
		NT distance = 0.0;
		if (p == 0.0)
		{
			for (int i = 0; i < b.dimension(); ++i) {
				if ((*The_weights)[i]*(b.lower(i) - Point[i]) > distance)
					distance = (*The_weights)[i] * (b.lower(i)-Point[i]);
				if ((*The_weights)[i] * (Point[i] - b.upper(i)) > distance)
					distance = (*The_weights)[i] * (Point[i]-b.upper(i));
			}
		}
		else
		{
			for (int i = 0; i < b.dimension(); ++i) {
				if (Point[i] < b.lower(i))
					distance += (*The_weights)[i] * pow(b.lower(i)-Point[i],p);
				if (Point[i] > b.upper(i))
					distance += (*The_weights)[i] * pow(Point[i]-b.upper(i),p);
			}
		};
		return distance;
	}

	inline NT upper_bound_distance_to_box(const Item& Point,
					      const Box<NT>& b) {
		NT distance=0.0;
		if (p == 0.0)
		{
			for (int i = 0; i < b.dimension(); ++i) {
				if (Point[i] >= (b.lower(i) + b.upper(i))/2.0)
					if ((*The_weights)[i] * (Point[i] - b.lower(i)) > distance)
						distance = (*The_weights)[i] * (Point[i]-b.lower(i));
				else
					if ((*The_weights)[i] * (b.upper(i) - Point[i]) > distance)
						distance = (*The_weights)[i] * ( b.upper(i)-Point[i]);
			}
		}
		else
		{
			for (int i = 0; i < b.dimension(); ++i) {
				if (Point[i] <= (b.lower(i)+b.upper(i))/2.0)
					distance += (*The_weights)[i] * pow(b.upper(i)-Point[i],p);
				else
					distance += (*The_weights)[i] * pow(Point[i]-b.lower(i),p);
			}
		};
		return distance;
	}

	inline NT new_distance(NT& dist, NT old_off, NT new_off,
			int cutting_dimension)  {
		NT new_dist;
		if (p == 0.0)
		{
			if ((*The_weights)[cutting_dimension]*fabs(new_off) > dist) 
				new_dist= (*The_weights)[cutting_dimension]*fabs(new_off);
			else new_dist=dist;
		}
		else
		{
			new_dist = dist + (*The_weights)[cutting_dimension] * 
					(pow(fabs(new_off),p)-pow(fabs(old_off),p));
		}
                return new_dist;
	}

  inline NT transformed_distance(NT d) {

		if (p <= 0.0) return d;
		else return pow(d,p);

	}

  inline NT inverse_of_transformed_distance(NT d) {

		if (p <= 0.0) return d;
		else return pow(d,1/p);

	}

  }; // class Weighted_Minkowski_distance

} // namespace CGAL
#endif // WEIGHTED_MINKOWSKI_DISTANCE_H
