//A method that computes the size of a mean neighborhood.
//Copyright (C) 2013  INRIA - Sophia Antipolis
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Author(s):      Thijs van Lankveld


#ifndef MEAN_NEIGHBORHOOD_BALL
#define MEAN_NEIGHBORHOOD_BALL

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Random.h>

// Compute the mean radius of the spheres containing the specified number of nearest neighbors from a point sample.
template < typename Kernel >
class Mean_neighborhood_ball {
	typedef typename Kernel::FT											Scalar;

	typedef typename CGAL::Search_traits_3<Kernel>						Tree_traits;

	typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits>	Neighbor_search;
	typedef typename Neighbor_search::iterator							Neighbor_iterator;
	typedef typename Neighbor_search::Tree								Neighbor_tree;

	typedef CGAL::Random												Random;

	unsigned int handled, checked;
	Scalar dist_far;

	Neighbor_tree tree;

	Random generator;
	
public:
	typedef typename Kernel::Point_3									Point;

	typedef Scalar														result_type;

private:
	unsigned int neighbors, samples;

public:
	Mean_neighborhood_ball(unsigned int nn = 30, unsigned int smp = 200): neighbors(nn), samples(smp) {}

    unsigned int get_neighbors() const { return neighbors; }
    void set_neighbors( unsigned int nn ) { neighbors = nn; }
    unsigned int get_samples() const { return samples; }
    void set_samples( unsigned int smp ) { samples = smp; }
	
	template < class InputIterator >
	void constructTree(InputIterator start, InputIterator end) {
		std::cout << "MNB: tree" << std::endl;
		tree.clear();
		tree.insert(start, end);
		if(!tree.is_built())
			tree.build();
	}

	void handlePoint(const Point& p) {
		Neighbor_search search(tree, p, neighbors+1);

		dist_far += CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(p, (search.end()-1)->first)));
	}

	result_type compute() {
		handled = 0;
		checked = 0;
		dist_far = 0;

		std::cout << "MNB: sample" << std::endl;
		for (typename Neighbor_tree::const_iterator it = tree.begin(); it != tree.end(); ++it) {
			if (generator.get_double() < double(samples - checked) / (tree.size() - handled)) {
				handlePoint(*it);
				++checked;
			}
			++handled;
		}

		dist_far /= double(checked);

		return dist_far;
	}

	template < class InputIterator >
	result_type operator()(InputIterator start, InputIterator end) {
		constructTree(start, end);
		return compute();
	}
}; // Mean_neighborhood_ball

#endif MEAN_NEIGHBORHOOD_BALL
