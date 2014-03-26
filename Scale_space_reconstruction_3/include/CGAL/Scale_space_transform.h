//A method to transform a point set.
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


#ifndef SCALE_SPACE_PERTURB
#define SCALE_SPACE_PERTURB

#include <omp.h>

#include <map>
#include <vector>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>

#include <Eigen/Dense>

#include "check3264.h"


// Transform a point set by projecting each point onto the best planar fit of the weighted local neighborhood.
template < typename Kernel >
class Scale_space_transform {
	typedef typename Kernel::FT													Scalar;

public:
	typedef typename Kernel::Point_3											Point;	
	typedef typename std::vector<Point>											Collection;

private:
	typedef CGAL::Search_traits_3<Kernel>										Tree_traits;
	
	typedef typename CGAL::Orthogonal_incremental_neighbor_search<Tree_traits>	Neighbor_search;
	typedef typename Neighbor_search::iterator									Neighbor_iterator;
	typedef typename Neighbor_search::Tree										Neighbor_Tree;
	
	typedef typename Kernel::Compare_squared_distance_3							Compare_squared_distance;

	Scalar _radius2;
	Collection _points;

public:
	Scale_space_transform(const Scalar& r): _radius2(r) {}

	void set_radius2(const Scalar& r) {_radius2 = r;}
	Scalar get_radius2() const {return _radius2;}

	template < class InputIterator >
	void assign(InputIterator start, InputIterator end) {
		std::cout << "SST: Assign" << std::endl;
		_points.assign(start, end);
	}

	Collection& get_collection() {return _points;}
	const Collection& get_collection() const {return _points;}

	void iterate() {
		std::cout << "SST: Iterate" << std::endl;
		typedef std::vector<unsigned int>					CountVec;
		typedef typename std::map<Point, size_t>			PIMap;

		typedef Eigen::Matrix<double, 3, Eigen::Dynamic>	Matrix3D;
		typedef Eigen::Array<double, 1, Eigen::Dynamic>		Array1D;
		typedef Eigen::Matrix3d								Matrix3;
		typedef Eigen::Vector3d								Vector3;
		typedef Eigen::SelfAdjointEigenSolver<Matrix3>		EigenSolver;

		typedef _ENV::s_ptr_type							p_size_t;
		
		// This method must be called after filling the point collection.
		CGAL_assertion(!_points.empty());
		Collection points;
		points.swap(_points);
		Scalar radius = _radius2;

		// Construct a search tree of the points.
		Neighbor_Tree tree(points.begin(), points.end());
		if(!tree.is_built())
			tree.build();

		// Collect the number of neighbors of each point.
		// This can be done parallel.
		CountVec neighbors(points.size(), 0);
		Compare_squared_distance compare;
		p_size_t count = points.size(); // openMP can only use signed variable
#pragma omp parallel for shared(count,tree,points,radius,neighbors) firstprivate(compare)
		for (p_size_t i = 0; i < count; ++i) {
			// Iterate over the neighbors until the first one is found that is too far.
			Neighbor_search search(tree, points[i]);
			for (Neighbor_iterator nit = search.begin(); nit != search.end() && compare(points[i], nit->first, radius) != CGAL::LARGER; ++nit)
				++neighbors[i];
		}
		
		// Construct a mapping from each point to its index.
		PIMap indices;
		for (p_size_t i = 0; i < count; ++i)
			indices[points[i]] = i;

		// Compute the tranformed point locations.
		// This can be done parallel.
#pragma omp parallel for shared(count,neighbors,tree,points,radius) firstprivate(compare)
		for (p_size_t i = 0; i < count; ++i) {
			// If the neighborhood is too small, the vertex is not moved.
			if (neighbors[i] < 4)
				continue;

			// Collect the vertices within the ball and their weights.
			Neighbor_search search(tree, points[i]);
			Matrix3D pts(3, neighbors[i]);
			Array1D wts(1, neighbors[i]);
			int column = 0;
			for (Neighbor_iterator nit = search.begin(); nit != search.end() && compare(points[i], nit->first, radius) != CGAL::LARGER; ++nit, ++column) {
				pts(0, column) = CGAL::to_double(nit->first[0]);
				pts(1, column) = CGAL::to_double(nit->first[1]);
				pts(2, column) = CGAL::to_double(nit->first[2]);
				wts(column) = 1.0 / neighbors[indices[nit->first]];
			}

			// Construct the barycenter.
			Vector3 bary = (pts.array().rowwise() * wts).rowwise().sum() / wts.sum();
			
			// Replace the points by their difference with the barycenter.
			pts = (pts.colwise() - bary).array().rowwise() * wts;

			// Construct the weighted covariance matrix.
			Matrix3 covariance = Matrix3::Zero();
			for (column = 0; column < pts.cols(); ++column)
				covariance += wts.matrix()(column) * pts.col(column) * pts.col(column).transpose();

			// Construct the Eigen system.
			EigenSolver solver(covariance);

			// If the Eigen system does not converge, the vertex is not moved.
			if (solver.info() != Eigen::Success)
				continue;

			// Find the Eigen vector with the smallest Eigen value.
			std::ptrdiff_t index;
			solver.eigenvalues().minCoeff(&index);
			if (solver.eigenvectors().col(index).imag() != Vector3::Zero()) {
				// This should actually never happen!
				CGAL_assertion(false);
				continue;
			}
			Vector3 n = solver.eigenvectors().col(index).real();

			// The vertex is moved by projecting it onto the plane
			// through the barycenter and orthogonal to the Eigen vector with smallest Eigen value.
			Vector3 bv = Vector3(CGAL::to_double(points[i][0]), CGAL::to_double(points[i][1]), CGAL::to_double(points[i][2])) - bary;
			Vector3 per = bary + bv - (n.dot(bv) * n);

			points[i] = Point(per(0), per(1), per(2));
		}
		
		_points.swap(points);
	}

	void iterate(unsigned int iterations) {
		for (unsigned int iter = 0; iter < iterations; ++iter)
			iterate();
	}
	
	template < class InputIterator >
	void operator()(InputIterator start, InputIterator end) {
		_points.assign(start, end);
		iterate();
	}

	template < class ForwardIterator >
	void operator()(ForwardIterator start, ForwardIterator end, unsigned int iterations) {
		// Perturb the points for a number of iterations.
		_points.assign(start, end);
		iterate(iterations);
	}
}; // class Scale_space_transform

#endif // SCALE_SPACE_PERTURB