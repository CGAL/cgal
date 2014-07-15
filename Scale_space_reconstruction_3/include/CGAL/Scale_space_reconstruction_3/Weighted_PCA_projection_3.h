//Projection for a point onto its weighted PCA plane.
//Copyright (C) 2014  INRIA - Sophia Antipolis
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


#ifndef CGAL_WEIGHTED_PCA_PROJECTION_3_H
#define CGAL_WEIGHTED_PCA_PROJECTION_3_H

#include <Eigen/Dense>


namespace CGAL {
    
/// projects a point orthogonally onto a weighted least-squares planar approximation of a point set.
/** \ingroup PkgScaleSpaceReconstruction3Classes
 *  
 *  This class uses the eigenvector solvers of \ref thirdpartyEigen.
 *
 *  Version 3.1.2 (or greater) of \ref thirdpartyEigen must be available on the
 *  system.
 *
 *  \tparam GeomTraits the geometric traits class of the input and output. It
 *  must be a model of Kernel with a RealEmbeddable FieldNumberType.
 *
 *  \note Irrespective of the geometric traits class, the projection is
 *  always estimated up to double precision.
 *
 *  \cgalModels `WeightedPCAProjection_3`
 */
template < class GeomTraits >
class Weighted_PCA_projection_3 {
public:
	typedef typename GeomTraits::FT         FT;
	typedef typename GeomTraits::Point_3    Point;

private:
    typedef Eigen::Matrix<double, 3, Eigen::Dynamic>	Matrix3D;       // 3-by-dynamic double-value matrix.
    typedef Eigen::Array<double, 1, Eigen::Dynamic>		Array1D;        // 1-by-dynamic double-value array.
    typedef Eigen::Matrix3d								Matrix3;        // 3-by-3 double-value matrix.
    typedef Eigen::Vector3d								Vector3;        // 3-by-1 double-value matrix.
    typedef Eigen::SelfAdjointEigenSolver< Matrix3 >    EigenSolver;    // Eigen value decomposition solver.

private:
    bool _comp;        // Whether the eigenvalues are computed.

    Matrix3D _pts;  // points.
    Array1D _wts;   // weights.

    Vector3 _bary;  // barycenter.
    Vector3 _norm;  // normal.

public:
/// \name Constructors
/// {
    /// constructs an empty projection to hold the points.
    /** \param size is the number of points that will be added.
     */
    Weighted_PCA_projection_3( unsigned int size ): _comp(false), _pts(3,size), _wts(1,size) {}

    /// constructs the weighted least-squares planar approximation of a point set.
    /** Similar to constructing an empty projection and calling
     *  <code>[set_points(points_begin, points_end, weights_begin, weights_end)](\ref WeightedPCAProjection_3::set_points )</code>
     *
     *  \param points_begin is an iterator to the first point.
     *  \param points_end is a past-the-end oterator for the points.
     *  \param weights_begin is an iterator to the weight of the first point.
     */
    template < typename PointIterator, typename WeightIterator >
    Weighted_PCA_projection_3( PointIterator points_begin, PointIterator points_end, WeightIterator weights_begin )
    : _comp(false) {
        std::size_t size = std::distance( points_begin, points_end );
        _pts = Matrix3D(3,size);
        _wts = Array1D(1,size);
        set_points( points_begin, points_end, weights_begin );
    }

/// \}

public:
    // sets a weighted point in the collection.
    void set_point( unsigned int i, const Point& p, const FT& w ) {
        _pts( 0, i ) = CGAL::to_double( p[0] );
        _pts( 1, i ) = CGAL::to_double( p[1] );
        _pts( 2, i ) = CGAL::to_double( p[2] );
        _wts( i ) = CGAL::to_double( w );
    }

    // sets the weighted points and compute PCA.
    template < typename PointIterator, typename WeightIterator >
    bool set_points( PointIterator points_begin, PointIterator points_end, WeightIterator weights_begin ) {
        for( unsigned int column = 0; points_begin != points_end; ++column, ++points_begin, ++weights_begin ) {
            _pts( 0, column ) = CGAL::to_double( (*points_begin)[0] );
            _pts( 1, column ) = CGAL::to_double( (*points_begin)[1] );
            _pts( 2, column ) = CGAL::to_double( (*points_begin)[2] );
            _wts( column ) = CGAL::to_double( *weights_begin );
        }
        approximate();
    }
    
    // compute weighted PCA.
    bool approximate() {
        // Construct the barycenter.
        _bary = ( _pts.array().rowwise() * _wts ).rowwise().sum() / _wts.sum();
			
        // Replace the points by their difference with the barycenter.
        _pts = ( _pts.colwise() - _bary ).array().rowwise() * _wts;

        // Construct the weighted covariance matrix.
        Matrix3 covar = Matrix3::Zero();
        for( unsigned int column = 0; column < _pts.cols(); ++column )
            covar += _wts.matrix()(column) * _pts.col(column) * _pts.col(column).transpose();

        // Construct the Eigen system.
        EigenSolver solver( covar );

        // If the Eigen system does not converge, the vertex is not moved.
        if( solver.info() != Eigen::Success )
            return false;

        // Find the Eigen vector with the smallest Eigen value.
        std::ptrdiff_t index;
        solver.eigenvalues().minCoeff( &index );
        if( solver.eigenvectors().col( index ).imag() != Vector3::Zero() ) {
            // This should actually never happen,
            // because the covariance matrix is symmetric!
            CGAL_assertion( false );
            return false;
        }

        // The normal is the Eigen vector with smallest Eigen value.
        _norm = solver.eigenvectors().col( index ).real();
        _comp = true;
        return true;
    }
    
    // checks whether the weighted least-squares approximating plane has been computed.
    bool is_approximated() const { return _comp; }
    
public:
    // projects a point onto the weighted PCA plane.
    Point project( const Point& p ) {
        CGAL_assertion( _comp );
        // The point is moved by projecting it onto the plane through the
        // barycenter and orthogonal to the normal.
        Vector3 to_p = Vector3( to_double( p[0] ),
                                to_double( p[1] ),
                                to_double( p[2] ) ) - _bary;
        Vector3 proj = _bary + to_p - ( _norm.dot(to_p) * _norm );
        return Point( proj(0), proj(1), proj(2) );
    }
}; // class Weighted_PCA_projection_3

} // namespace CGAL

#endif // CGAL_WEIGHTED_PCA_PROJECTION_3_H