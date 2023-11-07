#ifndef QEM_METRIC_H
#define QEM_METRIC_H

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
#include <CGAL/Aff_transformation_3.h>
#include "types.h"

typedef CGAL::Aff_transformation_3<Kernel>  Aff_transformation;
//template <class Kernel>
class QEM_metric
{
public:
    typedef Eigen::Matrix<double, 10, 1>        EVec10d;
    typedef Eigen::VectorXd                     VectorXd;
    typedef Eigen::MatrixXd                     MatrixXd;
    typedef Eigen::Matrix<double, 3, 3>         Matrix3d;
    typedef Eigen::Matrix<double, 4, 4>         Matrix4d;

private: // members

    EVec10d m_tensor;

    /*
    n n^T | -nnq^T
    -nnq    | (nq)^2

    aa ab ac ad
    ab bb bc bd
    ac bc cc cd
    ad bd cd dd

    m_tensor[0]: aa, m_tensor[1]: ab, m_tensor[2]: ac, m_tensor[3]: ad
    m_tensor[4]: bb, m_tensor[5]: bc, m_tensor[6]: bd, m_tensor[7]: cc
    m_tensor[8]: cd, m_tensor[9]: dd

    m_tensor[0]  m_tensor[1]  m_tensor[2]  m_tensor[3]
    m_tensor[1]  m_tensor[4]  m_tensor[5]  m_tensor[6]
    m_tensor[2]  m_tensor[5]  m_tensor[7]  m_tensor[8]
    m_tensor[3]  m_tensor[6]  m_tensor[8]  m_tensor[9]
    */

public: // functions

    QEM_metric(): m_tensor(EVec10d::Zero()) {}

    EVec10d& tensors() { return m_tensor; }
    const EVec10d& tensors() const { return m_tensor; }

    /*MatrixXd 4x4_matrix()
    {
        MatrixXd qem(4, 4);

        qem <<  m_tensor[0], m_tensor[1], m_tensor[2], m_tensor[3],
                m_tensor[1], m_tensor[4], m_tensor[5], m_tensor[6],
                m_tensor[2], m_tensor[5], m_tensor[7], m_tensor[8],
                m_tensor[3], m_tensor[6], m_tensor[8], m_tensor[9];

        return qem;
    }*/

    Matrix4d get_4x4_svd_matrix()
    {
        MatrixXd qem(4, 4);

        qem <<  m_tensor[0], m_tensor[1], m_tensor[2], m_tensor[3],
                m_tensor[1], m_tensor[4], m_tensor[5], m_tensor[6],
                m_tensor[2], m_tensor[5], m_tensor[7], m_tensor[8],
                0.         , 0.         , 0.         , 1.;

        return qem;
    }

    Matrix4d get_4x4_matrix()
    {
        MatrixXd qem(4, 4);

        qem <<  m_tensor[0], m_tensor[1], m_tensor[2], m_tensor[3],
                m_tensor[1], m_tensor[4], m_tensor[5], m_tensor[6],
                m_tensor[2], m_tensor[5], m_tensor[7], m_tensor[8],
                m_tensor[3], m_tensor[6], m_tensor[8], m_tensor[9];

        return qem;
    }

    Matrix3d get_3x3_matrix()
    {
        MatrixXd qem(3, 3);

        qem <<  m_tensor[0], m_tensor[1], m_tensor[2],
                m_tensor[1], m_tensor[4], m_tensor[5],
                m_tensor[2], m_tensor[5], m_tensor[7];

        return qem;
    }
    /// @brief Compute the tensor using a point and its normal and weighted by the area
    /// @param area 
    /// @param query 
    /// @param normal 
    void init_qem_metrics_face(const double& area, const Point& query, const Vector& normal)
    {
        double a = normal.x();
        double b = normal.y();
        double c = normal.z();
        double d = -1. * CGAL::scalar_product(normal, (query - CGAL::ORIGIN));

        m_tensor[0] = a * a;
        m_tensor[1] = a * b;
        m_tensor[2] = a * c;
        m_tensor[3] = a * d;
        m_tensor[4] = b * b;
        m_tensor[5] = b * c;
        m_tensor[6] = b * d;
        m_tensor[7] = c * c;
        m_tensor[8] = c * d;
        m_tensor[9] = d * d;

        m_tensor = m_tensor * area;
    }
    /// @brief Compute the tensor (probabilistic qem) using a point and its normal and weighted by the area
    /// @param area 
    /// @param mean_query 
    /// @param std_query 
    /// @param mean_normal 
    /// @param std_normal 
    void init_proba_qem_metrics_face(const double area,
    const Point& mean_query,
    const double std_query,
    const Vector& mean_normal,
    const double std_normal)
    {
        double sn2 = std_normal * std_normal;
        double sq2 = std_query * std_query;
        double dot_nq = CGAL::scalar_product(mean_normal, (mean_query - CGAL::ORIGIN));

        m_tensor[0] = mean_normal.x() * mean_normal.x() + sn2;
        m_tensor[1] = mean_normal.x() * mean_normal.y();
        m_tensor[2] = mean_normal.x() * mean_normal.z();
        m_tensor[4] = mean_normal.y() * mean_normal.y() + sn2;
        m_tensor[5] = mean_normal.y() * mean_normal.z();
        m_tensor[7] = mean_normal.z() * mean_normal.z() + sn2;

        Vector vec_nnq = mean_normal * dot_nq + (mean_query - CGAL::ORIGIN) * sn2;
        m_tensor[3] = -vec_nnq.x();
        m_tensor[6] = -vec_nnq.y();
        m_tensor[8] = -vec_nnq.z();

        m_tensor[9] = dot_nq * dot_nq + sn2 * CGAL::scalar_product(mean_query - CGAL::ORIGIN, mean_query - CGAL::ORIGIN) + 
                      sq2 * CGAL::scalar_product(mean_normal, mean_normal) + 3 * sq2 * sn2;

        m_tensor = m_tensor * area;
    }
    /// @brief Compute the qem tensor for a query point and a list of areas and normals
    /// @param query 
    /// @param areas 
    /// @param normals 
    void init_qem_metrics_vertex(Point& query,
        const std::vector<double>& areas, 
        const std::vector<Vector>& normals)
    {
        assert(areas.size() == normals.size());

        for(int i = 0; i < areas.size(); i++)
        {
            Vector normal = normals[i];
            double area = areas[i];

            double a = normal.x();
            double b = normal.y();
            double c = normal.z();
            double d = -1. * CGAL::scalar_product(normal, (query - CGAL::ORIGIN));

            m_tensor[0] += area * a * a;
            m_tensor[1] += area * a * b;
            m_tensor[2] += area * a * c;
            m_tensor[3] += area * a * d;
            m_tensor[4] += area * b * b;
            m_tensor[5] += area * b * c;
            m_tensor[6] += area * b * d;
            m_tensor[7] += area * c * c;
            m_tensor[8] += area * c * d;
            m_tensor[9] += area * d * d;
        }
    }
    /// @brief Compute the affine transformation on a sphere to display the qem error as en ellipsoid
    /// @return the affine transformation
    Aff_transformation aff_transform_sphere()
    {
        
        MatrixXd A(3, 3);
        A <<    m_tensor[0], m_tensor[1], m_tensor[2], 
                m_tensor[1], m_tensor[4], m_tensor[5], 
                m_tensor[2], m_tensor[5], m_tensor[7];

        Eigen::EigenSolver<MatrixXd> es(A);
        // eigenvalues
        VectorXd D = es.eigenvalues().real();
        D = D.cwiseMax(1e-5).cwiseInverse(); // take inverse
        MatrixXd D_inv = MatrixXd::Zero(3, 3);
        D_inv.diagonal() = D / D.sum(); // normalization
        // eigenvectors
        MatrixXd V = es.eigenvectors().real();
        // new matrix
        MatrixXd new_trans = V * D_inv * V.transpose();
        Aff_transformation aff( new_trans(0, 0), new_trans(0, 1), new_trans(0, 2),
                                new_trans(1, 0), new_trans(1, 1), new_trans(1, 2),
                                new_trans(2, 0), new_trans(2, 1), new_trans(2, 2));

        return aff;
    }

    QEM_metric& operator+(const QEM_metric& other)
    {
        this->tensors() = this->tensors() + other.tensors();
        return *this;
    }

    QEM_metric& operator-(const QEM_metric& other)
    {
        this->tensors() = this->tensors() - other.tensors();
        return *this;
    }

    QEM_metric& operator*(const double& scale)
    {
        this->tensors() = this->tensors() * scale;
        return *this;
    }
    /// @brief Compute optimal point using either SVD or the direct inverse
    /// @param cluster_qem 
    /// @param cluster_pole 
    /// @return the optimal point
    Point compute_optimal_point(QEM_metric& cluster_qem, Point& cluster_pole)
    {
        // solve Qx = b
        Eigen::MatrixXd qem_mat = cluster_qem.get_4x4_svd_matrix();
        Eigen::VectorXd qem_vec = qem_mat.row(3); // 0., 0., 0., 1.
        Eigen::VectorXd optim(4);

        // check rank
        Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(qem_mat);
        lu_decomp.setThreshold(1e-5);

        // full rank -> direct inverse
        if(lu_decomp.isInvertible())
        {
            optim = lu_decomp.inverse() * qem_vec;
        }
        else
        {   // low rank -> svd pseudo-inverse
            Eigen::JacobiSVD<Eigen::MatrixXd> svd_decomp(qem_mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
            svd_decomp.setThreshold(1e-5);

            optim(0) = cluster_pole.x();
            optim(1) = cluster_pole.y();
            optim(2) = cluster_pole.z();
            optim(3) = 1.;

            optim = optim + svd_decomp.solve(qem_vec - qem_mat * optim);
        }

        Point optim_point(optim(0), optim(1), optim(2));

        return optim_point;
    }

};




#endif