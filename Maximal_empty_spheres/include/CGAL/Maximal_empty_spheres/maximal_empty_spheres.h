// Copyright (c) 2025  TU Berlin
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Max Kohlbrenner

#ifndef CGAL_MAXIMAL_EMPTY_SPHERES_MAXIMAL_EMPTY_SPHERES_H
#define CGAL_MAXIMAL_EMPTY_SPHERES_MAXIMAL_EMPTY_SPHERES_H

// #include <CGAL/license/Maximal_empty_spheres.h>

#include <CGAL/Maximal_empty_spheres/internal/Lie_geometry.h>
#include <CGAL/Maximal_empty_spheres/internal/rotation_from_to.h>

#include <iostream>
#include <Eigen/Core>

#include <CGAL/Epick_d.h>
#include <CGAL/Triangulation.h>
#include <CGAL/algorithm.h>
#include <CGAL/Dimension.h>

namespace CGAL {

template<typename Dimension>
void maximal_empty_spheres(const Eigen::MatrixXd &G, Eigen::MatrixXd &result, Eigen::MatrixXi *contact_indices=NULL, double atol=1e-8, int debug_level=0) {

    bool full_simplices_only=true;
    const int D = Dimension::value;
    assert( D == G.cols()-1 && "Dimension does not match the number of columns in G");

    typedef CGAL::Epick_d< CGAL::Dimension_tag<D+2> >  Kernel;
    typedef CGAL::Triangulation_vertex<Kernel, std::ptrdiff_t>                                      Vertex;
    typedef CGAL::Triangulation_ds_full_cell<void,CGAL::TDS_full_cell_mirror_storage_policy>        DS_full_cell;
    typedef CGAL::Triangulation_full_cell<Kernel, int, DS_full_cell>                                Full_cell;
    typedef CGAL::Triangulation_data_structure<CGAL::Dimension_tag<D+2>,
                                               Vertex,
                                               Full_cell>                                           TDS;
    typedef CGAL::Triangulation<Kernel,TDS>                                                         Triangulation;
    typedef typename Triangulation::Vertex_handle                                                   Vertex_handle;

    using namespace CGAL::Maximal_empty_spheres::internal;

    std::cout << "Loaded SDF values of shape (" << G.rows() << ", " << G.cols() << ")"  << std::endl;

    Eigen::MatrixXd SE = G;
    SE.col(D) = G.col(D).array().abs();

    // std::cout << "SE: " << std::endl;
    // std::cout << SE << std::endl;

    std::cout << "Convert to Lie representation" << std::endl;
    Eigen::MatrixXd H;
    lie_ip_matrix(D,H);

    Eigen::MatrixXd SL;
    spheres_to_lie(SE,SL);
    // std::cout << "Lie Repr: " << std::endl << SL << std::endl;

    Eigen::MatrixXd NC = SL * H;

    // std::cout << "... on quadric? ";
    // std::cout << std::boolalpha << ( ((NC.array() * SL.array()).rowwise().sum()).matrix().norm() < atol) << std::endl;

    std::cout << std::boolalpha;

    Eigen::VectorXd s0  = Eigen::VectorXd::Zero(D+3);
    s0(D+2) = -1.;
    s0(D)   = -1.;
    Eigen::VectorXd dst = Eigen::VectorXd::Zero(D+3);
    dst[D+2] = 1.;
    Eigen::MatrixXd R = rotation_from_to(s0, dst);

    std::cout << "All normals in the same halfspace? : " <<  (((NC*s0).transpose()).array() >= atol).all() << std::endl;

    Eigen::MatrixXd NCR = NC * R.transpose(); // rotate
    std::cout << "Last component greater than zero?  : "  << (NCR.col(D+2).array()    >= atol).all() << std::endl;
    NCR.array().colwise() /= NCR.array().col(D+2);
    Eigen::MatrixXd Np = NCR.block(0,0,NCR.rows(),D+2);

    // ------------------------ CONVEX HULL COMPUTATION -------------------------

    // --> convert to point list --
    std::vector<typename Triangulation::Point> points;
    points.reserve(int(Np.rows()));
    std::vector<typename Kernel::FT> p_;
    p_.reserve(D+2);
    for (int i=0; i<Np.rows(); i++){
        p_.clear();
        for (int d=0; d<D+2; d++){
            p_.push_back(Np(i,d));
        }
        points.push_back(typename Triangulation::Point(p_.begin(),p_.end()));
    }

    Triangulation t(D+2);
    t.insert_and_index(points.begin(), points.end());
    // --> recover the infinite cells (contain the convex hull)
    std::vector<typename Triangulation::Full_cell_handle> infinite_cells;
    for(auto it = t.full_cells_begin(); it != t.full_cells_end(); ++it) {
        if(t.is_infinite(it)) {
            infinite_cells.push_back(it);
        }
    }

    // --> get a unique index per cell
    int ci=0;
    for(auto ch : infinite_cells) ch->data() = ci++;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    Eigen::MatrixXd Ks(infinite_cells.size(),D+3);
    Eigen::Vector<bool,Eigen::Dynamic> full_simplices(infinite_cells.size());

    Eigen::MatrixXd simplex(D+2,D+3);
    Eigen::RowVectorXd simplex_row(D+3);
    int ki=0;
    for(auto ch : infinite_cells) {

        int si=0;
        for (int i=0; i<D+3; i++) {
             if(ch->vertex(i) != t.infinite_vertex()) {
                 // std::cout << point_indices[ch->vertex(i)] << std::endl;
                 simplex.row(si++) = NC.row(ch->vertex(i)->data());
             }
        }
        // std::cout << "simplex: " << std::endl << simplex << std::endl;
        svd.compute(simplex, Eigen::ComputeFullV);
        Ks.row(ki)         = svd.matrixV().col(D+2).transpose();
        full_simplices(ki) = svd.singularValues().array().abs().minCoeff() > atol;
        ki++;

        /*
        if (svd.singularValues().array().abs().minCoeff() <= atol) {
            std::cout << "Not a full simplex" << std::endl;
        }
        */
    }

    // std::cout << "Ks.shape: (" << Ks.rows() << ", " << Ks.cols() << ")" << std::endl;
    // std::cout << "Ks: " << std::endl << Ks << std::endl;
    // std::cout << "Ks * NC.T" << std::endl << Ks * NC.transpose() << std::endl;

    int n_check_planes = (NC.rows()<=1000)? NC.rows(): 1000;
    // std::cout << "n_check_planes: " << n_check_planes << std::endl;
    std::cout << "NC.shape(): " << NC.rows() << ", " << NC.cols() << std::endl;
    Eigen::Vector<bool,Eigen::Dynamic> inv = ((Ks * NC.block(0,0,n_check_planes,D+3).transpose()).rowwise().maxCoeff().array() >= atol); // .rowwise().maxCoeff().array() >= 1e-5) << std::endl;
    for (int i=0; i<Ks.rows(); i++) if (inv(i)) Ks.row(i) *= -1;

    if (debug_level > 3) {
        // costly and not possible for large sets of spheres
        std::cout << "Debugging, remove later (matrices too large), all Ks in the cone?: ";
        std::cout << ((NC * Ks.transpose()).rowwise().maxCoeff().array() <= atol).array().all() << std::endl;
    }

    std::vector<Eigen::RowVectorXd> solutions_;
    std::vector<Eigen::RowVectorXi> contact_indices_;
    for(auto ch : infinite_cells) {
        int ci = ch->data();
        for (int i=0; i<D+3; i++) {
             if(ch->vertex(i) != t.infinite_vertex()) {
                int cj =  ch->neighbor(i)->data();
                if ((ci < cj) &&    // only treat each edge once
                    (!full_simplices_only || (full_simplices(ci) && full_simplices(cj)))  // only consider edges between full simplices
                    ) {
                    double l1,l2;
                    line_quadric_intersection(Ks.row(ci), Ks.row(cj), H, l1,l2);
                    if ((0.<=l1) && (l1<=1.)){
                        solutions_.push_back((1-l1)*Ks.row(ci)+l1*Ks.row(cj));
                        if (contact_indices){
                            Eigen::RowVectorXi st(D+1);
                            int ni=0;
                            for (int j=0; j<D+2; j++){
                                if (ch->vertex((i+1+j)%(D+3)) != t.infinite_vertex())
                                    st(ni++) = ch->vertex((i+1+j)%(D+3))->data();
                            }
                            contact_indices_.emplace_back(st);
                        }
                    }
                    if ((0.<=l2) && (l2<=1.)){
                        solutions_.push_back((1-l2)*Ks.row(ci)+l2*Ks.row(cj));

                        if (contact_indices){
                            Eigen::RowVectorXi st(D+1);
                            int ni=0;
                            for (int j=0; j<D+2; j++){
                                if (ch->vertex((i+1+j)%(D+3)) != t.infinite_vertex())
                                    st(ni++) = ch->vertex((i+1+j)%(D+3))->data();
                            }
                            contact_indices_.emplace_back(st);
                        }
                    }
                }
             }
        }
    }

    Eigen::MatrixXd solutions(solutions_.size(),D+3);
    for (int i=0; i<solutions_.size(); i++) solutions.row(i) = solutions_[i];
    std::cout << "Solutions.size: " << solutions.rows() << ", " << solutions.cols() << std::endl;

    if (debug_level > 3) {
        std::cout << "Debugging, remove later (matrices too large), all solutions in the cone?: ";
        std::cout << ((NC * solutions.transpose()).rowwise().maxCoeff().array() <= atol).array().all() << std::endl;
    }
    // ------- Post selection ------------
    Eigen::Vector<bool,Eigen::Dynamic> neg_convention = (solutions.col(D+2).array() < 0.);
    Eigen::Vector<bool,Eigen::Dynamic> neg_radius     = (solutions.col(D+1).array().sign() != solutions.col(D+2).array().sign());

    std::vector<Eigen::RowVectorXd> solutions_filtered_;
    std::vector<Eigen::RowVectorXi> contact_indices_filtered_;
    for (int i=0; i<solutions.rows(); i++){
        if (neg_convention(i) && neg_radius(i)){
            solutions_filtered_.push_back(solutions.row(i));
            if (contact_indices){
                contact_indices_filtered_.push_back(contact_indices_[i]);            
            }
        }
    }

    Eigen::MatrixXd solutions_filtered(solutions_filtered_.size(),D+3);
    for (int i=0; i<solutions_filtered_.size(); i++) solutions_filtered.row(i) = solutions_filtered_[i];
    std::cout << "solutions_filtered.shape: (" << solutions_filtered.rows() << ", " << solutions_filtered.cols() << ")"  << std::endl;

    if (contact_indices){
        contact_indices->resize(contact_indices_filtered_.size(),D+1);
        for (int i=0; i<contact_indices_filtered_.size(); i++) contact_indices->block(i,0,1,D+1) = contact_indices_filtered_[i];
    }

    lie_to_spheres(solutions_filtered, result);
    }

    template<typename SphereRange, typename OutputIterator, class DimensionTag>
    void maximal_empty_spheres(const SphereRange& /* input*/ , OutputIterator /* result*/, DimensionTag /* tag */ )
    {
       std::cout << "should not match any specialization" << std::endl;
    }

    /*!
     * \ingroup PkgMaximalEmptySpheresFunctions
     * \brief A dummy function to demonstrate the structure of the library.
     *
     * This function is a placeholder and does not perform any operations.
     * It is intended to be replaced with actual functionality in the future.
     */
    template<typename SphereRange, typename OutputIterator>
    void maximal_empty_spheres(const SphereRange& input, OutputIterator result, CGAL::Dimension_tag<3> /* tag */) {
        using Sphere_3 = typename SphereRange::value_type;
        using Kernel = typename Kernel_traits<Sphere_3>::Kernel;
        using Point_3 = typename Kernel::Point_3;

        Eigen::MatrixXd G(input.size(),4), res;
        for(auto it = input.begin(); it != input.end(); ++it) {
            G.row(it - input.begin()) = Eigen::RowVector4d(it->center().x(), it->center().y(), it->center().z(), sqrt(it->squared_radius()));
        }
        maximal_empty_spheres<Dimension_tag<3>>(G, res);
        for(int i=0; i<res.rows(); i++) {
          Point_3 p(res(i,0), res(i,1), res(i,2));
          *result++ = Sphere_3(p, res(i,3)*res(i,3)); // res(i,3) is the radius
        }

        const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        // write soution in lie representation:  AF: I think this is not Lie representation, but rather the sphere representation
        std::ofstream pointfile("solution_lie_3.csv");
        pointfile << res.format(CSVFormat);
    }

    template<typename SphereRange, typename OutputIterator>
    void maximal_empty_spheres(const SphereRange& input, OutputIterator result, CGAL::Dimension_tag<2> /* tag */) {
        using Circle_2 = typename SphereRange::value_type;
        using Kernel = typename Kernel_traits<Circle_2>::Kernel;
        using Point_2 = typename Kernel::Point_2;

        Eigen::MatrixXd G(input.size(),3), res;
        for(auto it = input.begin(); it != input.end(); ++it) {
            G.row(it - input.begin()) = Eigen::RowVector3d(it->center().x(), it->center().y(), sqrt(it->squared_radius()));
        }
        maximal_empty_spheres<Dimension_tag<2>>(G, res);
        for(int i=0; i<res.rows(); i++) {
          Point_2 p(res(i,0), res(i,1));
          *result++ = Circle_2(p, res(i,2)*res(i,2)); // res(i,2) is the radius
        }

        const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        // write soution in lie representation:
        std::ofstream pointfile("solution_lie_2.csv");
        pointfile << res.format(CSVFormat);
    }

    template<typename SphereRange, typename OutputIterator>
    void maximal_empty_spheres(const SphereRange& input, OutputIterator result) {
        using Sphere = typename SphereRange::value_type;
        constexpr int D = Ambient_dimension<Sphere>::value;
        maximal_empty_spheres(input, result, CGAL::Dimension_tag<D>());
    }
} // namespace CGAL

#endif // GAL_MAXIMAL_EMPTY_SPHERES_MAXIMAL_EMPTY_SPHERES_H
