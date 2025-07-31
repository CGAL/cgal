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

#ifndef CGAL_MAXIMAL_EMPTY_SPHERES_CONTACT_POINTS_FROM_SIGNED_DISTANCE_H
#define CGAL_MAXIMAL_EMPTY_SPHERES_CONTACT_POINTS_FROM_SIGNED_DISTANCE_H

#include <CGAL/license/Maximal_empty_spheres.h>


#include <CGAL/Maximal_empty_spheres/maximal_empty_spheres.h>
#include <CGAL/Dimension.h>
#include <CGAL/Bbox_3.h>

#include <Eigen/Core>

#include <utility>
#include <vector>
#include <iostream>

namespace CGAL {

template <typename OutputIterator>
OutputIterator contact_points(const Eigen::MatrixXd &G, OutputIterator Pwns, Eigen::MatrixXd *bbxl=nullptr, int debug_level=0)
{
    Eigen::MatrixXi contact_indices;
    Eigen::MatrixXd res;
    Eigen::MatrixXd G_ = G;
    G_.col(3) = G.col(3).array().abs();
    double tolerance = 1e-8;
    CGAL::maximal_empty_spheres<CGAL::Dimension_tag<3>>(G_, res, &contact_indices, tolerance, debug_level);

    // std::cout << res << std::endl;
    // std::cout << "contact_indices:" << std::endl;
    // std::cout << contact_indices << std::endl;

    // contact_indices contains the indices G and contains the spheres that each spheres of res is adjacent to.
    // calculate the contact sphere with largest (absolute) radius for each sphere in G:
    // (careful: results spheres all have negative radius)
    Eigen::VectorXi contact_point_indices = Eigen::VectorXi::Constant(G.rows(),-1);
    Eigen::VectorXd contact_point_radii   = Eigen::VectorXd::Constant(G.rows(),-1.);
    for (int i=0; i<contact_indices.rows(); i++){
        double r = fabs(res(i,3)); // contact spheres have negative radius, use abs value here
        if ((!bbxl)
            ||
            (  (res.block(i,0,1,3).array() >  bbxl->row(0).array()).all()
            && (res.block(i,0,1,3).array() <= bbxl->row(1).array()).all())){

            for (int j=0; j<contact_indices.cols(); j++){
                int n = contact_indices(i,j);
                if ( (contact_point_indices(n) < 0) || ((contact_point_indices(n) >= 0) && (contact_point_radii(n) < r)) ){
                    contact_point_indices(n) = i;
                    contact_point_radii(n)   = r;
                }
            }
        }
    }

    // std::cout << "CP Indices: " << std::endl;
    // std::cout << contact_point_indices.transpose() << std::endl;

    // std::cout << "CP Radii: " << std::endl;
    // std::cout << contact_point_radii.transpose() << std::endl;

    // calculate the contact points and normals
    for (int i=0; i<G.rows(); i++){
        if (contact_point_indices[i] >= 0){
            Point Csdf = Point(G(i,0),G(i,1),G(i,2));
            double rsdf = G(i,3);

            Point Ccontact = Point(
                        res(contact_point_indices[i],0),
                        res(contact_point_indices[i],1),
                        res(contact_point_indices[i],2)
                        );
            Vector D = Ccontact-Csdf;
            double vl = sqrt(D.squared_length());
            // std::cout << "vl-(r+rc): " << vl - (fabs(rsdf)+fabs(res(contact_point_indices[i],3))) << std::endl;
            D /= vl;
            // std::cout << "D.squared_lenght(): " << D.squared_length() << std::endl;
            Point  P = Csdf + fabs(rsdf)*D;
            Vector N = (rsdf >= 0)? -D: D;
            *Pwns++ = std::make_pair(P,N);
        }
    }
    return Pwns;
}


  template <typename Sphere, typename  OutputIterator>
  OutputIterator contact_points_from_signed_distances(const std::vector<std::pair<Sphere,int>>& input_spheres,
                                                      OutputIterator Pwns,
                                                      bool filter_contact_spheres_bbx = true,
                                                      int debug_level=0)
  {
    // Eigen::MatrixXd G(input_spheres.size(), 4);
    Bbox_3 bb;
    int np = 0, nn = 0;
    for (const auto& pair : input_spheres) {
      bb += pair.first.center().bbox();
      if(pair.second >= 0) { // positive radius means outside
        ++np;
      }else{
        ++nn;
      }
      }
    Eigen::MatrixXd bbxl(2,3);
    bbxl.row(0) = Eigen::RowVector3d(bb.xmin(), bb.ymin(), bb.zmin());
    bbxl.row(1) = Eigen::RowVector3d(bb.xmax(), bb.ymax(), bb.zmax());
    Eigen::MatrixXd Gp(np,4);
    Eigen::MatrixXd Gn(nn,4);
    np=0;
    nn=0;
    for (int i=0; i<input_spheres.size(); i++){
      const std::pair<Sphere,int>& sphere_and_orientation = input_spheres[i];
      const Sphere& sphere = sphere_and_orientation.first;
        if (sphere_and_orientation.second > 0){
            Gp.row(np++) = Eigen::RowVector4d(sphere.center().x(), sphere.center().y(), sphere.center().z(), sqrt(sphere.squared_radius()));
        } else {
            Gn.row(nn++) = Eigen::RowVector4d(sphere.center().x(), sphere.center().y(), sphere.center().z(), sqrt(sphere.squared_radius()));
        }
    }

    if (debug_level > 0) { std::cout << "... main calculation... " << std::endl; }
    if (true) {
        if (debug_level > 0) { std::cout << "Gp: " << Gp.rows() << std::endl;}
        Pwns = contact_points(Gp, Pwns, (filter_contact_spheres_bbx)?&bbxl:nullptr, debug_level);
    } else {
        std::cout << "WARNING: ignoring positive spheres" << std::endl;
    }

    if (true) {
        if (debug_level > 0) { std::cout << "Gn: " << Gn.rows() << std::endl; }
        contact_points(Gn, Pwns, (filter_contact_spheres_bbx)?&bbxl:nullptr, debug_level);
    } else {
        std::cout << "WARNING: ignoring negative spheres" << std::endl;
    }
    return Pwns;
  }

} // namespace CGAL

#endif // CGAL_MAXIMAL_EMPTY_SPHERES_CONTACT_POINTS_FROM_SIGNED_DISTANCE_H
