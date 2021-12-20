// Copyright (c) 2018  Liangliang Nan. All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Liangliang Nan

#ifndef CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_ALPHA_SHAPE_MESH_H
#define CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_ALPHA_SHAPE_MESH_H

#include <CGAL/license/Polygonal_surface_reconstruction.h>

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Surface_mesh.h>


// warn the user if he/she is using an unsupported version of CGAL
#if CGAL_VERSION_NR < 1041100000
#error CGAL 4.11 or above is required (due to the breaking change in CGAL 4.11). Please update your code.
#endif

/*!
\file alpha_shape_mesh.h
*/

namespace CGAL {

        namespace internal {

                /// \cond SKIP_IN_MANUA

                template <typename Ht>
                class Alpha_shape;

                /// \endcond

                /// \cond SKIP_IN_MANUA

                /* A vertex class with an additional member representing its index */
                template < class Gt, class VB = CGAL::Triangulation_hierarchy_vertex_base_2<Gt> >
                class AS_vertex_base : public VB
                {
                public:
                        typedef VB                                                                Base;
                        typedef typename VB::Vertex_handle      Vertex_handle;
                        typedef typename VB::Face_handle        Face_handle;
                        typedef typename VB::Point              Point;

                        template < typename TDS2 >
                        struct Rebind_TDS {
                                typedef typename VB::template Rebind_TDS<TDS2>::Other        VB2;
                                typedef AS_vertex_base<Gt, VB2>                         Other;
                        };

                public:
                        AS_vertex_base() : Base(), index_(-1) {}
                        AS_vertex_base(const Point & p) : Base(p), index_(-1) {}
                        AS_vertex_base(const Point & p, Face_handle f) : Base(f, p), index_(-1) {}
                        AS_vertex_base(Face_handle f) : Base(f), index_(-1) {}

                        void set_index(int idx) { index_ = idx; }
                        int  index() const { return index_; }

                private:
                        int index_;
                };


                template <typename Ht>
                class Alpha_shape : public Alpha_shape_2<Ht>
                {
                public:
                        typedef Alpha_shape_2<Ht>                                                Parent_class;
                        typedef typename Ht::Point_2                                        Point2;
                        typedef typename Parent_class::Vertex_handle        Vertex_handle;

                public:
                        // constructs alpha shapes from the input points
                        template <typename InputIterator>
                        Alpha_shape(InputIterator first, InputIterator beyond);
                };

                /// \endcond


                /**
                *        An Alpha Shape Mesh approximates the point covered region by a mesh representation.
                */

                template <typename Kernel>
                class Alpha_shape_mesh
                {
                        typedef typename Kernel::FT                                                FT;
                        typedef typename Kernel::Point_2                                Point2;
                        typedef typename Kernel::Point_3                                Point3;
                        typedef typename Kernel::Plane_3                                Plane3;
                        typedef CGAL::Surface_mesh<Point3>                                Mesh3;
                        typedef typename Mesh3::Vertex_index                        Vertex_descriptor;

                        typedef CGAL::Alpha_shape_vertex_base_2<Kernel>                        Avb;
                        typedef AS_vertex_base<Avb>                                                                Av;
                        typedef CGAL::Triangulation_face_base_2<Kernel>                        Tf;
                        typedef CGAL::Alpha_shape_face_base_2<Kernel, Tf>                Af;
                        typedef CGAL::Triangulation_default_data_structure_2<Kernel, Av, Af> Tds;
                        typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>                Dt;
                        typedef CGAL::Triangulation_hierarchy_2<Dt>                                Ht;

                public:
                        /// Given a set of 3D points lying on 'plane', constructs alpha shapes from the
                        /// the projection of the points onto 'plane'
                        template <typename InputIterator>
                        Alpha_shape_mesh(InputIterator first, InputIterator beyond, const Plane3& plane);

                        ~Alpha_shape_mesh() { delete alpha_shape_; }

                        /// Extracts the 3D mesh representation of the alpha shapes
                        bool extract_mesh(FT alpha_value, Mesh3& mesh);

                private:
                        Alpha_shape<Ht>*                        alpha_shape_;
                        std::vector<const Point3*>  original_points_;
                };


                //////////////////////////////////////////////////////////////////////////

                // implementation


                template <typename Traits>
                template <typename InputIterator>
                Alpha_shape<Traits>::Alpha_shape(InputIterator first, InputIterator beyond) {
                        InputIterator it = first;
                        for (int id = 0; it != beyond; ++it, ++id) {
                                const Point2& p = *it;
                                Vertex_handle vh = Traits::insert(p);
                                if (vh->index() == -1)
                                        vh->set_index(id);
                                else {
                                        // p was not inserted (there might be a duplicated point)
                                }
                        }

                        if (Parent_class::dimension() == 2) {
                                // Computes the associated _interval_face_map
                                Parent_class::initialize_interval_face_map();

                                // Computes the associated _interval_edge_map
                                Parent_class::initialize_interval_edge_map();

                                // Computes the associated _interval_vertex_map
                                Parent_class::initialize_interval_vertex_map();

                                // Merges the two maps
                                Parent_class::initialize_alpha_spectrum();
                        }
                }


                template <typename Kernel>
                template <typename InputIterator>
                Alpha_shape_mesh<Kernel>::Alpha_shape_mesh(InputIterator first, InputIterator beyond, const Plane3& plane) {
                        original_points_.clear();

                        std::vector<Point2> pts;
                        for (InputIterator it = first; it != beyond; ++it) {
                                const Point3& p = *it;
                                const Point2& q = plane.to_2d(p);
                                pts.push_back(q);
                                original_points_.push_back(&p);
                        }
                        alpha_shape_ = new Alpha_shape<Ht>(pts.begin(), pts.end());
                }


                template <typename Kernel>
                bool Alpha_shape_mesh<Kernel>::extract_mesh(FT alpha_value, Mesh3& mesh) {
                        alpha_shape_->set_alpha(alpha_value);

                        typedef std::vector<std::size_t> Triangle;
                        std::vector<Triangle>        faces;

                        typedef Alpha_shape<Ht>        Alpha_shape;

                        typename Alpha_shape::Finite_faces_iterator fit = alpha_shape_->finite_faces_begin();
                        for (; fit != alpha_shape_->finite_faces_end(); ++fit) {
                                if (alpha_shape_->classify(fit) == Alpha_shape::INTERIOR) {
                                        Triangle tri;
                                        for (int i = 0; i < 3; ++i) {
                                                typename Alpha_shape::Vertex_handle vh = fit->vertex(i);
                                                int idx = vh->index();
                                                tri.push_back(idx);
                                        }
                                        faces.push_back(tri);
                                }
                        }

                        if (faces.empty())
                                return false;

                        mesh.clear();

                        std::vector<Vertex_descriptor> descriptors(original_points_.size());
                        for (std::size_t i = 0; i < original_points_.size(); ++i) {
                                const Point3* p = original_points_[i];
                                descriptors[i] = mesh.add_vertex(*p);
                        }

                        for (std::size_t i = 0; i < faces.size(); ++i) {
                                std::vector<Vertex_descriptor> face;
                                const Triangle& tri = faces[i];
                                for (std::size_t j = 0; j < tri.size(); ++j) {
                                        std::size_t idx = tri[j];
                                        face.push_back(descriptors[idx]);
                                }
                                mesh.add_face(face);;
                        }

                        return true;
                }

        } //namespace internal

} //namespace CGAL

#endif        // CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_ALPHA_SHAPE_MESH_H
