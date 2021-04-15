// Copyright (c) 2018  Liangliang Nan. All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Liangliang Nan

#ifndef CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_HYPOTHESIS_H
#define CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_HYPOTHESIS_H

#include <CGAL/license/Polygonal_surface_reconstruction.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/intersections.h>
#include <CGAL/assertions.h>
#include <CGAL/internal/parameters.h>
#include <CGAL/internal/point_set_with_planes.h>

#include <set>
#include <unordered_map>


namespace CGAL {

        namespace internal {

                /**
                *         Generates candidate faces by pairwise intersecting of the supporting planes of
                *   the planar segments.
                */

                template <typename Kernel>
                class Hypothesis
                {
                private:
                        typedef typename Kernel::FT                                FT;
                        typedef typename Kernel::Point_3                        Point;
                        typedef typename Kernel::Point_2                        Point2;
                        typedef typename Kernel::Vector_3                        Vector;
                        typedef typename Kernel::Line_3                                Line;
                        typedef typename Kernel::Segment_3                        Segment;
                        typedef typename Kernel::Plane_3                        Plane;
                        typedef internal::Planar_segment<Kernel>                Planar_segment;
                        typedef internal::Point_set_with_planes<Kernel>                Point_set_with_planes;

                        typedef CGAL::Surface_mesh<Point>                        Polygon_mesh;
                        typedef typename Polygon_mesh::Face_index                Face_descriptor;
                        typedef typename Polygon_mesh::Edge_index                Edge_descriptor;
                        typedef typename Polygon_mesh::Vertex_index                Vertex_descriptor;

                        typedef typename Polygon_mesh::Halfedge_index        Halfedge_descriptor;

                public:
                        Hypothesis();
                        ~Hypothesis();

                        void generate(Point_set_with_planes& point_set, Polygon_mesh& candidate_faces);

                        /// 'Intersection' represents a set of faces intersecting at a common edge.
                        /// \note The faces are represented by their halfedges.
                        struct Intersection : public std::vector<Halfedge_descriptor> {
                                const Point* s;
                                const Point* t;
                        };
                        typedef typename std::vector<Intersection>                Adjacency;
                        /// Extracts the adjacency of the pairwise intersection.
                        /// The extracted adjacency will be used to formulate the hard constraints
                        /// in the face selection stage.
                        Adjacency extract_adjacency(const Polygon_mesh& candidate_faces);

                private:
                        // Merges near co-planar segments
                        void refine_planes();

                        // Constructs a mesh representing the bounding box of the point set
                        void construct_bbox_mesh(Polygon_mesh& bbox_mesh);

                        // Construct a mesh from the segments bounded by the bounding box mesh
                        void construct_proxy_mesh(const Polygon_mesh& bbox_mesh, Polygon_mesh& candidate_faces);

                        // Pairwise intersection
                        void pairwise_intersection(Polygon_mesh& candidate_faces);

                        // Counts the number of points that are with the dist_threshold to its supporting plane
                        std::size_t number_of_points_on_plane(const Planar_segment* s, const Plane* plane, FT dist_threshold);

                        // Merges two planar segments;
                        void merge(Planar_segment* s1, Planar_segment* s2);

                        // Pre-computes all potential intersections of plane triplets
                        void compute_triplet_intersections();

                        // Queries the intersecting point for a plane triplet
                        const Point* query_intersection(const Plane* plane1, const Plane* plane2, const Plane* plane3);

                        bool halfedge_exists(Vertex_descriptor v1, Vertex_descriptor v2, const Polygon_mesh& mesh);

                        // Tests if 'face' insects 'plane'
                        bool do_intersect(const Polygon_mesh& mesh, Face_descriptor face, const Plane* plane);

                        // Cuts face using the cutting_plane and returns the new faces
                        std::vector<Face_descriptor> cut(Face_descriptor face, const Plane* cutting_plane, Polygon_mesh& mesh);

                        // Collects all faces in 'mesh' that intersect 'face'. Store in std::set() so easier to erase
                        std::set<Face_descriptor> collect_intersecting_faces(Face_descriptor face, const Polygon_mesh& mesh);

                        // Represents an intersecting point at an edge
                        struct EdgePos {
                                EdgePos(Edge_descriptor e, const Point* p) : edge(e), pos(p) {}
                                Edge_descriptor edge;
                                const Point*        pos;
                        };

                        // Computes the intersecting points of face and cutting_plane. The intersecting points are returned
                        // by 'existing_vts' (if the plane intersects the face at its vertices) and 'new_vts' (if the plane
                        // intersects the face at its edges).
                        void compute_intersections(const Polygon_mesh& mesh,
                                Face_descriptor face, const Plane* cutting_plane,
                                std::vector<Vertex_descriptor>& existing_vts,
                                std::vector<EdgePos>& new_vts
                        );

                        // This function will
                        // - split an edge denoted by 'ep'
                        // - assign the new edges the supporting faces
                        // - return the halfedge pointing to the new vertex
                        // Internally it uses Euler split_edge().
                        Halfedge_descriptor split_edge(Polygon_mesh& mesh, const EdgePos& ep, const Plane* cutting_plane);

                        // Clears cached intermediate results
                        void clear();

                private:
                        // The input point cloud with planes
                        Point_set_with_planes * point_set_;

                        // The intersection of the planes can be unreliable when the planes are near parallel.
                        // Here are the tricks we use in our implementation:
                        //   - We first test if an intersection exists for every pair of planes. We then collect
                        //     plane triplets such that every pair in the plane triplet intersect. This is achieved
                        //     by testing each plane against the known intersecting pairs.
                        //         - The 3D vertices of the final faces are obtained by computing the intersections of
                        //     the plane triplets. To cope with limited floating point precision, each vertex is
                        //     identified by the pointers of (in an increasing order) of the three planes from
                        //     which it is computed. By doing so, two vertices with almost identical positions can
                        //     be distinguished. This turned out to be quite robust in handling very close and near
                        //     parallel planes.

                        // The supporting planes of all planar segments and the bounding box faces
                        std::vector<const Plane*>        supporting_planes_;

                        // Precomputed intersecting points of all plane triplets
                        std::vector<const Point*>        intersecting_points_;

                        typedef typename std::unordered_map<const Plane*, const Point*>                                Plane_to_point_map;
                        typedef typename std::unordered_map<const Plane*, Plane_to_point_map>                Two_planes_to_point_map;
                        typedef typename std::unordered_map<const Plane*, Two_planes_to_point_map>        Planes_to_point_map;
                        Planes_to_point_map                        triplet_intersections_;
                };


                //////////////////////////////////////////////////////////////////////////

                // implementation


                template <typename Kernel>
                Hypothesis<Kernel>::Hypothesis()
                {
                }


                template <typename Kernel>
                Hypothesis<Kernel>::~Hypothesis()
                {
                        clear();
                }


                template <typename Kernel>
                void Hypothesis<Kernel>::clear() {
                        for (std::size_t i = 0; i < supporting_planes_.size(); ++i)
                                delete supporting_planes_[i];
                        supporting_planes_.clear();

                        for (std::size_t i = 0; i < intersecting_points_.size(); ++i)
                                delete intersecting_points_[i];
                        intersecting_points_.clear();

                        triplet_intersections_.clear();
                }


                template <typename Kernel>
                void Hypothesis<Kernel>::generate(Point_set_with_planes& point_set, Polygon_mesh& candidate_faces) {
                        point_set_ = &point_set;

                        refine_planes();

                        Polygon_mesh bbox_mesh;
                        construct_bbox_mesh(bbox_mesh);

                        construct_proxy_mesh(bbox_mesh, candidate_faces);

                        pairwise_intersection(candidate_faces);
                }


                /// \cond SKIP_IN_MANUAL

                template <typename Planar_segment>
                class SegmentSizeIncreasing
                {
                public:
                        SegmentSizeIncreasing() {}
                        bool operator()(const Planar_segment* s0, const Planar_segment* s1) const {
                                return s0->size() < s1->size();
                        }
                };

                template <typename Planar_segment>
                class SegmentSizeDecreasing
                {
                public:
                        SegmentSizeDecreasing() {}
                        bool operator()(const Planar_segment* s0, const Planar_segment* s1) const {
                                return s0->size() > s1->size();
                        }
                };

                template <typename FT, typename Vector>
                void normalize(Vector& v) {
                        FT s = std::sqrt(v.squared_length());
                        if (s > 1e-30)
                                s = FT(1) / s;
                        v *= s;
                }

                template <typename BBox>
                typename BBox::FT bbox_radius(const BBox& box) {
                        typedef typename BBox::FT FT;
                        FT dx = box.xmax() - box.xmin();
                        FT dy = box.ymax() - box.ymin();
                        FT dz = box.zmax() - box.zmin();
                        return FT(0.5) * std::sqrt(dx * dx + dy * dy + dz * dz);
                }

                // Computes the intersection of a plane triplet
                // Returns true if the intersection exists (p returns the point)
                template <typename Plane, typename Point>
                bool intersect_plane_triplet(const Plane* plane1, const Plane* plane2, const Plane* plane3, Point& p) {
                        if (plane1 == plane2 || plane1 == plane3 || plane2 == plane3)
                                return false;

                        CGAL::Object obj = CGAL::intersection(*plane1, *plane2, *plane3);

                        // pt is the intersection point of the 3 planes
                        if (const Point* pt = CGAL::object_cast<Point>(&obj)) {
                                p = *pt;
                                return true;
                        }
                        else {
                                // If reached here, the reason might be:
                                //   (1) two or more are parallel;
                                //   (2) they intersect at the same line
                                // We can simply ignore these cases
                                return false;
                        }
                }


                template <typename VT>
                void sort_increasing(VT& v1, VT& v2, VT& v3) {
                        VT vmin = 0;
                        if (v1 < v2 && v1 < v3)
                                vmin = v1;
                        else if (v2 < v1 && v2 < v3)
                                vmin = v2;
                        else
                                vmin = v3;

                        VT vmid = 0;
                        if ((v1 > v2 && v1 < v3) || (v1 < v2 && v1 > v3))
                                vmid = v1;
                        else if ((v2 > v1 && v2 < v3) || (v2 < v1 && v2 > v3))
                                vmid = v2;
                        else
                                vmid = v3;

                        VT vmax = 0;
                        if (v1 > v2 && v1 > v3)
                                vmax = v1;
                        else if (v2 > v1 && v2 > v3)
                                vmax = v2;
                        else
                                vmax = v3;

                        v1 = vmin;
                        v2 = vmid;
                        v3 = vmax;
                }
                /// \endcond


                template <typename Kernel>
                std::size_t Hypothesis<Kernel>::number_of_points_on_plane(const Planar_segment* s, const Plane* plane, FT dist_threshold) {
                        CGAL_assertion(const_cast<Planar_segment*>(s)->point_set() == point_set_);

                        std::size_t count = 0;
                        const typename Point_set_with_planes::Point_map& points = point_set_->point_map();
                        for (std::size_t i = 0; i < s->size(); ++i) {
                                std::size_t idx = s->at(i);
                                const Point& p = points[idx];

                                FT sdist = CGAL::squared_distance(*plane, p);
                                FT dist = std::sqrt(sdist);
                                if (dist < dist_threshold)
                                        ++count;
                        }
                        return count;
                }

                template <typename Kernel>
                void Hypothesis<Kernel>::merge(Planar_segment* s1, Planar_segment* s2) {
                        CGAL_assertion(const_cast<Planar_segment*>(s1)->point_set() == point_set_);
                        CGAL_assertion(const_cast<Planar_segment*>(s2)->point_set() == point_set_);
                        std::vector< Planar_segment* >& segments = point_set_->planar_segments();

                        std::vector<std::size_t> points_indices;
                        points_indices.insert(points_indices.end(), s1->begin(), s1->end());
                        points_indices.insert(points_indices.end(), s2->begin(), s2->end());

                        Planar_segment* s = new Planar_segment(point_set_);
                        s->insert(s->end(), points_indices.begin(), points_indices.end());
                        s->fit_supporting_plane();
                        segments.push_back(s);

                        typename std::vector< Planar_segment* >::iterator pos = std::find(segments.begin(), segments.end(), s1);
                        if (pos != segments.end()) {
                                Planar_segment* tmp = *pos;
                                const Plane* plane = tmp->supporting_plane();
                                segments.erase(pos);
                                delete tmp;
                                delete plane;
                        }
                        else
                                std::cerr << "Fatal error: should not reach here" << std::endl;

                        pos = std::find(segments.begin(), segments.end(), s2);
                        if (pos != segments.end()) {
                                Planar_segment* tmp = *pos;
                                const Plane* plane = tmp->supporting_plane();
                                segments.erase(pos);
                                delete tmp;
                                delete plane;
                        }
                        else
                                std::cerr << "Fatal error: should not reach here" << std::endl;
                }

                template <typename Kernel>
                void Hypothesis<Kernel>::refine_planes() {
                        std::vector< Planar_segment* >& segments = point_set_->planar_segments();
                        const typename Point_set_with_planes::Point_map& points = point_set_->point_map();

                        FT avg_max_dist = 0;
                        for (std::size_t i = 0; i < segments.size(); ++i) {
                                Planar_segment* s = segments[i];
                                const Plane* plane = s->fit_supporting_plane(); // user may provide invalid plane fitting (we always fit)

                                FT max_dist = -(std::numeric_limits<FT>::max)();
                                for (std::size_t j = 0; j < s->size(); ++j) {
                                        std::size_t idx = s->at(j);
                                        const Point& p = points[idx];
                                        FT sdist = CGAL::squared_distance(*plane, p);
                                        max_dist = (std::max)(max_dist, std::sqrt(sdist));
                                }

                                avg_max_dist += max_dist;
                        }
                        avg_max_dist /= segments.size();
                        avg_max_dist /= FT(2.0);

                        FT theta = static_cast<FT>(CGAL_PI * 10.0 / FT(180.0));        // in radian
                        bool merged = false;
                        do {
                                merged = false;
                                // Segments with less points have less confidences and thus should be merged first.
                                // So we sort the segments according to their sizes.
                                std::sort(segments.begin(), segments.end(), internal::SegmentSizeIncreasing<Planar_segment>());

                                for (std::size_t i = 0; i < segments.size(); ++i) {
                                        Planar_segment* s1 = segments[i];
                                        const Plane* plane1 = s1->supporting_plane();
                                        Vector n1 = plane1->orthogonal_vector();
                                        internal::normalize<FT, Vector>(n1);

                                        FT num_threshold = s1->size() / FT(5.0);
                                        for (std::size_t j = i + 1; j < segments.size(); ++j) {
                                                Planar_segment* s2 = segments[j];
                                                const Plane* plane2 = s2->supporting_plane();
                                                Vector n2 = plane2->orthogonal_vector();
                                                internal::normalize<FT, Vector>(n2);

                                                if (std::abs(n1 * n2) > std::cos(theta)) {
                                                        std::size_t set1on2 = number_of_points_on_plane(s1, plane2, avg_max_dist);
                                                        std::size_t set2on1 = number_of_points_on_plane(s2, plane1, avg_max_dist);
                                                        if (set1on2 > num_threshold || set2on1 > num_threshold) {
                                                                merge(s1, s2);
                                                                merged = true;
                                                                break;
                                                        }
                                                }
                                        }
                                        if (merged)
                                                break;
                                }
                        } while (merged);

                        std::sort(segments.begin(), segments.end(), internal::SegmentSizeDecreasing<Planar_segment>());

                        // Stores all the supporting planes
                        for (std::size_t i = 0; i < segments.size(); ++i) {
                                Planar_segment* s = segments[i];
                                const Plane* plane = s->supporting_plane();
                                supporting_planes_.push_back(plane);
                        }
                }


                template <typename Kernel>
                void Hypothesis<Kernel>::construct_bbox_mesh(Polygon_mesh& mesh) {
                        const typename Point_set_with_planes::Point_map& points = point_set_->point_map();

                        typedef typename Kernel::Iso_cuboid_3 BBox;
                        const BBox& box = CGAL::bounding_box(points.begin(), points.end());

                        FT dx = box.xmax() - box.xmin();
                        FT dy = box.ymax() - box.ymin();
                        FT dz = box.zmax() - box.zmin();
                        FT radius = FT(0.5) * std::sqrt(dx * dx + dy * dy + dz * dz);
                        FT offset = radius * FT(0.05);

                        // make the box larger to ensure all points are enclosed.
                        FT xmin = box.xmin() - offset, xmax = box.xmax() + offset;
                        FT ymin = box.ymin() - offset, ymax = box.ymax() + offset;
                        FT zmin = box.zmin() - offset, zmax = box.zmax() + offset;

                        mesh.clear();

                        Vertex_descriptor v0 = mesh.add_vertex(Point(xmin, ymin, zmin));  // 0
                        Vertex_descriptor v1 = mesh.add_vertex(Point(xmax, ymin, zmin));  // 1
                        Vertex_descriptor v2 = mesh.add_vertex(Point(xmax, ymin, zmax));  // 2
                        Vertex_descriptor v3 = mesh.add_vertex(Point(xmin, ymin, zmax));  // 3
                        Vertex_descriptor v4 = mesh.add_vertex(Point(xmax, ymax, zmax));  // 4
                        Vertex_descriptor v5 = mesh.add_vertex(Point(xmax, ymax, zmin));  // 5
                        Vertex_descriptor v6 = mesh.add_vertex(Point(xmin, ymax, zmin));  // 6
                        Vertex_descriptor v7 = mesh.add_vertex(Point(xmin, ymax, zmax));  // 7

                        mesh.add_face(v0, v1, v2, v3);
                        mesh.add_face(v1, v5, v4, v2);
                        mesh.add_face(v1, v0, v6, v5);
                        mesh.add_face(v4, v5, v6, v7);
                        mesh.add_face(v0, v3, v7, v6);
                        mesh.add_face(v2, v4, v7, v3);

                        // The supporting plane of each face
                        typename Polygon_mesh::template Property_map<Face_descriptor, const Plane*> face_supporting_planes =
                                mesh.template add_property_map<Face_descriptor, const Plane*>("f:supp_plane").first;

                        // The supporting planes of each edge
                        typename Polygon_mesh::template Property_map<Edge_descriptor, std::set<const Plane*> > edge_supporting_planes
                                = mesh.template add_property_map<Edge_descriptor, std::set<const Plane*> >("e:supp_plane").first;

                        // The supporting planes of each vertex
                        typename Polygon_mesh::template Property_map<Vertex_descriptor, std::set<const Plane*> > vertex_supporting_planes
                                = mesh.template add_property_map<Vertex_descriptor, std::set<const Plane*> >("v:supp_plane").first;

                        // Assigns the original plane for each face
                        const typename  Polygon_mesh::template Property_map<Vertex_descriptor, Point>& coords = mesh.points();
                        for(auto fd : mesh.faces()) {
                                Halfedge_descriptor h = mesh.halfedge(fd);
                                Vertex_descriptor va = mesh.target(h);        const Point& pa = coords[va]; h = mesh.next(h);
                                Vertex_descriptor vb = mesh.target(h);        const Point& pb = coords[vb]; h = mesh.next(h);
                                Vertex_descriptor vc = mesh.target(h);        const Point& pc = coords[vc];
                                const Plane* plane = new Plane(pa, pb, pc);
                                supporting_planes_.push_back(plane);
                                face_supporting_planes[fd] = plane;
                        }

                        // Assigns the original planes for each edge
                        for( auto ed : mesh.edges()) {
                                Halfedge_descriptor h1 = mesh.halfedge(ed);
                                Halfedge_descriptor h2 = mesh.opposite(h1);

                                Face_descriptor f1 = mesh.face(h1);
                                Face_descriptor f2 = mesh.face(h2);
                                CGAL_assertion(f1 != Polygon_mesh::null_face()); // the bbox mesh is closed
                                CGAL_assertion(f2 != Polygon_mesh::null_face()); // the bbox mesh is closed

                                const Plane* plane1 = face_supporting_planes[f1];
                                const Plane* plane2 = face_supporting_planes[f2];
                                CGAL_assertion(plane1 && plane2 && plane1 != plane2);

                                edge_supporting_planes[ed].insert(plane1);
                                edge_supporting_planes[ed].insert(plane2);
                                CGAL_assertion(edge_supporting_planes[ed].size() == 2);
                        }

                        // Assigns the original planes for each vertex
                        for(auto vd : mesh.vertices()) {
                                CGAL_assertion(vertex_supporting_planes[vd].size() == 0);
                                CGAL::Halfedge_around_target_circulator<Polygon_mesh> hbegin(vd, mesh), done(hbegin);
                                do {
                                        Halfedge_descriptor h = *hbegin;
                                        Face_descriptor f = mesh.face(h);
                                        const Plane* plane = face_supporting_planes[f];
                                        vertex_supporting_planes[vd].insert(plane);
                                        ++hbegin;
                                } while (hbegin != done);
                                CGAL_assertion(vertex_supporting_planes[vd].size() == 3);
                        }

                        std::sort(supporting_planes_.begin(), supporting_planes_.end());

                        CGAL_assertion(mesh.is_valid());
                }


                template <typename Kernel>
                void Hypothesis<Kernel>::construct_proxy_mesh(const Polygon_mesh& bbox_mesh, Polygon_mesh& candidate_faces) {

                        // Properties of the bbox_mesh

                        typename Polygon_mesh::template Property_map<Edge_descriptor, std::set<const Plane*> > bbox_edge_supporting_planes
                                = bbox_mesh.template property_map<Edge_descriptor, std::set<const Plane*> >("e:supp_plane").first;
                        typename Polygon_mesh::template Property_map<Vertex_descriptor, std::set<const Plane*> > bbox_vertex_supporting_planes
                                = bbox_mesh.template property_map<Vertex_descriptor, std::set<const Plane*> >("v:supp_plane").first;

                        // The properties of the proxy mesh
                        candidate_faces.clear();

                        // The supporting plane of each face
                        typename Polygon_mesh::template Property_map<Face_descriptor, const Plane*> face_supporting_planes
                                = candidate_faces.template add_property_map<Face_descriptor, const Plane*>("f:supp_plane").first;

                        // The supporting planar segment of each face
                        typename Polygon_mesh::template Property_map<Face_descriptor, Planar_segment*> face_supporting_segments
                                = candidate_faces.template add_property_map<Face_descriptor, Planar_segment*>("f:supp_segment").first;

                        // The supporting planes of each edge
                        typename Polygon_mesh::template Property_map<Edge_descriptor, std::set<const Plane*> > edge_supporting_planes
                                = candidate_faces.template add_property_map<Edge_descriptor, std::set<const Plane*> >("e:supp_plane").first;

                        // The supporting planes of each vertex
                        typename Polygon_mesh::template Property_map<Vertex_descriptor, std::set<const Plane*> > vertex_supporting_planes
                                = candidate_faces.template add_property_map<Vertex_descriptor, std::set<const Plane*> >("v:supp_plane").first;

                        const std::vector<Planar_segment*>& segments = point_set_->planar_segments();
                        const typename Polygon_mesh::template Property_map<Vertex_descriptor, Point>& coords = bbox_mesh.points();

                        for (std::size_t i = 0; i < segments.size(); ++i) {
                                Planar_segment* g = segments[i];
                                const Plane* cutting_plane = g->supporting_plane();

                                std::vector<Point> intersecting_points;
                                std::vector< std::set<const Plane*> > intersecting_points_source_planes;

                                for(auto ed : bbox_mesh.edges()) {
                                        Vertex_descriptor sd = bbox_mesh.vertex(ed, 0);
                                        Vertex_descriptor td = bbox_mesh.vertex(ed, 1);
                                        const Point& s = coords[sd];
                                        const Point& t = coords[td];

                                        CGAL::Oriented_side ss = cutting_plane->oriented_side(s);
                                        CGAL::Oriented_side st = cutting_plane->oriented_side(t);

                                        if ((ss == CGAL::ON_POSITIVE_SIDE && st == CGAL::ON_NEGATIVE_SIDE) || (ss == ON_NEGATIVE_SIDE && st == CGAL::ON_POSITIVE_SIDE)) {
                                                CGAL::Object obj = CGAL::intersection(*cutting_plane, Line(s, t));
                                                if (const Point* p = CGAL::object_cast<Point>(&obj)) {
                                                        intersecting_points.push_back(*p);
                                                        std::set<const Plane*> planes = bbox_edge_supporting_planes[ed];

                                                        planes.insert(cutting_plane);
                                                        CGAL_assertion(planes.size() == 3);
                                                        intersecting_points_source_planes.push_back(planes);
                                                }
                                                else
                                                        std::cerr << "Fatal error: should have intersection" << std::endl;
                                        }
                                        else {
                                                if (ss == CGAL::ON_ORIENTED_BOUNDARY && st != CGAL::ON_ORIENTED_BOUNDARY) {
                                                        intersecting_points.push_back(s);
                                                        const std::set<const Plane*>& planes = bbox_vertex_supporting_planes[sd];
                                                        CGAL_assertion(planes.size() == 3);
                                                        intersecting_points_source_planes.push_back(planes);
                                                }
                                                else if (st == CGAL::ON_ORIENTED_BOUNDARY && ss != CGAL::ON_ORIENTED_BOUNDARY) {
                                                        intersecting_points.push_back(t);
                                                        const std::set<const Plane*>& planes = bbox_vertex_supporting_planes[td];
                                                        CGAL_assertion(planes.size() == 3);
                                                        intersecting_points_source_planes.push_back(planes);
                                                }
                                                else {
                                                        // The intersection is the current edge, nothing to do
                                                }
                                        }
                                }

                                // Decides the orientation of the points
                                if (intersecting_points.size() >= 3) {
                                        std::list<Point> pts;
                                        for (std::size_t i = 0; i < intersecting_points.size(); ++i) {
                                                const Point& p = intersecting_points[i];
                                                const Point2& q = cutting_plane->to_2d(p);
                                                pts.push_back(Point(q.x(), q.y(), FT(i))); // the z component stores the point index
                                        }

                                        typedef CGAL::Projection_traits_xy_3<Kernel>  Projection;
                                        std::list<Point> hull;
                                        CGAL::convex_hull_2(pts.begin(), pts.end(), std::back_inserter(hull), Projection());

                                        std::vector<Point> ch;
                                        std::vector< std::set<const Plane*> > ch_source_planes;
                                        for (typename std::list<Point>::iterator it = hull.begin(); it != hull.end(); ++it) {
                                                std::size_t idx = std::size_t(it->z());
                                                ch.push_back(intersecting_points[idx]);
                                                ch_source_planes.push_back(intersecting_points_source_planes[idx]);
                                        }

                                        if (ch.size() >= 3) {
                                                std::vector<Vertex_descriptor> descriptors;
                                                for (std::size_t j = 0; j < ch.size(); ++j) {
                                                        Vertex_descriptor vd = candidate_faces.add_vertex(ch[j]);
                                                        descriptors.push_back(vd);
                                                        vertex_supporting_planes[vd] = ch_source_planes[j];
                                                        CGAL_assertion(vertex_supporting_planes[vd].size() == 3);
                                                }

                                                Face_descriptor fd = candidate_faces.add_face(descriptors);
                                                face_supporting_segments[fd] = g;
                                                face_supporting_planes[fd] = cutting_plane;

                                                // Assigns each edge the supporting planes
                                                CGAL::Halfedge_around_face_circulator<Polygon_mesh> hbegin(candidate_faces.halfedge(fd), candidate_faces), done(hbegin);
                                                do {
                                                        Halfedge_descriptor hd = *hbegin;
                                                        Edge_descriptor ed = candidate_faces.edge(hd);

                                                        Vertex_descriptor s_vd = candidate_faces.source(hd);
                                                        Vertex_descriptor t_vd = candidate_faces.target(hd);
                                                        const std::set<const Plane*>& s_planes = vertex_supporting_planes[s_vd];
                                                        const std::set<const Plane*>& t_planes = vertex_supporting_planes[t_vd];
                                                        std::set<const Plane*> common_planes;
                                                        std::set_intersection(s_planes.begin(), s_planes.end(), t_planes.begin(), t_planes.end(), std::inserter(common_planes, common_planes.begin()));
                                                        if (common_planes.size() == 2) {
                                                                CGAL_assertion(edge_supporting_planes[ed].size() == 0);
                                                                edge_supporting_planes[ed] = common_planes;
                                                                CGAL_assertion(edge_supporting_planes[ed].size() == 2);
                                                        }
                                                        else // If reached here, there must be topological errors.
                                                                std::cerr << "topological error" << std::endl;

                                                        ++hbegin;
                                                } while (hbegin != done);
                                        }
                                }
                        }

                        CGAL_assertion(candidate_faces.is_valid());
                }


                template <typename Kernel>
                void Hypothesis<Kernel>::compute_triplet_intersections() {
                        triplet_intersections_.clear();
                        if (supporting_planes_.size() < 4) // no closed surface will be constructed from less than 4 planes
                                return;

                        for (std::size_t i = 0; i < supporting_planes_.size(); ++i) {
                                const Plane* plane1 = supporting_planes_[i];
                                for (std::size_t j = i + 1; j < supporting_planes_.size(); ++j) {
                                        const Plane* plane2 = supporting_planes_[j];
                                        for (std::size_t k = j + 1; k < supporting_planes_.size(); ++k) {
                                                const Plane* plane3 = supporting_planes_[k];
                                                CGAL_assertion(plane1 < plane2 && plane2 < plane3);
                                                Point p;
                                                if (internal::intersect_plane_triplet<Plane, Point>(plane1, plane2, plane3, p)) {
                                                        // Stores the intersection for future query
                                                        Point* new_point = new Point(p);
                                                        triplet_intersections_[plane1][plane2][plane3] = new_point;
                                                        intersecting_points_.push_back(new_point);
                                                }
                                        }
                                }
                        }
                }


                template <typename Kernel>
                const typename Hypothesis<Kernel>::Point*
                        Hypothesis<Kernel>::query_intersection(const Plane* min_plane, const Plane* mid_plane, const Plane* max_plane) {
                        CGAL_assertion(min_plane < mid_plane);
                        CGAL_assertion(mid_plane < max_plane);

                        if (triplet_intersections_.find(min_plane) == triplet_intersections_.end())
                                return nullptr;

                        Two_planes_to_point_map& map2 = triplet_intersections_[min_plane];
                        if (map2.find(mid_plane) == map2.end())
                                return nullptr;

                        Plane_to_point_map& map1 = map2[mid_plane];
                        if (map1.find(max_plane) == map1.end())
                                return nullptr;

                        return map1[max_plane];
                }


                template <typename Kernel>
                typename Hypothesis<Kernel>::Halfedge_descriptor
                        Hypothesis<Kernel>::split_edge(Polygon_mesh& mesh, const EdgePos& ep, const Plane* cutting_plane) {
                        // The supporting planes of each edge
                        typename Polygon_mesh::template Property_map<Edge_descriptor, std::set<const Plane*> > edge_supporting_planes =
                                mesh.template property_map<Edge_descriptor, std::set<const Plane*> >("e:supp_plane").first;

                        // The supporting planes of each vertex
                        typename Polygon_mesh::template Property_map<Vertex_descriptor, std::set<const Plane*> > vertex_supporting_planes
                                = mesh.template property_map<Vertex_descriptor, std::set<const Plane*> >("v:supp_plane").first;

                        // We cannot use const reference, because it will become invalid after splitting
                        std::set<const Plane*> sfs = edge_supporting_planes[ep.edge];
                        CGAL_assertion(sfs.size() == 2);

                        Halfedge_descriptor h = Euler::split_edge(mesh.halfedge(ep.edge), mesh);
                        if (h == Polygon_mesh::null_halfedge()) // failed splitting edge
                                return h;

                        Vertex_descriptor v = mesh.target(h);
                        if (v == Polygon_mesh::null_vertex())        // failed splitting edge
                                return Polygon_mesh::null_halfedge();

                        typename Polygon_mesh::template Property_map<Vertex_descriptor, Point>& coords = mesh.points();
                        coords[v] = *ep.pos;

                        Edge_descriptor e1 = mesh.edge(h);
                        edge_supporting_planes[e1] = sfs;
                        Edge_descriptor e2 = mesh.edge(mesh.next(h));
                        edge_supporting_planes[e2] = sfs;

                        vertex_supporting_planes[v] = sfs;
                        vertex_supporting_planes[v].insert(cutting_plane);
                        CGAL_assertion(vertex_supporting_planes[v].size() == 3);

                        return h;
                }


                // Cuts f using the cutter and returns the new faces
                template <typename Kernel>
                std::vector<typename Hypothesis<Kernel>::Face_descriptor>
                        Hypothesis<Kernel>::cut(Face_descriptor face, const Plane* cutting_plane, Polygon_mesh& mesh) {
                        std::vector<Face_descriptor> new_faces;

                        // The supporting plane of each face
                        typename Polygon_mesh::template Property_map<Face_descriptor, const Plane*> face_supporting_planes =
                                mesh.template property_map<Face_descriptor, const Plane*>("f:supp_plane").first;
                        const Plane* supporting_plane = face_supporting_planes[face];

                        if (supporting_plane == cutting_plane)
                                return new_faces;

                        // The supporting planar segment of each face
                        typename Polygon_mesh::template Property_map<Face_descriptor, Planar_segment*> face_supporting_segments =
                                mesh.template property_map<Face_descriptor, Planar_segment*>("f:supp_segment").first;

                        // The supporting planes of each edge
                        typename Polygon_mesh::template Property_map<Edge_descriptor, std::set<const Plane*> > edge_supporting_planes =
                                mesh.template property_map<Edge_descriptor, std::set<const Plane*> >("e:supp_plane").first;

                        Planar_segment* supporting_segment = face_supporting_segments[face];

                        std::vector<Vertex_descriptor> existing_vts;
                        std::vector<EdgePos> new_vts;
                        compute_intersections(mesh, face, cutting_plane, existing_vts, new_vts);

                        // We need to check here because new faces are emerging
                        if (existing_vts.size() + new_vts.size() != 2)
                                return new_faces;

                        else if (existing_vts.size() == 2) {
                                // Tests if the two intersecting points are both very close to an existing vertex.
                                // Since we allow snapping, we test if the two intersecting points are the same.
                                if (existing_vts[0] == existing_vts[1])
                                        return new_faces;

                                // Tests if an edge already exists, i.e., the plane cuts at this edge
                                if (halfedge_exists(existing_vts[0], existing_vts[1], mesh))
                                        return new_faces;
                        }

                        Halfedge_descriptor h0 = Polygon_mesh::null_halfedge();
                        Halfedge_descriptor h1 = Polygon_mesh::null_halfedge();

                        if (existing_vts.size() == 2) { // cutting_plane cuts the face at two existing vertices (not an edge)
                                h0 = mesh.halfedge(existing_vts[0]);
                                h1 = mesh.halfedge(existing_vts[1]);
                        }
                        else if (existing_vts.size() == 1) {
                                h0 = mesh.halfedge(existing_vts[0]);
                                h1 = split_edge(mesh, new_vts[0], cutting_plane);
                        }
                        else if (new_vts.size() == 2) {
                                h0 = split_edge(mesh, new_vts[0], cutting_plane);
                                h1 = split_edge(mesh, new_vts[1], cutting_plane);
                        }
                        CGAL_assertion(h0 != Polygon_mesh::null_halfedge());
                        CGAL_assertion(h1 != Polygon_mesh::null_halfedge());

                        // To split the face, `h0` and `h1` must be incident to the same face
                        if (mesh.face(h0) != face) {
                                Halfedge_descriptor end = h0;
                                do {
                                        h0 = mesh.opposite(mesh.next(h0));        // h0 = h0->next()->opposite();
                                        if (mesh.face(h0) == face)
                                                break;
                                } while (h0 != end);
                        }
                        CGAL_assertion(mesh.face(h0) == face);

                        if (mesh.face(h1) != face) {
                                Halfedge_descriptor end = h1;
                                do {
                                        h1 = mesh.opposite(mesh.next(h1)); // h1 = h1->next()->opposite();
                                        if (mesh.face(h1) == face)
                                                break;
                                } while (h1 != end);
                        }
                        CGAL_assertion(mesh.face(h1) == face);

                        Halfedge_descriptor h = Euler::split_face(h0, h1, mesh);
                        if (h == Polygon_mesh::null_halfedge() || mesh.face(h) == Polygon_mesh::null_face()) {
                                std::cerr << "Fatal error. could not split face" << std::endl;
                                return new_faces;
                        }

                        Edge_descriptor e = mesh.edge(h);
                        edge_supporting_planes[e].insert(supporting_plane);
                        edge_supporting_planes[e].insert(cutting_plane);
                        CGAL_assertion(edge_supporting_planes[e].size() == 2);

                        // Now the two faces
                        Face_descriptor f1 = mesh.face(h);
                        face_supporting_segments[f1] = supporting_segment;
                        face_supporting_planes[f1] = supporting_plane;
                        new_faces.push_back(f1);

                        Face_descriptor f2 = mesh.face(mesh.opposite(h));
                        face_supporting_segments[f2] = supporting_segment;
                        face_supporting_planes[f2] = supporting_plane;
                        new_faces.push_back(f2);

                        return new_faces;
                }


                template <typename Kernel>
                void Hypothesis<Kernel>::compute_intersections(const Polygon_mesh& mesh,
                        Face_descriptor face, const Plane* cutting_plane,
                        std::vector<Vertex_descriptor>& existing_vts,
                        std::vector<EdgePos>& new_vts)
                {
                        existing_vts.clear();
                        new_vts.clear();

                        // The supporting plane of each face
                        typename Polygon_mesh::template Property_map<Face_descriptor, const Plane*> face_supporting_planes =
                                mesh.template property_map<Face_descriptor, const Plane*>("f:supp_plane").first;
                        const Plane* supporting_plane = face_supporting_planes[face];
                        if (supporting_plane == cutting_plane)
                                return;

                        typename Polygon_mesh::template Property_map<Edge_descriptor, std::set<const Plane*> > edge_supporting_planes
                                = mesh.template property_map<Edge_descriptor, std::set<const Plane*> >("e:supp_plane").first;

                        const typename Polygon_mesh::template Property_map<Vertex_descriptor, Point>& coords = mesh.points();

                        Halfedge_descriptor cur = mesh.halfedge(face);
                        Halfedge_descriptor end = cur;
                        do {
                                Edge_descriptor ed = mesh.edge(cur);
                                const std::set<const Plane*>& supporting_planes = edge_supporting_planes[ed];
                                if (supporting_planes.find(cutting_plane) != supporting_planes.end()) // the edge lies on the cutting plane
                                        return;

                                Vertex_descriptor s_vd = mesh.source(cur);
                                Vertex_descriptor t_vd = mesh.target(cur);
                                const Point& s = coords[s_vd];
                                const Point& t = coords[t_vd];

                                CGAL::Oriented_side s_side = cutting_plane->oriented_side(s);
                                CGAL::Oriented_side t_side = cutting_plane->oriented_side(t);

                                if (t_side == CGAL::ON_ORIENTED_BOUNDARY) {
                                        if (s_side == CGAL::ON_ORIENTED_BOUNDARY)  // the edge lies on the cutting plane
                                                return;
                                        else
                                                existing_vts.push_back(t_vd);
                                }

                                else if (
                                        (s_side == CGAL::ON_POSITIVE_SIDE && t_side == CGAL::ON_NEGATIVE_SIDE) ||
                                        (s_side == CGAL::ON_NEGATIVE_SIDE && t_side == CGAL::ON_POSITIVE_SIDE)) // intersects at the interior of the edge
                                {
                                        FT s_sdist = CGAL::squared_distance(*cutting_plane, s);
                                        FT t_sdist = CGAL::squared_distance(*cutting_plane, t);

                                        if (s_sdist <= CGAL::snap_squared_distance_threshold<FT>())                // plane cuts at vertex 's'
                                                existing_vts.push_back(s_vd);
                                        else if (t_sdist <= CGAL::snap_squared_distance_threshold<FT>()) // plane cuts at vertex 't'
                                                existing_vts.push_back(t_vd);
                                        else {
                                                const Plane* plane1 = *(supporting_planes.begin());
                                                const Plane* plane2 = *(supporting_planes.rbegin());
                                                const Plane* plane3 = const_cast<const Plane*>(cutting_plane);

                                                if (plane3 != plane1 && plane3 != plane2) {
                                                        internal::sort_increasing(plane1, plane2, plane3);
                                                        const Point* p = query_intersection(plane1, plane2, plane3);
                                                        if (p) {
                                                                if (CGAL::squared_distance(*p, s) <= CGAL::snap_squared_distance_threshold<FT>())                // snap to 's'
                                                                        existing_vts.push_back(s_vd);
                                                                else if (CGAL::squared_distance(*p, t) <= CGAL::snap_squared_distance_threshold<FT>())        // snap to 't'
                                                                        existing_vts.push_back(t_vd);
                                                                else
                                                                        new_vts.push_back(EdgePos(ed, p));
                                                        }
                                                        else
                                                                std::cerr << "Fatal error: should have intersection" << std::endl;
                                                }
                                                else
                                                        std::cerr << "Fatal error: should not have duplicated planes." << std::endl;
                                        }
                                }

                                else {
                                        // Nothing needs to do here, we will test the next edge
                                }

                                cur = mesh.next(cur);
                        } while (cur != end);
                }


                template <typename Kernel>
                bool Hypothesis<Kernel>::halfedge_exists(Vertex_descriptor v1, Vertex_descriptor v2, const Polygon_mesh& mesh) {
                        Halfedge_descriptor h = mesh.halfedge(v1);
                        Halfedge_descriptor end = h;
                        do {
                                Halfedge_descriptor opp = mesh.opposite(h);
                                if (mesh.target(opp) == v2)
                                        return true;
                                h = mesh.prev(opp);
                        } while (h != end);
                        return false;
                }


                // Tests if face 'f' insects plane 'plane'
                template <typename Kernel>
                bool Hypothesis<Kernel>::do_intersect(const Polygon_mesh& mesh, Face_descriptor f, const Plane* plane) {
                        std::vector<Vertex_descriptor> existing_vts;
                        std::vector<EdgePos> new_vts;
                        compute_intersections(mesh, f, plane, existing_vts, new_vts);

                        if (existing_vts.size() == 2) {
                                if (!halfedge_exists(existing_vts[0], existing_vts[1], mesh))
                                        return true;
                        }
                        else if (existing_vts.size() + new_vts.size() == 2)
                                return true;

                        return false;
                }


                // Collects all faces in 'mesh' that intersect 'face', stored in std::set() so easier to erase
                template <typename Kernel>
                std::set<typename Hypothesis<Kernel>::Face_descriptor> Hypothesis<Kernel>::collect_intersecting_faces(Face_descriptor face, const Polygon_mesh& mesh)
                {
                        // The supporting plane of each face
                        typename Polygon_mesh::template Property_map<Face_descriptor, const Plane*> face_supporting_planes =
                                mesh.template property_map<Face_descriptor, const Plane*>("f:supp_plane").first;

                        // The supporting planar segment of each face
                        typename Polygon_mesh::template Property_map<Face_descriptor, Planar_segment*> face_supporting_segments =
                                mesh.template property_map<Face_descriptor, Planar_segment*>("f:supp_segment").first;

                        std::set<Face_descriptor> intersecting_faces;
                        for(auto f : mesh.faces()) {
                                if (f == face ||
                                        face_supporting_segments[f] == face_supporting_segments[face] ||
                                        face_supporting_planes[f] == face_supporting_planes[face])
                                {
                                        continue;
                                }

                                const Plane* plane = face_supporting_planes[face];
                                CGAL_assertion(plane != nullptr);
                                if (do_intersect(mesh, f, plane))
                                        intersecting_faces.insert(f);
                        }
                        return intersecting_faces;
                }


                template <typename Kernel>
                void Hypothesis<Kernel>::pairwise_intersection(Polygon_mesh& candidate_faces)
                {
                        // Pre-computes all potential intersection of plane triplets
                        compute_triplet_intersections();

                        // Since we are going to split faces, we can not use the reference.
                        // const Polygon_mesh::Face_range& all_faces = mesh.faces();
                        // So we make a local copy
                        std::vector<Face_descriptor> all_faces(candidate_faces.faces().begin(), candidate_faces.faces().end());

                        // The supporting plane of each face
                        typename Polygon_mesh::template Property_map<Face_descriptor, const Plane*> face_supporting_planes
                                = candidate_faces.template property_map<Face_descriptor, const Plane*>("f:supp_plane").first;

                        for (std::size_t i = 0; i < all_faces.size(); ++i) {
                                Face_descriptor face = all_faces[i];
                                const Plane* face_plane = face_supporting_planes[face];

                                std::set<Face_descriptor> intersecting_faces = collect_intersecting_faces(face, candidate_faces);
                                if (intersecting_faces.empty())
                                        continue;

                                std::vector<Face_descriptor> cutting_faces(intersecting_faces.begin(), intersecting_faces.end());

                                // 1. 'face' will be cut by all the intersecting faces
                                //    \note After each cut, the original face doesn't exist any more and it is replaced by multiple pieces.
                                //          then each piece will be cut by another face.
                                std::vector<Face_descriptor> faces_to_be_cut;
                                faces_to_be_cut.push_back(face);
                                while (!intersecting_faces.empty()) {
                                        Face_descriptor cutting_face = *(intersecting_faces.begin());
                                        const Plane* cutting_plane = face_supporting_planes[cutting_face];

                                        std::set<Face_descriptor> new_faces;                // stores the new faces
                                        std::set<Face_descriptor> remained_faces;        // faces that will be cut later
                                        for (std::size_t j = 0; j < faces_to_be_cut.size(); ++j) {
                                                Face_descriptor current_face = faces_to_be_cut[j];
                                                std::vector<Face_descriptor> tmp = cut(current_face, cutting_plane, candidate_faces);
                                                new_faces.insert(tmp.begin(), tmp.end());
                                                if (tmp.empty()) {  // no actual cut occurred. The face will be cut later
                                                        remained_faces.insert(current_face);
                                                }
                                        }

                                        // The new faces might be cut by other faces
                                        faces_to_be_cut = std::vector<Face_descriptor>(new_faces.begin(), new_faces.end());

                                        // Don't forget the remained faces
                                        faces_to_be_cut.insert(faces_to_be_cut.end(), remained_faces.begin(), remained_faces.end());

                                        // The job of cutting_face is done, remove it
                                        intersecting_faces.erase(cutting_face);
                                }

                                // 2. All the cutting_faces will be cut by f.
                                for (std::size_t j = 0; j < cutting_faces.size(); ++j) {
                                        cut(cutting_faces[j], face_plane, candidate_faces);
                                }
                        }

                        CGAL_assertion(candidate_faces.is_valid());
                }


                template <class Kernel>
                typename Hypothesis<Kernel>::Adjacency Hypothesis<Kernel>::extract_adjacency(const Polygon_mesh& candidate_faces)
                {
                        typename Polygon_mesh::template Property_map<Vertex_descriptor, std::set<const Plane*> > vertex_supporting_planes
                                = candidate_faces.template property_map<Vertex_descriptor, std::set<const Plane*> >("v:supp_plane").first;

                        // An edge is denoted by its two end points
                        typedef typename std::unordered_map<const Point*, std::set<Halfedge_descriptor> >        Edge_map;
                        typedef typename std::unordered_map<const Point*, Edge_map >                                                Face_pool;
                        Face_pool face_pool;

                        for(auto h : candidate_faces.halfedges()) {
                                Face_descriptor f = candidate_faces.face(h);
                                if (f == Polygon_mesh::null_face())
                                        continue;

                                Vertex_descriptor sd = candidate_faces.source(h);
                                Vertex_descriptor td = candidate_faces.target(h);
                                const std::set<const Plane*>& set_s = vertex_supporting_planes[sd];
                                const std::set<const Plane*>& set_t = vertex_supporting_planes[td];
                                CGAL_assertion(set_s.size() == 3);
                                CGAL_assertion(set_t.size() == 3);

                                std::vector<const Plane*> s_planes(set_s.begin(), set_s.end());
                                CGAL_assertion(s_planes[0] < s_planes[1]);
                                CGAL_assertion(s_planes[1] < s_planes[2]);
                                const Point* s = triplet_intersections_[s_planes[0]][s_planes[1]][s_planes[2]];

                                std::vector<const Plane*> t_planes(set_t.begin(), set_t.end());
                                CGAL_assertion(t_planes[0] < t_planes[1]);
                                CGAL_assertion(t_planes[1] < t_planes[2]);
                                const Point* t = triplet_intersections_[t_planes[0]][t_planes[1]][t_planes[2]];

                                if (s > t)
                                        std::swap(s, t);
                                face_pool[s][t].insert(candidate_faces.halfedge(f));
                        }

                        Adjacency fans;
                        typename Face_pool::const_iterator it = face_pool.begin();
                        for (; it != face_pool.end(); ++it) {
                                const Point* s = it->first;
                                const Edge_map& tmp = it->second;
                                typename Edge_map::const_iterator cur = tmp.begin();
                                for (; cur != tmp.end(); ++cur) {
                                        const Point* t = cur->first;
                                        const std::set<Halfedge_descriptor>& faces = cur->second;
                                        Intersection fan;
                                        fan.s = s;
                                        fan.t = t;
                                        fan.insert(fan.end(), faces.begin(), faces.end());
                                        fans.push_back(fan);
                                }
                        }

                        return fans;
                }

        } //namespace internal

} //namespace CGAL


#endif // CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_HYPOTHESIS_H
