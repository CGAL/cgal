// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



enum class CDT2_Filtering {
    ODD_EVEN = 0, // in domain := odd nesting levels, see mark_domain_in_triangulation()
    NOT_OUT // in domain := everything that is not in the connected component incident to the infinite faces
};

auto triangulate_facet_with_CDT2(FacetSPtr facet,
                                 CDT2_Filtering filtering_policy,
                                 std::map<Point3, std::size_t>& pids,
                                 std::vector<Point3>& points,
                                 std::vector<std::vector<std::size_t> >& triangles)
{
    if (facet->edges().size() < 3) {
        std::cerr << "Warning: face with < 3 edges" << std::endl;
        return;
    } else {
        // @todo factorize CDT usages, but mind the tags
        using Itag = CGAL::Exact_intersections_tag; // since we do shift + epsilon, we could have intersections within a face
        using PK = CGAL::Projection_traits_3<CGAL::K>;
        using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>;
        using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
        using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
        using PTDS = CGAL::Triangulation_data_structure_2<PVb,PFb>;
        using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
        using PCDT_VH = PCDT::Vertex_handle;
        using PCDT_FH = PCDT::Face_handle;

        Vector3SPtr n = KernelFactory::createVector3(facet->plane());
        CGAL_assertion(*n != CGAL::NULL_VECTOR);

        PK projection_traits(*n);
        PCDT pcdt(projection_traits);

        std::map<VertexSPtr, PCDT_VH> face_vhs; // might have multiple vertices at the same position

        std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
        while (it_v != facet->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            auto res = face_vhs.emplace(vertex, PCDT_VH());
            if(res.second) // first time seeing this point
            {
                PCDT_VH vh = pcdt.insert(*(vertex->getPoint()));
                res.first->second = vh;
            }
        }

        auto ne = 0;
        std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
        while (it_e != facet->edges().end()) {
            EdgeSPtr edge = *it_e++;
            VertexSPtr v0 = edge->src(facet);
            VertexSPtr v1 = edge->dst(facet);
            std::cout << "CDT2 constraint: " << *(v0->getPoint()) << " || " << *(v1->getPoint()) << std::endl;

            // degenerate edges can happen e.g. at t=0 for vertices with deg > 3
            if(*(v0->getPoint()) == *(v1->getPoint())) {
                std::cerr << "Warning: encountered degenerate edge" << std::endl;
            } else {
                PCDT_VH vh0 = face_vhs.at(v0);
                PCDT_VH vh1 = face_vhs.at(v1);

                try {
                    pcdt.insert_constraint(vh0, vh1);
                } catch(const typename PCDT::Intersection_of_constraints_exception&) {
                    std::cerr << "Error: Intersection of constraint w/ " << vh0->point() << " " << vh1->point() << std::endl;
                    DEBUG_VAR(facet->toString());
                    CGAL_assertion_msg(false, "Intersections are allowed, throw shouldn't happen");
                    std::exit(1);
                }
                ++ne;
            }
        }

        if(ne < 3) { // degenerate face
            std::cerr << "Warning: skipping degenerate face (ne < 3)" << std::endl;
            return;
        }

        // @fixme there can be self-intersections, so doing the domain thing is wrong?
#ifndef CGAL_SS3_FILTER_CDT2_FACES_WITH_WINDING_NUMBER
        std::unordered_map<PCDT_FH, bool> in_domain_map;
        boost::associative_property_map<std::unordered_map<PCDT_FH, bool> > in_domain(in_domain_map);

        if (filtering_policy == CDT2_Filtering::ODD_EVEN) {
            CGAL::mark_domain_in_triangulation(pcdt, in_domain);
        } else {
            for (PCDT_FH fh : pcdt.all_face_handles()) {
                put(in_domain, fh, true);
            }

            PCDT_FH seed_fh = pcdt.infinite_vertex()->face();
            std::stack<PCDT_FH> to_explore;
            to_explore.push(seed_fh);
            while (!to_explore.empty()) {
                PCDT_FH fh = to_explore.top();
                to_explore.pop();
                if (!get(in_domain, fh)) {
                    continue;
                } else {
                    put(in_domain, fh, false);
                }
                for (std::size_t j=0; j<3; ++j) {
                    if (!fh->is_constrained(j)) {
                        to_explore.push(fh->neighbor(j));
                    }
                }
            }
        }
#endif // CGAL_SS3_FILTER_CDT2_FACES_WITH_WINDING_NUMBER

        for(auto fh : pcdt.finite_face_handles()) {
        std::cout << "Face handle "
                  << fh->vertex(0)->point() << " "
                  << fh->vertex(1)->point() << " "
                  << fh->vertex(2)->point()
                  << " in domain? " << get(in_domain, fh) << std::endl;
#ifndef CGAL_SS3_FILTER_CDT2_FACES_WITH_WINDING_NUMBER
            if(!get(in_domain, fh)) {
                continue;
            }
#endif // CGAL_SS3_FILTER_CDT2_FACES_WITH_WINDING_NUMBER

            // purge degenerate faces as we will things conformal with autoref anyway
            if (CGAL::collinear(fh->vertex(0)->point(),
                                fh->vertex(1)->point(),
                                fh->vertex(2)->point())) {
                std::cerr << "Warning: degenerate face in autoref event handling" << std::endl;
                continue;
            }

            auto get_pid = [&pids, &points] (PCDT_VH vh) -> std::size_t {
                auto res = pids.emplace(vh->point(), points.size());
                if(res.second) { // first time seeing the vertex handle
                    points.push_back(vh->point());
                }
                return res.first->second;
            };

            // trying something to avoid inverted subparts of a face from surviving an event:
            // compute the winding number with angles (obviously terrible both in complexity AND robustness!)
#ifdef CGAL_SS3_FILTER_CDT2_FACES_WITH_WINDING_NUMBER

            // this doesn't work because we can't distinguish between doubly inverted faces pointing
            // in the wrong direction, and non-inverted faces
#               error

            Point3 centroid = CGAL::centroid(fh->vertex(0)->point(),
                                              fh->vertex(1)->point(),
                                              fh->vertex(2)->point());

            CGAL::FT cumulative_angle = 0;
            CGAL::internal::Evaluate<CGAL::FT> evaluate;

            std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
            while (it_e != facet->edges().end()) {
                EdgeSPtr edge = *it_e++;
                VertexSPtr v_src = edge->src(facet);
                VertexSPtr v_dst = edge->dst(facet);
                CGAL_assertion(*(v_src->getPoint()) != *(v_dst->getPoint())); // @todo handle this, if it can happen

                Vector3 ln = CGAL::cross_product(Vector3(centroid, *(v_src->getPoint())),
                                                  Vector3(centroid, *(v_dst->getPoint())));
                CGAL_assertion(ln != CGAL::NULL_VECTOR);
                const CGAL::FT s = (CGAL::scalar_product(*n, ln) > 0) ? 1 : -1;

                cumulative_angle += s * CGAL::approximate_angle(*(v_src->getPoint()), centroid, *(v_dst->getPoint()));
                evaluate(cumulative_angle);

                std::cout << centroid << " || " << *(v_src->getPoint()) << " || " << *(v_dst->getPoint()) << std::endl;
                std::cout << "sign: " << s << std::endl;
                std::cout << "angle: " << CGAL::approximate_angle(*(v_src->getPoint()), centroid, *(v_dst->getPoint())) << std::endl;
            }

            std::cout << "cumulative angle = " << cumulative_angle
                      << " (factor = " << cumulative_angle / 360 << ")" << std::endl;
            if (cumulative_angle > 180) { // && < 540
                std::cout << "1" << std::endl;
                triangles.push_back({get_pid(fh->vertex(0)),
                                      get_pid(fh->vertex(1)),
                                      get_pid(fh->vertex(2))});
            } else if (cumulative_angle < -180) {
                std::cout << "2" << std::endl;
                triangles.push_back({get_pid(fh->vertex(0)),
                                      get_pid(fh->vertex(2)), // invert
                                      get_pid(fh->vertex(1))});
            }
#else
            triangles.push_back({get_pid(fh->vertex(0)),
                                  get_pid(fh->vertex(1)),
                                  get_pid(fh->vertex(2))});
#endif // CGAL_SS3_FILTER_CDT2_FACES_WITH_WINDING_NUMBER
        }
    }
}; // triangulate_facet_with_CDT2()



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



template <typename ValueType,
          typename VisitorBase = CGAL::Polygon_mesh_processing::Autorefinement::Default_visitor>
struct Range_updating_autoref_visitor : public VisitorBase {
    Range_updating_autoref_visitor(const std::vector<ValueType>& old_range,
                                   std::vector<ValueType>& new_range,
                                   const VisitorBase& base = VisitorBase{})
        : VisitorBase(base), old_range_(old_range), new_range_(new_range) {
        new_range.reserve(old_range.size());
    }

    void verbatim_triangle_copy(std::size_t tgt_id, std::size_t src_id) {
        // std::cout << "verbatim_triangle_copy " << tgt_id << " from " << src_id << std::endl;
        VisitorBase::verbatim_triangle_copy(tgt_id, src_id);
        new_range_.resize(tgt_id + 1);
        new_range_[tgt_id] = old_range_[src_id];
    }

    void new_subtriangle (std::size_t tgt_id, std::size_t src_id) {
        // std::cout << "new_subtriangle " << tgt_id << " from " << src_id << std::endl;
        VisitorBase::new_subtriangle(tgt_id, src_id);
        new_range_.resize(tgt_id + 1);
        new_range_[tgt_id] = old_range_[src_id];
    }

private:
    const std::vector<ValueType>& old_range_;
    std::vector<ValueType>& new_range_;
};

template <typename ValueType,
          typename BaseVisitor = CGAL::Polygon_mesh_processing::internal::Default_repair_PS_visitor>
struct Range_updating_repair_PS_visitor : public BaseVisitor {
    Range_updating_repair_PS_visitor(std::vector<ValueType>& range,
                                     const BaseVisitor& base_visitor = BaseVisitor{})
        : BaseVisitor(base_visitor), range_(range) { }

    void swap(std::size_t pos_1, std::size_t pos_2) {
        BaseVisitor::swap(pos_1, pos_2);
        std::swap(range_[pos_1], range_[pos_2]);
    }
    void resize(std::size_t new_size) {
        BaseVisitor::resize(new_size);
        range_.resize(new_size);
    }
private:
    std::vector<ValueType>& range_;
};



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



std::pair<PolyhedronSPtr, CGAL::FT> SimpleStraightSkel::enablePerturbedMode(PolyhedronSPtr polyhedron,
                                                                            CGAL::FT currentOffset,
                                                                            CGAL::FT simultaneousOffset) {
    std::cout << "Enabling perturbed mode..." << std::endl;

    // @fixme the save offset should be taken into account in the limitation of the perturbation:
    // we need perturbed mode to end before the next save offset

    // shift a little bit the polyhedron inwards to avoid self intersections
    CGAL::FT shift = (currentOffset + simultaneousOffset) / 2 - currentOffset; // @fixme complete hack for now
    std::cout << "pre-enable shift = " << shift << std::endl;
    CGAL_assertion(!is_zero(shift));
    polyhedron = algo::_3d::PolyhedronTransformation::shiftFacets(polyhedron, shift);

    db::_3d::OBJFile::save("results/pertubation_pre_enable.obj", polyhedron, false /*do not triangulate*/);

    perturbationOffset_ = currentOffset + shift;

    // perturb
    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
#define CGAL_SS3_PERTURB_PLANE_COEFFICIENTS
#ifdef CGAL_SS3_PERTURB_PLANE_COEFFICIENTS
        facet->storePlaneCoefficients();
        facet->perturbPlaneCoefficients();
        facet->normalizePlaneCoefficients();
#else
        SkelFacetDataSPtr data;
        if (facet->hasData()) {
            data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
        } else {
            data = SkelFacetData::create(facet);
            data->setSpeed(1.0);
        }

        // store
        facet->cachedSpeed_ = data->getSpeed();
        std::cout << "caching speed: " << facet->cachedSpeed_ << std::endl;

        // perturb
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<> dist(0.0, 1e-10);

        data->setSpeed(facet->cachedSpeed_ * (1 + dist(gen)));
        std::cout << "speed is now: " << data->getSpeed() << std::endl;
#endif
    }

    polyhedron = algo::_3d::PolyhedronTransformation::shiftFacets(polyhedron, 0.0);
    db::_3d::OBJFile::save("results/pertubation_post_enable.obj", polyhedron, false /*do not triangulate*/);
    db::_3d::OBJFile::save("results/pertubation_post_enable_triangulated.obj", polyhedron, true /*do not triangulate*/);
    CGAL_assertion(bool(polyhedron));
    CGAL_assertion(polyhedron->isConsistent());

    usingTemporaryPerturbedMode_ = true;
    simultaneousOffset_ = simultaneousOffset;

    return { polyhedron, perturbationOffset_ };
}


// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


std::pair<PolyhedronSPtr, CGAL::FT> SimpleStraightSkel::disablePerturbedMode(PolyhedronSPtr polyhedron,
                                                                             CGAL::FT currentOffset,
                                                                             CGAL::FT nextEventOffset) {
    std::cout << "Disabling perturbed mode..." << std::endl;
    std::cout << "  Current offset: " << currentOffset << std::endl;
    std::cout << "  Next offset: " << nextEventOffset << std::endl;
    CGAL_precondition(usingTemporaryPerturbedMode_);

    // shift a little more to get some waylay
    CGAL::FT shift = (currentOffset + nextEventOffset) / 2 - currentOffset; // @fixme complete hack for now
    std::cout << "pre-disable shift = " << shift << std::endl;
    CGAL_assertion(!is_zero(shift));
    polyhedron = algo::_3d::PolyhedronTransformation::shiftFacets(polyhedron, shift);

    db::_3d::OBJFile::save("results/pertubation_pre_disable.obj", polyhedron, false /*do not triangulate*/);

    CGAL::FT perturbationEndOffset = (currentOffset + shift);
    std::cout << "  Perturbation Start: " << perturbationOffset_ << std::endl;
    std::cout << "  Perturbation End: " << perturbationEndOffset << std::endl;

    // un-perturb
    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
#ifdef CGAL_SS3_PERTURB_PLANE_COEFFICIENTS
        facet->restorePlaneCoefficients(perturbationOffset_, perturbationEndOffset);
#else
        CGAL_assertion(facet->hasData());
        SkelFacetDataSPtr data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
        data->setSpeed(facet->cachedSpeed_);
        std::cout << "speed is back to: " << data->getSpeed() << std::endl;
#endif
    }

    // @fixme what if faces incident to a vertex become a degenerate configuration and we can't
    // recompute the vertex? ==> just continue another iteration in perturbed mode?
    //
    // In this case, need to handle the case of the next offset event being the desired save offset

    polyhedron = algo::_3d::PolyhedronTransformation::shiftFacets(polyhedron, 0.0);
    db::_3d::OBJFile::save("results/pertubation_post_disable.obj", polyhedron, false /*do not triangulate*/);
    CGAL_assertion(bool(polyhedron));

    usingTemporaryPerturbedMode_ = false;
    simultaneousOffset_ = 0;

    return { polyhedron, perturbationEndOffset };
}



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



bool SimpleStraightSkel::handleSaveEventAtSimultaneity(PolyhedronSPtr polyhedron,
                                                       CGAL::FT current_offset,
                                                       CGAL::FT simultaneity_offset) {
    namespace PMP = CGAL::Polygon_mesh_processing;

    std::cout << "handleSaveEventAtSimultaneity()" << std::endl;

    CGAL_precondition(current_offset != simultaneity_offset);

    const CGAL::FT shift = simultaneity_offset - current_offset;
    polyhedron = PolyhedronTransformation::shiftFacets(polyhedron, shift);

    std::stringstream ss_filename;
    ss_filename << "results/offset_" << simultaneity_offset << ".obj";
    db::_3d::OBJFile::save(ss_filename.str(), polyhedron);

    polyhedron->initializeAllIDs();

    // convert to triangulated soup
    std::vector<Point3> points;
    std::vector<std::vector<std::size_t> > triangles;

    // @todo factorize CDT usages but mind the tags
    using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;
    using PK = CGAL::Projection_traits_3<CGAL::K>;
    using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>;
    using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
    using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
    using PTDS = CGAL::Triangulation_data_structure_2<PVb,PFb>;
    using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
    using PCDT_VH = PCDT::Vertex_handle;
    using PCDT_FH = PCDT::Face_handle;

    // not really required because the IDs are initialized in the order of vertex iteration
    std::map<VertexSPtr, std::size_t> v_ids;

    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        unsigned int id = vertex->getID();
        points.emplace_back(vertex->getX(), vertex->getY(), vertex->getZ());
        v_ids[vertex] = id;
    }

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        facet->makeFirstConvex();

        if (facet->edges().size() < 3) {
            std::cerr << "Warning: face with < 3 edges" << std::endl;
            continue;
        } else {
            Vector3SPtr n = KernelFactory::createVector3(facet->plane());
            CGAL_assertion(*n != CGAL::NULL_VECTOR);

            PK traits(*n);
            PCDT pcdt(traits);

            std::map<VertexSPtr, PCDT_VH> face_vhs; // might have multiple vertices at the same position

            std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
            while (it_v != facet->vertices().end()) {
                VertexSPtr vertex = *it_v++;
                auto res = face_vhs.emplace(vertex, PCDT_VH());
                if(res.second) // first time seeing this point
                {
                    PCDT_VH vh = pcdt.insert(*(vertex->getPoint()));
                    vh->info() = vertex;
                    res.first->second = vh;
                }
            }

            auto ne = 0;
            std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
            while (it_e != facet->edges().end()) {
                EdgeSPtr edge = *it_e++;
                VertexSPtr v0 = edge->src(facet);
                VertexSPtr v1 = edge->dst(facet);

                if(*(v0->getPoint()) == *(v1->getPoint()))
                {
                    std::cerr << "W: encountered degenerate edge @ " << *(v0->getPoint()) << std::endl;

                    CGAL_assertion(v0->degree() != 1); // @todo handle that...
                    VertexSPtr vm1 = edge->prev(facet)->src(facet);

                    // create a degenerate face, commented here because we will purge it anyway
                    // triangles.emplace_back({v_ids.at(vm1), v_ids.at(v0), v_ids.at(v1)});
                } else {
                    PCDT_VH vh0 = face_vhs.at(v0);
                    PCDT_VH vh1 = face_vhs.at(v1);

                    try {
                        pcdt.insert_constraint(vh0, vh1);
                    } catch(const typename PCDT::Intersection_of_constraints_exception&) {
                        std::cerr << "Error: Intersection of constraint w/ " << vh0->point() << " " << vh1->point() << std::endl;
                        DEBUG_VAR(facet->toString());
                        CGAL_warning_msg(false, "Intersections in CDT2 not allowed");
                        return false;
                    }
                    ++ne;
                }
            }

            if(ne < 3) { // degenerate face
                std::cerr << "Warning: skipping degenerate face" << std::endl;
                continue;
            }

            std::unordered_map<PCDT_FH, bool> in_domain_map;
            boost::associative_property_map<std::unordered_map<PCDT_FH, bool> > in_domain(in_domain_map);
            CGAL::mark_domain_in_triangulation(pcdt, in_domain);

            for(auto fh : pcdt.finite_face_handles()) {
                if(!get(in_domain, fh)) {
                    continue;
                }

                triangles.push_back({v_ids.at(fh->vertex(0)->info()),
                                     v_ids.at(fh->vertex(1)->info()),
                                     v_ids.at(fh->vertex(2)->info())});

            }
        }
    }

    CGAL::IO::write_OFF("results/simultaneous_save.off", points, triangles);

    // remove degenerate faces
    auto removal_predicate = [&](const std::vector<std::size_t>& polygon) {
        CGAL_precondition(polygon.size() == 3);
        return CGAL::collinear(points[polygon[0]], points[polygon[1]], points[polygon[2]]);
    };
    std::size_t size_pre = triangles.size();
    triangles.erase(std::remove_if(std::begin(triangles), std::end(triangles), removal_predicate),
                   std::end(triangles));
    std::cout << "Removed " << size_pre - triangles.size() << " degenerate faces" << std::endl;

    PMP::autorefine_triangle_soup(points, triangles,
                                  CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag()));
    CGAL::IO::write_OFF("results/simultaneous_save_autoref.off", points, triangles);

    PMP::repair_polygon_soup(points, triangles, CGAL::parameters::erase_all_duplicates(true)
                                                                .require_same_orientation(false)
                                                                .verbose(true));
    CGAL::IO::write_OFF("results/simultaneous_save_autoref_repaired.off", points, triangles);

    PMP::orient_polygon_soup(points, triangles);
    CGAL::IO::write_OFF("results/simultaneous_save_autoref_repaired_oriented.off", points, triangles);

    CGAL_assertion(PMP::does_triangle_soup_self_intersect(points, triangles));

    if (!PMP::is_polygon_soup_a_polygon_mesh(triangles)) {
        std::cerr << "Error: PS not a PM" << std::endl;
        return false;
    }

    return true;
}



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



std::pair<PolyhedronSPtr, CGAL::FT> SimpleStraightSkel::handleEventWithAutoref(AbstractEventSPtr event,
                                                                               CGAL::FT currentOffset,
                                                                               PolyhedronSPtr polyhedron) {
    namespace PMP = CGAL::Polygon_mesh_processing;

    using PID = std::size_t;
    using TID = std::size_t;
    using VID = std::size_t;

    std::cout << "handleEventWithAutoref()" << std::endl;

    CGAL_precondition(!usingTemporaryPerturbedMode_);

    // appendEventNode(event->getNode()); // @todo, if there is a point

    const CGAL::FT shift = event->getOffset() - currentOffset;
    const CGAL::FT delta = 0.2 * shift; // @fixme don't hardcode this value
    const CGAL::FT nudged_shift = shift + delta;
    std::cout << "currentOffset: " << currentOffset << std::endl;
    std::cout << "event->getOffset(): " << event->getOffset() << std::endl;
    std::cout << "shift: " << shift << std::endl;
    std::cout << "delta: " << delta << std::endl;
    std::cout << "total shift: " << nudged_shift << std::endl;

    static int autoref_event_id = -1;
    ++autoref_event_id;

    // {
    //     PolyhedronSPtr polyhedron_tmp = PolyhedronTransformation::shiftFacets(polyhedron, nudged_shift);
    //     db::_3d::OBJFile::save("results/autoref_event-shifted.obj", polyhedron_tmp, false /*do not triangulate*/);
    // }

    auto triangulate_quad_with_CDT2 = [](VertexSPtr vertex, VertexSPtr vertex_offset,
                                         VertexSPtr next_vertex, VertexSPtr next_vertex_offset,
                                         std::map<Point3, std::size_t>& pids,
                                         std::vector<Point3>& points,
                                         std::vector<std::vector<PID> >& triangles) {

        // @todo factorize CDT usages, but mind the tags
        using Itag = CGAL::Exact_intersections_tag;
        using PK = CGAL::Projection_traits_3<CGAL::K>;
        using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>;
        using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
        using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
        using PTDS = CGAL::Triangulation_data_structure_2<PVb,PFb>;
        using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
        using PCDT_VH = PCDT::Vertex_handle;
        using PCDT_FH = PCDT::Face_handle;

        std::cout << "CDT2 QUAD WITH\n";
        std::cout << "vertex = " << *(vertex->getPoint()) << std::endl;
        std::cout << "vertex offset = " << *(vertex_offset->getPoint()) << std::endl;
        std::cout << "next vertex = " << *(next_vertex->getPoint()) << std::endl;
        std::cout << "next vertex offset = " << *(next_vertex_offset->getPoint()) << std::endl;

        CGAL_precondition(*(vertex->getPoint()) != *(vertex_offset->getPoint()));
        CGAL_precondition(*(vertex_offset->getPoint()) != *(next_vertex_offset->getPoint()));
        CGAL_precondition(*(vertex->getPoint()) != *(next_vertex_offset->getPoint()));

        // currently the orientation doesn't matter because we re-order for compatibility
        // with base faces later
        Vector3 n = CGAL::cross_product(Vector3(*(vertex->getPoint()), *(vertex_offset->getPoint())),
                                        Vector3(*(vertex->getPoint()), *(next_vertex_offset->getPoint())));

        CGAL_assertion(n != CGAL::NULL_VECTOR);

        PK projection_traits(n);
        PCDT pcdt(projection_traits);

        Point3SPtr vertex_pt = vertex->getPoint();
        PCDT_VH vertex_vh = pcdt.insert(*vertex_pt);
        vertex_vh->info() = vertex;

        Point3SPtr next_vertex_pt = next_vertex->getPoint();
        PCDT_VH next_vertex_vh = pcdt.insert(*next_vertex_pt);
        next_vertex_vh->info() = next_vertex;

        Point3SPtr vertex_offset_pt = vertex_offset->getPoint();
        PCDT_VH vertex_offset_vh = pcdt.insert(*vertex_offset_pt);
        vertex_offset_vh->info() = vertex_offset;

        Point3SPtr next_vertex_offset_pt = next_vertex_offset->getPoint();
        PCDT_VH next_vertex_offset_vh = pcdt.insert(*next_vertex_offset_pt);
        next_vertex_offset_vh->info() = next_vertex_offset;

        CGAL_precondition(CGAL::orientation(*vertex_pt, *vertex_offset_pt,
                                            *next_vertex_pt, *next_vertex_offset_pt) == CGAL::COPLANAR);

        pcdt.insert_constraint(vertex_vh, next_vertex_vh);
        pcdt.insert_constraint(next_vertex_vh, next_vertex_offset_vh);
        pcdt.insert_constraint(next_vertex_offset_vh, vertex_offset_vh);
        pcdt.insert_constraint(vertex_offset_vh, vertex_vh);

        std::unordered_map<PCDT_FH, bool> in_domain_map;
        boost::associative_property_map<std::unordered_map<PCDT_FH, bool> > in_domain(in_domain_map);
        CGAL::mark_domain_in_triangulation(pcdt, in_domain);

        auto get_pid = [&pids, &points] (PCDT_VH vh) -> PID {
            auto res = pids.emplace(vh->point(), points.size());
            if(res.second) { // first time seeing the vertex handle
                points.push_back(vh->point());
            }
            return res.first->second;
        };

        std::cout << "new triangles from CDT2" << std::endl;

        for(auto fh : pcdt.finite_face_handles()) {
            if(!get(in_domain, fh)) {
                continue;
            }

            triangles.push_back({get_pid(fh->vertex(0)),
                                 get_pid(fh->vertex(1)),
                                 get_pid(fh->vertex(2))});

            std::cout << points[triangles.back()[0]] << " " << points[triangles.back()[1]] << " " << points[triangles.back()[2]] << std::endl;
        }
    }; // lambda 'triangulate_quad_with_CDT2'

    auto fill_edge_map = [](const std::vector<Point3>& points,
                            const std::vector<std::vector<PID> >& triangles,
                            std::vector<std::unordered_map<PID, std::vector<TID> > >& edge_map) {
        CGAL_precondition(edge_map.size() == points.size());

        // collect duplicated edges
        for (TID ti=0; ti<triangles.size(); ++ti) {
            for (std::size_t j=0; j<3; ++j) {
                std::pair<PID, PID> e_pids = CGAL::make_sorted_pair(triangles[ti][j],
                                                                    triangles[ti][(j+1)%3]);
                edge_map[e_pids.first][e_pids.second].push_back(ti);
            }
        }

        namespace pred = CGAL::Polygon_mesh_processing::Corefinement;

        for (std::size_t pid0=0; pid0<points.size(); ++pid0) {
            for (auto& pid1_and_edges : edge_map[pid0]) {
                std::vector<TID>& inc_triangles = pid1_and_edges.second;

                if (inc_triangles.size() == 2) { // edge is only incident to a single SS3 face
                    continue;
                }

                const PID pid1 = pid1_and_edges.first;
                CGAL_assertion(pid0 != pid1);

                std::cout << "processing a non-manifold edge: " << points[pid0] << " -- " << points[pid1] << std::endl;

                auto get_third_point_id = [&triangles, pid0, pid1](TID tid) -> PID
                {
                    std::size_t third;

                    // need to be careful that the orientation of the edge might not match the orientation of the triangle
                    if (triangles[tid][0] == pid0 || triangles[tid][0] == pid1) {
                        if (triangles[tid][1] == pid0 || triangles[tid][1] == pid1) {
                            third = triangles[tid][2];
                        } else {
                            third = triangles[tid][1];
                        }
                    } else {
                        third = triangles[tid][0];
                    }

                    CGAL_postcondition(third != pid0 && third != pid1);
                    return third;
                };

                const Point3& ref_pt = points.at(get_third_point_id(inc_triangles[0]));
                auto less = [&ref_pt, &points, pid0, pid1, get_third_point_id](TID tid1, TID tid2)
                {
                    return pred::sorted_around_edge<CGAL::K>(points.at(pid0), points.at(pid1),
                                                             ref_pt,
                                                             points.at(get_third_point_id(tid1)),
                                                             points.at(get_third_point_id(tid2)));
                };

                std::sort(inc_triangles.begin()+1, inc_triangles.end(), less);

                // std::cout << "Around edge [" << pid0 << " " << pid1 << "], faces are sorted: ";
                // for(TID tid : inc_triangles)
                //   std::cout << " " << tid;
                // std::cout << std::endl;
            }
        }
    }; // lambda 'fill_edge_map'

    enum class Volume_orientation {
        UNKNOWN = 0,
        UNREACHABLE, // 1
        INWARD, // 2
        OUTWARD, // 3
        INCONSISTENT // 4
    };

    auto build_volume_CC = [](const TID seed_tid,
                              const VID CC_ID,
                              const bool start_from_inverted_face,
                              const bool ignore_facet_with_incompatible_orientations,
                              const bool ignore_dangling_outside_facets,
                              const std::vector<Point3>& points,
                              const std::vector<std::vector<PID> >& triangles,
                              const auto& edge_map,
                              auto& volume_CCs,
                              auto& volume_orientations,
                              auto& face_volume_IDs) {

        std::cout << "Building volume #" << CC_ID << " from seed face " << seed_tid << std::endl;

        volume_CCs.emplace_back();
        volume_orientations.emplace_back();

        std::stack<std::pair<TID, bool> > to_visit;
        to_visit.emplace(seed_tid, start_from_inverted_face);

        while (!to_visit.empty())
        {
            TID current_tid;
            bool invert_face;
            std::tie(current_tid, invert_face) = to_visit.top();
            to_visit.pop();

            std::cout << "At face " << current_tid << " [" << triangles[current_tid][0]
                                                   << ", " << triangles[current_tid][1]
                                                   << ", " << triangles[current_tid][2] << "], ";
            std::cout << "invert: " << invert_face << ", ";
            std::cout << "VIDS: " << face_volume_IDs[current_tid][0] << " " << face_volume_IDs[current_tid][1] << std::endl;

            std::size_t pos = invert_face ? 0 : 1;
            if (face_volume_IDs[current_tid][pos] == CC_ID) {
                // already visited this facet during the flooding of this volume's boundary
                continue;
            }

            CGAL_assertion(face_volume_IDs[current_tid][pos] == VID(-1)); // triangle should only be encountered once
            CGAL_warning(face_volume_IDs[current_tid][(pos+1)%2] != CC_ID); // Moebius shenanigans should be an instance of a bug

            volume_CCs.back().push_back(current_tid);

            // mark face as visited
            face_volume_IDs[current_tid][pos] = CC_ID;

            // flood through the edges
            for (int j=0; j<3; ++j) {
                std::pair<PID, PID> e_pids = CGAL::make_sorted_pair(triangles[current_tid][j],
                                                                    triangles[current_tid][(j+1)%3]);
                const std::vector<TID>& inc_triangles = edge_map.at(e_pids.first).at(e_pids.second);
                CGAL_assertion(!inc_triangles.empty());

                std::cout << "  ~~ Crossing edge [" << e_pids.first << ", " << e_pids.second << "]" << std::endl;
                std::cout << "    pos: " << points[e_pids.first] << " " << points[e_pids.second] << std::endl;

                // The faces are ordered CCW while looking from pid0.
                // So the walking while looking from [j] depends on whether [j] is pid0 or not
                int iter_direction = (e_pids.first == triangles[current_tid][j]) ? 1 : -1;
                std::cout << "    iter_direction = " << iter_direction << std::endl;

                // and it also depends on whether we are walking above or below the face
                iter_direction *= invert_face ? 1 : -1;
                std::cout << "    invert_face = " << invert_face << std::endl;

                TID next_tid = current_tid;
                for (;;) {
                    if (inc_triangles.size() == 1) {
                        std::cerr << "Warning: dangling triangle..." << std::endl;
                        std::cout << "    over the edge, the triangle is ITSELF " << current_tid << " [" << triangles[next_tid][0] << ", " << triangles[next_tid][1] << ", " << triangles[next_tid][2] << "], ";
                        to_visit.emplace(current_tid, !invert_face);
                        break;
                    } else if (inc_triangles.size() == 2) {
                        // we should only be there once, meaning if we do not ignore orientations,
                        // then the faces MUST be compatible
                        CGAL_assertion(next_tid == current_tid);

                        next_tid = (inc_triangles[0] == current_tid) ? inc_triangles[1] : inc_triangles[0];
                        std::cout << "    over the edge, the triangle is TRIVIALLY " << next_tid << " [" << triangles[next_tid][0] << ", " << triangles[next_tid][1] << ", " << triangles[next_tid][2] << "], ";
                        std::cout << "VIDS " << face_volume_IDs[next_tid][0] << " " << face_volume_IDs[next_tid][1] << std::endl;
                        CGAL_assertion(next_tid != current_tid);
                    } else {
                        // tricky part, now
                        auto tid_it = std::find(std::begin(inc_triangles), std::end(inc_triangles), next_tid /*updates on every iteration*/);
                        CGAL_assertion(tid_it != inc_triangles.end());

                        if (iter_direction == 1) { // CCW
                            std::cout << "    CCW walk" << std::endl;
                            auto next_it = std::next(tid_it);
                            next_tid = (next_it == inc_triangles.end()) ? inc_triangles[0] : *next_it;
                        } else { // CW
                            std::cout << "    CW walk" << std::endl;
                            next_tid = (tid_it == inc_triangles.begin()) ? inc_triangles.back() : *(std::prev(tid_it));
                        }

                        std::cout << "    over the edge, the triangle is " << next_tid << " [" << triangles[next_tid][0] << ", " << triangles[next_tid][1] << ", " << triangles[next_tid][2] << "], ";
                        std::cout << "VIDS " << face_volume_IDs[next_tid][0] << " " << face_volume_IDs[next_tid][1] << std::endl;
                        CGAL_assertion(next_tid != current_tid);
                    }

                    // If the next face is incident to the outside (CC_ID == 0) on both sides, ignore it
                    if (ignore_dangling_outside_facets) {
                        // CC id #0 is the outside marker for volume building of facet prisms
                        bool is_dangling = (face_volume_IDs[next_tid][0] == 0 && face_volume_IDs[next_tid][1] == 0);
                        if (is_dangling) {
                            CGAL_assertion(inc_triangles.size() != 2);
                            std::cout << "  ignoring next TID because it is dangling" << std::endl;
                            continue;
                        }
                    }

                    // If the edge has the same direction in both faces (aka, the orientation changes),
                    // then we have to flip the direction of turning around the edge)
                    const auto j_it = std::find(std::begin(triangles[next_tid]),
                                              std::end(triangles[next_tid]),
                                              triangles[current_tid][j]);
                    CGAL_assertion(j_it != std::end(triangles[next_tid]));
                    const std::size_t pos = std::distance(std::begin(triangles[next_tid]), j_it);
                    CGAL_assertion(triangles[next_tid][pos] == triangles[current_tid][j]);

                    const bool flip_side = (triangles[next_tid][(pos+1)%3] == triangles[current_tid][(j+1)%3]);
                    std::cout << "    flipping? " << flip_side
                              << " (N: " << triangles[next_tid][(pos+1)%3] << " C: " << triangles[current_tid][(j+1)%3] << ")" << std::endl;
                    if (ignore_facet_with_incompatible_orientations && flip_side) {
                        if (inc_triangles.size() != 2) { // @tmp
                            std::cout << "  ignoring next TID because of incompatible orientation" << std::endl;
                            continue; // keep turning
                        } else {
                            std::cerr << "Warning: we should not have ignored incompatible face orientations, but there are only 2 incident triangles" << std::endl;
                        }
                    }

                    std::cout << "Final TID = " << next_tid << std::endl;
                    if (flip_side) {
                        CGAL_warning(!ignore_facet_with_incompatible_orientations);
                        volume_orientations[CC_ID] = Volume_orientation::INCONSISTENT;
                        to_visit.emplace(next_tid, !invert_face);
                    } else {
                        to_visit.emplace(next_tid, invert_face);
                    }

                    break;
                }
            }
        }
    }; // lambda 'build_volume_CC'

    // -- ACTUAL START

    // polyhedron offset, in a triangle soup form
    std::vector<Point3> points;
    std::vector<std::vector<PID> > triangles;

#ifdef CGAL_SS3_NO_INVERTED_FACE_FILTERING
    obsolete code, purge it
    polyhedron = PolyhedronTransformation::shiftFacets(polyhedron, nudged_shift);

    db::_3d::OBJFile::save("results/autoref_event-shifted.obj", polyhedron, false /*do not triangulate*/);

    polyhedron->initializeAllIDs();

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        triangulate_facet_with_CDT2(facet, triangles);
    }

    std::cout << points.size() << " points (pre auto-refine)" << std::endl;
    std::cout << triangles.size() << " triangles (pref auto-refine)" << std::endl;

    CGAL::IO::write_OFF("results/autoref_event-pre_autoref.off", points, triangles);

    PMP::autorefine_triangle_soup(points, triangles,
                                  CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag()));

    CGAL::IO::write_OFF("results/autoref_event-post_autoref.off", points, triangles);
#else // CGAL_SS3_NO_INVERTED_FACE_FILTERING

    // We need to offset each face, but we only output as correctly oriented the parts that
    // have never been inverted.
    //
    // To detect which parts are inverted, we first build the prism 3D polyhedron going from the base
    // face to the offset face (quadrangular faces to stitch those).
    // Then, we autorefine this triangle soup and identify the volume(s) incident to the base face.
    // Each top face that gets flagged is a face that hasn't been inverted.

    CGAL_assertion(triangles.empty());

    std::vector<int> is_unreachable_triangle; // because we want to swap, and vector<bool> is incompatible
    std::vector<Plane3SPtr> supporting_planes;
    std::vector<CGAL::FT> speeds;

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        std::cout << "\n-- Prism of Facet " << facet->getID() << std::endl;

        facet->makeFirstConvex();

        // Build the offset face
        FacetSPtr facet_offset = facet->clone();
        CGAL::FT speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
        Plane3SPtr offset_plane = KernelWrapper::offsetPlane(facet->plane(), nudged_shift*speed);
        facet_offset->setPlane(offset_plane);

        CGAL_assertion(facet->vertices().size() == facet_offset->vertices().size());

        // abusing the fact that vertices will have the same order in both facets
        std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
        std::list<VertexSPtr>::iterator it_v_offset = facet_offset->vertices().begin();
        while (it_v != facet->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            VertexSPtr offset_vertex = *it_v_offset++;
            CGAL_assertion(vertex->getPoint() == offset_vertex->getPoint()); // haven't offset yet
            Point3SPtr point_offset = PolyhedronTransformation::shiftPoint(vertex, nudged_shift);
            offset_vertex->setPoint(point_offset);
            std::cout << *(vertex->getPoint()) << " offsets to " << *point_offset << std::endl;
        }

        // Triangulate the base & offset faces
        std::vector<Point3> prism_points;
        std::map<Point3, PID> pids;

        std::vector<std::vector<PID> > prism_base_triangles, prism_offset_triangles;
        // odd even is good enough since we know the face is sane
        triangulate_facet_with_CDT2(facet, CDT2_Filtering::ODD_EVEN, pids, prism_points, prism_base_triangles);
        // here, we want all the faces initially, and volume filtering knows which are His own
        triangulate_facet_with_CDT2(facet_offset, CDT2_Filtering::NOT_OUT, pids, prism_points, prism_offset_triangles);

        // invert the offset triangles such that the orientation is opposite that of the base,
        // which will be used later on because we want to enable traversing faces with opposite
        // direction as to avoid erroneous clipping
        for (std::vector<PID>& t : prism_offset_triangles) {
            std::swap(t[0], t[1]);
        }

        std::cout << prism_base_triangles.size() << " triangles (base)" << std::endl;
        std::cout << prism_offset_triangles.size() << " triangles (offset)" << std::endl;

        // lateral faces link the base and the offset triangles
        std::vector<std::vector<PID> > prism_lateral_triangles;

        std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
        std::list<EdgeSPtr>::iterator it_e_offset = facet_offset->edges().begin();
        while (it_e != facet->edges().end()) {
            EdgeSPtr edge = *it_e++;
            EdgeSPtr edge_offset = *it_e_offset++;
            VertexSPtr vertex = edge->src(facet);
            VertexSPtr vertex_offset = edge_offset->src(facet_offset);
            VertexSPtr next_vertex = edge->dst(facet);
            VertexSPtr next_vertex_offset = edge_offset->dst(facet_offset);

            // whatever the diagonal, we might have a bow tie, and we don't want a fold
            // because the sort_around_edge() predicate cannot handle folds
            std::vector<std::vector<PID> > local_triangles;
            triangulate_quad_with_CDT2(vertex, vertex_offset, next_vertex, next_vertex_offset,
                                       pids, prism_points, local_triangles);

            // Re-orient triangles so that they match the boundaries of base & offset faces.
            // Orientations can be flipped because the edges of inner holes can be badly oriented.
            bool flip_orientations = false;

            // 'prism_base_triangles' are correctly oriented, and we need the lateral triangle
            // to have compatible orientations.
            // @todo do something less brute force
            for (const std::vector<PID>& lt : local_triangles) {
                for (std::size_t lj=0; lj<3; ++lj) {
                    for (const std::vector<PID>& bt : prism_base_triangles) {
                        for (std::size_t bj=0; bj<3; ++bj) {
                            if (prism_points[lt[lj]] == prism_points[bt[bj]] &&
                                prism_points[lt[(lj+1)%3]] == prism_points[bt[(bj+1)%3]]) {
                                std::cout << " incompatible orientations on " << prism_points[lt[lj]] << " " << prism_points[lt[(lj+1)%3]] << std::endl;
                                flip_orientations = true;
                                break;
                            }
                        }
                        if (flip_orientations) { break; }
                    }
                    if (flip_orientations) { break; }
                }
                if (flip_orientations) { break; }
            }

            if (flip_orientations) {
                for (std::vector<PID>& lt : local_triangles) {
                    std::swap(lt[0], lt[1]);
                }
            }

            prism_lateral_triangles.insert(std::end(prism_lateral_triangles),
                                           std::cbegin(local_triangles), std::cend(local_triangles));
        }

        std::cout << prism_lateral_triangles.size() << " triangles (lateral)" << std::endl;

        enum class Triangle_prism_location {
            UNKNOWN = 0,
            LATERAL, // = 1
            BASE, // = 2
            OFFSET, // = 3
        };

        // these are the ranges for the merge of base + offset + lateral
        std::vector<std::vector<PID> > prism_triangles;
        std::vector<Triangle_prism_location> triangle_locations;

        prism_triangles.insert(std::end(prism_triangles),
                               std::begin(prism_offset_triangles), std::end(prism_offset_triangles));
        triangle_locations.resize(prism_triangles.size(), Triangle_prism_location::OFFSET);

        CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-face-" + std::to_string(facet->getID()) + "-top.off",
                            prism_points, prism_triangles);

        prism_triangles.insert(std::end(prism_triangles),
                               std::begin(prism_base_triangles), std::end(prism_base_triangles));
        triangle_locations.resize(prism_triangles.size(), Triangle_prism_location::BASE);

        // ESSENTIAL: lateral faces are guaranted to be planar quads because the vertices
        // move in the bisecting plane of the two planes incident to the edge incident to the two vertices.
        prism_triangles.insert(std::end(prism_triangles),
                               std::begin(prism_lateral_triangles), std::end(prism_lateral_triangles));
        triangle_locations.resize(prism_triangles.size(), Triangle_prism_location::LATERAL);

        // @tmp debug
        {
            std::cout << prism_points.size() << " points in shifting prism (pre autorefine)" << std::endl;
            std::cout << prism_triangles.size() << " triangles in shifting prism (pre autorefine)" << std::endl;

            std::cout << prism_points.size() << " prism points DEBUG" << std::endl;
            std::cout << prism_triangles.size() << " prism triangles DEBUG" << std::endl;

            for (const auto& p : prism_points) {
                std::cout << p << std::endl;
            }
            for (const auto& t : prism_triangles) {
                std::cout << t[0] << " " << t[1] << " " << t[2] << std::endl;
            }

            for (std::size_t i=0; i<prism_triangles.size(); ++i) {
                CGAL_assertion(triangle_locations[i] != Triangle_prism_location::UNKNOWN);
                std::cout << "triangle location[" << i << "] = " << int(triangle_locations[i]) << std::endl;
            }

            CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-face-" + std::to_string(facet->getID()) + "-pre_autoref.off",
                                prism_points, prism_triangles);
        }

        // normally we shouldn't need to autoref the base & offset, but if autoref merges points (does it?),
        // it would mess up the other ranges
        //
        // And because we autorefine everything, we need a visitor because autoref might re-order
        // bottom and top faces even though they will not be refined (since we have already
        // computed a CDT2 and inserted intersections.)
        std::vector<Triangle_prism_location> updated_triangle_locations;
        Range_updating_autoref_visitor<Triangle_prism_location> autoref_visitor(triangle_locations, updated_triangle_locations);

        PMP::autorefine_triangle_soup(prism_points, prism_triangles,
                                      CGAL::parameters::visitor(autoref_visitor)
                                                       .concurrency_tag(CGAL::Parallel_if_available_tag()));

        triangle_locations = std::move(updated_triangle_locations);

        CGAL_assertion(triangle_locations.size() == prism_triangles.size());
        for (std::size_t i=0; i<prism_triangles.size(); ++i) {
            CGAL_assertion(triangle_locations[i] != Triangle_prism_location::UNKNOWN);
        }

        // because we could have overlaps and that's troublesome to determine volumes
        // @fixme does it need to be erase_all_duplicates(*TRUE*) sometimes?
        // @fixme the visitor needs to be smart about which face should be purged or not to have
        //        the correct value of reachable
        Range_updating_repair_PS_visitor<Triangle_prism_location> repair_ps_visitor(triangle_locations);
        PMP::merge_duplicate_polygons_in_polygon_soup(prism_points, prism_triangles,
                                                      CGAL::parameters::visitor(repair_ps_visitor)
                                                                       .erase_all_duplicates(false) /*keep one*/
                                                                       .require_same_orientation(false));

        std::cout << prism_points.size() << " points in shifting prism (post autorefine)" << std::endl;
        std::cout << prism_triangles.size() << " triangles in shifting prism (post autorefine)" << std::endl;

        CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-face-" + std::to_string(facet->getID()) + "-post_autoref.off",
                            prism_points, prism_triangles);

        // ----

        std::cout << "Walking into the prism offset of facet " << facet->getID() << std::endl;

        // extremity --> other extremity --> range of incident faces
        std::vector<std::unordered_map<PID, std::vector<TID> > > edge_map(prism_points.size());
        fill_edge_map(prism_points, prism_triangles, edge_map);

        // identify volumes in the shifting faces soup, and tag faces of the volumes
        // that are incident to the base face(s)
        std::vector<std::vector<TID> > unused_volume_CCs; // range of range (volume) of triangle IDs
        std::vector<Volume_orientation> unused_volume_orientations;
        std::vector<std::array<VID, 2> > face_volume_IDs(prism_triangles.size(),
                                                         // [0] is down, [1] is up
                                                         std::array<VID, 2>{VID(-1), VID(-1)});

        // need to flood from all possible base faces because the faces is not necessarily a single CC
        //
        // Two independent loops for clarity and safety:
        // - Flood outside
        // - Flood inner volumes

        // first loop, outside
        for(std::size_t i=0; i<prism_triangles.size(); ++i) {
            if (triangle_locations[i] != Triangle_prism_location::BASE ||
                face_volume_IDs[i][1] != VID(-1)) {
                continue;
            }

            // outer flooding starts 'up' because the base face points out
            build_volume_CC(i /*tid*/, 0 /*CC ID*/, false /*do not use inverted orientation*/,
                            false /*do not traverse boundaries with wrong orientations*/,
                            false /*ignore dangling faces while walking*/,
                            prism_points, prism_triangles, edge_map,
                            unused_volume_CCs, unused_volume_orientations, face_volume_IDs);
        }

        // second loop, inner volumes
        for(std::size_t i=0; i<prism_triangles.size(); ++i) {
            if (triangle_locations[i] != Triangle_prism_location::BASE ||
                face_volume_IDs[i][0] != VID(-1)) {
                continue;
            }

            // inner flooding starts 'down' because the base face points out
            build_volume_CC(i /*tid*/, 1 /*CC ID*/, true /*use inverted orientation*/,
                            true /*traverse boundaries with wrong orientations*/,
                            true /*ignore dangling faces while walking*/,
                            prism_points, prism_triangles, edge_map,
                            unused_volume_CCs, unused_volume_orientations, face_volume_IDs);

            CGAL_assertion(face_volume_IDs[i][0] == VID(1)); // inward
            CGAL_assertion(face_volume_IDs[i][1] == VID(0)); // outward
        }


        // just for debugging
        // {
        //     for(VID cc_id=0; cc_id<unused_volume_CCs.size(); ++cc_id) {
        //         std::cout << "Volume #" << cc_id
        //                   << ", size: " << unused_volume_CCs[cc_id].size() << std::endl;

        //         std::vector<std::vector<PID> > CC_triangles;
        //         for (TID tid : unused_volume_CCs[cc_id]) {
        //             std::cout << tid << " ";
        //             CC_triangles.push_back(prism_triangles[tid]);
        //         }
        //         std::cout << std::endl;
        //         CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-face-" + std::to_string(facet->getID()) +  "-volume_CC-" + std::to_string(cc_id) + ".off", prism_points, CC_triangles);
        //     }
        // }

        std::vector<std::vector<PID> > prism_triangle_contributions; // @tmp debug only
        std::vector<std::vector<PID> > prism_triangle_reachable_contributions; // @tmp debug only

        // Now, flag the offset faces that are not reachable
        for (std::size_t i=0; i<prism_triangles.size(); ++i) {
            CGAL_assertion(triangle_locations[i] != Triangle_prism_location::UNKNOWN);
            // std::cout << "triangle [" << i << "] position = " << prism_points[prism_triangles[i][0]] << " " << prism_points[prism_triangles[i][1]] << " " << prism_points[prism_triangles[i][2]] << std::endl;
            // std::cout << "triangle location[" << i << "] = " << int(triangle_locations[i]) << std::endl;

            if (triangle_locations[i] != Triangle_prism_location::OFFSET) {
                continue;
            }

            const bool is_dangling = (face_volume_IDs[i][0] == 0 /*exterior*/ &&
                                      face_volume_IDs[i][1] == 0 /*exterior*/);

            // we have walked inner volumes only, so check if the _lower_ side of offset faces is reached
            const bool is_unreached = (face_volume_IDs[i][0] == VID(-1));
            std::cout << "Triangle [" << i << "] reached? " << !is_unreached << std::endl;

#ifdef CGAL_SS3_DO_NOT_ADD_UNREACHED_TRIANGLES_TO_CONTRIBUTIONS
            this is too aggressive because we do not get volumes in the union of contributions
            and it is not clear if dangling faces must be extended or purged
            if (!is_unreached) {
#endif
                if(!is_dangling) {
                    // top's orientation was flipped to create a real prism volume, restore the real orientation
                    triangles.push_back({points.size() + prism_triangles[i][0],
                                         points.size() + prism_triangles[i][2],
                                         points.size() + prism_triangles[i][1]});

                    speeds.push_back(speed);
                    supporting_planes.push_back(offset_plane);
                    is_unreachable_triangle.push_back(is_unreached);
                    prism_triangle_contributions.push_back(prism_triangles[i]);
                }

                if (!is_unreached) {
                    prism_triangle_reachable_contributions.push_back(prism_triangles[i]);
                }

#ifdef CGAL_SS3_DO_NOT_ADD_UNREACHED_TRIANGLES_TO_CONTRIBUTIONS
            }
#endif
        }

        CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-face-" + std::to_string(facet->getID()) + "-contributions.off",
                            prism_points, prism_triangle_contributions);

        CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-face-" + std::to_string(facet->getID()) + "-reachable_contributions.off",
                            prism_points, prism_triangle_reachable_contributions);

        points.insert(std::end(points), std::cbegin(prism_points), std::cend(prism_points));

        std::cout << "Now, " << triangles.size() << " triangles in the polyhedron offset" << std::endl;
    }

    std::cout << " == MAIN POLYHEDRON OFFSET AUTOREF / VOLUME CHECKS ==" << std::endl;

    CGAL_postcondition(triangles.size() == speeds.size());
    CGAL_postcondition(triangles.size() == supporting_planes.size());
    CGAL_postcondition(triangles.size() == is_unreachable_triangle.size());
    for (std::size_t i=0; i<triangles.size(); ++i) {
        std::cout << "triangle [" << i << "] positions = " << points[triangles[i][0]] << " " << points[triangles[i][1]] << " " << points[triangles[i][2]] << std::endl;
        std::cout << "is_unreached[" << i << "] = " << is_unreachable_triangle[i] << " (before autoref)" << std::endl;
    }

    std::cout << points.size() << " points (pre auto-refine)" << std::endl;
    std::cout << triangles.size() << " triangles (pref auto-refine)" << std::endl;

    CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-pre_autoref.off", points, triangles);

    std::vector<int> updated_is_unreachable_triangle;
    std::vector<Plane3SPtr> updated_supporting_planes;
    std::vector<CGAL::FT> updated_speeds;
    Range_updating_autoref_visitor<int> bbvis(is_unreachable_triangle, updated_is_unreachable_triangle);
    Range_updating_autoref_visitor<Plane3SPtr, decltype(bbvis)> bvis(supporting_planes, updated_supporting_planes, bbvis);
    Range_updating_autoref_visitor<CGAL::FT, decltype(bvis)> vis(speeds, updated_speeds, bvis);

    PMP::autorefine_triangle_soup(points, triangles,
                                  CGAL::parameters::visitor(vis)
                                                   .concurrency_tag(CGAL::Parallel_if_available_tag()));

    is_unreachable_triangle = std::move(updated_is_unreachable_triangle);
    supporting_planes = std::move(updated_supporting_planes);
    speeds = std::move(updated_speeds);

    Range_updating_repair_PS_visitor<int> bbvis2(is_unreachable_triangle);
    Range_updating_repair_PS_visitor<Plane3SPtr, decltype(bbvis2)> bvis2(supporting_planes, bbvis2);
    Range_updating_repair_PS_visitor<CGAL::FT, decltype(bvis2)> vis2(speeds, bvis2);

    PMP::merge_duplicate_polygons_in_polygon_soup(points, triangles,
                                                  CGAL::parameters::visitor(vis2)
                                                                   .erase_all_duplicates(false) /*keep one*/
                                                                   .require_same_orientation(false));

    CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-post_autoref.off", points, triangles);
#endif // CGAL_SS3_NO_INVERTED_FACE_FILTERING

    std::cout << points.size() << " points (post auto-refine)" << std::endl;
    std::cout << triangles.size() << " triangles (post auto-refine)" << std::endl;

    CGAL_postcondition(triangles.size() == speeds.size());
    CGAL_postcondition(triangles.size() == supporting_planes.size());
    CGAL_postcondition(triangles.size() == is_unreachable_triangle.size());
    for (std::size_t i=0; i<triangles.size(); ++i) {
        std::cout << "triangle [" << i << "] position = " << points[triangles[i][0]] << " " << points[triangles[i][1]] << " " << points[triangles[i][2]] << std::endl;
        std::cout << "is_unreached[" << i << "] = " << is_unreachable_triangle[i] << " (after autoref)" << std::endl;
    }

    // ----------
    // Now, we have built a triangulation of all the offset faces of the polyhedron,
    // with proper orientation and flags
    //
    // What's next is identifying the volumes within this offset polyhedron, and getting rid
    // of the volumes that are not interesting

    // extremity --> other extremity --> range of incident faces
    std::vector<std::unordered_map<PID, std::vector<TID> > > edge_map(points.size());
    fill_edge_map(points, triangles, edge_map);

    // split the volumetric partition into volume connected components
    std::vector<std::vector<TID> > volume_CCs; // range of range (volume) of triangle IDs
    std::vector<Volume_orientation> volume_orientations;
    std::vector<std::array<VID, 2> > face_volume_IDs(triangles.size(),
                                                     // [0] is down, [1] is up
                                                     std::array<VID, 2>{VID(-1), VID(-1)});

#ifdef CGAL_SS3_FILTER_VOLUMES_WITH_ONLY_REACHABLE_FACES
    VID CC_ID = 0;
    for (std::size_t i=0; i<triangles.size(); ++i) {
        if (face_volume_IDs[i][0] == VID(-1) && // 'down' as we have inverted the offset triangles
            !is_unreachable_triangle[i]) {
            build_volume_CC(i /*seed ID*/, CC_ID++, true /*inverted orientation (down)*/,
                            false /*traverse boundaries with wrong orientations*/,
                            false /*ignore dangling faces while walking*/,
                            points, triangles, edge_map,
                            volume_CCs, volume_orientations, face_volume_IDs);
        }
    }

    for(VID cc_id=0; cc_id<volume_CCs.size(); ++cc_id) {
        for (TID tid : volume_CCs[cc_id]) {
            if (is_unreachable_triangle[tid]) {
                volume_orientations[cc_id] = Volume_orientation::UNREACHABLE;
                break;
            }
        }
    }

    // debug-only loop
    for(VID cc_id=0; cc_id<volume_CCs.size(); ++cc_id) {
        std::cout << "Volume #" << cc_id
                  << ", size: " << volume_CCs[cc_id].size()
                  << ", flag: " << int(volume_orientations[cc_id]) << std::endl;

        std::vector<std::vector<PID> > CC_triangles;
        for (TID tid : volume_CCs[cc_id]) {
            CC_triangles.push_back(triangles[tid]);
        }

        CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-volume_CC-" + std::to_string(cc_id) + ".off", points, CC_triangles);
    }

    if (true || event->getType() == AbstractEvent::SURFACE_EVENT) {
        for (std::size_t pid0=0; pid0<points.size(); ++pid0) {
            for (auto& pid1_and_edges : edge_map[pid0]) {
                const PID pid1 = pid1_and_edges.first;

                std::vector<TID>& inc_triangles = pid1_and_edges.second;
                if (inc_triangles.size() != 2) {
                    continue;
                }

                // convex | concave angle
                TID tid0 = inc_triangles[0], tid1 = inc_triangles[1];

                auto third_pos = [&triangles](const TID tid, const TID other) -> std::size_t {
                    for (std::size_t i=0; i<3; ++i)
                        if (std::find(std::cbegin(triangles[other]), std::cend(triangles[other]),
                                      triangles[tid][i]) == std::cend(triangles[other]))
                            return i;
                    return -1; // In case no such element is found
                };

                std::size_t third_pos0 = third_pos(tid0, tid1), third_pos1 = third_pos(tid1, tid0);
                PID pid = triangles[tid0][(third_pos0 + 1)%3];
                PID qid = triangles[tid0][(third_pos0 + 2)%3];
                PID rid = triangles[tid0][third_pos0];
                PID sid = triangles[tid1][third_pos1];
                CGAL_assertion((pid != qid) && (pid != rid) && (pid != sid) && (qid != rid) && (qid != sid) && (rid != sid));
                CGAL::FT dh = CGAL::approximate_dihedral_angle(points[pid], points[qid], points[rid], points[sid]);

                std::cout << "checking checking" << std::endl;
                if (is_unreachable_triangle[tid0] != is_unreachable_triangle[tid1]) {
                    std::cerr << "deg 2 edge with reachable/unreachable triangles @ edge ";
                    std::cerr << points[pid0] << " -- " << points[pid1] << std::endl;
                }
            }
        }
    }

    // end debug

    // Purge faces incident to volumes either unreachable or inconsistent
    auto removal_predicate = [&](const TID tid) {
        Volume_orientation down_flag = (face_volume_IDs[tid][0] != VID(-1)) ? volume_orientations.at(face_volume_IDs[tid][0]) : Volume_orientation::UNREACHABLE;
        Volume_orientation up_flag = (face_volume_IDs[tid][1] != VID(-1)) ? volume_orientations.at(face_volume_IDs[tid][1]) : Volume_orientation::UNREACHABLE;

        const bool down_rejected = (down_flag == Volume_orientation::UNREACHABLE || down_flag == Volume_orientation::INCONSISTENT);
        const bool up_rejected = (up_flag == Volume_orientation::UNREACHABLE || up_flag == Volume_orientation::INCONSISTENT);

        return (down_rejected && up_rejected);
    };
#else
    this is bad because if there are holes, the volume will be negative
    but we do not want to discard it

    VID CC_ID = 0;
    for (std::size_t i=0; i<triangles.size(); ++i) {
        std::cout << "face " << i << " VIDS: " << face_volume_IDs[i][0] << " " << face_volume_IDs[i][1] << std::endl;
        if (face_volume_IDs[i][0] == VID(-1)) {
            build_volume_CC(i /*seed ID*/, CC_ID++, true /*inverted orientation (down)*/,
                            false /*traverse boundaries with wrong orientations*/,
                            false /*ignore dangling faces while walking*/,
                            points, triangles, edge_map,
                            volume_CCs, volume_orientations, face_volume_IDs);
        }

        if (face_volume_IDs[i][1] == VID(-1)) {
            build_volume_CC(i /*seed ID*/, CC_ID++, false /*normal orientation (up)*/,
                            false /*traverse boundaries with wrong orientations*/,
                            false /*ignore dangling faces while walking*/,
                            points, triangles, edge_map,
                            volume_CCs, volume_orientations, face_volume_IDs);
        }


        CGAL_postcondition(face_volume_IDs[i][0] != VID(-1) && face_volume_IDs[i][1] != VID(-1));
    }

    // @fixme some volumes are computed twice (take, e.g., a sole tetrahedron)
    std::cout << CC_ID << " volumes in the offset polyhedron" << std::endl;

    // ----------

    // determine the orientation of consistent CCs
    for(VID cc_id=0; cc_id<volume_CCs.size(); ++cc_id) {
        // just for debugging
        {
            std::cout << "Volume #" << cc_id
                      << ", size: " << volume_CCs[cc_id].size()
                      << ", orientation: " << int(volume_orientations[cc_id]) << std::endl;

            std::vector<std::vector<PID> > CC_triangles;
            for (TID tid : volume_CCs[cc_id]) {
                CC_triangles.push_back(triangles[tid]);
            }
            CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-volume_CC-" + std::to_string(cc_id) + ".off", points, CC_triangles);
        }

        if (volume_orientations[cc_id] != Volume_orientation::UNKNOWN) {
            continue;
        }

#ifndef CGAL_SS3_NO_INVERTED_FACE_FILTERING
        bool reject_due_to_incident_unreachable_triangle = false;
        for (TID tid : volume_CCs[cc_id]) {
#ifdef CGAL_SS3_DO_NOT_ADD_UNREACHED_TRIANGLES_TO_CONTRIBUTIONS
            CGAL_assertion(!is_unreachable_triangle[tid]);
# else
            if (is_unreachable_triangle[tid]) {
                std::cout << "Rejecting because " << tid << " is NOT reachable" << std::endl;
                volume_orientations[cc_id] = Volume_orientation::UNREACHABLE;
                reject_due_to_incident_unreachable_triangle = true;
            }

            if (reject_due_to_incident_unreachable_triangle) {
                break;
            }
# endif
        }

        if (reject_due_to_incident_unreachable_triangle) {
            continue;
        }
#endif

        // the volume CC has consistent orientations, but we need to know if it's inward or outward oriented

        // @fixme for now, evaluating volumes because finding the equivalent of:
        //    target(next(min_slope_he, pmesh), pmesh))
        //    target(next(opposite(min_slope_he, pmesh)
        // in:
        //    PMP::internal::is_outward_oriented
        // is tedious
#if 0
        // find the extremum point
        PID ext_pid = -1;

        for (TID tid : volume_CCs[cc_id]) {
            for (std::size_t j=0; j<3; ++j) {
                PID pid = triangles[tid][j];
                if (ext_pid == PID(-1) ||
                    points[ext_pid].z() < points[pid].z()) {
                    ext_pid = pid;
                }
            }
        }

        CGAL_postcondition(ext_pid != PID(-1));

        // find the edges adjacent to the extrem point
        std::unordered_set<PID> incident_vertices;
        for (TID tid : volume_CCs[cc_id]) {
            for (std::size_t j=0; j<3; ++j) {
                if (triangles[tid][j] == ext_pid) {
                  incident_vertices.insert(triangles[tid][(j+1)%3]);
                  incident_vertices.insert(triangles[tid][(j+2)%3]); // should be sufficient to only do one
                }
            }
        }

        // now we compare slopes, but we need the third points on either side of the edge.......
        ... @todo
#else
        // outward = positive volume, inward = negative volume
        Point3 origin(0,0,0);
        CGAL::FT volume = 0;
        CGAL::internal::Evaluate<CGAL::FT> evaluate;
        for (TID tid : volume_CCs[cc_id]) {
            volume += CGAL::volume(origin,
                                   points[triangles[tid][0]],
                                   points[triangles[tid][1]],
                                   points[triangles[tid][2]]);
            evaluate(volume);
        }

        volume_orientations[cc_id] = (volume > 0) ? Volume_orientation::OUTWARD
                                                  : Volume_orientation::INWARD;

        std::cout << "Volume #" << cc_id
                  << ", volume: " << volume
                  << ", orientation: " << int(volume_orientations[cc_id]) << std::endl;
#endif
    }

    // ----------

    // Purge faces that are not incident to an "interior" volume
    auto removal_predicate = [&](const TID tid) {
        const Volume_orientation first_orientation = volume_orientations[face_volume_IDs[tid][0]];
        const Volume_orientation second_orientation = volume_orientations[face_volume_IDs[tid][1]];

        CGAL_warning(!(volume_orientations[face_volume_IDs[tid][0]] == Volume_orientation::OUTWARD &&
                       volume_orientations[face_volume_IDs[tid][1]] == Volume_orientation::OUTWARD));
        CGAL_warning(!(volume_orientations[face_volume_IDs[tid][0]] == Volume_orientation::INWARD &&
                       volume_orientations[face_volume_IDs[tid][1]] == Volume_orientation::INWARD));

        // discard any face that is not incident to at least one non-rejected volume
        return (first_orientation != Volume_orientation::OUTWARD && second_orientation != Volume_orientation::OUTWARD);
    };
#endif // CGAL_SS3_FILTER_VOLUMES_WITH_ONLY_REACHABLE_FACES

    for (TID tid = 0; tid < triangles.size();) {
        if (removal_predicate(tid)) {
            std::size_t swap_pos = triangles.size() - 1;
            if (tid != swap_pos) {
                std::swap(face_volume_IDs[tid], face_volume_IDs[swap_pos]);
                std::swap(triangles[tid], triangles[swap_pos]);
                std::swap(supporting_planes[tid], supporting_planes[swap_pos]);
                std::swap(speeds[tid], speeds[swap_pos]);
            }
            triangles.pop_back();
            supporting_planes.pop_back();
            speeds.pop_back();
        } else {
            ++tid;
        }
    }

    std::cout << "final has " << triangles.size() << " triangles" << std::endl;

    PMP::remove_isolated_points_in_polygon_soup(points, triangles);

    CGAL::IO::write_OFF("results/autoref_event-" + std::to_string(autoref_event_id) + "-final_soup.off", points, triangles);

    // ----------

    polyhedron = soup_to_polyhedron(points, triangles, supporting_planes, speeds);
    CGAL_assertion(polyhedron && polyhedron->isConsistent());

    db::_3d::AbstractFile::mergeCoplanarFacets(polyhedron);
    db::_3d::AbstractFile::removeVerticesDegLt3(polyhedron);
    CGAL_postcondition(polyhedron && polyhedron->isConsistent());

    db::_3d::OBJFile::save("results/autoref_event-almost_final.obj", polyhedron, false /*do not triangulate*/);

    // @fixme this is very dangerous because now we will get a different normalization for the planes.
    // It would be better to re-use the planes
    PolyhedronTransformation::harmonizeFacetPlanes(polyhedron);
    polyhedron = algo::_3d::PolyhedronTransformation::shiftFacets(polyhedron, 0.0);
    CGAL_postcondition(polyhedron && polyhedron->isConsistent());

    polyhedron->initializeAllIDs();

    db::_3d::OBJFile::save("results/autoref_event-final.obj", polyhedron, false /*do not triangulate*/);

#ifndef CGAL_SS3_NO_SKELETON_DS
    event->setPolyhedronResult(polyhedron);
#endif
    skel_result_->addEvent(event);

    return { polyhedron, event->getOffset() + delta };
}




// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //




/*
  vanishAt:

      \                      /
       \       F0           /
        \                  /
   F1    ------------------   F3
        /                  \
       /       F2           \
      /                      \

  crashAt:

      \                      /                  \                      /
       \       F0           /                    \       F2           /
        \                  /                      \                  /
         ------------------                          ----------------
        /                  \                      /                  \
       /       F1           \                    /       F3           \
      /                      \                  /                      \
*/
enum Quadplane_parallelism {
    PLANES_PARALLELISM_NONE = 0, // No planes are parallel
    PLANES_PARALLEL_01,          // Planes 0 and 1 are parallel (1)
    PLANES_PARALLEL_02,          // Planes 0 and 2 are parallel (2)
    PLANES_PARALLEL_03,          // Planes 0 and 3 are parallel (3)
    PLANES_PARALLEL_12,          // Planes 1 and 2 are parallel (4)
    PLANES_PARALLEL_13,          // Planes 1 and 3 are parallel (5)
    PLANES_PARALLEL_23,          // Planes 2 and 3 are parallel (6)
    PLANES_PARALLEL_012,         // Planes 0, 1, and 2 are parallel (7)
    PLANES_PARALLEL_013,         // Planes 0, 1, and 3 are parallel (8)
    PLANES_PARALLEL_023,         // Planes 0, 2, and 3 are parallel (9)
    PLANES_PARALLEL_123,         // Planes 1, 2, and 3 are parallel (10)
    PLANES_PARALLEL_01_23,       // Planes 0 and 1 are parallel, and planes 2 and 3 are parallel (11)
    PLANES_PARALLEL_02_13,       // Planes 0 and 2 are parallel, and planes 1 and 3 are parallel (12)
    PLANES_PARALLEL_03_12,       // Planes 0 and 3 are parallel, and planes 1 and 2 are parallel (13)
    PLANES_PARALLEL_0123         // All planes are parallel (14)
};

Quadplane_parallelism quadplane_parallelism(FacetSPtr facet_0,
                                            FacetSPtr facet_1,
                                            FacetSPtr facet_2,
                                            FacetSPtr facet_3)
{
    Plane3SPtr plane_0 = facet_0->getPlane();
    Plane3SPtr plane_1 = facet_1->getPlane();
    Plane3SPtr plane_2 = facet_2->getPlane();
    Plane3SPtr plane_3 = facet_3->getPlane();

    // Function to check if two planes are parallel
    auto are_planes_parallel = [](const Plane3SPtr pl1, const Plane3SPtr pl2) {
        // plane coefficients are normalized
        return ((pl1->a() == pl2->a() && pl1->b() == pl2->b() && pl1->c() == pl2->c()) ||
                (pl1->a() == - pl2->a() && pl1->b() == - pl2->b() && pl1->c() == - pl2->c()));
    };

    // @todo could avoid some computations: 1&2 + 2&3 ==> no need to do 1&3,
    // but it would be a lot more code and branches
    int mask = (are_planes_parallel(plane_0, plane_1) << 5) |
               (are_planes_parallel(plane_0, plane_2) << 4) |
               (are_planes_parallel(plane_0, plane_3) << 3) |
               (are_planes_parallel(plane_1, plane_2) << 2) |
               (are_planes_parallel(plane_1, plane_3) << 1) |
               (are_planes_parallel(plane_2, plane_3));

    // Switch based on the bitmask to return the appropriate enum
    switch (mask) {
        case 0b111111: return PLANES_PARALLEL_0123;    // All planes are parallel
        case 0b111000: return PLANES_PARALLEL_012;     // 0, 1, 2 are parallel
        case 0b101100: return PLANES_PARALLEL_013;     // 0, 1, 3 are parallel
        case 0b011010: return PLANES_PARALLEL_023;     // 0, 2, 3 are parallel
        case 0b000111: return PLANES_PARALLEL_123;     // 1, 2, 3 are parallel
        case 0b100001: return PLANES_PARALLEL_01_23;   // 0, 1 are parallel, and 2, 3 are parallel
        case 0b001001: return PLANES_PARALLEL_02_13;   // 0, 2 are parallel, and 1, 3 are parallel
        case 0b010010: return PLANES_PARALLEL_03_12;   // 0, 3 are parallel, and 1, 2 are parallel
        case 0b100000: return PLANES_PARALLEL_01;      // 0, 1 are parallel
        case 0b010000: return PLANES_PARALLEL_02;      // 0, 2 are parallel
        case 0b001000: return PLANES_PARALLEL_03;      // 0, 3 are parallel
        case 0b000100: return PLANES_PARALLEL_12;      // 1, 2 are parallel
        case 0b000010: return PLANES_PARALLEL_13;      // 1, 3 are parallel
        case 0b000001: return PLANES_PARALLEL_23;      // 2, 3 are parallel
        default:       return PLANES_PARALLELISM_NONE; // No planes are parallel
    }
}



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



std::pair<Point3SPtr, CGAL::FT> SimpleStraightSkel::vanishesAtOnePairOpposite(FacetSPtr facet_0,
                                                                              FacetSPtr facet_1,
                                                                              FacetSPtr facet_2,
                                                                              FacetSPtr facet_3)
{
    Point3SPtr point = Point3SPtr();

    Plane3SPtr plane_0 = facet_0->plane();
    Plane3SPtr plane_1 = facet_1->plane();
    Plane3SPtr plane_2 = facet_2->plane();
    Plane3SPtr plane_3 = facet_3->plane();

    CGAL_precondition(CGAL::parallel(*plane_1, *plane_3));

    CGAL::FT speed_0 = std::dynamic_pointer_cast<SkelFacetData>(facet_0->getData())->getSpeed();
    CGAL::FT speed_1 = std::dynamic_pointer_cast<SkelFacetData>(facet_1->getData())->getSpeed();
    CGAL::FT speed_2 = std::dynamic_pointer_cast<SkelFacetData>(facet_2->getData())->getSpeed();
    CGAL::FT speed_3 = std::dynamic_pointer_cast<SkelFacetData>(facet_3->getData())->getSpeed();

    // The facet #1 and #3 are parallel
    // - Get the time of intersection of #1 and #3 as the intersection is at that time
    // - Intersec the three non planar planes at that time

    Vector3SPtr n_1 = KernelFactory::createVector3(plane_1);
    Vector3SPtr n_3 = KernelFactory::createVector3(plane_3);
    Point3 seed_1 = plane_1->point();
    Line3 orth_line (seed_1, *n_1);

#ifdef CGAL_SS3_DEBUG_QUAD_PLANE_INTERSECTIONS
    std::cout << "facet1 = " << facet_1->getID() << std::endl;
    std::cout << "facet3 = " << facet_3->getID() << std::endl;
    std::cout << "plane1 = " << *plane_1 << std::endl;
    std::cout << "plane3 = " << *plane_3 << std::endl;
    std::cout << "n1 = " << *n_1 << std::endl;
    std::cout << "n3 = " << *n_3 << std::endl;
    std::cout << "speed1 = " << speed_1 << std::endl;
    std::cout << "speed3 = " << speed_3 << std::endl;
    std::cout << "seed_1 = " << seed_1 << std::endl;
#endif

    CGAL::Object obj = CGAL::intersection(orth_line, *plane_3);
    if (const CGAL::Point3* seed_3_ptr = CGAL::object_cast<CGAL::Point3>(&obj)) {
#ifdef CGAL_SS3_DEBUG_QUAD_PLANE_INTERSECTIONS
        std::cout << "seed_3 = " << *seed_3_ptr << std::endl;
#endif

        // (1): p = seed_1 + t * speed_1 * n_1
        // (2): p = seed_3 + t * speed_3 * n_3
        // t = (seed_3.i - seed_1.i) / (speed_1 * n_1.i - speed_3 * n_3.i)

        CGAL::FT t;
        for (int i=0; i<3; ++i) {
            if(is_zero(n_1->operator[](i)) && is_zero(n_3->operator[](i)))
                continue;

            const CGAL::FT den = speed_1 * n_1->operator[](i) - speed_3 * n_3->operator[](i);
            if (is_zero(den)) {
                std::cerr << "Moving in the same direction, no event" << std::endl;
                return { };
            }

            t = (seed_3_ptr->operator[](i) - seed_1[i]) / den;
            break;
        }

        if (t > 0) {
            std::cerr << "Event in the past" << std::endl;
            return { };
        }

        Plane3SPtr shifted_planed_0 = KernelWrapper::offsetPlane(plane_0, t*speed_0);
        Plane3SPtr shifted_planed_1 = KernelWrapper::offsetPlane(plane_1, t*speed_1);
        Plane3SPtr shifted_planed_2 = KernelWrapper::offsetPlane(plane_2, t*speed_2);

        point = KernelWrapper::intersection(shifted_planed_0, shifted_planed_1, shifted_planed_2);
        if (!point) {
            std::cerr << "No intersection of shifted planes" << std::endl;
            return { };
        } else {
            return { point, t };
        }
    } else {
        CGAL_assertion_msg(false, "Unknown second point...?");
    }

    return { };
}

std::pair<Point3SPtr, CGAL::FT> SimpleStraightSkel::vanishesAtOnePairContiguous(FacetSPtr,
                                                                                FacetSPtr,
                                                                                FacetSPtr,
                                                                                FacetSPtr)
{
    CGAL_assertion(false);
    return {};
}

std::pair<Point3SPtr, CGAL::FT> SimpleStraightSkel::vanishesAtTwoPairs(FacetSPtr,
                                                                       FacetSPtr,
                                                                       FacetSPtr,
                                                                       FacetSPtr)
{
    // the only possibility is if the middle edge is degenerate?
    // otherwise, the bisecting planes going through edges O2 and 13 are parallel and not intersecting?

    CGAL_assertion(false);
    return {};
}




// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //





        Quadplane_parallelism parallelism;
        if (usePerturbations) {
            parallelism = PLANES_PARALLELISM_NONE;
        } else {
            // @todo get rid of all the degenerate code
            std::cerr << "You are not supposed to encounter degenerate configurations" << std::endl;
            std::exit(1); // @tmp
            parallelism = quadplane_parallelism(facetL, facetP, facetR, facetN);
        }

        if(parallelism != PLANES_PARALLELISM_NONE)
          std::cerr << "parallelism = " << parallelism << std::endl;

        if (parallelism == PLANES_PARALLELISM_NONE) {
            // generic case
            std::tie(point, offset_event) = vanishesAtGeneric(facetL, facetP, facetR, facetN, current_offset);
        } else if (parallelism == PLANES_PARALLEL_01) {
            CGAL_assertion(false);
            // std::tie(point, offset_event) = vanishesAtOnePairContiguous();
        } else if (parallelism == PLANES_PARALLEL_02) {
            // if the two main faces are parallel, everything should be parallel
            CGAL_assertion_msg(false, "This configuration shouldn't be possible with vanish events");
        } else if (parallelism == PLANES_PARALLEL_03) {
            CGAL_assertion(false);
            // std::tie(point, offset_event) = vanishesAtOnePairContiguous();
        } else if (parallelism == PLANES_PARALLEL_12) {
            CGAL_assertion(false);
            // std::tie(point, offset_event) = vanishesAtOnePairContiguous();
        } else if (parallelism == PLANES_PARALLEL_13) {
            std::tie(point, offset_event) = vanishesAtOnePairOpposite(facetL, facetP, facetR, facetN);
        } else if (parallelism == PLANES_PARALLEL_23) {
            CGAL_assertion(false);
            // std::tie(point, offset_event) = vanishesAtOnePairContiguous();
        } else if (parallelism == PLANES_PARALLEL_012 ||
                   parallelism == PLANES_PARALLEL_013 ||
                   parallelism == PLANES_PARALLEL_023 ||
                   parallelism == PLANES_PARALLEL_123) {
         // if 3 are parallel, it should be all 4 since the faces are contiguous
          CGAL_assertion_msg(false, "This configuration shouldn't be possible with vanish events");
        } else if (parallelism == PLANES_PARALLEL_02_13) {
            CGAL_assertion_msg(false, "This configuration shouldn't be possible with vanish events");
        } else if (parallelism == PLANES_PARALLEL_03_12) {
            CGAL_assertion(false);
            // std::tie(point, offset_event) = vanishesAtTwoPairs();
        } else if (parallelism == PLANES_PARALLEL_01_23) {
            CGAL_assertion(false);
            // std::tie(point, offset_event) = vanishesAtTwoPairs();
        } else if (parallelism == PLANES_PARALLEL_0123) {
            // nothing happens
        } else {
            CGAL_assertion_msg(false, "Unknown parallelism enum value");
        }




// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //






bool SimpleStraightSkel::run() {
    if (controller_) {
        controller_->wait();
        controller_->setDispPolyhedron(polyhedron_);
        controller_->setDispSkel3d(skel_result_);
    }

    DEBUG_PRINT("== Straight Skeleton 3D started ==");

    util::ConfigurationSPtr config = util::Configuration::getInstance();

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    PolyhedronSPtr polyhedron = polyhedron_->clone();

    basePlanes_.reserve(polyhedron->facets().size());
    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        CGAL_assertion(bool(facet->getPlane()));
        facet->setBasePlaneID(basePlanes_.size());
        basePlanes_.push_back(facet->getPlane());
    }

    db::_3d::OBJFile::save("results/init_pre.obj", polyhedron, false /*do not triangulate*/);

// #define CGAL_SS3_ACUTE_WEIGHTS
// #define CGAL_SS3_MERGING_WEIGHTS
// #define CGAL_SS3_PERFORMANCE_WEIGHTS

#if defined(CGAL_SS3_ACUTE_WEIGHTS) || defined(CGAL_SS3_MERGING_WEIGHTS) || defined(CGAL_SS3_PERFORMANCE_WEIGHTS)
# ifdef CGAL_SS3_ACUTE_WEIGHTS
    const CGAL::FT x_speed = 20;
    const CGAL::FT y_speed = 20;
    const CGAL::FT z_speed = 20;
    const CGAL::FT other_speed = 18.7939;
# elif defined(CGAL_SS3_MERGING_WEIGHTS)
    const CGAL::FT x_speed = 20;
    const CGAL::FT y_speed = 20;
    const CGAL::FT z_speed = 20;
    const CGAL::FT other_speed = 19.8777;
# elif defined(CGAL_SS3_PERFORMANCE_WEIGHTS)
    const CGAL::FT x_speed = 5;
    const CGAL::FT y_speed = 5;
    const CGAL::FT z_speed = 2;
    const CGAL::FT other_speed = 5;
# else
#  error
# endif

    std::size_t fi = 0;
    it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        CGAL::FT speed = other_speed;
        const auto pl = facet->plane();
        const auto normal = KernelFactory::createVector3(pl);
        std::cout << "SP X " << CGAL::scalar_product(*normal, Vector3(1,0,0)) << std::endl;
        std::cout << "SP Y " << CGAL::scalar_product(*normal, Vector3(0,1,0)) << std::endl;
        std::cout << "SP Z " << CGAL::scalar_product(*normal, Vector3(0,0,1)) << std::endl;
        if(CGAL::abs(CGAL::abs(CGAL::scalar_product(*normal, Vector3(1,0,0))) - 1) < 1e-3)
          speed = x_speed;
        if(CGAL::abs(CGAL::abs(CGAL::scalar_product(*normal, Vector3(0,1,0))) - 1) < 1e-3)
          speed = y_speed;
        if(CGAL::abs(CGAL::abs(CGAL::scalar_product(*normal, Vector3(0,0,1))) - 1) < 1e-3)
          speed = z_speed;

        SkelFacetDataSPtr data = SkelFacetData::create(facet);
        data->setSpeed(speed);
        std::cout << "speed to " << speed << std::endl;

        // for visualization, use color_ply_inputs.cpp to compile it into a single colored ply file
        std::vector<CGAL::EPICK::Point_3> points;
        std::vector<std::vector<std::size_t> > faces;

        CGAL::Cartesian_converter<CGAL::K, CGAL::EPICK> to_epick;

        std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
        while (it_v != facet->vertices().end()) {
          points.push_back(to_epick(*((*it_v++)->getPoint())));
        }

        std::vector<std::size_t> f(points.size());
        std::iota(f.begin(), f.end(), 0);
        faces.push_back(f);

        std::cout << points.size() << std::endl;
        for(auto p : points)
          std::cout << "  " << p << std::endl;

        std::cout << faces.size() << std::endl;
        std::cout << faces[0].size() << std::endl;
        for(std::size_t i : faces[0])
          std::cout << "  " << i << std::endl;

        ++fi;
    }
#endif

    DEBUG_VAL("Using " << vertex_splitter_->toString() << " to initialize polyhedron.");
    if (init(polyhedron)) {
        if (controller_) {
            controller_->wait();
        }

        CGAL::FT offset = 0.0;
        CGAL::FT offset_prev = 0.0;
        CGAL::FT offset_next = 0.0;

        db::_3d::OBJFile::save("results/init_post.obj", polyhedron, false /*do not triangulate*/);

        for(;;) {
            static int event_id = 0;

            DEBUG_PRINT(" =========== ITERATION #" << event_id << " AT OFFSET " << offset);
            DEBUG_PRINT(polyhedron->vertices().size() << " NV " << polyhedron->facets().size() << " NF");

            // std::stringstream ss_filename;
            // ss_filename << "results/" << save_path_.string() << "/face_count.txt";
            // std::ofstream out(ss_filename.str(), std::ios::app);
            // if (out) {
            //     out.precision(17);
            //     out << polyhedron->facets().size() << "\n";
            //     out.close();
            // }

            // @debug +
            DEBUG_PRINT(intersectionCache_.size() << " elements in crash cache");

            std::list<FacetSPtr>::iterator it_f_tmp = polyhedron->facets().begin();
            while (it_f_tmp != polyhedron->facets().end()) {
                FacetSPtr facet = *it_f_tmp++;
                CGAL_assertion(facet->getPlane()->a() == basePlanes_.at(facet->getBasePlaneID())->a());
                CGAL_assertion(facet->getPlane()->b() == basePlanes_.at(facet->getBasePlaneID())->b());
                CGAL_assertion(facet->getPlane()->c() == basePlanes_.at(facet->getBasePlaneID())->c());
                CGAL_assertion_code(CGAL::FT speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();)
                CGAL_assertion(facet->getPlane()->d() == basePlanes_.at(facet->getBasePlaneID())->d() - speed * offset);
            }
            // @debug -

            PQ queue;
            collectEvents(polyhedron, offset, queue);

#ifdef CGAL_SS3_DEBUG_PRINT_QUEUE
            {
                std::cout << "------------------------------" << std::endl;
                std::cout << "--- Event queue (" << event_id << ") ---" << std::endl;
                std::cout << "------------------------------" << std::endl;
                PQ duplicate_queue = queue;
                while (!duplicate_queue.empty()) {
                    AbstractEventSPtr event = duplicate_queue.top();
                    std::cout << event->toString() << std::endl;
                    duplicate_queue.pop();
                }
                std::cout << "Saves:";
                for (CGAL::FT save_offset : save_offsets_) {
                    std::cout << " " << save_offset;
                }
                std::cout << std::endl;
                std::cout << "-------------------" << std::endl;
                std::cout << "-------------------" << std::endl;
            }
#endif

            if (queue.empty()) {
                break; // we are done
            }

            // treat the next event
            AbstractEventSPtr event;
            bool doSave;
            bool simultaneousEvents;

            if (!save_offsets_.empty()) {
                CGAL::FT next_save_offset = save_offsets_.front();
                if (next_save_offset > queue.top()->getOffset()) { // save is strictly earlier
                    simultaneousEvents = false;
                    doSave = true;
                    event = SaveOffsetEvent::create(save_offsets_.front());
                } else {
                    std::tie(event, simultaneousEvents) = nextEvent(queue);
                    doSave = (next_save_offset == event->getOffset());
                }
            } else {
                std::tie(event, simultaneousEvents) = nextEvent(queue);
                doSave = false;
            }

            // @tmp +
            if (simultaneousEvents) {
              std::cerr << "Error: there should not be any simultaneous events these days" << std::endl;
              std::cerr << "Forgot to enable perturbations in the config file?" << std::endl;
              std::exit(1);
            }
            // @tmp -

            offset_next = event->getOffset();

            DEBUG_PRINT(" current offset: " << offset << "\n"
                     << " next offset: " << offset_next << " (type " << event->getType() << ")\n"
                     << " simultaneous? " << simultaneousEvents << "\n"
                     << " save? " << doSave << "\n"
                     << " in perturbed mode? " << usingTemporaryPerturbedMode_
                     << " @ " << simultaneousOffset_);
#ifdef CGAL_SS3_RUN_TIMERS
            std::cout << "current elapsed time: " << timer.time(); << std::endl;
#endif

#ifdef CGAL_SS3_USE_AUTOREF_FOR_ALL_EVENTS
            std::cout << "============== USING AUTOREF APPROACH =============" << std::endl;
            simultaneousEvents = false;
#endif

            // switch ON perturbation if needed
            if (simultaneousEvents) {
                if (doSave) {
                    std::cout << "Save event at simultaneity" << std::endl;
                    return handleSaveEventAtSimultaneity(polyhedron, offset, offset_next);
                }

#ifndef CGAL_SS3_USE_AUTOREF_FOR_SIMULTANEOUS_EVENTS
                if (offset_next == 0) {
                    std::cerr << "NYI: simultaneous event at offset=0" << std::endl;
                    std::exit(1);
                }

                if (usingTemporaryPerturbedMode_) {
                    std::cerr << "Error: met a simultaneous event in perturbed mode!" << std::endl;
                    std::exit(1); // nuke it
                } else {
                    std::cout << "simultaneous events ahead, at offset: " << offset_next << std::endl;
                    std::tie(polyhedron, offset) = enablePerturbedMode(polyhedron, offset, offset_next);
                    CGAL_assertion(offset > offset_next);
                    continue; // recompute events
                }
#endif
            }

            // switch OFF perturbation if needed
            if (usingTemporaryPerturbedMode_) {
                CGAL::FT delta = 1e-5; // @fixme don't hardcode this value
                // if the next event is far and we are in perturbed mode, disable perturbation
                if (simultaneousOffset_ - delta > offset_next) { // offsets are negative
                    std::tie(polyhedron, offset) = disablePerturbedMode(polyhedron, offset, offset_next);
                    CGAL_assertion(offset > offset_next);
                    continue; // recompute the polyhedron and its events
                }
            }

            DEBUG_PRINT("\n-----------------------------------------------------");
            DEBUG_PRINT("-- Event #" << event_id << " " << event->toString() << " --");
            DEBUG_PRINT("-----------------------------------------------------\n");

            if (controller_) {
                controller_->wait();
            }

            Point3SPtr p_box_min = PolyhedronTransformation::boundingBoxMin(polyhedron);
            Point3SPtr p_box_max = PolyhedronTransformation::boundingBoxMax(polyhedron);

#ifdef CGAL_SS3_USE_AUTOREF_FOR_ALL_EVENTS
            // @fixme delta should depend on the next save and next constant events
            // @fixme delta should be small enough so that there are no other events
            CGAL::FT nudged_offset;
            std::tie(polyhedron, nudged_offset) = handleEventWithAutoref(event, offset, polyhedron);

            offset_prev = offset;
            offset = nudged_offset;
#else
# ifdef CGAL_SS3_USE_AUTOREF_FOR_SIMULTANEOUS_EVENTS
            if (simultaneousEvents) {
                if (doSave) {
                    return handleSaveEventAtSimultaneity(polyhedron, offset, offset_next);
                } else {
                    CGAL::FT nudged_offset;
                    std::tie(polyhedron, nudged_offset) = handleEventWithAutoref(event, offset, polyhedron);

                    offset_prev = offset;
                    offset = nudged_offset;
                }
            } else
# endif
            {
                offset_prev = offset;
                offset = event->getOffset();
                const CGAL::FT shift = offset - offset_prev;
                CGAL_warning(!usingTemporaryPerturbedMode_ || !is_zero(shift));
                const bool recompute_positions = (shift != 0);
                polyhedron = PolyhedronTransformation::shiftFacets(polyhedron, shift, recompute_positions);

                // below will have degeneracies since we haven't treated the event yet
                // db::_3d::OBJFile::save("results/shift_" + std::to_string(event_id) + ".obj",
                //                        polyhedron,
                //                        false /*do not triangulate*/);

                if (event->getType() == AbstractEvent::CONST_OFFSET_EVENT) {
#ifndef CGAL_SS3_NO_SKELETON_DS
                    event->setPolyhedronResult(polyhedron);
#endif
                    skel_result_->addEvent(event);

                    bool screenshot_on_const_offset_event =
                            util::Configuration::getInstance()->getBool(
                            "algo_3d_SimpleStraightSkel", "screenshot_on_const_offset_event");
                    if (controller_ && screenshot_on_const_offset_event) {
                        controller_->screenshot();
                    }
                } else if (event->getType() == AbstractEvent::EDGE_EVENT) {
                    handleEdgeEvent(std::dynamic_pointer_cast<EdgeEvent>(event), polyhedron);
                } else if (event->getType() == AbstractEvent::EDGE_MERGE_EVENT) {
                    handleEdgeMergeEvent(std::dynamic_pointer_cast<EdgeMergeEvent>(event), polyhedron);
                } else if (event->getType() == AbstractEvent::TRIANGLE_EVENT) {
                    handleTriangleEvent(std::dynamic_pointer_cast<TriangleEvent>(event), polyhedron);
                } else if (event->getType() == AbstractEvent::DBL_EDGE_MERGE_EVENT) {
                    handleDblEdgeMergeEvent(std::dynamic_pointer_cast<DblEdgeMergeEvent>(event), polyhedron);
                } else if (event->getType() == AbstractEvent::DBL_TRIANGLE_EVENT) {
                    handleDblTriangleEvent(std::dynamic_pointer_cast<DblTriangleEvent>(event), polyhedron);
                } else if (event->getType() == AbstractEvent::TETRAHEDRON_EVENT) {
                    handleTetrahedronEvent(std::dynamic_pointer_cast<TetrahedronEvent>(event), polyhedron);
                } else if (event->getType() == AbstractEvent::VERTEX_EVENT) {
                    handleVertexEvent(std::dynamic_pointer_cast<VertexEvent>(event), polyhedron);
                } else if (event->getType() == AbstractEvent::FLIP_VERTEX_EVENT) {
                    handleFlipVertexEvent(std::dynamic_pointer_cast<FlipVertexEvent>(event), polyhedron);
                } else if (event->getType() == AbstractEvent::SURFACE_EVENT) {
                    handleSurfaceEvent(std::dynamic_pointer_cast<SurfaceEvent>(event), polyhedron);
                } else if (event->getType() == AbstractEvent::POLYHEDRON_SPLIT_EVENT) {
                    handlePolyhedronSplitEvent(std::dynamic_pointer_cast<PolyhedronSplitEvent>(event), polyhedron);
                } else if (event->getType() == AbstractEvent::SPLIT_MERGE_EVENT) {
                    handleSplitMergeEvent(std::dynamic_pointer_cast<SplitMergeEvent>(event), polyhedron);
                } else if (event->getType() == AbstractEvent::EDGE_SPLIT_EVENT) {
                    handleEdgeSplitEvent(std::dynamic_pointer_cast<EdgeSplitEvent>(event), polyhedron);
                } else if (event->getType() == AbstractEvent::PIERCE_EVENT) {
                    handlePierceEvent(std::dynamic_pointer_cast<PierceEvent>(event), polyhedron);
                }
            }
#endif // CGAL_SS3_USE_AUTOREF_FOR_ALL_EVENTS

            DEBUG_PRINT("-- Finished handling Event --");

            db::_3d::OBJFile::save("results/iter_" + std::to_string(event_id) + ".obj", polyhedron, false /*do triangulate*/);
            db::_3d::OBJFile::save("results/iter_" + std::to_string(event_id) + "_triangulated.obj", polyhedron);

            // std::cout << "-- Degen count --" << std::endl;
            // std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
            // while (it_f != polyhedron->facets().end()) {
            //     FacetSPtr facet = *it_f++;
            //     auto it_v = facet->vertices().begin();
            //     Point3SPtr p0 = (*(it_v++))->getPoint();
            //     Point3SPtr p1 = (*(it_v++))->getPoint();
            //     Point3SPtr p2 = (*it_v)->getPoint();
            //     if (CGAL::collinear(*p0, *p1, *p2)) {
            //         std::cout << *p0 << " " << *p1 << " " << *p2 << " is degen" << std::endl;
            //     }
            // }

            CGAL_assertion(polyhedron->isConsistent());
#ifndef CGAL_SS3_NO_SKELETON_DS
            CGAL_assertion(skel_result_->isConsistent());
#endif
            CGAL_assertion(p_box_min && p_box_max);
            CGAL_assertion(PolyhedronTransformation::isInsideBox(polyhedron, p_box_min, p_box_max));

            // this is tempting, but the mesh is usually not in a nice state here
            // CGAL_assertion(!SelfIntersection::hasSelfIntersectingSurface(polyhedron));

            if (controller_) {
                controller_->wait();
            }

            if (doSave) {
                savePolyhedron(polyhedron, offset,
                               true /*triangulate*/,
                               true /*convert to double*/,
                               false /*attempt untilting*/);

                save_offsets_.pop_front();

                if (save_offsets_.empty()) {
                    if (config->isLoaded()) {
                        if ((config->contains("main", "stop_after_last_save_event") &&
                             config->getBool("main", "stop_after_last_save_event"))) {
                            break;
                        }
                    }
                }
            }

            // Can't perturb to treat simultaneous events AND get a meaningful "SAVE" result
            if (simultaneousEvents && event->getType() == AbstractEvent::SAVE_OFFSET_EVENT) {
                std::cout << "WARNING: simultaneous event @ save time, not proceeding farther" << std::endl;
                break;
            }

            ++event_id;
        }


        DEBUG_PRINT("== Straight Skeleton 3D finished ==");

#ifdef CGAL_SS3_RUN_TIMERS
        timer.stop();
        skel_result_->appendDescription("time=" + timer.time() + "; ");
#endif

        //skel_result_->appendDescription("controller=" +
        //        util::StringFactory::fromBoolean(controller_) + "; ");
        std::cout << skel_result_->toString() << std::endl;
    } else {
        DEBUG_PRINT("Error: Failed to initialize");
        return false;
    }

    return true;
}



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



void SimpleStraightSkel::registerAllowedInteractions(PolyhedronSPtr polyhedron,
                                                     double maxInteractionDistance)
{
    std::cout << "MEGA HACK" << std::endl;

    allowedInteractions_.resize(polyhedron->facets().size());

    // Create bounding boxes for each facet
    using Box = CGAL::Box_intersection_d::Box_with_info_d<double, 3, FacetSPtr>;

    std::vector<Box> boxes;
    boxes.reserve(polyhedron->facets().size());

    // triangulated faces... should be a mesh and correspondence, really...
    std::vector<Point3> points;
    std::map<Point3, std::size_t> pids;
    std::vector<std::vector<std::vector<std::size_t> > > faces_triangles;
#ifdef CGAL_SS3_USE_BOX_D_HACK_WITH_FACE_DISTANCE_FILTERING
    faces_triangles.reserve(polyhedron->facets().size());
#endif

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;

        std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
        CGAL::Bbox_3 bbox = (*it_v++)->getPoint()->bbox(); // don't initialize with 0!!
        while (it_v != facet->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            bbox += vertex->getPoint()->bbox();
        }

        bbox = CGAL::Bbox_3(bbox.xmin() - maxInteractionDistance,
                            bbox.ymin() - maxInteractionDistance,
                            bbox.zmin() - maxInteractionDistance,
                            bbox.xmax() + maxInteractionDistance,
                            bbox.ymax() + maxInteractionDistance,
                            bbox.zmax() + maxInteractionDistance);
        boxes.emplace_back(bbox, facet);

#ifdef CGAL_SS3_USE_BOX_D_HACK_WITH_FACE_DISTANCE_FILTERING
        std::vector<std::vector<std::size_t> > triangles;
        triangulate_facet_with_CDT2(facet, CDT2_Filtering::ODD_EVEN, pids, points, triangles);
        faces_triangles.emplace_back(std::move(triangles));
#endif
    }

    // Define the callback for box intersection
    auto callback = [this, maxInteractionDistance, &points, &faces_triangles](const Box& a, const Box& b) {
        FacetSPtr fa = a.info(), fb = b.info();
        if (fa->getID() == fb->getID()) {
            return;
        }

#ifdef CGAL_SS3_USE_BOX_D_HACK_WITH_FACE_DISTANCE_FILTERING
        // real brute force...
        CGAL::FT min_sq_d = (std::numeric_limits<double>::max)();
        for (const auto& tri_a : faces_triangles[fa->getID()]) {
            for (const auto& tri_b : faces_triangles[fb->getID()]) {
                Triangle3 t_a(points[tri_a[0]], points[tri_a[1]], points[tri_a[2]]);
                Triangle3 t_b(points[tri_b[0]], points[tri_b[1]], points[tri_b[2]]);
                CGAL::FT sq_d = CGAL::squared_distance(t_a, t_b);
                if (sq_d < min_sq_d) {
                    min_sq_d = sq_d;
                }
            }
        }

        const bool is_close_enough = (min_sq_d < CGAL::square(maxInteractionDistance));
        if (!is_close_enough) {
            return;
        }
#endif

        allowedInteractions_[fa->getID()].push_back(fb->getID());
        allowedInteractions_[fb->getID()].push_back(fa->getID());
    };

    // Perform box intersection
    CGAL::box_self_intersection_d(boxes.begin(), boxes.end(), callback);
}

std::vector<std::vector<std::size_t> > allowedInteractions_; // CGAL_SS3_USE_BOX_D_HACK


// ----
// Macros to hack a reduction of eligible face-face interactions based on some max distance

// #define CGAL_SS3_USE_BOX_D_HACK
// #define CGAL_SS3_USE_BOX_D_HACK_WITH_FACE_DISTANCE_FILTERING
#ifdef CGAL_SS3_USE_BOX_D_HACK_WITH_FACE_DISTANCE_FILTERING
# ifndef CGAL_SS3_USE_BOX_D_HACK
#  define CGAL_SS3_USE_BOX_D_HACK
# endif
#endif


#ifdef CGAL_SS3_USE_BOX_D_HACK
        std::cout << "Compute allowed interactions between faces" << std::endl;
        registerAllowedInteractions(polyhedron, 1500);
#endif


#ifdef CGAL_SS3_USE_BOX_D_HACK
    std::array<FacetSPtr, 4> facets = { facet_0, facet_1, facet_2, facet_3 };
    for (std::size_t i=0; i<4; ++i) {
        for (std::size_t j=i+1; j<4; ++j) {
            const auto& ai_i = allowedInteractions_[facets[i]->getID()];
            if (std::find(std::cbegin(ai_i), std::cend(ai_i), facets[j]->getID()) == std::cend(ai_i)) {
                return { };
            }
        }
    }
#endif



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


#ifdef CGAL_SS3_FILTER_PIERCE_EVENTS_AT_POP_TIME
    // Check if it is a valid pierce event: the piercing point must be in the face
    bool reject_event = false;

# if 1 // not sure which one is better
    // with ray shooting
    reject_event = !(SelfIntersection::isInsideWithRayShooting(node->getPoint(), facet_offset));
# else
    using Itag = CGAL::No_constraint_intersection_tag;
    using PK = CGAL::Projection_traits_3<CGAL::K>;
    using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>; // @todo not needed
    using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
    using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
    using PTDS = CGAL::Triangulation_data_structure_2<PVb,PFb>;
    using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
    using PCDT_VH = PCDT::Vertex_handle;
    using PCDT_FH = PCDT::Face_handle;

    Vector3SPtr n = KernelFactory::createVector3(facet_offset->plane());
    CGAL_assertion(*n != CGAL::NULL_VECTOR);
    PK traits(*n);
    PCDT pcdt(traits);

    std::map<VertexSPtr, PCDT_VH> face_vhs;
    std::list<VertexSPtr>::iterator it_v = facet_offset->vertices().begin();
    while (it_v != facet_offset->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        auto res = face_vhs.emplace(vertex, PCDT_VH());
        if(res.second) // first time seeing this point
        {
            PCDT_VH vh = pcdt.insert(*(vertex->getPoint()));
            res.first->second = vh;
        }
    }

    std::list<EdgeSPtr>::iterator it_e = facet_offset->edges().begin();
    while (it_e != facet_offset->edges().end()) {
        EdgeSPtr edge = *it_e++;
        VertexSPtr v0 = edge->src(facet_offset);
        VertexSPtr v1 = edge->dst(facet_offset);

        if(*(v0->getPoint()) == *(v1->getPoint()))
        {
            std::cerr << "W: encountered degenerate edge @ " << *(v0->getPoint()) << std::endl;
        }
        else
        {
            PCDT_VH vh0 = face_vhs.at(v0);
            PCDT_VH vh1 = face_vhs.at(v1);

            try
            {
                pcdt.insert_constraint(vh0, vh1);
            }
            catch(const typename PCDT::Intersection_of_constraints_exception&)
            {
                std::cerr << "Error: Intersection of constraints" << std::endl;
                DEBUG_VAR(facet_offset->toString());
                CGAL_assertion_msg(false, "Intersections in CDT2 not allowed");
                break;
            }
        }
    }

    std::unordered_map<PCDT_FH, bool> in_domain_map;
    boost::associative_property_map<std::unordered_map<PCDT_FH, bool> > in_domain(in_domain_map);
    CGAL::mark_domain_in_triangulation(pcdt, in_domain);

    int li;
    PCDT::Locate_type lt;
    PCDT_FH fh = pcdt.locate(*(node->getPoint()), lt, li);

    if (lt == PCDT::VERTEX || lt == PCDT::EDGE) {
        CGAL_assertion(false); // @todo handle this by looping over incident faces
    } else if (lt == PCDT::FACE) {
        reject_event = !(get(in_domain, fh));
    } else { // outside of convex hull
        reject_event = true;
    }
# endif // ray shooting or CDT
    if (reject_event) {
        std::cout << "Pierce Event rejected" << std::endl;
        event->setPolyhedronResult(polyhedron);
        skel_result_->addEvent(event);
        return;
    }
#endif // CGAL_SS3_FILTER_PIERCE_EVENTS_AT_POP_TIME


// #define CGAL_SS3_FILTER_PIERCE_EVENTS_AT_POP_TIME
# ifdef CGAL_SS3_FILTER_PIERCE_EVENTS_AT_POP_TIME
                    // Now, we might naively wish to filter using bisectors like in 2D SLS code,
                    // but unlike a segment, a face in the 3D SS code has no reason to be convex,
                    // which changes everything and can result in false positives.
                    //
                    // The bisector filter in 2D is equivalent to checking if the point is on the offset
                    // face. We could check this here, but determining what is the offset face at this
                    // point (i.e., while searching for events) is rough: plenty of other events
                    // might modify the face before this particular pierce event appears, and so
                    // we can't just do shift(facet) because the result might be a self-intersecting
                    // polygon with holes.
                    //
                    // Instead, we do not filter here, but simply put it in the queue. When the event
                    // will be popped, then we know it's the next event globally and nothing else
                    // can mess up the face, and we can do the in-test check then.
                    //
                    // See HandlePierceEvent()
                    //
                    // One exception: if offset_event is 0 (meaning, pierce event at the current time),
                    // then we can and should do the filtering at this particular point for two reasons:
                    // - if we put it in the queue and filter at pop, we might get an endless loop
                    // - it is safe: there can't be topological changes like described above
                    //
                    // Note that it'll avoid the infinite loop as such:
                    // - pierce event at t = t_0 is put in the queue
                    // - pierce event is top of the queue
                    // - polyhedron shifts
                    // - the pierce is rejected
                    // - pierce event is detected (again) and rejected (again)
                    // - pierce event is not put in the queue
                    // @todo needless computations on the penultimate step
                    if (offset_event == 0)
                    {
                        if (!SelfIntersection::isInsideWithRayShooting(point, facet)) {
                            std::cout << "Filtered T=0 pierce event" << std::endl;
                            continue;
                        }
                    }
# else

#  if 0
                    // @fixme actually bugged? (construction of the offset of the facet was broken when I put the other function...)
                    if (!IsLineInFacet(facet_offset, arc->line())) {
                        // std::cout << "IsLineInFacet rejects" << std::endl;
                        continue;
                    }
#  else
                    // Note that this result could be meaningless if the offset face
                    // is not a simple polygon. However, if it's not simple, then some event
                    // has happened before the pierce, and the pierce event - if whitelisted -
                    // would be checked again later, thus it's safe to call.
                    if (!SelfIntersection::isInsideWithRayShooting(point, facet_offset)) {
                        // std::cout << "isInsideWithRayShooting rejects" << std::endl;
                        continue;
                    }
#  endif
# endif






// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //







#if 0 // individual events of a given type
EdgeEventSPtr SimpleStraightSkel::nextEdgeEvent(PolyhedronSPtr polyhedron,
                                                const CGAL::FT current_offset,
                                                CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Edge Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    EdgeEventSPtr result = EdgeEventSPtr();

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        VertexSPtr vertex_src = edge->getVertexSrc();
        VertexSPtr vertex_dst = edge->getVertexDst();
        if (vertex_src->getPoint() == vertex_dst->getPoint()) {
            continue;
        }

        if (isLocked(edge)) { // one of the vertices cannot be rebuilt from its incident faces
            continue;
        }

        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        if (isTriangle(facet_l, edge) || isTriangle(facet_r, edge)) {
            // triangle event
            continue;
        }

        Point3SPtr point = Point3SPtr();
        CGAL::FT offset_event;
        std::tie(point, offset_event) = vanishesAt(edge, current_offset);
        if (!point) {
            continue;
        }

        std::cout << "Tentative edge event @ " << edge->toString() << std::endl;
        std::cout << "Would give: " << *point << " (at t=" << offset_event << ")" << std::endl;

        FacetSPtr facet_src = getFacetSrc(edge);
        FacetSPtr facet_dst = getFacetDst(edge);

        // This does not work when there is more than one edge between both facets.
        // EdgeSPtr edge_2 = facet_src->findEdge(facet_dst);
        std::list<EdgeSPtr> edges_2 = facet_src->findEdges(facet_dst);

        bool split_event = false;
        std::list<EdgeSPtr>::iterator it_e2 = edges_2.begin();
        while (it_e2 != edges_2.end()) {
            EdgeSPtr edge_2 = *it_e2++;

            if (isLocked(edge_2)) {
                continue;
            }

            // @todo exit as soon as there is a single rejection (to be done once the old code is safe to remove)
#if defined(CGAL_SS3_OLD_CODE_BOUND_CHECKS) || defined(CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS)
            bool split_event_current_1 = true;
            bool split_event_current_2 = true;
            bool split_event_current_3 = true;

            SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(edge_2->getData());
            Vector3SPtr normal_2 = KernelFactory::createVector3(data_2->getSheet()->getPlane());
            Line3SPtr line_normal_2 = KernelFactory::createLine3(point, normal_2);
            if (KernelWrapper::orientation(line(edge_2), line_normal_2) < 0) {
                // out of bounded area
                split_event_current_1 = false;
            }
            SkelVertexDataSPtr data_2_src = std::dynamic_pointer_cast<SkelVertexData>(
                edge_2->getVertexSrc()->getData());
            ArcSPtr arc_2_src = data_2_src->getArc();
            if (KernelWrapper::orientation(arc_2_src->line(), line_normal_2) > 0) {
                // out of bounded area
                split_event_current_2 = false;
            }
            SkelVertexDataSPtr data_2_dst = std::dynamic_pointer_cast<SkelVertexData>(
                edge_2->getVertexDst()->getData());
            ArcSPtr arc_2_dst = data_2_dst->getArc();
            if (KernelWrapper::orientation(arc_2_dst->line(), line_normal_2) < 0) {
                // out of bounded area
                split_event_current_3 = false;
            }

            const bool split_event_current = (split_event_current_1 &&
                                              split_event_current_2 &&
                                              split_event_current_3);
#endif

#if !defined(CGAL_SS3_OLD_CODE_BOUND_CHECKS) || defined(CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS)
            bool split_event_current_1_b = true;
            bool split_event_current_2_b = true;
            bool split_event_current_3_b = true;

            FacetSPtr facet_l2 = edge_2->getFacetL();
            FacetSPtr facet_r2 = edge_2->getFacetR();
            FacetSPtr facet_2_src = getFacetSrc(edge_2);
            FacetSPtr facet_2_dst = getFacetDst(edge_2);

            Plane3SPtr plane_l2 = facet_l2->plane();
            Plane3SPtr plane_r2 = facet_r2->plane();
            CGAL::FT speed_l2 = std::dynamic_pointer_cast<SkelFacetData>(facet_l2->getData())->getSpeed();
            CGAL::FT speed_r2 = std::dynamic_pointer_cast<SkelFacetData>(facet_r2->getData())->getSpeed();

            CGAL::FT l2a = plane_l2->a();
            CGAL::FT l2b = plane_l2->b();
            CGAL::FT l2c = plane_l2->c();
            CGAL::FT l2d = plane_l2->d();
            CGAL::FT r2a = plane_r2->a();
            CGAL::FT r2b = plane_r2->b();
            CGAL::FT r2c = plane_r2->c();
            CGAL::FT r2d = plane_r2->d();

            CGAL::FT lt2 = (l2a * point->x() + l2b * point->y() + l2c * point->z() + l2d) / speed_l2;
            CGAL::FT rt2 = (r2a * point->x() + r2b * point->y() + r2c * point->z() + r2d) / speed_r2;
            CGAL_assertion(lt2 == rt2);

            if ((lt2 > 0) || (rt2 > 0)) {
                split_event_current_1_b = false;
            }

            if (!check_bisector(edge_2, facet_r2, rt2, facet_2_src, point)) {
                split_event_current_2_b = false;
            }

            if (!check_bisector(edge_2, facet_l2, lt2, facet_2_dst, point)) {
                split_event_current_3_b = false;
            }

            const bool split_event_current_b = (split_event_current_1_b &&
                                                split_event_current_2_b &&
                                                split_event_current_3_b);
#endif

#ifdef CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS
            CGAL_assertion(split_event_current_1 == split_event_current_1_b);
            CGAL_assertion(split_event_current_2 == split_event_current_2_b);
            CGAL_assertion(split_event_current_3 == split_event_current_3_b);

            CGAL_assertion(split_event_current == split_event_current_b);
#endif

            if (split_event_current_b) {
                split_event = true;
                break;
            }
        }
        if (split_event) {
            continue;
        }
        // edge merge event
        EdgeSPtr edge_prev = edge->prev(facet_l);
        EdgeSPtr edge_next = edge->next(facet_l)->next(facet_l);
        if (edge_prev->hasSameFacets(edge_next)) {
            continue;
        }
        edge_prev = edge->prev(facet_l)->prev(facet_l);
        edge_next = edge->next(facet_l);
        if (edge_prev->hasSameFacets(edge_next)) {
            continue;
        }
        edge_prev = edge->prev(facet_r);
        edge_next = edge->next(facet_r)->next(facet_r);
        if (edge_prev->hasSameFacets(edge_next)) {
            continue;
        }
        edge_prev = edge->prev(facet_r)->prev(facet_r);
        edge_next = edge->next(facet_r);
        if (edge_prev->hasSameFacets(edge_next)) {
            continue;
        }
        // @todo could move up this check
        if (offset_event < 0 && offset_event > current_offset_to_nearest_event)
        {
            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = EdgeEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setEdge(edge);

#ifndef CGAL_SS3_NO_SKELETON_DS
            SkelVertexDataSPtr data_src = std::dynamic_pointer_cast<SkelVertexData>(
                    edge->getVertexSrc()->getData());
            SkelVertexDataSPtr data_dst = std::dynamic_pointer_cast<SkelVertexData>(
                    edge->getVertexDst()->getData());
            node->addArc(data_src->getArc());
            node->addArc(data_dst->getArc());
            SkelEdgeDataSPtr data_edge = std::dynamic_pointer_cast<SkelEdgeData>(
                    edge->getData());
            node->addSheet(data_edge->getSheet());
#endif

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }

    return result;
}

EdgeMergeEventSPtr SimpleStraightSkel::nextEdgeMergeEvent(PolyhedronSPtr polyhedron,
                                                          const CGAL::FT current_offset,
                                                          CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Edge Merge Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    EdgeMergeEventSPtr result = EdgeMergeEventSPtr();

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;

        VertexSPtr vertex_src = edge->getVertexSrc();
        VertexSPtr vertex_dst = edge->getVertexDst();
        if (vertex_src->getPoint() == vertex_dst->getPoint()) {
            continue;
        }

        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        if (isTriangle(facet_l, edge) || isTriangle(facet_r, edge)) {
            // triangle event
            continue;
        }

        FacetSPtr facet_other = edge->getFacetL();
        EdgeSPtr edge_next = edge->next(facet_other);
        facet_other = edge_next->other(facet_other);
        edge_next = edge_next->prev(facet_other);
        facet_other = edge_next->other(facet_other);
        edge_next = edge_next->next(facet_other);
        facet_other = edge_next->other(facet_other);
        edge_next = edge_next->prev(facet_other);
        if (edge_next == edge) {
            // dbl edge merge event
            continue;
        }

        facet_other = edge->getFacetR();
        edge_next = edge->prev(facet_other);
        facet_other = edge_next->other(facet_other);
        edge_next = edge_next->next(facet_other);
        facet_other = edge_next->other(facet_other);
        edge_next = edge_next->prev(facet_other);
        facet_other = edge_next->other(facet_other);
        edge_next = edge_next->next(facet_other);
        if (edge_next == edge) {
            // dbl edge merge event
            continue;
        }

        FacetSPtr facet = FacetSPtr();
        EdgeSPtr edge_1 = EdgeSPtr();
        EdgeSPtr edge_2 = EdgeSPtr();
        EdgeSPtr edge_prev = edge->prev(facet_l);
        edge_next = edge->next(facet_l)->next(facet_l);
        if (edge_prev->hasSameFacets(edge_next) && edge_prev != edge_next) {
            facet = facet_l;
            edge_1 = edge_prev;
            edge_2 = edge_next;
        }
        // @todo do we need to test these other combinations if above matched?
        edge_prev = edge->prev(facet_l)->prev(facet_l);
        edge_next = edge->next(facet_l);
        if (edge_prev->hasSameFacets(edge_next) && edge_prev != edge_next) {
            facet = facet_l;
            edge_1 = edge_prev;
            edge_2 = edge_next;
        }
        edge_prev = edge->prev(facet_r);
        edge_next = edge->next(facet_r)->next(facet_r);
        if (edge_prev->hasSameFacets(edge_next) && edge_prev != edge_next) {
            facet = facet_r;
            edge_1 = edge_prev;
            edge_2 = edge_next;
        }
        edge_prev = edge->prev(facet_r)->prev(facet_r);
        edge_next = edge->next(facet_r);
        if (edge_prev->hasSameFacets(edge_next) && edge_prev != edge_next) {
            facet = facet_r;
            edge_1 = edge_prev;
            edge_2 = edge_next;
        }
        if (!(facet && edge_1 && edge_2)) {
            continue;
        }

        Point3SPtr point = Point3SPtr();
        CGAL::FT offset_event;
        std::tie(point, offset_event) = vanishesAt(edge, current_offset);
        if (!point) {
            continue;
        }
        if (offset_event < 0 && offset_event > current_offset_to_nearest_event) {
            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = EdgeMergeEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setFacet(facet);
            result->setEdge1(edge_1);
            result->setEdge2(edge_2);

#ifndef CGAL_SS3_NO_SKELETON_DS
            EdgeSPtr edge_toremove_1 = edge_1->next(facet);
            EdgeSPtr edge_toremove_2 = edge_toremove_1->next(facet);
            SkelVertexDataSPtr data_vertex = std::dynamic_pointer_cast<SkelVertexData>(
                    edge_toremove_1->src(facet)->getData());
            node->addArc(data_vertex->getArc());
            data_vertex = std::dynamic_pointer_cast<SkelVertexData>(
                    edge_toremove_1->dst(facet)->getData());
            node->addArc(data_vertex->getArc());
            data_vertex = std::dynamic_pointer_cast<SkelVertexData>(
                    edge_toremove_2->dst(facet)->getData());
            node->addArc(data_vertex->getArc());
            SkelEdgeDataSPtr data_edge = std::dynamic_pointer_cast<SkelEdgeData>(
                    edge_toremove_1->getData());
            node->addSheet(data_edge->getSheet());
            data_edge = std::dynamic_pointer_cast<SkelEdgeData>(
                    edge_toremove_2->getData());
            node->addSheet(data_edge->getSheet());
#endif

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }

    if(result)
    {
      std::cout << "Best Edge Merge Event with time increment of: " << result->getOffset() << std::endl;
    }

    return result;
}

TriangleEventSPtr SimpleStraightSkel::nextTriangleEvent(PolyhedronSPtr polyhedron,
                                                        const CGAL::FT current_offset,
                                                        CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Triangle Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    TriangleEventSPtr result = TriangleEventSPtr();

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;

        if (edge->getVertexSrc()->getPoint() == edge->getVertexDst()->getPoint()) {
            continue;
        }

        if (isTetrahedron(edge)) {
            // tetrahedron event
            continue;
        }

        FacetSPtr facet;
        if (isTriangle(edge->getFacetL(), edge)) {
            facet = edge->getFacetL();
        } else if (isTriangle(edge->getFacetR(), edge)) {
            facet = edge->getFacetR();
        } else {
            continue;
        }

        bool dbl_triangle_event = false;
        EdgeSPtr edge_tmp = edge;
        for (unsigned int i = 0; i < 3; i++) {
            FacetSPtr facet_tmp_l = edge_tmp->getFacetL();
            FacetSPtr facet_tmp_r = edge_tmp->getFacetR();
            if (facet_tmp_l && facet_tmp_r) {
                if (isTriangle(facet_tmp_l, edge_tmp) &&
                        isTriangle(facet_tmp_r, edge_tmp)) {
                    dbl_triangle_event = true;
                    break;
                }
            }
            edge_tmp = edge_tmp->next(facet);
        }
        if (dbl_triangle_event) {
            continue;
        }

        Point3SPtr point = Point3SPtr();
        CGAL::FT offset_event;
        std::tie(point, offset_event) = vanishesAt(edge, current_offset);
        if (!point) {
            continue;
        }

        if ((KernelWrapper::side(edge->getFacetL()->plane(), point) > 0) ||
                KernelWrapper::side(edge->getFacetR()->plane(), point) > 0) {
            // triangle may not be a hole
            // after pierce event
            continue;
        }

        if (offset_event < 0 && offset_event > current_offset_to_nearest_event) {
            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = TriangleEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setFacet(facet);
            result->setEdgeBegin(edge);

#ifndef CGAL_SS3_NO_SKELETON_DS
            VertexSPtr vertices[3];
            result->getVertices(vertices);
            for (unsigned int i = 0; i < 3; i++) {
                SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(
                        vertices[i]->getData());
                ArcSPtr arc = data->getArc();
                node->addArc(arc);
            }
            EdgeSPtr edges[3];
            result->getEdges(edges);
            for (unsigned int i = 0; i < 3; i++) {
                SkelEdgeDataSPtr data = std::dynamic_pointer_cast<SkelEdgeData>(
                        edges[i]->getData());
                SheetSPtr sheet = data->getSheet();
                node->addSheet(sheet);
            }
#endif

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }
    return result;
}

DblEdgeMergeEventSPtr SimpleStraightSkel::nextDblEdgeMergeEvent(PolyhedronSPtr polyhedron,
                                                                const CGAL::FT current_offset,
                                                                CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Dbl Edge Merge Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    DblEdgeMergeEventSPtr result = DblEdgeMergeEventSPtr();

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;

        if (isLocked(edge)) {
            continue;
        }

        if (!isReflex(edge)) {
            continue;
        }

        bool is_dbl_edge_merge_event = false;
        FacetSPtr facet_1;
        EdgeSPtr edge_11;
        EdgeSPtr edge_12;
        FacetSPtr facet_2;
        EdgeSPtr edge_21;
        EdgeSPtr edge_22;
        FacetSPtr facet_other = edge->getFacetL();
        EdgeSPtr edge_next = edge->next(facet_other);
        facet_other = edge_next->other(facet_other);
        edge_next = edge_next->prev(facet_other);
        facet_other = edge_next->other(facet_other);
        edge_next = edge_next->next(facet_other);
        facet_other = edge_next->other(facet_other);
        edge_next = edge_next->prev(facet_other);
        if (edge_next == edge) {
            is_dbl_edge_merge_event = true;
            facet_1 = edge->getFacetL();
            edge_11 = edge->prev(facet_1);
            edge_12 = edge->next(facet_1)->next(facet_1);
            facet_2 = edge->getFacetR();
            edge_21 = edge->prev(facet_2);
            edge_22 = edge->next(facet_2)->next(facet_2);
        }

        facet_other = edge->getFacetR();
        edge_next = edge->prev(facet_other);
        facet_other = edge_next->other(facet_other);
        edge_next = edge_next->next(facet_other);
        facet_other = edge_next->other(facet_other);
        edge_next = edge_next->prev(facet_other);
        facet_other = edge_next->other(facet_other);
        edge_next = edge_next->next(facet_other);
        if (edge_next == edge) {
            is_dbl_edge_merge_event = true;
            facet_1 = edge->getFacetR();
            edge_11 = edge->prev(facet_1)->prev(facet_1);
            edge_12 = edge->next(facet_1);
            facet_2 = edge->getFacetL();
            edge_21 = edge->prev(facet_2)->prev(facet_2);
            edge_22 = edge->next(facet_2);
        }

        if (edge_11 == edge_12 || edge_21 == edge_22) {
            // double triangle event
            continue;
        }

        if (!is_dbl_edge_merge_event) {
            continue;
        }

        Point3SPtr point = Point3SPtr();
        CGAL::FT offset_event;
        std::tie(point, offset_event) = vanishesAt(edge, current_offset);
        if (!point) {
            continue;
        }
        if (offset_event < 0 && offset_event > current_offset_to_nearest_event) {
            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = DblEdgeMergeEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setFacet1(facet_1);
            result->setEdge11(edge_11);
            result->setEdge12(edge_12);
            result->setFacet2(facet_2);
            result->setEdge21(edge_21);
            result->setEdge22(edge_22);
            VertexSPtr vertices[4];
            result->getVertices(vertices);
            for (unsigned int i = 0; i < 4; i++) {
                SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
                        vertices[i]->getData());
                ArcSPtr arc = vertex_data->getArc();
                node->addArc(arc);
            }
            EdgeSPtr edges[4];
            result->getEdges(edges);
            for (unsigned int i = 0; i < 4; i++) {
                SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
                        edges[i]->getData());
                SheetSPtr sheet = edge_data->getSheet();
                node->addSheet(sheet);
            }

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }
    return result;
}

DblTriangleEventSPtr SimpleStraightSkel::nextDblTriangleEvent(PolyhedronSPtr polyhedron,
                                                              const CGAL::FT current_offset,
                                                              CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Dbl Triangle Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    DblTriangleEventSPtr result = DblTriangleEventSPtr();

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (isTetrahedron(edge)) {
            continue;
        }
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        if (!facet_l || !facet_r) {
            continue;
        }
        if (!(isTriangle(facet_l, edge) &&
                isTriangle(facet_r, edge))) {
            continue;
        }

        Point3SPtr point = Point3SPtr();
        CGAL::FT offset_event;
        std::tie(point, offset_event) = vanishesAt(edge, current_offset);
        if (!point) {
            continue;
        }
        if (offset_event < 0 && offset_event > current_offset_to_nearest_event) {
            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = DblTriangleEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setEdge(edge);

#ifndef CGAL_SS3_NO_SKELETON_DS
            VertexSPtr vertices[4];
            result->getVertices(vertices);
            for (unsigned int i = 0; i < 4; i++) {
                SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(
                        vertices[i]->getData());
                ArcSPtr arc = data->getArc();
                node->addArc(arc);
            }
            EdgeSPtr edges[5];
            result->getEdges(edges);
            for (unsigned int i = 0; i < 5; i++) {
                SkelEdgeDataSPtr data = std::dynamic_pointer_cast<SkelEdgeData>(
                        edges[i]->getData());
                SheetSPtr sheet = data->getSheet();
                node->addSheet(sheet);
            }
#endif

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }
    return result;
}

TetrahedronEventSPtr SimpleStraightSkel::nextTetrahedronEvent(PolyhedronSPtr polyhedron,
                                                              const CGAL::FT current_offset,
                                                              CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Tetrahedron Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    TetrahedronEventSPtr result = TetrahedronEventSPtr();

    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (isTetrahedron(edge)) {
            FacetSPtr facet = edge->getFacetL();
            Point3SPtr point = Point3SPtr();
            CGAL::FT offset_event;
            std::tie(point, offset_event) = vanishesAt(edge, current_offset);
            if (!point) {
                continue;
            }
            if (offset_event < 0 && offset_event > current_offset_to_nearest_event) {
                NodeSPtr node;
                if (!result) {
                    node = Node::create(point);
                    result = TetrahedronEvent::create(polyhedron);
                    result->setNode(node);
                }
                node = result->getNode();
                node->clear();
                node->setOffset(current_offset + offset_event);
                node->setPoint(point);
                result->setEdgeBegin(edge);

#ifndef CGAL_SS3_NO_SKELETON_DS
                VertexSPtr vertices[4];
                result->getVertices(vertices);
                for (unsigned int i = 0; i < 4; i++) {
                    SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
                            vertices[i]->getData());
                    ArcSPtr arc = vertex_data->getArc();
                    node->addArc(arc);
                }
                EdgeSPtr edges[6];
                result->getEdges(edges);
                for (unsigned int i = 0; i < 6; i++) {
                    SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
                            edges[i]->getData());
                    SheetSPtr sheet = edge_data->getSheet();
                    node->addSheet(sheet);
                }
#endif

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
                current_offset_to_nearest_event = offset_event;
#endif
            }
        }
    }
    return result;
}

VertexEventSPtr SimpleStraightSkel::nextVertexEvent(PolyhedronSPtr polyhedron,
                                                    const CGAL::FT current_offset,
                                                    CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Vertex Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    VertexEventSPtr result = VertexEventSPtr();

    std::list<VertexSPtr>::iterator it_v1 = polyhedron->vertices().begin();
    while (it_v1 != polyhedron->vertices().end()) {
        VertexSPtr vertex_1 = *it_v1++;

        if (isLocked(vertex_1)) {
            continue;
        }

        if (isConvex(vertex_1)) {
            continue;
        }

        std::list<VertexSPtr> vertices_2;
        std::list<FacetWPtr>::iterator it_f = vertex_1->facets().begin();
        while (it_f != vertex_1->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (!facet_wptr.expired()) {
                FacetSPtr facet(facet_wptr);
                vertices_2.insert(vertices_2.end(),
                        facet->vertices().begin(), facet->vertices().end());
            }
        }
        std::list<VertexSPtr>::iterator it_v2 = vertices_2.begin();
        while (it_v2 != vertices_2.end()) {
            VertexSPtr vertex_2 = *it_v2++;
            if (vertex_1 == vertex_2) {
                continue;
            }
            // @fixme is this really meant to be a shared ptr comparison?
            if (vertex_1->getPoint() == vertex_2->getPoint()) {
                continue;
            }
            if (isLocked(vertex_2)) {
                continue;
            }
            if (isConvex(vertex_2)) {
                continue;
            }

            if (vertex_1->findEdge(vertex_2)) {
                // edge event
                continue;
            }

            FacetSPtr facet_1;
            FacetSPtr facet_2;
            int num_equal_facets = 0;
            std::list<FacetWPtr>::iterator it_f1 = vertex_1->facets().begin();
            while (it_f1 != vertex_1->facets().end()) {
                FacetWPtr facet_1_wptr = *it_f1++;
                if (!facet_1_wptr.expired()) {
                    std::list<FacetWPtr>::iterator it_f2 = vertex_2->facets().begin();
                    while (it_f2 != vertex_2->facets().end()) {
                        FacetWPtr facet_2_wptr = *it_f2++;
                        if (facet_1_wptr == facet_2_wptr) {
                            if (num_equal_facets == 0) {
                                facet_1 = FacetSPtr(facet_1_wptr);
                            } else {
                                facet_2 = FacetSPtr(facet_2_wptr);
                            }
                            num_equal_facets++;
                        }
                    }
                }
            }
            if (num_equal_facets != 2) {
                continue;
            }
            if (facet_1->next(vertex_1) != facet_2) {
                FacetSPtr facet_tmp = facet_1;
                facet_1 = facet_2;
                facet_2 = facet_tmp;
            }
            if (vertex_1->next(facet_1)->next(facet_1) == vertex_2 ||
                    vertex_1->next(facet_2)->next(facet_2) == vertex_2 ||
                    vertex_1->prev(facet_1)->prev(facet_1) == vertex_2 ||
                    vertex_1->prev(facet_2)->prev(facet_2) == vertex_2) {
                // edge merge event
                continue;
            }

            EdgeSPtr edge_11 = EdgeSPtr();
            EdgeSPtr edge_12 = EdgeSPtr();
            std::list<EdgeWPtr>::iterator it_e1 = vertex_1->edges().begin();
            while (it_e1 != vertex_1->edges().end()) {
                EdgeWPtr edge_1_wptr = *it_e1++;
                if (!edge_1_wptr.expired()) {
                    EdgeSPtr edge_1 = EdgeSPtr(edge_1_wptr);
                    FacetSPtr facet_1l = edge_1->getFacetL();
                    FacetSPtr facet_1r = edge_1->getFacetR();
                    if ((facet_1l == facet_1 && facet_1r != facet_2) ||
                            (facet_1r == facet_1 && facet_1l != facet_2)) {
                        edge_11 = edge_1;
                    } else if ((facet_1l == facet_2 && facet_1r != facet_1) ||
                            (facet_1r == facet_2 && facet_1l != facet_1)) {
                        edge_12 = edge_1;
                    }
                }
            }
            EdgeSPtr edge_21 = EdgeSPtr();
            EdgeSPtr edge_22 = EdgeSPtr();
            std::list<EdgeWPtr>::iterator it_e2 = vertex_2->edges().begin();
            while (it_e2 != vertex_2->edges().end()) {
                EdgeWPtr edge_2_wptr = *it_e2++;
                if (!edge_2_wptr.expired()) {
                    EdgeSPtr edge_2 = EdgeSPtr(edge_2_wptr);
                    FacetSPtr facet_2l = edge_2->getFacetL();
                    FacetSPtr facet_2r = edge_2->getFacetR();
                    if ((facet_2l == facet_1 && facet_2r != facet_2) ||
                            (facet_2r == facet_1 && facet_2l != facet_2)) {
                        edge_21 = edge_2;
                    } else if ((facet_2l == facet_2 && facet_2r != facet_1) ||
                            (facet_2r == facet_2 && facet_2l != facet_1)) {
                        edge_22 = edge_2;
                    }
                }
            }
            if (!((edge_11->next(vertex_1) == edge_12 && edge_22->next(vertex_2) == edge_21) ||
                    (edge_12->next(vertex_1) == edge_11 && edge_21->next(vertex_2) == edge_22))) {
                // flip vertex event
                continue;
            }
            bool conv_split_event = false;
            FacetSPtr facet_1b = facet_2->next(vertex_1);
            FacetSPtr facet_2b = facet_1->next(vertex_2);
            EdgeSPtr edge_cur = edge_11->next(facet_1b);
            while (edge_cur != edge_11) {
                if ((edge_cur->getFacetL() == facet_1b && edge_cur->getFacetR() == facet_2b) ||
                        (edge_cur->getFacetR() == facet_1b && edge_cur->getFacetL() == facet_2b)) {
                    conv_split_event = true;
                    break;
                }
                edge_cur = edge_cur->next(facet_1b);
            }
            if (conv_split_event) {
                continue;
            }

            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_11, edge_22, current_offset, current_offset_to_nearest_event);
            if (!point)
                continue;
            CGAL_assertion(offset_event > current_offset_to_nearest_event);

            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = VertexEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            SkelVertexDataSPtr data_1 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_1->getData());
            node->addArc(data_1->getArc());
            SkelVertexDataSPtr data_2 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_2->getData());
            node->addArc(data_2->getArc());
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setVertex1(vertex_1);
            result->setVertex2(vertex_2);
            result->setFacet1(facet_1);
            result->setFacet2(facet_2);

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }
    return result;
}

FlipVertexEventSPtr SimpleStraightSkel::nextFlipVertexEvent(PolyhedronSPtr polyhedron,
                                                            const CGAL::FT current_offset,
                                                            CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Flip Vertex Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    FlipVertexEventSPtr result = FlipVertexEventSPtr();

    std::list<VertexSPtr>::iterator it_v1 = polyhedron->vertices().begin();
    while (it_v1 != polyhedron->vertices().end()) {
        VertexSPtr vertex_1 = *it_v1++;
        if (isConvex(vertex_1)) {
            continue;
        }

        std::list<VertexSPtr> vertices_2;
        std::list<FacetWPtr>::iterator it_f = vertex_1->facets().begin();
        while (it_f != vertex_1->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (!facet_wptr.expired()) {
                FacetSPtr facet(facet_wptr);
                vertices_2.insert(vertices_2.end(),
                        facet->vertices().begin(), facet->vertices().end());
            }
        }
        std::list<VertexSPtr>::iterator it_v2 = vertices_2.begin();
        while (it_v2 != vertices_2.end()) {
            VertexSPtr vertex_2 = *it_v2++;
            if (vertex_1 == vertex_2) {
                continue;
            }
            if (vertex_1->getPoint() == vertex_2->getPoint()) {
                continue;
            }
            if (isConvex(vertex_2)) {
                continue;
            }

            if (vertex_1->findEdge(vertex_2)) {
                // edge event
                continue;
            }

            FacetSPtr facet_1;
            FacetSPtr facet_2;
            int num_equal_facets = 0;
            std::list<FacetWPtr>::iterator it_f1 = vertex_1->facets().begin();
            while (it_f1 != vertex_1->facets().end()) {
                FacetWPtr facet_1_wptr = *it_f1++;
                if (!facet_1_wptr.expired()) {
                    std::list<FacetWPtr>::iterator it_f2 = vertex_2->facets().begin();
                    while (it_f2 != vertex_2->facets().end()) {
                        FacetWPtr facet_2_wptr = *it_f2++;
                        if (facet_1_wptr == facet_2_wptr) {
                            if (num_equal_facets == 0) {
                                facet_1 = FacetSPtr(facet_1_wptr);
                            } else {
                                facet_2 = FacetSPtr(facet_2_wptr);
                            }
                            num_equal_facets++;
                        }
                    }
                }
            }
            if (num_equal_facets != 2) {
                continue;
            }
            if (facet_1->next(vertex_1) != facet_2) {
                FacetSPtr facet_tmp = facet_1;
                facet_1 = facet_2;
                facet_2 = facet_tmp;
            }
            if (vertex_1->next(facet_1)->next(facet_1) == vertex_2 ||
                    vertex_1->next(facet_2)->next(facet_2) == vertex_2 ||
                    vertex_1->prev(facet_1)->prev(facet_1) == vertex_2 ||
                    vertex_1->prev(facet_2)->prev(facet_2) == vertex_2) {
                // edge merge event
                continue;
            }

            EdgeSPtr edge_11 = EdgeSPtr();
            EdgeSPtr edge_12 = EdgeSPtr();
            std::list<EdgeWPtr>::iterator it_e1 = vertex_1->edges().begin();
            while (it_e1 != vertex_1->edges().end()) {
                EdgeWPtr edge_1_wptr = *it_e1++;
                if (!edge_1_wptr.expired()) {
                    EdgeSPtr edge_1 = EdgeSPtr(edge_1_wptr);
                    FacetSPtr facet_1l = edge_1->getFacetL();
                    FacetSPtr facet_1r = edge_1->getFacetR();
                    if ((facet_1l == facet_1 && facet_1r != facet_2) ||
                            (facet_1r == facet_1 && facet_1l != facet_2)) {
                        edge_11 = edge_1;
                    } else if ((facet_1l == facet_2 && facet_1r != facet_1) ||
                            (facet_1r == facet_2 && facet_1l != facet_1)) {
                        edge_12 = edge_1;
                    }
                }
            }
            EdgeSPtr edge_21 = EdgeSPtr();
            EdgeSPtr edge_22 = EdgeSPtr();
            std::list<EdgeWPtr>::iterator it_e2 = vertex_2->edges().begin();
            while (it_e2 != vertex_2->edges().end()) {
                EdgeWPtr edge_2_wptr = *it_e2++;
                if (!edge_2_wptr.expired()) {
                    EdgeSPtr edge_2 = EdgeSPtr(edge_2_wptr);
                    FacetSPtr facet_2l = edge_2->getFacetL();
                    FacetSPtr facet_2r = edge_2->getFacetR();
                    if ((facet_2l == facet_1 && facet_2r != facet_2) ||
                            (facet_2r == facet_1 && facet_2l != facet_2)) {
                        edge_21 = edge_2;
                    } else if ((facet_2l == facet_2 && facet_2r != facet_1) ||
                            (facet_2r == facet_2 && facet_2l != facet_1)) {
                        edge_22 = edge_2;
                    }
                }
            }
            if (!(edge_12->next(vertex_1) == edge_11 && edge_22->next(vertex_2) == edge_21)) {
                // vertex event
                continue;
            }
            bool conv_split_event = false;
            FacetSPtr facet_1b = facet_2->next(vertex_1);
            FacetSPtr facet_2b = facet_2->next(vertex_2);
            EdgeSPtr edge_cur = edge_11->next(facet_1b);
            while (edge_cur != edge_11) {
                if ((edge_cur->getFacetL() == facet_1b && edge_cur->getFacetR() == facet_2b) ||
                        (edge_cur->getFacetR() == facet_1b && edge_cur->getFacetL() == facet_2b)) {
                    conv_split_event = true;
                    break;
                }
                edge_cur = edge_cur->next(facet_1b);
            }
            if (conv_split_event) {
                continue;
            }

            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_11, edge_22, current_offset, current_offset_to_nearest_event);
            if (!point)
                continue;
            CGAL_assertion(offset_event > current_offset_to_nearest_event);

            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = FlipVertexEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            SkelVertexDataSPtr data_1 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_1->getData());
            node->addArc(data_1->getArc());
            SkelVertexDataSPtr data_2 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_2->getData());
            node->addArc(data_2->getArc());
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setVertex1(vertex_1);
            result->setVertex2(vertex_2);
            result->setFacet1(facet_1);
            result->setFacet2(facet_2);

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }
    return result;
}

SurfaceEventSPtr SimpleStraightSkel::nextSurfaceEvent(PolyhedronSPtr polyhedron,
                                                      const CGAL::FT current_offset,
                                                      CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Surface Event" << std::endl;
    CGAL::Real_timer timer;
    timer.start();

    ReadLock l(polyhedron->mutex());
    SurfaceEventSPtr result = SurfaceEventSPtr();

    std::list<EdgeSPtr>::iterator it_e1 = polyhedron->edges().begin();
    while (it_e1 != polyhedron->edges().end()) {
        EdgeSPtr edge_1 = *it_e1++;
        FacetSPtr facet_1_src = getFacetSrc(edge_1);
        FacetSPtr facet_1_dst = getFacetDst(edge_1);
        std::list<EdgeSPtr> edges_2;
        edges_2.insert(edges_2.end(),
                facet_1_src->edges().begin(), facet_1_src->edges().end());
        edges_2.insert(edges_2.end(),
                facet_1_dst->edges().begin(), facet_1_dst->edges().end());
        std::list<EdgeSPtr>::iterator it_e2 = edges_2.begin();
        while (it_e2 != edges_2.end()) {
            EdgeSPtr edge_2 = *it_e2++;
            if (edge_1 == edge_2) {
                continue;
            }
            if (edge_1->getFacetL() == edge_2->getFacetL() ||
                    edge_1->getFacetL() == edge_2->getFacetR() ||
                    edge_1->getFacetR() == edge_2->getFacetL() ||
                    edge_1->getFacetR() == edge_2->getFacetR()) {
                // on same facet
                continue;
            }
            if (edge_1->getVertexSrc()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
                    edge_1->getVertexSrc()->getPoint() == edge_2->getVertexDst()->getPoint() ||
                    edge_1->getVertexDst()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
                    edge_1->getVertexDst()->getPoint() == edge_2->getVertexDst()->getPoint()) {
                // share a vertex
                continue;
            }
            // vertex of edge_1 splits edge_2
            if (!((edge_2->getFacetL() == facet_1_src && edge_2->getFacetR() != facet_1_dst) ||
                    (edge_2->getFacetL() == facet_1_dst && edge_2->getFacetR() != facet_1_src) ||
                    (edge_2->getFacetR() == facet_1_src && edge_2->getFacetL() != facet_1_dst) ||
                    (edge_2->getFacetR() == facet_1_dst && edge_2->getFacetL() != facet_1_src))) {
                // no surface event
                continue;
            }
            FacetSPtr facet_2_src = getFacetSrc(edge_2);
            FacetSPtr facet_2_dst = getFacetDst(edge_2);
            if ((edge_1->getFacetL() == facet_2_src && edge_1->getFacetR() != facet_2_dst) ||
                    (edge_1->getFacetL() == facet_2_dst && edge_1->getFacetR() != facet_2_src) ||
                    (edge_1->getFacetR() == facet_2_src && edge_1->getFacetL() != facet_2_dst) ||
                    (edge_1->getFacetR() == facet_2_dst && edge_1->getFacetL() != facet_2_src)) {
                // flip vertex event
                continue;
            }
            if (edge_1->getVertexSrc()->findEdge(edge_2->getVertexSrc()) ||
                    edge_1->getVertexSrc()->findEdge(edge_2->getVertexDst()) ||
                    edge_1->getVertexDst()->findEdge(edge_2->getVertexSrc()) ||
                    edge_1->getVertexDst()->findEdge(edge_2->getVertexDst()) ) {
                // edge event (when a pyramid grows outwards)
                // a surface split is not possible with only one edge in between
                continue;
            }
            if ((edge_1->getFacetL() == facet_2_src && facet_1_src == edge_2->getFacetL()) ||
                    (edge_1->getFacetL() == facet_2_dst && facet_1_src == edge_2->getFacetR()) ||
                    (edge_1->getFacetR() == facet_2_src && facet_1_dst == edge_2->getFacetL()) ||
                    (edge_1->getFacetR() == facet_2_dst && facet_1_dst == edge_2->getFacetR()) ||
                    (edge_1->getFacetR() == facet_2_src && facet_1_src == edge_2->getFacetR()) ||
                    (edge_1->getFacetR() == facet_2_dst && facet_1_src == edge_2->getFacetL()) ||
                    (edge_1->getFacetL() == facet_2_src && facet_1_dst == edge_2->getFacetR()) ||
                    (edge_1->getFacetL() == facet_2_dst && facet_1_dst == edge_2->getFacetL())) {
                // vertex event
                continue;
            }
            bool is_conv_split_event = false;
            std::list<EdgeSPtr> edges = edge_1->getFacetL()->findEdges(edge_1->getFacetR());
            std::list<EdgeSPtr>::iterator it_e = edges.begin();
            while (it_e != edges.end()) {
                EdgeSPtr edge = *it_e++;
                if (edge == edge_1) {
                    continue;
                }
                FacetSPtr facet_src = getFacetSrc(edge);
                FacetSPtr facet_dst = getFacetDst(edge);
                if (facet_1_src == edge_2->getFacetL() ||
                        facet_1_dst == edge_2->getFacetL()) {
                    if (facet_src == edge_2->getFacetR() ||
                            facet_dst == edge_2->getFacetR()) {
                        is_conv_split_event = true;
                        break;
                    }
                } else if (facet_1_src == edge_2->getFacetR() ||
                        facet_1_dst == edge_2->getFacetR()) {
                    if (facet_src == edge_2->getFacetL() ||
                            facet_dst == edge_2->getFacetL()) {
                        is_conv_split_event = true;
                        break;
                    }
                }
            }
            if (is_conv_split_event) {
                continue;
            }

            // calculate intersection point
            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_1, edge_2, current_offset, current_offset_to_nearest_event);
            if (!point)
                continue;

            CGAL_assertion(offset_event > current_offset_to_nearest_event);

            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = SurfaceEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setEdge1(edge_1);
            result->setEdge2(edge_2);

            SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(
                    edge_1->getData());
            SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(
                    edge_2->getData());
            node->addSheet(data_1->getSheet());
            node->addSheet(data_2->getSheet());

            if (facet_1_src == edge_2->getFacetL() ||
                    facet_1_src == edge_2->getFacetR()) {
                SkelVertexDataSPtr data_1_src = std::dynamic_pointer_cast<SkelVertexData>(
                    edge_1->getVertexSrc()->getData());
                node->addArc(data_1_src->getArc());
            }
            if (facet_1_dst == edge_2->getFacetL() ||
                    facet_1_dst == edge_2->getFacetR()) {
                SkelVertexDataSPtr data_1_dst = std::dynamic_pointer_cast<SkelVertexData>(
                    edge_1->getVertexDst()->getData());
                node->addArc(data_1_dst->getArc());
            }

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }

    timer.stop();
    std::cout << "Sought Surface Event in: " << timer.time() << std::endl;

    return result;
}

PolyhedronSplitEventSPtr SimpleStraightSkel::nextPolyhedronSplitEvent(PolyhedronSPtr polyhedron,
                                                                      const CGAL::FT current_offset,
                                                                      CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Polyhedron Split Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    PolyhedronSplitEventSPtr result = PolyhedronSplitEventSPtr();

    std::list<EdgeSPtr>::iterator it_e1 = polyhedron->edges().begin();
    while (it_e1 != polyhedron->edges().end()) {
        EdgeSPtr edge_1 = *it_e1++;
        std::cout << "consider edge " << edge_1->toString() << std::endl;
        if (isLocked(edge_1)) {
            continue;
        }
        if (!isReflex(edge_1)) {
            continue;
        }
        FacetSPtr facet_1_src = getFacetSrc(edge_1);
        FacetSPtr facet_1_dst = getFacetDst(edge_1);
        std::list<EdgeSPtr>::iterator it_e2 = facet_1_src->edges().begin();
        while (it_e2 != facet_1_src->edges().end()) {
            EdgeSPtr edge_2 = *it_e2++;
            if (edge_1->getVertexSrc()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
                    edge_1->getVertexSrc()->getPoint() == edge_2->getVertexDst()->getPoint() ||
                    edge_1->getVertexDst()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
                    edge_1->getVertexDst()->getPoint() == edge_2->getVertexDst()->getPoint()) {
                // share a vertex
                continue;
            }
            if (!((edge_2->getFacetL() == facet_1_src && edge_2->getFacetR() == facet_1_dst) ||
                    (edge_2->getFacetL() == facet_1_dst && edge_2->getFacetR() == facet_1_src))) {
                // no polyhedron split event
                continue;
            }
            if (edge_1->getVertexSrc()->findEdge(edge_2->getVertexSrc()) ||
                    edge_1->getVertexSrc()->findEdge(edge_2->getVertexDst()) ||
                    edge_1->getVertexDst()->findEdge(edge_2->getVertexSrc()) ||
                    edge_1->getVertexDst()->findEdge(edge_2->getVertexDst())) {
                // does not work when there is only one edge in between
                continue;
            }

            // calculate intersection point
            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_1, edge_2, current_offset, current_offset_to_nearest_event);
            if (!point)
                continue;

            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = PolyhedronSplitEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setEdge1(edge_1);
            result->setEdge2(edge_2);

            SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(
                    edge_1->getData());
            SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(
                    edge_2->getData());
            node->addSheet(data_1->getSheet());
            node->addSheet(data_2->getSheet());

            if (facet_1_src == edge_2->getFacetL() ||
                    facet_1_src == edge_2->getFacetR()) {
                SkelVertexDataSPtr data_1_src = std::dynamic_pointer_cast<SkelVertexData>(
                    edge_1->getVertexSrc()->getData());
                node->addArc(data_1_src->getArc());
            }
            if (facet_1_dst == edge_2->getFacetL() ||
                    facet_1_dst == edge_2->getFacetR()) {
                SkelVertexDataSPtr data_1_dst = std::dynamic_pointer_cast<SkelVertexData>(
                    edge_1->getVertexDst()->getData());
                node->addArc(data_1_dst->getArc());
            }

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }
    return result;
}

SplitMergeEventSPtr SimpleStraightSkel::nextSplitMergeEvent(PolyhedronSPtr polyhedron,
                                                            const CGAL::FT current_offset,
                                                            CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Split Merge Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    SplitMergeEventSPtr result = SplitMergeEventSPtr();

    std::list<VertexSPtr>::iterator it_v1 = polyhedron->vertices().begin();
    while (it_v1 != polyhedron->vertices().end()) {
        VertexSPtr vertex_1 = *it_v1++;
        if (isConvex(vertex_1)) {
            continue;
        }

        std::list<VertexSPtr> vertices_2;
        std::list<FacetWPtr>::iterator it_f = vertex_1->facets().begin();
        while (it_f != vertex_1->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (!facet_wptr.expired()) {
                FacetSPtr facet(facet_wptr);
                vertices_2.insert(vertices_2.end(),
                        facet->vertices().begin(), facet->vertices().end());
            }
        }
        std::list<VertexSPtr>::iterator it_v2 = vertices_2.begin();
        while (it_v2 != vertices_2.end()) {
            VertexSPtr vertex_2 = *it_v2++;
            if (vertex_1 == vertex_2) {
                continue;
            }
            if (vertex_1->getPoint() == vertex_2->getPoint()) {
                continue;
            }
            if (isConvex(vertex_2)) {
                continue;
            }

            if (vertex_1->findEdge(vertex_2)) {
                // edge event
                continue;
            }

            FacetSPtr facet_1;
            FacetSPtr facet_2;
            int num_equal_facets = 0;
            std::list<FacetWPtr>::iterator it_f1 = vertex_1->facets().begin();
            while (it_f1 != vertex_1->facets().end()) {
                FacetWPtr facet_1_wptr = *it_f1++;
                if (!facet_1_wptr.expired()) {
                    std::list<FacetWPtr>::iterator it_f2 = vertex_2->facets().begin();
                    while (it_f2 != vertex_2->facets().end()) {
                        FacetWPtr facet_2_wptr = *it_f2++;
                        if (facet_1_wptr == facet_2_wptr) {
                            if (num_equal_facets == 0) {
                                facet_1 = FacetSPtr(facet_1_wptr);
                            } else {
                                facet_2 = FacetSPtr(facet_2_wptr);
                            }
                            num_equal_facets++;
                        }
                    }
                }
            }
            if (num_equal_facets != 2) {
                continue;
            }
            if (facet_1->next(vertex_1) != facet_2) {
                FacetSPtr facet_tmp = facet_1;
                facet_1 = facet_2;
                facet_2 = facet_tmp;
            }
            if (vertex_1->next(facet_1)->next(facet_1) == vertex_2 ||
                    vertex_1->next(facet_2)->next(facet_2) == vertex_2 ||
                    vertex_1->prev(facet_1)->prev(facet_1) == vertex_2 ||
                    vertex_1->prev(facet_2)->prev(facet_2) == vertex_2) {
                // edge merge event
                continue;
            }

            EdgeSPtr edge_11 = EdgeSPtr();
            EdgeSPtr edge_12 = EdgeSPtr();
            std::list<EdgeWPtr>::iterator it_e1 = vertex_1->edges().begin();
            while (it_e1 != vertex_1->edges().end()) {
                EdgeWPtr edge_1_wptr = *it_e1++;
                if (!edge_1_wptr.expired()) {
                    EdgeSPtr edge_1 = EdgeSPtr(edge_1_wptr);
                    FacetSPtr facet_1l = edge_1->getFacetL();
                    FacetSPtr facet_1r = edge_1->getFacetR();
                    if ((facet_1l == facet_1 && facet_1r != facet_2) ||
                            (facet_1r == facet_1 && facet_1l != facet_2)) {
                        edge_11 = edge_1;
                    } else if ((facet_1l == facet_2 && facet_1r != facet_1) ||
                            (facet_1r == facet_2 && facet_1l != facet_1)) {
                        edge_12 = edge_1;
                    }
                }
            }
            EdgeSPtr edge_21 = EdgeSPtr();
            EdgeSPtr edge_22 = EdgeSPtr();
            std::list<EdgeWPtr>::iterator it_e2 = vertex_2->edges().begin();
            while (it_e2 != vertex_2->edges().end()) {
                EdgeWPtr edge_2_wptr = *it_e2++;
                if (!edge_2_wptr.expired()) {
                    EdgeSPtr edge_2 = EdgeSPtr(edge_2_wptr);
                    FacetSPtr facet_2l = edge_2->getFacetL();
                    FacetSPtr facet_2r = edge_2->getFacetR();
                    if ((facet_2l == facet_1 && facet_2r != facet_2) ||
                            (facet_2r == facet_1 && facet_2l != facet_2)) {
                        edge_21 = edge_2;
                    } else if ((facet_2l == facet_2 && facet_2r != facet_1) ||
                            (facet_2r == facet_2 && facet_2l != facet_1)) {
                        edge_22 = edge_2;
                    }
                }
            }
            bool conv_split_event = false;
            FacetSPtr facet_1b = facet_2->next(vertex_1);
            FacetSPtr facet_2b = facet_1->next(vertex_2);
            if (facet_2b == facet_2) {
                // flip vertex event
                facet_2b = facet_2b->next(vertex_2);
            }
            EdgeSPtr edge_cur = edge_11->next(facet_1b);
            while (edge_cur != edge_11) {
                if ((edge_cur->getFacetL() == facet_1b && edge_cur->getFacetR() == facet_2b) ||
                        (edge_cur->getFacetR() == facet_1b && edge_cur->getFacetL() == facet_2b)) {
                    conv_split_event = true;
                    break;
                }
                edge_cur = edge_cur->next(facet_1b);
            }
            if (!conv_split_event) {
                continue;
            }

            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_11, edge_22, current_offset, current_offset_to_nearest_event);
            if (!point)
                continue;

            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = SplitMergeEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            SkelVertexDataSPtr data_1 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_1->getData());
            node->addArc(data_1->getArc());
            SkelVertexDataSPtr data_2 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex_2->getData());
            node->addArc(data_2->getArc());
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setVertex1(vertex_1);
            result->setVertex2(vertex_2);
            result->setFacet1(facet_1);
            result->setFacet2(facet_2);

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }
    return result;
}

EdgeSplitEventSPtr SimpleStraightSkel::nextEdgeSplitEvent(PolyhedronSPtr polyhedron,
                                                          const CGAL::FT current_offset,
                                                          CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Edge Split Event" << std::endl;
    CGAL::Real_timer timer;
    timer.start();

    ReadLock l(polyhedron->mutex());
    EdgeSplitEventSPtr result = EdgeSplitEventSPtr();

    std::list<EdgeSPtr> edges_reflex;
    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (isLocked(edge)) {
            continue; // I have no idea what this does
        }
        if (isReflex(edge)) {
            edges_reflex.push_back(edge);
        }
    }

    std::list<EdgeSPtr>::iterator it_e1 = edges_reflex.begin();
    while (it_e1 != edges_reflex.end()) {
        EdgeSPtr edge_1 = *it_e1++;
        FacetSPtr facet_1_src = getFacetSrc(edge_1);
        FacetSPtr facet_1_dst = getFacetDst(edge_1);

        std::list<EdgeSPtr>::iterator it_e2 = it_e1;
        while (it_e2 != edges_reflex.end()) {
            EdgeSPtr edge_2 = *it_e2++;
            if (edge_1->getFacetL() == edge_2->getFacetL() ||
                    edge_1->getFacetL() == edge_2->getFacetR() ||
                    edge_1->getFacetR() == edge_2->getFacetL() ||
                    edge_1->getFacetR() == edge_2->getFacetR()) {
                // on same facet
                continue;
            }
            if (edge_1->getVertexSrc()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
                    edge_1->getVertexSrc()->getPoint() == edge_2->getVertexDst()->getPoint() ||
                    edge_1->getVertexDst()->getPoint() == edge_2->getVertexSrc()->getPoint() ||
                    edge_1->getVertexDst()->getPoint() == edge_2->getVertexDst()->getPoint()) {
                // share a vertex
                continue;
            }
            if (((edge_2->getFacetL() == facet_1_src && edge_2->getFacetR() == facet_1_dst) ||
                    (edge_2->getFacetL() == facet_1_dst && edge_2->getFacetR() == facet_1_src))) {
                // polyhedron split event
                continue;
            }
            FacetSPtr facet_2_src = getFacetSrc(edge_2);
            FacetSPtr facet_2_dst = getFacetDst(edge_2);
            if ((edge_2->getFacetL() == facet_1_src && edge_2->getFacetR() != facet_1_dst) ||
                    (edge_2->getFacetL() == facet_1_dst && edge_2->getFacetR() != facet_1_src) ||
                    (edge_2->getFacetR() == facet_1_src && edge_2->getFacetL() != facet_1_dst) ||
                    (edge_2->getFacetR() == facet_1_dst && edge_2->getFacetL() != facet_1_src)) {
                // surface event
                continue;
            }
            if ((edge_1->getFacetL() == facet_2_src && edge_1->getFacetR() != facet_2_dst) ||
                    (edge_1->getFacetL() == facet_2_dst && edge_1->getFacetR() != facet_2_src) ||
                    (edge_1->getFacetR() == facet_2_src && edge_1->getFacetL() != facet_2_dst) ||
                    (edge_1->getFacetR() == facet_2_dst && edge_1->getFacetL() != facet_2_src)) {
                // surface event
                continue;
            }

            // calculate intersection point
            Point3SPtr point;
            CGAL::FT offset_event;
            std::tie(point, offset_event) = crashAt(edge_1, edge_2, current_offset, current_offset_to_nearest_event);
            if (!point)
                continue;

            NodeSPtr node;
            if (!result) {
                node = Node::create(point);
                result = EdgeSplitEvent::create(polyhedron);
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(current_offset + offset_event);
            node->setPoint(point);
            result->setEdge1(edge_1);
            result->setEdge2(edge_2);

            SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(
                    edge_1->getData());
            SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(
                    edge_2->getData());
            node->addSheet(data_1->getSheet());
            node->addSheet(data_2->getSheet());

            if (facet_1_src == edge_2->getFacetL() ||
                    facet_1_src == edge_2->getFacetR()) {
                SkelVertexDataSPtr data_1_src = std::dynamic_pointer_cast<SkelVertexData>(
                    edge_1->getVertexSrc()->getData());
                node->addArc(data_1_src->getArc());
            }
            if (facet_1_dst == edge_2->getFacetL() ||
                    facet_1_dst == edge_2->getFacetR()) {
                SkelVertexDataSPtr data_1_dst = std::dynamic_pointer_cast<SkelVertexData>(
                    edge_1->getVertexDst()->getData());
                node->addArc(data_1_dst->getArc());
            }

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
            current_offset_to_nearest_event = offset_event;
#endif
        }
    }

    timer.stop();
    std::cout << "Sought Edge Split Event in: " << timer.time() << std::endl;

    return result;
}

PierceEventSPtr SimpleStraightSkel::nextPierceEvent(PolyhedronSPtr polyhedron,
                                                    const CGAL::FT current_offset,
                                                    CGAL::FT& current_offset_to_nearest_event)
{
    std::cout << ">>> Seek -- Next Pierce Event" << std::endl;

    ReadLock l(polyhedron->mutex());
    PierceEventSPtr result = PierceEventSPtr();

    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        if (isLocked(vertex)) {
            continue;
        }

        if (isReflex(vertex)) {
            SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
            ArcSPtr arc = data->getArc();

            std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
            while (it_f != polyhedron->facets().end()) {
                FacetSPtr facet = *it_f++;
                bool contains_vertex = false;
                std::list<VertexSPtr>::iterator it_v2 = facet->vertices().begin();
                while (it_v2 != facet->vertices().end()) {
                    VertexSPtr vertex_2 = *it_v2++;
                    if (vertex_2->getPoint() == vertex->getPoint()) {
                        contains_vertex = true;
                        break;
                    }
                }
                if (contains_vertex) {
                    continue;
                }

                bool has_edge_to_facet = false;
                std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
                while (it_e != vertex->edges().end()) {
                    EdgeWPtr edge_wptr = *it_e++;
                    if (!edge_wptr.expired()) {
                        EdgeSPtr edge(edge_wptr);
                        // @todo checking both because we don't know if vertex is the src or the dst of the edge?
                        FacetSPtr facet_src = edge->getFacetL()->next(edge->getVertexSrc());
                        FacetSPtr facet_dst = edge->getFacetR()->next(edge->getVertexDst());
                        if (facet == facet_src || facet == facet_dst) {
                            has_edge_to_facet = true;
                            break;
                        }
                    }
                }
                if (has_edge_to_facet) {
                    continue;
                }

                // @fixme assumes positive weights
                if (KernelWrapper::side(facet->plane(), vertex->getPoint()) > 0) {
                    continue;
                }

                Point3SPtr point;
                CGAL::FT offset_event;

// #define CGAL_SS3_OLD_CODE_PIERCE_EVENT
#ifdef CGAL_SS3_OLD_CODE_PIERCE_EVENT
                // @fixme this filter seems overly restrictive?
                // Could not we have an intersection between the OFFSET face
                // and the arc's supporting line (but we don't know at what offset yet)?
                if (!IsLineInFacet(facet, arc->line())) {
                    continue;
                }

                FacetSPtr facet_vertex = FacetSPtr(vertex->facets().front());
                CGAL::FT facet_speed_vertex = 1.0;
                if (facet_vertex->hasData()) {
                    facet_speed_vertex = std::dynamic_pointer_cast<SkelFacetData>(
                            facet_vertex->getData())->getSpeed();
                }
                Plane3SPtr plane_vertex_offset = KernelWrapper::offsetPlane(facet_vertex->plane(), -facet_speed_vertex);
                Point3SPtr point_vertex_offset = KernelWrapper::intersection(plane_vertex_offset, arc->line());
                CGAL::FT speed_vertex = KernelWrapper::distance(vertex->getPoint(), point_vertex_offset);

                Point3SPtr point_facet = KernelWrapper::intersection(facet->plane(), arc->line());
                CGAL::FT facet_speed = 1.0;
                if (facet->hasData()) {
                    facet_speed = std::dynamic_pointer_cast<SkelFacetData>(
                            facet->getData())->getSpeed();
                }
                Plane3SPtr plane_facet_offset = KernelWrapper::offsetPlane(facet->plane(), -facet_speed);
                Point3SPtr point_facet_offset = KernelWrapper::intersection(plane_facet_offset, arc->line());
                CGAL::FT speed_facet = KernelWrapper::distance(point_facet, point_facet_offset);

                CGAL::FT distance = KernelWrapper::distance(vertex->getPoint(), point_facet);
                CGAL::FT dist_vertex = (distance * speed_vertex) / (speed_vertex + speed_facet);
                if (KernelWrapper::comparePoints(arc->getDirection(),
                        point_facet, point_facet_offset) < 0) {
                    // for weighted straight skeleton
                    // reflex vertex and facet move into same direction
                    if (speed_facet < speed_vertex) {
                        // facet to slow
                        continue;
                    }
                    dist_vertex = (distance * speed_vertex) / (speed_facet - speed_vertex);
                }

                offset_event = -dist_vertex / speed_vertex;
                point = KernelWrapper::offsetPoint(vertex->getPoint(), arc->getDirection(), dist_vertex);
                // std::cout << "old result: " << *point << " time: " << offset_event << std::endl;
#else
                CGAL_assertion(vertex->facets().size() == 3);

                FacetWPtr wf0 = *(std::next(vertex->facets().begin(), 0)); // @todo can be factorized
                FacetWPtr wf1 = *(std::next(vertex->facets().begin(), 1));
                FacetWPtr wf2 = *(std::next(vertex->facets().begin(), 2));

                CGAL_assertion(!wf0.expired());
                CGAL_assertion(!wf1.expired());
                CGAL_assertion(!wf2.expired());

                FacetSPtr f0(wf0);
                FacetSPtr f1(wf1);
                FacetSPtr f2(wf2);

                Plane3SPtr plane = facet->plane();
                Plane3SPtr plane_0 = f0->plane();
                Plane3SPtr plane_1 = f1->plane();
                Plane3SPtr plane_2 = f2->plane();

                CGAL::FT speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
                CGAL::FT speed_0 = std::dynamic_pointer_cast<SkelFacetData>(f0->getData())->getSpeed();
                CGAL::FT speed_1 = std::dynamic_pointer_cast<SkelFacetData>(f1->getData())->getSpeed();
                CGAL::FT speed_2 = std::dynamic_pointer_cast<SkelFacetData>(f2->getData())->getSpeed();

                std::tie(point, offset_event) = KernelWrapper::intersectionAndTimeOffsetPlanes(
                    plane, speed, plane_0, speed_0, plane_1, speed_1, plane_2, speed_2);

                // std::cout << " **" << std::endl;
                // std::cout << "facet = " << facet->getID() << std::endl;
                // std::cout << "plane = " << *plane << std::endl;
                // std::cout << "f0 = " << f0->getID() << std::endl;
                // std::cout << "f1 = " << f1->getID() << std::endl;
                // std::cout << "f2 = " << f2->getID() << std::endl;
                // std::cout << "ipoint: " << *point << " time: " << offset_event << std::endl;

                // Now, we might naively wish to filter using bisectors like in 2D SLS code,
                // but unlike a segment, a face in the 3D SS code has no reason to be convex,
                // which changes everything and can result in false positives.
                //
                // The bisector filter in 2D is equivalent to checking if the point is on the offset
                // face. We could check this here, but determining what is the offset face at this
                // point (i.e., while searching for events) is rough: plenty of other events
                // might modify the face before this particular pierce event appears, and so
                // we can't just do shift(facet) because the result might be a self-intersecting
                // polygon with holes.
                //
                // Instead, we do not filter here, but simply put it in the queue. When the event
                // will be popped, then we know it's the next event globally and nothing else
                // can mess up the face, and we can do the in-test check then.
                //
                // See HandlePierceEvent()

#endif // CGAL_SS3_OLD_CODE_PIERCE_EVENT
                std::cout << "Accepted event" << std::endl;

                if (offset_event == 0)
                {
                  // @fixme this should be <= like the other nextSomethingEvent
                  // but since we handle filtering at pop time, it will just loop indefinitely.
                  CGAL_assertion_msg(false, "this should have been accepted!!");
                }

                if (offset_event < 0 && offset_event > current_offset_to_nearest_event) {
                    NodeSPtr node;
                    if (!result) {
                        node = Node::create(point);
                        result = PierceEvent::create(polyhedron);
                        result->setNode(node);
                    }
                    node = result->getNode();
                    node->clear();
                    node->addArc(arc);
                    node->setOffset(current_offset + offset_event);
                    node->setPoint(point);
                    result->setFacet(facet);
                    result->setVertex(vertex);

#ifndef CGAL_SS3_DO_NOT_FILTER_FUTURE_EVENTS
                    current_offset_to_nearest_event = offset_event;
#endif
                }
            }
        }
    }
    return result;
}
#endif // individual events of a certain type





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //





Point3SPtr KernelWrapper::intersectionOffsetPlanes(Plane3SPtr plane_0,
                                                   const CGAL::FT& w0,
                                                   Plane3SPtr plane_1,
                                                   const CGAL::FT& w1,
                                                   Plane3SPtr plane_2,
                                                   const CGAL::FT& w2,
                                                   Plane3SPtr plane_3,
                                                   const CGAL::FT& w3)
{
    if(is_zero(w0) && is_zero(w1) && is_zero(w2) && is_zero(w3)) {
        return { };
    }

    const CGAL::FT& a0 = plane_0->a();
    const CGAL::FT& b0 = plane_0->b();
    const CGAL::FT& c0 = plane_0->c();
    const CGAL::FT& d0 = plane_0->d();
    const CGAL::FT& a1 = plane_1->a();
    const CGAL::FT& b1 = plane_1->b();
    const CGAL::FT& c1 = plane_1->c();
    const CGAL::FT& d1 = plane_1->d();
    const CGAL::FT& a2 = plane_2->a();
    const CGAL::FT& b2 = plane_2->b();
    const CGAL::FT& c2 = plane_2->c();
    const CGAL::FT& d2 = plane_2->d();
    const CGAL::FT& a3 = plane_3->a();
    const CGAL::FT& b3 = plane_3->b();
    const CGAL::FT& c3 = plane_3->c();
    const CGAL::FT& d3 = plane_3->d();

// #define CGAL_SS3_DEBUG_PLANES_INTERSECTION
#ifdef CGAL_SS3_DEBUG_PLANES_INTERSECTION
    std::cout << "Coefficients\n" << a0 << " " << b0 << " " << c0 << " " << d0 << "\n"
                                  << a1 << " " << b1 << " " << c1 << " " << d1 << "\n"
                                  << a2 << " " << b2 << " " << c2 << " " << d2 << "\n"
                                  << a3 << " " << b3 << " " << c3 << " " << d3 << std::endl;
    std::cout << "Weights\n" << w0 << " " << w1 << " " << w2 << " " << w3 << std::endl;

    std::cout << "CHECK det " << CGAL::determinant(a0, b0, c0, d0,
                                                   a1, b1, c1, d1,
                                                   a2, b2, c2, d2,
                                                   a3, b3, c3, d3) << std::endl;
#endif

    CGAL_assertion((a0*a0 + b0*b0 + c0*c0 - 1) <= 1e-5);
    CGAL_assertion((a1*a1 + b1*b1 + c1*c1 - 1) <= 1e-5);
    CGAL_assertion((a2*a2 + b2*b2 + c2*c2 - 1) <= 1e-5);
    CGAL_assertion((a3*a3 + b3*b3 + c3*c3 - 1) <= 1e-5);

    CGAL::FT den = (-a0*b1*c2*w3 + a0*b1*c3*w2 + a0*b2*c1*w3 - a0*b2*c3*w1 - a0*b3*c1*w2 + a0*b3*c2*w1 + a1*b0*c2*w3 - a1*b0*c3*w2 - a1*b2*c0*w3 + a1*b2*c3*w0 + a1*b3*c0*w2 - a1*b3*c2*w0 - a2*b0*c1*w3 + a2*b0*c3*w1 + a2*b1*c0*w3 - a2*b1*c3*w0 - a2*b3*c0*w1 + a2*b3*c1*w0 + a3*b0*c1*w2 - a3*b0*c2*w1 - a3*b1*c0*w2 + a3*b1*c2*w0 + a3*b2*c0*w1 - a3*b2*c1*w0);

    util::ConfigurationSPtr config = util::Configuration::getInstance();
    bool usePerturbations = false;
    if (config->isLoaded()) {
        if ((config->contains("main", "rand_move_points") &&
            config->getBool("main", "rand_move_points")) ||
            (config->contains("main", "rand_move_points_when_degenerated") &&
            config->getBool("main", "rand_move_points_when_degenerated"))) {
            usePerturbations = true;
        }
    }

    if (!usePerturbations) {
        std::exit(1); // @tmp
        if(CGAL::is_zero(den))
        {
            std::cerr << "Warning: no solution in 4 shifted plane system" << std::endl;
            return { };
        }
    }

    // warning: only valid for normalized coefficients!!
    CGAL::FT x = (b0*c1*d2*w3 - b0*c1*d3*w2 - b0*c2*d1*w3 + b0*c2*d3*w1 + b0*c3*d1*w2 - b0*c3*d2*w1 - b1*c0*d2*w3 + b1*c0*d3*w2 + b1*c2*d0*w3 - b1*c2*d3*w0 - b1*c3*d0*w2 + b1*c3*d2*w0 + b2*c0*d1*w3 - b2*c0*d3*w1 - b2*c1*d0*w3 + b2*c1*d3*w0 + b2*c3*d0*w1 - b2*c3*d1*w0 - b3*c0*d1*w2 + b3*c0*d2*w1 + b3*c1*d0*w2 - b3*c1*d2*w0 - b3*c2*d0*w1 + b3*c2*d1*w0) / den;

    CGAL::FT y = (-a0*c1*d2*w3 + a0*c1*d3*w2 + a0*c2*d1*w3 - a0*c2*d3*w1 - a0*c3*d1*w2 + a0*c3*d2*w1 + a1*c0*d2*w3 - a1*c0*d3*w2 - a1*c2*d0*w3 + a1*c2*d3*w0 + a1*c3*d0*w2 - a1*c3*d2*w0 - a2*c0*d1*w3 + a2*c0*d3*w1 + a2*c1*d0*w3 - a2*c1*d3*w0 - a2*c3*d0*w1 + a2*c3*d1*w0 + a3*c0*d1*w2 - a3*c0*d2*w1 - a3*c1*d0*w2 + a3*c1*d2*w0 + a3*c2*d0*w1 - a3*c2*d1*w0) / den;

    CGAL::FT z = (a0*b1*d2*w3 - a0*b1*d3*w2 - a0*b2*d1*w3 + a0*b2*d3*w1 + a0*b3*d1*w2 - a0*b3*d2*w1 - a1*b0*d2*w3 + a1*b0*d3*w2 + a1*b2*d0*w3 - a1*b2*d3*w0 - a1*b3*d0*w2 + a1*b3*d2*w0 + a2*b0*d1*w3 - a2*b0*d3*w1 - a2*b1*d0*w3 + a2*b1*d3*w0 + a2*b3*d0*w1 - a2*b3*d1*w0 - a3*b0*d1*w2 + a3*b0*d2*w1 + a3*b1*d0*w2 - a3*b1*d2*w0 - a3*b2*d0*w1 + a3*b2*d1*w0) / den;

#ifdef CGAL_SS3_DEBUG_PLANES_INTERSECTION
    std::cout << "CHECK x|y|z " << x << " " << y << " " << z << std::endl;
#endif

    Point3SPtr result = KernelFactory::createPoint3(x, y, z);
    return result;
}





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //




    // Filter if edges are too far
    //
    // @speed the best filtering would be to first have a loop over all combinations
    // to tighten the bound, and only then compute the crashAt, but below is
    // to get an easy speed-up

    this does not filter enough because edges can still see very far
    also broken after update of vanish/crash functions not returning local offsets

    // check shifted planes and edges orientations: if a shifted edge is still in the negative
    // side of the planes, we can 'continue'
    auto is_shifted_edge_definitely_on_negative_side_of_planes = [&offset_of_nearest_event](EdgeSPtr lhs, EdgeSPtr rhs) {
        // @speed note that here and in other places, we could speed up shifting computations
        // because there is a lot of redundant computations: for example here, the two
        // points are the intersection of 2 planes with the same shifted edge line
        // but we compute it from scratch for both
        Point3SPtr offset_e1so = PolyhedronTransformation::shiftPoint(lhs->getVertexSrc(), offset_of_nearest_event);
        Point3SPtr offset_e1to = PolyhedronTransformation::shiftPoint(lhs->getVertexDst(), offset_of_nearest_event);

        Plane3SPtr offset_pl0 = PolyhedronTransformation::shiftPlane(rhs->getFacetL(), offset_of_nearest_event);
        Plane3SPtr offset_pl1 = PolyhedronTransformation::shiftPlane(rhs->getFacetR(), offset_of_nearest_event);

        auto is_shifted_edge_definitely_on_negative_side_of_plane = [](Point3SPtr os, Point3SPtr od, Plane3SPtr opl, Plane3SPtr other_opl) {
            auto orient_s = KernelWrapper::side(opl, os);
            auto orient_d = KernelWrapper::side(opl, od);
            if (orient_s * orient_d < 0) {
                // the offset edge crosses opl; is it entirely on the negative side of the other plane?
                //
                // @todo we can still have the segment be on the union of the negative sides
                // even if it crosses both planes
                return (KernelWrapper::side(other_opl, os) < 0 && KernelWrapper::side(other_opl, od) < 0);
            } else {
                return (orient_s < 0);
            }
        };

        return (is_shifted_edge_definitely_on_negative_side_of_plane(offset_e1so, offset_e1to, offset_pl0, offset_pl1) ||
                is_shifted_edge_definitely_on_negative_side_of_plane(offset_e1so, offset_e1to, offset_pl1, offset_pl0));
    };

    if (is_shifted_edge_definitely_on_negative_side_of_planes(edge_1, edge_2) ||
        is_shifted_edge_definitely_on_negative_side_of_planes(edge_2, edge_1)) {
        ++filtered_candidates;
        continue;
    } else {
        // std::cout << "Checking possible edge split event\n\t"
        //           << edge_1->toString() << "\n\t"
        //           << edge_2->toString() << std::endl;
    }


// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //













// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //