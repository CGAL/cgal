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







// individual events of a given type
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





edge merge handling
#if 1 // using EPECK for now so we are exact
    if (crashAt(edge_1, edge_b).first || crashAt(edge_b, edge_1).first) {
        // crashAt(...) should be commutative,
        // but because of rounding errors it's not in certain cases.
        // TODO: this may need an improvement
        // fix coordinates because point is exactly on merged edge
        DEBUG_VAL("Fixing coordinates of " << vertex->toString());
        Point3SPtr p = vertex->getPoint();
        Vector3SPtr dir = KernelWrapper::normalize(arc->getDirection());
        Point3SPtr p_fixed = KernelFactory::createPoint3((*p) + (*dir)*0.0001);
        event->getNode()->setPoint(p_fixed);
        vertex->setPoint(p_fixed);
    }
#endif





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



void SimpleStraightSkel::collectLocalEvents(const std::set<FacetSPtr>& base_facet_set,
  PolyhedronSPtr polyhedron,
  const CGAL::FT current_offset,
  PQ& queue)
{
std::cout << "collectLocalEvents(" << base_facet_set.size() << "," << current_offset << ")" << std::endl;

// return collectEvents(polyhedron, current_offset, queue);

#ifdef CGAL_SS3_RUN_TIMERS
CGAL::Real_timer timer;
timer.start();
#endif

std::set<VertexSPtr> vertex_set;
std::set<EdgeSPtr> edge_set;
std::set<FacetSPtr> facet_set;

// @todo this is for the Surface Event: the edge might be in the tagged zone, but we also need
// to grab edges of facets that are incident to extremities of edges of the tagged zone.
// @speed event-specific local sets...?
std::set<FacetSPtr> ring_1_facet_set;
for (FacetSPtr facet : base_facet_set) {
for (VertexSPtr vertex : facet->vertices()) {
for (FacetWPtr extra_wf : vertex->facets()) {
if (FacetSPtr extra_f = extra_wf.lock()) {
ring_1_facet_set.insert(extra_f);
}
}
}
}

// in case the event modifies a face incident to a face that is incident to edge2 of Surface event...
std::set<FacetSPtr> ring_2_facet_set;
for (FacetSPtr facet : ring_1_facet_set) {
for (VertexSPtr vertex : facet->vertices()) {
for (FacetWPtr extra_wf : vertex->facets()) {
if (FacetSPtr extra_f = extra_wf.lock()) {
ring_2_facet_set.insert(extra_f);
}
}
}
}

facet_set.insert(base_facet_set.begin(), base_facet_set.end());
facet_set.insert(ring_2_facet_set.begin(), ring_2_facet_set.end());

for (FacetSPtr facet : facet_set) {
for (VertexSPtr vertex : facet->vertices()) {
vertex_set.insert(vertex);
}
for (EdgeSPtr edge : facet->edges()) {
edge_set.insert(edge);
}
}

// @todo template the collectors or something
#if 1
std::list<VertexSPtr> local_vertices(vertex_set.begin(), vertex_set.end());
std::list<EdgeSPtr> local_edges(edge_set.begin(), edge_set.end());
std::list<FacetSPtr> local_facets(facet_set.begin(), facet_set.end());
#else
// just to check if other mechanisms are OK
std::list<VertexSPtr> local_vertices = polyhedron->vertices();
std::list<EdgeSPtr> local_edges = polyhedron->edges();
std::list<FacetSPtr> local_facets = polyhedron->facets();
#endif

std::cout << "Base Facets [" << base_facet_set.size() << "]:\t";
for(FacetSPtr f : base_facet_set) {
std::cout << " " << f->getID();
}
std::cout << "\nVertices [" << local_vertices.size() << "]:\t";
for(VertexSPtr v : local_vertices) {
std::cout << " " << v->getID();
}
std::cout << "\nEdges [" << local_edges.size() << "]:\t";
for(EdgeSPtr e : local_edges) {
std::cout << " " << e->getID();
}
std::cout << "\nFacets [" << local_facets.size() << "]:\t";
for(FacetSPtr f : local_facets) {
std::cout << " " << f->getID();
}
std::cout << std::endl;

// two types of useless events:
// - events that are in the past:
//     offset > current_offset <--- values are negative and decreasing!
// - events that are stricly later than the current next tentative offset:
//     offset < curr_earliest_next_offset
CGAL::FT offset_of_nearest_event = - std::numeric_limits<double>::max();

// if we stop immediately after the last save event, there is no point registering events
// that are farther away
if (!save_offsets_.empty()) {
util::ConfigurationSPtr config = util::Configuration::getInstance();
if (config->isLoaded()) {
if ((config->contains("main", "stop_after_last_save_event") &&
config->getBool("main", "stop_after_last_save_event"))) {
offset_of_nearest_event = (std::max)(offset_of_nearest_event, save_offsets_.back());
}
}
}

DEBUG_PRINT("Past bound = " << current_offset);
DEBUG_PRINT("Initial future bound = " << offset_of_nearest_event);

#ifdef CGAL_SS3_USE_GENERIC_VANISH_EVENT
collectVanishEvents(local_edges, polyhedron, current_offset, offset_of_nearest_event, queue);
#else
collectEdgeEvents(local_edges, polyhedron, current_offset, offset_of_nearest_event, queue);
collectEdgeMergeEvents(local_edges, polyhedron, current_offset, offset_of_nearest_event, queue);
collectTriangleEvents(local_edges, polyhedron, current_offset, offset_of_nearest_event, queue);
collectDblEdgeMergeEvents(local_edges, polyhedron, current_offset, offset_of_nearest_event, queue);
collectDblTriangleEvents(local_edges, polyhedron, current_offset, offset_of_nearest_event, queue);
collectTetrahedronEvents(local_edges, polyhedron, current_offset, offset_of_nearest_event, queue);
#endif

collectVertexEvents(local_vertices, polyhedron, current_offset, offset_of_nearest_event, queue);
collectFlipVertexEvents(local_vertices, polyhedron, current_offset, offset_of_nearest_event, queue);
collectPolyhedronSplitEvents(local_edges, polyhedron, current_offset, offset_of_nearest_event, queue);
collectSplitMergeEvents(local_vertices, polyhedron, current_offset, offset_of_nearest_event, queue);

collectSurfaceEvents(local_edges, polyhedron, current_offset, offset_of_nearest_event, queue);
collectEdgeSplitEvents(local_edges, polyhedron->edges(), polyhedron, current_offset, offset_of_nearest_event, queue);

// currently need to do it both ways because if the facet is part of the modified zone,
// and the concave vertex is not, the event is still marked as obsolete.
// obviously it would be best not to have to do it both ways...
collectPierceEvents(local_vertices, polyhedron->facets(), polyhedron, current_offset, offset_of_nearest_event, queue);
collectPierceEvents(polyhedron->vertices(), local_facets, polyhedron, current_offset, offset_of_nearest_event, queue);

#ifdef CGAL_SS3_RUN_TIMERS
timer.stop();
std::cout << "Sought All Local Events in: " << timer.time() << std::endl;
#endif

#ifdef CGAL_SS3_DEBUG_PRINT_QUEUE
printQueue(queue);
#endif

// checkQueueCorrectness(queue, polyhedron, current_offset);
}



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //




// #define CGAL_SS3_OLD_CODE_PIERCE_EVENT
#ifdef CGAL_SS3_OLD_CODE_PIERCE_EVENT
                SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
                ArcSPtr arc = data->getArc();

                this branch is bad because of the code below which is overly restrictive:
                the line might not intersect the face yet at T=t0 the current offset,
                but might at T=t1 the offset of intersection
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
                        // facet too slow
                        continue;
                    }
                    dist_vertex = (distance * speed_vertex) / (speed_facet - speed_vertex);
                }

                offset_event = -dist_vertex / speed_vertex;
                point = KernelWrapper::offsetPoint(vertex->getPoint(), arc->getDirection(), dist_vertex);
                // std::cout << "old result: " << *point << " time: " << offset_event << std::endl;
#else
                CGAL_assertion(vertex->facets().size() >= 3);

                FacetSPtr fs[3];
                for (int i = 0; i < 3; ++i) {
                    FacetWPtr wf = *(std::next(vertex->facets().begin(), i));
                    CGAL_assertion(!wf.expired());
                    fs[i] = wf.lock();
                }

                std::tie(point, offset_event) = intersectionPointAndTimeOffsetPlanes(facet, fs[0], fs[1], fs[2], current_offset, offset_of_nearest_event);
                if (!point) {
                    continue;
                }

                // std::cout << "  Filter E" << std::endl;

                CGAL_assertion(offset_event < current_offset && offset_event > offset_of_nearest_event);
#endif // CGAL_SS3_OLD_CODE_PIERCE_EVENT





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



void Polyhedron::rebuild(const std::list<FacetSPtr>& facets) {
  std::map<VertexWPtr, VertexSPtr, std::owner_less<VertexWPtr> > old_to_new_vertex;
  std::map<EdgeWPtr, EdgeSPtr, std::owner_less<EdgeWPtr> > old_to_new_edge;
  std::map<FacetWPtr, FacetSPtr, std::owner_less<FacetWPtr> > old_to_new_facet;

  // First pass: create new elements and build replacement maps
  // Create new elements for input facets
  for(const FacetSPtr& f : facets) {
      FacetSPtr new_f = Facet::create();
      new_f->setID(f->getID());
      new_f->setPlane(f->getPlane());
      new_f->setBasePlaneID(f->getBasePlaneID());
      new_f->cachedPlane_ = f->cachedPlane_;
      new_f->cachedSpeed_ = f->cachedSpeed_;
      if(f->hasData())
          new_f->setData(f->getData());

      old_to_new_facet[f] = new_f;
  }

  for(const FacetSPtr& f : facets) {
      FacetSPtr new_f = old_to_new_facet[f];

      // Process vertices
      for(const VertexSPtr& v : f->vertices()) {
          VertexSPtr new_v;

          auto it = old_to_new_vertex.find(v);
          if(it == old_to_new_vertex.end()) {
              new_v = Vertex::create(v->getPoint());
              new_v->setID(v->getID());
              if(v->hasData())
                  new_v->setData(v->getData());

              old_to_new_vertex[v] = new_v;
              addVertex(new_v);
          } else {
              new_v = it->second;
          }

          new_f->addVertex(new_v);
      }

      // Process edges
      for(const EdgeSPtr& e : f->edges()) {
          EdgeSPtr new_e;

          auto it = old_to_new_edge.find(e);
          if(it == old_to_new_edge.end()) {
              new_e = Edge::create(old_to_new_vertex[e->getVertexSrc()],
                                   old_to_new_vertex[e->getVertexDst()]);
              new_e->setID(e->getID());
              if(e->hasData())
                  new_e->setData(e->getData());

              old_to_new_edge[e] = new_e;
              addEdge(new_e);
          } else {
              new_e = it->second;
          }

          if(e->getFacetL() == f) {
              new_e->setFacetL(new_f);
          }
          if(e->getFacetR() == f) {
              new_e->setFacetR(new_f);
          }
          new_f->addEdge(new_e);
      }

      addFacet(new_f);
  }

  std::cout << "Invalided the following elements:" << std::endl;
  std::cout << "Vertices:";
  for(const auto& [old_v, new_v] : old_to_new_vertex) {
      if(!old_v.expired()) {
        std::cout << " " << old_v.lock()->getID();
      }
  }
  std::cout << std::endl;
  std::cout << "Edges:";
  for(const auto& [old_e, new_e] : old_to_new_edge) {
      if(!old_e.expired()) {
        std::cout << " " << old_e.lock()->getID();
      }
  }
  std::cout << std::endl;
  std::cout << "Facets:";
  for(const auto& [old_f, new_f] : old_to_new_facet) {
      if(!old_f.expired()) {
        std::cout << " " << old_f.lock()->getID();
      }
  }
  std::cout << std::endl;

  // Update references in non-replaced elements
  for(const auto& [old_v, new_v] : old_to_new_vertex) {
      if(VertexSPtr v = old_v.lock()) {
          // Look at incident facets to find those that need updates
          for(const FacetWPtr& fw : v->facets()) {
              if(!fw.expired()) {
                  FacetSPtr f = FacetSPtr(fw);
                  if(old_to_new_facet.find(f) == old_to_new_facet.end()) {
                      // This facet is kept but might need vertex and edge updates
                      std::list<VertexSPtr> new_vertices;
                      std::list<EdgeSPtr> new_edges;

                      // Update vertices
                      for(const VertexSPtr& v : f->vertices()) {
                          auto it = old_to_new_vertex.find(v);
                          new_vertices.push_back(it != old_to_new_vertex.end() ? it->second : v);
                      }

                      // Update edges
                      for(const EdgeSPtr& e : f->edges()) {
                          auto it = old_to_new_edge.find(e);
                          if(it != old_to_new_edge.end()) {
                              new_edges.push_back(it->second);
                              if(e->getFacetL() == f)
                                  it->second->setFacetL(f);
                              if(e->getFacetR() == f)
                                  it->second->setFacetR(f);
                          }
                          else {
                              new_edges.push_back(e);
                              // Update vertices of kept edges
                              auto src_it = old_to_new_vertex.find(e->getVertexSrc());
                              if(src_it != old_to_new_vertex.end())
                                  e->replaceVertexSrc(src_it->second);

                              auto dst_it = old_to_new_vertex.find(e->getVertexDst());
                              if(dst_it != old_to_new_vertex.end())
                                  e->replaceVertexDst(dst_it->second);
                          }
                      }

                      // Replace the lists
                      f->vertices_.clear();
                      f->edges_.clear();
                      for(const auto& v : new_vertices) {
                          f->addVertex(v);
                      }
                      for(const auto& e : new_edges) {
                          f->addEdge(e);
                      }
                  }
              }
          }

          // Update edges not in replaced facets
          for(const EdgeWPtr& ew : v->edges()) {
              if(!ew.expired()) {
                  EdgeSPtr e = EdgeSPtr(ew);
                  if(old_to_new_edge.find(e) == old_to_new_edge.end()) {
                      auto src_it = old_to_new_vertex.find(e->getVertexSrc());
                      if(src_it != old_to_new_vertex.end())
                          e->replaceVertexSrc(src_it->second);

                      auto dst_it = old_to_new_vertex.find(e->getVertexDst());
                      if(dst_it != old_to_new_vertex.end())
                          e->replaceVertexDst(dst_it->second);
                  }
              }
          }
      }
  }

  // Remove old elements
  for(const auto& [old_v, _] : old_to_new_vertex) {
      CGAL_assertion(!old_v.expired());
      vertices_.erase(old_v.lock()->getPolyhedronListIt());
  }
  for(const auto& [old_e, _] : old_to_new_edge) {
      CGAL_assertion(!old_e.expired());
      edges_.erase(old_e.lock()->getPolyhedronListIt());
  }
  for(const auto& f : facets) {
      facets_.erase(f->getPolyhedronListIt());
  }

  // Verify expiration
  for(const auto& [old_v, _] : old_to_new_vertex)
      CGAL_assertion(!old_v.expired());
  for(const auto& [old_e, _] : old_to_new_edge)
      CGAL_assertion(!old_e.expired());
  for(const auto& [old_f, _] : old_to_new_facet)
      CGAL_assertion(!old_f.expired());

  std::cout << "Rebuilt finished" << std::endl;
  CGAL_postcondition(isConsistent());
}




// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //






bool Facet::initPlane() {
  bool result = false;
  std::list<TriangleSPtr>::iterator it_t = triangles_.begin();
  while (it_t != triangles_.end()) {
      TriangleSPtr triangle = *it_t++;
      Plane3SPtr plane = triangle->plane();
      if (plane) {
          plane_ = plane;
          result = true;
          break;
      }
  }
  if (!result) {
      Point3SPtr point_prev;
      Point3SPtr points[3];
      unsigned int i = 0;
      std::list<VertexSPtr>::iterator it_v = vertices_.begin();
      while (i < 3 && it_v != vertices_.end()) {
          VertexSPtr vertex = *it_v++;
          Point3SPtr point = vertex->getPoint();
          if (point_prev != point) {
              points[i] = point;
              i++;
          }
          point_prev = point;
      }
      if (i >= 3) {
          Plane3SPtr plane = KernelFactory::createPlane3(
                  points[0], points[1], points[2]);
          if (plane) {
              plane_ = plane;
              result = true;
          }
      }
  }

  CGAL_assertion(result);
  return result;
}





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //






static PolyhedronSPtr soup_to_polyhedron(const std::vector<Point3>& points,
  const std::vector<std::vector<std::size_t> >& triangles,
  const std::vector<Plane3SPtr>& planes,
  const std::vector<CGAL::FT>& speeds);




// belongs somewhere else, probably
PolyhedronSPtr SimpleStraightSkel::soup_to_polyhedron(const std::vector<Point3>& points,
  const std::vector<std::vector<std::size_t> >& triangles,
  const std::vector<Plane3SPtr>& planes,
  const std::vector<CGAL::FT>& speeds) {
PolyhedronSPtr result = Polyhedron::create();

unsigned int vertex_id_new = 0;
for (const Point3& p : points) {
Point3SPtr point = KernelFactory::createPoint3(p);
VertexSPtr vertex = Vertex::create(point);
vertex->setID(vertex_id_new++);
result->addVertex(vertex);
}

unsigned int facet_id_new = 0;
for (std::size_t tid=0; tid<triangles.size(); ++tid) {
const std::vector<std::size_t>& t = triangles[tid];
CGAL_assertion(t.size() == 3); // @tmp could handle more, but no need right now
VertexSPtr poly_vertices[3];
for (unsigned int i = 0; i < 3; i++) {
poly_vertices[i] = *(std::next(result->vertices().begin(), t[i]));
}

FacetSPtr facet = Facet::create(3, poly_vertices);
facet->setID(facet_id_new++);
facet->setPlane(planes[tid]);
Triangle::create(facet, poly_vertices);

data::_3d::skel::SkelFacetDataSPtr data;
if (facet->hasData()) {
data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
} else {
data = data::_3d::skel::SkelFacetData::create(facet);
}
data->setSpeed(speeds[tid]);

result->addFacet(facet);
}

CGAL_postcondition(result->vertices().size() == points.size());
CGAL_postcondition(result->facets().size() == triangles.size());

return result;
}







// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



void tagFacetsWithStepID(const std::set<FacetSPtr>& facets);

void SimpleStraightSkel::tagFacetsWithStepID(const std::set<FacetSPtr>& facets)
{
    std::cout << "Tag facets:";
    for (FacetSPtr f : facets) {
        std::cout << " " << f->getID();
    }
    std::cout << " with step ID " << step_id_ << std::endl;

    for(FacetSPtr facet : facets) {
        if (facet->hasData()) {
            std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->setStepID(step_id_);
        }
    }
}




// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //




    // don't modify this 'define' without modifying the one in SimpleStraightSkel.cpp
// #define CGAL_SS3_ALWAYS_PERTURB_WITH_PLANE_TILTS
#ifndef CGAL_SS3_ALWAYS_PERTURB_WITH_PLANE_TILTS
    canUsePlaneTilts = all_degree_3;
#endif

// @todo something cleaner (return type from the preprocess function?)
bool forced_tilts = false;
// don't modify this 'define' without modifying the one in PolyhedronTransformation.cpp
// #define CGAL_SS3_ALWAYS_PERTURB_WITH_PLANE_TILTS
#ifdef CGAL_SS3_ALWAYS_PERTURB_WITH_PLANE_TILTS
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
      VertexSPtr vertex = *it_v++;
      if (vertex->degree() != 3) { // splitter is in init()
        DEBUG_PRINT("Tilts were used, but there are degree > 3 vertices like: " << vertex->getID())
        forced_tilts = true;
        break;
      }
    }
#endif


#ifdef CGAL_SS3_ALWAYS_PERTURB_WITH_PLANE_TILTS
        // If we have forced tilts despite having vertices with degree > 3 (pre split),
        // then we have a bit of an awkward situation: the points are not on the planes
        // but if we recompute them, we might have (small) self-intersections because of the
        // tilts. So, we recompute points AND shift forward a bit: we know the future
        // is without self-intersections because the splitter picked valid combinations
        if (forced_tilts) {
            current_offset = -1e-5; // @fixme hardcoded value that is obviously not sound
            DEBUG_PRINT("Dummy tilt to " << current_offset);
            PolyhedronTransformation::shiftFacetsInPlace(polyhedron, current_offset);
            db::_3d::OBJFile::save("results/shift_post_forced_tilt.obj", polyhedron, false /*do not triangulate*/);
            CGAL_assertion(polyhedron->isConsistent());
            CGAL_assertion(!SelfIntersection::hasSelfIntersectingSurface(polyhedron));
        }
#endif





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



/**
 * Modify the geometry of the polyhedron as to avoid degenerate configurations.
 */
static PolyhedronSPtr preprocess(PolyhedronSPtr polyhedron);

PolyhedronSPtr SimpleStraightSkel::preprocess(PolyhedronSPtr polyhedron) {
  bool rand_move_points = false;
  bool rand_move_points_when_degenerated = false;
  CGAL::FT rand_move_points_range = 0.01;

  util::ConfigurationSPtr config = util::Configuration::getInstance();
  if (config->isLoaded()) {
      rand_move_points = config->getBool("main", "rand_move_points");
      rand_move_points_when_degenerated = config->getBool("main", "rand_move_points_when_degenerated");
      rand_move_points_range = config->getFT("main", "rand_move_points_range");
  }

  // unperturbed is not currently supported
  CGAL_warning_msg(rand_move_points_when_degenerated || rand_move_points,
    "Not perturbing...?");

  if (rand_move_points_when_degenerated && !rand_move_points) {
      DEBUG_PRINT("Checking if all combinations of 3 facet supporting planes intersect in a point.");
      DEBUG_PRINT("If this takes too long, disable 'rand_move_points_when_degenerated'.");
      if (!PolyhedronTransformation::doAll3PlanesIntersect(polyhedron)) {
          DEBUG_PRINT("Not all combinations of 3 planes intersect.");
          rand_move_points = true;
      } else {
          DEBUG_PRINT("No need to perturb");
      }
  }

  if (rand_move_points) {
      polyhedron = PolyhedronTransformation::perturb(polyhedron);
  }

  // Failure here is likely from a bad perturbation (after all, there is an epsilon probability
  // that we create something that is degenerate).
  // @todo try again with another perturbation, or smarter (iterative) perturbation
  CGAL_assertion_msg(polyhedron && polyhedron->isConsistent(),
                     "Error: invalid polyhedron (bad perturbation?)");

  return polyhedron;
}







// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //




PolyhedronSPtr PolyhedronTransformation::perturb(PolyhedronSPtr polyhedron) {
  // Check if we can tilt facets' planes (i.e., nudge plane coefficients) directly.
  // A sufficient condition is that all vertices have degree 3: in that case, a small tilt
  // of the plane will still yield a single intersection point.
  // That's not the case (in general) for degree > 3 vertices as there would no longer be
  // a single intersection point for the tilted planes.
  //
  // The advantage is that we then manipulate smaller meshes since the faces are polygonal.
  bool canUsePlaneTilts = true;

  // copy the polyhedron because we will merge (almost) coplanar faces and check if the result
  // is a mesh with only degree 3 vertices.
  PolyhedronSPtr polyhedron_cpy = polyhedron->clone();

  // @todo?
  // could we merge non-connected input faces as to assign them the same (tilted) plane?
  // Often in inputs we have many faces that correspond to the same plane, but vertical faces
  // split it into separate connected components.
  // The important thing is that we don't want to create degenerate conditions so the CCs
  // should NOT interact with each other; how to prevent that?...
  db::_3d::AbstractFile::mergeCoplanarFacets(polyhedron_cpy);
  db::_3d::AbstractFile::removeVerticesDegLt3(polyhedron_cpy);

  db::_3d::OBJFile::save("results/pre-tilt_merged.obj", polyhedron_cpy,
                          false /*do_triangulate*/,
                          true /*convert_to_double*/);

  CGAL_assertion(polyhedron_cpy && polyhedron_cpy->isConsistent());

  bool all_degree_3 = true;

  std::list<VertexSPtr>::iterator it_v = polyhedron_cpy->vertices().begin();
  while (it_v != polyhedron_cpy->vertices().end()) {
    VertexSPtr vertex = *it_v++;
    if (vertex->degree() != 3) {
      DEBUG_PRINT("Can't use plane tilts because of " << vertex->toString());
      all_degree_3 = false;
      break;
    }
  }

  if (canUsePlaneTilts) {
      DEBUG_PRINT("Tilting the polyhedron's facets");

      std::list<FacetSPtr>::iterator it_f = polyhedron_cpy->facets().begin();
      while (it_f != polyhedron_cpy->facets().end()) {
          FacetSPtr facet = *it_f++;
          facet->perturbPlaneCoefficients();
      }

      // If we have forced tilts with degree > 3 vertices, we cannot yet recompute point
      // positions: the incident facets will no longer intersect in a single point,
      if (all_degree_3) {
          resetPoints(polyhedron_cpy);
      }

      polyhedron = polyhedron_cpy;
  } else {
      randMovePoints(polyhedron);
  }

  DEBUG_PRINT("Done with perturbation");

  CGAL_postcondition(polyhedron && polyhedron->isConsistent());

  return polyhedron;
}





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //






void
OutwardMeshOffset::
nudge_points(Mesh& sm)
{
    DEBUG_PRINT("Apply random point perturbations...");
    CGAL_precondition(CGAL::is_triangle_mesh(sm));
    CGAL::Polygon_mesh_processing::random_perturbation(sm, 1e-10, CGAL::parameters::do_project(false).random_seed(0));
}




// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //




bool
OutwardMeshOffset::
merge_coplanar_faces(Mesh& sm)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  remesh almost planar is no good because Surface_mesh cannot represent faces with multiple boundaries

  // run the remeshing algorithm using filled properties
  std::vector<Vector3> normal_map(num_faces(sm));
  for(face_descriptor f : faces(sm)) {
      normal_map[f] = plane_map[f].orthogonal_vector();
  }

  Mesh out;
  PMP::remesh_almost_planar_patches(sm,
                                    out,
                                    nb_regions, nb_corners,
                                    CGAL::make_random_access_property_map(region_ids),
                                    CGAL::make_random_access_property_map(corner_id_map),
                                    CGAL::make_random_access_property_map(ecm),
                                    CGAL::parameters::patch_normal_map(CGAL::make_random_access_property_map(normal_map)),
                                    CGAL::parameters::do_not_triangulate_faces(true));

  // carry over the weight information from 'sm' to 'out'
  auto fwm = sm.property_map<face_descriptor, double>("f:weight");
  CGAL_assertion(bool(fwm));
  auto fwm_out = out.add_property_map<face_descriptor, double>("f:weight").first;

  std::map<std::size_t, double> patch_weights;
  for (face_descriptor f : faces(sm)) {
      std::size_t region_id = region_ids[f];
      std::cout << "face " << f << " in patch " << region_id << " wants to assign its weight " << get(*fwm, f);
      if (patch_weights.count(region_id) != 0) {
          std::cout << ", and patch already has weight " << patch_weights[region_id] << std::endl;

          // "warning" only because we can have numerical errors (note that these are only because
          // of possible approximations in the input: there is nothing approximate with EPECK
          // in the normal computation or weight interpolation).
          CGAL_warning(patch_weights[region_id] == get(*fwm, f));
      } else {
          std::cout << std::endl;
      }
      patch_weights[region_id] = get(*fwm, f);
  }

  for (face_descriptor f : faces(out)) {
      std::size_t region_id = region_ids[f];
      CGAL_assertion(region_id < patch_weights.size());
      put(fwm_out, f, patch_weights[region_id]);
  }

  sm = out;

  return is_valid(sm) && CGAL::is_valid_face_graph(sm);
}





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //




static bool chamfer_vertices(const std::vector<vertex_descriptor>& vs, Mesh& sm);
static bool chamfer_high_degree_vertices(Mesh& sm);

bool
OutwardMeshOffset::
chamfer_vertices(const std::vector<vertex_descriptor>& vertices_to_chamfer,
                 Mesh& sm)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;

  DEBUG_PRINT("Chamfer " << vertices_to_chamfer.size() << " vertices with high degree...");

  auto fwm = sm.property_map<face_descriptor, double>("f:weight");
  CGAL_assertion(bool(fwm));
  double min_weight = std::numeric_limits<double>::max(); // 'double' on purpose
  for(face_descriptor f : faces(sm)) {
      min_weight = (std::min)(min_weight, get(*fwm, f));
  }

  // the weight of chamfering faces is very small, such that they do not live long
  // @fixme chamfering at a vertex that has a face with weight 0...?
  // @fixme weight being typed double
  // @fixme currently the min_weight at that point is that of the bbox faces, which is already
  // 1e-10 times the min of zero weight faces, which is already 1e-10 the "real" min weight
  double chamfer_weight = 1e-10 * min_weight;

  // For all vertices with degree > 3, compute the minimal distance to the nearest neighboring point
  // and construct a regular tetrahedron centered at the vertex and inscribed in a sphere 1000
  // times smaller than the distance to the nearest neighboring point
  std::vector<Tetrahedron3> tetrahedra;

  for (vertex_descriptor vertex : vertices_to_chamfer) {
      DEBUG_PRINT("Chamfer " << vertex << "...");

      // Compute the distance to the nearest neighboring point
      CGAL::FT min_sq_distance = std::numeric_limits<double>::max();
      for (vertex_descriptor nv : CGAL::vertices_around_target(halfedge(vertex, sm), sm)) {
          CGAL::FT sq_distance = CGAL::squared_distance(sm.point(vertex), sm.point(nv));
          if (sq_distance < min_sq_distance) {
              min_sq_distance = sq_distance;
          }
      }

      // Create the tetrahedron centered at the vertex
      CGAL::FT radius = CGAL::approximate_sqrt(min_sq_distance) / 10; // @tmp put something much smaller

      auto random_unit_axis = []() {
          static std::random_device rd;
          unsigned int s = 0; // rd()
          static std::mt19937 gen(s);
          static std::uniform_real_distribution<> angle_dist(0, 2 * CGAL_PI);
          static std::uniform_real_distribution<> axis_dist(-1, 1);
          double angle = angle_dist(gen);
          double z = axis_dist(gen);
          double x = std::sqrt(1 - z * z) * std::cos(angle);
          double y = std::sqrt(1 - z * z) * std::sin(angle);
          return Vector3(x, y, z);
      };

      // Rotate the tetrahedron around a random axis
      Vector3 first_axis = random_unit_axis();
      Vector3 second_axis = random_unit_axis();
      second_axis = second_axis - CGAL::scalar_product(first_axis, second_axis) * first_axis;
      second_axis = second_axis / CGAL::approximate_sqrt(second_axis.squared_length()); // handle null length (generate another 2nd)
      Vector3 third_axis = CGAL::cross_product(first_axis, second_axis);

      // std::cout << "rotation matrix" << std::endl;
      // std::cout << first_axis << " N = " << first_axis.squared_length() << std::endl;
      // std::cout << second_axis << " N = " << second_axis.squared_length() << std::endl;
      // std::cout << third_axis << " N = " << third_axis.squared_length() << std::endl;

      // the rotation matrix
      CGAL::Aff_transformation_3<CGAL::K> rotation_matrix(
          first_axis.x(), first_axis.y(), first_axis.z(),
          second_axis.x(), second_axis.y(), second_axis.z(),
          third_axis.x(), third_axis.y(), third_axis.z());

      // the rotated tetrahedron vertices
      Vector3 v0 = rotation_matrix.transform(Vector3(1, 1, 1));
      Vector3 v1 = rotation_matrix.transform(Vector3(-1, -1, 1));
      Vector3 v2 = rotation_matrix.transform(Vector3(-1, 1, -1));
      Vector3 v3 = rotation_matrix.transform(Vector3(1, -1, -1));

      // std::cout << "v0 = " << v0 << std::endl;
      // std::cout << "v1 = " << v1 << std::endl;
      // std::cout << "v2 = " << v1 << std::endl;
      // std::cout << "v3 = " << v2 << std::endl;

      tetrahedra.emplace_back(sm.point(vertex) + radius * v0,
                              sm.point(vertex) + radius * v1,
                              sm.point(vertex) + radius * v2,
                              sm.point(vertex) + radius * v3);
  }

  // Boolean operation to remove all tetrahedra from the mesh
  Mesh chamfered_mesh;
  Mesh tetrahedra_mesh;

  for (const auto& tetrahedron : tetrahedra) {
      Mesh tetrahedron_mesh;
      CGAL::make_tetrahedron(tetrahedron[0], tetrahedron[1], tetrahedron[2], tetrahedron[3],
                             tetrahedron_mesh);
      CGAL_assertion(!PMP::does_self_intersect(tetrahedron_mesh));

      static int tet_id = -1;
      CGAL::IO::write_polygon_mesh("results/tet_" + std::to_string(++tet_id) + ".obj", tetrahedron_mesh,
                                   CGAL::parameters::stream_precision(17));

      CGAL::copy_face_graph(tetrahedron_mesh, tetrahedra_mesh);
  }

  CGAL_assertion(!PMP::does_self_intersect(tetrahedra_mesh));

  auto tet_fwm = tetrahedra_mesh.add_property_map<face_descriptor, double>("f:weight").first;
  for (face_descriptor f : faces(tetrahedra_mesh)) {
      put(tet_fwm, f, chamfer_weight);
  }

  Mesh result;
  auto res_fwm = result.add_property_map<face_descriptor, double>("f:weight").first;

  // the visitor keeps the weight from faces coming from 'sm', and sets a very small weight
  // for faces coming from the tetrahedron. The point is that the faces of the tetrahedron
  // should disappear very quickly, "eaten" by other faces that shift at much higher speed
  // such that in the result, we do not see that we had chamfered the vertices
  struct Weight_setter_visitor
    : public CGAL::Polygon_mesh_processing::Corefinement::Default_visitor<Mesh>
  {
      // @fixme surface mesh
      boost::container::flat_map<const Mesh*, Mesh::Property_map<Mesh::Face_index, double> > properties;
      double weight = -1;

      Weight_setter_visitor() { properties.reserve(3); }

      void before_subface_creations(face_descriptor f_split, Mesh& tm)
      {
        weight = properties[&tm][f_split];
      }

      void after_subface_created(face_descriptor f_new, Mesh& tm)
      {
        properties[&tm][f_new] = weight;
      }

      void after_face_copy(face_descriptor f_src, Mesh& tm_src,
                            face_descriptor f_tgt, Mesh& tm_tgt)
      {
        properties[&tm_tgt][f_tgt] = properties[&tm_src][f_src];
      }
  };

  Weight_setter_visitor visitor;
  visitor.properties[&sm] = *fwm;
  visitor.properties[&tetrahedra_mesh] = tet_fwm;
  visitor.properties[&result] = res_fwm;

  CGAL_assertion(PMP::does_bound_a_volume(sm));
  CGAL_assertion(PMP::does_bound_a_volume(tetrahedra_mesh));

  PMP::corefine_and_compute_difference(sm, tetrahedra_mesh, result,
                                        CGAL::parameters::visitor(visitor));

  utils::save_colored_mesh(result, res_fwm, "results/coref_regions.ply");

  sm = result;
  fwm = sm.property_map<face_descriptor, double>("f:weight");

  CGAL::IO::write_polygon_mesh("results/chamfered.obj", sm, CGAL::parameters::stream_precision(17));

  // test if the chamfered mesh has only degree 3 vertices
  // preprocess should return the PolyhedronSPtr
  // most complicated part: simplify the mesh according to the shape detection regions
  //   need some kind of SM_edge SS3Polyhedron_edge correspodence (build it in the load())
  // weights should be CGAL::FT
  // performance acute merging -> write the PLY
  // write the code to get rid of main.cpp (inward offsetting)

  return true;
}

bool
OutwardMeshOffset::
chamfer_high_degree_vertices(Mesh& sm)
{
    namespace PMP = CGAL::Polygon_mesh_processing;

    DEBUG_PRINT("Chamfer mesh...");

    CGAL::IO::write_polygon_mesh("results/pre-perturbation.obj", sm, CGAL::parameters::stream_precision(17));
    CGAL_assertion(CGAL::is_triangle_mesh(sm));

    namespace PMP = CGAL::Polygon_mesh_processing;

    DEBUG_PRINT("Attempting to tilt the polyhedron's facets...");

    // Check if we can tilt facets' planes (i.e., nudge plane coefficients) directly.
    // A sufficient condition is that all vertices have degree 3: in that case, a small tilt
    // of the plane will still yield a single intersection point.
    // That's not the case (in general) for degree > 3 vertices as there would no longer be
    // a single intersection point for the tilted planes.
    //
    // The advantage is that we then manipulate smaller meshes since the faces are polygonal.
    bool can_use_plane_tilts = true;

    // @todo?
    // could we merge non-connected input faces as to assign them the same (tilted) plane?
    // Often in inputs we have many faces that correspond to the same plane, but vertical faces
    // split it into separate connected components.
    // The important thing is that we don't want to create degenerate conditions so the CCs
    // should NOT interact with each other; how to prevent that?...

    // declare vectors to store mesh properties
    std::vector<std::size_t> region_ids(num_faces(sm));
    boost::vector_property_map<Plane3> plane_map; // supporting planes of the regions detected

    // detect planar regions in the mesh
    // @fixme growing should stop if it merges faces with different weights
    // and give an error for true coplanar faces with different weights
    std::size_t nb_regions =
        PMP::region_growing_of_planes_on_faces(sm,
                                               CGAL::make_random_access_property_map(region_ids),
                                               CGAL::parameters::cosine_of_maximum_angle(0.98)
                                                                .region_primitive_map(plane_map)
                                                                .maximum_distance(0.001));

    utils::save_colored_mesh(sm, region_ids, "results/regions.ply");

    // detect corner vertices on the boundary of planar regions
    std::vector<std::size_t> corner_id_map(num_vertices(sm), -1); // corner status of vertices
    std::vector<bool> ecm(num_edges(sm), false); // mark edges at the boundary of regions

    PMP::detect_corners_of_regions(sm,
                                   CGAL::make_random_access_property_map(region_ids),
                                   nb_regions,
                                   CGAL::make_random_access_property_map(corner_id_map),
                                   CGAL::parameters::cosine_of_maximum_angle(0.98).
                                                     maximum_distance(0.001).
                                                     edge_is_constrained_map(CGAL::make_random_access_property_map(ecm)));

    std::vector<vertex_descriptor> vertices_to_chamfer;

    vertex_iterator vit = vertices(sm).begin(), vend = vertices(sm).end();
    std::size_t cid = 0;
    for (; vit!=vend; ++vit, ++cid) {
        std::set<std::size_t> incident_regions;
        for (face_descriptor f : CGAL::faces_around_target(halfedge(*vit, sm), sm)) {
            auto res = incident_regions.insert(region_ids[f]);
            if (res.second && incident_regions.size() > 3) {
                DEBUG_PRINT("Region corner with degree > 3: " << sm.point(*vit));
                can_use_plane_tilts = false;
                vertices_to_chamfer.push_back(*vit);
                break;
            }
        }
    }

    // @todo with incident face _IDS_ instead of incident faces,
    // maybe we can avoid having to triangulate EVERYTHING as soon as one vertex
    // does not have degree 3 by marking the vertex that needs to be moved manually
    // and tilts are possible for the incident faces as long as they have less than 3 vertices
    // that need to be manually moved (if 3 out of n vertices must be moved manually, we can still
    // tilt the face such that it matches the 3 moved vertices and other vertices will move
    // with the tilted plane)
    if (!can_use_plane_tilts) {
        can_use_plane_tilts = chamfer_vertices(vertices_to_chamfer, sm);
        CGAL_assertion(can_use_plane_tilts);
    }

    return is_valid(sm) && is_valid_face_graph(sm) && can_use_plane_tilts;
}



PolyhedronSPtr
OutwardMeshOffset::
preprocess(Mesh& sm)
{
    DEBUG_PRINT("Preprocessing mesh...");

    CGAL_precondition(CGAL::is_triangle_mesh(sm));

#ifdef CGAL_SS3_APPLY_CHAMFERING
    bool apply_perturbation = false;
    bool apply_perturbation_when_degenerated = false;
    CGAL::FT rand_move_points_range = 0.01;

    util::ConfigurationSPtr config = util::Configuration::getInstance();
    if (config->isLoaded()) {
        apply_perturbation = config->getBool("main", "rand_move_points");
        apply_perturbation_when_degenerated = config->getBool("main", "rand_move_points_when_degenerated");
        rand_move_points_range = config->getFT("main", "rand_move_points_range");
    }

    // unperturbed is not currently supported
    CGAL_warning_msg(apply_perturbation_when_degenerated || apply_perturbation,
      "Not perturbing...?");

    if (apply_perturbation_when_degenerated && !apply_perturbation) {
        DEBUG_PRINT("Checking if all combinations of 3 facet supporting planes intersect in a point.");
        DEBUG_PRINT("If this takes too long, disable 'rand_move_points_when_degenerated'.");
        // @todo implement a smarter doAll3PlanesIntersect() for Surface_mesh, careful about normalization of planes
        if (true) {
            DEBUG_PRINT("Not all combinations of 3 planes intersect.");
            apply_perturbation = true;
          } else {
            DEBUG_PRINT("No need to perturb");
          }
        }

        if(apply_perturbation) {
        // as to force every plane to have degree 3
        // @todo make this a config option
        chamfer_high_degree_vertices(sm);
    }
#endif

    PolyhedronSPtr polyhedron = convert(sm);
    CGAL_assertion(polyhedron && polyhedron->isConsistent());

    polyhedron = PolyhedronTransformation::perturb(polyhedron);
    CGAL_assertion(polyhedron && polyhedron->isConsistent());

    if (!polyhedron || !polyhedron->isConsistent()) {
        // Failure here is likely from a bad perturbation (after all, there is a probabily
        // epsilon that we create something that is degenerate).
        // @todo try again with another perturbation, or smarter (iterative) perturbation
        std::cerr << "Error: invalid polyhedron (bad perturbation?)" << std::endl;
        return {};
    }

    return polyhedron;
}




// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //





VertexSPtr Polyhedron::ith_vertex(const std::size_t i)
{
    CGAL_assertion(i < this->vertices_.size());
    return *(std::next(std::cbegin(this->vertices_), i));
}
EdgeSPtr Polyhedron::ith_edge(const std::size_t i)
{
    CGAL_assertion(i < this->edges_.size());
    return *(std::next(std::cbegin(this->edges_), i));
}
FacetSPtr Polyhedron::ith_facet(const std::size_t i)
{
    CGAL_assertion(i < this->facets_.size());
    return *(std::next(std::cbegin(this->facets_), i));
}





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



useless because the quintessence of vanish events is that some edge length(s) go to 0

struct Stop_on_long_enough_edge
    : public Base_mesh_offset_visitor
{
    Stop_on_long_enough_edge(CGAL::FT min_edge_length) : sq_min_edge_length_(CGAL::square(min_edge_length)) { }

    bool go_further(int step_id, PolyhedronSPtr polyhedron, CGAL::FT offset) override {
        for (EdgeSPtr edge : polyhedron->edges()) {
            CGAL::FT sq_l = CGAL::squared_distance(*(edge->getVertexSrc()->getPoint()),
                                                   *(edge->getVertexDst()->getPoint()));
            if (sq_l < sq_min_edge_length_) {
                std::cout << edge->getID() << " with length " << sq_l << std::endl;
                return true;
            }
        }

        db::_3d::OBJFile::save("results/intermediate-exact.obj", polyhedron,
                                true /*triangulate*/,
                                false /*do not convert to double*/);
        return false;
    }

    void before_offset_event(PolyhedronSPtr polyhedron, CGAL::FT offset) override { }

    void on_save_offset_event(PolyhedronSPtr polyhedron, CGAL::FT offset) override { }

    void after_offset_event(PolyhedronSPtr polyhedron, CGAL::FT offset) override { }

private:
    CGAL::FT sq_min_edge_length_;
};






// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //





#if 1
        DEBUG_PRINT("Running first pass");

        SimpleStraightSkelSPtr algoskel3d_f =
            SimpleStraightSkel::create(polyhedron, controller, save_offsets, save_path);

        // @fixme this doesn't guarantee that all edges are big enough to have stability through remeshing...
        Mesh tmp_sm;
        utils::Stop_on_far_enough_event visitor(1e-1);
        algoskel3d_f->setVisitor(&visitor);

        // @todo handle the (unlikely) case where run() handled the whole offsetting process
        try {
          algoskel3d_f->run();
        } catch(utils::Far_enough_event) {
          std::cout << "caught the throw" << std::endl;
          bool success = db::_3d::Surface_meshIO::save(visitor.polyhedron_, tmp_sm,
                                                       true /*triangulate*/, false /*no doubles*/);
          CGAL_assertion(success);
        }
#else
        // This one isn't great because we don't know when the first event is, and it can be very early

        // @todo give that function a nicer name, here we are just interested in the split
        // @todo use the splitter from the config file
        SimpleStraightSkel::init(polyhedron, algo::_3d::ConvexVertexSplitter::create());

        // We move forward in time by epsilon
        // @fixme, technically we should guarantee that we are below the first event
        //         - we could scan and do a const event with lower time
        //         - when scanning, can do a linear pass
        //         - the shift must also be below the earliest save offset value
        //
        // @fixme shift the save offsets by the shift value, or initialize the shiftting algorithm at "shift"
        CGAL::FT shift = -1e-10;
        PolyhedronTransformation::shiftFacetsInPlace(polyhedron, shift);
        db::_3d::OBJFile::save("results/intermediate-polyhedron.obj", polyhedron, false /*do not triangulate*/);

        Mesh tmp_sm;
        bool success = db::_3d::Surface_meshIO::save(polyhedron, tmp_sm,
                                                     true /*triangulate*/,
                                                     false /*do not convert to double*/);
        CGAL_assertion(success);
#endif





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //






PolyhedronSPtr
OutwardMeshOffset::
convert(Mesh& sm,
        const bool force_simplification)
{
    namespace PMP = CGAL::Polygon_mesh_processing;

    DEBUG_PRINT("Converting mesh...");

    bool merge_faces = false;

    util::ConfigurationSPtr config = util::Configuration::getInstance();
    std::string section("main");
    if (config->isLoaded() &&
        config->contains(section, "merge_coplanar_faces") &&
        config->getBool(section, "merge_coplanar_faces")) {
        merge_faces = true;
    }

    if (!merge_faces) {
        return db::_3d::Surface_meshIO::load(sm);
    }

    CGAL::Bbox_3 bbox = PMP::bbox(sm);
    const CGAL::FT diag_length = CGAL::approximate_sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                                        CGAL::square(bbox.ymax() - bbox.ymin()) +
                                                        CGAL::square(bbox.zmax() - bbox.zmin()));

    // Use shape detection to analyze the mesh
    std::vector<std::size_t> region_ids(num_faces(sm));
    boost::vector_property_map<Plane3> plane_map; // supporting planes of the regions detected

    const CGAL::FT cos_of_max_angle = 0.98;
    const CGAL::FT max_distance = 0.0001 * diag_length;

    // detect planar regions in the mesh
    // @todo growing should:
    // - use the .ini value of 'epsilon_coplanarity'
    // - stop if it merges faces with different weights
    // - give an error for adjacent coplanar faces that have different weights
    std::size_t nb_regions =
        PMP::region_growing_of_planes_on_faces(sm,
                                               CGAL::make_random_access_property_map(region_ids),
                                               CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle)
                                                                .region_primitive_map(plane_map)
                                                                .maximum_distance(max_distance));

    static int region_dump_id = -1;
    utils::save_colored_mesh(sm, region_ids, "results/regions_" + std::to_string(++region_dump_id) + ".ply");

    // detect corner vertices on the boundary of planar regions
    std::vector<std::size_t> corner_id_map(num_vertices(sm), -1); // corner status of vertices
    std::vector<bool> ecm(num_edges(sm), false); // mark edges at the boundary of regions

    PMP::detect_corners_of_regions(sm,
                                   CGAL::make_random_access_property_map(region_ids),
                                   nb_regions,
                                   CGAL::make_random_access_property_map(corner_id_map),
                                   CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle).
                                                     maximum_distance(max_distance).
                                                     edge_is_constrained_map(CGAL::make_random_access_property_map(ecm)));

    // @debug
    {
        for (face_descriptor f : faces(sm)) {
            std::cout << "face " << f << " is in region " << region_ids[f] << std::endl;
        }
    }

    // the almost-coplanar merge is performed after the conversion to the Polyhedron
    // data structure because we want to be able to create faces that have holes,
    // which the CGAL::Surface_mesh class does not support
    std::map<edge_descriptor, EdgeWPtr> e2e;
    PolyhedronSPtr polyhedron = db::_3d::Surface_meshIO::load(sm, {}, e2e);

    PolyhedronTransformation::normalizeFacetPlanes(polyhedron);

    // std::cout << "Converted polyhedron (out): " << polyhedron->toString() << std::endl;
    db::_3d::OBJFile::save("results/convert-base_polyhedron.obj", polyhedron, false /*do not triangulate*/);

    // If everything is degree 3 in the region growing, merge the facets
    bool all_degree_3 = true;
    vertex_iterator vit = vertices(sm).begin(), vend = vertices(sm).end();
    std::size_t cid = 0;
    for (; vit!=vend; ++vit, ++cid) {
        std::set<std::size_t> incident_regions;
        for (face_descriptor f : CGAL::faces_around_target(halfedge(*vit, sm), sm)) {
            // could break early but it's useful to know how many high-degree corners were detected
            incident_regions.insert(region_ids[f]);
        }

        if (incident_regions.size() > 3) {
            DEBUG_PRINT("Region corner with degree " << incident_regions.size() << " at " << sm.point(*vit));
            all_degree_3 = false;
            break;
        }
    }

    DEBUG_PRINT("all_degree_3 = " << all_degree_3);
    if (all_degree_3 || force_simplification) {
        // merge the facets incident to an unconstrained edge (i.e., the edge is interior to a region)
        for (edge_descriptor e: edges(sm)) {
            if (ecm[e]) {
                continue;
            }

            EdgeSPtr edge = e2e[e].lock();
            if (!edge) {
                continue;
            }

            DEBUG_PRINT("Merging facets " << edge->getFacetL()->getID() << " and " << edge->getFacetR()->getID());
            CGAL_assertion(sm.point(source(e, sm)) == *(edge->getVertexSrc()->getPoint()));
            CGAL_assertion(sm.point(target(e, sm)) == *(edge->getVertexDst()->getPoint()));

            // @todo it seems like intermediate states are somewhat unsound during edge merging
            db::_3d::AbstractFile::mergeFacets(edge, polyhedron);
        }
    }

    db::_3d::OBJFile::save("results/convert-before_sanitize.obj", polyhedron, false /*do not triangulate*/);

    DEBUG_PRINT("Sanitizing..");
    db::_3d::AbstractFile::sanitize(polyhedron);

    if (all_degree_3 || force_simplification) {
        bool success = PolyhedronTransformation::resetPoints(polyhedron);
        if (!success) {
            std::cerr << "Warning: invalid polyhedron" << std::endl;
            return {};
        }
    }

    std::cout << "Converted, " << polyhedron->facets().size() << " facets" << std::endl;
    db::_3d::OBJFile::save("results/convert-final.obj", polyhedron, false /*do not triangulate*/);

    return polyhedron;
}




// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //




// dump intermediate offsets every -0.01
for (double i = -0.01; i > -1.0; i -= 0.01) {
    PolyhedronSPtr poly_c_offset_i = PolyhedronTransformation::shiftFacets(poly_c, i);
    if (poly_c_offset_i) {
        db::_3d::OBJFile::save("results/split_" + std::to_string(test_id) +
                                "_offset_" + std::to_string(i) + ".obj",
                                poly_c_offset_i, false);
    }
}








// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //





#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>


// this is bad because the constraints should be that planes incident to the same vertex
// intersect in a single point, not that planes go through the perturbed positions
void PolyhedronTransformation::randTiltPlanesSystem(PolyhedronSPtr polyhedron) {

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        facet->perturbPlaneCoefficients();
    }

    resetPoints(polyhedron);

    using Program = CGAL::Quadratic_program<CGAL::FT>;
    using Solution = CGAL::Quadratic_program_solution<CGAL::FT>;

    Program qp(CGAL::EQUAL, true, -1, true, 1); // Minimize, no lower/upper bounds

    // D & c
    int pl_idx = 0;
    CGAL::FT c0 = 0;
    for (FacetSPtr facet : polyhedron->facets()) {
        // Minimize (x - x0)^2 for a,b,c,d
        qp.set_d(pl_idx + 0, pl_idx + 0, 2.);
        qp.set_d(pl_idx + 1, pl_idx + 1, 2.);
        qp.set_d(pl_idx + 2, pl_idx + 2, 2.);
        qp.set_d(pl_idx + 3, pl_idx + 3, 2.);
        qp.set_c(pl_idx + 0, - 2. * facet->getPlane()->a());
        qp.set_c(pl_idx + 1, - 2. * facet->getPlane()->b());
        qp.set_c(pl_idx + 2, - 2. * facet->getPlane()->c());
        qp.set_c(pl_idx + 3, - 2. * facet->getPlane()->d());
        c0 += CGAL::square(facet->getPlane()->a());
        c0 += CGAL::square(facet->getPlane()->b());
        c0 += CGAL::square(facet->getPlane()->c());
        c0 += CGAL::square(facet->getPlane()->d());
        pl_idx += 4;
    }
    qp.set_c0(c0);

    // A & b (planarity constraints)
    int cst_idx = 0;
    pl_idx = 0;
    for (FacetSPtr facet : polyhedron->facets()) {
        for (VertexSPtr vertex : facet->vertices()) {
            CGAL_assertion(facet->getPlane()->has_on(*(vertex->getPoint())));
            qp.set_a(cst_idx, pl_idx + 0, vertex->getPoint()->x());
            qp.set_a(cst_idx, pl_idx + 1, vertex->getPoint()->y());
            qp.set_a(cst_idx, pl_idx + 2, vertex->getPoint()->z());
            qp.set_a(cst_idx, pl_idx + 3, 1.0);
            qp.set_b(cst_idx, 0.0); // Planarity constraint
            ++cst_idx;
        }
        pl_idx += 4;
    }

    // Solve the quadratic program
    Solution solution = CGAL::solve_quadratic_program(qp, CGAL::FT());
    if (solution.is_optimal()) {
        std::cout << "yay, solution is:\n" << solution << std::endl;

        Solution::Variable_value_iterator s_it = solution.variable_values_begin();
        std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
        while (it_f != polyhedron->facets().end()) {
            FacetSPtr facet = *it_f++;
            std::cout << "  From coefficients ["
              << facet->getPlane()->a() << " " << facet->getPlane()->b() << " "
              << facet->getPlane()->c() << " " << facet->getPlane()->d() << "]" << std::endl;
            CGAL::FT na = CGAL::quotient_truncation(*s_it++);
            CGAL::FT nb = CGAL::quotient_truncation(*s_it++);
            CGAL::FT nc = CGAL::quotient_truncation(*s_it++);
            CGAL::FT nd = CGAL::quotient_truncation(*s_it++);
            facet->setPlane(KernelFactory::createPlane3(na, nb, nc, nd));
            std::cout << "  To coefficients ["
              << facet->getPlane()->a() << " " << facet->getPlane()->b() << " "
              << facet->getPlane()->c() << " " << facet->getPlane()->d() << "]" << std::endl;
        }
    } else {
        std::cerr << "Optimization failed!" << std::endl;
    }

    std::exit(1);
}





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //





// 0 fixed points: nudge all coefficients independently
if (fixed_points.size() == 0) {
    perturbPlaneCoefficientsNudge(range);
} else if (fixed_points.size() == 1) {
    // 1 fixed point: nudge (a, b, c), recompute d so the plane passes through the point
    const Point3& point = *(fixed_points[0]);

    CGAL::FT na = nudge(plane_->a());
    CGAL::FT nb = nudge(plane_->b());
    CGAL::FT nc = nudge(plane_->c());
    const CGAL::FT& x0 = point.x();
    const CGAL::FT& y0 = point.y();
    const CGAL::FT& z0 = point.z();
    CGAL::FT nd = - (na * x0 + nb * y0 + nc * z0);
    plane_ = KernelFactory::createPlane3(na, nb, nc, nd);
    CGAL_postcondition(plane_->has_on(point));
} else if (fixed_points.size() == 2) {
    // 2 fixed points: nudge (a, b, c) with constraint that plane passes through both points
    // Nudge a and b, solve for c and d so plane passes through both points
    const Point3& p0 = *(fixed_points[0]);
    const Point3& p1 = *(fixed_points[1]);

    // The plane goes through p0 and p1:
    // na*x0 + nb*y0 + nc*z0 + nd = 0
    // na*x1 + nb*y1 + nc*z1 + nd = 0
    // => na*(x1-x0) + nb*(y1-y0) = (z0-z1)*nc
    // => nc = -(na*(x0-x1) + nb*(y0-y1)) / (z0-z1)
    const CGAL::FT& a0 = plane_->a();
    const CGAL::FT& b0 = plane_->b();
    const CGAL::FT& c0 = plane_->c();
    const CGAL::FT& d0 = plane_->d();
    const CGAL::FT& x0 = p0.x();
    const CGAL::FT& y0 = p0.y();
    const CGAL::FT& z0 = p0.z();
    const CGAL::FT& x1 = p1.x();
    const CGAL::FT& y1 = p1.y();
    const CGAL::FT& z1 = p1.z();

    CGAL::FT na, nb, nc, nd;

    CGAL::FT dz = z0 - z1;
    if (!CGAL::is_zero(dz)) {
        na = nudge(a0);
        nb = nudge(b0);
        nc = - (na*(x0-x1) + nb*(y0-y1)) / dz;
        nd = - (na*x0 + nb*y0 + nc*z0);
    } else {
        // If dz == 0, try to use nb & nd as variables
        CGAL::FT dy = y0 - y1;
        if (!CGAL::is_zero(dy)) {
            na = nudge(a0);
            nc = nudge(c0);
            nb = -(na*(x0-x1) + nc*(z0-z1)) / dy;
            nd = - (na*x0 + nb*y0 + nc*z0);
        } else {
            // If both dz and dy are zero, try na & nd as variables...
            CGAL::FT dx = x0 - x1;
            if (!CGAL::is_zero(dx)) {
                nb = nudge(b0);
                nc = nudge(c0);
                na = -(nb*(y0-y1) + nc*(z0-z1)) / dx;
                nd = - (na*x0 + nb*y0 + nc*z0);
            } else {
                // Points are identical, fallback to 1-point case
                na = nudge(a0);
                nb = nudge(b0);
                nc = nudge(c0);
                nd = - (na*x0 + nb*y0 + nc*z0);
            }
        }
    }
    plane_ = KernelFactory::createPlane3(na, nb, nc, nd);

    // If we have 2, at least one of `a`, `b`, or `c` coefficients are affected
    // and could change significantly
    // normalizePlaneCoefficients();

    CGAL_postcondition(plane_->has_on(p0));
    CGAL_postcondition(plane_->has_on(p1));
} else {
    DEBUG_PRINT("Error: called fixed point facet perturbation with > 2 fixed points");
}





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //





void PolyhedronTransformation::iterativeTiltPlanes(PolyhedronSPtr polyhedron) {
    // Assume polyhedron is already simplified and faces are polygonal.

    std::cout << "Iterative tilting of planes" << std::endl;

    enum class FacetStatus { FREE = 0, NUDGED, OVER_CONSTRAINED };
    int max_facet_id = -1;
    for (auto it = polyhedron->facets().begin(); it != polyhedron->facets().end(); ++it) {
        if ((*it)->getID() > max_facet_id) {
            max_facet_id = (*it)->getID();
        }
    }
    std::vector<FacetStatus> facet_status(max_facet_id + 1, FacetStatus::FREE);

    // Helper: for a vertex, count incident nudged faces
    auto vertex_nudged_face_count = [&](VertexSPtr v) {
        int n = 0;
        for (auto it = v->facets().begin(); it != v->facets().end(); ++it) {
            FacetWPtr fw = *it;
            if (!fw.expired()) {
                FacetSPtr f = FacetSPtr(fw);
                if (facet_status[f->getID()] == FacetStatus::NUDGED)
                    ++n;
            }
        }
        return n;
    };

    // Main loop: iteratively nudge faces with the largest polygonal face incident to a high-degree vertex
    for (;;) {
        VertexSPtr v_max;
        FacetSPtr f_max;
        int max_f_degree = 0;
        for (auto vit = polyhedron->vertices().begin(); vit != polyhedron->vertices().end(); ++vit) {
            VertexSPtr v = *vit;
            if (v->degree() <= 3) {
                continue;
            }
            for (auto fit = v->facets().begin(); fit != v->facets().end(); ++fit) {
                FacetWPtr fw = *fit;
                if (!fw.expired()) {
                    FacetSPtr f = FacetSPtr(fw);
                    if (facet_status[f->getID()] == FacetStatus::FREE) {
                        int deg = f->vertices().size();
                        if (deg > max_f_degree) {
                            max_f_degree = deg;
                            f_max = f;
                            v_max = v;
                        }
                    }
                }
            }
        }

        if (!f_max) {
            break;
        }

        std::cout << "  Processing face #" << f_max->getID() << " with " << f_max->vertices().size()
                  << " vertices" << std::endl;

        // Check if face is a triangle (not interesting)
        if (f_max->vertices().size() == 3) {
            facet_status[f_max->getID()] = FacetStatus::NUDGED;
            continue;
        }

        // Count high-degree vertices with >= 2 incident nudged faces
        std::vector<Point3SPtr> pivots;
        for (auto it = f_max->vertices().begin(); it != f_max->vertices().end(); ++it) {
            VertexSPtr v = *it;
            if (v->degree() > 3) {
                int fc = vertex_nudged_face_count(v);
                if (fc >= 2) {
                    pivots.push_back(v->getPoint());
                }
            }
        }

        if (pivots.size() > 3) {
            facet_status[f_max->getID()] = FacetStatus::OVER_CONSTRAINED;
        } else {
            facet_status[f_max->getID()] = FacetStatus::NUDGED;

            // dump the face into a single OFF file
            std::map<Point3, std::size_t> point_ids;
            std::vector<Point3> points;
            std::vector<std::vector<std::size_t> > polygons;
            triangulate_facet_with_CDT2(f_max, CDT2_Filtering::ODD_EVEN, point_ids, points, polygons);

            static int nudged_face_id = -1;
            CGAL::IO::write_OFF("results/nudged_face_" + std::to_string(++nudged_face_id) + ".OFF", points, polygons);

            f_max->perturbPlaneCoefficientsFixedPoints(1e-10, pivots);

            // reset the point positions of all vertices that are incident to 3 nudged facets
            // (do it now because some faces will be perturbed with constraints, and we need
            // the constraints to be correct geometrically)
            for (auto it_v = f_max->vertices().begin(); it_v != f_max->vertices().end(); ++it_v) {
                VertexSPtr v = *it_v;
                std::vector<Plane3SPtr> planes;
                for (auto it = v->facets().begin(); it != v->facets().end(); ++it) {
                    FacetWPtr fw = *it;
                    if (!fw.expired()) {
                        FacetSPtr f = FacetSPtr(fw);
                        if (facet_status[f->getID()] == FacetStatus::NUDGED) {
                           std::cout << "    " << f->getID() << " is nudged for vertex #" << v->getID() << std::endl;
                           planes.push_back(f->getPlane());
                        }
                    }
                }

                std::cout << planes.size() << " incident nudged facets for vertex #" << v->getID() << " @ " << *(v->getPoint()) << std::endl;
                CGAL_assertion(planes.size() <= 3);
                if (planes.size() == 3) {
                    resetPoint(v, { planes[0], planes[1], planes[2] });
                }
            }
        }
    }

    std::cout << "Statuses:" << std::endl;
    unsigned int nudged_faces_n = 0, oc_faces_n = 0, unexplored_n = 0;
    for (auto it_f = polyhedron->facets().begin(); it_f != polyhedron->facets().end(); ++it_f) {
        FacetSPtr f = *it_f;
        if (facet_status[f->getID()] == FacetStatus::NUDGED) {
            ++nudged_faces_n;
        } else if(facet_status[f->getID()] == FacetStatus::OVER_CONSTRAINED) {
            ++oc_faces_n;
        } else {
            ++unexplored_n;
        }
    }

    std::cout << nudged_faces_n << " nudged faces" << std::endl;
    std::cout << oc_faces_n << " over constrained faces" << std::endl;
    std::cout << unexplored_n << " unexplored faces" << std::endl;

    // nudge everything else
    for (auto it = polyhedron->facets().begin(); it != polyhedron->facets().end(); ++it) {
        FacetSPtr f = *it;
        if (facet_status[f->getID()] == FacetStatus::FREE) {
            f->perturbPlaneCoefficients();
        }
    }
}







// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //






void PolyhedronTransformation::randTiltPlanesv3(PolyhedronSPtr polyhedron) {
    // Identify high-degree vertices
    // std::set<VertexSPtr> high_degree_vertices;
    // for (VertexSPtr v : polyhedron->vertices()) {
    //     if (v->degree() > 3) {
    //         high_degree_vertices.insert(v);
    //     }
    // }

    // high degree vertex --> incident faces used to compute the vertex position
    std::map<VertexSPtr, std::vector<FacetSPtr>> canonical_faces;

    // Sort by size as to avoid triangulating the largest faces
    std::list<FacetSPtr> facets_to_process;
    for (FacetSPtr facet : polyhedron->facets()) {
        facets_to_process.push_back(facet);
    }
    facets_to_process.sort([](const FacetSPtr& a, const FacetSPtr& b) {
        // Sort by number of vertices (descending)
        return a->vertices().size() > b->vertices().size();
    });

    while (!facets_to_process.empty()) {
        FacetSPtr facet = facets_to_process.front();
        facets_to_process.pop_front();

        // Count high-degree vertices in this facet
        std::vector<VertexSPtr> determined_hd_vertices;
        for (VertexSPtr v : facet->vertices()) {
            if (v->degree() > 3 && canonical_faces[v].size() == 3) {
                determined_hd_vertices.push_back(v);
            }
        }

        std::size_t dhd_count = determined_hd_vertices.size();

        if(dhd_count >= 3) {
            // Red: triangulate and add resulting triangles to the end of the list
            std::list<FacetSPtr> triangles = PolyhedronTransformation::triangulate(facet, polyhedron);
            for (FacetSPtr tri : triangles) {
                facets_to_process.push_back(tri);
            }
            continue; // skip the rest for this facet
        }

        std::vector<Point3SPtr> dhd_points;
        for (VertexSPtr v : determined_hd_vertices) {
            dhd_points.push_back(v->getPoint());
        }

        facet->perturbPlaneCoefficientsFixedPoints(1e-15, dhd_points);

        // Update all high degree vertices of the facet
        for (VertexSPtr v : facet->vertices()) {
            if (v->degree() > 3 && canonical_faces[v].size() < 3) {
                canonical_faces[v].push_back(facet);
                if (canonical_faces[v].size() == 3) {
                    // Recompute its position as intersection of its 3 incident planes
                    std::array<Plane3SPtr, 3> planes;
                    for (int j = 0; j < 3; ++j) {
                        planes[j] = canonical_faces[v][j]->plane();
                    }

                    unstable stuff, not how it should be done
                    vvvvvvvvvvv
                    Point3SPtr new_point = KernelWrapper::intersection(planes[0], planes[1], planes[2]);
                    if (new_point) {
                        v->setPoint(new_point);
                    }
                }
            }
        }
    }


}






// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //






std::vector<VertexSPtr> dhd_vertices;
bool predetermined = true;
for (VertexSPtr v : facet->vertices()) {
    if (v->degree() > 3 && determining_faces[v].size() >= 2) {
        dhd_vertices.push_back(v);
        std::cout << "  V" << v->getID() << " determined by";
        for (FacetSPtr inc_f : determining_faces[v]) {
            std::cout << " F" << inc_f->getID();
        }
        std::cout << std::endl;
        predetermined = predetermined && (determining_faces[v].size() == 3);
    }
}
predetermined = predetermined && (dhd_vertices.size() >= 2);

std::size_t dhd_count = dhd_vertices.size();
std::cout << "  predetermined: " << predetermined << std::endl;






// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //







// For triangles, we can have more than 2 determined vertices.
// The case of 3 predetermined vertices is handled before the main loop
// The case of 2 predetermined vertices and 1 determined (during visit) is handled in the normal loop
// We cannot have a triangle with two predetermined vertices and no new determined vertices:
//
if (facet->isTriangle() && predetermined_hdv_n == 2 && determined_n == 0) {
    std::vector<Point3SPtr> fps;
    for (VertexSPtr fv : fixing_vertices[facet]) {
        fps.push_back(fv->getPoint());
    }
    CGAL_assertion(fps.size() == 2);
    facet->perturbPlaneCoefficientsFixedPoints(range, fps);

    std::cout << "Nudge and fix F" << facet->getID() << " predetermined (2)" << std::endl;
    dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_predetermined.OFF", facet);
}



 else {
    // triangle face

    // A special case: if we have 2 fully predetermined vertices from start,
    // then no additional vertex can get determined by the addition of 'facet' (otherwise we would
    // triangulate), except if it's a triangle
    unsigned int predetermined_hdv_n = 0;
    for (VertexSPtr fv : fixing_vertices[facet]) {
        CGAL_assertion(determining_facets[fv].size() == 3);
        if (std::find(determining_facets[fv].begin(), determining_facets[fv].end(), facet) == determining_facets[fv].end()) {
            ++predetermined_hdv_n;
        }
    }

    CGAL_assertion(predetermined_hdv_n < 3);
}


if (predetermined_hdv_n == 3) {
    CGAL_assertion(facet->vertices().size() == 3); // otherwise we would have triangulated
    // vertices are determined, so they were nudged, so the plane is nudged
    // but can there be instability...?
    facet->initPlane();
    facet->normalizePlaneCoefficients();

    std::cout << "Fix F" << facet->getID() << " predetermined (3)" << std::endl;
    dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_fixed_triangle.OFF", facet);

    continue;
}





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //






static void randTiltPlanesv2(PolyhedronSPtr polyhedron);

void PolyhedronTransformation::randTiltPlanesv2(PolyhedronSPtr polyhedron) {
    double range = 1e-15; // Small nudge

    // Step 1: Identify high-degree vertices
    std::set<VertexSPtr> high_degree_vertices;
    for (VertexSPtr v : polyhedron->vertices()) {
        if (v->degree() > 3) {
            high_degree_vertices.insert(v);
        }
    }

    // Step 2: Classify facets as green (no high-degree vertex) or blue (at least one high-degree vertex)
    std::vector<FacetSPtr> green_facets;
    std::vector<FacetSPtr> blue_facets;
    for (FacetSPtr facet : polyhedron->facets()) {
        bool is_blue = false;
        for (VertexSPtr v : facet->vertices()) {
            if (high_degree_vertices.count(v)) {
                std::cout << "Facet " << facet->getID() << " is blue because of " << *(v->getPoint()) << std::endl;
                is_blue = true;
                break;
            }
        }
        if (is_blue) {
            blue_facets.push_back(facet);
        } else {
            green_facets.push_back(facet);
        }

        std::cout << "Facet " << facet->getID() << " is " << (is_blue ? "blue" : "green") << std::endl;
    }

    std::cout << green_facets.size() << " green facets" << std::endl;
    std::cout << blue_facets.size() << " blue facets" << std::endl;

    // Step 3: Nudge all green facets (randomly perturb their plane coefficients)
    for (FacetSPtr facet : green_facets) {
        facet->perturbPlaneCoefficientsNudge(range);
    }

    // Step 4: Triangulate all blue facets, creating red facets
    std::vector<FacetSPtr> red_facets;
    for (FacetSPtr facet : blue_facets) {
        std::list<FacetSPtr> triangles = PolyhedronTransformation::triangulate(facet, polyhedron);
        red_facets.insert(red_facets.end(), triangles.begin(), triangles.end());
    }

    polyhedron->initializeAllIDs();

    db::_3d::OBJFile::save("results/triangulated.obj", polyhedron,
                           false /*do_triangulate*/,
                           true /*convert_to_double*/);

    // Step 5 & 6: For all vertices, nudge and project according to incident green facets
    for (VertexSPtr v : polyhedron->vertices()) {
        // Find all incident green facets
        std::vector<FacetSPtr> incident_green_facets;
        for (FacetWPtr fw : v->facets()) {
            if (!fw.expired()) {
                FacetSPtr f = FacetSPtr(fw);
                if (std::find(green_facets.begin(), green_facets.end(), f) != green_facets.end()) {
                    incident_green_facets.push_back(f);
                }
            }
        }
        size_t n_green = incident_green_facets.size();
        CGAL_assertion(n_green <= 3);
        if (n_green == 3) {
            // All incident facets are green, just recompute position (green planes have been nudged)
            resetPoint(v);
        } else {
            Point3SPtr p = v->getPoint();
            std::array<double, 3> v_r = randVec(-range/2.0, range/2.0);
            Point3SPtr p_nudged = KernelFactory::createPoint3(CGAL::to_double(p->x()) + v_r[0],
                                                              CGAL::to_double(p->y()) + v_r[1],
                                                              CGAL::to_double(p->z()) + v_r[2]);

            Point3SPtr p_new;
            if (n_green == 0) {
                p_new = p_nudged;
            } else if (n_green == 1) {
                Plane3SPtr plane = incident_green_facets[0]->plane();
                p_new = KernelWrapper::projection(plane, p_nudged);
            } else if (n_green == 2) {
                Plane3SPtr plane1 = incident_green_facets[0]->plane();
                Plane3SPtr plane2 = incident_green_facets[1]->plane();
                Line3SPtr line = KernelWrapper::intersection(plane1, plane2);
                p_new = KernelWrapper::projection(line, p_nudged);
            }
            v->setPoint(p_new);
        }
    }

    // Recompute planes for all red facets
    for (FacetSPtr facet : polyhedron->facets()) {
        if (std::find(green_facets.begin(), green_facets.end(), facet) == green_facets.end()) {
            facet->initPlane();
        }
    }

    db::_3d::OBJFile::save("results/tilt_v2.obj", polyhedron,
                           false /*do_triangulate*/,
                           true /*convert_to_double*/);
    db::_3d::OBJFile::save("results/tilt_v2-triangulated.obj", polyhedron,
                           true /*do_triangulate*/,
                           true /*convert_to_double*/);

    // @debug
    for (FacetSPtr facet : polyhedron->facets()) {
        std::cout << "check facet " << facet->getID() << std::endl;
        std::cout << "green? " << (std::find(green_facets.begin(), green_facets.end(), facet) != green_facets.end()) << std::endl;
        for (VertexSPtr v : facet->vertices()) {
            CGAL_assertion(facet->getPlane()->has_on(*(v->getPoint())));
        }
    }
}






// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //






} else if (fixed_points.size() == 2) {
        // 2 fixed points: construct a plane through both points, nudge the normal within the allowed family
        Point3SPtr p0 = fixed_points[0];
        Point3SPtr p1 = fixed_points[1];
        CGAL_assertion(*p0 != *p1);

        const CGAL::FT& p0x = p0->x();
        const CGAL::FT& p0y = p0->y();
        const CGAL::FT& p0z = p0->z();
        const CGAL::FT& p1x = p1->x();
        const CGAL::FT& p1y = p1->y();
        const CGAL::FT& p1z = p1->z();

        // Step 1: Direction vector between points
        CGAL::FT ux = p1x - p0x;
        CGAL::FT uy = p1y - p0y;
        CGAL::FT uz = p1z - p0z;
        CGAL::FT uu = ux*ux + uy*uy + uz*uz;

        // Step 2: Original normal
        const CGAL::FT& a0 = plane_->a();
        const CGAL::FT& b0 = plane_->b();
        const CGAL::FT& c0 = plane_->c();

        // Step 2: Project original normal onto plane orthogonal to u
        CGAL::FT dot = a0*ux + b0*uy + c0*uz;
        CGAL::FT ab = a0 - dot * ux / uu;
        CGAL::FT bb = b0 - dot * uy / uu;
        CGAL::FT cb = c0 - dot * uz / uu;

        // std::cout << "a0 = " << exact(a0) << std::endl;
        // std::cout << "ux = " << exact(ux) << std::endl;
        // std::cout << "dot = " << exact(dot) << std::endl;
        // std::cout << "uu = " << exact(uu) << std::endl;

#ifdef CGAL_SS3_CHECK_ALMOST_DEGENERATE_VECTORS_IN_FACET_PERTURBATIONS // these shouldn't be needed because in general the original normal is sane and very far from 'u'
        const double eps = 1e-12;

        // If the projection is (almost) zero, pick any orthogonal normal
        if (CGAL::abs(ab) < eps && CGAL::abs(bb) < eps && CGAL::abs(cb) < eps) {
            ab = 0;
            bb = uz;
            cb = -uy;
            if (CGAL::abs(bb) < eps && CGAL::abs(cb) < eps) {
                ab = -uz;
                bb = 0;
                cb = ux;
            }
        }
#endif

#if 0
        std::cout << "uy = " << exact(uy) << std::endl;
        std::cout << "cb = " << exact(cb) << std::endl;
        std::cout << "uz = " << exact(uz) << std::endl;
        std::cout << "bb = " << exact(bb) << std::endl;
#endif

        // Step 3: Find a direction to nudge (cross product)
        CGAL::FT vx = uy * cb - uz * bb;
        CGAL::FT vy = uz * ab - ux * cb;
        CGAL::FT vz = ux * bb - uy * ab;

#ifdef CGAL_SS3_CHECK_ALMOST_DEGENERATE_VECTORS_IN_FACET_PERTURBATIONS // these shouldn't be needed because in general the original normal is sane and very far from 'u'
        // If v is zero, pick another orthogonal direction
        if (CGAL::abs(vx) < eps && CGAL::abs(vy) < eps && CGAL::abs(vz) < eps) {
            // Use cross with (1,0,0)
            vx = 0;
            vy = uz;
            vz = -uy;
            if (CGAL::abs(vy) < eps && CGAL::abs(vz) < eps) {
                vx = -uz;
                vy = 0;
                vz = ux;
            }
        }
#endif

        // Step 4: Nudge the normal
        CGAL::FT epsilon = nudge_to_simplest_rational_in_interval(rdist(gen));

        // std::cout << "ab = " << exact(ab) << std::endl;
        // std::cout << "vx = " << exact(vx) << std::endl;

        CGAL::FT a1 = ab + epsilon * vx;
        CGAL::FT b1 = bb + epsilon * vy;
        CGAL::FT c1 = cb + epsilon * vz;

        // std::cout << "a1 = " << exact(a1) << std::endl;

        // Step 5: Compute d so plane passes through p0
        CGAL::FT d1 = - (a1 * p0x + b1 * p0y + c1 * p0z);
        plane_ = KernelFactory::createPlane3(a1, b1, c1, d1);

        CGAL_postcondition(plane_->has_on(*p0));
        CGAL_postcondition(plane_->has_on(*p1));
    } else {
        DEBUG_PRINT("Error: called fixed point facet perturbation with > 2 fixed points");
    }






// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //






        } else if (fixed_points.size() == 2) {
        // 2 fixed points: construct a plane through both points, nudge the normal within the allowed family
        Point3SPtr p0 = fixed_points[0];
        Point3SPtr p1 = fixed_points[1];
        CGAL_assertion(*p0 != *p1);

        // Step 1: Direction vector between points
        CGAL::FT ux = p1->x() - p0->x();
        CGAL::FT uy = p1->y() - p0->y();
        CGAL::FT uz = p1->z() - p0->z();
        CGAL::FT uu = ux*ux + uy*uy + uz*uz;

        // Step 2: Project original normal onto plane orthogonal to u
        const CGAL::FT& a0 = plane_->a();
        const CGAL::FT& b0 = plane_->b();
        const CGAL::FT& c0 = plane_->c();
        CGAL::FT dot = a0*ux + b0*uy + c0*uz;
        CGAL::FT ab = a0 - dot * ux / uu;
        CGAL::FT bb = b0 - dot * uy / uu;
        CGAL::FT cb = c0 - dot * uz / uu;

        // Step 3: Pick a simple seed vector not parallel to u
        CGAL::FT sx = 1, sy = 0, sz = 0;
        if (CGAL::abs(ux) > 0.99) { // if u ~= (1,0,0), let's pick something else
            sx = 0; sy = 1; sz = 0;
        }

        // Step 4: Compute a vector orthogonal to u (seed^u)
        // No normalization needed: all vectors are already normalized or nearly so
        CGAL::FT vx = sy * uz - sz * uy;
        CGAL::FT vy = sz * ux - sx * uz;
        CGAL::FT vz = sx * uy - sy * ux;

        // Step 5: Nudge projected normal by a small random multiple of v
        CGAL::FT lambda = nudge_to_simplest_rational_in_interval(rdist(gen));
        CGAL::FT a1 = ab + lambda * vx;
        CGAL::FT b1 = bb + lambda * vy;
        CGAL::FT c1 = cb + lambda * vz;

        // Step 6: Compute d so plane passes through p0
        // Since the normal is orthogonal to u, it will pass through both p0 and p1
        CGAL::FT d1 = - (a1 * p0->x() + b1 * p0->y() + c1 * p0->z());
        plane_ = KernelFactory::createPlane3(a1, b1, c1, d1);
        CGAL_postcondition(plane_->has_on(*p0));
        CGAL_postcondition(plane_->has_on(*p1));
    } else {
        DEBUG_PRINT("Error: called fixed point facet perturbation with > 2 fixed points");
    }





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //





// Triangular facets always have lowest priority
bool a_is_triangle = a->isTriangle();
bool b_is_triangle = b->isTriangle();
if (a_is_triangle && !b_is_triangle) { return false; }
if (!a_is_triangle && b_is_triangle) { return true; }







// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


void PolyhedronTransformation::randTiltPlanesv3(PolyhedronSPtr polyhedron) {
    const double range = 1e-10; // Small nudge

    unsigned int had_to_triangulate_n = 0;

    // high degree vertex --> first 3 incident facets determining the vertex
    std::map<VertexSPtr, std::set<FacetSPtr> > determining_facets;

    // facet --> first 2 determined high degree vertices
    //
    // The facet becomes fixed at 2 vertices and not 3 vertices despite the vertices being perturbed
    // because if we do 3 random perturbations of vertices, the normal can vary wildly.
    //
    // Ideally, it could be fixed with 3 high degree vertices and a smarter perturbation (which needs
    // to take into account all incident facets of these 3 fixing vertices...)
    std::map<FacetSPtr, std::set<VertexSPtr> > fixing_vertices;

    auto dump_facet = [](std::string filename, FacetSPtr f) {
        auto pcdt = triangulate_facet_with_CDT2(f);

        using PCDT = decltype(pcdt);
        using PCDT_VH = typename PCDT::Vertex_handle;
        using PCDT_FH = typename PCDT::Face_handle;

        std::unordered_map<PCDT_FH, bool> in_domain_map;
        boost::associative_property_map<std::unordered_map<PCDT_FH, bool>> in_domain(in_domain_map);
        CGAL::mark_domain_in_triangulation(pcdt, in_domain);

        std::map<PCDT_VH, std::size_t> point_to_id;
        std::vector<Point3> points;
        std::vector<std::vector<std::size_t> > triangles;
        for (PCDT_VH vh : pcdt.finite_vertex_handles()) {
            point_to_id[vh] = points.size();
            points.push_back(vh->point());
        }

        for (PCDT_FH fh : pcdt.finite_face_handles()) {
            if(!get(in_domain, fh)) {
                continue;
            }

            triangles.push_back({point_to_id[fh->vertex(0)],
                                 point_to_id[fh->vertex(1)],
                                 point_to_id[fh->vertex(2)]});
        }

        CGAL::IO::write_OFF(filename, points, triangles, CGAL::parameters::stream_precision(17));
    };

    auto has_high_degree_vertices = [](FacetSPtr f) -> bool {
        for (VertexSPtr v : f->vertices()) {
            if (v->degree() > 3) {
                return true;
            }
        }
        return false;
    };

    // @debug
    unsigned int visited_face_id = 0;
    unsigned int nudged_face_id = 0;

    std::list<FacetSPtr> facets_to_process;

    // Sort by number of high degree vertices as to avoid triangulating as much as possible
    auto facet_sorter = [&](FacetSPtr a, FacetSPtr b)
    {
        auto hdv_count = [](FacetSPtr f) -> unsigned int {
            unsigned int hdv_n = 0;
            for (VertexSPtr v : f->vertices()) {
                if (v->degree() > 3) {
                    ++hdv_n;
                }
            }
            return hdv_n;
        };

#if 0
        // Triangular facets always have lowest priority
        bool a_is_triangle = a->isTriangle();
        bool b_is_triangle = b->isTriangle();
        if (a_is_triangle && !b_is_triangle) { return false; }
        if (!a_is_triangle && b_is_triangle) { return true; }
        return hdv_count(a) > hdv_count(b);
#elif 0
        // Give priority to the largest facets in term of vertices
        return a->vertices().size() > b->vertices().size();
#else
        // Give priority to facets with no determined vertices (using determining_facets)
        // If both or neither have constrained vertices, give priority to the largest hdv count
        //
        // The point is to avoid cascading exact number types, even if we have to triangulate a little more
        auto get_determined_count = [&](FacetSPtr f) -> unsigned int
        {
          unsigned int res = 0;
          for (VertexSPtr v : f->vertices()) {
              auto it = determining_facets.find(v);
              if (it != determining_facets.end() && it->second.size() == 3) {
                  ++res;
              }
          }
          return res;
        };

        unsigned int a_determined_n = get_determined_count(a);
        unsigned int b_determined_n = get_determined_count(b);
        if (a_determined_n != b_determined_n) {
            // Give priority to the one with the least determined vertices
            return a_determined_n < b_determined_n;
        }

        // same number of determined vertices, give priority to the facet with the most high degree vertices
        unsigned int a_hdv_n = hdv_count(a);
        unsigned int b_hdv_n = hdv_count(b);
        if (a_hdv_n != b_hdv_n) {
            // Give priority to the one with the most high degree vertices
            return a_hdv_n > b_hdv_n;
        }

        // same number of determined vertices and high degree vertices, give priority to the largest facet
        return a->vertices().size() > b->vertices().size();
#endif
    };

    for (FacetSPtr facet : polyhedron->facets()) {
        if (has_high_degree_vertices(facet)) {
            facets_to_process.push_back(facet);
        } else {
            // If the facet has no high degree vertices, we can just tilt it randomly and it will
            // be fine because by definition all of its vertices are degree 3 and will be stable
            // because an unstable configuration results from almost coplanar facets, which have been
            // merged ahead of randomization.
            // UNLESS we have to triangulate a face incident to one vertex of this face without
            // high degree vertices and then the face now has a high degree vertex. If that happens,
            // we want the high degree vertex to constrain to the (up to 2) facets with no high degree
            // vertices which we are constraining here ahead of the flooding process.
            // Hence, we mark this as fixed with dummy vertices and add 'v' (a non high degree vertex)
            // to the 'determining_facets' map.
            std::cout << "Nudge and fix F" << facet->getID() << std::endl;
            facet->perturbPlaneCoefficientsNudge(range);

            // @debug dump the face into a single OFF file
            dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_low_degree.OFF", facet);

            // the point of this is that the facet with low degree facet is still used as a constraining place
            // if nudging a vertex incident to it that would become high degree after triangulation
            fixing_vertices[facet].insert(facet->vertices().front());
            fixing_vertices[facet].insert(facet->vertices().back());

            // the point of this is that if the vertex becomes high degree after triangulation,
            // one (or two) facet with low degree vertices will appear in the determining facets
            for (VertexSPtr v : facet->vertices()) {
                determining_facets[v].insert(facet);
            }
        }
    }

    // Forward declarations for mutually recursive lambdas
    std::function<void(FacetSPtr, VertexSPtr)> add_fixing_vertex;
    std::function<void(VertexSPtr)> determine_vertex;

    auto nudge_determined_vertex = [&](VertexSPtr v) {
        std::vector<Plane3SPtr> constraining_planes;
        for (FacetSPtr df : determining_facets[v]) {
            if (fixing_vertices[df].size() >= 2) {
                constraining_planes.push_back(df->getPlane());
                std::cout << "  F" << df->getID() << " constrains the nudge" << std::endl;
            }
        }

        CGAL_assertion(constraining_planes.size() <= 3);

        const size_t n_fixed = constraining_planes.size();
        if (n_fixed == 3) {
            resetPoint(v, { constraining_planes[0], constraining_planes[1], constraining_planes[2] });
            std::cout << "V" << v->getID() << " reset to " << *v->getPoint() << std::endl;
            return;
        }

        Point3SPtr p = v->getPoint();

        static std::random_device rd;
        unsigned int s = 0; // rd()
        // std::cout << "seed = " << s << std::endl;
        static std::mt19937 gen(s);
        static std::uniform_real_distribution<> rdist(-range, range);

        auto nudge = [&](const CGAL::FT& v) {
            // Since we are perturbing, we might as well collapse the DAG of 'v'.
            // the point is also that once 'nv' is a double, its interval will be a singleton,
            // and we will have access to static filters
            double step = rdist(gen);
            double nv = CGAL::to_double(v) + step;
            return nv;
        };

        auto nudge_to_simplest_rational_in_interval = [&](const CGAL::FT& v) {
            double d1 = nudge(v);
            double d2 = nudge(v);
            if (d2 < d1) {
                std::swap(d1, d2);
            }
            CGAL::FT nv = CGAL::simplest_rational_in_interval<CGAL::K::Exact_kernel::FT>(d1, d2);
            return nv;
        };

#if 0
        std::array<double, 3> v_r = randVec(-range/2.0, range/2.0);
        double x = CGAL::to_double(p->x()) + v_r[0];
        double y = CGAL::to_double(p->y()) + v_r[1];
        double z = CGAL::to_double(p->z()) + v_r[2];
        Point3SPtr p_nudged = KernelFactory::createPoint3(Point3(x, y, z));
#else
        CGAL::FT x = nudge_to_simplest_rational_in_interval(p->x());
        CGAL::FT y = nudge_to_simplest_rational_in_interval(p->y());
        CGAL::FT z = nudge_to_simplest_rational_in_interval(p->z());
        Point3SPtr p_nudged = KernelFactory::createPoint3(x, y, z);
#endif

        Point3SPtr p_new;

        if (n_fixed == 0) {
            p_new = p_nudged;
        } else if (n_fixed == 1) {
            Plane3SPtr plane = constraining_planes[0];
#if 0
            p_new = KernelWrapper::projection(plane, p_nudged);
#else
            // something similar but a little more subtle:
            // 1. project the point onto the plane
            // 2. express the point as a linear combination of the plane's origin and basis: pp = o + l1 * b1 + l2 * b2
            // 3. nudge l1 and l2 to l1' and l2' with a random interval around l1 and l2, and
            //    simplest_rational_in_interval
            // 4. recompute the point as pp = o + l1' * b1 + l2' * b2
            Point3SPtr pp = KernelWrapper::projection(plane, p_nudged);
            const Point3& o = plane->point();
            const Vector3& b1 = plane->base1();
            const Vector3& b2 = plane->base2();
            CGAL::FT l1 = CGAL::scalar_product(*pp - o, b1);
            CGAL::FT l2 = CGAL::scalar_product(*pp - o, b2);
            CGAL::FT nl1 = nudge_to_simplest_rational_in_interval(l1);
            CGAL::FT nl2 = nudge_to_simplest_rational_in_interval(l2);
            p_new = KernelFactory::createPoint3(o.x() + nl1 * b1.x() + nl2 * b2.x(),
                                                o.y() + nl1 * b1.y() + nl2 * b2.y(),
                                                o.z() + nl1 * b1.z() + nl2 * b2.z());
#endif
        } else if (n_fixed == 2) {
            Plane3SPtr plane1 = constraining_planes[0];
            Plane3SPtr plane2 = constraining_planes[1];
            Line3SPtr line = KernelWrapper::intersection(plane1, plane2);
#if 0
            p_new = KernelWrapper::projection(line, p_nudged);
#else
            // something similar but a little more subtle:
            // 1. project the point onto the line
            // 2. express the point as a linear combination of the line's origin and basis: pp = o + l * v
            // 3. nudge l to l' with a random interval around l, and simplest_rational_in_interval
            // 4. recompute the point as pp = o + l' * v
            Point3SPtr pp = KernelWrapper::projection(line, p_nudged);
            const Point3& o = line->point();
            const Vector3& d = line->to_vector();
            CGAL::FT l = CGAL::scalar_product(*pp - o, d);
            CGAL::FT nl = nudge_to_simplest_rational_in_interval(l);
            p_new = KernelFactory::createPoint3(o.x() + nl * d.x(),
                                                o.y() + nl * d.y(),
                                                o.z() + nl * d.z());
#endif
        }

        std::cout << "V" << v->getID() << " nudged to " << *p_new << std::endl;
        v->setPoint(p_new);
    };

    auto is_facet_fixed = [&](FacetSPtr f) -> bool {
        CGAL_assertion(fixing_vertices[f].size() <= 3);
        return (f->isTriangle() && fixing_vertices[f].size() == 3) ||
               (!f->isTriangle() && fixing_vertices[f].size() == 2);
    };

    // The face has enough determined vertices to be fixed. Compute its random perturbation,
    // and add the facet ID to its vertices
    add_fixing_vertex = [&](FacetSPtr f, VertexSPtr v) {
        CGAL_precondition(fixing_vertices[f].size() <= 3);

        if (is_facet_fixed(f)) {
            std::cout << "  F" << f->getID() << " is already fully fixed" << std::endl;
            return;
        }

        std::cout << "  Fix F" << f->getID() << " with V" << v->getID() << std::endl;
        fixing_vertices[f].insert(v);

        if (!is_facet_fixed(f)) {
            // nothing to do yet, there are still degrees of freedom in the facet
            return;
        }

        std::cout << "F" << f->getID() << " is now fully fixed by";
        for (VertexSPtr fv : fixing_vertices[f]) {
            std::cout << " V" << fv->getID();
        }
        std::cout << std::endl;

        if (f->isTriangle()) {
            CGAL_assertion(fixing_vertices[f].size() == 3); // just to be clear

            // for triangles, all vertices are determined, and there is nothing to nudge
            // (note that vertices were themselves nudged so the facet is nudged).
            f->initPlane();
            f->normalizePlaneCoefficients();

            dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_fixed_3.OFF", f);
            return;
        }

        std::vector<Point3SPtr> fixed_points;
        for (VertexSPtr fixed_v : fixing_vertices[f]) {
            fixed_points.push_back(fixed_v->getPoint());
        }

        std::cout << "Nudge and fix F" << f->getID() << " with " << fixed_points.size() << " fixed point(s)" << std::endl;
        f->perturbPlaneCoefficientsFixedPoints(range, fixed_points);

        dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_fixed_" + std::to_string(fixed_points.size()) + ".OFF", f);

        // Need to now tag the vertices of the facet
        for (VertexSPtr v : f->vertices()) {
            if (v->degree() <= 3) {
                continue;
            }

            if (determining_facets[v].size() < 3) {
                determining_facets[v].insert(f);
                if (determining_facets[v].size() == 3) {
                    // When the vertex becomes fixed (its 3 determining facets become known), we need:
                    // - to perturb the position of the vertex
                    // - to update all incident facets to check if they are now fixed and in that case,
                    //   compute their plane coefficients
                    determine_vertex(v);
                }
            }
        }
    };

    determine_vertex = [&](VertexSPtr v) {
        CGAL_precondition(determining_facets[v].size() == 3);

        std::cout << "V" << v->getID() << " is now fully determined by F";
        auto it = determining_facets[v].begin();
        std::cout << (*it++)->getID() << " F";
        std::cout << (*it++)->getID() << " F";
        std::cout << (*it)->getID() << std::endl;

        // set the nudged position for the vertex: a nudge constrained by already fixed incident facets
        nudge_determined_vertex(v);

        // compute the plane coefficients of any incident facet that becomes fixed
        // by this vertex becoming determined
        for (FacetWPtr wf : v->facets()) {
            if (FacetSPtr f = wf.lock()) {
                add_fixing_vertex(f, v);
            }
        }
    };

    auto is_facet_overconstrained = [&](FacetSPtr f) -> bool {
        if (f->vertices().size() == 3) {
            return false;
        }

        // If the facet has too many fixed vertices after triangulation, it is overconstrained
        if (fixing_vertices[f].size() > 2) {
            return true;
        }

        // we cannot handle that facet without triangulation if adding this facet to high degree
        // vertices would create too many determined vertices (> 2) in any facet incident to
        // the determined high degree vertices of this facet
        std::map<FacetSPtr, unsigned int> facets_to_test; // facets + number of appearances
        for (VertexSPtr hdv : f->vertices()) {
            if (hdv->degree() > 3) {
                for (FacetWPtr inc_f : hdv->facets()) {
                    FacetSPtr f = inc_f.lock();
                    if (f) {
                        ++facets_to_test[f];
                    }
                }
            }
        }

        for (const auto& [ft, count] : facets_to_test) {
            // If the two facets appear more than twice in high degree vertices,
            // then it means we have at least 4 high degree vertices along the common line
            // between two facets (whichever the perturbation). This can be tricky because
            // if we - for example - have a third facet that is sharing two of these
            // high degree vertices and then there's no way to nudge this.
            //
            // @todo we can definitely do better than just triangulate everything: for example,
            // we could just split the facet in two (pick the edge among the Delaunay edges
            // of the triangulated facet that would split the facet in a nice manner).
            //
            // @todo we should split the smaller of the two facets
            if (ft != f && count > 2) {
                std::cout << "Sharing too many vertices with F" << ft->getID() << " (" << count << ")" << std::endl;
                return true;
            }

            // Count the number of high-degree vertices with either:
            // - 3 determining facets
            // - 2 determining facets and incident to 'facet'
            // These are vertices that are determined, or would be determined once we 'add'
            // the facet to its high degree vertices.
            unsigned int constrain_n = 0;
            for (VertexSPtr v : f->vertices()) {
                if (v->degree() > 3 ) {
                    if (determining_facets[v].size() == 3) {
                        ++constrain_n;
                    } else if (determining_facets[v].size() == 2 && ft->hasVertex(v)) {
                        ++constrain_n;
                    }
                }

                if(constrain_n > 2) {
                    std::cout << "F" << ft->getID() << " would be over constrained by fixing of F" << f->getID() << std::endl;
                    return true;
                }
            }
        }

        return false;
    };

    while (!facets_to_process.empty()) {
        facets_to_process.sort(facet_sorter);
        FacetSPtr facet = facets_to_process.front();
        facets_to_process.pop_front();

        std::cout << "Pop F" << facet->getID() << std::endl;
        CGAL_assertion(fixing_vertices[facet].size() <= 2);
        std::cout << "  Fixing vertices:";
        for (VertexSPtr fv : fixing_vertices[facet]) {
            std::cout << " V" << fv->getID();
        }
        std::cout << std::endl;

        dump_facet("results/visited_face_" + std::to_string(visited_face_id++) + ".OFF", facet);

        CGAL_assertion(facet->vertices().size() >= 3);

        if (is_facet_overconstrained(facet)) {
            std::list<FacetSPtr> facets_to_triangulate = { facet };
            std::list<FacetSPtr> final_facets;
            while (!facets_to_triangulate.empty())
            {
                FacetSPtr facet_tt = facets_to_triangulate.front();
                facets_to_triangulate.pop_front();

                // @debug
                // the facet is not yet fixed, so no vertex can have it as determining facet
                CGAL_assertion(fixing_vertices[facet_tt].size() < 2);
                for (VertexSPtr v : facet_tt->vertices()) {
                    CGAL_assertion(determining_facets[v].size() <= 3);
                    CGAL_assertion(determining_facets[v].count(facet_tt) == 0);
                }

                std::cout << "Triangulate F" << facet_tt->getID() << std::endl;
                ++had_to_triangulate_n;

                const Triangulation_strategy strategy = Triangulation_strategy::DEFAULT;
                std::list<FacetSPtr> new_facets = PolyhedronTransformation::triangulate(facet_tt, polyhedron, strategy);

                db::_3d::OBJFile::save("results/post-triangulate.obj", polyhedron,
                                       false /*do_triangulate*/,
                                       true /*convert_to_double*/);

                VertexSPtr va, vb;
                if (strategy == Triangulation_strategy::MID_CUT) {
                    FacetSPtr f1 = new_facets.front(), f2 = new_facets.back();

                    // Identify the common edge, i.e. the two vertices that are shared by both facets.
                    // @todo A little unefficient because we had that information during triangulating.
                    for (VertexSPtr v : f1->vertices()) {
                        if (f2->hasVertex(v)) {
                            if (!va) {
                                va = v;
                            } else {
                                vb = v;
                                break;
                            }
                        }
                    }

                    // even if we have to further subdivide, those vertices must be determined
                    for (VertexSPtr sv : {va, vb}) {
                        if (determining_facets[sv].size() == 3) {
                            continue;
                        }

                        determine_vertex(sv);
                    }
                }

                // existing already-determined vertices are fixed points for the new facets
                for (FacetSPtr nf : new_facets) {
                    std::cout << "spawned F" << nf->getID() << std::endl;

                    for (VertexSPtr v : nf->vertices()) {
                        if (determining_facets[v].size() == 3) {
                            fixing_vertices[nf].insert(v);
                        }
                    }

                    bool is_overconstrained = is_facet_overconstrained(nf);

                    if (!is_overconstrained && strategy == Triangulation_strategy::MID_CUT) {
                        // When using the mid cut strategy, we have to fix the cutting edge for both
                        // (non-triangle) new facets, otherwise any perturbation not involving these
                        // vertices will be very unstable and send the split vertices far.
                        if (!nf->isTriangle()) {
                            for (VertexSPtr v : nf->vertices()) {
                                if (v != va && v != vb && determining_facets[v].size() == 3) {
                                    is_overconstrained = true;
                                    break;
                                }
                            }
                        }
                    }

                    if (is_overconstrained) {
                        std::cout << "Newborn face F" << nf << " is over constrained" << std::endl;
                        facets_to_triangulate.push_back(nf); // further triangulate it
                        continue;
                    }


                    final_facets.push_back(nf);
                }

                for (FacetSPtr nf : final_facets) {
                    for (VertexSPtr v : nf->vertices()) {
                        if (determining_facets[v].size() == 3) {
                            std::cout << "newborn F" << nf->getID() << " is constrained by V" << v->getID() << std::endl;
                            fixing_vertices[nf].insert(v);
                        }
                    }

                    facets_to_process.push_back(nf);
                }
            }

            continue;
        }

        // Now, adding the facet to the high degree vertices will not over constrain the facet, so do it:
        for (VertexSPtr v : facet->vertices()) {
            if (v->degree() <= 3) {
                continue;
            }

            if (determining_facets[v].size() < 3) {
                determining_facets[v].insert(facet);
                if (determining_facets[v].size() == 3) {
                    // When the vertex becomes fixed (its 3 determining facets become known), we need:
                    // - to perturb the position of the vertex
                    // - to update all incident facets to check if they are now fixed and in that case,
                    //   compute their plane coefficients
                    determine_vertex(v);
                }
            }
        }
    }

    // Some facets might have high degree, but still some freedom of movement after the flooding, fix them
    for (FacetSPtr facet : polyhedron->facets()) {
        // if there were 0, it would be a facet without high degree vertices and those are nudged first
        // if there were 2, it would have been nudged when the 2nd vertex got determined
        if (fixing_vertices[facet].size() != 1) {
            continue;
        }

        VertexSPtr v = *(fixing_vertices[facet].begin());
        std::vector<Point3SPtr> fixed_points = { v->getPoint() };
        std::cout << "Nudge and fix F" << facet->getID() << " fixed by V" << v->getID() << " [remaining]" << std::endl;
        facet->perturbPlaneCoefficientsFixedPoints(range, fixed_points);

        dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_remaining.OFF", facet);
    }

    for (VertexSPtr v : polyhedron->vertices()) {
        CGAL_assertion(determining_facets[v].size() <= 3);
        CGAL_assertion(v->degree() == 3 || determining_facets[v].size() == 3);
    }

    for (FacetSPtr f : polyhedron->facets()) {
        std::cout << "check fixing_vertices[" << f->getID() << "].size() = " << fixing_vertices[f].size() << std::endl;
        CGAL_assertion(fixing_vertices[f].size() <= 3);
    }

    db::_3d::OBJFile::save("results/tilt_v3-pre_reset.obj", polyhedron,
                           false /*do_triangulate*/,
                           true /*convert_to_double*/);

    // Recompute all points which were not fixed
    for (VertexSPtr v : polyhedron->vertices()) {
       // @fixme high degree vertices that are not entirely determined are not yet recomputed
        if (v->degree() == 3) {
            std::cout << "Reset V" << v->getID() << std::endl;
            resetPoint(v);
        }
    }

    std::cout << "All facets processed" << std::endl;

    db::_3d::OBJFile::save("results/tilt_v3.obj", polyhedron,
                           false /*do_triangulate*/,
                           true /*convert_to_double*/);
    // db::_3d::OBJFile::save("results/tilt_v3-triangulated.obj", polyhedron,
    //                        true /*do_triangulate*/,
    //                        true /*convert_to_double*/);

    // @debug
    for (FacetSPtr facet : polyhedron->facets()) {
        for (VertexSPtr v : facet->vertices()) {
            std::cout << "Is V" << v->getID() << " on F" << facet->getID() << std::endl;
            CGAL_assertion(facet->getPlane()->has_on(*(v->getPoint())));
        }
    }

    // @debug
    std::cout << "Had to triangulate " << had_to_triangulate_n << " facets" << std::endl;

    unsigned int tr_n = 0;
    for (FacetSPtr facet : polyhedron->facets()) {
        if (facet->vertices().size() == 3) {
            ++tr_n;
        }
    }
    std::cout << tr_n << " triangle facets" << std::endl;

    // @debug
    for (VertexSPtr v : polyhedron->vertices()) {
        std::cout << "V" << v->getID() << " has depth " << CGAL::depth(*(v->getPoint())) << std::endl;
        // std::cout << "  " << exact(*(v->getPoint())) << std::endl;
    }

    for (FacetSPtr f : polyhedron->facets()) {
        std::cout << "F" << f->getID() << " has depth " << CGAL::depth(*(f->getPlane())) << std::endl;
    }
}






// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //





{
  namespace PMP = CGAL::Polygon_mesh_processing;

  DEBUG_PRINT("Inverting and adding a Bbox...");
  CGAL_precondition(!CGAL::is_empty(sm));

  auto fwm = sm.property_map<face_descriptor, double>("f:weight");
  CGAL_assertion(bool(fwm));

  struct Weight_setter_visitor
    : public CGAL::Polygon_mesh_processing::Triangulate_faces::Default_visitor<Mesh>
  {
    Mesh::Property_map<Mesh::Face_index, double> property;
    double weight = 1.0;
    void before_subface_creations(face_descriptor f_old) { weight = get(property, f_old); }
    void after_subface_created(face_descriptor f_new) { put(property, f_new, weight); }
  };

  Weight_setter_visitor visitor;
  visitor.property = *fwm;

  PMP::triangulate_faces(sm, CGAL::parameters::visitor(visitor));

  // check the sanity of the input
  bool has_SI = PMP::does_self_intersect(sm);
  if(has_SI) {
    std::cerr << "Error: input has self intersections" << std::endl;
    return false;
  }

  auto vol_id_map = sm.add_property_map<face_descriptor, std::size_t>().first;
  std::size_t vccn = PMP::volume_connected_components(sm, vol_id_map,
                                                      CGAL::parameters::do_orientation_tests(false));
  std::size_t ccn = PMP::internal::number_of_connected_components(sm);
  if(vccn != ccn) {
    std::cerr << "Error: input has nested connected components" << std::endl;
    return false;
  }

  PMP::orient_to_bound_a_volume(sm);
  PMP::reverse_face_orientations(sm);

#define CGAL_SS3_DO_NOT_USE_ENCLOSING_BBOX
#ifndef CGAL_SS3_DO_NOT_USE_ENCLOSING_BBOX
  const CGAL::Bbox_3 bb = PMP::bbox(sm, CGAL::parameters::bbox_scaling(100));
  Mesh bbox_mesh;
  CGAL::make_hexahedron(Iso_cuboid3(Point3(bb.xmin(), bb.ymin(), bb.zmin()),
                                    Point3(bb.xmax(), bb.ymax(), bb.zmax())),
                        bbox_mesh,
                        CGAL::parameters::do_not_triangulate_faces(false));

  // if the face:weight pmap exists, get the smallest value
  // as to assign an even smaller value to the bounding box's faces
  double min_weight = std::numeric_limits<double>::max(); // 'double' on purpose
  for(face_descriptor f : faces(sm)) {
    min_weight = (std::min)(min_weight, get(*fwm, f));
  }
  DEBUG_PRINT("min weight: " << min_weight)

  std::unordered_map<face_descriptor, face_descriptor> f2f;
  CGAL::copy_face_graph(bbox_mesh, sm,
                        CGAL::parameters::face_to_face_output_iterator(
                          std::inserter(f2f, f2f.end())));


  if(fwm) {
    for(const auto& e : f2f)
      put(*fwm, e.second, 1e-10 * min_weight);
  }
#endif

  if(CGAL::is_empty(sm) || !CGAL::is_closed(sm) || !CGAL::is_triangle_mesh(sm)) {
    std::cerr << "Error: empty or open output" << std::endl;
    return false;
  }

  return true;
}


bool
OutwardMeshOffset::
remove_bbox_and_invert(Mesh& sm)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  DEBUG_PRINT("Removing Bbox and inverting...");

  std::cout << vertices(sm).size() << " NV " << faces(sm).size() << " NF (with BBOX)" << std::endl;

  CGAL_precondition(!CGAL::is_empty(sm));
  CGAL_precondition(CGAL::is_valid_face_graph(sm) && sm.is_valid());

  if(CGAL::is_empty(sm)) {
    std::cerr << "Error: empty output" << std::endl;
    return false;
  } else if(!CGAL::is_closed(sm)) {
    std::cerr << "Error: open output" << std::endl;
    return false;
  }

#ifndef CGAL_SS3_DO_NOT_USE_ENCLOSING_BBOX
  auto vpm = get(CGAL::vertex_point, sm);

  vertex_descriptor extreme_v = *(vertices(sm).begin());
  for(vertex_descriptor v : vertices(sm)) {
    if(get(vpm, v).z() > get(vpm, extreme_v).z())
      extreme_v = v;
  }

  face_descriptor extreme_f = face(halfedge(extreme_v, sm), sm);
  CGAL_assertion(extreme_f != boost::graph_traits<Mesh>::null_face());

  std::vector<face_descriptor> fs { extreme_f };
  PMP::remove_connected_components(sm, fs);

  std::cout << vertices(sm).size() << " NV " << faces(sm).size() << " NF" << std::endl;

  if(CGAL::is_empty(sm)) {
    std::cerr << "Error: empty output" << std::endl;
    return false;
  }
#endif

  PMP::reverse_face_orientations(sm);

  return true;
}


// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //




// just to verify that we cannot just simply ignore the enforcement of vertices on plane supports
// 2 FAILURES on 639 test cases
PolyhedronSPtr p = convert(sm, true /*force merge*/);
p->initializeAllIDs();
std::cout << "IN: " << p->vertices().size() << " NV " << p->facets().size() << " NF" << std::endl;
PolyhedronTransformation::normalizeFacetPlanes(p);
PolyhedronTransformation::randTiltPlanes(p);
std::cout << "OUT: " << p->vertices().size() << " NV " << p->facets().size() << " NF" << std::endl;

PolyhedronTransformation::resetPoints(p);

CGAL_assertion(PolyhedronTransformation::doAll3PlanesIntersect(p));
std::cout << "generic position" << std::endl;

// below checks for vertices on plane, which we obviously won't have here
// CGAL_assertion(!SelfIntersection::hasSelfIntersectingSurface(p));
// std::cout << "no self intersections" << std::endl;

// run the skeleton code
algo::ControllerSPtr controller = { };
SimpleStraightSkelSPtr algoskel3d =
    SimpleStraightSkel::create(p, controller, save_offsets, save_path);
bool success = algoskel3d->run();
if (!success) {
    return false;
}




// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //




// From CGAL::Surface_mesh to the custom polyhedron data structure
// @tmp force simplification
PolyhedronSPtr polyhedron = convert(sm);
CGAL_assertion(polyhedron && polyhedron->isConsistent());

// Apply perturbations to ensure generic configuration
DEBUG_PRINT("Perturbing mesh...");

// Check if we can tilt facets' planes (i.e., nudge plane coefficients) directly.
// The advantage is that we then manipulate smaller meshes since faces are polygonal.
// @todo convert already knows this...
bool use_plane_tilts = PolyhedronTransformation::isTiltCompatible(polyhedron);
if (use_plane_tilts) {
    DEBUG_PRINT("Tilting the polyhedron's facets...");
    PolyhedronTransformation::normalizeFacetPlanes(polyhedron);
    PolyhedronTransformation::randTiltPlanes(polyhedron);
    PolyhedronTransformation::resetPoints(polyhedron);
} else {
    DEBUG_PRINT("Moving vertices randomly...");
    PolyhedronTransformation::normalizeFacetPlanes(polyhedron);
    PolyhedronTransformation::randMovePoints(polyhedron);
}

CGAL_assertion(polyhedron && polyhedron->isConsistent());
db::_3d::OBJFile::save("results/first_input.obj", polyhedron, false /*do not triangulate*/);

CGAL_expensive_assertion(PolyhedronTransformation::doAll3PlanesIntersect(polyhedron));
CGAL_expensive_assertion(!SelfIntersection::hasSelfIntersectingSurface(polyhedron));
// @todo In safe mode, this should be:
// if (bad) -> reduce amplitude of perturbation and try again"

DEBUG_PRINT("Constructing offset...");

algo::ControllerSPtr controller = { };

if (use_plane_tilts) {
    // run the skeleton code
    SimpleStraightSkelSPtr algoskel3d =
        SimpleStraightSkel::create(polyhedron, controller, save_offsets, save_path);
    bool success = algoskel3d->run();
    if (!success) {
        return false;
    }
} else {
    // this one is a little fancier: we split vertices, move forward a bit in time,
    // then try to remesh it to go back to a simple surface

    // Run the main offset loop with a visitor that breaks once the minimal edge length is above
    // a certain bound (or we have reached the desired value).

    DEBUG_PRINT("Running first pass");

    SimpleStraightSkelSPtr algoskel3d_f =
        SimpleStraightSkel::create(polyhedron, controller, save_offsets, save_path);

    // @todo this doesn't guarantee that all edges are big enough
    Mesh tmp_sm;
    utils::Stop_on_far_enough_event visitor(1e-1);
    algoskel3d_f->setVisitor(&visitor);

    // @todo handle the (unlikely) case where run() handled the whole offsetting process
    try {
        bool success = algoskel3d_f->run();
        if (!success) {
            return false;
        }
    } catch(utils::Far_enough_event) {
        std::cout << "caught the throw" << std::endl;
        bool success = db::_3d::Surface_meshIO::save(visitor.polyhedron_, tmp_sm,
                                                    true /*triangulate*/, false /*no doubles*/);
        CGAL_assertion(success);
    }

      CGAL::IO::write_polygon_mesh("results/second_input-surface_mesh.off", tmp_sm, CGAL::parameters::stream_precision(17));

      CGAL_assertion(tmp_sm.is_valid());
      CGAL_assertion(is_valid_face_graph(tmp_sm));
      CGAL_assertion(!CGAL::is_empty(tmp_sm));
      CGAL_assertion(CGAL::is_closed(tmp_sm));
      CGAL_assertion(CGAL::is_triangle_mesh(tmp_sm));
      CGAL_assertion(!PMP::has_degenerate_faces(tmp_sm));
      CGAL_assertion(!PMP::does_self_intersect(tmp_sm));

    // @todo instead of region growing, could we re-use the first coplanar partition and
    // edge merge? The split cannot separate two regions?
    polyhedron = convert(tmp_sm); // performs a coplanar retriangulation of the mesh

    CGAL_assertion(polyhedron && polyhedron->isConsistent());
    db::_3d::OBJFile::save("results/second_input-converted.obj", polyhedron, false /*do not triangulate*/);

    // Perturbing here is just for safety in the unlikely event that fusion of triangular facets
    // recreated a degenerate configuration
    use_plane_tilts = PolyhedronTransformation::isTiltCompatible(polyhedron);
    if (use_plane_tilts) {
        DEBUG_PRINT("Tilting the polyhedron's facets (2nd pass)...");
        PolyhedronTransformation::normalizeFacetPlanes(polyhedron);
        PolyhedronTransformation::randTiltPlanes(polyhedron);
        PolyhedronTransformation::resetPoints(polyhedron);
    } else {
        // We really should not be there
        DEBUG_PRINT("Warning: moving vertices randomly (again?!)...");
        PolyhedronTransformation::normalizeFacetPlanes(polyhedron);
        PolyhedronTransformation::randMovePoints(polyhedron);
    }

    CGAL_assertion(polyhedron && polyhedron->isConsistent());
    db::_3d::OBJFile::save("results/second_input.obj", polyhedron, false /*do not triangulate*/);
    CGAL_expensive_assertion(PolyhedronTransformation::doAll3PlanesIntersect(polyhedron));
    CGAL_expensive_assertion(!SelfIntersection::hasSelfIntersectingSurface(polyhedron));

    SimpleStraightSkelSPtr algoskel3d =
        SimpleStraightSkel::create(polyhedron, controller, save_offsets, save_path);
    bool success = algoskel3d->run();
    if (!success) {
        return false;
    }
}



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



struct Stop_after_n_events
    : public utils::Base_mesh_offset_visitor
{
    Stop_after_n_events(int n) : n_(n) { }

    bool go_further(int step_id, PolyhedronSPtr polyhedron, CGAL::FT offset) override {
        bool stop = (event_count_ >= n_);
        return !stop;
    }

    void on_save_offset_event(PolyhedronSPtr polyhedron, CGAL::FT offset) override { }

    void after_offset_event(PolyhedronSPtr polyhedron, CGAL::FT offset) override {
        ++event_count_;
    }

private:
    int n_;
    int event_count_ = 0;
};

class Far_enough_event
    : public std::exception
{
  const char* what() const throw () {
      return "Unauthorized intersections of constraints";
  }
};

// The point is that once the edges are long enough, we can merge faces and recompute planes
// and the small error is safe.
struct Stop_on_far_enough_event
    : public Base_mesh_offset_visitor
{
    Stop_on_far_enough_event(CGAL::FT min_event_distance) : min_event_distance_(min_event_distance) { }

    bool go_further(int step_id, PolyhedronSPtr polyhedron, CGAL::FT offset) override {
        return true;
    }

    void before_offset_event(PolyhedronSPtr polyhedron,
                             CGAL::FT current_offset,
                             AbstractEventSPtr event) override {

        // @speed this can be put beneath the delta filter, I'm putting it here for debugging purposes
        {
            EdgeSPtr min_length_edge;
            CGAL::FT sq_min_edge_length = (std::numeric_limits<double>::max)();
            for (EdgeSPtr edge : polyhedron->edges()) {
                CGAL::FT sq_l = CGAL::squared_distance(*(edge->getVertexSrc()->getPoint()),
                                                       *(edge->getVertexDst()->getPoint()));
                if (sq_l < sq_min_edge_length) {
                    min_length_edge = edge;
                    sq_min_edge_length = sq_l;
                }
            }
            std::cout << "min edge length @ " << current_offset << " = " << CGAL::approximate_sqrt(sq_min_edge_length) << std::endl;
            std::cout << min_length_edge->toString() << std::endl;
        }

        CGAL::FT event_offset = event->getOffset();
        std::cout << "event delta = " << CGAL::abs(event_offset - current_offset) << std::endl;
        if (CGAL::abs(event_offset - current_offset) < min_event_distance_) {
            return;
        }

        db::_3d::OBJFile::save("results/interrupted.obj", polyhedron, false /*do not triangulate*/);

        std::cout << "Event @ " << event_offset << " is far enough from " << current_offset << std::endl;

        CGAL::FT shift = current_offset + (event_offset - current_offset) / 2;
        std::cout << "safety shift by: " << shift << std::endl;
        PolyhedronTransformation::shiftFacetsInPlace(polyhedron, shift);

        // @debug
        {
            EdgeSPtr min_length_edge;
            CGAL::FT sq_min_edge_length = (std::numeric_limits<double>::max)();
            for (EdgeSPtr edge : polyhedron->edges()) {
                CGAL::FT sq_l = CGAL::squared_distance(*(edge->getVertexSrc()->getPoint()),
                                                       *(edge->getVertexDst()->getPoint()));
                if (sq_l < sq_min_edge_length) {
                    min_length_edge = edge;
                    sq_min_edge_length = sq_l;
                }
            }
            std::cout << "min edge length @ " << current_offset + shift << " = " << CGAL::approximate_sqrt(sq_min_edge_length) << std::endl;
            std::cout << min_length_edge->toString() << std::endl;
        }

        db::_3d::OBJFile::save("results/interrupted-shifted.obj", polyhedron, false /*do not triangulate*/);

        polyhedron_ = polyhedron;
        throw Far_enough_event();
    }

    void on_save_offset_event(PolyhedronSPtr polyhedron, CGAL::FT offset) override { }

    void after_offset_event(PolyhedronSPtr polyhedron, CGAL::FT offset) override { }

public:
    const CGAL::FT min_event_distance_;
    PolyhedronSPtr polyhedron_;
};





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //





std::stringstream ss_filename;
ss_filename << "results/" << save_path_.string() << "/face_count.txt";
std::ofstream out(ss_filename.str(), std::ios::app);
if (out) {
    out.precision(17);
    out << polyhedron->facets().size() << "\n";
    out.close();
}





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //






function process_file {
  FILE=$1
  BASE_NAME=$(basename "$FILE" | sed 's/\.[^.]*$//')

  # Extract the number between "input_" and the next underscore
  FULL_ID=$(echo ${BASE_NAME} | sed -E 's/^[^_]*_(.*\..*|.*)/\1/') # input_nnnnnn_n.ply to nnnnnn_n
  ID=$(echo "${BASE_NAME}" | sed 's/^[^_]*_\([^_]*\)_.*$/\1/') # input_nnnnnn_n.ply to nnnnnn, for the offset file

  mkdir -p $OUTPUT_DIRECTORY/${FULL_ID}

  # to write the result, concatenated later into the arrays
  RESULT_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/result.txt

  LOCAL_LOG=${OUTPUT_DIRECTORY}/${FULL_ID}/local_log.txt
  exec > ${LOCAL_LOG} 2>&1

  echo "FILE: $FILE" | tee -a "${LOCAL_LOG}"

  echo "BASE_NAME: ${BASE_NAME}"
  echo "Extracted FULL_ID: $FULL_ID"
  echo "Extracted ID: $ID"

  # ----------------------------------------------------------------
  # Convert to weighted PLY

  WEIGHT_FILE="${DATA_PATH}/offsets_${ID}.txt"
  echo "WEIGHT_FILE: ${WEIGHT_FILE}"

  WEIGHTED_INPUT="$OUTPUT_DIRECTORY/${FULL_ID}/input.ply"
  echo "WEIGHTED_INPUT: ${WEIGHTED_INPUT}"

  CMD="./convert_to_weighted_PLY ${FILE} ${WEIGHTED_INPUT} ${WEIGHT_FILE}"
  LOG_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/log_conversion.txt

  if [ ! -f $WEIGHT_FILE ]; then
    echo "====== [ERROR]: missing offset file?! ======"
    echo "$CMD" > $RESULT_FILE
    echo "PREPROCESS FAILURE" >> $RESULT_FILE
    return
  fi

  echo "---- Calling:"
  echo "  $CMD"
  echo "  $LOG_FILE"

  # Run conversion code
  timeout --preserve-status $TIMEOUT_VALUE $CMD > "$LOG_FILE" 2>&1
  RES=$?
  if [ ! "$RES" -eq 0 ]; then
    echo "====== [ERROR]: failed to create weighted PLY?! ======"
    echo "$CMD" > $RESULT_FILE
    echo "PREPROCESS FAILURE" >> $RESULT_FILE
    return
  fi

  # ----------------------------------------------------------------
  # If outward offset, add bbox and invert
  if [ "${OFFSET_DIRECTION}" == "out" ]; then
    WEIGHTED_INPUT_WITH_BBOX="$OUTPUT_DIRECTORY/${FULL_ID}/input_inverted_with_bbox.ply"
    echo "WEIGHTED_INPUT_WITH_BBOX: ${WEIGHTED_INPUT_WITH_BBOX}"

    CMD="./add_or_remove_bbox ${WEIGHTED_INPUT} add ${WEIGHTED_INPUT_WITH_BBOX}"
    LOG_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/log_add_bbox.txt

    echo "---- Calling:"
    echo "  $CMD"
    echo "  $LOG_FILE"

    # Run code
    timeout --preserve-status $TIMEOUT_VALUE $CMD > "$LOG_FILE" 2>&1
    RES=$?
    if [ ! "$RES" -eq 0 ]; then
      echo "====== [ERROR]: failed to preprocess ======"
      echo "$CMD" > $RESULT_FILE
      echo "PREPROCESS FAILURE" >> $RESULT_FILE
      return
    fi

    WEIGHTED_INPUT=${WEIGHTED_INPUT_WITH_BBOX}
  fi

  # ----------------------------------------------------------------
  # Run skeleton and compare results

  # below will write offset_-1.obj and offset_-1_exact.obj into the item's folder
  CMD="./StraightSkel 3d load ${WEIGHTED_INPUT} --no-window --save-offsets -1 --save-path ${OUTPUT_DIRECTORY}/${FULL_ID}"
  LOG_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/log_offset.txt

  echo "---- Calling:"
  echo "  $CMD"
  echo "  $LOG_FILE"

  TIMING_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/timing.txt
  echo "$TIMING_FILE"

  # Run offset code
  RES_CODE_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/res_code.txt
  time ( timeout --preserve-status $TIMEOUT_VALUE $CMD > "$LOG_FILE" 2>&1; echo $? > "$RES_CODE_FILE" ) 2> "$TIMING_FILE"

  # Extract runtime from the output of timeout
  runtime=$(cat "$TIMING_FILE" | grep "real" | awk '{print $2}')

  # Convert runtime to seconds
  if [[ $runtime =~ m ]]; then
    minutes=$(echo "$runtime" | sed 's/m.*//')
    seconds=$(echo "$runtime" | sed 's/.*m//;s/s//')
    runtime=$(echo "$minutes * 60 + $seconds" | bc)
  else
    runtime=$(echo "$runtime" | sed 's/s//')
  fi

  RUNTIME_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/runtime.txt
  echo $runtime > ${RUNTIME_FILE}

  RES=$(cat "$RES_CODE_FILE")
  if [ "$RES" -eq 0 ]; then
    # ----------------------------------------------------------------
    # Got a result, compare it

    OC_OUTPUT=${FILE/input/output}
    echo "OC_OUTPUT: ${OC_OUTPUT}"

    if [ ! -f $OC_OUTPUT ]; then
      echo "[ERROR]: missing OC output file?!"
      echo "$CMD" > $RESULT_FILE
      echo "PREPROCESS FAILURE" >> $RESULT_FILE
      return
    fi

    cp ${OC_OUTPUT} ${OUTPUT_DIRECTORY}/${FULL_ID}

    OUTPUT=${OUTPUT_DIRECTORY}/${FULL_ID}/result.obj

    # ----------------------------------------------------------------
    # If outward offset, add Bbox and invert

    if [ "${OFFSET_DIRECTION}" == "out" ]; then
      OUTPUT_WITH_BBOX=${OUTPUT_DIRECTORY}/${FULL_ID}/offset_-1_exact.obj

      CMD="./add_or_remove_bbox ${OUTPUT_WITH_BBOX} remove ${OUTPUT}"
      LOG_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/log_remove_bbox.txt

      echo "---- Calling:"
      echo "  $CMD"
      echo "  $LOG_FILE"

      # Run code
      timeout --preserve-status $TIMEOUT_VALUE $CMD > "$LOG_FILE" 2>&1
      RES=$?
      if [ ! "$RES" -eq 0 ]; then
        echo "====== [ERROR]: failed postprocess ======"
        echo "$CMD" > $RESULT_FILE
        echo "POSTPROCESS FAILURE" >> $RESULT_FILE
        return
      fi
    fi

    # ----------------------------------------------------------------
    # Compare

    CMD="./compare_outputs ${OUTPUT} ${OC_OUTPUT}"
    LOG_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/log_comparison.txt

    echo "---- Calling:"
    echo "  $CMD"
    echo "  $LOG_FILE"

    timeout --preserve-status $TIMEOUT_VALUE $CMD > "$LOG_FILE" 2>&1
    RES=$?
    if [ "$RES" -eq 0 ]; then
      echo "====== [SUCCESS] ======"
      echo "$CMD" > $RESULT_FILE
      echo "SUCCESS" >> $RESULT_FILE
    else
      echo "====== [ERROR]: outputs differ ======"
      echo "$CMD" > $RESULT_FILE
      echo "BAD OUTPUT" >> $RESULT_FILE
    fi
  elif [ "$RES" -eq 143 ]; then
    echo "====== [ERROR]: time out! ======"
    echo "$CMD" > $RESULT_FILE
    echo "TIMEOUT" >> $RESULT_FILE
  elif [ "$RES" -eq 137 ]; then
    echo "====== [ERROR]: process killed! ======"
    echo "$CMD" > $RESULT_FILE
    echo "OTHER FAILURE" >> $RESULT_FILE
  else
    echo "====== [ERROR]: offset failure! ======"
    echo "$CMD" > $RESULT_FILE
    echo "OFFSET FAILURE" >> $RESULT_FILE
  fi
  echo ""
}





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //





void PolyhedronTransformation::randTiltPlanes(PolyhedronSPtr polyhedron) {
    // Nudge vertices.
    // If we only nudged planes with fixed point constraints, we might not ensure generic position,
    // for example if two pairs of constraints are along the same line.
    //
    // @todo could restrict to only high degree vertices in facets that have 2 high degree vertices
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr p = vertex->getPoint();
        std::array<double, 3> v_r = randVec(-1e-15, 1e-15); // @fixme hardcoded values
        Point3SPtr p_t = KernelFactory::createPoint3(CGAL::to_double(p->x()) + v_r[0],
                                                     CGAL::to_double(p->y()) + v_r[1],
                                                     CGAL::to_double(p->z()) + v_r[2]);
        vertex->setPoint(p_t);
    }

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        facet->perturbPlaneCoefficientsNudge(1e-10);
    }
}



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


for (FacetSPtr f : polyhedron->facets()) {
    if (!f->isTriangle() || !has_high_degree_vertices(f)) {
        continue;
    }

    std::cout << "Fix triangle F" << f->getID() << std::endl;

    // Triangle facets are great because we have freedom to move points and compute planes afterwards
    // Check at all (3) vertices if we still have freedom of movement (i.e., they have fewer
    // than 3 determining facets)
    for (VertexSPtr v : f->vertices()) {
        std::cout << "Recompute position of triangle vertex V" << v->getID() << std::endl;
        nudge_constrained_vertex(v);
        fixing_vertices[f].insert(v);
    }

    f->initPlane();
    f->normalizePlaneCoefficients();

    dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_triangle.OFF", f);

    // Here we need to also need to update the determining facets because some neighboring
    // facets could be unfixed high degree triangles
    for (VertexSPtr v : f->vertices()) {
        if (determining_facets[v].size() < 3) {
            determining_facets[v].insert(f);
            std::cout << "  V" << v->getID() << " is determined by " << f->getID() << std::endl;
            // no need to cascade here, we know only triangles are left
        }
    }
}





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //





/**
 * Collects the facets that are affected by the event.
 */
static std::list<FacetSPtr> eventFacets(AbstractEventSPtr event);


// these are facets that will be modified by the event
// @todo should this list be stored in the event? and then -> isvalid for all + check is_same
// but then: « sure, your faces are the valid & the same, but what about edges|vertices? »...
// so, should this be "eventPolyhedron" that extracts the whole subset of the polyhedron
// involved in the event?
std::list<FacetSPtr> SimpleStraightSkel::eventFacets(AbstractEventSPtr event)
{
    // std::cout << "eventFacets()" << std::endl;
    CGAL_precondition(event->isValid());

    std::list<FacetSPtr> result;

    if (event->getType() == AbstractEvent::EDGE_EVENT) {
        EdgeEventSPtr edge_event = std::dynamic_pointer_cast<EdgeEvent>(event);
        EdgeSPtr edge = edge_event->getEdge();
        DEBUG_SPTR(edge);
        result.push_back(edge->getFacetR());
        result.push_back(edge->getFacetL());
        result.push_back(edge->getFacetSrc());
        result.push_back(edge->getFacetDst());
        CGAL_postcondition(std::distance(result.begin(), std::unique(result.begin(), result.end())) == 4);
    } else if (event->getType() == AbstractEvent::EDGE_MERGE_EVENT) {
        EdgeMergeEventSPtr edge_merge_event = std::dynamic_pointer_cast<EdgeMergeEvent>(event);
        EdgeSPtr edge = edge_merge_event->getEdge1()->next(edge_merge_event->getEdge1()->getFacetL());
        DEBUG_SPTR(edge);
        result.push_back(edge->getFacetR());
        result.push_back(edge->getFacetL());
        result.push_back(edge->getFacetSrc());
        result.push_back(edge->getFacetDst());
        CGAL_postcondition(std::distance(result.begin(), std::unique(result.begin(), result.end())) == 4);
    } else if (event->getType() == AbstractEvent::TRIANGLE_EVENT) {
        TriangleEventSPtr triangle_event = std::dynamic_pointer_cast<TriangleEvent>(event);
        EdgeSPtr edge = triangle_event->getEdgeBegin();
        DEBUG_SPTR(edge);
        result.push_back(edge->getFacetR());
        result.push_back(edge->getFacetL());
        result.push_back(edge->getFacetSrc());
        result.push_back(edge->getFacetDst());
        CGAL_postcondition(std::distance(result.begin(), std::unique(result.begin(), result.end())) == 4);
    } else if (event->getType() == AbstractEvent::DBL_EDGE_MERGE_EVENT) {
        DblEdgeMergeEventSPtr dbl_edge_merge_event = std::dynamic_pointer_cast<DblEdgeMergeEvent>(event);
        // could be a L-R-Src-Dst too, of e.g. Edge11->next(Edge11->getFacetL)
        result.push_back(dbl_edge_merge_event->getEdge11()->getFacetL());
        result.push_back(dbl_edge_merge_event->getEdge11()->getFacetR());
        result.push_back(dbl_edge_merge_event->getEdge22()->getFacetL());
        result.push_back(dbl_edge_merge_event->getEdge22()->getFacetR());
        CGAL_postcondition(std::set<FacetSPtr>(result.begin(), result.end()).size() == 4);
    } else if (event->getType() == AbstractEvent::DBL_TRIANGLE_EVENT) {
        DblTriangleEventSPtr dlb_triangle_event = std::dynamic_pointer_cast<DblTriangleEvent>(event);
        EdgeSPtr edge = dlb_triangle_event->getEdge();
        result.push_back(edge->getFacetR());
        result.push_back(edge->getFacetL());
        result.push_back(edge->getFacetSrc());
        result.push_back(edge->getFacetDst());
        CGAL_postcondition(std::set<FacetSPtr>(result.begin(), result.end()).size() == 4);
    } else if (event->getType() == AbstractEvent::TETRAHEDRON_EVENT) {
        TetrahedronEventSPtr tetrahedron_event = std::dynamic_pointer_cast<TetrahedronEvent>(event);
        EdgeSPtr edge = tetrahedron_event->getEdgeBegin();
        result.push_back(edge->getFacetL());
        result.push_back(edge->getFacetR());
        result.push_back(edge->getFacetSrc());
        result.push_back(edge->getFacetDst());
        CGAL_postcondition(std::set<FacetSPtr>(result.begin(), result.end()).size() == 4);
    } else if (event->getType() == AbstractEvent::VERTEX_EVENT) {
        VertexEventSPtr vertex_event = std::dynamic_pointer_cast<VertexEvent>(event);
        for (FacetWPtr wf : vertex_event->getVertex1()->facets()) {
            result.push_back(wf.lock());
        }
        for (FacetWPtr wf : vertex_event->getVertex2()->facets()) {
            result.push_back(wf.lock());
        }
        // CGAL_postcondition(std::set<FacetSPtr>(result.begin(), result.end()).size() == 4);
    } else if (event->getType() == AbstractEvent::FLIP_VERTEX_EVENT) {
        FlipVertexEventSPtr flip_vertex_event = std::dynamic_pointer_cast<FlipVertexEvent>(event);
        for (FacetWPtr wf : flip_vertex_event->getVertex1()->facets()) {
            result.push_back(wf.lock());
        }
        for (FacetWPtr wf : flip_vertex_event->getVertex2()->facets()) {
            result.push_back(wf.lock());
        }
        // CGAL_postcondition(std::set<FacetSPtr>(result.begin(), result.end()).size() == 4);
    } else if (event->getType() == AbstractEvent::SURFACE_EVENT) {
        SurfaceEventSPtr surface_event = std::dynamic_pointer_cast<SurfaceEvent>(event);
        EdgeSPtr edge_1 = surface_event->getEdge1();
        EdgeSPtr edge_2 = surface_event->getEdge2();
        result.push_back(edge_1->getFacetL());
        result.push_back(edge_1->getFacetR());
        result.push_back(edge_2->getFacetL());
        result.push_back(edge_2->getFacetR());
        // CGAL_postcondition(std::set<FacetSPtr>(result.begin(), result.end()).size() == 4);
    } else if (event->getType() == AbstractEvent::POLYHEDRON_SPLIT_EVENT) {
        PolyhedronSplitEventSPtr polyhedron_split_event =
            std::dynamic_pointer_cast<PolyhedronSplitEvent>(event);
        EdgeSPtr edge_1 = polyhedron_split_event->getEdge1();
        EdgeSPtr edge_2 = polyhedron_split_event->getEdge2();
        result.push_back(edge_1->getFacetL());
        result.push_back(edge_1->getFacetR());
        result.push_back(edge_2->getFacetL());
        result.push_back(edge_2->getFacetR());
        // CGAL_postcondition(std::set<FacetSPtr>(result.begin(), result.end()).size() == 4);
    } else if (event->getType() == AbstractEvent::SPLIT_MERGE_EVENT) {
        SplitMergeEventSPtr split_merge_event = std::dynamic_pointer_cast<SplitMergeEvent>(event);
        for (FacetWPtr wf : split_merge_event->getVertex1()->facets()) {
            result.push_back(wf.lock());
        }
        for (FacetWPtr wf : split_merge_event->getVertex2()->facets()) {
            result.push_back(wf.lock());
        }
        // CGAL_postcondition(std::set<FacetSPtr>(result.begin(), result.end()).size() == 4);
    } else if (event->getType() == AbstractEvent::EDGE_SPLIT_EVENT) {
        EdgeSplitEventSPtr edge_split_event = std::dynamic_pointer_cast<EdgeSplitEvent>(event);
        EdgeSPtr edge_1 = edge_split_event->getEdge1();
        EdgeSPtr edge_2 = edge_split_event->getEdge2();
        result.push_back(edge_1->getFacetL());
        result.push_back(edge_1->getFacetR());
        result.push_back(edge_2->getFacetL());
        result.push_back(edge_2->getFacetR());
        // CGAL_postcondition(std::set<FacetSPtr>(result.begin(), result.end()).size() == 4);
    } else if (event->getType() == AbstractEvent::PIERCE_EVENT) {
        PierceEventSPtr pierce_event = std::dynamic_pointer_cast<PierceEvent>(event);
        result.push_back(pierce_event->getFacet());
        for (FacetWPtr wf : pierce_event->getVertex()->facets()) {
            result.push_back(wf.lock());
        }
        // CGAL_postcondition(std::set<FacetSPtr>(result.begin(), result.end()).size() == 4);
    } else {
        std::cerr << "Error: unknown event type" << std::endl;
        CGAL_assertion(false);
    }

    return result;
}

/**
 * Returns `true` if the event has facets that have a "step_id" that is
 * greater than the event's "step_id".
 */
static bool isEventPotentiallyObsolete(AbstractEventSPtr event);

// First events are initialized with step_ID '-1'.  At step #i, events' step_ID is 'i'.
// When not refreshing the queue, an event is potentially obsolete
// if its step ID is (strictly) smaller than that of any of its faces
bool SimpleStraightSkel::isEventPotentiallyObsolete(AbstractEventSPtr event)
{
    CGAL_precondition(event->isValid());

    // if any of the facets involved in the event has been "touched"
    // after the event has been created, then the event is tagged as obsolete
    std::list<FacetSPtr> facets = eventFacets(event);
    for (FacetSPtr facet : facets) {
        int event_step_id = event->getStepID();
        int facet_step_id = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getStepID();
        if (facet_step_id > event_step_id) {
            // std::cout << "Event is potentially obsolete:\n";
            // std::cout << " --> event " << event->getID() << " has stamp " << event_step_id << "\n";
            // std::cout << " --> facet " << facet->getID() << " has stamp " << facet_step_id << std::endl;
            return true;
        }
    }

    return false;
}





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //




enum class Triangulation_strategy
{
    DEFAULT = 0,
    MID_CUT,
    EQUALIZE_HDV
};


std::list<FacetSPtr> PolyhedronTransformation::triangulate(FacetSPtr facet,
                                                           PolyhedronSPtr polyhedron,
                                                           Triangulation_strategy strategy) {
    std::list<FacetSPtr> created_facets;
    CGAL_precondition(facet && polyhedron && facet->vertices().size() >= 3);

    if (facet->vertices().size() == 3) {
        return { facet };
    }

    facet->sortVertices();

    // Prepare containers for triangulation
    auto pcdt = triangulate_facet_with_CDT2(facet);

    using PCDT = decltype(pcdt);
    using PCDT_VH = typename PCDT::Vertex_handle;
    using PCDT_FH = typename PCDT::Face_handle;

    std::unordered_map<PCDT_FH, bool> in_domain_map;
    boost::associative_property_map<std::unordered_map<PCDT_FH, bool>> in_domain(in_domain_map);
    CGAL::mark_domain_in_triangulation(pcdt, in_domain);

    // @debug
    unsigned int tr_n = 0;
    for (auto f : pcdt.finite_face_handles()) {
        if (get(in_domain, f)) {
            ++tr_n;
        }
    }
    std::cout << tr_n << " triangulation faces to partition" << std::endl;

   // Get the speed from the parent facet (if any)
    CGAL::FT parent_speed = 1.0;
    if (facet->hasData()) {
        SkelFacetDataSPtr parent_data = std::dynamic_pointer_cast<SkelFacetData>(facet->getData());
        if (parent_data) {
            parent_speed = parent_data->getSpeed();
        }
    }

    if (strategy == Triangulation_strategy::DEFAULT)
    {
        polyhedron->removeFacet(facet);

        // Create new facets and edges for each triangle
        for(auto fh : pcdt.finite_face_handles()) {
            if(!get(in_domain, fh)) {
                continue;
            }

            VertexSPtr v0 = fh->vertex(0)->info();
            VertexSPtr v1 = fh->vertex(1)->info();
            VertexSPtr v2 = fh->vertex(2)->info();
            VertexSPtr verts[3] = {v0, v1, v2};
            FacetSPtr new_facet = Facet::create(3, verts);
            Plane3SPtr plane = KernelFactory::createPlane3(v0->getPoint(),
                                                           v1->getPoint(),
                                                           v2->getPoint());
            new_facet->setPlane(plane);
            new_facet->normalizePlaneCoefficients();
            SkelFacetDataSPtr new_data = SkelFacetData::create(new_facet);
            new_data->setSpeed(parent_speed);
            polyhedron->addFacet(new_facet);
            created_facets.push_back(new_facet);
        }
    } else if (strategy == Triangulation_strategy::MID_CUT) {
        // Find all candidate internal edges (not constraints)

        struct Candidate {
            PCDT_VH a, b;
            double min_area;
        };

        Candidate best;
        best.min_area = -1.0;
        for (auto eit = pcdt.finite_edges_begin(); eit != pcdt.finite_edges_end(); ++eit) {
            if (pcdt.is_constrained(*eit)) {
                continue;
            }
            auto fh = eit->first;
            int i = eit->second;
            if (!get(in_domain, fh)) {
                continue;
            }

            PCDT_VH va = fh->vertex(pcdt.cw(i));
            PCDT_VH vb = fh->vertex(pcdt.ccw(i));
            std::cout << "test " << va->point() << " " << vb->point() << std::endl;

            // BFS on faces to split triangles into two groups
            std::set<PCDT_FH> group1, group2;
            std::map<PCDT_FH, bool> visited;
            // Mark the edge as cut
            std::queue<PCDT_FH> q;
            q.push(fh);
            visited[fh] = true;
            while (!q.empty()) {
                PCDT_FH cur = q.front();
                q.pop();
                CGAL_assertion(!pcdt.is_infinite(cur) && get(in_domain, cur));
                group1.insert(cur);
                for (int j = 0; j < 3; ++j) {
                    if ((cur == fh && (j == i)) ||
                         pcdt.is_constrained(std::make_pair(cur, j))) {
                        continue;
                    }

                    PCDT_FH neigh = cur->neighbor(j);
                    CGAL_assertion(!pcdt.is_infinite(neigh));

                    if (!visited[neigh]) {
                        visited[neigh] = true;
                        q.push(neigh);
                    }
                }
            }
            // The rest go to group2
            for (auto f : pcdt.finite_face_handles()) {
                if (get(in_domain, f) && !visited[f]) {
                    group2.insert(f);
                }
            }

            std::cout << "Group sizes: " << group1.size() << " " << group2.size() << std::endl;
            CGAL_assertion(group1.size() + group2.size() == tr_n);
            CGAL_assertion(!group1.empty());
            CGAL_assertion(!group2.empty());

            // Compute area for each group (sum triangle areas)
            auto group_area = [&](const std::set<PCDT_FH>& group) {
                double area = 0.0;
                for (auto f : group) {
                    const auto& p0 = *(f->vertex(0)->info()->getPoint());
                    const auto& p1 = *(f->vertex(1)->info()->getPoint());
                    const auto& p2 = *(f->vertex(2)->info()->getPoint());
                    area += std::abs(CGAL::to_double(CGAL::approximate_sqrt(CGAL::squared_area(p0, p1, p2))));
                }
                return area;
            };

            double area1 = group_area(group1);
            double area2 = group_area(group2);
            double min_area = (std::min)(area1, area2);

            // @todo give absolute priority if the edge cut uses already fixed vertices
            // @todo if priority if the edge uses high degree vertices
            // give priority to the edge that equalizes areas
            if (min_area > best.min_area) {
                best = {va, vb, min_area};
                std::cout << "new best " << va->point() << " " << vb->point() << std::endl;
            }
        }

        // Collect unique vertices for each subpart:
        // walk facet->vertices() and fill verts1 and verts2
        // e.g. '0 1 2 3 4 5 6' with '3 5' being the cut edge,
        // then the facets should be split into
        // verts1 = {0, 1, 2, 3, 5, 6}
        // verts2 = {5, 3, 4}
        std::vector<VertexSPtr> verts1, verts2;

        // Loop over all vertices in the facet
        bool inserting_into_verts1 = true;

        std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
        while (it_v != facet->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            if (vertex == best.a->info() || vertex == best.b->info()) {
                verts1.push_back(vertex);
                verts2.push_back(vertex);
                std::cout << "inserting into both: " << vertex->getID() << std::endl;
                inserting_into_verts1 = !inserting_into_verts1;
            } else {
                if (inserting_into_verts1) {
                    std::cout << "inserting into verts1: " << vertex->getID() << std::endl;
                    verts1.push_back(vertex);
                } else {
                    std::cout << "inserting into verts2: " << vertex->getID() << std::endl;
                    verts2.push_back(vertex);
                }
            }
        }

        std::cout << "Verts sizes: " << verts1.size() << " " << verts2.size() << std::endl;
        CGAL_assertion(verts1.size() >= 3 && verts2.size() >= 3);
        CGAL_assertion(verts1.size() + verts2.size() == facet->vertices().size() + 2);

        polyhedron->removeFacet(facet);

        FacetSPtr new_facet1 = Facet::create(verts1.size(), verts1.data());
        FacetSPtr new_facet2 = Facet::create(verts2.size(), verts2.data());
        new_facet1->setPlane(facet->plane());
        new_facet2->setPlane(facet->plane());
        new_facet1->normalizePlaneCoefficients();
        new_facet2->normalizePlaneCoefficients();
        SkelFacetDataSPtr new_data1 = SkelFacetData::create(new_facet1);
        SkelFacetDataSPtr new_data2 = SkelFacetData::create(new_facet2);
        new_data1->setSpeed(parent_speed);
        new_data2->setSpeed(parent_speed);

        polyhedron->addFacet(new_facet1);
        polyhedron->addFacet(new_facet2);

        created_facets.push_back(new_facet1);
        created_facets.push_back(new_facet2);
    } else if (strategy == Triangulation_strategy::EQUALIZE_HDV) {

    }

    polyhedron->initializeAllIDs();

    return created_facets;
}




// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


// This particular reader can read weights if they are stored as a property in the PLY file
PolyhedronSPtr PLYFile::load(const std::string& filename) {
    PolyhedronSPtr result = PolyhedronSPtr();

    std::ifstream ifs(filename.c_str());
    if (ifs.is_open()) {
        // @todo avoid building the intermediate mesh (but still use CGAL readers)
        CGAL::Surface_mesh<Point3> sm;
        bool success = CGAL::IO::read_PLY(ifs, sm);
        if (!success) {
          return {};
        }

        result = Surface_meshIO::load(sm);
    }

    return result;
}


// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


static PolyhedronSPtr merge_and_perturb(PolyhedronSPtr polyhedron);


PolyhedronSPtr PolyhedronTransformation::merge_and_perturb(PolyhedronSPtr polyhedron) {
  // Check if we can tilt facets' planes (i.e., nudge plane coefficients) directly.
  // A sufficient condition is that all vertices have degree 3: in that case, a small tilt
  // of the plane will still yield a single intersection point.
  // That's not the case (in general) for degree > 3 vertices as there would no longer be
  // a single intersection point for the tilted planes.
  //
  // The advantage is that we can manipulate much smaller meshes since the facets are polygonal.
  bool canUsePlaneTilts;

  // copy the polyhedron because we will merge (almost) coplanar facets and check if the result
  // is a mesh with only degree 3 vertices.
  PolyhedronSPtr polyhedron_cpy = polyhedron->clone();

  // @todo?
  // could we merge non-connected input facets as to assign them the same (tilted) plane?
  // Often in inputs we have many facets that correspond to the same plane, but vertical facets
  // split it into separate connected components.
  // The important thing is that we don't want to create degenerate conditions so the CCs
  // should NOT interact with each other; how to prevent that?...
  db::_3d::AbstractFile::mergeCoplanarFacets(polyhedron_cpy);

  CGAL_assertion(polyhedron_cpy && polyhedron_cpy->isConsistent());

  PolyhedronTransformation::normalizeFacetPlanes(polyhedron_cpy);

  canUsePlaneTilts = PolyhedronTransformation::isTiltCompatible(polyhedron_cpy);
  DEBUG_PRINT("Tiltability: " << canUsePlaneTilts);

  if (canUsePlaneTilts) {
      polyhedron = polyhedron_cpy;
      PolyhedronTransformation::randTiltPlanes(polyhedron);
      resetPoints(polyhedron);
  } else {
      // this is not 'polyhedron_cpy' because the polyhedron must be triangulated
      // for vertices to remain on the planes of their incident facets
      randMovePoints(polyhedron);
  }

  DEBUG_PRINT("Done with perturbation");

  return polyhedron;
}




// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



namespace PMP = CGAL::Polygon_mesh_processing;

CGAL::Bbox_3 bbox = PMP::bbox(sm);
const CGAL::FT diag_length = CGAL::approximate_sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                                    CGAL::square(bbox.ymax() - bbox.ymin()) +
                                                    CGAL::square(bbox.zmax() - bbox.zmin()));

// Use shape detection to analyze the mesh
std::vector<std::size_t> region_ids(num_faces(sm));
boost::vector_property_map<Plane3> plane_map; // supporting planes of the regions detected

const CGAL::FT cos_of_max_angle = 0.98;
const CGAL::FT max_distance = 0.0001 * diag_length;

// detect planar regions in the mesh
// @todo growing should:
// - use the .ini value of 'epsilon_coplanarity'
// - stop if it merges faces with different weights
// - give an error for adjacent coplanar faces that have different weights
std::size_t nb_regions =
    PMP::region_growing_of_planes_on_faces(sm,
                                            CGAL::make_random_access_property_map(region_ids),
                                            CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle)
                                                            .region_primitive_map(plane_map)
                                                            .maximum_distance(max_distance));

static int region_dump_id = -1;
utils::save_colored_mesh(sm, region_ids, "results/regions_" + std::to_string(++region_dump_id) + ".ply");

// detect corner vertices on the boundary of planar regions
std::vector<std::size_t> corner_ids(num_vertices(sm), -1); // corner status of vertices
std::vector<bool> ecm(num_edges(sm), false); // mark edges at the boundary of regions

std::size_t nb_corners =
    PMP::detect_corners_of_regions(sm,
                                  CGAL::make_random_access_property_map(region_ids),
                                  nb_regions,
                                  CGAL::make_random_access_property_map(corner_ids),
                                  CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle).
                                                    maximum_distance(max_distance).
                                                    edge_is_constrained_map(CGAL::make_random_access_property_map(ecm)));

// @debug
{
    for (face_descriptor f : faces(sm)) {
        std::cout << "face " << f << " is in region " << region_ids[f] << std::endl;
    }
}

// the almost-coplanar merge is performed after the conversion to the Polyhedron
// data structure because we want to be able to create faces that have holes,
// which the CGAL::Surface_mesh class does not support
std::map<edge_descriptor, EdgeWPtr> e2e;
PolyhedronSPtr polyhedron = db::_3d::Surface_meshIO::load(sm, e2e, np);

// If regions are tilt-compatible, i.e. they have at most 2 high degree vertices, merge the facets
bool should_merge = true;
std::vector<unsigned int> high_degree_corners_n(nb_regions, 0); // region ID -> number of high degree corners
vertex_iterator vit = vertices(sm).begin(), vend = vertices(sm).end();
for (; vit!=vend; ++vit) {
    if (corner_ids[*vit] == static_cast<std::size_t>(-1)) {
        continue;
    }

    std::set<std::size_t> incident_regions;
    for (face_descriptor f : CGAL::faces_around_target(halfedge(*vit, sm), sm)) {
        incident_regions.insert(region_ids[f]);
    }

    // more than 3 incident regions ==> high degree vertex
    if (incident_regions.size() > 3) {
        DEBUG_PRINT("Corner with high degree " << incident_regions.size() << " at " << sm.point(*vit));
        for (std::size_t ri : incident_regions) {
            ++(high_degree_corners_n[ri]);
        }
    }
}

// more than 2 ==> we cannot constrain
for (unsigned int n : high_degree_corners_n) {
    if (n > 2) {
        should_merge = false;
    }
}

DEBUG_PRINT("should_merge = " << should_merge);
if (should_merge) {
    // merge the facets incident to an unconstrained edge (i.e., the edge is interior to a region)
    for (edge_descriptor e: edges(sm)) {
        if (ecm[e]) {
            continue;
        }

        EdgeSPtr edge = e2e[e].lock();
        if (!edge) {
            continue;
        }

        DEBUG_PRINT("Merging facets " << edge->getFacetL()->getID() << " and " << edge->getFacetR()->getID());
        CGAL_assertion(sm.point(source(e, sm)) == *(edge->getVertexSrc()->getPoint()));
        CGAL_assertion(sm.point(target(e, sm)) == *(edge->getVertexDst()->getPoint()));

        // @todo it seems like intermediate states are somewhat unsound during edge merging
        db::_3d::AbstractFile::mergeFacets(edge, polyhedron);
    }

    polyhedron->initializeAllIDs();
}

DEBUG_PRINT("Sanitizing...");
db::_3d::AbstractFile::sanitize(polyhedron);



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



template <typename NamedParameters>
bool Surface_meshIO::save(const PolyhedronSPtr& polyhedron,
                          CGAL::Surface_mesh<Point3>& sm,
                          const NamedParameters& np)
{
    using CGAL::parameters::choose_parameter;
    using CGAL::parameters::get_parameter;

   // @fixme factorize with OBJ/PLY save()
    using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;
    using PK = CGAL::Projection_traits_3<CGAL::K>;
    using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>;
    using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
    using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
    using PTDS = CGAL::Triangulation_data_structure_2<PVb, PFb>;
    using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
    using PCDT_VH = PCDT::Vertex_handle;
    using PCDT_FH = PCDT::Face_handle;

    std::cout << "--> SM of polyhedron with " << polyhedron->vertices().size() << " vertices and "
              << polyhedron->facets().size() << " facets" << std::endl;

    bool do_triangulate = !choose_parameter(get_parameter(np, CGAL::internal_np::do_not_triangulate_faces), false);

    // Vertices
    std::unordered_map<VertexSPtr, vertex_descriptor> v_map;
    for (VertexSPtr v : polyhedron->vertices()) {
        vertex_descriptor idx = sm.add_vertex(*(v->getPoint()));
        v_map[v] = idx;
    }

    // Write facets
    auto weight_pmap = choose_parameter(get_parameter(np, CGAL::internal_np::face_weight),
                                        CGAL::Constant_property_map<std::size_t, double>(1.0));

    while (FacetSPtr facet : polyhedron->facets()) {
        CGAL_assertion(facet->getID() != -1);
        // std::cout << "handle facet " << facet->getID() << std::endl;

        double speed = 1.;
        if (facet->hasData()) {
            auto skel_data = std::dynamic_pointer_cast<data::_3d::skel::SkelFacetData>(facet->getData());
            if (skel_data) {
                speed = CGAL::to_double(skel_data->getSpeed());
            }
        }

        bool do_triangulate_face = do_triangulate;
        if (facet->edges().size() < 3)
            do_triangulate_face = false;

        if (do_triangulate_face)
        {
            Vector3SPtr n = KernelFactory::createVector3(facet->plane());

            // @todo might have to do something fancier than staring from a single vertex
            // for degenerate faces with zigzagging edges...
            CGAL_assertion(*n != CGAL::NULL_VECTOR);

            PK traits(*n);
            PCDT pcdt(traits);

            std::map<VertexSPtr, PCDT_VH> face_vhs;

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
                    // std::cerr << "W: encountered degenerate edge @ " << *(v0->getPoint()) << std::endl;

                    CGAL_assertion(v0->degree() != 1);
                    VertexSPtr vm1 = edge->prev(facet)->src(facet);

                    face_descriptor sm_f = sm.add_face(v_map[vm1], v_map[v0], v_map[v1]);
                    if (sm_f == boost::graph_traits<CGAL::Surface_mesh<Point3>>::null_face()) {
                        std::cerr << "Error: failed to add face to surface mesh (1)" << std::endl;
                        return false;
                    }

                    put(weight_pmap, sm_f, speed);
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
                        DEBUG_PRINT("Error: Intersection of constraints");
                        DEBUG_PRINT("While inserting " << *(v0->getPoint()) << " || " << *(v1->getPoint()));
                        DEBUG_PRINT(facet->toString());
                        CGAL_warning_msg(false, "Intersections in CDT2 not allowed");
                        do_triangulate_face = false;
                        break;
                    }
                    ++ne;
                }
            }

            if(ne < 3) // degenerate face
            {
                std::cerr << "Warning: skipping degenerate face" << std::endl;
                continue;
            }

            if(do_triangulate_face)
            {
                std::unordered_map<PCDT_FH, bool> in_domain_map;
                boost::associative_property_map<std::unordered_map<PCDT_FH, bool>> in_domain(in_domain_map);

                CGAL::mark_domain_in_triangulation(pcdt, in_domain);

                for(auto fh : pcdt.finite_face_handles())
                {
                    if(!get(in_domain, fh))
                      continue;

                    bool canAddFace = CGAL::Euler::can_add_face(CGAL::make_array(v_map[fh->vertex(0)->info()],
                                                               v_map[fh->vertex(1)->info()],
                                                               v_map[fh->vertex(2)->info()]),
                                                               sm);
                    std::cout << "canAddFace: " << canAddFace << std::endl;

                    face_descriptor sm_f =
                        CGAL::Euler::add_face(CGAL::make_array(v_map[fh->vertex(0)->info()],
                                                               v_map[fh->vertex(1)->info()],
                                                               v_map[fh->vertex(2)->info()]),
                                                               sm);
                    if (sm_f == boost::graph_traits<CGAL::Surface_mesh<Point3>>::null_face()) {
                        std::cerr << "Error: failed to add face to surface mesh (2)" << std::endl;
                        std::cerr << "Face:\n" << fh->vertex(0)->point() << " "
                                  << fh->vertex(1)->point() << " "
                                  << fh->vertex(2)->point() << std::endl;
                        CGAL::IO::write_polygon_mesh("results/failed.off", sm, CGAL::parameters::stream_precision(17));
                        return false;
                    }

                    put(weight_pmap, sm_f, speed);
                }
            }
        }

        if(!do_triangulate_face)
        {
            std::set<EdgeSPtr> visited_edges;

            std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
            while (it_e != facet->edges().end()) {
                EdgeSPtr edge = *it_e++;
                if (visited_edges.find(edge) != visited_edges.end()) {
                    continue; // already visited
                }

                std::vector<vertex_descriptor> boundary_vertices;
                EdgeSPtr start_edge = edge;
                // std::cout << "start a walk @ " << start_edge->src(facet)->getID() << " " << start_edge->dst(facet)->getID() << std::endl;
                bool is_open = false;

                // Walk forward to collect boundary vertices
                do {
                    visited_edges.insert(edge);
                    boundary_vertices.push_back(v_map[edge->src(facet)]);
                    if (edge->dst(facet)->degree() == 1) {
                        is_open = true;
                        break;
                    }
                    edge = edge->next(facet);
                } while (edge != start_edge);

                // If open, also walk backward to collect remaining boundary vertices
                if (is_open) {
                    boundary_vertices.push_back(v_map[edge->dst(facet)]);
                    edge = start_edge;
                    while (edge->src(facet)->degree() != 1) {
                        edge = edge->prev(facet);
                        visited_edges.insert(edge);
                        boundary_vertices.insert(boundary_vertices.begin(), v_map[edge->src(facet)]);
                    }
                }

                // Write the boundary as a face
                face_descriptor sm_f = sm.add_face(boundary_vertices);
                if (sm_f == boost::graph_traits<CGAL::Surface_mesh<Point3>>::null_face()) {
                    std::cerr << "Error: failed to add face to surface mesh (3)" << std::endl;
                    return false;
                }
                put(weight_pmap, sm_f, speed);
            }
        }
    }

    return true;
}




// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


static CGAL::Sign SideOfBisector(FacetSPtr facet, VertexSPtr vertex, Point3SPtr point);


// Returns:
// - ON_POSITIVE_SIDE  if the query is on the positive side (see below for definition)
//   of the planar bisector of the two edges incident to 'vertex'.
// - ON_ORIENTED_BOUNDARY if the query is on the bisector
// - ON_NEGATIVE_SIDE if the query is on the negative side of the bisector
//
// The positive side of the bisector is the side that has p_out
// in the following (facet seen from above):
//
//                         | <-- planar bisector (positive side is here)
//                         |
//           edge_in       |          edge_out
//  p_in --------------> vertex -------------------> p_out
//                         |
//                         |
//
// \pre `point` is coplanar with the edges incident to `vertex`
CGAL::Sign SelfIntersection::SideOfBisector(FacetSPtr facet, VertexSPtr vertex, Point3SPtr point) {
#ifndef USE_CGAL
# error "This function is not compatible with the old kernel"
#endif

    CGAL_precondition(vertex->degree() > 1);

    CGAL::Sign result;
    EdgeSPtr edge_in;
    EdgeSPtr edge_out;

    std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
    while (it_e != vertex->edges().end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge(edge_wptr);
            if (edge->src(facet) == vertex) {
                edge_out = edge;
            } else if (edge->dst(facet) == vertex) {
                edge_in = edge;
            }
        }
    }

    if (edge_in && edge_out) {
        // @fixme plenty of constructions below
        Point3SPtr p_mid = vertex->getPoint();
        Point3SPtr p_in = edge_in->src(facet)->getPoint();
        Point3SPtr p_out = edge_out->dst(facet)->getPoint();
        Vector3SPtr normal = KernelFactory::createVector3(facet->plane());
        Point3SPtr p_normal = KernelFactory::createPoint3(*point + *normal);

        CGAL::Orientation or_in = CGAL::orientation(*p_mid, *p_normal, *p_in, *point);
        CGAL::Orientation or_out = CGAL::orientation(*p_mid, *p_out, *p_normal, *point);

        // positive = right turn (reflex), null = collinear, negative = left turn (convex)
        CGAL::Orientation or_edge = CGAL::orientation(*p_mid, *p_normal, *p_in, *p_out);

        // @todo only compute this when necessary and it should be a predicate
        // Note that there are no weights involved here: we are checking for self-intersections
        // in a fixed polyhedron
        CGAL::FT sq_dist_in = KernelWrapper::squared_distance(edge_in->line(), point);
        CGAL::FT sq_dist_out = KernelWrapper::squared_distance(edge_out->line(), point);

        // DEBUG_PRINT("p_in = " << *p_in);
        // DEBUG_PRINT("p_mid = " << *p_mid);
        // DEBUG_PRINT("p_out = " << *p_out);
        // DEBUG_PRINT("query = " << *point);
        // DEBUG_PRINT("or_in = " << or_in);
        // DEBUG_PRINT("or_out = " << or_out);
        // DEBUG_PRINT("or_edge = " << or_edge);
        // DEBUG_PRINT("sq_dist_in = " << sq_dist_in);
        // DEBUG_PRINT("sq_dist_out = " << sq_dist_out);

        // that's for the left turn (convex); for reflex, it's the opposite
        if(or_in == CGAL::POSITIVE) {
            if(or_out == CGAL::POSITIVE) {
                // upper quadrant
                result = CGAL::compare(sq_dist_out, sq_dist_in);
            } else { // or_out != CGAL::POSITIVE
                // right quadrant
                if(*point == *p_mid) {
                  result = CGAL::ON_ORIENTED_BOUNDARY;
                } else {
                  result = CGAL::ON_NEGATIVE_SIDE;
                }
            }
        } else { // or_in != CGAL::POSITIVE
            if(or_out == CGAL::POSITIVE) {
                // left quadrant
                if(*point == *p_mid) {
                    result = CGAL::ON_ORIENTED_BOUNDARY;
                } else {
                    result = CGAL::ON_POSITIVE_SIDE;
                }
            } else { // or_out != CGAL::POSITIVE
                // bottom quadrant
                result = CGAL::compare(sq_dist_in, sq_dist_out);
            }
        }

        // opposite for reflex
        if (or_edge == CGAL::POSITIVE /*reflex*/) {
          result = (result == CGAL::POSITIVE) ? CGAL::NEGATIVE : CGAL::POSITIVE;
        }

    } else {
        CGAL_unreachable();
    }

    return result;
}


// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


static Plane3SPtr bisector(FacetSPtr facet, VertexSPtr vertex);


Plane3SPtr SelfIntersection::bisector(FacetSPtr facet, VertexSPtr vertex) {
    Plane3SPtr result;
    EdgeSPtr edge_in;
    EdgeSPtr edge_out;
    std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
    while (it_e != vertex->edges().end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge(edge_wptr);
            if (edge->src(facet) == vertex) {
                edge_out = edge;
            } else if (edge->dst(facet) == vertex) {
                edge_in = edge;
            }
        }
    }
    if (edge_in && edge_out) {
        Point3SPtr point = vertex->getPoint();
        Point3SPtr p_in = edge_in->src(facet)->getPoint();
        Point3SPtr p_out = edge_out->dst(facet)->getPoint();
        Vector3SPtr normal = KernelFactory::createVector3(facet->plane());
        Point3SPtr p_normal = KernelFactory::createPoint3(*point + *normal);
        Plane3SPtr plane_in = KernelFactory::createPlane3(p_in, point, p_normal);
        Plane3SPtr plane_out = KernelFactory::createPlane3(point, p_out, p_normal);
        bool reflex = false;
        if (KernelWrapper::side(plane_in, p_out) > 0) {
            reflex = true;
        }

        // p_out is on the positive side of the bisector
        if (!reflex) {
            result = KernelWrapper::bisector(
                    KernelWrapper::opposite(plane_in), plane_out);
        } else {
            result = KernelWrapper::bisector(
                    plane_in, KernelWrapper::opposite(plane_out));
        }

        CGAL_postcondition(KernelWrapper::side(result, p_out) >= 0);
    }
    return result;
}




// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


this does not work because the closest point is usually a vertex, and then the wrong edge can be picked

static bool isInsideWithClosestEdge(const Point3& point,
                                    FacetSPtr facet,
                                    const bool handle_deg1_as_ray);

bool SelfIntersection::isInsideWithClosestEdge(const Point3& point,
                                               FacetSPtr facet,
                                               const bool handle_degree_1_as_ray)
{
    Plane3SPtr pl = facet->plane();
    Vector3SPtr normal = KernelFactory::createVector3(pl);

    CGAL::FT sq_dist_to_closest = std::numeric_limits<double>::max();
    EdgeSPtr closest_edge = nullptr;

    for (EdgeSPtr edge : facet->edges()) {
        // check collinearity
        if (handle_degree_1_as_ray) {
            CGAL_assertion(edge->getVertexSrc() != edge->getVertexDst());
            CGAL_precondition(edge->getVertexSrc()->degree() != 1 || edge->getVertexDst()->degree() != 1);

            VertexSPtr r_src = nullptr;
            VertexSPtr r_dst = nullptr;
            if (edge->getVertexSrc()->degree() == 1) {
                r_src = edge->getVertexDst();
                r_dst = edge->getVertexSrc();
            } else if (edge->getVertexDst()->degree() == 1) {
                r_src = edge->getVertexSrc();
                r_dst = edge->getVertexDst();
            }

            if (r_src && r_dst) {
                Ray3 r { *r_src->getPoint(), *r_dst->getPoint() };
                if (r.has_on(point)) {
                    // intersection if it's on the ray except if its the source
                    return (point != *(r_src->getPoint()));
                }
            } else {
                Segment3 s { *edge->getVertexSrc()->getPoint(),
                             *edge->getVertexDst()->getPoint() };
                if (s.has_on(point)) {
                    // intersection if it's on the ray except if its the source or target
                    return (point != *(edge->getVertexSrc()->getPoint()) &&
                            point != *(edge->getVertexDst()->getPoint()));
                }
            }
        } else {
            Segment3 s { *edge->getVertexSrc()->getPoint(),
                          *edge->getVertexDst()->getPoint() };
            if (s.has_on(point)) {
                // intersection if it's on the ray except if its the source or target
                return (point != *(edge->getVertexSrc()->getPoint()) &&
                        point != *(edge->getVertexDst()->getPoint()));
            }
        }

        // now we know the point is not exactly on the edge

        // We'll use a (coplanar) orientation check with the edge to determine
        // whether we are on the 'outside' or 'inside'
        VertexSPtr v_src = edge->src(facet);
        VertexSPtr v_dst = edge->dst(facet);
        Point3SPtr p_src = v_src->getPoint();
        Point3SPtr p_dst = v_dst->getPoint();

        // this only stands for EPECK and flat faces
        CGAL_assertion(pl->has_on(point));
        CGAL_assertion(pl->has_on(*p_src));
        CGAL_assertion(pl->has_on(*p_dst));
        CGAL_assertion(CGAL::scalar_product(Vector3(*p_src, *p_dst), *normal) == 0);

        // this means collinear but not on the geometry
        if (CGAL::collinear(*p_src, *p_dst, point)) {
            // std::cout << " skipping because collinear" << std::endl;
            continue;
        }

        auto treat_edge = [&](EdgeSPtr target_edge, const auto& edge_geometry) -> void
        {
            CGAL::FT sqd = CGAL::squared_distance(point, edge_geometry);
            // std::cout << "intersects & sq_dst: " << sqd << std::endl;
            CGAL::Comparison_result res = CGAL::compare(sqd, sq_dist_to_closest);
            if (res == CGAL::SMALLER) {
                sq_dist_to_closest = sqd;
                closest_edge = target_edge;
            }
        };

        if (handle_degree_1_as_ray) {
            if (v_src->degree() == 1) {
                if (v_dst->degree() == 1) {
                    // std::cout << "L3" << std::endl;
                    treat_edge(edge, Line3(*p_src, *p_dst));
                } else {
                    // std::cout << "R3 from DST " << *p_dst << std::endl;
                    treat_edge(edge, Ray3(*p_dst, *p_src));
                }
            } else if (v_dst->degree() == 1) {
                // std::cout << "R3 from SRC " << *p_src << std::endl;
                treat_edge(edge, Ray3(*p_src, *p_dst));
            } else {
                // std::cout << "S3" << std::endl;
                treat_edge(edge, Segment3(*p_src, *p_dst));
            }
        } else {
            treat_edge(edge, Segment3(*p_src, *p_dst));
        }
    }

    CGAL_postcondition (closest_edge != EdgeSPtr());

    // std::cout << "closest_edge = " << closest_edge->toString() << std::endl;
    Point3SPtr p_src = closest_edge->src(facet)->getPoint();
    Point3SPtr p_dst = closest_edge->dst(facet)->getPoint();

    DEBUG_PRINT("p_src = " << *p_src);
    DEBUG_PRINT("p_src + normal = " << *p_src + *normal);
    DEBUG_PRINT("p_dst = " << *p_dst);
    DEBUG_PRINT("point = " << point);
    CGAL_assertion(!CGAL::collinear(*p_src, *p_src + *normal, *p_dst));
    CGAL_assertion(CGAL::scalar_product(Vector3(*p_src, *p_src + *normal), Vector3(*p_src, *p_dst)) == 0);

    CGAL::Orientation o = CGAL::orientation(*p_src, *p_src + *normal, *p_dst, point);
    std::cout << "Orientation = " << o << std::endl;

    return (o != CGAL::NEGATIVE);
}




// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



CGAL::FT cx = 0, cy = 0, cz = 0;
std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
while (it_v != vertices_.end()) {
    VertexSPtr vertex = *it_v++;
    std::cout << "facet v " << *(vertex->getPoint()) << std::endl;
    cx += vertex->getPoint()->x();
    cy += vertex->getPoint()->y();
    cz += vertex->getPoint()->z();
}
cx /= CGAL::FT(vertices_.size());
cy /= CGAL::FT(vertices_.size());
cz /= CGAL::FT(vertices_.size());

// we want the plane such that at time 0 we are going through the centroid
CGAL::FT d = 0 /*speed*time */ - (cachedPlane_->a() * cx + cachedPlane_->b() * cy + cachedPlane_->c() * cz);



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


        // Step 3: Project original normal onto plane orthogonal to u
        CGAL::FT dot = a0*ux + b0*uy + c0*uz;
        CGAL::FT ab = a0 - dot * ux / uu;
        CGAL::FT bb = b0 - dot * uy / uu;
        CGAL::FT cb = c0 - dot * uz / uu;

#ifdef CGAL_SS3_CHECK_ALMOST_DEGENERATE_VECTORS_IN_FACET_PERTURBATIONS // these shouldn't be needed because in general the original normal is sane and very far from 'u'
        const double eps = 1e-12;

        // If the projection is (almost) zero, pick any orthogonal normal
        if (CGAL::abs(ab) < eps && CGAL::abs(bb) < eps && CGAL::abs(cb) < eps) {
            ab = 0;
            bb = uz;
            cb = -uy;
            if (CGAL::abs(bb) < eps && CGAL::abs(cb) < eps) {
                ab = -uz;
                bb = 0;
                cb = ux;
            }
        }
#endif

        // Step 4: Find a direction to nudge (cross product)
        CGAL::FT vx = uy * cb - uz * bb;
        CGAL::FT vy = uz * ab - ux * cb;
        CGAL::FT vz = ux * bb - uy * ab;

#ifdef CGAL_SS3_CHECK_ALMOST_DEGENERATE_VECTORS_IN_FACET_PERTURBATIONS // these shouldn't be needed because in general the original normal is sane and very far from 'u'
        // If v is zero, pick another orthogonal direction
        if (CGAL::abs(vx) < eps && CGAL::abs(vy) < eps && CGAL::abs(vz) < eps) {
            // Use cross with (1,0,0)
            vx = 0;
            vy = uz;
            vz = -uy;
            if (CGAL::abs(vy) < eps && CGAL::abs(vz) < eps) {
                vx = -uz;
                vy = 0;
                vz = ux;
            }
        }
#endif



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


static bool hasParallelPlanes(PolyhedronSPtr polyhedron);

bool PolyhedronTransformation::hasParallelPlanes(PolyhedronSPtr polyhedron) {
    bool result = false;
    std::list<FacetSPtr>::iterator it_f1 = polyhedron->facets().begin();
    while (it_f1 != polyhedron->facets().end()) {
        FacetSPtr facet1 = *it_f1++;
        std::list<FacetSPtr>::iterator it_f2 = it_f1;
        while (it_f2 != polyhedron->facets().end()) {
            FacetSPtr facet2 = *it_f2++;
            if (!KernelWrapper::intersection(
                    facet1->plane(), facet2->plane())) {
                result = true;
                break;
            }
        }
        if (result) {
            break;
        }
    }
    return result;
}

// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


/**
  * Normalize facet planes, ensuring parallel facets receive the same plane coefficients.
*/
static void harmonizeFacetPlanes(PolyhedronSPtr polyhedron);

// normalize coefficients and ensure that coplanar facets have the same coefficients
//
// @todo no point doing this with EPECK since we have SQRT errors
// and if the planes were singletons, we would benefit from static filters
void PolyhedronTransformation::harmonizeFacetPlanes(PolyhedronSPtr polyhedron)
{
    // @todo also harmonize parallel but non equal planes (i.e., opposite normal directions)

    // @todo this must be rewritten to use kernel predicates to sort the planes
    // Order all the supporting planes in a global order
    // The order is completely arbitrary, the only thing that we care about
    // is to quickly find which planes are identical
    //
    // for now, abuse unordered containers because it's less code to write than above
    struct FEqual
    {
        bool operator()(FacetSPtr facet_1, FacetSPtr facet_2) const {
            Plane3SPtr plane_1 = facet_1->plane(); // calls initPlane() if needed
            Plane3SPtr plane_2 = facet_2->plane();
            return CGAL::parallel(*plane_1, *plane_2);
        }
    };

    struct FHash
    {
        std::size_t operator()(FacetSPtr) const {
            return 1; // equality is only checked if the hash is the same
        }
    };

    std::unordered_map<FacetSPtr, Plane3SPtr, FHash, FEqual> facet_coefficients;

    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;

        auto res = facet_coefficients.emplace(facet, Plane3SPtr{});
        if(res.second) { // never seen this direction of planes before
            facet->normalizePlaneCoefficients();
            res.first->second = facet->getPlane();
        } else {
            // std::cout << "Facet #" << facet->getID() << " is reusing coefficients from Facet #" << res.first->first->getID() << std::endl;

            // a parallel plane already exists, so we re-use its a, b, c coordinates
            // but need the 'd' to be shifted so the plane goes through the facet's vertices
            Plane3SPtr plane = facet->plane(); // calls initPlane() if needed
            Plane3SPtr parallel_plane = res.first->second;

            Vector3SPtr n = KernelFactory::createVector3(plane);
            Vector3SPtr n_p = KernelFactory::createVector3(parallel_plane);
            CGAL::FT sign = (CGAL::scalar_product(*n, *n_p) > 0) ? 1 : -1;

            VertexSPtr v = facet->vertices().front();
            Point3SPtr p = v->getPoint();
            const CGAL::FT a = sign * parallel_plane->a();
            const CGAL::FT b = sign * parallel_plane->b();
            const CGAL::FT c = sign * parallel_plane->c();
            const CGAL::FT d = - (a * p->x() + b * p->y() + c * p->z());
            Plane3SPtr normalized_plane = KernelFactory::createPlane3(a, b, c, d);
            CGAL_assertion(normalized_plane->has_on(*p));
            facet->setPlane(normalized_plane);
        }

        const CGAL::FT a = facet->plane()->a();
        const CGAL::FT b = facet->plane()->b();
        const CGAL::FT c = facet->plane()->c();
        const CGAL::FT d = facet->plane()->d();
        CGAL_assertion_code(CGAL::FT sq_n = a*a + b*b + c*c;)
        // std::cout << "normalization check (post harmonize): " << sq_n
        //           << " (" << CGAL::to_interval(sq_n).first << "; "
        //           << CGAL::to_interval(sq_n).second << ")" << std::endl;
        // std::cout << "a: " << a << std::endl;
        // std::cout << "b: " << b << std::endl;
        // std::cout << "c: " << c << std::endl;
        // std::cout << "d: " << d << std::endl;

        // inaccuracies during normalization since the sqrt is (usually) not exact
        CGAL_assertion((sq_n - 1) < 1e-5);
    }
}



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //




std::pair<Point3SPtr, CGAL::FT>
SimpleStraightSkel::crashAt(EdgeSPtr edge_1, EdgeSPtr edge_2,
                            const CGAL::FT& offset_past_bound,
                            const CGAL::FT& offset_future_bound)
{
    FacetSPtr facet_l1 = edge_1->getFacetL();
    FacetSPtr facet_r1 = edge_1->getFacetR();
    FacetSPtr facet_l2 = edge_2->getFacetL();
    FacetSPtr facet_r2 = edge_2->getFacetR();
    CGAL::FT speed_l1 = std::dynamic_pointer_cast<SkelFacetData>(facet_l1->getData())->getSpeed();
    CGAL::FT speed_r1 = std::dynamic_pointer_cast<SkelFacetData>(facet_r1->getData())->getSpeed();
    CGAL::FT speed_l2 = std::dynamic_pointer_cast<SkelFacetData>(facet_l2->getData())->getSpeed();
    CGAL::FT speed_r2 = std::dynamic_pointer_cast<SkelFacetData>(facet_r2->getData())->getSpeed();

    CGAL_SS3_CORE_TRACE("-- Crash At\n    " << edge_1->toString() << "\n    " << edge_2->toString());

    CGAL_SS3_CORE_TRACE("Facet L1 = " << facet_l1->getID());
    CGAL_SS3_CORE_TRACE("Facet R1 = " << facet_r1->getID());
    CGAL_SS3_CORE_TRACE("Facet L2 = " << facet_l2->getID());
    CGAL_SS3_CORE_TRACE("Facet R2 = " << facet_r2->getID());

    Point3SPtr point;
    CGAL::FT offset_event;
    std::tie(point, offset_event) = intersectionPointAndTimeOffsetPlanes(facet_l1, facet_r1, facet_l2, facet_r2,
                                                                         offset_past_bound, offset_future_bound);

    if (!point) {
        return { };
    }

    // @speed we could delay point computation and orientation checks below to queue pop time,
    // but it probably does not gain much: 99.99% of the time, we compute an intersection time
    // and it is filtered, so the times where we compute a point and apply filters is somewhat
    // negligible (for now).

    CGAL_SS3_CORE_TRACE("Intersection: " << *point << " @ " << offset_event);

    // Check that the point is inside bounds
    FacetSPtr facet_1_src = edge_1->getFacetSrc();
    FacetSPtr facet_1_dst = edge_1->getFacetDst();
    FacetSPtr facet_2_src = edge_2->getFacetSrc();
    FacetSPtr facet_2_dst = edge_2->getFacetDst();

    CGAL_SS3_CORE_TRACE("Facet 1 SRC = " << facet_1_src->getID());
    CGAL_SS3_CORE_TRACE("Facet 1 DST = " << facet_1_dst->getID());
    CGAL_SS3_CORE_TRACE("Facet 2 SRC = " << facet_2_src->getID());
    CGAL_SS3_CORE_TRACE("Facet 2 DST = " << facet_2_dst->getID());

    CGAL_SS3_TRACE_CODE(CGAL::FT current_offset = (basePlanes_.at(facet_l1->getBasePlaneID())->d()
                                                   - facet_l1->getPlane()->d()) / speed_l1);
    CGAL_SS3_TRACE_CODE(CGAL::FT shift_offset = offset_event - current_offset);
    DEBUG_PRINT("current offset " << current_offset);
    DEBUG_PRINT("shift offset " << shift_offset);
    CGAL_SS3_TRACE_CODE(Segment3SPtr offset_e1 = PolyhedronTransformation::shiftEdge(edge_1, shift_offset);)
    CGAL_SS3_TRACE_CODE(Segment3SPtr offset_e2 = PolyhedronTransformation::shiftEdge(edge_2, shift_offset);)
    DEBUG_PRINT("Offset edge 1: " << *offset_e1);
    DEBUG_PRINT("Offset edge 2: " << *offset_e2);

#if defined(CGAL_SS3_OLD_CODE_BOUND_CHECKS) || defined(CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS)
    SkelEdgeDataSPtr data_1 = std::dynamic_pointer_cast<SkelEdgeData>(edge_1->getData());
    SkelEdgeDataSPtr data_2 = std::dynamic_pointer_cast<SkelEdgeData>(edge_2->getData());

    bool reject_1 = false;
    bool reject_2 = false;
    bool reject_3 = false;
    bool reject_4 = false;
    bool reject_5 = false;
    bool reject_6 = false;
    Vector3SPtr normal_1 = KernelFactory::createVector3(data_1->getSheet()->getPlane());
    Line3SPtr line_normal_1 = KernelFactory::createLine3(point, normal_1);
    if (KernelWrapper::orientation(line(edge_1), line_normal_1) <= 0) {
        DEBUG_PRINT("reject #1");
        reject_1 = true;
    }
    if (!(facet_1_src == facet_l2 ||
            facet_1_src == facet_r2)) {
        SkelVertexDataSPtr data_1_src = std::dynamic_pointer_cast<SkelVertexData>(
            edge_1->getVertexSrc()->getData());
        ArcSPtr arc_1_src = data_1_src->getArc();
        if (KernelWrapper::orientation(arc_1_src->line(), line_normal_1) > 0) {
            DEBUG_PRINT("reject #2");
            reject_2 = true;
        }
    }
    if (!(facet_1_dst == facet_l2 ||
            facet_1_dst == facet_r2)) {
        SkelVertexDataSPtr data_1_dst = std::dynamic_pointer_cast<SkelVertexData>(
            edge_1->getVertexDst()->getData());
        ArcSPtr arc_1_dst = data_1_dst->getArc();
        if (KernelWrapper::orientation(arc_1_dst->line(), line_normal_1) < 0) {
            DEBUG_PRINT("reject #3");
            reject_3 = true;
        }
    }
    Vector3SPtr normal_2 = KernelFactory::createVector3(data_2->getSheet()->getPlane());
    Line3SPtr line_normal_2 = KernelFactory::createLine3(point, normal_2);
    if (KernelWrapper::orientation(line(edge_2), line_normal_2) <= 0) {
            DEBUG_PRINT("reject #4");
            reject_4 = true;
    }
    if (!(facet_2_src == facet_l1 ||
            facet_2_src == facet_r1)) {
        SkelVertexDataSPtr data_2_src = std::dynamic_pointer_cast<SkelVertexData>(
            edge_2->getVertexSrc()->getData());
        ArcSPtr arc_2_src = data_2_src->getArc();
        if (KernelWrapper::orientation(arc_2_src->line(), line_normal_2) > 0) {
            DEBUG_PRINT("reject #5");
            reject_5 = true;
        }
    }
    if (!(facet_2_dst == facet_l1 ||
            facet_2_dst == facet_r1)) {
        SkelVertexDataSPtr data_2_dst = std::dynamic_pointer_cast<SkelVertexData>(
            edge_2->getVertexDst()->getData());
        ArcSPtr arc_2_dst = data_2_dst->getArc();
        if (KernelWrapper::orientation(arc_2_dst->line(), line_normal_2) < 0) {
            DEBUG_PRINT("reject #6");
            reject_6 = true;
        }
    }

    const bool reject_old = (reject_1 || reject_2 || reject_3 || reject_4 || reject_5 || reject_6);

#ifndef CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS
    if (reject_old) {
        return { };
    }
#endif

#endif

#if !defined(CGAL_SS3_OLD_CODE_BOUND_CHECKS) || defined(GAL_SS3_COMPARE_BOTH_BOUND_CHECKS)
    bool reject_2b = false;
    bool reject_3b = false;
    bool reject_5b = false;
    bool reject_6b = false;

    // Not to be in the past is equivalent to be on the negative side (we are shrinking)
    // of the planes of the facets incident to the edge
    //
    // Since we have filtered positive times, this shouldn't be needed (assuming positive weights)

    Plane3SPtr plane_l1 = facet_l1->getPlane();
    Plane3SPtr plane_r1 = facet_r1->getPlane();
    Plane3SPtr plane_l2 = facet_l2->getPlane();
    Plane3SPtr plane_r2 = facet_r2->getPlane();

    // this assumes positive weights
    CGAL_assertion(!(KernelWrapper::side(plane_l1, point) > 0 ||
                     KernelWrapper::side(plane_r1, point) > 0));
    CGAL_assertion(!(KernelWrapper::side(plane_l2, point) > 0 ||
                     KernelWrapper::side(plane_r2, point) > 0));

    // @todo this should be a predicate (oriented_side_of_event_point_wrt_bisectorC2)
    CGAL::FT l1a = plane_l1->a();
    CGAL::FT l1b = plane_l1->b();
    CGAL::FT l1c = plane_l1->c();
    CGAL::FT l1d = plane_l1->d();
    CGAL::FT r1a = plane_r1->a();
    CGAL::FT r1b = plane_r1->b();
    CGAL::FT r1c = plane_r1->c();
    CGAL::FT r1d = plane_r1->d();
    CGAL::FT l2a = plane_l2->a();
    CGAL::FT l2b = plane_l2->b();
    CGAL::FT l2c = plane_l2->c();
    CGAL::FT l2d = plane_l2->d();
    CGAL::FT r2a = plane_r2->a();
    CGAL::FT r2b = plane_r2->b();
    CGAL::FT r2c = plane_r2->c();
    CGAL::FT r2d = plane_r2->d();

    CGAL::FT lt1 = (l1a * point->x() + l1b * point->y() + l1c * point->z() + l1d) / speed_l1;
    CGAL::FT rt1 = (r1a * point->x() + r1b * point->y() + r1c * point->z() + r1d) / speed_r1;
    CGAL::FT lt2 = (l2a * point->x() + l2b * point->y() + l2c * point->z() + l2d) / speed_l2;
    CGAL::FT rt2 = (r2a * point->x() + r2b * point->y() + r2c * point->z() + r2d) / speed_r2;
    // DEBUG_PRINT("time from l1 " << lt1);
    // DEBUG_PRINT("time from r1 " << rt1);
    // DEBUG_PRINT("time from l2 " << lt2);
    // DEBUG_PRINT("time from r2 " << rt2);
    CGAL_assertion(lt1 == rt1 && lt1 == lt2 && lt1 == rt2);

    // We want the point to be left of the right arc, and right of the left arc
    //
    // We can just check the position of the query 'point' with respect
    // to one of the other bisector. A few tricky parts:
    // - The side of bisector is easy to know: it's comparing times, but that
    //   depends on whether the edge is reflex or convex
    // - Once we know on which side of the bisector it is, it's not done yet:
    //   we need to know which side of the bisector is the correct one and since
    //   faces can be non-convex polygons, we need to check if we are in a concave
    //   (within the facet) vertex to know if the clipping bisector is inverted

    // @todo lots of duplicate computations

    if (!(facet_1_src == facet_l2 ||
            facet_1_src == facet_r2)) {
        // DEBUG_PRINT("-- check 1_src");
        // src is the target of the edge when in the right facet
        if (!check_bisector(edge_1, facet_r1, rt1, facet_1_src, point)) {
            // DEBUG_PRINT("reject #2b");
#ifdef CGAL_SS3_EXIT_ASAP
            return { };
#else
            reject_2b = true;
#endif
        }
    }

    if (!(facet_1_dst == facet_l2 ||
            facet_1_dst == facet_r2)) {
        // DEBUG_PRINT("-- check 1_dst");
        if (!check_bisector(edge_1, facet_l1, lt1, facet_1_dst, point)) {
            // DEBUG_PRINT("reject #3b");
#ifdef CGAL_SS3_EXIT_ASAP
            return { };
#else
            reject_3b = true;
#endif
        }
    }

    if (!(facet_2_src == facet_l1 ||
            facet_2_src == facet_r1)) {
        // DEBUG_PRINT("-- check 2_src");
        if (!check_bisector(edge_2, facet_r2, rt2, facet_2_src, point)) {
            // DEBUG_PRINT("reject #5b");
#ifdef CGAL_SS3_EXIT_ASAP
            return { };
#else
            reject_5b = true;
#endif
        }
    }

    if (!(facet_2_dst == facet_l1 ||
            facet_2_dst == facet_r1)) {
          // DEBUG_PRINT("-- check 2_dst");
          if (!check_bisector(edge_2, facet_l2, lt2, facet_2_dst, point)) {
              // DEBUG_PRINT("reject #6b");
#ifdef CGAL_SS3_EXIT_ASAP
              return { };
#else
              reject_6b = true;
#endif
        }
    }

    const bool reject_new = (reject_2b || reject_3b || reject_5b || reject_6b);
# ifndef CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS
    if(reject_new) {
        return { };
    }
# endif
#endif

#ifdef CGAL_SS3_COMPARE_BOTH_BOUND_CHECKS
# ifdef CGAL_SS3_EXIT_ASAP
#  error // can't use ASAP-exits if we want to compare
# endif
    CGAL_assertion(reject_2 == reject_2b);
    CGAL_assertion(reject_3 == reject_3b);
    CGAL_assertion(reject_5 == reject_5b);
    CGAL_assertion(reject_6 == reject_6b);

    CGAL_assertion(reject_old == reject_new);

    if(reject_new) {
        return { };
    }
#endif

    return { point, offset_event };
}





// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //





bool SimpleStraightSkel::isReflex(EdgeSPtr edge,
                                  const bool future_facing) {
    bool result = false;

    VertexSPtr vertex_src = edge->getVertexSrc();
    VertexSPtr vertex_dst = edge->getVertexDst();

    if (*(vertex_src->getPoint()) == *(vertex_dst->getPoint())) {
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        FacetSPtr facet_src = edge->getFacetSrc();
        FacetSPtr facet_dst = edge->getFacetDst();

        CGAL::FT speed_l = 1.0;
        if (facet_l->hasData()) {
            speed_l = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_l->getData())->getSpeed();
        }
        CGAL::FT speed_r = 1.0;
        if (facet_r->hasData()) {
            speed_r = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_r->getData())->getSpeed();
        }
        CGAL::FT speed_src = 1.0;
        if (facet_src->hasData()) {
            speed_src = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_src->getData())->getSpeed();
        }
        CGAL::FT speed_dst = 1.0;
        if (facet_dst->hasData()) {
            speed_dst = std::dynamic_pointer_cast<SkelFacetData>(
                    facet_dst->getData())->getSpeed();
        }

        CGAL::FT od = future_facing ? -1 : 1;
        Plane3SPtr offset_plane_l = KernelWrapper::offsetPlane(facet_l->getPlane(), od * speed_l);
        Plane3SPtr offset_plane_r = KernelWrapper::offsetPlane(facet_r->getPlane(), od * speed_r);
        Plane3SPtr offset_plane_src = KernelWrapper::offsetPlane(facet_src->getPlane(), od * speed_src);
        Plane3SPtr offset_plane_dst = KernelWrapper::offsetPlane(facet_dst->getPlane(), od * speed_dst);

        Point3SPtr p_src = KernelWrapper::intersection(offset_plane_src, offset_plane_l, offset_plane_r);
        Point3SPtr p_dst = KernelWrapper::intersection(offset_plane_dst, offset_plane_l, offset_plane_r);
        CGAL_assertion(p_src && p_dst);

        Vector3SPtr v_dir = KernelFactory::createVector3((*p_dst) - (*p_src));
        CGAL_assertion(*v_dir != CGAL::NULL_VECTOR);

        Vector3SPtr normal_l = KernelFactory::createVector3(offset_plane_l);
        CGAL_assertion(*normal_l != CGAL::NULL_VECTOR);

        Vector3SPtr v_cross = KernelWrapper::cross(normal_l, v_dir);
        Point3SPtr p = KernelFactory::createPoint3((*p_src) + (*v_cross));
        if (KernelWrapper::side(offset_plane_r, p) > 0) {
            result = true;
        }
    } else {
        result = edge->isReflex();
    }
    return result;
}






// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


if ((src_y == CGAL::SMALLER && dst_y == CGAL::LARGER) ||
    (src_y == CGAL::LARGER && dst_y == CGAL::SMALLER)) {
    CGAL::Orientation o = orientation_2(*p_src, *p_dst, *point);
    if (o == CGAL::LEFT_TURN) is_inside = !is_inside;
    else if (o == CGAL::COLLINEAR) return CGAL::ON_BOUNDARY;
} else if (src_y == CGAL::EQUAL) {
    switch (compare_x_2(*point, *p_src)) {
        case CGAL::SMALLER:
            if (dst_y == CGAL::LARGER) is_inside = !is_inside;
            break;
        case CGAL::EQUAL:
            return CGAL::ON_BOUNDARY;
        case CGAL::LARGER:
            break;
    }
} else if (dst_y == CGAL::EQUAL) {
    switch (compare_x_2(*point, *p_dst)) {
        case CGAL::SMALLER:
            if (src_y == CGAL::LARGER) is_inside = !is_inside;
            break;
        case CGAL::EQUAL:
            return CGAL::ON_BOUNDARY;
        case CGAL::LARGER:
            break;
    }
}




// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



    VertexSPtr third_l = this->next(facet_l)->dst(facet_l);
    VertexSPtr third_r = this->next(facet_r)->dst(facet_r);
    Point3SPtr p_s = vertex_src_->getPoint();
    Point3SPtr p_d = vertex_dst_->getPoint();
    Point3SPtr p_l = third_l->getPoint();
    Point3SPtr p_r = third_r->getPoint();
    std::cout << "L NV " << facet_l->vertices().size() << std::endl;
    std::cout << "R NV " << facet_r->vertices().size() << std::endl;
    std::cout << "ps: " << *p_s << std::endl;
    std::cout << "pd: " << *p_d << std::endl;
    std::cout << "pl: " << *p_l << std::endl;
    std::cout << "pr: " << *p_r << std::endl;
    CGAL_assertion(!CGAL::collinear(*p_s, *p_d, *p_l));
    CGAL_assertion(!CGAL::collinear(*p_s, *p_d, *p_r));
    if (CGAL::is_positive(CGAL::scalar_product(*p_r - *p_l, *normal_l))) {
        result = true;
    }


// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

#ifndef CGAL_SS3_NO_CACHING
    std::unordered_map<std::array<int, 4>,
                       std::pair<Point3SPtr, CGAL::FT> > intersectionCache_;
#endif


// Sometimes should be disabled not to run out of memory.
// @todo Does not seem very effective anymore, need to bench again to check
// where the real cost is in events
// @todo could store only a Boolean value to improve memory footprint at the expense
// of computation.
// #define CGAL_SS3_NO_CACHING

#ifndef CGAL_SS3_NO_CACHING
            CGAL_SS3_CORE_TRACE_V(0, intersectionCache_.size() << " elements in crash cache");
#endif

    auto compute_point_and_time = [&]() -> std::pair<Point3SPtr, CGAL::FT>
    {
        Plane3SPtr plane_0 = facet_0->getBasePlane();
        Plane3SPtr plane_1 = facet_1->getBasePlane();
        Plane3SPtr plane_2 = facet_2->getBasePlane();
        Plane3SPtr plane_3 = facet_3->getBasePlane();

        CGAL::FT speed_0 = std::dynamic_pointer_cast<SkelFacetData>(facet_0->getData())->getSpeed();
        CGAL::FT speed_1 = std::dynamic_pointer_cast<SkelFacetData>(facet_1->getData())->getSpeed();
        CGAL::FT speed_2 = std::dynamic_pointer_cast<SkelFacetData>(facet_2->getData())->getSpeed();
        CGAL::FT speed_3 = std::dynamic_pointer_cast<SkelFacetData>(facet_3->getData())->getSpeed();

        return KernelWrapper::intersectionPointAndTimeOffsetPlanes(plane_0, speed_0, plane_1, speed_1,
                                                                   plane_2, speed_2, plane_3, speed_3,
                                                                   offset_past_bound, offset_future_bound);
    };

#ifdef CGAL_SS3_NO_CACHING
    return compute_point_and_time();
#else
    int ids[] = { facet_0->getID(),
                  facet_1->getID(),
                  facet_2->getID(),
                  facet_3->getID() };
    std::sort(std::begin(ids), std::end(ids));
    std::array<int, 4> canonical_ids = CGAL::make_array(ids[0], ids[1], ids[2], ids[3]);

    std::pair<Point3SPtr, CGAL::FT> dummy;
    auto res = intersectionCache_.emplace(canonical_ids, dummy);
    if (res.second) { // successful insertion, first time seeing it so actual computation is required
        CGAL_SS3_CORE_TRACE_V(0, "compute needed: " << canonical_ids[0] << " " << canonical_ids[1] << " "
                                                     << canonical_ids[2] << " " << canonical_ids[3]);

        res.first->second = compute_point_and_time();
    } else {
        CGAL_SS3_CORE_TRACE_V(0, "used cache value: " << canonical_ids[0] << " " << canonical_ids[1] << " "
                                                       << canonical_ids[2] << " " << canonical_ids[3]);
    }

    return res.first->second;
#endif // CGAL_SS3_NO_CACHING

// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //