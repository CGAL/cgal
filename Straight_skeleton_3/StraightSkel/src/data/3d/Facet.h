// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * @file   data/3d/Facet.h
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#ifndef DATA_3D_FACET_H
#define DATA_3D_FACET_H

#include "data/2d/ptrs.h"
#include "data/3d/ptrs.h"

#include <list>
#include <string>
#include <vector>

namespace data { namespace _3d {

class Facet : public std::enable_shared_from_this<Facet> {
public:
    virtual ~Facet();

    Facet();

    static FacetSPtr create();

    static FacetSPtr create(const std::vector<VertexSPtr>& vertices);

    /**
     * Each edge may have a facet to the left and a facet to the right.
     * Left and right is defined by a view outside the polyhedron.
     * In case an edge is used the first time, it is assumed that the left
     * side of the edge points to the interior of the facet.
     * In case the edge is already part of a facet,
     * it is assumed that the right side of the edge points to
     * the interior of the facet that is created.
     */
    static FacetSPtr create(const std::vector<EdgeSPtr>& edges);

    FacetSPtr clone() const;

    void addVertex(VertexSPtr vertex);
    bool removeVertex(VertexSPtr vertex);
    bool hasVertex(VertexSPtr vertex);

    bool isTriangle() const;

    void addEdge(EdgeSPtr edge);
    bool removeEdge(EdgeSPtr edge);

    /**
     * Searches for an edge to the given facet.
     * Does not work when there is more than one edge adjacent to both facets.
     * @deprecated
     */
    EdgeSPtr findEdge(FacetSPtr facet) const;

    std::list<EdgeSPtr> findEdges(FacetSPtr facet) const;

    bool containsVertex(VertexSPtr vertex) const;
    bool containsEdge(EdgeSPtr edge) const;

    void sortVertices();
    void sortEdges();

    PolyhedronSPtr getPolyhedron() const;
    void setPolyhedron(PolyhedronSPtr polyhedron);
    std::list<FacetSPtr>::iterator getPolyhedronListIt() const;
    void setPolyhedronListIt(std::list<FacetSPtr>::iterator list_it);

    FacetDataSPtr getData() const;
    void setData(FacetDataSPtr data);
    bool hasData() const;

    std::list<VertexSPtr>& vertices();
    std::list<EdgeSPtr>& edges();

    /**
     * counter clockwise from outside
     */
    FacetSPtr next(VertexSPtr vertex) const;
    FacetSPtr prev(VertexSPtr vertex) const;

    /**
     * merge 'facet' into this facet.
     */
    void merge(FacetSPtr facet);

    int getID() const;
    void setID(int id);

    /**
     * The direction of the normal points to the outside.
     */
    Plane3SPtr getPlane() const;
    void setPlane(Plane3SPtr plane);

    Plane3SPtr getBasePlane() const;
    void setBasePlane(Plane3SPtr plane);

    bool hasFinalPlane() const;
    Plane3SPtr getFinalPlane() const;
    void setFinalPlane(Plane3SPtr plane);

    /**
     * First vertices have to form a triangle that is inside.
     */
    bool initPlane();
    Plane3SPtr plane();

    /**
     * Normalize the plane coefficients to obtain a canonical plane representation
     */
    void normalizePlaneCoefficients();

    /**
     * Check if the plane is normalized
     */
    bool isNormalizedPlane();

    /**
     * Store the current plane ahead of perturbation
     */
    void storePlaneCoefficients();

    enum class PerturbationType {
        NUDGE,
        STEPS,
        EXACT,
        HIGH_DEGREES
    };

    /**
     * Nudge the plane coefficients by a random value in the range [low, high].
     */
    void perturbPlaneCoefficientsNudge(const double range);

    /**
     * Nudge the plane coefficients by a number of nextafter steps
     * in a random direction.
     */
    void perturbPlaneCoefficientsSteps(int steps);

    /**
     * Nudge the plane coefficients by a random value in the range [0, 1] / den.
     */
    void perturbPlaneCoefficientsExact(const CGAL::FT& den);

    /**
     * Nudge the plane coefficients but ensure that the perturbed plane goes through 0, 1, or 2 fixed points.
     * If 0 points: nudge all coefficients independently.
     * If 1 point: nudge (a, b, c), recompute d so the plane passes through the point.
     * If 2 points: nudge (a, b, c) with the constraint that the new plane passes through both points.
     */
    void perturbPlaneCoefficientsFixedPoints(const double range,
                                             const std::vector<Point3SPtr>& fixed_points);

    void perturbPlaneCoefficientsHighDegrees(const double range);

    /**
     * Nudge the plane coefficients to get rid of simultaneous events
     */
    void perturbPlaneCoefficients(PerturbationType type = PerturbationType::HIGH_DEGREES);

    /**
     * Restore the plane coefficients to the previous value, updating 'd' so that the plane
     * matches the desired offset.
     */
    void restorePlaneCoefficients(const CGAL::FT& perturbationOffset,
                                  const CGAL::FT& perturbationEndOffset);

    bool makeFirstConvex();

    /**
     * returns the number of vertices whose degree is higher than 3
     */
    int numHighDegreeVertices() const;

    data::_2d::PolygonSPtr toPolygon();

    std::string toString() const;

public:
    std::list<VertexSPtr> vertices_;
    std::list<EdgeSPtr> edges_;
    PolyhedronWPtr polyhedron_;
    std::list<FacetSPtr>::iterator polyhedron_list_it_;
    FacetDataSPtr data_;

    Plane3SPtr plane_;
    Plane3SPtr base_plane_;
    Plane3SPtr final_plane_;

    int id_;

    Plane3SPtr cachedPlane_;
    CGAL::FT cachedSpeed_;
};

} }

#endif /* DATA_3D_FACET_H */
