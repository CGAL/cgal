/**
 * @file   data/3d/Polyhedron.h
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#ifndef DATA_3D_POLYHEDRON_H
#define DATA_3D_POLYHEDRON_H

#include "data/3d/ptrs.h"
#include "typedefs_thread.h"
#include <list>
#include <string>

namespace data { namespace _3d {

class Polyhedron : public std::enable_shared_from_this<Polyhedron> {
public:
    virtual ~Polyhedron();

    static PolyhedronSPtr create();
    static PolyhedronSPtr create(unsigned int num_facets, FacetSPtr facets[]);

    PolyhedronSPtr clone() const;

    void addVertex(VertexSPtr vertex);
    bool removeVertex(VertexSPtr vertex);

    /**
     * Searches for a vertex with the same coordinates as the given vertex.
     */
    VertexSPtr findVertex(VertexSPtr needle);

    void addEdge(EdgeSPtr edge);
    bool removeEdge(EdgeSPtr edge);

    /**
     * Searches for an edge with the same coordinates as the given edge.
     * The orientation of the edge is ignored.
     */
    EdgeSPtr findEdge(EdgeSPtr needle);

    void addFacet(FacetSPtr facet);
    bool removeFacet(FacetSPtr facet);

    void initPlanes();
    void clearData();

    SharedMutex& mutex();

    std::list<VertexSPtr>& vertices();
    std::list<EdgeSPtr>& edges();
    std::list<FacetSPtr>& facets();

    bool isConsistent() const;
    void clear();

    int getID() const;
    void setID(int id);

    void resetAllIDs();

    std::string getDescription() const;
    void setDescription(const std::string& desc);
    void appendDescription(const std::string& desc);

    std::string toString() const;

protected:
    Polyhedron();
    mutable SharedMutex mutex_;
    std::list<VertexSPtr> vertices_;
    std::list<EdgeSPtr> edges_;
    std::list<FacetSPtr> facets_;
    int id_;
    std::string description_;
};

} }

#endif /* DATA_3D_POLYHEDRON_H */
