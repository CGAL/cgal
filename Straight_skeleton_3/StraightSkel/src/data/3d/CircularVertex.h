/**
 * @file   data/3d/CircularVertex.h
 * @author Gernot Walzl
 * @date   2012-11-28
 */

#ifndef DATA_3D_CIRCULARVERTEX_H
#define DATA_3D_CIRCULARVERTEX_H

#include "data/3d/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _3d {

class CircularVertex {
public:
    virtual ~CircularVertex();

    static CircularVertexSPtr create(Point3SPtr point);

    CircularVertexSPtr clone() const;

    Point3SPtr getPoint() const;
    void setPoint(Point3SPtr point);
    bool isPointValid() const;
    void setPointValid(bool point_valid);
    CircularEdgeSPtr getEdgeIn() const;
    void setEdgeIn(CircularEdgeSPtr edge);
    CircularEdgeSPtr getEdgeOut() const;
    void setEdgeOut(CircularEdgeSPtr edge);
    SphericalPolygonSPtr getPolygon() const;
    void setPolygon(SphericalPolygonSPtr polygon);
    std::list<CircularVertexSPtr>::iterator getListIt() const;
    void setListIt(std::list<CircularVertexSPtr>::iterator list_it);
    CircularVertexDataSPtr getData() const;
    void setData(CircularVertexDataSPtr data);
    bool hasData() const;

    CircularVertexSPtr next() const;
    CircularVertexSPtr prev() const;

    double getX() const;
    double getY() const;
    double getZ() const;

    std::string toString() const;

protected:
    CircularVertex(Point3SPtr point);
    CircularVertex(const CircularVertex& vertex);
    Point3SPtr point_;
    bool point_valid_;
    CircularEdgeWPtr edge_in_;
    CircularEdgeWPtr edge_out_;
    SphericalPolygonWPtr polygon_;
    std::list<CircularVertexSPtr>::iterator list_it_;
    CircularVertexDataSPtr data_;
};

} }

#endif /* DATA_3D_CIRCULARVERTEX_H */

