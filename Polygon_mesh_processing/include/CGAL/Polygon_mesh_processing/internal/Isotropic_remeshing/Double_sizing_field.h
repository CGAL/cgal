//
// Created by kabir on 30/03/21.
//

#ifndef CGAL_DOUBLE_SIZING_FIELD_H
#define CGAL_DOUBLE_SIZING_FIELD_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>
#include "Sizing_field.h"
#include <CGAL/number_utils.h>
#include "Sizing_field.h"

namespace CGAL

    template<class PolygonMesh>
    {
    class Double_sizing_field : public CGAL::Sizing_feild<PolygonMesh> {
    private:
        typedef CGAL::Sizing_field<PolygonMesh> Base;
    public:
        typedef typename Base::FT         FT;
        typedef typename Base::Point_3    Point_3;
        typedef typename Base::halfedge_descriptor halfedge_descriptor;
        typedef typename Base::vertex_descriptor   vertex_descriptor;
    private:
        const PolygonMesh& polygonMesh;
        FT m_sq_short;
        FT m_sq_long;
    public:
        FT sizing_field;

        Double_sizing_field(const Base::FT &size, const PolygonMesh &mesh) :
                polygonMesh(mesh),
                sizing_field(2*size),
                m_sq_short(CGAL::square(8./5.*size)),
                m_sq_long(CGAL::square(8./3.*size))
                {}

    private:
        FT sqlength(const vertex_descriptor& va,
                    const vertex_descriptor& vb) const
        {
            typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type
                    vpmap = get(CGAL::vertex_point, m_pmesh);
            return FT(CGAL::squared_distance(get(vpmap, va), get(vpmap, vb)));
        }

        FT sqlength(const halfedge_descriptor& h) const
        {
            return sqlength(target(h, m_pmesh), source(h, m_pmesh));
        }

    public:
        void print() {
            return sizing_field;
        }
        boost::optional<FT> is_too_long(const halfedge_descriptor& h) const
        {
            const FT sqlen = sqlength(h);
            if(sqlen > m_sq_long)
                return sqlen;
            else
                return boost::none;
        }

        boost::optional<FT> is_too_long(const vertex_descriptor& va,
                                        const vertex_descriptor& vb) const
        {
            const FT sqlen = sqlength(va, vb);
            if (sqlen > m_sq_long)
                return sqlen;
            else
                return boost::none;
        }

        boost::optional<FT> is_too_short(const halfedge_descriptor& h) const
        {
            const FT sqlen = sqlength(h);
            if (sqlen < m_sq_long)
                return sqlen;
            else
                return boost::none;
        }

        virtual Point_3 split_placement(const halfedge_descriptor& h) const
        {
            typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type
                    vpmap = get(CGAL::vertex_point, m_pmesh);
            return CGAL::midpoint(get(vpmap, target(h, m_pmesh)),
                                  get(vpmap, source(h, m_pmesh)));
        }
    };
}

}// ending namespace CGAL


#endif //CGAL_DOUBLE_SIZING_FIELD_H
