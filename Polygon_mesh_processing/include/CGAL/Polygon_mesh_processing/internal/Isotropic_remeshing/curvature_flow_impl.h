#ifndef CURVATURE_FLOW_IMPL_H
#define CURVATURE_FLOW_IMPL_H

#include <CGAL/Polygon_mesh_processing/Weights.h>
#include <CGAL/Monge_via_jet_fitting.h>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {



template<typename PolygonMesh, typename VertexPointMap, typename GeomTraits>
class Curvature_flow
{

    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;


    typedef typename CGAL::Monge_via_jet_fitting<GeomTraits> Monge_via_jet_fitting;
    typedef typename Monge_via_jet_fitting::Monge_form Monge_form;

    typedef typename GeomTraits::Point_3  Point;
    typedef typename GeomTraits::Vector_3 Vector;
    typedef typename GeomTraits::FT       FT;



public:


    Curvature_flow(PolygonMesh& pmesh, VertexPointMap& vpmap) : mesh_(pmesh), vpmap_(vpmap),
                                                                cot_calculator_(mesh_, vpmap_)
    {

        //std::vector<vertex_descriptor> points
        std::size_t num_points = vertices(mesh_).size();
        std::size_t min_num_of_points = 6; // (d+1)(d+2)/2, for d=2
        if(num_points < min_num_of_points)
        {
            CGAL_error_msg("Not enough points in the mesh.");
        }

    }




    void curvature_smoothing()
    {
        for(vertex_descriptor v : vertices(mesh_))
        {

            if(!is_border(v, mesh_)) // add && !is_constrained
            {


                double k_mean = compute_mean_curvature();
                std::cout<<"ok"<<std::endl;




            }


        }




    }






private:


    double compute_mean_curvature()
    {
        // gather points around v
        std::vector<Point> incident_points = gather_all_points();

        Monge_form monge_form;
        Monge_via_jet_fitting monge_fit;

        std::size_t d_fit = 2; // d_fit >= d_monge
        std::size_t d_monge = 2; // need 2 principal coeeficients
        std::size_t Nd = (d_fit + 1)*(d_fit + 1) / 2.0;
        CGAL_assertion(incident_points.size() >= Nd);

        monge_form = monge_fit(incident_points.begin(), incident_points.end(), d_fit, d_monge);
        const FT k1 = monge_form.principal_curvatures(0);
        const FT k2 = monge_form.principal_curvatures(1);

        return (k1 + k2) / 2.0;
    }


    std::vector<Point> points_around_vertex(vertex_descriptor v)
    {
        std::vector<Point> incident_vertices;
        for(halfedge_descriptor h : halfedges_around_target(v, mesh_))
        {
            vertex_descriptor vs = source(h, mesh_);
            incident_vertices.push_back(get(vpmap_, vs));
        }

        // temp assertion
        std::vector<Point> incident_vertices2;
        for(vertex_descriptor vi : vertices_around_target(v, mesh_))
        {
           incident_vertices2.push_back(get(vpmap_, vi));
        }
        CGAL_assertion(incident_vertices.size() == incident_vertices2.size());

        return incident_vertices;
    }


    std::vector<Point> gather_all_points()
    {
        std::vector<Point> points; // todo: preallocate it and fill it with pointing
        for(vertex_descriptor v : vertices(mesh_))
        {
            points.push_back(get(vpmap_, v));
        }

        return points;
    }



    // data members
    // ------------
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;
    // Cotagent calculator class
    CGAL::internal::Cotangent_value_Meyer<
        PolygonMesh,
        typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type >   cot_calculator_;

};







} // internal
} // Polygon_mesh_processing
} // CGAL





#endif // CURVATURE_FLOW_IMPL_H
