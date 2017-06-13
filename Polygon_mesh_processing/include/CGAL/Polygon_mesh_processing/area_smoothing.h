#ifndef AREA_SMOOTHING_H
#define AREA_SMOOTHING_H

#include <fstream>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>


#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
//#include <CGAL/AABB_face_graph_triangle_primitive.h>





namespace CGAL {

namespace Polygon_mesh_processing {








template<typename PolygonMesh, typename VertexPointMap, typename GeomTraits>
class Area_remesher
{

    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

    typedef typename GeomTraits::Point_3 Point;
    typedef typename GeomTraits::Vector_3 Vector;
    typedef typename GeomTraits::Triangle_3 Triangle;

    typedef std::vector<Triangle> Triangle_list;


    typedef CGAL::AABB_triangle_primitive<GeomTraits, typename Triangle_list::iterator> Primitive;
    typedef CGAL::AABB_traits<GeomTraits, Primitive> AABB_Traits;
    typedef CGAL::AABB_tree<AABB_Traits> Tree;






public:
    Area_remesher(PolygonMesh& pmesh, VertexPointMap& vpmap) : mesh_(pmesh), vpmap_(vpmap)
    {}

    ~Area_remesher()
    {
        delete tree_ptr_;
    }


    template<typename FaceRange>
    void init_remeshing(const FaceRange& face_range)
    {
        for(face_descriptor f : face_range)
            input_triangles_.push_back(triangle(f));


        tree_ptr_ = new Tree(input_triangles_.begin(), input_triangles_.end());
        tree_ptr_->accelerate_distance_queries();

    }



    void area_relaxation()
    {

        for(vertex_descriptor v : vertices(mesh_))
        {

             if(!is_border(v, mesh_))
             {
                 std::cout<<"processing vertex: "<< v << std::endl;

                 std::cout<<"point moved from: "<<get(vpmap_, v);

                 gradient_descent(v);

                 Point p_query = get(vpmap_, v);
                 Point projected = tree_ptr_->closest_point(p_query);

                 //do the projection
                 put(vpmap_, v, projected);

                 std::cout<<" after projection: "<<get(vpmap_, v)<<std::endl;


             } // not on border


        }


    }







private:


    Triangle triangle(face_descriptor f) const
    {
        halfedge_descriptor h = halfedge(f, mesh_);
        vertex_descriptor v1 = target(h, mesh_);
        vertex_descriptor v2 = target(next(h, mesh_), mesh_);
        vertex_descriptor v3 = target(next(next(h, mesh_), mesh_), mesh_);
        return Triangle(get(vpmap_, v1), get(vpmap_, v2), get(vpmap_, v3));
    }



    double element_area(const vertex_descriptor& p1,
                        const vertex_descriptor& p2,
                        const vertex_descriptor& p3) const
    {
        return to_double(CGAL::approximate_sqrt(
                             GeomTraits().compute_squared_area_3_object()(
                                   get(vpmap_, p1),
                                   get(vpmap_, p2),
                                   get(vpmap_, p3))));

    }

    double compute_average_area_around(const vertex_descriptor& v)
    {

        double sum_areas = 0;
        unsigned int number_of_edges = 0;

        for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
        {
            // opposite vertices
            vertex_descriptor pi = source(next(h, mesh_), mesh_);
            vertex_descriptor pi1 = target(next(h, mesh_), mesh_);
            //std::cout<<"("<<v<<","<<pi<<","<<pi1<<")";

            double S = element_area(v, pi, pi1);
            //std::cout<<"- Area= "<<S<<std::endl;
            sum_areas += S;
            number_of_edges++;
        }

        return sum_areas / number_of_edges;

    }


    double measure_energy(const vertex_descriptor& v, const double& S_av)
    {
        double energy = 0;
        unsigned int number_of_edges = 0;


        for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
        {
            vertex_descriptor pi = source(next(h, mesh_), mesh_);
            vertex_descriptor pi1 = target(next(h, mesh_), mesh_);
            double S = element_area(v, pi, pi1);

            energy += (S - S_av)*(S - S_av);
            number_of_edges++;

        }

        return (double) (energy / (2* number_of_edges) );


    }




    void compute_derivatives(double& dFdx, double& dFdy, double& dFdz, const vertex_descriptor& v, const double& S_av)
    {

        for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
        {

            vertex_descriptor pi = source(next(h, mesh_), mesh_);
            vertex_descriptor pi1 = target(next(h, mesh_), mesh_);
            double S = element_area(v, pi, pi1);
            //std::cout<<"h="<<h<<" - A= "<<S<<std::endl;


            Vector vec(get(vpmap_, pi), get(vpmap_, pi1));


            //dFdx += (S - S_av) * (yi - yi1);
            //dFdx += (S - S_av) * (-vec.y());

            //dFdy += (S - S_av) * (xi1 - xi);
            //dFdy += (S - S_av) * (vec.x());


            dFdx += (S - S_av) * (vec.z() - vec.y());
            dFdy += (S - S_av) * (vec.x() - vec.z());
            dFdz += (S - S_av) * (vec.y() - vec.x());


        }
        //std::cout<<std::endl;

        dFdx *= 2;
        dFdy *= 2;
        dFdz *= 2;

    }




    //void gradient_descent(vertex_descriptor& v, bool do_project)
    void gradient_descent(vertex_descriptor& v)
    {

        double eta = 0.01; //learning rate
        double precision = 0.001;
        double energy;
        double x, y, z, x_new, y_new, z_new, dFdx, dFdy, dFdz;
        double previous_step_size_x = precision + 1; // temp; to be improved
        double previous_step_size_y = precision + 1;
        double previous_step_size_z = precision + 1;


        double S_av = compute_average_area_around(v);

        x = get(vpmap_, v).x();
        y = get(vpmap_, v).y();
        z = get(vpmap_, v).z();

        std::ofstream out("data/energy.txt");


        while(previous_step_size_x > precision || previous_step_size_y > precision || previous_step_size_z > precision)
        {

            dFdx=0, dFdy=0, dFdz=0;
            compute_derivatives(dFdx, dFdy, dFdz, v, S_av);

            x_new = x - eta * dFdx;
            y_new = y - eta * dFdy;
            z_new = z - eta * dFdz;

            Point moved(x_new, y_new, z_new);
            put(vpmap_, v, moved);

            previous_step_size_x = CGAL::abs(x_new - x);
            previous_step_size_y = CGAL::abs(y_new - y);
            previous_step_size_z = CGAL::abs(z_new - z);
            //std::cout<<"previous_step_size_x: "<<previous_step_size_x<<std::endl;

            x = x_new;
            y = y_new;
            z = z_new;

            energy = measure_energy(v, S_av);
            out << energy << std::endl;


        }



    } //gradient descent





private:
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;
    const Tree* tree_ptr_;
    Triangle_list input_triangles_;









};














template<typename PolygonMesh, typename NamedParameters, typename FaceRange>
void area_remeshing(PolygonMesh& pmesh, const NamedParameters& np, const FaceRange& faces)
{


    using boost::choose_param;
    using boost::get_param;

    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VertexPointMap;
    VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                 get_const_property_map(CGAL::vertex_point, pmesh));


    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GeomTraits;


    //bool do_project = choose_param(get_param(np, internal_np::do_project), true);


    CGAL::Polygon_mesh_processing::Area_remesher<PolygonMesh, VertexPointMap, GeomTraits> remesher(pmesh, vpmap);
    remesher.init_remeshing(faces);
    remesher.area_relaxation();






}










} // namespace Polygon_mesh_processing
} //namespace CGAL








#endif // AREA_SMOOTHING_H


