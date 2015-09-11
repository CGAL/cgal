#include "Scene_polyhedron_item.h"
#include <CGAL/AABB_intersections.h>
#include "Kernel_type.h"
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>



#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_2_filtered_projection_traits_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <QVariant>
#include <list>
#include <queue>
#include <iostream>
#include <QDebug>


#include <limits>

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_traits;
typedef CGAL::AABB_tree<AABB_traits> Input_facets_AABB_tree;
const char* aabb_property_name = "Scene_polyhedron_item aabb tree";

Input_facets_AABB_tree* get_aabb_tree(Scene_polyhedron_item* item)
{
    QVariant aabb_tree_property = item->property(aabb_property_name);
    if(aabb_tree_property.isValid()) {
        void* ptr = aabb_tree_property.value<void*>();
        return static_cast<Input_facets_AABB_tree*>(ptr);
    }
    else {
        Polyhedron* poly = item->polyhedron();
        if(poly) {
            Input_facets_AABB_tree* tree =
                    new Input_facets_AABB_tree(faces(*poly).first,
                                               faces(*poly).second,
                                               *poly);
            item->setProperty(aabb_property_name,
                              QVariant::fromValue<void*>(tree));
            return tree;
        }
        else return 0;
    }
}

void delete_aabb_tree(Scene_polyhedron_item* item)
{
    QVariant aabb_tree_property = item->property(aabb_property_name);
    if(aabb_tree_property.isValid()) {
        void* ptr = aabb_tree_property.value<void*>();
        Input_facets_AABB_tree* tree = static_cast<Input_facets_AABB_tree*>(ptr);
        if(tree) {
            delete tree;
            tree = 0;
        }
        item->setProperty(aabb_property_name, QVariant());
    }
}

typedef Polyhedron::Traits Traits;
typedef Polyhedron::Facet Facet;
typedef CGAL::Triangulation_2_filtered_projection_traits_3<Traits>   P_traits;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
struct Face_info {
    Polyhedron::Halfedge_handle e[3];
    bool is_external;
};
typedef CGAL::Triangulation_vertex_base_with_info_2<Halfedge_handle,
P_traits>        Vb;
typedef CGAL::Triangulation_face_base_with_info_2<Face_info,
P_traits>          Fb1;
typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>   Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                  TDS;
typedef CGAL::No_intersection_tag                                    Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits,
TDS,
Itag>             CDTbase;
typedef CGAL::Constrained_triangulation_plus_2<CDTbase>              CDT;

//Make sure all the facets are triangles
void
Scene_polyhedron_item::is_Triangulated() const
{
    typedef Polyhedron::Halfedge_around_facet_circulator HF_circulator;
    Facet_iterator f = poly->facets_begin();
    int nb_points_per_facet =0;

    for(f = poly->facets_begin();
        f != poly->facets_end();
        f++)
    {
        HF_circulator he = f->facet_begin();
        HF_circulator end = he;
        CGAL_For_all(he,end)
        {
            nb_points_per_facet++;
        }

        if(nb_points_per_facet !=3)
        {
            is_Triangle = false;
            break;
        }

        nb_points_per_facet = 0;
    }
}
void
Scene_polyhedron_item::triangulate_facet(Facet_iterator fit) const
{
    //Computes the normal of the facet
    Traits::Vector_3 normal =
            CGAL::Polygon_mesh_processing::compute_face_normal(fit,*poly);

    P_traits cdt_traits(normal);
    CDT cdt(cdt_traits);

    Facet::Halfedge_around_facet_circulator
            he_circ = fit->facet_begin(),
            he_circ_end(he_circ);

    // Iterates on the vector of facet handles
    CDT::Vertex_handle previous, first;
    do {
        CDT::Vertex_handle vh = cdt.insert(he_circ->vertex()->point());
        if(first == 0) {
            first = vh;
        }
        vh->info() = he_circ;
        if(previous != 0 && previous != vh) {
            cdt.insert_constraint(previous, vh);
        }
        previous = vh;
    } while( ++he_circ != he_circ_end );
    cdt.insert_constraint(previous, first);

    // sets mark is_external
    for(CDT::All_faces_iterator
        fit2 = cdt.all_faces_begin(),
        end = cdt.all_faces_end();
        fit2 != end; ++fit2)
    {
        fit2->info().is_external = false;
    }
    //check if the facet is external or internal
    std::queue<CDT::Face_handle> face_queue;
    face_queue.push(cdt.infinite_vertex()->face());
    while(! face_queue.empty() ) {
        CDT::Face_handle fh = face_queue.front();
        face_queue.pop();
        if(fh->info().is_external) continue;
        fh->info().is_external = true;
        for(int i = 0; i <3; ++i) {
            if(!cdt.is_constrained(std::make_pair(fh, i)))
            {
                face_queue.push(fh->neighbor(i));
            }
        }
    }

    //iterates on the internal faces to add the vertices to the positions
    //and the normals to the appropriate vectors
    for(CDT::Finite_faces_iterator
        ffit = cdt.finite_faces_begin(),
        end = cdt.finite_faces_end();
        ffit != end; ++ffit)
    {
        if(ffit->info().is_external)
            continue;

        positions_facets.push_back(ffit->vertex(0)->point().x());
        positions_facets.push_back(ffit->vertex(0)->point().y());
        positions_facets.push_back(ffit->vertex(0)->point().z());
        positions_facets.push_back(1.0);

        positions_facets.push_back(ffit->vertex(1)->point().x());
        positions_facets.push_back(ffit->vertex(1)->point().y());
        positions_facets.push_back(ffit->vertex(1)->point().z());
        positions_facets.push_back(1.0);

        positions_facets.push_back(ffit->vertex(2)->point().x());
        positions_facets.push_back(ffit->vertex(2)->point().y());
        positions_facets.push_back(ffit->vertex(2)->point().z());
        positions_facets.push_back(1.0);


        typedef Kernel::Vector_3	    Vector;
        Vector n = CGAL::Polygon_mesh_processing::compute_face_normal(fit, *poly);
        normals_flat.push_back(n.x());
        normals_flat.push_back(n.y());
        normals_flat.push_back(n.z());

        normals_flat.push_back(n.x());
        normals_flat.push_back(n.y());
        normals_flat.push_back(n.z());

        normals_flat.push_back(n.x());
        normals_flat.push_back(n.y());
        normals_flat.push_back(n.z());

        normals_gouraud.push_back(n.x());
        normals_gouraud.push_back(n.y());
        normals_gouraud.push_back(n.z());

        normals_gouraud.push_back(n.x());
        normals_gouraud.push_back(n.y());
        normals_gouraud.push_back(n.z());

        normals_gouraud.push_back(n.x());
        normals_gouraud.push_back(n.y());
        normals_gouraud.push_back(n.z());

    }
}

void
Scene_polyhedron_item::triangulate_facet_color(Facet_iterator fit) const
{
    Traits::Vector_3 normal =
            CGAL::Polygon_mesh_processing::compute_face_normal(fit, *poly);

    P_traits cdt_traits(normal);
    CDT cdt(cdt_traits);

    Facet::Halfedge_around_facet_circulator
            he_circ = fit->facet_begin(),
            he_circ_end(he_circ);

    // Iterates on the vector of facet handles
    CDT::Vertex_handle previous, first;
    do {
        CDT::Vertex_handle vh = cdt.insert(he_circ->vertex()->point());
        if(first == 0) {
            first = vh;
        }
        vh->info() = he_circ;
        if(previous != 0 && previous != vh) {
            cdt.insert_constraint(previous, vh);
        }
        previous = vh;
    } while( ++he_circ != he_circ_end );
    cdt.insert_constraint(previous, first);

    // sets mark is_external
    for(CDT::All_faces_iterator
        afit = cdt.all_faces_begin(),
        end = cdt.all_faces_end();
        afit != end; ++afit)
    {
        afit->info().is_external = false;
    }
    //check if the facet is external or internal
    std::queue<CDT::Face_handle> face_queue;
    face_queue.push(cdt.infinite_vertex()->face());
    while(! face_queue.empty() ) {
        CDT::Face_handle fh = face_queue.front();
        face_queue.pop();
        if(fh->info().is_external) continue;
        fh->info().is_external = true;
        for(int i = 0; i <3; ++i) {
            if(!cdt.is_constrained(std::make_pair(fh, i)))
            {
                face_queue.push(fh->neighbor(i));
            }
        }
    }

    //iterates on the internal faces to add the vertices to the positions vector
    for(CDT::Finite_faces_iterator
        ffit = cdt.finite_faces_begin(),
        end = cdt.finite_faces_end();
        ffit != end; ++ffit)
    {
        if(ffit->info().is_external)
            continue;
        //Add Colors
        for(int i = 0; i<3; ++i)
        {
            const int this_patch_id = fit->patch_id();

                color_facets_selected.push_back(colors_[this_patch_id].lighter(120).redF());
                color_facets_selected.push_back(colors_[this_patch_id].lighter(120).greenF());
                color_facets_selected.push_back(colors_[this_patch_id].lighter(120).blueF());

                color_facets_selected.push_back(colors_[this_patch_id].lighter(120).redF());
                color_facets_selected.push_back(colors_[this_patch_id].lighter(120).greenF());
                color_facets_selected.push_back(colors_[this_patch_id].lighter(120).blueF());

                color_facets.push_back(colors_[this_patch_id].redF());
                color_facets.push_back(colors_[this_patch_id].greenF());
                color_facets.push_back(colors_[this_patch_id].blueF());

                color_facets.push_back(colors_[this_patch_id].redF());
                color_facets.push_back(colors_[this_patch_id].greenF());
                color_facets.push_back(colors_[this_patch_id].blueF());


        }
    }
}

#include <QObject>
#include <QMenu>
#include <QAction>
#include <CGAL/gl_render.h>


void
Scene_polyhedron_item::initialize_buffers(Viewer_interface* viewer) const
{
    //vao containing the data for the unselected facets
    {
        program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
        program->bind();
        //flat
        vaos[0]->bind();
        buffers[0].bind();
        buffers[0].allocate(positions_facets.data(),
                            static_cast<int>(positions_facets.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        buffers[0].release();



        buffers[1].bind();
        buffers[1].allocate(normals_flat.data(),
                            static_cast<int>(normals_flat.size()*sizeof(float)));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_FLOAT,0,3);
        buffers[1].release();

        buffers[2].bind();
        buffers[2].allocate(color_facets.data(),
                            static_cast<int>(color_facets.size()*sizeof(float)));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_FLOAT,0,3);
        buffers[2].release();
        vaos[0]->release();

        //gouraud
        vaos[4]->bind();
        buffers[0].bind();
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        buffers[0].release();

        buffers[10].bind();
        buffers[10].allocate(normals_gouraud.data(),
                            static_cast<int>(normals_gouraud.size()*sizeof(float)));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_FLOAT,0,3);
        buffers[10].release();

        buffers[2].bind();
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_FLOAT,0,3);
        buffers[2].release();
        vaos[4]->release();

        program->release();

    }
    //vao containing the data for the unselected lines
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        vaos[1]->bind();

        buffers[3].bind();
        buffers[3].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        buffers[3].release();

        buffers[4].bind();
        buffers[4].allocate(color_lines.data(),
                            static_cast<int>(color_lines.size()*sizeof(float)));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_FLOAT,0,3);
        buffers[4].release();
        program->release();

        vaos[1]->release();

    }
    //vao containing the data for the selected facets
    {

        program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
        program->bind();

        //flat
        vaos[2]->bind();
        buffers[5].bind();
        buffers[5].allocate(positions_facets.data(),
                            static_cast<int>(positions_facets.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        buffers[5].release();

        buffers[6].bind();
        buffers[6].allocate(normals_flat.data(),
                            static_cast<int>(normals_flat.size()*sizeof(float)));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_FLOAT,0,3);
        buffers[6].release();


        buffers[7].bind();
        buffers[7].allocate(color_facets_selected.data(),
                            static_cast<int>(color_facets_selected.size()*sizeof(float)));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_FLOAT,0,3);
        buffers[7].release();
        vaos[2]->release();

        //gouraud
        vaos[5]->bind();
        buffers[5].bind();
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        buffers[5].release();

        buffers[10].bind();
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_FLOAT,0,3);
        buffers[10].release();


        buffers[7].bind();
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_FLOAT,0,3);
        buffers[7].release();
        vaos[5]->release();

        program->release();
    }
    //vao containing the data for the selected lines
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        vaos[3]->bind();
        program->bind();
        buffers[8].bind();
        buffers[8].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        buffers[8].release();

        buffers[9].bind();
        buffers[9].allocate(color_lines_selected.data(),
                            static_cast<int>(color_lines_selected.size()*sizeof(float)));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_FLOAT,0,3);
        buffers[9].release();
        program->release();

        vaos[3]->release();
    }
    nb_lines = positions_lines.size();
    positions_lines.resize(0);
    std::vector<float>(positions_lines).swap(positions_lines);
    nb_facets = positions_facets.size();
    positions_facets.resize(0);
    std::vector<float>(positions_facets).swap(positions_facets);


    color_lines_selected.resize(0);
    std::vector<float>(color_lines_selected).swap(color_lines_selected);
    color_facets_selected.resize(0);
    std::vector<float>(color_facets_selected).swap(color_facets_selected);
    color_lines.resize(0);
    std::vector<float>(color_lines).swap(color_lines);
    color_facets.resize(0);
    std::vector<float>(color_facets).swap(color_facets);
    normals_flat.resize(0);
    std::vector<float>(normals_flat).swap(normals_flat);
    normals_gouraud.resize(0);
    std::vector<float>(normals_gouraud).swap(normals_gouraud);
    are_buffers_filled = true;
}

void
Scene_polyhedron_item::compute_normals_and_vertices(void) const
{
    positions_facets.resize(0);
    positions_lines.resize(0);
    normals_flat.resize(0);
    normals_gouraud.resize(0);


    //Facets
    typedef Polyhedron::Traits	    Kernel;
    typedef Kernel::Point_3	    Point;
    typedef Kernel::Vector_3	    Vector;
    typedef Polyhedron::Facet_iterator Facet_iterator;
    typedef Polyhedron::Halfedge_around_facet_circulator HF_circulator;



    Facet_iterator f = poly->facets_begin();

    for(f = poly->facets_begin();
        f != poly->facets_end();
        f++)
    {

        if(!is_Triangle)
            triangulate_facet(f);
        else
        {
            HF_circulator he = f->facet_begin();
            HF_circulator end = he;
            CGAL_For_all(he,end)
            {

                // // If Flat shading:1 normal per polygon added once per vertex

                Vector n = CGAL::Polygon_mesh_processing::compute_face_normal(f, *poly);
                normals_flat.push_back(n.x());
                normals_flat.push_back(n.y());
                normals_flat.push_back(n.z());


                //// If Gouraud shading: 1 normal per vertex

                n = CGAL::Polygon_mesh_processing::compute_vertex_normal(he->vertex(), *poly);
                normals_gouraud.push_back(n.x());
                normals_gouraud.push_back(n.y());
                normals_gouraud.push_back(n.z());

                //position
                const Point& p = he->vertex()->point();
                positions_facets.push_back(p.x());
                positions_facets.push_back(p.y());
                positions_facets.push_back(p.z());
                positions_facets.push_back(1.0);
            }
        }
    }
    //Lines
    typedef Kernel::Point_3		Point;
    typedef Polyhedron::Edge_iterator	Edge_iterator;

    Edge_iterator he;
    if(!show_only_feature_edges_m) {
        for(he = poly->edges_begin();
            he != poly->edges_end();
            he++)
        {
            if(he->is_feature_edge()) continue;
            const Point& a = he->vertex()->point();
            const Point& b = he->opposite()->vertex()->point();
            positions_lines.push_back(a.x());
            positions_lines.push_back(a.y());
            positions_lines.push_back(a.z());
            positions_lines.push_back(1.0);

            positions_lines.push_back(b.x());
            positions_lines.push_back(b.y());
            positions_lines.push_back(b.z());
            positions_lines.push_back(1.0);

        }
    }
    for(he = poly->edges_begin();
        he != poly->edges_end();
        he++)
    {
        if(!he->is_feature_edge()) continue;
        const Point& a = he->vertex()->point();
        const Point& b = he->opposite()->vertex()->point();

        positions_lines.push_back(a.x());
        positions_lines.push_back(a.y());
        positions_lines.push_back(a.z());
        positions_lines.push_back(1.0);

        positions_lines.push_back(b.x());
        positions_lines.push_back(b.y());
        positions_lines.push_back(b.z());
        positions_lines.push_back(1.0);


    }


    //set the colors
    compute_colors();
}

void
Scene_polyhedron_item::compute_colors() const
{
    color_lines.resize(0);
    color_facets.resize(0);
    color_lines_selected.resize(0);
    color_facets_selected.resize(0);
    //Facets
    typedef Polyhedron::Facet_iterator Facet_iterator;
    typedef Polyhedron::Halfedge_around_facet_circulator HF_circulator;



    // int patch_id = -1;
   // Facet_iterator f = poly->facets_begin();

    for(Facet_iterator f = poly->facets_begin();
        f != poly->facets_end();
        f++)
    {
        if(!is_Triangle)
            triangulate_facet_color(f);
        else
        {
            HF_circulator he = f->facet_begin();
            HF_circulator end = he;
            CGAL_For_all(he,end)
            {

                const int this_patch_id = f->patch_id();

                    color_facets_selected.push_back(colors_[this_patch_id].lighter(120).redF());
                    color_facets_selected.push_back(colors_[this_patch_id].lighter(120).greenF());
                    color_facets_selected.push_back(colors_[this_patch_id].lighter(120).blueF());

                    color_facets.push_back(colors_[this_patch_id].redF());
                    color_facets.push_back(colors_[this_patch_id].greenF());
                    color_facets.push_back(colors_[this_patch_id].blueF());

            }

        }
    }
    //Lines
    typedef Polyhedron::Edge_iterator	Edge_iterator;

    Edge_iterator he;
    if(!show_only_feature_edges_m) {
        for(he = poly->edges_begin();
            he != poly->edges_end();
            he++)
        {
            if(he->is_feature_edge()) continue;

                color_lines_selected.push_back(0.0);
                color_lines_selected.push_back(0.0);
                color_lines_selected.push_back(0.0);

                color_lines_selected.push_back(0.0);
                color_lines_selected.push_back(0.0);
                color_lines_selected.push_back(0.0);

                color_lines.push_back(this->color().lighter(50).redF());
                color_lines.push_back(this->color().lighter(50).greenF());
                color_lines.push_back(this->color().lighter(50).blueF());

                color_lines.push_back(this->color().lighter(50).redF());
                color_lines.push_back(this->color().lighter(50).greenF());
                color_lines.push_back(this->color().lighter(50).blueF());




        }
    }
    for(he = poly->edges_begin();
        he != poly->edges_end();
        he++)
    {
        if(!he->is_feature_edge()) continue;
        color_lines.push_back(1.0);
        color_lines.push_back(0.0);
        color_lines.push_back(0.0);

        color_lines.push_back(1.0);
        color_lines.push_back(0.0);
        color_lines.push_back(0.0);

        color_lines_selected.push_back(1.0);
        color_lines_selected.push_back(0.0);
        color_lines_selected.push_back(0.0);

        color_lines_selected.push_back(1.0);
        color_lines_selected.push_back(0.0);
        color_lines_selected.push_back(0.0);
    }
}

Scene_polyhedron_item::Scene_polyhedron_item()
    : Scene_item(11,6),
      is_Triangle(true),
      poly(new Polyhedron),
      show_only_feature_edges_m(false),
      facet_picking_m(false),
      erase_next_picked_facet_m(false),
      plugin_has_set_color_vector_m(false)
{
    cur_shading=FlatPlusEdges;
    is_selected = true;
    nb_facets = 0;
    nb_lines = 0;
    init();

}

Scene_polyhedron_item::Scene_polyhedron_item(Polyhedron* const p)
    : Scene_item(11,6),
      is_Triangle(true),
      poly(p),
      show_only_feature_edges_m(false),
      facet_picking_m(false),
      erase_next_picked_facet_m(false),
      plugin_has_set_color_vector_m(false)
{
    cur_shading=FlatPlusEdges;
    is_selected = true;
    nb_facets = 0;
    nb_lines = 0;
    init();
    invalidate_buffers();
}

Scene_polyhedron_item::Scene_polyhedron_item(const Polyhedron& p)
    : Scene_item(11,6),
      is_Triangle(true),
      poly(new Polyhedron(p)),
      show_only_feature_edges_m(false),
      facet_picking_m(false),
      erase_next_picked_facet_m(false),
      plugin_has_set_color_vector_m(false)
{
    cur_shading=FlatPlusEdges;
    is_selected=true;
    init();
    nb_facets = 0;
    nb_lines = 0;
    invalidate_buffers();
}

Scene_polyhedron_item::~Scene_polyhedron_item()
{
    delete_aabb_tree(this);
    delete poly;
}

#include "Color_map.h"

void
Scene_polyhedron_item::
init()
{
    typedef Polyhedron::Facet_iterator Facet_iterator;

    if ( !plugin_has_set_color_vector_m )
    {
        // Fill indices map and get max subdomain value
        int max = 0;
        for(Facet_iterator fit = poly->facets_begin(), end = poly->facets_end() ;
            fit != end; ++fit)
        {
            max = (std::max)(max, fit->patch_id());
        }

        colors_.resize(0);
        compute_color_map(this->color(), max + 1,
                          std::back_inserter(colors_));
    }

  volume=-std::numeric_limits<double>::infinity();
  area=-std::numeric_limits<double>::infinity();
  if (poly->is_pure_triangle())
  {
    // compute the volume if the polyhedron is closed
    if (poly->is_closed())
    {
      volume=0;
      Polyhedron::Vertex::Point p(0,0,0);
      Q_FOREACH(Polyhedron::Face_handle fh, faces(*poly))
      {
        volume+=CGAL::volume( p,
                    fh->halfedge()->vertex()->point(),
                    fh->halfedge()->next()->vertex()->point(),
                    fh->halfedge()->prev()->vertex()->point() );
      }
    }

    // compute the surface area
    area=0;
    Q_FOREACH(Polyhedron::Face_handle fh, faces(*poly))
    {
      area+=std::sqrt( CGAL::squared_area(
              fh->halfedge()->vertex()->point(),
              fh->halfedge()->next()->vertex()->point(),
              fh->halfedge()->prev()->vertex()->point() )
            );
    }
  }
}


Scene_polyhedron_item*
Scene_polyhedron_item::clone() const {
    return new Scene_polyhedron_item(*poly);}

// Load polyhedron from .OFF file
bool
Scene_polyhedron_item::load(std::istream& in)
{


    in >> *poly;

    if ( in && !isEmpty() )
    {
        invalidate_buffers();
        return true;
    }
    return false;
}

// Write polyhedron to .OFF file
bool
Scene_polyhedron_item::save(std::ostream& out) const
{
  out.precision(17);
    out << *poly;
    return (bool) out;
}

QString
Scene_polyhedron_item::toolTip() const
{
    if(!poly)
        return QString();

  QString str =
         QObject::tr("<p>Polyhedron <b>%1</b> (mode: %5, color: %6)</p>"
                       "<p>Number of vertices: %2<br />"
                       "Number of edges: %3<br />"
                     "Number of facets: %4")
            .arg(this->name())
            .arg(poly->size_of_vertices())
            .arg(poly->size_of_halfedges()/2)
            .arg(poly->size_of_facets())
            .arg(this->renderingModeName())
            .arg(this->color().name());
  if (volume!=-std::numeric_limits<double>::infinity())
    str+=QObject::tr("<br />Volume: %1").arg(volume);
  if (area!=-std::numeric_limits<double>::infinity())
    str+=QObject::tr("<br />Area: %1").arg(area);
  str+="</p>";

  return str;
}

QMenu* Scene_polyhedron_item::contextMenu()
{
    const char* prop_name = "Menu modified by Scene_polyhedron_item.";

    QMenu* menu = Scene_item::contextMenu();

    // Use dynamic properties:
    // http://doc.qt.io/qt-5/qobject.html#property
    bool menuChanged = menu->property(prop_name).toBool();

    if(!menuChanged) {

        QAction* actionShowOnlyFeatureEdges =
                menu->addAction(tr("Show only &feature edges"));
        actionShowOnlyFeatureEdges->setCheckable(true);
        actionShowOnlyFeatureEdges->setObjectName("actionShowOnlyFeatureEdges");
        connect(actionShowOnlyFeatureEdges, SIGNAL(toggled(bool)),
                this, SLOT(show_only_feature_edges(bool)));

        QAction* actionPickFacets =
                menu->addAction(tr("Facets picking"));
        actionPickFacets->setCheckable(true);
        actionPickFacets->setObjectName("actionPickFacets");
        connect(actionPickFacets, SIGNAL(toggled(bool)),
                this, SLOT(enable_facets_picking(bool)));

        QAction* actionEraseNextFacet =
                menu->addAction(tr("Erase next picked facet"));
        actionEraseNextFacet->setCheckable(true);
        actionEraseNextFacet->setObjectName("actionEraseNextFacet");
        connect(actionEraseNextFacet, SIGNAL(toggled(bool)),
                this, SLOT(set_erase_next_picked_facet(bool)));

        menu->setProperty(prop_name, true);
    }
    QAction* action = menu->findChild<QAction*>("actionPickFacets");
    if(action) action->setChecked(facet_picking_m);
    action = menu->findChild<QAction*>("actionEraseNextFacet");
    if(action) action->setChecked(erase_next_picked_facet_m);
    return menu;
}

void Scene_polyhedron_item::show_only_feature_edges(bool b)
{
    show_only_feature_edges_m = b;
  Q_EMIT itemChanged();
}

void Scene_polyhedron_item::enable_facets_picking(bool b)
{
    facet_picking_m = b;
}

void Scene_polyhedron_item::set_erase_next_picked_facet(bool b)
{
    if(b) { facet_picking_m = true; } // automatically activate facet_picking
    erase_next_picked_facet_m = b;
}

void Scene_polyhedron_item::draw(Viewer_interface* viewer) const {
    if(!are_buffers_filled)
    {
        is_Triangulated();
        compute_normals_and_vertices();
        initialize_buffers(viewer);
    }


    if(!is_selected && renderingMode() == Flat)
        vaos[0]->bind();
    else if(!is_selected && renderingMode() == Gouraud)
    {
        vaos[4]->bind();
    }
    else if (is_selected && renderingMode() == Gouraud)
    {
        vaos[5]->bind();
    }
    else
    {
        vaos[2]->bind();
    }
    attrib_buffers(viewer, PROGRAM_WITH_LIGHT);
    program = getShaderProgram(PROGRAM_WITH_LIGHT);
    program->bind();
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(nb_facets/4));
    program->release();
    if(!is_selected && renderingMode() == Flat)
        vaos[0]->release();
    else if(!is_selected && renderingMode() == Gouraud)
        vaos[4]->release();
    else if (is_selected && renderingMode() == Gouraud)
        vaos[5]->release();
    else
        vaos[2]->release();
}

// Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
void Scene_polyhedron_item::draw_edges(Viewer_interface* viewer) const {
    if(!are_buffers_filled)
    {
        is_Triangulated();
        compute_normals_and_vertices();
        initialize_buffers(viewer);
    }

    if(!is_selected)
    {
        vaos[1]->bind();
    }
    else
    {
        vaos[3]->bind();
    }
    attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    program->bind();
    //draw the edges
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(nb_lines/4));
    program->release();
    if(!is_selected)
        vaos[1]->release();
    else
        vaos[3]->release();
}

void
Scene_polyhedron_item::draw_points(Viewer_interface* viewer) const {
    if(!are_buffers_filled)
    {
        is_Triangulated();
        compute_normals_and_vertices();
        initialize_buffers(viewer);
    }

    vaos[1]->bind();
    attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    program->bind();
    //draw the points
    viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(nb_lines/4));
    // Clean-up
    program->release();
    vaos[1]->release();
}

Polyhedron*
Scene_polyhedron_item::polyhedron()       { return poly; }
const Polyhedron*
Scene_polyhedron_item::polyhedron() const { return poly; }

bool
Scene_polyhedron_item::isEmpty() const {
    return (poly == 0) || poly->empty();
}

Scene_polyhedron_item::Bbox
Scene_polyhedron_item::bbox() const {
    const Kernel::Point_3& p = *(poly->points_begin());
    CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    for(Polyhedron::Point_iterator it = poly->points_begin();
        it != poly->points_end();
        ++it) {
        bbox = bbox + it->bbox();
    }
    return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
                bbox.xmax(),bbox.ymax(),bbox.zmax());
}


void
Scene_polyhedron_item::
invalidate_buffers()
{
  Q_EMIT item_is_about_to_be_changed();
    delete_aabb_tree(this);
    init();
    Base::invalidate_buffers();
    is_Triangulated();
    compute_normals_and_vertices();
    are_buffers_filled = false;

}
void
Scene_polyhedron_item::
contextual_changed()
{
    prev_shading = cur_shading;
    cur_shading = renderingMode();
    if(prev_shading != cur_shading)
        if(cur_shading == Flat || cur_shading == FlatPlusEdges ||cur_shading == Gouraud)
        {
            //Change the normals
           // invalidate_buffers();
        }

}
void
Scene_polyhedron_item::selection_changed(bool p_is_selected)
{
    if(p_is_selected != is_selected)
    {
        is_selected = p_is_selected;
         //are_buffers_filled = false;
    }

}

void
Scene_polyhedron_item::setColor(QColor c)
{
  // reset patch ids
  if (colors_.size()>2 || plugin_has_set_color_vector_m)
  {
    BOOST_FOREACH(Polyhedron::Facet_handle fh, faces(*poly))
      fh->set_patch_id(1);
    colors_[1]=c;
  }
  Scene_item::setColor(c);
}

void
Scene_polyhedron_item::select(double orig_x,
                              double orig_y,
                              double orig_z,
                              double dir_x,
                              double dir_y,
                              double dir_z)
{
    if(facet_picking_m) {
        typedef Input_facets_AABB_tree Tree;
        typedef Tree::Object_and_primitive_id Object_and_primitive_id;

        Tree* aabb_tree = get_aabb_tree(this);
        if(aabb_tree) {
            const Kernel::Point_3 ray_origin(orig_x, orig_y, orig_z);
            const Kernel::Vector_3 ray_dir(dir_x, dir_y, dir_z);
            const Kernel::Ray_3 ray(ray_origin, ray_dir);

            typedef std::list<Object_and_primitive_id> Intersections;
            Intersections intersections;

            aabb_tree->all_intersections(ray, std::back_inserter(intersections));

            Intersections::iterator closest = intersections.begin();
            if(closest != intersections.end()) {
                const Kernel::Point_3* closest_point =
                        CGAL::object_cast<Kernel::Point_3>(&closest->first);

                for(Intersections::iterator
                    it = boost::next(intersections.begin()),
                    end = intersections.end();
                    it != end; ++it)
                {
                    if(! closest_point) {
                        closest = it;
                    }
                    else {
                        const Kernel::Point_3* it_point =
                                CGAL::object_cast<Kernel::Point_3>(&it->first);
                        if(it_point &&
                                (ray_dir * (*it_point - *closest_point)) < 0)
                        {
                            closest = it;
                            closest_point = it_point;
                        }
                    }
                }
                if(closest_point) {
                    Polyhedron::Facet_handle selected_fh = closest->second;

                    // The computation of the nearest vertex may be costly.  Only
                    // do it if some objects are connected to the signal
                    // 'selected_vertex'.
                    if(QObject::receivers(SIGNAL(selected_vertex(void*))) > 0)
                    {
                        Polyhedron::Halfedge_around_facet_circulator
                                he_it = selected_fh->facet_begin(),
                                around_end = he_it;

                        Polyhedron::Vertex_handle v = he_it->vertex(), nearest_v = v;

                        Kernel::FT sq_dist = CGAL::squared_distance(*closest_point,
                                                                    v->point());

                        while(++he_it != around_end) {
                            v = he_it->vertex();
                            Kernel::FT new_sq_dist = CGAL::squared_distance(*closest_point,
                                                                            v->point());
                            if(new_sq_dist < sq_dist) {
                                sq_dist = new_sq_dist;
                                nearest_v = v;
                            }
                        }

            Q_EMIT selected_vertex((void*)(&*nearest_v));
                    }

                    if(QObject::receivers(SIGNAL(selected_edge(void*))) > 0
                            || QObject::receivers(SIGNAL(selected_halfedge(void*))) > 0)
                    {
                        Polyhedron::Halfedge_around_facet_circulator
                                he_it = selected_fh->facet_begin(),
                                around_end = he_it;

                        Polyhedron::Halfedge_handle nearest_h = he_it;
                        Kernel::FT sq_dist = CGAL::squared_distance(*closest_point,
                                                                    Kernel::Segment_3(he_it->vertex()->point(), he_it->opposite()->vertex()->point()));

                        while(++he_it != around_end) {
                            Kernel::FT new_sq_dist = CGAL::squared_distance(*closest_point,
                                                                            Kernel::Segment_3(he_it->vertex()->point(), he_it->opposite()->vertex()->point()));
                            if(new_sq_dist < sq_dist) {
                                sq_dist = new_sq_dist;
                                nearest_h = he_it;
                            }
                        }

            Q_EMIT selected_halfedge((void*)(&*nearest_h));
            Q_EMIT selected_edge((void*)(std::min)(&*nearest_h, &*nearest_h->opposite()));
                    }

          Q_EMIT selected_facet((void*)(&*selected_fh));
                    if(erase_next_picked_facet_m) {
                        polyhedron()->erase_facet(selected_fh->halfedge());
                        polyhedron()->normalize_border();
                        //set_erase_next_picked_facet(false);
                        invalidate_buffers();
            Q_EMIT itemChanged();
                    }
                }
            }
        }
    }
    Base::select(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z);
}

void Scene_polyhedron_item::update_vertex_indices()
{
    std::size_t id=0;
    for (Polyhedron::Vertex_iterator vit = polyhedron()->vertices_begin(),
         vit_end = polyhedron()->vertices_end(); vit != vit_end; ++vit)
    {
        vit->id()=id++;
    }
}
void Scene_polyhedron_item::update_facet_indices()
{
    std::size_t id=0;
    for (Polyhedron::Facet_iterator  fit = polyhedron()->facets_begin(),
         fit_end = polyhedron()->facets_end(); fit != fit_end; ++fit)
    {
        fit->id()=id++;
    }
}
void Scene_polyhedron_item::update_halfedge_indices()
{
    std::size_t id=0;
    for (Polyhedron::Halfedge_iterator hit = polyhedron()->halfedges_begin(),
         hit_end = polyhedron()->halfedges_end(); hit != hit_end; ++hit)
    {
        hit->id()=id++;
    }
}
void Scene_polyhedron_item::invalidate_aabb_tree()
{
  delete_aabb_tree(this);
}


