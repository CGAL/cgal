#include "Scene_combinatorial_map_item.h"
#include "Scene_polyhedron_item.h"
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Viewer_interface.h>

#include <QObject>
#include <QMenu>
#include <QAction>
#include <QtDebug>
#include <QApplication>
#include <QDebug>
#include <QKeyEvent>
#include <CGAL/corefinement_operations.h>

struct Scene_combinatorial_map_item_priv
{
  Scene_combinatorial_map_item_priv(CGAL::Three::Scene_interface* scene, void* ad_A, Scene_combinatorial_map_item* parent)
  :last_known_scene(scene),volume_to_display(0),exportSelectedVolume(NULL),address_of_A(ad_A)
  {
    item = parent;
    item->are_buffers_filled = false;
    nb_points = 0;
    nb_lines =0;
    nb_facets =0;
  }

  Kernel::Vector_3 compute_face_normal(Combinatorial_map_3::Dart_const_handle adart) const;

  template <class Predicate>
  void export_as_polyhedron(Predicate,const QString&) const;
  void initialize_buffers(CGAL::Three::Viewer_interface *viewer) const;
  void compute_elements(void) const;

  Scene_combinatorial_map_item *item;
  enum VAOs {
      Edges = 0,
      Points,
      Facets,
      NbOfVaos
  };
  enum VBOs {
      Edges_vertices = 0,
      Points_vertices,
      Facets_vertices,
      Facets_normals,
      NbOfVbos
  };

  CGAL::Three::Scene_interface* last_known_scene;
  std::size_t volume_to_display;
  QAction* exportSelectedVolume;
  void* address_of_A;

  mutable QOpenGLShaderProgram *program;
  mutable std::vector<double> positions_lines;
  mutable std::vector<double> positions_points;
  mutable std::vector<double> positions_facets;
  mutable std::vector<double> normals;
  mutable std::size_t nb_lines;
  mutable std::size_t nb_points;
  mutable std::size_t nb_facets;
};
Scene_combinatorial_map_item::Scene_combinatorial_map_item(CGAL::Three::Scene_interface* scene,void* address){m_combinatorial_map=NULL; d = new Scene_combinatorial_map_item_priv(scene, address, this);}
Scene_combinatorial_map_item::~Scene_combinatorial_map_item(){if (m_combinatorial_map!=NULL) delete m_combinatorial_map; delete d;}

Scene_combinatorial_map_item* Scene_combinatorial_map_item::clone() const{return NULL;}

Kernel::Vector_3 Scene_combinatorial_map_item_priv::compute_face_normal(Combinatorial_map_3::Dart_const_handle adart) const
{
    typedef Combinatorial_map_3::Dart_of_orbit_const_range<1> Dart_in_facet_range;
    typedef Kernel::Vector_3 Vector_3;
    Vector_3 normal = CGAL::NULL_VECTOR;

    Dart_in_facet_range vertices=item->combinatorial_map().darts_of_orbit<1>(adart);
    Kernel::Point_3 points[3];
    int index=0;
    Dart_in_facet_range::const_iterator pit=vertices.begin();
    for (;pit!=vertices.end() && index!=3;++pit,++index ){
        points[index]=pit->attribute<0>()->point();
    }

    if (index!=3) return normal;

    do{
        Vector_3 n = CGAL::cross_product(points[2]-points[1],points[0]-points[1]);
        if (n != Vector_3(0,0,0) )
            normal = normal + (n / std::sqrt(n*n));
        points[0]=points[1];
        points[1]=points[2];
        if ( pit==vertices.end() ) break;

        points[2]=pit->attribute<0>()->point();
        ++pit;
    }while(true);

    return normal == Vector_3(0,0,0)? normal : normal / std::sqrt(normal * normal);
}  

void Scene_combinatorial_map_item::set_next_volume(){
    //Update des vectors faits ici
    ++d->volume_to_display;
    d->volume_to_display=d->volume_to_display%(combinatorial_map().attributes<3>().size()+1);
  are_buffers_filled = false;
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();

    if (d->exportSelectedVolume!=NULL && ( d->volume_to_display==1 || d->volume_to_display==0 ) )
        d->exportSelectedVolume->setEnabled(!d->exportSelectedVolume->isEnabled());
}


template <class Predicate> 
void Scene_combinatorial_map_item_priv::export_as_polyhedron(Predicate pred,const QString& name) const {
    typedef Combinatorial_map_3::Dart_const_handle Dart_handle;
    typedef Combinatorial_map_3::One_dart_per_cell_range<3> One_dart_per_vol_range;
    typedef CGAL::internal::Import_volume_as_polyhedron<Polyhedron::HalfedgeDS> Volume_import_modifier;

    std::vector<Dart_handle> darts;
    One_dart_per_vol_range cell_range=item->combinatorial_map().template one_dart_per_cell<3>();


    for (One_dart_per_vol_range::const_iterator it = cell_range.begin();it!= cell_range.end() ; ++it )
        if ( pred(it) ){
            darts.push_back(it);
            if (Predicate::only_one_run) break;
        }

    if (!darts.empty())
    {
        Volume_import_modifier modifier=Predicate::swap_orientation?
                    Volume_import_modifier(item->combinatorial_map(),darts.begin(),darts.end(),Predicate::swap_orientation):
                    Volume_import_modifier(item->combinatorial_map(),darts.begin(),darts.end());

        Polyhedron* new_poly=new Polyhedron();
        new_poly->delegate(modifier);
        Scene_polyhedron_item* new_item = new Scene_polyhedron_item(new_poly);
        new_item->setName(name);
        last_known_scene->addItem(new_item);
    }
}

struct Select_volume{
    static const bool only_one_run=true;
    static const bool swap_orientation=false;
    Select_volume(std::size_t i):volume_to_select(i),index(0){}
    template <class Dart_handle>
    bool operator() (Dart_handle){
        return ++index==volume_to_select;
    }
private:
    std::size_t volume_to_select;
    std::size_t index;
};

void Scene_combinatorial_map_item::export_current_volume_as_polyhedron() const {
    if (d->volume_to_display==0) return; //no volume selected

    Select_volume predicate(d->volume_to_display);
    d->export_as_polyhedron(predicate,QString("%1_%2").arg(this->name()).arg(d->volume_to_display-1));
}

struct Select_union{
    static const bool only_one_run=false;
    static const bool swap_orientation=true;
    template <class Dart_handle>
    bool operator() (Dart_handle d){ return d->template attribute<3>()->info().outside.size()==2; }
};

struct Select_inter{
    static const bool only_one_run=false;
    static const bool swap_orientation=false;
    template <class Dart_handle>
    bool operator() (Dart_handle d){ return d->template attribute<3>()->info().inside.size()==2; }
};

struct Select_A_minus_B{
    static const bool only_one_run=false;
    static const bool swap_orientation=false;
    Select_A_minus_B(void* address):address_of_A(address){}
    template <class Dart_handle>
    bool operator() (Dart_handle d){
        return d->template attribute<3>()->info().inside.size()==1 &&
                static_cast<void*>(*d->template attribute<3>()->info().inside.begin())==address_of_A;
    }
private:
    void* address_of_A;
};

struct Select_B_minus_A{
    static const bool only_one_run=false;
    static const bool swap_orientation=false;
    Select_B_minus_A(void* address):address_of_A(address){}
    template <class Dart_handle>
    bool operator() (Dart_handle d){
        return d->template attribute<3>()->info().inside.size()==1 &&
                static_cast<void*>(*d->template attribute<3>()->info().inside.begin())!=address_of_A;
    }
private:
    void* address_of_A;
};

void Scene_combinatorial_map_item::export_union_as_polyhedron() const {
    d->export_as_polyhedron(Select_union(),QString("%1_union_%2").arg("A").arg("B"));
}
void Scene_combinatorial_map_item::export_intersection_as_polyhedron() const{
    d->export_as_polyhedron(Select_inter(),QString("%1_inter_%2").arg("A").arg("B"));
}
void Scene_combinatorial_map_item::export_A_minus_B_as_polyhedron() const{
    Select_A_minus_B predicate(d->address_of_A);
    d->export_as_polyhedron(predicate,QString("%1_minus_%2").arg("A").arg("B"));
}
void Scene_combinatorial_map_item::export_B_minus_A_as_polyhedron() const{
    Select_B_minus_A predicate(d->address_of_A);
    d->export_as_polyhedron(predicate,QString("%1_minus_%2").arg("B").arg("A"));
}

QMenu* Scene_combinatorial_map_item::contextMenu()
{
    const char* prop_name = "Menu modified by Scene_combinatorial_map_item.";

    QMenu* menu = Scene_item::contextMenu();

    // Use dynamic properties:
    // http://doc.qt.io/qt-5/qobject.html#property
    bool menuChanged = menu->property(prop_name).toBool();

    if(!menuChanged) {
        QAction* actionSelectNextVolume =
                menu->addAction(tr("Iterate over volumes"));
        actionSelectNextVolume->setObjectName("actionSelectNextVolume");
        connect(actionSelectNextVolume, SIGNAL(triggered()),this, SLOT(set_next_volume()));

        d->exportSelectedVolume =
                menu->addAction(tr("Export current volume as polyhedron"));
        d->exportSelectedVolume->setObjectName("exportSelectedVolume");
        connect(d->exportSelectedVolume, SIGNAL(triggered()),this, SLOT(export_current_volume_as_polyhedron()));
        d->exportSelectedVolume->setEnabled(d->volume_to_display!=0);
        menu->setProperty(prop_name, true);

        if(is_from_corefinement()){
            //Export union as polyhedron
            QAction* exportUnion =
                    menu->addAction(tr("Export union as polyhedron"));
            exportUnion->setObjectName("exportUnion");
            connect(exportUnion, SIGNAL(triggered()),this, SLOT(export_union_as_polyhedron()));

            //Export intersection as polyhedron
            QAction* exportIntersection =
                    menu->addAction(tr("Export intersection as polyhedron"));
            exportIntersection->setObjectName("exportIntersection");
            connect(exportIntersection, SIGNAL(triggered()),this, SLOT(export_intersection_as_polyhedron()));

            //Export A minus B as polyhedron
            QAction* exportAMinusB =
                    menu->addAction(tr("Export A minus B as polyhedron"));
            exportAMinusB->setObjectName("exportAMinusB");
            connect(exportAMinusB, SIGNAL(triggered()),this, SLOT(export_A_minus_B_as_polyhedron()));

            //Export B minus A as polyhedron
            QAction* exportBMinusA =
                    menu->addAction(tr("Export B minus A as polyhedron"));
            exportBMinusA->setObjectName("exportBMinusA");
            connect(exportBMinusA, SIGNAL(triggered()),this, SLOT(export_B_minus_A_as_polyhedron()));

        }
    }
    return menu;
}

bool Scene_combinatorial_map_item::keyPressEvent(QKeyEvent* e){
    if (e->key()==Qt::Key_N){
        set_next_volume();
        return true;
    }
    return false;
}

void Scene_combinatorial_map_item_priv::compute_elements(void) const{
    QApplication::setOverrideCursor(Qt::WaitCursor);

    positions_facets.resize(0);
    normals.resize(0);
    positions_lines.resize(0);
    positions_points.resize(0);

    //Facets
    {
    std::size_t index = 0;
    Combinatorial_map_3::size_type voltreated
      = item->combinatorial_map().get_new_mark();
    Combinatorial_map_3::size_type facetreated
      = item->combinatorial_map().get_new_mark();
    Combinatorial_map_3::Dart_const_range::const_iterator
            darts_it=item->combinatorial_map().darts().begin(), darts_end=item->combinatorial_map().darts().end();
    for( ; darts_it!=darts_end; ++darts_it)
    {
        if ( !item->combinatorial_map().is_marked(darts_it,voltreated) )
        {
            ++index;
            //iterate over all the darts of the volume
            Combinatorial_map_3::Dart_of_cell_const_range<3>::const_iterator
                    vol_it=item->combinatorial_map().darts_of_cell<3>(darts_it).begin(),
                    vol_end=item->combinatorial_map().darts_of_cell<3>(darts_it).end();
            if ( volume_to_display!=0 && index!=volume_to_display )
            {
                //only mark darts if the volume is not the one to display
                for ( ;vol_it!=vol_end; ++vol_it )
                {
                    item->combinatorial_map().mark(vol_it,facetreated);
                    item->combinatorial_map().mark(vol_it, voltreated);
                }
            }
            else
            {
                for ( ;vol_it!=vol_end; ++vol_it )
                {
                    if ( !item->combinatorial_map().is_marked(vol_it,facetreated) )
                    {
                        Kernel::Vector_3 normal = compute_face_normal(vol_it);
                        for(int i=0; i<3; i++)
                        {
                            normals.push_back(normal.x());
                            normals.push_back(normal.y());
                            normals.push_back(normal.z());
                        }

                        //iterate over all darts of facets
                        for ( Combinatorial_map_3::Dart_of_orbit_const_range<1>::const_iterator
                              face_it=item->combinatorial_map().darts_of_orbit<1>(vol_it).begin(),
                              face_end=item->combinatorial_map().darts_of_orbit<1>(vol_it).end();
                              face_it!=face_end; ++face_it)
                        {
                            const Kernel::Point_3& p= face_it->attribute<0>()->point();
                            positions_facets.push_back(p.x());
                            positions_facets.push_back(p.y());
                            positions_facets.push_back(p.z());
                            item->combinatorial_map().mark(face_it,facetreated);
                            item->combinatorial_map().mark(face_it, voltreated);
                        }
                    }
                }
            }
            if ( index==volume_to_display ) break;
        }
    }
    //mark remaining darts to have an O(1) free_mark
    for( ;  darts_it!=darts_end; ++darts_it)
    {
        item->combinatorial_map().mark(darts_it, facetreated);
        item->combinatorial_map().mark(darts_it, voltreated);
    }

    item->combinatorial_map().free_mark(facetreated);
    item->combinatorial_map().free_mark(voltreated);
    }

    //edges
    {

        typedef Combinatorial_map_3::One_dart_per_cell_range<1,3> Edge_darts;
        Edge_darts darts=item->combinatorial_map().one_dart_per_cell<1>();
        for (Edge_darts::const_iterator dit=darts.begin();dit!=darts.end();++dit){
            CGAL_assertion(!item->combinatorial_map().is_free(dit,1));
            const Kernel::Point_3& a = dit->attribute<0>()->point();
            const Kernel::Point_3& b = dit->beta(1)->attribute<0>()->point();
            positions_lines.push_back(a.x());
            positions_lines.push_back(a.y());
            positions_lines.push_back(a.z());

            positions_lines.push_back(b.x());
            positions_lines.push_back(b.y());
            positions_lines.push_back(b.z());

        }
    }

    //points
    {
        typedef Combinatorial_map_3::Attribute_const_range<0>::type Point_range;
        const Point_range& points=item->combinatorial_map().attributes<0>();
        for(Point_range::const_iterator pit=boost::next(points.begin());pit!=points.end();++pit){
            const Kernel::Point_3& p=pit->point();
            positions_points.push_back(p.x());
            positions_points.push_back(p.y());
            positions_points.push_back(p.z());
        }

    }
    QApplication::restoreOverrideCursor();
}


void Scene_combinatorial_map_item_priv::initialize_buffers(CGAL::Three::Viewer_interface *viewer) const
{
    //vao for the edges
    {
        program = item->getShaderProgram(Scene_combinatorial_map_item::PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();

        item->vaos[Scene_combinatorial_map_item_priv::Edges]->bind();
        item->buffers[Scene_combinatorial_map_item_priv::Edges_vertices].bind();
        item->buffers[Scene_combinatorial_map_item_priv::Edges_vertices].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        item->buffers[Scene_combinatorial_map_item_priv::Edges_vertices].release();
        nb_lines = positions_lines.size();
        positions_lines.resize(0);
        std::vector<double>(positions_lines).swap(positions_lines);
        item->vaos[Scene_combinatorial_map_item_priv::Edges]->release();
        program->release();
    }
    //vao for the points
    {
        program = item->getShaderProgram(Scene_combinatorial_map_item::PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();

        item->vaos[Scene_combinatorial_map_item_priv::Points]->bind();
        item->buffers[Scene_combinatorial_map_item_priv::Points_vertices].bind();
        item->buffers[Scene_combinatorial_map_item_priv::Points_vertices].allocate(positions_points.data(),
                            static_cast<int>(positions_points.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        item->buffers[Scene_combinatorial_map_item_priv::Points_vertices].release();
        item->vaos[Scene_combinatorial_map_item_priv::Points]->release();
        nb_points = positions_points.size();
        positions_points.resize(0);
        std::vector<double>(positions_points).swap(positions_points);
        program->release();
    }
    //vao for the facets
    {
        program = item->getShaderProgram(Scene_combinatorial_map_item::PROGRAM_WITH_LIGHT, viewer);
        program->bind();

        item->vaos[Scene_combinatorial_map_item_priv::Facets]->bind();
        item->buffers[Scene_combinatorial_map_item_priv::Facets_vertices].bind();
        item->buffers[Scene_combinatorial_map_item_priv::Facets_vertices].allocate(positions_facets.data(),
                            static_cast<int>(positions_facets.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        item->buffers[Scene_combinatorial_map_item_priv::Facets_vertices].release();

        item->buffers[Scene_combinatorial_map_item_priv::Facets_normals].bind();
        item->buffers[Scene_combinatorial_map_item_priv::Facets_normals].allocate(normals.data(),
                            static_cast<int>(normals.size()*sizeof(double)));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_DOUBLE,0,3);
        item->buffers[Scene_combinatorial_map_item_priv::Facets_normals].release();
        nb_facets = positions_facets.size();
        positions_facets.resize(0);
        std::vector<double>(positions_facets).swap(positions_facets);
        normals.resize(0);
        std::vector<double>(normals).swap(normals);
        item->vaos[Scene_combinatorial_map_item_priv::Facets]->release();
        program->release();
    }
    item->are_buffers_filled = true;


}

bool Scene_combinatorial_map_item::isEmpty() const {return combinatorial_map().number_of_darts()==0;}

void
Scene_combinatorial_map_item::compute_bbox() const {
    typedef Combinatorial_map_3::Attribute_const_range<0>::type Point_range;
    const Point_range& points=combinatorial_map().attributes<0>();
    CGAL::Bbox_3 bbox=points.begin()->point().bbox();
    for(Point_range::const_iterator pit=boost::next(points.begin());pit!=points.end();++pit)
        bbox=bbox+pit->point().bbox();
    _bbox = Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
                bbox.xmax(),bbox.ymax(),bbox.zmax());
}


QString Scene_combinatorial_map_item::toolTip() const{ 
    if(!m_combinatorial_map)
        return QString();

    std::vector<unsigned int> cells(5);
    for (unsigned int i=0; i<=4; ++i)
        cells[i]=i;
    std::vector<unsigned int> res = combinatorial_map().count_cells(cells);
    if (d->volume_to_display==0)
        return QObject::tr("<p>Combinatorial_map_3 <b>%1</b> (mode: %8, color: %9)</p>"
                           "<p>Number of darts: %2<br />"
                           "Number of vertices: %3<br />"
                           "Number of edges: %4<br />"
                           "Number of facets: %5<br />"
                           "Number of volumes: %6<br />"
                           "Number of connected components: %7</p>")
                .arg(this->name())
                .arg(combinatorial_map().number_of_darts())
                .arg(res[0])
                .arg(res[1])
                .arg(res[2])
                .arg(res[3])
                .arg(res[4])
                .arg(this->renderingModeName())
                .arg(this->color().name());
    return QObject::tr("<p>Combinatorial_map_3 <b>%1</b> (mode: %8, color: %9)</p>"
                       "<p>Number of darts: %2<br />"
                       "Number of vertices: %3<br />"
                       "Number of edges: %4<br />"
                       "Number of facets: %5<br />"
                       "Number of volumes: %6<br />"
                       "Number of connected components: %7 <br />"
                       "Currently Displaying facets of volume: %10 </p>")
            .arg(this->name())
            .arg(combinatorial_map().number_of_darts())
            .arg(res[0])
            .arg(res[1])
            .arg(res[2])
            .arg(res[3])
            .arg(res[4])
            .arg(this->renderingModeName())
            .arg(this->color().name())
            .arg(d->volume_to_display-1);
}


void Scene_combinatorial_map_item::draw(CGAL::Three::Viewer_interface* viewer) const
{
    if(!are_buffers_filled)
    {
        d->compute_elements();
        d->initialize_buffers(viewer);
    }
    vaos[Scene_combinatorial_map_item_priv::Facets]->bind();
    d->program=getShaderProgram(PROGRAM_WITH_LIGHT);
    attribBuffers(viewer,PROGRAM_WITH_LIGHT);
    d->program->bind();
    d->program->setAttributeValue("colors", this->color());
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(d->nb_facets/3));
    vaos[Scene_combinatorial_map_item_priv::Facets]->release();
    d->program->release();
}
 void Scene_combinatorial_map_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
{
     if(!are_buffers_filled)
     {
         d->compute_elements();
         d->initialize_buffers(viewer);
     }
     vaos[Scene_combinatorial_map_item_priv::Edges]->bind();
     d->program=getShaderProgram(PROGRAM_WITHOUT_LIGHT);
     attribBuffers(viewer,PROGRAM_WITHOUT_LIGHT);
     d->program->bind();
     d->program->setAttributeValue("colors", this->color());
     viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(d->nb_lines/3));
     vaos[Scene_combinatorial_map_item_priv::Edges]->release();
     d->program->release();

}
 void Scene_combinatorial_map_item::drawPoints(CGAL::Three::Viewer_interface* viewer) const
{
     if(!are_buffers_filled)
     {
         d->compute_elements();
         d->initialize_buffers(viewer);
     }
     vaos[Scene_combinatorial_map_item_priv::Points]->bind();
     d->program=getShaderProgram(PROGRAM_WITHOUT_LIGHT);
     attribBuffers(viewer,PROGRAM_WITHOUT_LIGHT);
     d->program->bind();
     d->program->setAttributeValue("colors", this->color());
     viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(d->nb_points/3));
     vaos[Scene_combinatorial_map_item_priv::Points]->release();
     d->program->release();
}

 bool Scene_combinatorial_map_item::is_from_corefinement() const{return d->address_of_A!=NULL;}
