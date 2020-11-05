#include <CGAL/basic.h>
#include <CGAL/Nef_3/SNC_io_parser.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <sstream>

namespace CGAL {

template <typename Nef_>
class grid_generator : public CGAL::SNC_decorator<typename Nef_::SNC_structure> {
  typedef Nef_                                          Nef_polyhedron;
  typedef typename Nef_polyhedron::Aff_transformation_3 Aff_transformation_3;
  typedef typename Nef_polyhedron::SNC_structure        SNC_structure;
  typedef typename Nef_polyhedron::Kernel               Kernel;
  typedef typename Kernel::RT                           RT;
  typedef typename CGAL::SNC_decorator<SNC_structure>   Base;
  typedef typename SNC_structure::Sphere_map  Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>      SM_decorator;
  typedef typename SNC_structure::Infi_box    Infi_box;

  typedef typename SNC_structure::Vertex_iterator Vertex_iterator;
  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Vertex_const_handle Vertex_const_handle;
  typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::Halfedge_const_handle Halfedge_const_handle;
  typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator;
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;
  typedef typename SNC_structure::Halffacet_const_handle Halffacet_const_handle;
  typedef typename SNC_structure::Volume_iterator Volume_iterator;
  typedef typename SNC_structure::Volume_handle Volume_handle;
  typedef typename SNC_structure::Volume_const_handle Volume_const_handle;
  typedef typename SNC_structure::SVertex_iterator SVertex_iterator;
  typedef typename SNC_structure::SVertex_handle SVertex_handle;
  typedef typename SNC_structure::SVertex_const_handle SVertex_const_handle;
  typedef typename SNC_structure::SHalfedge_iterator SHalfedge_iterator;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SHalfedge_const_handle SHalfedge_const_handle;
  typedef typename SNC_structure::SFace_iterator SFace_iterator;
  typedef typename SNC_structure::SFace_handle SFace_handle;
  typedef typename SNC_structure::SFace_const_handle SFace_const_handle;
  typedef typename SNC_structure::SHalfloop_iterator SHalfloop_iterator;
  typedef typename SNC_structure::SHalfloop_handle SHalfloop_handle;
  typedef typename SNC_structure::SHalfloop_const_handle SHalfloop_const_handle;
  typedef typename SNC_structure::Object_iterator Object_iterator;
  typedef typename SNC_structure::Object_handle Object_handle;
  typedef typename SNC_structure::SFace_cycle_iterator SFace_cycle_iterator;
  typedef typename SNC_structure::Halffacet_cycle_iterator Halffacet_cycle_iterator;
  typedef typename SNC_structure::Shell_entry_iterator Shell_entry_iterator;
  typedef typename SNC_structure::SHalfedge_around_svertex_circulator
                                  SHalfedge_around_svertex_circulator;
  typedef typename SNC_structure::SHalfedge_around_sface_circulator
                                  SHalfedge_around_sface_circulator;

  typedef typename Nef_polyhedron::Point_3  Point_3;
  typedef typename Nef_polyhedron::Vector_3 Vector_3;
  typedef typename Nef_polyhedron::Plane_3  Plane_3;

  typedef typename Infi_box::Standard_kernel  Standard_kernel;
  typedef typename Standard_kernel::Point_3   Standard_point;
  typedef typename Standard_kernel::Vector_3  Standard_vector;
  typedef typename Standard_kernel::Plane_3   Standard_plane;

 private:
  std::ostream& out;
  bool verbose;
  bool reduce;
  bool sorted;
  bool equilateral;

  int vertex_offset;
  int edge_offset;
  int facet_offset;
  int volume_offset;
  int sedge_offset;
  int sloop_offset;
  int sface_offset;

  RT x_min, x_max, y_min, y_max, z_min, z_max;

  CGAL::Object_index<Vertex_iterator> VI;
  CGAL::Object_index<Halfedge_iterator> EI;
  CGAL::Object_index<Halffacet_iterator>    FI;
  CGAL::Object_index<Volume_iterator>   CI;
  CGAL::Object_index<SHalfedge_iterator> SEI;
  CGAL::Object_index<SHalfloop_iterator>   SLI;
  CGAL::Object_index<SFace_iterator>     SFI;
  std::list<Vertex_iterator> VL;
  std::list<Halfedge_iterator> EL;
  std::list<Halffacet_iterator> FL;
  std::list<Volume_iterator> CL;
  std::list<SHalfedge_iterator> SEL;
  std::list<SHalfloop_iterator> SLL;
  std::list<SFace_iterator> SFL;
  std::vector<Vertex_iterator>   Vertex_of;
  std::vector<Halfedge_iterator> Edge_of;
  std::vector<Halffacet_iterator>    Halffacet_of;
  std::vector<Volume_iterator>   Volume_of;
  std::vector<SVertex_iterator>  SVertex_of;
  std::vector<SHalfedge_iterator> SEdge_of;
  std::vector<SHalfloop_iterator> SLoop_of;
  std::vector<SFace_iterator>     SFace_of;
  long i,vn,en,fn,cn,sen,sln,sfn;

public:
  std::string add_string_and_int(std::string s, int offset) const {
    int i=atoi(s.c_str());
    i+=offset;
    std::ostringstream out;
    out << i;
    return out.str();
  }

  std::string index(Vertex_iterator v) const
  { return VI(v,verbose); }
  std::string index(Halfedge_iterator e) const
  { return EI(e,verbose); }
  std::string index(Halffacet_iterator f) const
  { return FI(f,verbose); }
  std::string index(Volume_iterator c) const
  { return CI(c,verbose); }
  std::string index(SHalfedge_iterator e) const
  { return SEI(e,verbose); }
  std::string index(SHalfloop_iterator l) const
  { return SLI(l,verbose); }
  std::string index(SFace_iterator f) const
  { return SFI(f,verbose); }


  std::string index(Vertex_iterator v, int offset) const {
    return add_string_and_int(VI(v,verbose),offset);
  }
  std::string index(Halfedge_iterator v, int offset) const {
    return add_string_and_int(EI(v,verbose),offset);
  }
  std::string index(Halffacet_iterator v, int offset) const {
    return add_string_and_int(FI(v,verbose),offset);
  }
  std::string index(Volume_iterator v, int offset) const {
    return add_string_and_int(CI(v,verbose),offset);
  }
  std::string index(SHalfedge_iterator v, int offset) const {
    return add_string_and_int(SEI(v,verbose),offset);
  }
  std::string index(SHalfloop_iterator v, int offset) const {
    return add_string_and_int(SLI(v,verbose),offset);
  }
  std::string index(SFace_iterator v, int offset) const {
    return add_string_and_int(SFI(v,verbose),offset);
  }

  bool is_zero(Volume_iterator c) const {
    std::string s(CI(c,verbose));
    int i = atoi(s.c_str());
    return (i== 0);
  }

 template <typename Iter, typename Index>
    void output_sorted_indexes(Iter begin, Iter end, Index i, int offset) const {
   std::string s(i(begin,verbose));
   int low = atoi(s.c_str());
   int high = low;
   for(Iter it=begin; it != end; it++) {
     s = i(it,verbose);
     int t = atoi(s.c_str());
     if(t < low) low = t;
     if(t > high) high = t;
   }
   out << low+offset << " " << high+offset << ", ";
 }

  void print_vertex(Vertex_handle v, const Aff_transformation_3& aff) const {
    // syntax: index { svs sve, ses see, sfs sfe, sl | point } mark
    SM_decorator SD(&*v);
    out << index(v,vertex_offset) << " { ";
    if(sorted) {
      output_sorted_indexes(v->svertices_begin(),
                            v->svertices_end(), EI, edge_offset);
      output_sorted_indexes(v->shalfedges_begin(),
                            v->shalfedges_end(), SEI, sedge_offset);
      output_sorted_indexes(v->sfaces_begin(),
                            v->sfaces_end(), SFI, sface_offset);
      out << index(SD.shalfloop(),sloop_offset) << " | ";
    }
    else {
      out
        << index(v->svertices_begin(),edge_offset) << " "
        << index(v->svertices_last(),edge_offset) << ", "
        << index(v->shalfedges_begin(),sedge_offset) << " "
        << index(v->shalfedges_last(),sedge_offset) << ", "
        << index(v->sfaces_begin(),sface_offset) << " "
        << index(v->sfaces_last(),sface_offset) << ", "
        << index(SD.shalfloop(),sloop_offset) << " | ";
    }
    Point_3 p(v->point());
    p=p.transform(aff);
    if(reduce) {
      Standard_point sp(Infi_box::standard_point(p));
      out << sp.hx() << " " << sp.hy() << " " << sp.hz() << " " << sp.hw();
    }
    else
      out << p.hx() << " " << p.hy() << " " << p.hz() << " " << p.hw();

    out << " } "  << v->mark() << std::endl;
  }

  void print_edge(Halfedge_handle e) const {
    // syntax: index { twin, source, isolated incident_object | spoint } mark
    SM_decorator D(&*e->source());
    out << index(e,edge_offset) << " { "
        << index(e->twin(),edge_offset) << ", "
        << index(e->source(),vertex_offset) << ", ";
    if ( D.is_isolated(e) ) out << "1 " << index(e->incident_sface(),sface_offset);
    else out << "0 " << index(D.first_out_edge(e),sedge_offset);
    out << " | ";
    if(reduce) {
      Standard_vector p(Infi_box::standard_vector(e->vector()));
      out << p.hx() << " " << p.hy() << " " << p.hz() << " " << p.hw();
    }
    else {
      Vector_3 p(e->vector());
      out << p.hx() << " " << p.hy() << " " << p.hz() << " " << p.hw();
    }
    out << " } "<< e->mark() << std::endl;
  }

  void print_facet(Halffacet_handle f, const Aff_transformation_3& aff) const {
    // syntax: index { twin, fclist, ivlist, volume | plane } mark
    out << index(f,facet_offset) << " { ";
    out << index(f->twin(),facet_offset) << ", ";
    Halffacet_cycle_iterator it;
    CGAL_forall_facet_cycles_of(it,f)
      if ( it.is_shalfedge() )
        out << index(SHalfedge_handle(it),sedge_offset) << ' ';
    out << ", ";
    CGAL_forall_facet_cycles_of(it,f)
      if ( it.is_shalfloop() )
        out << index(SHalfloop_handle(it),sloop_offset) << ' ';
    out << ", "
        << (is_zero(f->incident_volume())?"0":index(f->incident_volume(),volume_offset))
        << " | ";
    Plane_3 p(f->plane());
    p=p.transform(aff);
    if(reduce) {
      Standard_plane sp(Infi_box::standard_plane(p));
      out << sp.a() << " " << sp.b() << " " << sp.c() << " " << sp.d();
    }
    else {
      out << p.a() << " " << p.b() << " " << p.c() << " " << p.d();
    }
    out << " } " << f->mark() << std::endl;
  }

  void print_outer_volume(int n) const {
    Volume_handle c = *CL.begin();
    out << index(c,volume_offset) << " { ";
    Shell_entry_iterator it;
    int offset=0;
    for(int i=0; i<n; ++i) {
      CGAL_forall_shells_of(it,c) {
        out << index(SFace_handle(it),offset) << ' ';
      }
      offset+=SFL.size();
    }
    out << "} " << c->mark() << std::endl;
  }

  void print_volume(Volume_handle c) const {
    // syntax: index { shlist } mark
    out << index(c,volume_offset) << " { ";
    Shell_entry_iterator it;
    CGAL_forall_shells_of(it,c)
      if(!reduce || Infi_box::is_standard(SFace_handle(it)->center_vertex()->point())) {
        out << index(SFace_handle(it),sface_offset) << ' ';
      }
    out << "} " << c->mark() << std::endl;
  }

  void print_sedge(SHalfedge_handle e) const {
    //index { twin, sprev, snext, source, sface, prev, next, facet | circle } mark
    SM_decorator D(&*e->source()->source());
    out << index(e,sedge_offset) << " { "
        << index(e->twin(),sedge_offset) << ", "
        << index(e->sprev(),sedge_offset) << ", "
        << index(e->snext(),sedge_offset) << ", "
        << index(e->source(),edge_offset) << ", "
        << index(e->incident_sface(),sface_offset) << ", "
        << index(e->prev(),sedge_offset) << ", "
        << index(e->next(),sedge_offset) << ", "
        << index(e->facet(),facet_offset)
        << " | ";
    if(reduce) {
      Standard_plane p(Infi_box::standard_plane(e->circle()));
      out << p.a() << " " << p.b() << " " << p.c() << " " << p.d();
    }
    else {
      Plane_3 p(e->circle());
      out << p.a() << " " << p.b() << " " << p.c() << " " << p.d();
    }
    out << " } " << e->mark() << "\n";
  }

  void print_sloop(SHalfloop_handle l) const {
    // syntax: index { twin, sface, facet | circle } mark
    SM_decorator D(&*l->incident_sface()->center_vertex());
    out << index(l,sloop_offset) << " { "
        << index(l->twin(),sloop_offset) << ", "
        << index(l->incident_sface(),sface_offset) << ", "
        << index(l->facet(),facet_offset)
        << " | ";
    if(reduce) {
      Standard_plane p(Infi_box::standard_plane(l->circle()));
      out << p.a() << " " << p.b() << " " << p.c() << " " << p.d();
    }
    else {
      Plane_3 p(l->circle());
      out << p.a() << " " << p.b() << " " << p.c() << " " << p.d();
    }
    out << " } " << l->mark() << "\n";
  }

  void print_sface(SFace_handle f) const {
    // syntax: index { vertex, fclist, ivlist, sloop, volume }
    SM_decorator D(&*f->center_vertex());
    out << index(f,sface_offset) << " { "
        << index(f->center_vertex(),vertex_offset) << ", ";
    SFace_cycle_iterator it;
    CGAL_forall_sface_cycles_of(it,f)
      if ( it.is_shalfedge() )
        out << index(SHalfedge_handle(it),sedge_offset) << ' ';
    out << ", ";
    CGAL_forall_sface_cycles_of(it,f)
      if ( it.is_svertex() )
        out << index(SVertex_handle(it),edge_offset) << ' ';
    out << ", ";
    CGAL_forall_sface_cycles_of(it,f)
      if ( it.is_shalfloop() )
        out << index(SHalfloop_handle(it),sloop_offset);
    out << ", "
        << (is_zero(f->volume())?"0":index(f->volume(),volume_offset))
        << " } "
        << f->mark() <<"\n";
  }

  void print_items(const Aff_transformation_3& aff, bool first) {

    typename std::list<Vertex_iterator>::const_iterator v;
    for(v=VL.begin();v!=VL.end();v++)
      print_vertex(*v, aff);
    vertex_offset+=VL.size();

    typename std::list<Halfedge_iterator>::const_iterator e;
    for(e=EL.begin();e!=EL.end();e++)
      print_edge(*e);
    edge_offset+=EL.size();

    typename std::list<Halffacet_iterator>::const_iterator f;
    for(f=FL.begin();f!=FL.end();f++)
      print_facet(*f,aff);
    facet_offset+=FL.size();

    typename std::list<Volume_iterator>::const_iterator c;
    c=CL.begin();
    if(first) c++;
    for(;c!=CL.end();c++)
      print_volume(*c);
    volume_offset+=CL.size();
    if(!first) volume_offset--;

    typename std::list<SHalfedge_iterator>::const_iterator se;
    for(se=SEL.begin();se!=SEL.end();se++)
      print_sedge(*se);
    sedge_offset+=SEL.size();

    typename std::list<SHalfloop_iterator>::const_iterator sl;
    for(sl=SLL.begin();sl!=SLL.end();sl++)
      print_sloop(*sl);
    sloop_offset+=SLL.size();

    typename std::list<SFace_iterator>::const_iterator sf;
    for(sf=SFL.begin();sf!=SFL.end();sf++)
      print_sface(*sf);
    sface_offset+=SFL.size();
  }

  Aff_transformation_3 compute_translation() {
    Vertex_iterator vi(this->vertices_begin());
    x_min = vi->point().hx()/vi->point().hw();
    x_max = vi->point().hx()/vi->point().hw();
    y_min = vi->point().hy()/vi->point().hw();
    y_max = vi->point().hy()/vi->point().hw();
    z_min = vi->point().hz()/vi->point().hw();
    z_max = vi->point().hz()/vi->point().hw();
    for(;vi != this->vertices_end();++vi) {
      if(vi->point().hx()/vi->point().hw() < x_min)
        x_min = vi->point().hx()/vi->point().hw();
      if(vi->point().hx()/vi->point().hw() > x_max)
        x_max = vi->point().hx()/vi->point().hw();
      if(vi->point().hy()/vi->point().hw() < y_min)
        y_min = vi->point().hy()/vi->point().hw();
      if(vi->point().hy()/vi->point().hw() > y_max)
        y_max = vi->point().hy()/vi->point().hw();
      if(vi->point().hz()/vi->point().hw() < z_min)
        z_min = vi->point().hz()/vi->point().hw();
      if(vi->point().hz()/vi->point().hw() > z_max)
        z_max = vi->point().hz()/vi->point().hw();
    }

    x_min-= (x_min<=0?1:0);
    y_min-= (y_min<=0?1:0);
    z_min-= (z_min<=0?1:0);
    x_max+=1;
    y_max+=1;
    z_max+=1;

    if(equilateral) {
      RT max = x_max-x_min;
      if(y_max-y_min > max) max=y_max-y_min;
      if(z_max-z_min > max) max=z_max-z_min;
      return Aff_transformation_3(max, 0, 0, -x_min,
                                  0, max, 0, -y_min,
                                  0, 0, max, -z_min,
                                  1);
    }

    return Aff_transformation_3(x_max-x_min, 0, 0, -x_min,
                                0, y_max-y_min, 0, -y_min,
                                0, 0, z_max-z_min, -z_min,
                                1);
  }

 public:
  void print(int nx, int ny, int nz) {

    out << "Selective Nef Complex" << std::endl;
    if(this->is_extended_kernel() && (!reduce || !this->is_bounded()))
      out << "extended" << std::endl;
    else
      out << "standard" << std::endl;
    out << "vertices   " << VL.size()*nx*ny*nz << std::endl;
    out << "halfedges  " << EL.size()*nx*ny*nz << std::endl;
    out << "facets     " << FL.size()*nx*ny*nz << std::endl;
    out << "volumes    " << (CL.size()-1)*nx*ny*nz+1 << std::endl;
    out << "shalfedges " << SEL.size()*nx*ny*nz << std::endl;
    out << "shalfloops " << SLL.size()*nx*ny*nz << std::endl;
    out << "sfaces     " << SFL.size()*nx*ny*nz << std::endl;

    Aff_transformation_3 trans = compute_translation();

    for(int dx=0; dx < nx; ++dx)
      for(int dy=0; dy < ny; ++dy)
        for(int dz=0; dz < nz; ++dz) {
          Point_3 vp(dx,dy,dz);
          vp = vp.transform(trans);
          Aff_transformation_3 aff(CGAL::TRANSLATION,vp-CGAL::ORIGIN);
          typename std::list<Vertex_iterator>::const_iterator v;
          for(v=VL.begin();v!=VL.end();v++)
            print_vertex(*v, aff);
          vertex_offset+=VL.size();
          edge_offset+=EL.size();
          sedge_offset+=SEL.size();
          sloop_offset+=SLL.size();
          sface_offset+=SFL.size();
        }

    vertex_offset=0;
    edge_offset=0;
    sedge_offset=0;
    sloop_offset=0;
    sface_offset=0;

    for(int dx=0; dx < nx; ++dx)
      for(int dy=0; dy < ny; ++dy)
        for(int dz=0; dz < nz; ++dz) {
          typename std::list<Halfedge_iterator>::const_iterator e;
          for(e=EL.begin();e!=EL.end();e++)
            print_edge(*e);
          vertex_offset+=VL.size();
          edge_offset+=EL.size();
          sedge_offset+=SEL.size();
          sface_offset+=SFL.size();
        }

    vertex_offset=0;
    edge_offset=0;
    sedge_offset=0;
    sface_offset=0;

    bool first=true;
    for(int dx=0; dx < nx; ++dx)
      for(int dy=0; dy < ny; ++dy)
        for(int dz=0; dz < nz; ++dz) {
          Point_3 vp(dx,dy,dz);
          vp = vp.transform(trans);
          Aff_transformation_3 aff(CGAL::TRANSLATION,vp-CGAL::ORIGIN);
          typename std::list<Halffacet_iterator>::const_iterator f;
          for(f=FL.begin();f!=FL.end();f++)
            print_facet(*f,aff);
          facet_offset+=FL.size();
          sedge_offset+=SEL.size();
          sloop_offset+=SLL.size();
          volume_offset+=CL.size()-1;
          first = false;
        }

    facet_offset=0;
    sedge_offset=0;
    sloop_offset=0;
    volume_offset=0;

    first=true;
    for(int dx=0; dx < nx; ++dx)
      for(int dy=0; dy < ny; ++dy)
        for(int dz=0; dz < nz; ++dz) {
          typename std::list<Volume_iterator>::const_iterator c;
          c=CL.begin();
          if(first)
            print_outer_volume(nx*ny*nz);
          c++;
          for(;c!=CL.end();c++)
            print_volume(*c);
          volume_offset+=CL.size()-1;
          sface_offset+=SFL.size();
          first = false;
        }

    sface_offset=0;
    volume_offset=0;

    for(int dx=0; dx < nx; ++dx)
      for(int dy=0; dy < ny; ++dy)
        for(int dz=0; dz < nz; ++dz) {
          typename std::list<SHalfedge_iterator>::const_iterator se;
          for(se=SEL.begin();se!=SEL.end();se++)
            print_sedge(*se);
          sedge_offset+=SEL.size();
          sface_offset+=SFL.size();
          facet_offset+=FL.size();
          edge_offset+=EL.size();
        }

    facet_offset=0;
    sedge_offset=0;
    edge_offset=0;
    sface_offset=0;

    for(int dx=0; dx < nx; ++dx)
      for(int dy=0; dy < ny; ++dy)
        for(int dz=0; dz < nz; ++dz) {
          typename std::list<SHalfloop_iterator>::const_iterator sl;
          for(sl=SLL.begin();sl!=SLL.end();sl++)
            print_sloop(*sl);
          sloop_offset+=SLL.size();
          sface_offset+=SFL.size();
          facet_offset+=FL.size();
        }

    facet_offset=0;
    sloop_offset=0;
    sface_offset=0;

    for(int dx=0; dx < nx; ++dx)
      for(int dy=0; dy < ny; ++dy)
        for(int dz=0; dz < nz; ++dz) {
          typename std::list<SFace_iterator>::const_iterator sf;
          for(sf=SFL.begin();sf!=SFL.end();sf++)
            print_sface(*sf);
          sface_offset+=SFL.size();
          sedge_offset+=SEL.size();
          sloop_offset+=SLL.size();
          vertex_offset+=VL.size();
          volume_offset+=CL.size()-1;
        }

    out << "end Selective Nef complex" << std::endl;
  }

  RT size_x() { return x_max-x_min; }
  RT size_y() { return y_max-y_min; }
  RT size_z() { return z_max-z_min; }

  grid_generator(std::ostream& os, const Nef_polyhedron& N, bool equil = true) :
    Base(*const_cast<SNC_structure*>(N.sncp())), out(os), sorted(true), equilateral(equil),
    FI(halffacets_begin(),halffacets_end(),'F'),
    CI(volumes_begin(),volumes_end(),'C'),
    SEI(shalfedges_begin(),shalfedges_end(),'e'),
    SLI(shalfloops_begin(),shalfloops_end(),'l'),
    SFI(sfaces_begin(),sfaces_end(),'f'),
    vn(N.sncp()->number_of_vertices()),
    en(N.sncp()->number_of_halfedges()),
    fn(N.sncp()->number_of_halffacets()),
    cn(N.sncp()->number_of_volumes()),
    sen(N.sncp()->number_of_shalfedges()),
    sln(N.sncp()->number_of_shalfloops()),
    sfn(N.sncp()->number_of_sfaces()) {

    vertex_offset=0;
    edge_offset=0;
    facet_offset=0;
    volume_offset=0;
    sedge_offset=0;
    sloop_offset=0;
    sface_offset=0;

    verbose = (out.iword(CGAL::IO::mode) != CGAL::IO::ASCII &&
               out.iword(CGAL::IO::mode) != CGAL::IO::BINARY);
    reduce = this->is_extended_kernel() && sorted && this->is_bounded();

    Vertex_iterator vi;
    CGAL_forall_vertices(vi, *this->sncp()) {
      VL.push_back(vi);
      if(sorted) {
        vi->point() = normalized(vi->point());
        if(vi->has_shalfloop() &&
           sort_sloops<SNC_structure>(*this->sncp())(vi->shalfloop()->twin(),
                                                     vi->shalfloop()))
          vi->shalfloop() = vi->shalfloop()->twin();
      }
    }
    if(sorted) {
      VL.sort(sort_vertices<SNC_structure>(*this->sncp()));
    }
    if(reduce)
      for(int k=0; k<4; k++){
        VL.pop_front(); VL.pop_back();
      }
    int i = 0;
    typename std::list<Vertex_iterator>::iterator vl;
    for(vl = VL.begin(); vl != VL.end(); vl++)
      VI[*vl] = i++;

    SM_decorator SD;
    Halfedge_iterator ei;
    CGAL_forall_halfedges(ei, *this->sncp()) {
      EL.push_back(ei);
      if(sorted) {
        //      std::cerr << point(ei) << " | " << normalized(point(ei)) << " |";
        ei->point() = normalized(ei->point());
        //      std::cerr << point(ei) << std::endl;
        sort_sedges<SNC_structure> sortSE(*this->sncp());
        SHalfedge_handle new_outedge = ei->out_sedge();
        SHalfedge_around_svertex_circulator cb(new_outedge), ce(cb);
        CGAL_For_all(cb,ce) {
          if(cb != new_outedge && sortSE(cb,new_outedge))
            new_outedge = cb;
        }
      ei->out_sedge() = new_outedge;
      }
    }
    if(sorted) EL.sort(sort_edges<SNC_structure>(*this->sncp()));
    if(reduce)
      for(int k=0; k<12; k++){
        EL.pop_front(); EL.pop_back();
      }
    i = 0;
    typename std::list<Halfedge_iterator>::iterator el;
    for(el = EL.begin(); el != EL.end(); el++)
      EI[*el] = i++;

    Halffacet_iterator fi;
    CGAL_forall_halffacets(fi, *this->sncp()){
      if(sorted) {
        fi->plane() = normalized(fi->plane());
        fi->boundary_entry_objects().sort(sort_facet_cycle_entries<Base>((Base) *this));
      }
      FL.push_back(fi);
    }
    if(sorted) FL.sort(sort_facets<SNC_structure>(*this->sncp()));
    if(reduce) {
      for(int k=0; k<6; k++){
        FL.pop_front();
        FL.pop_back();
      }
    }
    i = 0;
    typename std::list<Halffacet_iterator>::iterator fl;
    for(fl = FL.begin(); fl != FL.end(); fl++)
      FI[*fl] = i++;

    SHalfedge_iterator sei;
    CGAL_forall_shalfedges(sei, *this->sncp()) {
      SEL.push_back(sei);
      if(sorted)
        sei->circle() = normalized(sei->circle());
    }
    if(sorted) SEL.sort(sort_sedges<SNC_structure>(*this->sncp()));
    if(reduce)
      for(int k=0; k<24; k++){
        SEL.pop_front(); SEL.pop_back();
      }
    i = 0;
    typename std::list<SHalfedge_iterator>::iterator sel;
    for(sel = SEL.begin(); sel != SEL.end(); sel++)
      SEI[*sel] = i++;

    SHalfloop_iterator sli;
    CGAL_forall_shalfloops(sli, *this->sncp()) {
      SLL.push_back(sli);
      if(sorted)
        sli->circle() = normalized(sli->circle());
    }
    if(sorted) SLL.sort(sort_sloops<SNC_structure>(*this->sncp()));
    i = 0;
    typename std::list<SHalfloop_iterator>::iterator sll;
    for(sll = SLL.begin(); sll != SLL.end(); sll++)
      SLI[*sll] = i++;

    SFace_iterator sfi;
    CGAL_forall_sfaces(sfi, *this->sncp()) {
      if(sorted) {
        SFace_cycle_iterator fc;
        CGAL_forall_sface_cycles_of(fc, sfi) {
          if(fc.is_shalfedge()) {
            SHalfedge_handle se(fc);
            SHalfedge_around_sface_circulator cb(se), ce(cb);
            CGAL_For_all(cb,ce) {
              if(cb->source() != se->source()) {
                if(lexicographically_xyz_smaller(cb->source()->twin()->source()->point(),
                                                 se->source()->twin()->source()->point()))
                  se = cb;
              }
              else
                if(lexicographically_xyz_smaller(cb->source()->twin()->source()->point(),
                                                 se->source()->twin()->source()->point()))
                  se = cb;
            }
            *fc = se;
          }
        }
        sfi->boundary_entry_objects().sort(sort_sface_cycle_entries<Base>((Base) *this));
      }
      SFL.push_back(sfi);
    }
    if(sorted) SFL.sort(sort_sfaces<SNC_structure>(*this->sncp()));
    if(reduce)
      for(int k=0; k<8; k++){
        SFL.pop_front(); SFL.pop_back();
      }
    i = 0;
    typename std::list<SFace_iterator>::iterator sfl;
    for(sfl = SFL.begin(); sfl != SFL.end(); sfl++)
      SFI[*sfl] = i++;

    Volume_iterator ci;
    CGAL::Unique_hash_map<SFace_handle,bool> Done(false);
    find_minimal_sface_of_shell<SNC_structure> findMinSF(*this->sncp(),Done);
    CGAL_forall_volumes(ci, *this->sncp()) {
      if(sorted) {
        Shell_entry_iterator it;
        CGAL_forall_shells_of(it,ci) {
          findMinSF.minimal_sface() = SFace_handle(it);
          visit_shell_objects(SFace_handle(it),findMinSF);
          *it = findMinSF.minimal_sface();
        }
        ci->shell_entry_objects().sort(sort_shell_entries<Base>((Base)*this));
      }
      CL.push_back(ci);
    }

    if(sorted) CL.sort(sort_volumes<SNC_structure>(*this->sncp()));
    if(reduce)
      CL.pop_front();
    i = 0;
    typename std::list<Volume_iterator>::iterator cl;
    for(cl = CL.begin(); cl != CL.end(); cl++)
      CI[*cl] = i++;

    VI[vertices_end()]=-2;
    EI[halfedges_end()]=-2;
    FI[halffacets_end()]=-2;
    CI[volumes_end()]=-2;
    SEI[shalfedges_end()]=-2;
    SLI[shalfloops_end()]=-2;
    SFI[sfaces_end()]=-2;
  }
};

} //namespace CGAL
