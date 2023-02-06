#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <fstream>

#define TEST_ALL_ORIENTATIONS 1

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> SMesh;

template<class T>
struct deque_wrapper : public std::deque<T>
{
  template<class U>
  void reserve(const U&){}
};

template<class TriangleMesh, class NamedParameters>
bool test_orientation(const TriangleMesh& tm, bool is_positive, const NamedParameters& np)
{
  typedef boost::graph_traits<TriangleMesh> Graph_traits;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  typedef typename CGAL::GetVertexPointMap<TriangleMesh, NamedParameters>::const_type Vpm;
  typedef typename CGAL::GetInitializedFaceIndexMap<TriangleMesh, NamedParameters>::const_type Fid_map;

  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  Vpm vpm = choose_parameter(get_parameter(np, CGAL::internal_np::vertex_point),
                             CGAL::get_const_property_map(boost::vertex_point, tm));

  Fid_map fid_map = CGAL::get_initialized_face_index_map(tm, np);

  std::vector<std::size_t> face_cc(num_faces(tm), std::size_t(-1));

  // set the connected component id of each face
  std::size_t nb_cc = PMP::connected_components(tm,
                                           CGAL::make_compose_property_map(fid_map,CGAL::make_property_map(face_cc)),
                                           CGAL::parameters::face_index_map(fid_map));

  // extract a vertex with max z coordinate for each connected component
  std::vector<vertex_descriptor> xtrm_vertices(nb_cc, Graph_traits::null_vertex());
  for(vertex_descriptor vd : vertices(tm))
  {
    face_descriptor test_face = face(halfedge(vd, tm), tm);
    if(test_face == Graph_traits::null_face())
      test_face = face(opposite(halfedge(vd, tm), tm), tm);
    std::size_t cc_id = face_cc[get(fid_map,test_face )];
    if (xtrm_vertices[cc_id]==Graph_traits::null_vertex())
      xtrm_vertices[cc_id]=vd;
    else
      if (get(vpm, vd).z()>get(vpm,xtrm_vertices[cc_id]).z())
        xtrm_vertices[cc_id]=vd;
  }
  std::vector<std::vector<face_descriptor> > ccs(nb_cc);
  for(face_descriptor fd : faces(tm))
  {
    ccs[face_cc[get(fid_map,fd)]].push_back(fd);
  }

  //test ccs orientation
  for(std::size_t id=0; id<nb_cc; ++id)
  {
    if((!PMP::internal::is_outward_oriented(xtrm_vertices[id], tm, np)
         &&  is_positive)
       || (PMP::internal::is_outward_oriented(xtrm_vertices[id], tm, np)
           &&  !is_positive))
    {
      std::cerr<<" the orientation failed"<<std::endl;
      return false;
    }
  }
  return true;
}


int main()
{

  std::ifstream input("data-coref/nested_cubes_invalid_volume.off");
  assert(input);
  SMesh sm1, sm2, sm3, sm4, volume, volume_copy;
  input >> sm1;
  sm2 = sm1;
  sm3 = sm1;
  sm4 = sm1;
  volume = sm1;
  volume_copy = volume;
  PMP::orient(sm1);
  if(!test_orientation(sm1, true, CGAL::parameters::default_values()))
    return 1;
  typedef boost::property_map<SMesh, CGAL::vertex_point_t>::type Ppmap;
  typedef boost::property_map<SMesh, CGAL::face_index_t>::type Fidmap;
  Ppmap vpmap2 = get(CGAL::vertex_point, sm2);
  Fidmap fidmap2 = get(CGAL::face_index, sm2);

  PMP::orient(sm2, CGAL::parameters::vertex_point_map(vpmap2)
                                    .face_index_map(fidmap2));
  if(!test_orientation(sm2, true, CGAL::parameters::vertex_point_map(vpmap2)
                                                   .face_index_map(fidmap2)))
  {
    std::cerr << "ERROR for test1\n";
    return 1;
  }

  PMP::orient(sm3, CGAL::parameters::outward_orientation(false));
  if(!test_orientation(sm3, false, CGAL::parameters::default_values()))
  {
    std::cerr << "ERROR for test2\n";
    return 1;
  }

  Ppmap vpmap4 = get(CGAL::vertex_point, sm4);
  Fidmap fidmap4 = get(CGAL::face_index, sm4);

  PMP::orient(sm4, CGAL::parameters::vertex_point_map(vpmap4)
                                    .face_index_map(fidmap4)
                                    .outward_orientation(false));
  if(!test_orientation(sm4, false, CGAL::parameters::vertex_point_map(vpmap4)
                                        .face_index_map(fidmap4)))
  {
    std::cerr << "ERROR for test3\n";
    return 1;
  }

  assert(!PMP::does_bound_a_volume(volume));
  PMP::orient_to_bound_a_volume(volume);
  if( !PMP::does_bound_a_volume(volume))
  {
    std::cerr << "ERROR for test4\n";
    return 1;
  }

  PMP::orient_to_bound_a_volume(volume);
  if( !PMP::does_bound_a_volume(volume))
  {
    std::cerr << "ERROR for test5\n";
    return 1;
  }

  PMP::orient_to_bound_a_volume(volume_copy, CGAL::parameters::outward_orientation(false));
  if( !PMP::does_bound_a_volume(volume_copy) )
  {
    std::cerr << "ERROR for test6\n";
    return 1;
  }

#ifdef TEST_ALL_ORIENTATIONS //takes around 2 hours
  std::cout<<"testing ALL orientations..."<<std::endl;
  SMesh::Property_map<SMesh::Face_index, std::size_t> fccmap =
      sm1.add_property_map<SMesh::Face_index, std::size_t>("f:CC").first;
  std::vector<bool> is_cc_o_or;
  PMP::orient_to_bound_a_volume(sm1);
  PMP::does_bound_a_volume(sm1,
                           CGAL::parameters::is_cc_outward_oriented(std::ref(is_cc_o_or)));

  std::size_t nb_ccs = PMP::connected_components(sm1, fccmap);

  std::vector< std::vector<SMesh::Face_index> > faces_per_cc(nb_ccs);
  for(SMesh::Face_index fd : faces(sm1))
  {
    std::size_t cc_id = get(fccmap, fd);
    faces_per_cc[cc_id].push_back(fd);
  }
  //double total_length = 1<<20;
  CGAL::Timer timer;
  timer.start();
  //for(std::size_t i=1; i<total_length; ++i)//0 is initial state, already tested
  int i = 1;
  {
    SMesh loop_m = sm1;
    int i_bis = i;
    int cc = 0;
    while(i_bis)
    {
      if(i_bis&1)
      {
        //reverse_orientation(cc)
        PMP::reverse_face_orientations(faces_per_cc[cc], loop_m);
      }
      //test volume
      i_bis>>=1;
      ++cc;
    }
    deque_wrapper<bool> is_loop_o_o;
    PMP::orient_to_bound_a_volume(loop_m);

    if( !PMP::does_bound_a_volume(loop_m,
                                  CGAL::parameters::is_cc_outward_oriented(std::ref(is_loop_o_o))) )
    {
      std::cerr << "ERROR for test7\n";
      return 1;
    }
    if(is_loop_o_o.empty())
      return 1;
    for(std::size_t k = 0; k< is_loop_o_o.size(); ++k){
      if(is_loop_o_o[k]!= is_cc_o_or[k])
      {
        std::cerr << "ERROR for test7\n";
        return 1;
      }
    }
    if(i%1000 == 0){
      timer.stop();
      double remaining = ((1<<20) -i)* timer.time()/1000.0 ;
      std::cout<<remaining/60.0<<"min remaining."<<std::endl;
      timer.reset();
      timer.start();
    }
  }
  std::cout<<"finished ! "<<std::endl;
#endif
  return 0;
}
