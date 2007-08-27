#include<CGAL/Shortest_path_2.h>
typedef CGAL::Shortest_path_2<Gt> SP;
typedef SP::Visibility_complex_2 VC;
typedef VC::Constraint_input Constraint_input;

#include<sstream>


int main() {
  std::istringstream di(input_disks);
  std::istringstream ci(input_constraints);

  std::istream_iterator<Gt::Disk> disk_it(di),disk_end;
  std::istream_iterator<Constraint_input> constraint_it(ci),constraint_end;

  VC vc(disk_it,disk_end,constraint_it,constraint_end);
  SP sp(vc);
  std::vector<Gt::Bitangent_2> vb;
  sp.compute_shortest_paths(vc.disks_begin()[0]);
  sp.get_path_bitangents(vc.disks_begin()[1],std::back_inserter(vb));

  std::istringstream pi(path);
  std::istream_iterator<Constraint_input> path_it(pi),path_end;
  std::vector<Constraint_input> pv CGAL_make_vector(path_it,path_end);
  if (vb.size()!=pv.size()) return 1;
  for (std::vector<Gt::Bitangent_2>::size_type i=0;i<vb.size();++i) {
    Gt::Bitangent_2& a=vb.begin()[i];
    Gt::Bitangent_2 b(pv.begin()[i].type(),
                      &(vc.disks_begin()[pv.begin()[i].source()]),
                      &(vc.disks_begin()[pv.begin()[i].target()]));
    if (i==0) {
      if (a.target_object()!=b.target_object()||
          a.is_xx_left()!=a.is_xx_left()) return 1;
    } else if (i==vb.size()-1) {
      if (a.source_object()!=b.source_object()||
          a.is_left_xx()!=a.is_left_xx()) return 1;      
    } else if (a!=b) return 1;
  }
  return 0;
}
