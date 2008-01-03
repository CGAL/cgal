#include<CGAL/Visibility_complex_2.h>
typedef CGAL::Compute_free_bitangents_2<Gt> VC;
typedef VC::Constraint_input Constraint_input;

#include<sstream>


int main() {
  std::istringstream di(input_disks);
  std::istringstream ci(input_constraints);

  std::istream_iterator<Gt::Disk> disk_it(di),disk_end;
  std::istream_iterator<Constraint_input> constraint_it(ci),constraint_end;
  std::vector<Gt::Disk> disks (disk_it,disk_end);
  std::vector<VC::Bitangent_2> result;
  std::vector<VC::Bitangent_2> resultc;

  VC()(disks.begin(),disks.end(),std::back_inserter(result));
  if (result.size()!=
      static_cast<std::vector<VC::Bitangent_2>::size_type>(nbit)) {
    std::cout<<"Wrong number of bitangent for unconstrained scene:"
             <<result.size()<<"\n"<<std::flush;
    return 1;
  }
  VC()(disks.begin(),disks.end(),constraint_it,constraint_end,
       std::back_inserter(resultc));
 
  if (resultc.size()!=
      static_cast<std::vector<VC::Bitangent_2>::size_type>(ncbit)) {
    std::cout<<"Wrong number of bitangent for constrained scene:"
             <<resultc.size()<<"\n"<<std::flush;
    return 1;
  }
  return 0;
}
