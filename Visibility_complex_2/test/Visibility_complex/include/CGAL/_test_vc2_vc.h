#include<CGAL/Visibility_complex_2.h>
#include<CGAL/Shortest_path_2.h>
// typedef CGAL::Shortest_path_2<Gt> SP;
// typedef SP::Visibility_complex_2 VC;
typedef CGAL::Visibility_complex_2<Gt> VC;
typedef VC::Constraint_input Constraint_input;

#include<sstream>


int main() {
  std::istringstream di(input_disks);
  std::istringstream ci(input_constraints);

  std::istream_iterator<Gt::Disk> disk_it(di),disk_end;
  std::istream_iterator<Constraint_input> constraint_it(ci),constraint_end;

  VC vc(disk_it,disk_end);
  VC vcc(vc.disks_begin(),vc.disks_end(),constraint_it,constraint_end);
 
  if (vc.size()!=2*nbit) return 1;
  if (vcc.size()!=2*ncbit) return 1;
  return 0;
}
