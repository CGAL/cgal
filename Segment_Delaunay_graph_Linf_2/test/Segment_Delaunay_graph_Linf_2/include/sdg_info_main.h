#ifndef CGAL_SDG_INFO_MAIN_H
#define CGAL_SDG_INFO_MAIN_H 1


int main()
{
  bool b;
  SDG2 sdg;

  print_separator();
  b = test_info(sdg, "data/sitesx.cin");
  assert( b );

  print_separator();
  b = test_info(sdg, "data/sitesxx.cin");
  assert( b );

  print_separator();

  return 0;
}

#endif // CGAL_SDG_INFO_MAIN_H
