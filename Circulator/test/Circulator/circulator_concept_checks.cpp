#include <vector>
#include <list>

#include <boost/concept/assert.hpp>

#include <CGAL/circulator.h>
#include <CGAL/Circulator/Circulator_concepts.h>

int main()
{
  // test Circulator_from_container
  typedef CGAL::Circulator_from_container< std::vector<int> > Circulator_from_vec;
  typedef CGAL::Circulator_from_container< std::list<int> > Circulator_from_list;

  BOOST_CONCEPT_ASSERT((CGAL::Concepts::RandomAccessCirculator<Circulator_from_vec>)) CGAL_UNUSED;
  BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<Circulator_from_list>)) CGAL_UNUSED;

  typedef CGAL::Circulator_from_iterator<int*> Circulator_from_intp;
  typedef CGAL::Circulator_from_iterator<std::vector<int>::iterator> Circulator_from_veci;
  typedef CGAL::Circulator_from_iterator<std::list<int>::iterator> Circulator_from_listi;
  
  BOOST_CONCEPT_ASSERT((CGAL::Concepts::RandomAccessCirculator<Circulator_from_intp>)) CGAL_UNUSED;
  BOOST_CONCEPT_ASSERT((CGAL::Concepts::RandomAccessCirculator<Circulator_from_veci>)) CGAL_UNUSED;
  BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<Circulator_from_listi>)) CGAL_UNUSED;

  return 0;
}
