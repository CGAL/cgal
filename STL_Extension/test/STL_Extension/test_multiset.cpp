// ============================================================================
//
// Copyright (c) 2005 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : src/test_multiset.C
// package       : $CGAL_Package: STL_Extension $
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
//
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Ron Wein <wein@post.tau.ac.il>
// maintainer    : Ron Wein <wein@post.tau.ac.il>
// coordinator   : ETH Zurich
//
// A test for the CGAL::Multiset container.
// ============================================================================

#include <CGAL/Multiset.h>
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cassert>

void test_massive_insert_and_erase();
void test_iterators();
void test_lookup();
void test_comparison();
void test_catenate();
void test_split();
void test_swap_and_replace();

int main ()
{
  std::cout << std::endl
            << "Testing insert and erase:" << std::endl
            << "=========================" << std::endl << std::endl;
  test_massive_insert_and_erase();

  std::cout << std::endl
            << "Testing the iterators:" << std::endl
            << "======================" << std::endl << std::endl;
  test_iterators();

  std::cout << std::endl
            << "Testing the lookup methods:" << std::endl
            << "===========================" << std::endl << std::endl;
  test_lookup();

  std::cout << std::endl
            << "Testing the comparison functions:" << std::endl
            << "=================================" << std::endl << std::endl;
  test_comparison();

  std::cout << std::endl
            << "Testing catenate:" << std::endl
            << "=================" << std::endl << std::endl;
  test_catenate();

  std::cout << std::endl
            << "Testing split:" << std::endl
            << "==============" << std::endl << std::endl;
  test_split();

  std::cout << std::endl
            << "Testing swap and replace:" << std::endl
            << "=========================" << std::endl << std::endl;
  test_swap_and_replace();

  return (0);
}

// ---------------------------------------------------------------------------
void test_massive_insert_and_erase ()
{
  typedef CGAL::Multiset<int>                Set;
  typedef Set::iterator                      Set_iter;

  Set         set;
  Set_iter    iter;
  const int   n = 2000;
  int         k;
  int         val;

  std::cout << "Inserting " << n << " numbers into the set... " << std::flush;
  for (k = 0; k < n; k++)
  {
    val = static_cast<int> (n * static_cast<double>(std::rand()) / RAND_MAX);
    set.insert (val);
    assert (set.is_valid());
  }
  std::cout << "Done." << std::endl;

  // Print the set characteristics.
  std::cout << set.size() << " numbers in the set "
            << "(tree height is " << set.height()
            << ", black-height is " << set.black_height() << ")."
            << std::endl;

  // Delete some numbers.
  size_t         n_del = 0;

  std::cout << "Deleting numbers from the set... " << std::flush;
  for (k = 0; k < 3*n; k++)
  {
    val = static_cast<int> (n * static_cast<double>(std::rand()) / RAND_MAX);
    n_del += set.erase (val);
    assert (set.is_valid());
  }
  std::cout << "Done." << std::endl;
  std::cout << n_del << " numbers have been deleted." << std::endl;

  // Print the set.
  std::cout << set.size() << " numbers in the set "
            << "(tree height is " << set.height()
            << ", black-height is " << set.black_height()
            << ")." << std::endl;

  for (iter = set.begin(); iter != set.end(); ++iter)
    std::cout << *iter << ' ';
  std::cout << std::endl;

  return;
}

// ---------------------------------------------------------------------------
void test_iterators ()
{
  typedef CGAL::Multiset<int>                Set;
  typedef Set::iterator                      Set_iter;
  typedef Set::reverse_iterator              Set_rev_iter;

  Set         set;
  const int   range = 101;           // A prime number.
  int         k;
  int         val;

  std::cout << "Inserting numbers:" << std::endl;
  for (k = 0; k < range; k++)
  {
    val = (k * (range /2)) % range;
    std::cout << val << ' ';
    set.insert (val);
  }
  std::cout << std::endl << std::endl;

  // Print the set.
  Set_iter   iter;

  std::cout << set.size() << " numbers in the set (tree height is "
            << set.height() << "):" << std::endl;
  for (iter = set.begin(); iter != set.end(); ++iter)
    std::cout << *iter << ' ';
  std::cout << std::endl;

  // Insert a few negative numbers.
  Set_iter   pos;

  pos = set.insert (-20);
  set.insert_before (pos, -30);
  set.insert_after (pos, -10);

  // Erase a few numbers.
  set.erase (pos);
  set.erase (0);
  set.erase (17);
  set.erase (19);

  // Replace 18 by 19. This does not violate the set properties:
  pos = set.find (18);
  set.replace (pos, 19);

  // Insert back the number 17, but using an inexact hint.
  pos++;
  set.insert (pos, 17);

  // Print the set in a reversed order.
  Set_rev_iter   rev_iter;

  std::cout << set.size() << " numbers in the set (tree height is "
            << set.height() << "):" << std::endl;
  for (rev_iter = set.rbegin(); rev_iter != set.rend(); ++rev_iter)
    std::cout << *rev_iter << ' ' << std::flush;
  std::cout << std::endl;

  return;
}

// ---------------------------------------------------------------------------

struct Word_index
{
  std::string     word;
  int             index;
};

struct Word_index_compare
{
  // Compare two Word_index objects.
  CGAL::Comparison_result operator() (const Word_index& w1,
                                      const Word_index& w2) const
  {
    int    res = std::strcmp (w1.word.c_str(), w2.word.c_str());

    if (res == 0)
      return (CGAL::EQUAL);
    else if (res < 0)
      return (CGAL::SMALLER);
    else
      return (CGAL::LARGER);
  }
};

struct String_word_index_compare
{
  // Compare a string to a Word_index object.
  CGAL::Comparison_result operator() (const std::string& str,
                                      const Word_index& w2) const
  {
    int    res = std::strcmp (str.c_str(), w2.word.c_str());

    if (res == 0)
      return (CGAL::EQUAL);
    else if (res < 0)
      return (CGAL::SMALLER);
    else
      return (CGAL::LARGER);
  }
};


void test_lookup ()
{
  typedef CGAL::Multiset<Word_index, Word_index_compare>   Set;
  typedef Set::const_iterator                              Set_iter;

  const char *words[] =
    {"There's", "a", "lady", "who's", "sure",
     "all", "that", "glitters", "is", "gold",
     "And", "she's", "buying", "a", "stairway", "to", "heaven",
     "And", "when", "she", "gets", "there", "she", "knows",
     "if", "the", "stores", "are", "closed",
     "With", "a", "word", "she", "can", "get", "what", "she", "came", "for",
     "Woe", "oh", "oh", "oh", "oh", "oh",
     "And", "she's", "buying", "a", "stairway", "to", "heaven"};
  const int   n_words = sizeof(words) / sizeof(char*);

  // Insert all indexed words into the set.
  Set         set;
  Word_index  obj;
  int         k;

  for (k = 0; k < n_words; k++)
  {
    obj.word = words[k];
    obj.index = k;
    set.insert (obj);
  }

  // Print the indices of the words that occur more than once.
  CGAL::Multiset<std::string>   used_words;
  const Set   cset = set;
  Set_iter    lower, upper;
  std::pair<Set_iter, Set_iter> range;
  std::pair<Set_iter, bool>     res;

  for (k = 0; k < n_words; k++)
  {
    if (cset.count (words[k], String_word_index_compare()) > 1 &&
        used_words.find (words[k]) == used_words.end())
    {
      used_words.insert (words[k]);

      std::cout << words[k] << ':';
      lower = cset.lower_bound (words[k], String_word_index_compare());
      upper = cset.upper_bound (words[k], String_word_index_compare());
      range = cset.equal_range (words[k], String_word_index_compare());

      res = cset.find_lower (words[k], String_word_index_compare());
      assert ((res.second && lower != upper) ||
                      (! res.second && lower == upper));

      assert (lower == range.first && upper == range.second);

      while (lower != upper)
      {
        std::cout << ' ' << lower->index;
        lower++;
      }
      std::cout << std::endl;
    }
  }

  return;
}

// ---------------------------------------------------------------------------

void test_comparison ()
{
  typedef CGAL::Multiset<int>                Set;

  // Construct a random set.
  const int   n = 10;
  Set         set1;
  int         k;

  for (k = 1; k <= n; k++)
    set1.insert (10*k);

  Set          set2 = set1;

  assert (set1 == set2);

  // Add elements, then compare.
  set2.insert (200);

  std::cout << "After first insertion: ";
  if (set1 == set2)
    std::cout << "S1 == S2" << std::endl;
  else if (set1 < set2)
    std::cout << "S1 < S2" << std::endl;
  else
    std::cout << "S1 > S2" << std::endl;

  set1.insert (25);

  std::cout << "After second insertion: ";
  if (set1 == set2)
    std::cout << "S1 == S2" << std::endl;
  else if (set1 < set2)
    std::cout << "S1 < S2" << std::endl;
  else
    std::cout << "S1 > S2" << std::endl;

  // Swap the two sets.
  set1.swap (set2);

  assert (set1.is_valid());
  assert (set2.is_valid());

  std::cout << "After swapping the sets: ";
  if (set1 == set2)
    std::cout << "S1 == S2" << std::endl;
  else if (set1 < set2)
    std::cout << "S1 < S2" << std::endl;
  else
    std::cout << "S1 > S2" << std::endl;

  return;
}

// ---------------------------------------------------------------------------

void test_catenate ()
{
  typedef CGAL::Multiset<int>                Set;
  typedef Set::iterator                      Set_iter;

  // Test catenation of small sets.
  Set         s1, s1_rev;
  Set         s2, s2_rev;
  const int   max_size = 4;
  int         size1, size2;
  int         k;

  for (size1 = 1; size1 <= max_size; size1++)
  {
    for (size2 = 1; size2 <= max_size; size2++)
    {
      s1.clear();
      s1_rev.clear();
      for (k = 1; k <= size1; k++)
      {
        s1.insert (k);
        s1_rev.insert (size1 - k + 1);
      }

      s2.clear();
      s2_rev.clear();
      for (k = 1; k <= size2; k++)
      {
        s2.insert (size1 + k);
        s2_rev.insert (size1 + size2 - k + 1);
      }

      s1.catenate (s2);
      assert (s1.is_valid());
      s1_rev.catenate (s2_rev);
      assert (s1_rev.is_valid());

      s1.clear();
      s1_rev.clear();
      for (k = 1; k <= size1; k++)
      {
        s1.insert (k);
        s1_rev.insert (size1 - k + 1);
      }

      s2.clear();
      s2_rev.clear();
      for (k = 1; k <= size2; k++)
      {
        s2.insert (size1 + k);
        s2_rev.insert (size1 + size2 - k + 1);
      }

      s1.catenate (s2_rev);
      assert (s1.is_valid());
      s1_rev.catenate (s2);
      assert (s1_rev.is_valid());
    }
  }

  // Construct two random sets.
  const int   n1 = 1000;
  Set         set1;
  const int   n2 = 100;
  Set         set2;
  int         val;

  for (k = 0; k < n1; k++)
  {
    val = static_cast<int> (1000 * static_cast<double>(std::rand()) /
                                   static_cast<double>(RAND_MAX));
    set1.insert (val);
  }
  assert (set1.is_valid());

  std::cout << set1.size() << " numbers in the first set "
            << "(tree height is " << set1.height()
            << ", black-height is " << set1.black_height()
            << ")." << std::endl;

  for (k = 0; k < n2; k++)
  {
    val = static_cast<int> (2000 + 1000 * static_cast<double>(std::rand()) /
                                          static_cast<double>(RAND_MAX));
    set2.insert (val);
  }
  assert (set2.is_valid());

  std::cout << set2.size() << " numbers in the second set "
            << "(tree height is " << set2.height()
            << ", black-height is " << set2.black_height()
            << ")." << std::endl;

  // Merge the sets.
  set1.catenate (set2);

  // Print the resulting set.
  Set_iter   iter;

  std::cout << set1.size() << " numbers in the merged set "
            << "(tree height is " << set1.height()
            << ", black-height is " << set1.black_height()
            << ")." << std::endl;
  for (iter = set1.begin(); iter != set1.end(); ++iter)
    std::cout << *iter << ' ';
  std::cout << std::endl;

  if (set2.size() == 0)
    std::cout << "The other set is now empty." << std::endl;

  return;
}

// ---------------------------------------------------------------------------

void test_split ()
{
  typedef CGAL::Multiset<int>                Set;
  typedef Set::iterator                      Set_iter;

  Set         set;
  const int   range = 101;           // A prime number.
  int         k;
  int         val;

  std::cout << "Inserting numbers:" << std::endl;
  for (k = 0; k < range; k++)
  {
    val = (k * (range /2)) % range;
    set.insert (val);
  }
  assert (set.is_valid());

  // Print the set.
  Set_iter    iter;

  std::cout << set.size() << " numbers in the set "
            << "(tree height is " << set.height()
            << ", black-height is " << set.black_height()
            << ")." << std::endl;
  for (iter = set.begin(); iter != set.end(); ++iter)
    std::cout << *iter << ' ';
  std::cout << std::endl;

  // Split the set at all possible positions.
  for (k = 0; k < range; k++)
  {
    Set      set1 = set;
    Set      set2;

    std::cout << "    Splitting at " << k << " ... " << std::flush;
    set1.split (k, set2);
    assert (set1.is_valid());
    assert (set2.is_valid());

    assert (set.size() == set1.size() + set2.size());
    std::cout << "OK." << std::endl;
  }

  return;
}

// ---------------------------------------------------------------------------

struct Compare_floor
{
  CGAL::Comparison_result operator() (const double& x1,
                                      const double& x2) const
  {
    const int  ix1 = static_cast<int> (x1);
    const int  ix2 = static_cast<int> (x2);

    if (ix1 < ix2)
      return (CGAL::SMALLER);
    else if (ix1 == ix2)
      return (CGAL::EQUAL);
    else
      return (CGAL::LARGER);
  }
};


void test_swap_and_replace ()
{
  typedef CGAL::Multiset<double, Compare_floor>  Set;
  typedef Set::iterator                          Set_iter;
  typedef std::list<Set_iter>                    Iter_list;

  Set                     set;
  Set_iter                iter;
  const int               M = 10;
  std::vector<Iter_list>  handles (M);
  const int               n = 30;
  double                  val;
  int                     k;

  std::cout << "Inserting " << n << " numbers into the set... " << std::flush;
  for (k = 0; k < n; k++)
  {
    val = M * (static_cast<double>(std::rand()) / RAND_MAX);
    handles[static_cast<int>(val)].push_back (set.insert (val));
  }
  assert (set.is_valid());
  std::cout << "Done." << std::endl;

  // Print the set.
  std::cout << set.size() << " numbers in the set "
            << "(tree height is " << set.height()
            << ", black-height is " << set.black_height() << ")."
            << std::endl;

  for (iter = set.begin(); iter != set.end(); ++iter)
    std::cout << *iter << ' ';
  std::cout << std::endl;

  // Go over the handle buckets and swap elements.
  std::cout << "Swapping elements in the set... " << std::flush;
  for (k = 0; k < M; k ++)
  {
    if (handles[k].size() >= 2)
    {
      set.swap (handles[k].front(), handles[k].back());
      assert (set.is_valid());
    }
    else if (handles[k].size() == 1)
    {
      set.replace (handles[k].front(), k);
    }
  }
  std::cout << "Done." << std::endl;

  for (iter = set.begin(); iter != set.end(); ++iter)
    std::cout << *iter << ' ';
  std::cout << std::endl;

  // Check the case of swapping two sibling nodes.
  typedef CGAL::Multiset<int>  Int_set;
  typedef Int_set::iterator    Int_set_iter;

  Int_set                     iset;
  Int_set_iter                pos1, pos2, iitr;

  iset.insert (0);
  iset.insert (2000);
  iset.insert (1000);
  iset.insert (1500);
  pos1 = iset.insert (1250);
  iset.insert (1375);
  pos2 = iset.insert (1437);

  for (iitr = iset.begin(); iitr != iset.end(); ++iitr)
    std::cout << *iitr << ' ';
  std::cout << std::endl;
  assert (iset.is_valid());

  *pos1 = 1438;
  *pos2 = 1251;

  for (iitr = iset.begin(); iitr != iset.end(); ++iitr)
    std::cout << *iitr << ' ';
  std::cout << std::endl;
  // Note that the set is temporarily invalid!
  assert (! iset.is_valid());

  // Perform the swap to make it valid again.
  iset.swap (pos1, pos2);

  for (iitr = iset.begin(); iitr != iset.end(); ++iitr)
    std::cout << *iitr << ' ';
  std::cout << std::endl;
  assert (iset.is_valid());

  return;
}
