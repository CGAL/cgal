#include <CGAL/assertions.h>
#include <CGAL/use.h>

#include <bitset>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// 0, nothing
// > 0, print RT_sufficient/FT_necessary errors and successes
// > 1, same as above + predicate being tested
// > 2, same as above + some general indications on what is going on
// > 4, same as above + even more indications on what is going on
// > 8, everything
#define CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY 2

std::vector<std::string> predicates_types = { };

// @todo, technically somebody could create predicates with non-kernel objects (nor FT/Origin), e.g. `int`.
// In that case, these arguments would have to be added to the lists below since there is no scrapping
// of the predicate arguments, but simply trying all combinations of objects from these lists.
std::vector<std::string> object_types_2 = { "FT", "Origin" };
std::vector<std::string> object_types_3 = { "FT", "Origin" };

// @todo potential operator()s with fewer than MIN_ARITY and more than MAX_ARITY are not tested
constexpr std::size_t MIN_ARITY = 0;
constexpr std::size_t MAX_ARITY = 12;

const std::string kernel_name = "Simple_cartesian";
const std::string FT_div = "double";
const std::string RT_no_div = "CGAL::Mpzf";

enum Needs_FT_checks
{
  NO_CHECK = 0,
  CHECK_NEEDS_FT,
  CHECK_NO_NEEDS_FT
};

enum Compilation_result
{
  SUCCESSFUL = 0, // if it got to linking, it is also a successful compilation
  FAILED_NO_MATCH,
  FAILED_AMBIGUOUS_CALL, // ambiguous calls means the arity is valid
  FAILED_NO_DIVISION_OPERATOR, // used to detect if a valid compilation can be done with RT
  FAILED_STATIC_ASSERTION, // used to detect failures in the result type checks
  UNKNOWN
};

inline const char* get_error_message(int error_code)
{
  // Messages corresponding to Error_code list above. Must be kept in sync!
  static const char* error_message[UNKNOWN+1] =
  {
    "Success!",
    "Failed: no match!",
    "Failed: ambiguous call!",
    "Failed: called division operator!",
    "Failed: static assertion violated!",
    "Unexpected error!"
  };

  if(error_code > UNKNOWN || error_code < 0)
    return "Doubly unexpected error!!";
  else
    return error_message[error_code];
}

std::string kernel_with_FT(const std::string& FT_name)
{
  return "CGAL::" + kernel_name + "<" + FT_name + ">";
}

// convert from e.g. Point_2 to CGAL::Point_2<Kernel::FT>
std::string parameter_with_namespace(const std::string& FT_name,
                                     const std::string& o)
{
  if(o == "Any")
    return "CGAL::Kernel_23_tests::Any";
  else if(o == "FT")
    return "K::FT";
  else if(o == "Origin")
    return "CGAL::Origin";
  else
    return "CGAL::" + o + "<" + kernel_with_FT(FT_name) + " >";
}

int compile()
{
#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 4)
  std::cout << "====== Compiling atomic file... ======" << std::endl;
#endif

  return std::system("cmake --build " CGAL_STRINGIZE(CMAKE_BINARY_DIR) " -t atomic_compilation_test > log.txt 2>&1");
}

Compilation_result parse_output(const std::string& predicate_name,
                                const std::string& FT_name = {},
                                const std::vector<std::string>& parameters = {})
{
#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 4)
  std::cout << "====== Parsing compilation log... ======" << std::endl;
#endif
  Compilation_result res = UNKNOWN;

  std::ifstream in("log.txt");
  if(!in)
  {
#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 0)
    std::cerr << "Error: failed to open log file" << std::endl;
#endif
    return UNKNOWN;
  }

#ifdef CGAL_KERNEL_23_TEST_RT_FT_PREDICATES_TEST_PREDICATES_WITH_TEMPLATED_OPERATORS
  // Compare_(squared)_distance_23 have templated operator()s, which are a lot of combinations to test.
  // In templated operator()s, the compare is simply a call to squared_distance()s and a CGAL::compare().
  // Below prunes some exploration branches in case the first squared_distance() call does not even compile.
  bool prune_compare_distance_branches = false;
  if(predicate_name == "Compare_distance_2" || predicate_name == "Compare_distance_3" ||
     predicate_name == "Compare_squared_distance_2" || predicate_name == "Compare_squared_distance_3")
  {
    prune_compare_distance_branches = true;
  }
#else
  CGAL_USE(predicate_name);
  CGAL_USE(FT_name);
  CGAL_USE(parameters);
#endif

  std::string line;
  while(getline(in, line))
  {
#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 8)
    std::cout << line << std::endl;
#endif

    if(line.find("no match for call") != std::string::npos) {
      res = FAILED_NO_MATCH;
      break;
    } else if(line.find("use of deleted function") != std::string::npos) {
      res = FAILED_NO_MATCH;
      break;
#ifdef CGAL_KERNEL_23_TEST_RT_FT_PREDICATES_TEST_PREDICATES_WITH_TEMPLATED_OPERATORS
    } else if(prune_compare_distance_branches && parameters.size() > 1 &&
              parameters[0] != "Any" && parameters[1] != "Any" &&
              line.find(std::string{"no matching function for call to ‘squared_distance(const " +
                                    parameter_with_namespace(FT_name, parameters[0]) + "&, const " +
                                    parameter_with_namespace(FT_name, parameters[1])}) != std::string::npos) {
      res = FAILED_NO_MATCH;
      break;
#endif
    } else if(line.find("ambiguous") != std::string::npos) {
      res = FAILED_AMBIGUOUS_CALL;
      break;
    } else if(line.find("candidate") != std::string::npos) {
      res = FAILED_AMBIGUOUS_CALL;
      break;
    } else if(line.find("no match for ‘operator/’") != std::string::npos) {
      res = FAILED_NO_DIVISION_OPERATOR;
      break;
    } else if(line.find("no match for ‘operator/=’") != std::string::npos) {
      res = FAILED_NO_DIVISION_OPERATOR;
      break;
    } else if(line.find("static assertion failed") != std::string::npos) {
      res = FAILED_STATIC_ASSERTION;
      break;
    } else if(line.find("Built") != std::string::npos) {
      res = SUCCESSFUL;
      break;
    } else if(line.find("undefined reference") != std::string::npos) {
      // Can happen because the conversion Any -> kernel object is not implemented
      res = SUCCESSFUL;
      break;
    }
  }

#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 4)
  std::cout << "Result of atomic test file is: " << get_error_message(res) << std::endl;
#endif
  CGAL_postcondition(res != UNKNOWN);

  return res;
}

void generate_atomic_compilation_test(const std::string& FT_name,
                                      const std::string& predicate_name,
                                      const std::vector<std::string>& parameters,
                                      const Needs_FT_checks check = NO_CHECK)
{
#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 4)
  std::cout << "\n====== Generate atomic compilation test... ======" << std::endl;
#endif

#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 4)
  std::cout << "\t" << predicate_name << "(";
  for(std::size_t j=0, i=parameters.size(); j<i; ++j)
    std::cout << ((j != 0) ? ", " : "") << parameter_with_namespace(FT_name, parameters[j]);
  std::cout << ")" << std::endl;
#endif

  std::ofstream out(CGAL_STRINGIZE(CMAKE_CURRENT_SOURCE_DIR) "/atomic_compilation_test.cpp");
  if(!out)
  {
    std::cerr << "Error: could not write into atomic compilation test" << std::endl;
    std::exit(1);
  }

  out << "#include \"atomic_RT_FT_predicate_headers.h\"\n";

  out << "using K = " + kernel_with_FT(FT_name) + ";\n";
  out << "using P = K::" << predicate_name << ";\n";
  out << "using B = P::result_type;\n";
  out << "using NFT_B = CGAL::Needs_FT<B>;\n";

  out << "int main(int, char**)\n";
  out << "{\n";

  out << "  P p{};\n";
  for(std::size_t j=0, i=parameters.size(); j<i; ++j)
    out << "  " << parameter_with_namespace(FT_name, parameters[j]) << " o" << j << "{};\n";

  out << "  p(";
  for(std::size_t j=0, i=parameters.size(); j<i; ++j)
    out << ((j != 0) ? ", o" : "o") << j;
  out << ");\n";

  if(check != NO_CHECK)
  {
    out << "  static_assert(std::is_same<decltype(";

    out << "p(";
    for(std::size_t j=0, i=parameters.size(); j<i; ++j)
      out << ((j != 0) ? ", o" : "o") << j;
    out << ")";

    if(check == CHECK_NO_NEEDS_FT)
      out << ", B>::value));\n";
    else if(check == CHECK_NEEDS_FT)
      out << ", NFT_B>::value));\n";
  }

  out << "  return EXIT_SUCCESS;\n";
  out << "}\n";
  out.close();
}

void ensure_NO_Needs_FT(const std::string& predicate_name,
                        const std::vector<std::string>& parameters)
{
#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 4)
  std::cout << predicate_name << "(";
  for(std::size_t j=0, i=parameters.size(); j<i; ++j)
    std::cout << ((j != 0) ? ", " : "") << parameters[j];
  std::cout << ") is RT-sufficient; check that the return type is *NOT* wrapped..." << std::endl;
#endif

  // RT is sufficient, check that `Needs_FT` is not in the operator()'s return type
  generate_atomic_compilation_test(RT_no_div, predicate_name, parameters, CHECK_NO_NEEDS_FT);
  auto compilation_result = compile();
  Compilation_result res = compilation_result == 0 ?
                           SUCCESSFUL :
                           parse_output(predicate_name);

  if(res == SUCCESSFUL)
  {
#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 0)
    std::cout << predicate_name << "(";
    for(std::size_t j=0, i=parameters.size(); j<i; ++j)
      std::cout << ((j != 0) ? ", " : "") << parameters[j];
    std::cout << ") is RT-sufficient, and the wrap `Needs_FT` is (correctly) absent!" << std::endl;
#endif
  }
  else if(res == FAILED_STATIC_ASSERTION)
  {
#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 0)
    std::cout << "Error: " << predicate_name << "(";
    for(std::size_t j=0, i=parameters.size(); j<i; ++j)
      std::cout << ((j != 0) ? ", " : "") << parameters[j];
    std::cout << ") is RT-sufficient, but the return type is wrong (superfluous `Needs_FT`?)!" << std::endl;
#endif
  }
  else
  {
#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 0)
    std::cerr << "Unexpected error during Needs_FT checks" << std::endl;
#endif
    assert(false);
  }
}

void ensure_Needs_FT(const std::string& predicate_name,
                     const std::vector<std::string>& parameters)
{
#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 4)
  std::cout << predicate_name << "(";
  for(std::size_t j=0, i=parameters.size(); j<i; ++j)
    std::cout << ((j != 0) ? ", " : "") << parameters[j];
  std::cout << ") is FT-necessary; check that the return type is wrapped..." << std::endl;
#endif

  // The predicate requires a FT with division, ensure that Needs_FT is present in the operator()'s return type
  generate_atomic_compilation_test(FT_div, predicate_name, parameters, CHECK_NEEDS_FT);
  auto compilation_result = compile();
  Compilation_result res = compilation_result == 0 ?
                           SUCCESSFUL :
                           parse_output(predicate_name);

  if(res == SUCCESSFUL)
  {
#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 0)
    std::cout << predicate_name << "(";
    for(std::size_t j=0, i=parameters.size(); j<i; ++j)
      std::cout << ((j != 0) ? ", " : "") << parameters[j];
    std::cout << ") is FT-necessary, and the wrap `Needs_FT` is (correctly) present!" << std::endl;
#endif
  }
  else if(res == FAILED_STATIC_ASSERTION)
  {
#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 0)
    std::cout << "Error: " << predicate_name << "(";
    for(std::size_t j=0, i=parameters.size(); j<i; ++j)
      std::cout << ((j != 0) ? ", " : "") << parameters[j];
    std::cout << ") is FT-necessary, but the wrap `Needs_FT` is missing!" << std::endl;
#endif
  }
  else
  {
#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 0)
    std::cerr << "Unexpected error during Needs_FT checks" << std::endl;
#endif
    assert(false);
  }
}

void test_predicate(const std::string& predicate_name,
                    const std::size_t object_pos,
                    const std::size_t arity,
                    // intentional copy, each sub-branch gets its own parameter list
                    std::vector<std::string> parameters)
{
  const std::size_t last = arity - 1;
  CGAL_precondition(object_pos <= last);

  CGAL_precondition(predicate_name.back() == '2' || predicate_name.back() == '3');
  const auto& object_types = (predicate_name.back() == '2') ? object_types_2 : object_types_3;

  for(const std::string& object_type : object_types)
  {
#ifdef CGAL_KERNEL_23_TEST_RT_FT_PREDICATES_TEST_PREDICATES_WITH_TEMPLATED_OPERATORS
    // This pruning could be done for other predicates, but they're not as expensive so it doesn't matter
    if((predicate_name == "Compare_distance_2" || predicate_name == "Compare_distance_3" ||
        predicate_name == "Compare_squared_distance_2" || predicate_name == "Compare_squared_distance_3" ||
        predicate_name == "Do_intersect_2" || predicate_name == "Do_intersect_3") &&
       object_type == "FT")
    {
      continue;
    }
#endif

    parameters[object_pos] = object_type;
    generate_atomic_compilation_test(RT_no_div, predicate_name, parameters);
    auto compilation_result = compile();
    Compilation_result res = compilation_result == 0 ?
                             SUCCESSFUL :
                             parse_output(predicate_name, RT_no_div, parameters);

    // See if we can already (i.e., possibly with `Any`s) conclude on the current parameter list
    if(res == FAILED_NO_MATCH)
    {
      // The object at the current position yields a compilation error, do not explore any further
      continue;
    }
    else if(object_pos == last)
    {
      if(res == SUCCESSFUL)
      {
        ensure_NO_Needs_FT(predicate_name, parameters);
      }
      else if(res == FAILED_NO_DIVISION_OPERATOR)
      {
        ensure_Needs_FT(predicate_name, parameters);
      }
    }
    else
    {
      // The object at the current position does not invalid the call, explore further this list
      test_predicate(predicate_name, object_pos + 1, arity, parameters);
    }
  }
}

void test_predicate(const std::string& predicate_name,
                    const std::size_t arity)
{
#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 2)
  std::cout << "\n\n==== Test predicate with arity " << arity << "... ====" << std::endl;
#endif

  // Use "Any" to prune early:
  // 1st try "Object_1, Any, ..., Any" (i - 1 "Any")
  // -> if that doesn't compile, we're done with Object_1 and try "Object_2, Any, ..., Any" (i-1 "Any")
  // -> if that compiles, try "Object_1, Object_1, Any, ..., Any" (i-2 "Any")
  // etc.

  // the position of the object being changed/tested, when object_pos == arity - 1,
  // then this is a call with full objects on which we can do the RT test
  std::vector<std::string> parameters(arity, "Any");

  // Quick try to see if it even matches anything
  generate_atomic_compilation_test(RT_no_div, predicate_name, parameters);
  auto compilation_result = compile();
  Compilation_result res = compilation_result == 0 ?
                           SUCCESSFUL :
                           parse_output(predicate_name);
  if(res == FAILED_NO_MATCH) // No point with this current arity
    return;

  std::size_t object_pos = 0;
  test_predicate(predicate_name, object_pos, arity, parameters);
}

void test_predicate(const std::string& predicate_name)
{
#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 1)
  std::cout << "\n\n\n== Test predicate: " << predicate_name << "... ==" << std::endl;
#endif

#ifndef CGAL_KERNEL_23_TEST_RT_FT_PREDICATES_TEST_PREDICATES_WITH_TEMPLATED_OPERATORS
  if(predicate_name == "Compare_distance_2" || predicate_name == "Compare_distance_3" ||
     predicate_name == "Compare_squared_distance_2" || predicate_name == "Compare_squared_distance_3" ||
     predicate_name == "Do_intersect_2" || predicate_name == "Do_intersect_3")
  {
#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 1)
    std::cout << "Skipping because 'CGAL_KERNEL_23_TEST_RT_FT_PREDICATES_TEST_PREDICATES_WITH_TEMPLATED_OPERATORS' is not defined!" << std::endl;
#endif
    return;
  }
#endif

  for(std::size_t i=MIN_ARITY; i<=MAX_ARITY; ++i)
    test_predicate(predicate_name, i);
}

// Just to not get a diff at the end of the test
void restore_atomic_file()
{
  std::ofstream out("../atomic_compilation_test.cpp");
  if(!out)
  {
    std::cerr << "Error: could not write into atomic compilation test" << std::endl;
    std::exit(1);
  }

  out << "// This executable is used by test_RT_or_FT_predicates.cpp\n";
  out << "int main(int, char**) { }\n";
  out.close();
}

int main(int , char**)
{
  // Get the predicates
  #define CGAL_Kernel_pred(X, Y) predicates_types.push_back(#X);
  #include <CGAL/Kernel/interface_macros.h>

  // Get the objects
  #define CGAL_Kernel_obj(X) { const std::string O = #X; \
                               CGAL_precondition(O.back() == '2' || O.back() == '3'); \
                               if(O.back() == ('2')) \
                                 object_types_2.push_back(#X); \
                               else \
                                 object_types_3.push_back(#X); }
  #include <CGAL/Kernel/interface_macros.h>

#if (CGAL_KERNEL_23_TEST_RT_FT_VERBOSITY > 1)
  std::cout << predicates_types.size() << " predicates:" << std::endl;
  for(const std::string& s : predicates_types)
    std::cout << s << "\n";
  std::cout << std::endl;

  std::cout << object_types_2.size() << " 2D objects:" << std::endl;
  for(const std::string& o : object_types_2)
    std::cout << o << "\n";
  std::cout << std::endl;

  std::cout << object_types_3.size() << " 3D objects:" << std::endl;
  for(const std::string& o : object_types_3)
    std::cout << o << "\n";
  std::cout << std::endl;
#endif

  // Actual tests
  for(const std::string& predicate_name : predicates_types)
  {
    test_predicate(predicate_name);
  }

  restore_atomic_file();

  return EXIT_SUCCESS;
}
