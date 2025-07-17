cmake_minimum_required(VERSION 3.24...3.31)

include_guard(GLOBAL)

include(FetchContent)

FetchContent_Declare(
  tl-expected
  DOWNLOAD_EXTRACT_TIMESTAMP NO
  GIT_REPOSITORY https://github.com/TartanLlama/expected.git
  GIT_TAG        292eff8bd8ee230a7df1d6a1c00c4ea0eb2f0362 # https://github.com/TartanLlama/expected/releases/tag/v1.1.0
  FIND_PACKAGE_ARGS CONFIG
)

set(EXPECTED_BUILD_TESTS FALSE CACHE BOOL "Enable tl::expected tests (defaults to OFF for CGAL)")

FetchContent_MakeAvailable(tl-expected)
