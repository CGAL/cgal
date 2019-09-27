if(NOT TARGET benchmark::benchmark)
  include(ExternalProject)
  ExternalProject_Add(google-benchmark
    URL               https://github.com/google/benchmark/archive/v1.5.0.tar.gz
    URL_HASH          SHA512=a0df9aa3d03f676e302c76d83b436de36eea0a8517ab50a8f5a11c74ccc68a1f5128fa02474901002d8e6b5a4d290ef0272a798ff4670eab3e2d78dc86bb6cd3
    CMAKE_CACHE_ARGS  -DBENCHMARK_ENABLE_GTEST_TESTS:BOOL=FALSE -DCMAKE_INSTALL_PREFIX:STRING=<INSTALL_DIR>
    EXCLUDE_FROM_ALL  true
    LOG_INSTALL       true
    LOG_CONFIGURE     true
    )
  ExternalProject_Get_Property(google-benchmark INSTALL_DIR)
  find_package(benchmark HINTS ${INSTALL_DIR})
endif()
