cmake_minimum_required(VERSION 2.9)

project(StraightSkel)

include_directories(src)
include_directories(test)

# set default options
# to configure it, use:
# $ ccmake .
option(BUILD_OWN_KERNEL "Build our own Kernel (double precision)" ON)
option(USE_CGAL "Use the CGAL Kernel" OFF)
option(BUILD_TESTS "Use Boost's Unit Test Framework to test the implementation" OFF)
option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)

# it may be useful to set CGAL_DONT_OVERRIDE_CMAKE_FLAGS
# /usr/lib/CGAL/CGAL_SetupFlags.cmake
# set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS ON)

if (CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "== GCC C++ Compiler detected ==")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()


# create config.h
if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(DEBUG ON)
endif()
configure_file(${CMAKE_SOURCE_DIR}/src/config.h.in ${CMAKE_SOURCE_DIR}/src/config.h)


# to get help about a cmake module run:
# $ cmake --help-module FindBoost
if(BUILD_TESTS)
    set(BOOST_COMPONENTS ${BOOST_COMPONENTS} unit_test_framework)
    find_package(Boost REQUIRED ${BOOST_COMPONENTS})
endif()


# build util
set(util_SOURCES
        src/util/StringFactory.cpp
        src/util/StringFuncs.cpp
        src/util/StackTrace.cpp
        src/util/Configuration.cpp
        src/util/Timer.cpp)
add_library(skelutil SHARED ${util_SOURCES})


if(BUILD_OWN_KERNEL)
    set(kernel_SOURCES
            src/kernel/Point2.cpp
            src/kernel/Vector2.cpp
            src/kernel/Line2.cpp
            src/kernel/Segment2.cpp
            src/kernel/Ray2.cpp
            src/kernel/Point3.cpp
            src/kernel/Vector3.cpp
            src/kernel/Line3.cpp
            src/kernel/Segment3.cpp
            src/kernel/Plane3.cpp
            src/kernel/Sphere3.cpp
            src/kernel/intersection.cpp
            src/kernel/distance.cpp
            src/kernel/bisector.cpp
            src/kernel/projection.cpp)
    add_library(skelkernel SHARED ${kernel_SOURCES})
    target_link_libraries(skelkernel m)
endif()


# build data structures
set(data_2d_SOURCES
        src/data/2d/KernelFactory.cpp
        src/data/2d/Vertex.cpp
        src/data/2d/VertexData.cpp
        src/data/2d/Edge.cpp
        src/data/2d/EdgeData.cpp
        src/data/2d/Polygon.cpp
        src/data/2d/skel/StraightSkeleton.cpp
        src/data/2d/skel/Node.cpp
        src/data/2d/skel/Arc.cpp
        src/data/2d/skel/AbstractEvent.cpp
        src/data/2d/skel/ConstOffsetEvent.cpp
        src/data/2d/skel/EdgeEvent.cpp
        src/data/2d/skel/SplitEvent.cpp
        src/data/2d/skel/TriangleEvent.cpp
        src/data/2d/skel/SkelVertexData.cpp
        src/data/2d/skel/SkelEdgeData.cpp
        src/data/2d/mesh/Mesh.cpp
        src/data/2d/mesh/MeshVertex.cpp
        src/data/2d/mesh/MeshCell.cpp
        src/data/2d/mesh/MeshRay.cpp
        src/data/2d/mesh/MeshEdgeData.cpp)
set(data_3d_SOURCES
        src/data/3d/KernelFactory.cpp
        src/data/3d/Vertex.cpp
        src/data/3d/VertexData.cpp
        src/data/3d/Edge.cpp
        src/data/3d/EdgeData.cpp
        src/data/3d/Facet.cpp
        src/data/3d/FacetData.cpp
        src/data/3d/Triangle.cpp
        src/data/3d/Polyhedron.cpp
        src/data/3d/CircularVertex.cpp
        src/data/3d/CircularVertexData.cpp
        src/data/3d/CircularEdge.cpp
        src/data/3d/CircularEdgeData.cpp
        src/data/3d/SphericalPolygon.cpp
        src/data/3d/skel/StraightSkeleton.cpp
        src/data/3d/skel/Node.cpp
        src/data/3d/skel/Arc.cpp
        src/data/3d/skel/Sheet.cpp
        src/data/3d/skel/AbstractEvent.cpp
        src/data/3d/skel/ConstOffsetEvent.cpp
        src/data/3d/skel/SaveOffsetEvent.cpp
        src/data/3d/skel/EdgeEvent.cpp
        src/data/3d/skel/EdgeMergeEvent.cpp
        src/data/3d/skel/TriangleEvent.cpp
        src/data/3d/skel/DblEdgeMergeEvent.cpp
        src/data/3d/skel/DblTriangleEvent.cpp
        src/data/3d/skel/TetrahedronEvent.cpp
        src/data/3d/skel/VertexEvent.cpp
        src/data/3d/skel/FlipVertexEvent.cpp
        src/data/3d/skel/SurfaceEvent.cpp
        src/data/3d/skel/PolyhedronSplitEvent.cpp
        src/data/3d/skel/SplitMergeEvent.cpp
        src/data/3d/skel/EdgeSplitEvent.cpp
        src/data/3d/skel/PierceEvent.cpp
        src/data/3d/skel/SkelVertexData.cpp
        src/data/3d/skel/SkelEdgeData.cpp
        src/data/3d/skel/SkelFacetData.cpp
        src/data/3d/skel/SphericalSkeleton.cpp
        src/data/3d/skel/CircularNode.cpp
        src/data/3d/skel/CircularArc.cpp
        src/data/3d/skel/SphericalAbstractEvent.cpp
        src/data/3d/skel/SphericalConstOffsetEvent.cpp
        src/data/3d/skel/SphericalEdgeEvent.cpp
        src/data/3d/skel/SphericalSplitEvent.cpp
        src/data/3d/skel/SphericalTriangleEvent.cpp
        src/data/3d/skel/SphericalDblEdgeEvent.cpp
        src/data/3d/skel/SphericalLeaveEvent.cpp
        src/data/3d/skel/SphericalReturnEvent.cpp
        src/data/3d/skel/SphericalDblLeaveEvent.cpp
        src/data/3d/skel/SphericalDblReturnEvent.cpp
        src/data/3d/skel/SphericalVertexEvent.cpp
        src/data/3d/skel/SphericalEdgeMergeEvent.cpp
        src/data/3d/skel/SphericalInversionEvent.cpp
        src/data/3d/skel/SphericalSkelVertexData.cpp
        src/data/3d/skel/SphericalSkelEdgeData.cpp
        src/data/3d/skel/SphericalOffset.cpp)
add_library(skeldata2d SHARED ${data_2d_SOURCES})
add_library(skeldata3d SHARED ${data_3d_SOURCES})
if(USE_CGAL)
    find_package(CGAL REQUIRED Core)
    if(CGAL_FOUND)
        include(${CGAL_USE_FILE})
        target_link_libraries(skeldata2d skelutil ${Boost_LIBRARIES} ${CGAL_LIBRARIES})
        target_link_libraries(skeldata3d skeldata2d skelutil ${Boost_LIBRARIES} ${CGAL_LIBRARIES})
    endif()
else()
    target_link_libraries(skeldata2d skelutil skelkernel ${Boost_LIBRARIES})
    target_link_libraries(skeldata3d skeldata2d skelutil skelkernel ${Boost_LIBRARIES})
endif()


# build database connector
set(db_SOURCES
        src/db/SQLiteDatabase.cpp
        src/db/SQLiteStmt.cpp)
set(db_2d_SOURCES
        src/db/2d/DAOFactory.cpp
        src/db/2d/PointDAO.cpp
        src/db/2d/VertexDAO.cpp
        src/db/2d/EdgeDAO.cpp
        src/db/2d/PolygonDAO.cpp
        src/db/2d/NodeDAO.cpp
        src/db/2d/ArcDAO.cpp
        src/db/2d/EventDAO.cpp
        src/db/2d/StraightSkeletonDAO.cpp
        src/db/2d/AbstractFile.cpp
        src/db/2d/FLMAFile.cpp)
set(db_3d_SOURCES
        src/db/3d/DAOFactory.cpp
        src/db/3d/PointDAO.cpp
        src/db/3d/PlaneDAO.cpp
        src/db/3d/VertexDAO.cpp
        src/db/3d/EdgeDAO.cpp
        src/db/3d/TriangleDAO.cpp
        src/db/3d/FacetDAO.cpp
        src/db/3d/PolyhedronDAO.cpp
        src/db/3d/NodeDAO.cpp
        src/db/3d/ArcDAO.cpp
        src/db/3d/SheetDAO.cpp
        src/db/3d/EventDAO.cpp
        src/db/3d/StraightSkeletonDAO.cpp
        src/db/3d/AbstractFile.cpp
        src/db/3d/OBJFile.cpp
        src/db/3d/FLMAFile.cpp)
add_library(skeldb SHARED ${db_SOURCES} ${db_2d_SOURCES} ${db_3d_SOURCES})
target_link_libraries(skeldb skeldata2d skeldata3d skelutil sqlite3 dl)


# build algorithm
set(algo_SOURCES
        src/algo/Controller.cpp)
set(algo_2d_SOURCES
        src/algo/2d/KernelWrapper.cpp
        src/algo/2d/PolygonTransformation.cpp
        src/algo/2d/SimpleStraightSkel.cpp
        src/algo/2d/FastStraightSkel.cpp
        src/algo/2d/MeshModifier.cpp
        src/algo/2d/SkelMeshGenerator.cpp)
set(algo_3d_SOURCES
        src/algo/3d/KernelWrapper.cpp
        src/algo/3d/SimpleStraightSkel.cpp
        src/algo/3d/PolyhedronBuilder.cpp
        src/algo/3d/LineInFacet.cpp
        src/algo/3d/SelfIntersection.cpp
        src/algo/3d/AbstractVertexSplitter.cpp
        src/algo/3d/AngleVertexSplitter.cpp
        src/algo/3d/CombiVertexSplitter.cpp
        src/algo/3d/ConvexVertexSplitter.cpp
        src/algo/3d/VolumeVertexSplitter.cpp
        src/algo/3d/WeightVertexSplitter.cpp
        src/algo/3d/SphereVertexSplitter.cpp
        src/algo/3d/PolyhedronTransformation.cpp
        src/algo/3d/PolyhedronIntersection.cpp
        src/algo/3d/AbstractSimpleSphericalSkel.cpp
        src/algo/3d/ProjSimpleSphericalSkel.cpp
        src/algo/3d/RotSimpleSphericalSkel.cpp
        src/algo/3d/TransSimpleSphericalSkel.cpp
        src/algo/3d/SpeedSimpleSphericalSkel.cpp
        src/algo/3d/GraphChecker.cpp)
add_library(skelalgo2d SHARED ${algo_SOURCES} ${algo_2d_SOURCES})
add_library(skelalgo3d SHARED ${algo_3d_SOURCES})
target_link_libraries(skelalgo2d skeldata2d)
target_link_libraries(skelalgo3d skelalgo2d skeldata3d)


# build user interface
set(ui_SOURCES
        src/ui/vecmath.cpp)
set(ui_gl_SOURCES
        src/ui/gl/Camera.cpp
        src/ui/gl/OpenGLWindow.cpp
        src/ui/gl/MainOpenGLWindow.cpp
        src/ui/gl/KeyboardAdapter.cpp
        src/ui/gl/MouseAdapter.cpp)
set(ui_ps_SOURCES
        src/ui/ps/PSPrinter.cpp
        src/ui/ps/SpacePSPrinter.cpp
        src/ui/ps/PlanePSPrinter.cpp
        src/ui/ps/CutPatternPrinter.cpp)
add_library(skelui SHARED ${ui_SOURCES} ${ui_gl_SOURCES} ${ui_ps_SOURCES})
target_link_libraries(skelui skelutil m GL GLU glut ${Boost_LIBRARIES})


# build executable
set(StraightSkel_SOURCES
        src/main.cpp)
add_executable(StraightSkel ${StraightSkel_SOURCES})
target_link_libraries(StraightSkel skeldata2d skelalgo2d skeldata3d skelalgo3d skeldb skelui)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/res/StraightSkel.ini
    ${CMAKE_CURRENT_BINARY_DIR}/StraightSkel.ini COPYONLY)


if(BUILD_TESTS)
    # run the tests:
    # $ make test
    enable_testing()
    if(Boost_FOUND)
        include_directories(${Boost_INCLUDE_DIRS})
        set(UtilTest_SOURCES
                test/util/UtilTestRunner.cpp
                test/util/StringFuncsTest.cpp)
        add_executable(UtilTestRunner ${UtilTest_SOURCES})
        target_link_libraries(UtilTestRunner skelutil ${Boost_LIBRARIES})
        add_test(UtilTestRunner UtilTestRunner)
        if(BUILD_OWN_KERNEL)
            set(KernelTest_SOURCES
                    test/kernel/KernelTestRunner.cpp
                    test/kernel/Vector2Test.cpp
                    test/kernel/Point2Test.cpp
                    test/kernel/Line2Test.cpp
                    test/kernel/Segment2Test.cpp
                    test/kernel/Vector3Test.cpp
                    test/kernel/Point3Test.cpp
                    test/kernel/Line3Test.cpp
                    test/kernel/Segment3Test.cpp
                    test/kernel/Plane3Test.cpp
                    test/kernel/IntersectionTest.cpp
                    test/kernel/DistanceTest.cpp
                    test/kernel/BisectorTest.cpp)
            add_executable(KernelTestRunner ${KernelTest_SOURCES})
            target_link_libraries(KernelTestRunner skelkernel ${Boost_LIBRARIES})
            add_test(KernelTestRunner KernelTestRunner)
        endif()
        set(Data2DTest_SOURCES
                test/data/2d/Data2DTestRunner.cpp
                test/data/2d/KernelFactoryTest.cpp
                test/data/2d/VertexTest.cpp
                test/data/2d/EdgeTest.cpp
                test/data/2d/PolygonTest.cpp)
        add_executable(Data2DTestRunner ${Data2DTest_SOURCES})
        target_link_libraries(Data2DTestRunner skeldata2d ${Boost_LIBRARIES})
        add_test(Data2DTestRunner Data2DTestRunner)
        set(Data3DTest_SOURCES
                test/data/3d/Data3DTestRunner.cpp
                test/data/3d/EdgeTest.cpp
                test/data/3d/FacetTest.cpp
                test/data/3d/PolyhedronTest.cpp)
        add_executable(Data3DTestRunner ${Data3DTest_SOURCES})
        target_link_libraries(Data3DTestRunner skeldata3d ${Boost_LIBRARIES})
        add_test(Data3DTestRunner Data3DTestRunner)
        set(DBTest_SOURCES
                test/db/DBTestRunner.cpp
                test/db/2d/VertexDAOTest.cpp
                test/db/2d/EdgeDAOTest.cpp
                test/db/2d/PolygonDAOTest.cpp
                test/db/3d/PolyhedronDAOTest.cpp)
        add_executable(DBTestRunner ${DBTest_SOURCES})
        target_link_libraries(DBTestRunner skeldb ${Boost_LIBRARIES})
        add_test(DBTestRunner DBTestRunner)
        set(Algo2DTest_SOURCES
                test/algo/2d/Algo2DTestRunner.cpp
                test/algo/2d/KernelWrapperTest.cpp)
        add_executable(Algo2DTestRunner ${Algo2DTest_SOURCES})
        target_link_libraries(Algo2DTestRunner skelalgo2d ${Boost_LIBRARIES})
        add_test(Algo2DTestRunner Algo2DTestRunner)
        set(Algo3DTest_SOURCES
                test/algo/3d/Algo3DTestRunner.cpp
                test/algo/3d/KernelWrapperTest.cpp
                test/algo/3d/PolyhedronBuilderTest.cpp
                test/algo/3d/LineInFacetTest.cpp
                test/algo/3d/CombiVertexSplitterTest.cpp
                test/algo/3d/VolumeVertexSplitterTest.cpp)
        add_executable(Algo3DTestRunner ${Algo3DTest_SOURCES})
        target_link_libraries(Algo3DTestRunner skelalgo3d ${Boost_LIBRARIES})
        add_test(Algo3DTestRunner Algo3DTestRunner)
    endif()
endif()


if(BUILD_DOCUMENTATION)
    # build documentation:
    # $ make doc
    find_package(Doxygen REQUIRED)
    if(DOXYGEN_FOUND)
        set(DOXYGEN_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/doc)
        set(DOXYGEN_INPUT ${PROJECT_SOURCE_DIR}/src)
        configure_file(Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile @ONLY)
        add_custom_target(doc
                          COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
                          SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
    endif()
endif()


install(TARGETS skelutil DESTINATION lib${LIB_SUFFIX})
if(BUILD_OWN_KERNEL)
    install(TARGETS skelkernel DESTINATION lib${LIB_SUFFIX})
endif()
install(TARGETS skeldata2d DESTINATION lib${LIB_SUFFIX})
install(TARGETS skeldata3d DESTINATION lib${LIB_SUFFIX})
install(TARGETS skeldb DESTINATION lib${LIB_SUFFIX})
install(TARGETS skelalgo2d DESTINATION lib${LIB_SUFFIX})
install(TARGETS skelalgo3d DESTINATION lib${LIB_SUFFIX})
install(TARGETS skelui DESTINATION lib${LIB_SUFFIX})
install(TARGETS StraightSkel DESTINATION bin)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/StraightSkel.ini DESTINATION etc/StraightSkel)
