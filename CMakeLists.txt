cmake_minimum_required(VERSION 3.20)
project(ALOD VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
option(ALOD_BUILD_TESTS "Build migration regression tests" ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED COMPONENTS CXX)
find_path(UMFPACK_INCLUDE_DIR umfpack.h PATH_SUFFIXES suitesparse REQUIRED)
find_library(UMFPACK_LIBRARY NAMES umfpack REQUIRED)

add_library(alod_core STATIC
    src/mesh/edges.cpp
    src/mesh/refine.cpp
    src/mesh/initial.cpp
    src/fem/boundary.cpp
    src/fem/quadrature.cpp
    src/fem/manufactured.cpp
    src/fem/operators.cpp
    src/problems/legacy_cases.cpp
    src/problems/problems.cpp)
target_include_directories(alod_core PUBLIC include PRIVATE ${UMFPACK_INCLUDE_DIR})
target_compile_definitions(alod_core PRIVATE LOD2D_HAVE_UMFPACK=1)
target_link_libraries(alod_core PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX PRIVATE ${UMFPACK_LIBRARY})

add_executable(alod_fem_smoke apps/fem_smoke.cpp)
target_link_libraries(alod_fem_smoke PRIVATE alod_core)

if(ALOD_BUILD_TESTS)
    enable_testing()
    find_package(Python3 3.10 REQUIRED COMPONENTS Interpreter)
    add_test(NAME repository_integrity
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_repository.py)
    add_test(NAME tooling_regression
        COMMAND ${Python3_EXECUTABLE} -m unittest discover -s tests -p test_*.py)
    set_tests_properties(tooling_regression PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    foreach(problem E1 E2)
        add_test(NAME fem_${problem}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_oracle.py
                --executable $<TARGET_FILE:alod_fem_smoke>
                --fixture ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/fem_${problem}.json)
        set_tests_properties(fem_${problem} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
    endforeach()
    set_tests_properties(repository_integrity tooling_regression fem_E1 fem_E2 PROPERTIES LABELS smoke TIMEOUT 120)
endif()
