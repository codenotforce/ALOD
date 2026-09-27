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
    src/problems/problems.cpp
    src/fem/afem.cpp
    src/lod/quasi_interp.cpp
    src/lod/patches.cpp
    src/lod/slod.cpp
    src/lod/lod.cpp
    src/estimator/riesz.cpp
    src/estimator/marking.cpp
    src/localization/theta.cpp
    src/localization/dense.cpp
    src/mesh/nested.cpp
    src/adaptive/schedule.cpp
    src/adaptive/refine.cpp
    src/enrichment/aot.cpp
    src/enrichment/regional.cpp
    src/io/checkpoint.cpp
    src/fem/batch.cpp
    src/helmholtz/patch_system.cpp
    src/helmholtz/patch_solver.cpp
    src/helmholtz/corrector.cpp)
target_include_directories(alod_core PUBLIC include PRIVATE ${UMFPACK_INCLUDE_DIR})
target_compile_definitions(alod_core PRIVATE LOD2D_HAVE_UMFPACK=1)
target_link_libraries(alod_core PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX PRIVATE ${UMFPACK_LIBRARY})

add_executable(alod_fem_smoke apps/fem_smoke.cpp)
target_link_libraries(alod_fem_smoke PRIVATE alod_core)
add_executable(alod_run apps/alod_run.cpp apps/adaptive_run.cpp apps/audit_run.cpp apps/audit_worker.cpp apps/shared_audit.cpp)
target_link_libraries(alod_run PRIVATE alod_core)
add_executable(alod_fixed apps/alod_fixed.cpp)
target_link_libraries(alod_fixed PRIVATE alod_core)
add_executable(alod_audit apps/alod_audit.cpp apps/audit_run.cpp)
target_link_libraries(alod_audit PRIVATE alod_core)
add_executable(alod_geometry apps/alod_geometry.cpp)
target_link_libraries(alod_geometry PRIVATE alod_core)
add_executable(alod_profile apps/alod_profile.cpp)
target_link_libraries(alod_profile PRIVATE alod_core)

if(ALOD_BUILD_TESTS)
    enable_testing()
    add_executable(benchmark_spaces EXCLUDE_FROM_ALL tests/benchmark_spaces.cpp)
    target_link_libraries(benchmark_spaces PRIVATE alod_core)
    add_executable(stress_checkpoint EXCLUDE_FROM_ALL tests/stress_checkpoint.cpp)
    target_link_libraries(stress_checkpoint PRIVATE alod_core)
    foreach(name cross_state bordered_solve kernel_lift zero_kernel patch_recovery relative_theta integration space_assembly)
        add_executable(test_${name} tests/test_${name}.cpp)
        target_link_libraries(test_${name} PRIVATE alod_core)
        add_test(NAME v4_${name} COMMAND test_${name})
        set_tests_properties(v4_${name} PROPERTIES TIMEOUT 240 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
    endforeach()
    add_test(NAME campaign_resume COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_campaign.py $<TARGET_FILE:alod_run>)
    set_tests_properties(campaign_resume PROPERTIES TIMEOUT 180 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
    add_test(NAME asynchronous_restart COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_async_restart.py $<TARGET_FILE:alod_run>)
    add_test(NAME shared_audit_transport COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_shared_audit.py $<TARGET_FILE:alod_run>)
    set_tests_properties(shared_audit_transport PROPERTIES TIMEOUT 180 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
    add_test(NAME physical_training_coordinates COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_large_training.py $<TARGET_FILE:alod_run>)
    set_tests_properties(physical_training_coordinates PROPERTIES TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
    set_tests_properties(asynchronous_restart PROPERTIES TIMEOUT 360 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
    add_test(NAME v4_geometry COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_v4_geometry.py $<TARGET_FILE:alod_run>)
    set_tests_properties(v4_geometry PROPERTIES TIMEOUT 180 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
    add_test(NAME v4_initial_evidence COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_v4_initial.py $<TARGET_FILE:alod_run>)
    set_tests_properties(v4_initial_evidence PROPERTIES TIMEOUT 600 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
    add_test(NAME v4_normalized_ell COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_normalized_ell.py $<TARGET_FILE:alod_run>)
    set_tests_properties(v4_normalized_ell PROPERTIES TIMEOUT 600 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
    add_executable(test_p34 tests/test_p34.cpp)
    target_link_libraries(test_p34 PRIVATE alod_core)
    add_test(NAME p34_mathematics COMMAND test_p34)
    add_test(NAME p34_cli COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_adaptive_cli.py $<TARGET_FILE:alod_run>)
    set_tests_properties(p34_mathematics p34_cli PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1" TIMEOUT 240)
    add_executable(test_p2 tests/test_p2.cpp)
    target_link_libraries(test_p2 PRIVATE alod_core)
    add_test(NAME p2_mathematics COMMAND test_p2)
    set_tests_properties(p2_mathematics PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1" TIMEOUT 180)
    find_package(Python3 3.10 REQUIRED COMPONENTS Interpreter)
    add_executable(test_optimization_reuse tests/test_optimization_reuse.cpp)
    target_link_libraries(test_optimization_reuse PRIVATE alod_core)
    add_test(NAME optimization_reuse COMMAND test_optimization_reuse)
    set_tests_properties(optimization_reuse PROPERTIES TIMEOUT 180 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
    add_test(NAME optimization_workflow COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_optimization_workflow.py $<TARGET_FILE:alod_run>)
    set_tests_properties(optimization_workflow PROPERTIES TIMEOUT 240 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
    add_test(NAME workflow_basis_reuse COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_workflow_reuse.py $<TARGET_FILE:alod_run>)
    set_tests_properties(workflow_basis_reuse PROPERTIES TIMEOUT 180 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
    add_test(NAME variable_wavenumber COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_wavenumber.py $<TARGET_FILE:alod_run>)
    set_tests_properties(variable_wavenumber PROPERTIES TIMEOUT 360 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
    add_test(NAME p56_checkpoint_audit COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_checkpoint_cli.py $<TARGET_FILE:alod_run>)
    set_tests_properties(p56_checkpoint_audit PROPERTIES TIMEOUT 360 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
    foreach(fixture e1_prefix e2_initial)
        add_test(NAME p34_${fixture} COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_adaptive.py
            --executable $<TARGET_FILE:alod_run> --fixture ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/p34/${fixture}.json)
        set_tests_properties(p34_${fixture} PROPERTIES TIMEOUT 240 LABELS smoke)
    endforeach()
    foreach(problem E1 E2)
        foreach(state uniform graded ritz)
            add_test(NAME p2_${problem}_${state}
                COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_p2.py
                    --executable $<TARGET_FILE:alod_fixed>
                    --fixture ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/p2/${problem}_${state}.json)
            set_tests_properties(p2_${problem}_${state} PROPERTIES TIMEOUT 180 LABELS smoke)
        endforeach()
    endforeach()
    add_test(NAME p2_thread_reproducibility
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_p2.py
            --executable $<TARGET_FILE:alod_fixed>
            --fixture ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/p2/E2_graded.json --threads 2)
    set_tests_properties(p2_thread_reproducibility PROPERTIES TIMEOUT 180 LABELS smoke)
    add_test(NAME p2_cli COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_fixed_cli.py $<TARGET_FILE:alod_fixed>)
    set_tests_properties(p2_cli PROPERTIES TIMEOUT 120 LABELS smoke)
    add_executable(quasi_average_probe tests/quasi_average_probe.cpp)
    target_link_libraries(quasi_average_probe PRIVATE alod_core)
    add_test(NAME quasi_average_discrepancy COMMAND quasi_average_probe)
    add_executable(test_structural_reuse tests/test_structural_reuse.cpp)
    target_link_libraries(test_structural_reuse PRIVATE alod_core)
    add_test(NAME structural_reuse COMMAND test_structural_reuse)
    set_tests_properties(structural_reuse PROPERTIES TIMEOUT 180 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
    foreach(k 8 128)
        add_test(NAME structural_reuse_k${k} COMMAND test_structural_reuse ${k})
        set_tests_properties(structural_reuse_k${k} PROPERTIES TIMEOUT 180 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
    endforeach()
    add_executable(test_mathematics tests/test_mathematics.cpp)
    target_link_libraries(test_mathematics PRIVATE alod_core)
    add_test(NAME mathematical_regression COMMAND test_mathematics)
    set_tests_properties(mathematical_regression PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1" TIMEOUT 120)
    add_test(NAME repository_integrity
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_repository.py)
    add_test(NAME tooling_regression
        COMMAND ${Python3_EXECUTABLE} -m unittest discover -s tests -p test_*.py)
    set_tests_properties(tooling_regression PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    foreach(problem E1 E2)
        foreach(method AFEM UFEM SLOD)
            add_test(NAME baseline_${problem}_${method}
                COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_baselines.py
                    --executable $<TARGET_FILE:alod_run>
                    --fixture ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/baselines/${problem}_${method}.json)
            set_tests_properties(baseline_${problem}_${method} PROPERTIES TIMEOUT 180 LABELS smoke)
        endforeach()
        add_test(NAME fem_${problem}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_oracle.py
                --executable $<TARGET_FILE:alod_fem_smoke>
                --fixture ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/fem_${problem}.json)
        set_tests_properties(fem_${problem} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
    endforeach()
    add_test(NAME baseline_thread_reproducibility
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_baselines.py
            --executable $<TARGET_FILE:alod_run>
            --fixture ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/baselines/E2_SLOD.json --threads 2)
    set_tests_properties(baseline_thread_reproducibility PROPERTIES TIMEOUT 180 LABELS smoke)
    add_test(NAME baseline_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_cli.py $<TARGET_FILE:alod_run>)
    set_tests_properties(baseline_cli PROPERTIES TIMEOUT 60 LABELS smoke)
    set_tests_properties(repository_integrity tooling_regression fem_E1 fem_E2 PROPERTIES LABELS smoke TIMEOUT 120)
endif()
