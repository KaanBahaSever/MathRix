cmake_minimum_required(VERSION 3.16)
project(MathRix VERSION 0.1.0 DESCRIPTION "Header-only C++17 math library" LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(MATHRIX_TOP_LEVEL ON)
else()
    set(MATHRIX_TOP_LEVEL OFF)
endif()
option(MATHRIX_BUILD_TESTS "Build MathRix tests" ${MATHRIX_TOP_LEVEL})
option(MATHRIX_BUILD_EXAMPLES "Build MathRix examples" ${MATHRIX_TOP_LEVEL})

add_library(MathRix INTERFACE)
add_library(MathRix::MathRix ALIAS MathRix)
target_include_directories(MathRix INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_features(MathRix INTERFACE cxx_std_17)

if(MATHRIX_BUILD_TESTS)
    enable_testing()
    add_executable(mathrix_tests tests/main.cpp tests/test_geometry.cpp tests/test_linear_algebra.cpp
                                 tests/test_interpolation.cpp tests/test_ode.cpp)
    target_link_libraries(mathrix_tests PRIVATE MathRix::MathRix)
    if(MSVC)
        target_compile_options(mathrix_tests PRIVATE /W4 /permissive-)
    else()
        target_compile_options(mathrix_tests PRIVATE -Wall -Wextra -Wpedantic -Wshadow)
    endif()
    add_test(NAME mathrix_tests COMMAND mathrix_tests)
endif()

if(MATHRIX_BUILD_EXAMPLES)
    foreach(example example linear_system)
        add_executable(mathrix_${example} examples/${example}.cpp)
        target_link_libraries(mathrix_${example} PRIVATE MathRix::MathRix)
    endforeach()
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS MathRix EXPORT MathRixTargets)
install(EXPORT MathRixTargets NAMESPACE MathRix:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MathRix)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/MathRixConfig.cmake "include(\${CMAKE_CURRENT_LIST_DIR}/MathRixTargets.cmake)\n")
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/MathRixConfigVersion.cmake COMPATIBILITY SameMajorVersion ARCH_INDEPENDENT)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/MathRixConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/MathRixConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MathRix)
