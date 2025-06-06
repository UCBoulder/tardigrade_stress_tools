#===================================================================================================== PROJECT SETUP ===
cmake_minimum_required(VERSION 3.18)
project(tardigrade_stress_tools)

# Set common project paths relative to project root directory
set(CPP_SRC_PATH "src/cpp")
set(CPP_TEST_PATH "${CPP_SRC_PATH}/tests")
set(CMAKE_SRC_PATH "src/cmake")
set(CMAKE_SRC_PATH "src/cmake")
set(DOXYGEN_SRC_PATH "docs/doxygen")
set(SPHINX_SRC_PATH "docs/sphinx")
# Set the internal support libraries
set(INTERNAL_SUPPORT_LIBRARIES "linear_elasticity")
set(ADDITIONAL_HEADER_ONLY_LIBRARIES "tardigrade_mass_change_deformation")
set(PROJECT_SOURCE_FILES ${PROJECT_NAME}.cpp ${PROJECT_NAME}.h)
set(PROJECT_PRIVATE_HEADERS "")
foreach(package ${INTERNAL_SUPPORT_LIBRARIES})
    set(PROJECT_SOURCE_FILES ${PROJECT_SOURCE_FILES} ${package}.cpp ${package}.h)
    set(PROJECT_PRIVATE_HEADERS ${PROJECT_PRIVATE_HEADERS} ${package}.h)
endforeach(package)
set(PROJECT_LINK_LIBRARIES tardigrade_error_tools tardigrade_constitutive_tools)

# Added a flag for whether the python bindings for stress tools get built
set(TARDIGRADE_STRESS_TOOLS_BUILD_PYTHON_BINDINGS ON CACHE BOOL "Flag for whether the python bindings should be built")

# Add a flag for if a full build of all tardigrade repositories should be performed
set(TARDIGRADE_FULL_BUILD OFF CACHE BOOL "Flag for whether a full build of Tardigrade should be performed (i.e., all repos pulled from git and built)")

# Add the cmake folder to locate project CMake module(s)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}" ${CMAKE_MODULE_PATH})

# Set build type checks
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_lower)
set(upstream_required "")
set(not_conda_test "true")
if(cmake_build_type_lower STREQUAL "release")
    set(upstream_required "REQUIRED")
elseif(cmake_build_type_lower STREQUAL "conda-test")
    set(upstream_required "REQUIRED")
    set(not_conda_test "false")
    # Find the installed project libraries
    find_package(${PROJECT_NAME} ${upstream_required} CONFIG)
    # Find the installed project umat
endif()

if(${not_conda_test} STREQUAL "true")
    foreach(package ${ADDITIONAL_HEADER_ONLY_LIBRARIES})
        add_library(${package} INTERFACE "${CPP_SRC_PATH}/${package}.cpp" "${CPP_SRC_PATH}/${package}.h")
    endforeach(package)
endif()

# Get version number from Git
set(VERSION_UPDATE_FROM_GIT True)
if(${not_conda_test} STREQUAL "true")
    # TODO: On osx-arm64 some c++ projects struggle to find $PREFIX/bin/python during conda-builds and others don't.
    # Trace the source and remove the hint when CMake configures with the correct python version during conda-build.
    if (DEFINED ENV{PREFIX})
        set(Python_ROOT_DIR "$ENV{PREFIX}/bin")
    endif()
    set(Python_FIND_STRATEGY LOCATION)
    find_package(Python COMPONENTS Interpreter REQUIRED)
    execute_process(COMMAND ${Python_EXECUTABLE} -m setuptools_scm
                            OUTPUT_VARIABLE ${PROJECT_NAME}_VERSION_STRING_FULL
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX MATCH "[0-9]+\.[0-9]+\.[0-9]+" ${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_VERSION_STRING_FULL})
    project(${PROJECT_NAME} VERSION ${${PROJECT_NAME}_VERSION})
endif()

# Add installation directory variables
include(GNUInstallDirs)

# Make the code position independent
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set the c++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lrt -O3")
endif()

# Enable CTest
enable_testing()

#================================================================================================= FIND DEPENDENCIES ===
# Find eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
if(EIGEN3_FOUND)
    message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
endif()

# Save the eigen directory for use in the python interface
set(EIGEN_DIR ${EIGEN3_INCLUDE_DIR} CACHE PATH "The path to the eigen include directory")

# Find related, but external, projects in installed environment
include(FetchContent)
set(upstream_packages "tardigrade_error_tools" "tardigrade_vector_tools" "tardigrade_constitutive_tools")
foreach(package ${upstream_packages})
    string(TOUPPER "${package}" package_upper)
    set(${package_upper}_BUILD_PYTHON_BINDINGS ${TARDIGRADE_STRESS_TOOLS_BUILD_PYTHON_BINDINGS} CACHE INTERNAL "Setting ${package}'s python binding flag to the global value")
    if(TARDIGRADE_FULL_BUILD)
        message("Not attempting to find ${package}")
    else()
        find_package(${package} ${upstream_required} CONFIG)
    endif()
    if(${package}_FOUND)
        message(STATUS "Found ${package}: ${${package}_DIR}")
    else()
        # Find related, but external, projects using FetchContent and building locally
        message(WARNING "Did not find an installed ${package} package. Attempting local build with FetchContent.")
        if(NOT DEFINED ${package_upper}_FETCHCONTENT_VERSION)
            set(${package_upper}_FETCHCONTENT_VERSION "origin/dev")
        endif()
        message("${package_upper} is being built with version ${${package_upper}_FETCHCONTENT_VERSION}")
        FetchContent_Declare(
            ${package}
            GIT_REPOSITORY http://github.com/UCBoulder/${package}.git
            GIT_TAG        ${${package_upper}_FETCHCONTENT_VERSION}
        )
        FetchContent_MakeAvailable(${package})
    endif()
endforeach(package)

#=============================================================================================== ADD PROJECT TARGETS ===
# MUST COME AFTER DEPENDENCY LOCATING
# Add project source directories
if(${not_conda_test} STREQUAL "true")
    include_directories("${CPP_SRC_PATH}")
    add_subdirectory(${CPP_SRC_PATH})
endif()

# Only add tests and documentation for current project builds. Protects downstream project builds.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Find Boost. Required for tests
    find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework)
    # Add tests and docs
    add_subdirectory(${CPP_TEST_PATH})
    if(${not_conda_test} STREQUAL "true")
        add_subdirectory(${DOXYGEN_SRC_PATH})
        add_subdirectory(${SPHINX_SRC_PATH})
    endif()
endif()

#==================================================================================== SETUP INSTALLATION CMAKE FILES ===
if(${not_conda_test} STREQUAL "true")
    foreach(package ${PROJECT_NAME} ${UMAT})
        include(CMakePackageConfigHelpers)
        write_basic_package_version_file("${package}ConfigVersion.cmake"
                                         VERSION ${PROJECT_VERSION}
                                         COMPATIBILITY SameMajorVersion)
        configure_package_config_file(
          "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}/Config.cmake.in"
          "${PROJECT_BINARY_DIR}/${package}Config.cmake"
          INSTALL_DESTINATION
          ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)

        # CMake won't build the targets for local builds of upstream projects
        if(cmake_build_type_lower STREQUAL release)
            install(EXPORT ${package}_Targets
                    FILE ${package}Targets.cmake
                    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)
        endif()

        install(FILES "${PROJECT_BINARY_DIR}/${package}Config.cmake"
                      "${PROJECT_BINARY_DIR}/${package}ConfigVersion.cmake"
                DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)


    endforeach(package)
    foreach(package ${INTERNAL_SUPPORT_LIBRARIES})
        install(FILES ${CPP_SRC_PATH}/${package}.cpp
                DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    endforeach(package)
endif()
