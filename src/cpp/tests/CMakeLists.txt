set(TEST_NAME "test_${PROJECT_NAME}")
add_executable(${TEST_NAME} "${TEST_NAME}.cpp")
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
target_compile_options(${TEST_NAME} PRIVATE)
target_link_libraries(
    ${TEST_NAME}
    PUBLIC ${project_link_string} ${PROJECT_LINK_LIBRARIES}
)

# Local builds of upstream projects require local include paths
if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(
        ${TEST_NAME}
        PUBLIC
            ${Boost_INCLUDE_DIRS} # Required for MacOSX CMake builds. Not sure if it's a CMake or clang issue.
            ${tardigrade_vector_tools_SOURCE_DIR}/src/cpp
            ${tardigrade_error_tools_SOURCE_DIR}/src/cpp
            ${tardigrade_constitutive_tools_SOURCE_DIR}/src/cpp
    )
endif()

get_target_property(target_sources ${TEST_NAME} SOURCES)
foreach(source ${target_sources})
    set_property(
        GLOBAL
        APPEND
        PROPERTY
            CLANG_FORMAT_SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${source}"
    )
endforeach(source)

foreach(support_module ${INTERNAL_SUPPORT_LIBRARIES})
    set(TEST_NAME "test_${support_module}")
    add_executable(${TEST_NAME} "${TEST_NAME}.cpp")
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    target_link_libraries(${TEST_NAME} PUBLIC ${PROJECT_NAME})

    # Local builds of upstream projects require local include paths
    if(NOT cmake_build_type_lower STREQUAL "release")
        target_include_directories(
            ${TEST_NAME}
            PUBLIC
                ${Boost_INCLUDE_DIRS} # Required for MacOSX CMake builds. Not sure if it's a CMake or clang issue.
                ${tardigrade_vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                ${tardigrade_error_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                ${tardigrade_constitutive_tools_SOURCE_DIR}/${CPP_SRC_PATH}
        )
    endif()

    get_target_property(target_sources ${TEST_NAME} SOURCES)
    foreach(source ${target_sources})
        set_property(
            GLOBAL
            APPEND
            PROPERTY
                CLANG_FORMAT_SOURCE_FILES
                    "${CMAKE_CURRENT_SOURCE_DIR}/${source}"
        )
    endforeach(source)
endforeach(support_module)

get_target_property(target_sources ${TEST_NAME} SOURCES)
foreach(source ${target_sources})
    set_property(
        GLOBAL
        APPEND
        PROPERTY
            CLANG_FORMAT_SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${source}"
    )
endforeach(source)

foreach(support_module ${ADDITIONAL_HEADER_ONLY_LIBRARIES})
    set(TEST_NAME "test_${support_module}")
    add_executable(${TEST_NAME} "${TEST_NAME}.cpp")
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    target_link_libraries(${TEST_NAME} PUBLIC ${PROJECT_NAME})

    # Local builds of upstream projects require local include paths
    if(NOT cmake_build_type_lower STREQUAL "release")
        target_include_directories(
            ${TEST_NAME}
            PUBLIC
                ${Boost_INCLUDE_DIRS} # Required for MacOSX CMake builds. Not sure if it's a CMake or clang issue.
                ${tardigrade_vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                ${tardigrade_error_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                ${tardigrade_constitutive_tools_SOURCE_DIR}/${CPP_SRC_PATH}
        )
    endif()

    get_target_property(target_sources ${TEST_NAME} SOURCES)
    foreach(source ${target_sources})
        set_property(
            GLOBAL
            APPEND
            PROPERTY
                CLANG_FORMAT_SOURCE_FILES
                    "${CMAKE_CURRENT_SOURCE_DIR}/${source}"
        )
    endforeach(source)
endforeach(support_module)
