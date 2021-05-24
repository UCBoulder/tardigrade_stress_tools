set(TEST_NAME "test_${PROJECT_NAME}")
add_executable(${TEST_NAME} "${TEST_NAME}.cpp")
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
set_target_properties(${TEST_NAME} PROPERTIES CXX_STANDARD 11)
target_link_libraries(${TEST_NAME} PUBLIC ${PROJECT_NAME} error_tools constitutive_tools)

# Local builds of upstream projects require local include paths
if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(${PROJECT_NAME} PUBLIC
                               ${Boost_INCLUDE_DIRS}
                               ${vector_tools_SOURCE_DIR}/src/cpp
                               ${error_tools_SOURCE_DIR}/src/cpp
                               ${constitutive_SOURCE_DIR}/src/cpp)
endif()
