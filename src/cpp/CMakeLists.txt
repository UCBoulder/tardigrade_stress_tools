add_library (${PROJECT_NAME} SHARED "${PROJECT_NAME}.cpp" "${PROJECT_NAME}.h")
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11
                                                 PUBLIC_HEADER ${PROJECT_NAME}.h)
target_link_libraries(${PROJECT_NAME} error_tools constitutive_tools)
target_compile_options(${PROJECT_NAME} PUBLIC)

# Support libraries
foreach(support_library ${INTERNAL_SUPPORT_LIBRARIES})
    string(TOUPPER ${support_library} library_uppercase)
    add_library (${support_library} SHARED "${support_library}.cpp" "${support_library}.h")
    set_target_properties(${support_library} PROPERTIES CXX_STANDARD 11
                                                        PUBLIC_HEADER ${support_library}.h)
    target_link_libraries(${support_library} error_tools "${${library_uppercase}_UPSTREAM_LIBRARIES}")
    target_compile_options(${support_library} PUBLIC)
endforeach(support_library)

# Local builds of upstream projects require local include paths
if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(${PROJECT_NAME} PUBLIC
                               "${vector_tools_SOURCE_DIR}/src/cpp"
                               "${error_tools_SOURCE_DIR}/src/cpp"
                               "${constitutive_tools_SOURCE_DIR}/src/cpp")
endif()

install(TARGETS ${PROJECT_NAME} ${INTERNAL_SUPPORT_LIBRARIES}
        EXPORT ${PROJECT_NAME}_Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
