cmake_minimum_required(VERSION 3.17)
project(geftools)

set(CMAKE_CXX_STANDARD 11)

find_package(OpenCV REQUIRED COMPONENTS opencv_core opencv_highgui opencv_imgcodecs opencv_imgproc)
include_directories(${OpenCV_INCLUDE_DIRS})

find_package(HDF5 REQUIRED COMPONENTS C)
include_directories(${HDF5_INCLUDE_DIRS})

set(Headers cell_bin.h common_bin.h polygon.h cell_bin_writer.h  mask.h  utils.h)

aux_source_directory(. SRC_LIST)
add_library(objlib OBJECT ${SRC_LIST})
add_library(gef-static STATIC $<TARGET_OBJECTS:objlib>)
add_library(gef SHARED $<TARGET_OBJECTS:objlib>)
add_executable(geftools $<TARGET_OBJECTS:objlib>)
target_link_libraries(geftools ${HDF5_LIBRARIES} ${OpenCV_LIBRARIES})
target_link_libraries(gef ${HDF5_LIBRARIES} ${OpenCV_LIBRARIES})

set_target_properties(gef-static PROPERTIES OUTPUT_NAME "gef")
set_target_properties(gef PROPERTIES CLEAN_DIRECT_OUTPUT 1)
set_target_properties(gef-static PROPERTIES CLEAN_DIRECT_OUTPUT 1)


install(TARGETS gef
        LIBRARY DESTINATION
        ${CMAKE_INSTALL_PREFIX}/lib)

install(TARGETS gef-static
        ARCHIVE DESTINATION
        ${CMAKE_INSTALL_PREFIX}/lib)

install(TARGETS geftools
        RUNTIME DESTINATION
        ${CMAKE_INSTALL_PREFIX}/bin)

install(FILES ${Headers}
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include)


#-----------------------------------------------------------------------------
# Option to build documentation
#-----------------------------------------------------------------------------
set (GEFTOOLS_DOXYGEN_DIR  ${CMAKE_SOURCE_DIR}/docs)
option (GEFTOOLS_BUILD_DOC  "Build documentation" OFF)
if (GEFTOOLS_BUILD_DOC AND EXISTS "${GEFTOOLS_DOXYGEN_DIR}" AND IS_DIRECTORY "${GEFTOOLS_DOXYGEN_DIR}")
    # check if Doxygen is installed
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        message(STATUS "Doxygen version: ${DOXYGEN_VERSION}")
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile)

        # request to configure the file
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        add_custom_target( docs ALL
                COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating API documentation with Doxygen"
                VERBATIM )
    else ()
        message(STATUS "Doxygen needs to be installed to generate the doxygen documentation")
    endif ()
endif ()

