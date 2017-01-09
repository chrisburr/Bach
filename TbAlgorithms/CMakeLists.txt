include_directories(src ${DD4hep_INCLUDE_DIRS})

file(GLOB SOURCES "src/*.cpp")
add_library(TbAlgorithms SHARED ${SOURCES})
target_link_libraries(TbAlgorithms Millepede TbKernel ${ROOT_LIBRARIES})
target_include_directories(TbAlgorithms PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)
