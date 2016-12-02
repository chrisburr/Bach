include_directories(src)

file(GLOB SOURCES "src/*.cpp")
add_library(TbKernel SHARED ${SOURCES})
target_link_libraries(TbKernel ${ROOT_LIBRARIES})
target_include_directories(TbKernel PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)
