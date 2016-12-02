include_directories(.)

add_library(Millepede SHARED Millepede.cpp)
target_link_libraries(Millepede ${ROOT_LIBRARIES})
target_include_directories(Millepede PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/.>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)

