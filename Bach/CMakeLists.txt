include_directories(src)

add_executable(Bach src/bach.cpp)
target_link_libraries(Bach TbAlgorithms TbKernel ${ROOT_LIBRARIES})
