cmake_minimum_required(VERSION 2.6)

project(ANN CXX)
set(major 1)
set(minor 1)
set(revision 2)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${ANN_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ANN_SOURCE_DIR}/bin)

# the ANN library
include_directories(include)
include_directories(src)
file(GLOB sources src/*.cpp)
add_library(ANN SHARED ${sources})
set_target_properties(ANN PROPERTIES VERSION ${major}.${minor}.${revision}
	SOVERSION ${major}
	COMPILE_FLAGS 
	"-O3")

# ANN test utility
include_directories(test)
file(GLOB test_sources test/*.cpp)
add_executable(ann_test ${test_sources})
target_link_libraries(ann_test ANN)

# ANN sample utility
add_executable(ann_sample sample/ann_sample.cpp)
target_link_libraries(ann_sample ANN)

# ann2fig utility
add_executable(ann2fig ann2fig/ann2fig.cpp)
target_link_libraries(ann2fig ANN)

# install
install(TARGETS ANN ann_test ann_sample ann2fig 
        RUNTIME DESTINATION bin
              LIBRARY DESTINATION lib)
install(DIRECTORY sample doc DESTINATION share/ANN)
