cmake_minimum_required(VERSION 2.8.3)
project(alined)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -march=native -Wfatal-errors")
## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  roscpp
  cmake_modules
  )

## System dependencies are found with CMake's conventions
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(Eigen REQUIRED)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME}
  CATKIN_DEPENDS roscpp cmake_modules
  DEPENDS Ceres

)

file(GLOB ${PROJECT_NAME}_CODE_BASE
  "src/code_base/*.cpp"
)


include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
  )

add_executable(alined_test
  src/ros_nodes/alined_test_node.cpp
  ${${PROJECT_NAME}_CODE_BASE})
target_link_libraries(alined_test ${catkin_LIBRARIES})



add_library(${PROJECT_NAME} STATIC
  ${${PROJECT_NAME}_CODE_BASE}
)

target_link_libraries(${PROJECT_NAME}
   ${catkin_LIBRARIES}
   ${OpenCV_LIBRARIES}
)
