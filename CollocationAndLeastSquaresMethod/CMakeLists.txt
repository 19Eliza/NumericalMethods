cmake_minimum_required (VERSION 3.8)

project("Eigen" VERSION 0.0.1 LANGUAGES CXX)
find_package(Eigen3 CONFIG REQUIRED)

file(GLOB PROJECT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_executable(${PROJECT_NAME})

target_sources(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCES})

set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)