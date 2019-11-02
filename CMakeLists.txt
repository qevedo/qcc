cmake_minimum_required(VERSION 2.8)
project(qcc)

# Find files
file(GLOB HEADERS include/*.h)
file(GLOB SOURCES src/*.cpp)

# Include header files
include_directories(include)

# Create shared library
add_library(${PROJECT_NAME} SHARED ${SOURCES})

# Install library
#install(TARGETS ${PROJECT_NAME} DESTINATION lib/${PROJECT_NAME})

# Install library headers
#install(FILES ${HEADERS} DESTINATION include/${PROJECT_NAME})

# Tests
enable_testing()
add_subdirectory(tests)
