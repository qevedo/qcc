add_executable(test_basic basic.cpp)
target_link_libraries(test_basic ${PROJECT_NAME})

add_test(basic test_basic)
