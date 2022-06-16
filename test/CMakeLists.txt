project(test)

enable_testing(true)

set(HDRS network.hpp
         ../include/airflownetwork/properties.hpp)

set(SRCS opening_tests.cpp
         powerlaw_tests.cpp
         two_way_flow_tests.cpp
         specified_flow_tests.cpp
)

#include_directories(../)
#include_directories(../dependencies/include)

add_executable(afn_tests afn_tests.cpp ${HDRS} ${SRCS})
add_test(afn_tests afn_tests)
#target_link_libraries(afn_tests libafn)
