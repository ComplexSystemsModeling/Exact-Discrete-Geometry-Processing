cmake_minimum_required(VERSION 2.8)

project( PointInCircle )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

add_executable( IsPointInCircle 
  IsPointInCircle.cxx
  )
target_link_libraries( IsPointInCircle
  ${ITK_LIBRARIES}
  )

enable_testing()

add_test( TestIsPointInCircle_1
  IsPointInCircle -0.2 0.2
  )
set_tests_properties(
    TestIsPointInCircle_1 
    PROPERTIES WILL_FAIL "TRUE"
    )

add_test( TestIsPointInCircle_2
  IsPointInCircle 0.2 -0.2
  )
