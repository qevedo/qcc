project(qcc)

get_filename_component(PARENT_DIR ${PROJECT_SOURCE_DIR} DIRECTORY)

add_library(qcc
  platform.cc
  program.cc
  exception.cc
  arch/cc/eqasm_backend_cc.cc
  arch/cc/codegen_cc.cc
  arch/cc/vcd.cc
)

# support 'include <ql/foo.h>'
#target_include_directories(ql SYSTEM
#  PUBLIC ${PARENT_DIR}
#  )

set_target_properties(qcc PROPERTIES POSITION_INDEPENDENT_CODE ON)

unset(PARENT_DIR)
