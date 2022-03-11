set(QUICC_WORK_DIR
  "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks")

# Base dir for reference data archive
set(QUICC_REF_ARCHIVE_DIR "${CMAKE_BINARY_DIR}/Models/_refdata")

set(_impl "Explicit")
set(_benchmarkExe "BoussinesqShellRTC${_impl}Model")
set(_benchmarkTarget "Benchmark${_benchmarkExe}")

set(_refdir "${QUICC_WORK_DIR}/_refdata/${_impl}")
set(_rundir "${QUICC_WORK_DIR}/_data/${_impl}")
add_custom_target(${_benchmarkTarget} ALL
  COMMAND ${CMAKE_COMMAND} -E create_symlink
    "${_refdir}/parameters.cfg"
    "${_rundir}/parameters.cfg"
  COMMAND ${CMAKE_COMMAND} -E create_symlink
    "${_refdir}/state_initial.hdf5"
    "${_rundir}/state_initial.hdf5"
  COMMAND ${CMAKE_COMMAND} -E copy
    "${CMAKE_CURRENT_SOURCE_DIR}/validate_benchmark_${_impl}.py"
    "${_rundir}/validate_benchmark.py"
  COMMAND ${CMAKE_COMMAND} -E copy
    "${CMAKE_CURRENT_SOURCE_DIR}/colorcodes.py"
    "${_rundir}/"
  COMMAND ${CMAKE_COMMAND} -E copy
    "${CMAKE_CURRENT_SOURCE_DIR}/validation_tools.py"
    "${_rundir}/"
  )
add_dependencies(${_benchmarkTarget} ${_benchmarkExe})

# Fetch reference data
include(FetchBenchmarkReference)
quicc_fetch_benchmark_reference(
  ${_benchmarkTarget}
  MODEL BoussinesqShellRTC
  FILENAME "${_impl}.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)

set(_testname "Benchmark${_benchmarkExe}")
add_test(
  NAME ${_testname}
  COMMAND "${_benchmarkExe}"
  WORKING_DIRECTORY
  "${_rundir}"
  )

set(_testname "ValidateBenchmark${_benchmarkExe}")
add_test(
  NAME ${_testname}
  COMMAND "${Python_EXECUTABLE}" validate_benchmark.py
    -d "${_rundir}"
    -r "${_refdir}"
  WORKING_DIRECTORY
  "${_rundir}"
  )
set_tests_properties(ValidateBenchmark${_benchmarkExe} PROPERTIES
  PASS_REGULAR_EXPRESSION "All benchmark validation tests passed!"
  TIMEOUT 60
  )

set_tests_properties(ValidateBenchmark${_benchmarkExe} PROPERTIES DEPENDS "Benchmark${_benchmarkExe}")
