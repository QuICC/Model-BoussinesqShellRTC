include(ConfigureBenchmarks)
quicc_add_benchmark(Explicit
  MODEL "BoussinesqShellRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:serial"
  MPIRANKS 6
  GITTAG "rratio"
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqShellRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:single1d"
  MPIRANKS 6
  GITTAG "rratio"
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqShellRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:single2d"
  MPIRANKS 6
  GITTAG "rratio"
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqShellRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:tubular"
  MPIRANKS 6
  GITTAG "rratio"
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqShellRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:serial"
  MPIRANKS 6
  GITTAG "rratio"
  )

quicc_add_benchmark(Implicit
  MODEL "BoussinesqShellRTC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:single2d"
  MPIRANKS 6
  GITTAG "rratio"
  )
