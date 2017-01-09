include_directories(src)

include ( ${DD4hep_DIR}/cmake/DD4hep.cmake )

dd4hep_configure_output()
dd4hep_package(   Bach MAJOR 0 MINOR 0 PATCH 1
  USES            [ROOT   REQUIRED COMPONENTS Geom GenVector MathCore]
                  [DD4hep REQUIRED COMPONENTS DDCore DDCond DDAlign]
  LINK_LIBRARIES  TbAlgorithms TbKernel
  OPTIONAL        XERCESC
  )

dd4hep_add_plugin( AIDABach SOURCES src/bach.cpp )
dd4hep_install_dir( compact DESTINATION examples/Bach )
dd4hep_configure_scripts( Bach DEFAULT_SETUP)
