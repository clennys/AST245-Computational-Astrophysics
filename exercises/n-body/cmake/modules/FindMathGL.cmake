FIND_PATH(MathGL_INCLUDE_DIRS NAMES mgl2/mgl.h
  PATHS
  /opt/local/include
  /usr/include
  /usr/local/include
)

FIND_LIBRARY(MathGL_LIB NAMES mgl
  PATHS
  /opt/local/lib
  /usr/local/lib
  /usr/lib
  /usr/lib64
)
FIND_LIBRARY(MathGL_QT_LIB NAMES mgl-qt
  PATHS
  /opt/local/lib
  /usr/local/lib
  /usr/lib
  /usr/lib64
)


SET(MathGL_LIBRARIES ${MathGL_LIB} ${MathGL_QT_LIB})

IF (MathGL_INCLUDE_DIRS AND MathGL_LIBRARIES)
  SET(MathGL_FOUND TRUE)
  MESSAGE(STATUS "MathGL found")
  MESSAGE(STATUS "MathGL Include dirs: " ${MathGL_INCLUDE_DIRS})
  MESSAGE(STATUS "MathGL Library: " ${MathGL_LIB})
  MESSAGE(STATUS "MathGL-Qt Library: " ${MathGL_QT_LIB})


  # Add the MathGL library
  add_library(MathGL::mgl SHARED IMPORTED)
  set_target_properties(MathGL::mgl PROPERTIES
    IMPORTED_LOCATION ${MathGL_LIB}
    # INTERFACE_INCLUDE_DIRECTORIES ${MathGL_INCLUDE_DIR}
  )

  # Add the MathGL-Qt library
  add_library(MathGL::mgl-qt SHARED IMPORTED)
  set_target_properties(MathGL::mgl-qt PROPERTIES
    IMPORTED_LOCATION ${MathGL_QT_LIB}
    # INTERFACE_INCLUDE_DIRECTORIES ${MathGL_INCLUDE_DIR}
  )

ELSE (MathGL_INCLUDE_DIRS AND MathGL_LIBRARIES)
  MESSAGE(STATUS "MathGL was not found")
ENDIF(MathGL_INCLUDE_DIRS AND MathGL_LIBRARIES)
