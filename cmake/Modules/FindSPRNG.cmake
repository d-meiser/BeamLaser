# Use SPRNG_ROOT_DIR to specify location of cgreen
# When successful this defines
# SPRNG_FOUND
# SPRNG_LIB
# SPRNG_INCLUDE_DIR

set(SPRNG_GUESS_DIR ${PROJECT_SOURCE_DIR}/sprng2.0)
find_path(SPRNG_INCLUDE_DIR sprng.h
    HINTS ${SPRNG_ROOT_DIR}/include;${SPRNG_GUESS_DIR}/include;/usr/include;/usr/local/include)
find_library(SPRNG_LIB sprng
    HINTS ${SPRNG_ROOT_DIR}/lib;${SPRNG_GUESS_DIR}/lib;/usr/local/lib;/usr/lib;/usr/lib64)
if (SPRNG_INCLUDE_DIR AND SPRNG_LIB)
  set(SPRNG_FOUND TRUE)
endif ()

