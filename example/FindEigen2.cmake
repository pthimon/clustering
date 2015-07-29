# - Try to find Eigen2
# Once done, this will define
#
#  Eigen2_FOUND - system has Eigen2
#  Eigen2_INCLUDE_DIR - the Eigen2 include directory

# Include dir
find_path(Eigen2_INCLUDE_DIR
  NAMES Eigen/Core
  PATH_SUFFIXES eigen2
)

set(Eigen2_FOUND "NO")
if(Eigen2_INCLUDE_DIR)
  set(Eigen2_FOUND "YES")
endif(Eigen2_INCLUDE_DIR)

