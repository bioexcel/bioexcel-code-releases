# - Try to find IB Verbs
# Once done this will define
#  IB_VERBS_FOUND - System has IB Verbs
#  IB_VERBS_INCLUDE_DIRS - The IB Verbs include directories
#  IB_VERBS_LIBRARIES - The libraries needed to use IB Verbs

find_path(IB_VERBS_INCLUDE_DIR verbs.h
        HINTS /usr/local/include /usr/include/infiniband)
find_library(IB_VERBS_LIBRARY NAMES ibverbs
        PATHS /usr/local/lib /usr/lib /usr/lib64)
find_path(RDMA_CM_INCLUDE_DIR rdma_cma.h
        HINTS /usr/local/include /usr/include/rdma)
find_library(RDMA_CM_LIBRARY NAMES rdmacm
        PATHS /usr/local/lib /usr/lib /usr/lib64)

set(IB_VERBS_INCLUDE_DIRS ${IB_VERBS_INCLUDE_DIR})
set(IB_VERBS_LIBRARIES ${IB_VERBS_LIBRARY})
set(RDMA_CM_INCLUDE_DIRS ${RDMA_CM_INCLUDE_DIR})
set(RDMA_CM_LIBRARIES ${RDMA_CM_LIBRARY})
set(IB_INCLUDE_DIRS ${IB_VERBS_INCLUDE_DIR} ${IB_VERBS_LIBRARY} ${RDMA_CM_INCLUDE_DIR} ${RDMA_CM_LIBRARY})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set IB_VERBS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(IB_VERBS DEFAULT_MSG
        IB_VERBS_INCLUDE_DIR IB_VERBS_LIBRARY)
find_package_handle_standard_args(RDMA_CM DEFAULT_MSG
        RDMA_CM_INCLUDE_DIR RDMA_CM_LIBRARY)

set(IB_LIBRARIES ${IB_VERBS_LIBRARY} ${RDMA_CM_LIBRARY})

mark_as_advanced(IB_VERBS_INCLUDE_DIR IB_VERBS_LIBRARY)
mark_as_advanced(RDMA_CM_INCLUDE_DIR RDMA_CM_LIBRARY)
mark_as_advanced(IB_INCLUDE_DIRS IB_LIBRARIES)