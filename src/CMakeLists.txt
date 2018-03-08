ADD_EXECUTABLE(cortex_walk cortex_walk.cpp minc_wrappers.cpp)
TARGET_LINK_LIBRARIES(cortex_walk bicpl ${ITK_LIBRARIES}  ${VOLUME_IO_LIBRARIES} ${LIBMINC_LIBRARIES})

