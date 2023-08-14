#include_directories()
#link_directories()

ADD_EXECUTABLE( gco_classify
                gco_classify.cpp
                strtok.h
                utils.h
                utils.cpp
                gco/block.h
                gco/energy.h
                gco/GCoptimization.cpp
                gco/GCoptimization.h
                gco/graph.cpp
                gco/graph.h
                gco/LinkedBlockList.cpp
                gco/LinkedBlockList.h
                gco/maxflow.cpp
                )

TARGET_LINK_LIBRARIES(gco_classify ${EZMINC_LIBRARIES})

# ADD_EXECUTABLE( gco_segment
#                 gco_segment.cpp
#                 ../gco/graph.h
#                 )

#TARGET_LINK_LIBRARIES(gco_segment ${EZMINC_LIBRARIES})

ADD_EXECUTABLE( label_interaction_estimate 
                label_interaction_estimate.cpp)

TARGET_LINK_LIBRARIES(label_interaction_estimate ${EZMINC_LIBRARIES})


INSTALL(TARGETS 
  gco_classify
  label_interaction_estimate
  RUNTIME DESTINATION bin )


# kate: space-indent on; indent-width 2; mixedindent off; indent-mode cstyle;replace-tabs on;hl cmake
