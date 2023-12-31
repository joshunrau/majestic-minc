IF(GSL_FOUND AND HAVE_MINC1)
  ADD_EXECUTABLE(volume_pol
      volume_pol.cpp )

  TARGET_LINK_LIBRARIES( volume_pol 
    minc_image_proc ${EZMINC_LIBRARIES}  ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
    
  ADD_EXECUTABLE( random_volume
                  random_volume.cpp )

  TARGET_LINK_LIBRARIES(random_volume ${EZMINC_LIBRARIES}  ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})

  ADD_EXECUTABLE( fdr_threshold
                fdr_threshold.cpp )

  TARGET_LINK_LIBRARIES(fdr_threshold ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}  ${EZMINC_LIBRARIES} )
  
  ADD_EXECUTABLE(t2_fit
      t2_fit.cpp )

  TARGET_LINK_LIBRARIES(t2_fit ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EZMINC_LIBRARIES} )
  
  INSTALL(TARGETS 
        random_volume
        fdr_threshold
        volume_pol
        t2_fit
        RUNTIME DESTINATION bin )
        
  # PCA utils
  
  ADD_LIBRARY(pca_utils pca_utils.cpp  pca_utils.h  utils.cpp  utils.h )
  TARGET_LINK_LIBRARIES(pca_utils ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EZMINC_LIBRARIES} )
  
  ADD_EXECUTABLE(volumes_pca volumes_pca.cpp)
  TARGET_LINK_LIBRARIES(volumes_pca pca_utils)

  ADD_EXECUTABLE(grids_pca grids_pca.cpp)
  TARGET_LINK_LIBRARIES(grids_pca pca_utils)
  
  ADD_EXECUTABLE(volumes_pca_preprocess volumes_pca_preprocess.cpp)
  TARGET_LINK_LIBRARIES(volumes_pca_preprocess pca_utils)
  
  ADD_EXECUTABLE(volumes_lm volumes_lm.cpp)
  TARGET_LINK_LIBRARIES(volumes_lm pca_utils)
  
  ADD_EXECUTABLE(volumes_lsq volumes_lsq.cpp)
  TARGET_LINK_LIBRARIES(volumes_lsq pca_utils)
  
  ADD_EXECUTABLE(em_classify em_classify.cpp)
  TARGET_LINK_LIBRARIES(em_classify pca_utils)
  
  INSTALL(TARGETS 
        volumes_pca_preprocess
        volumes_pca
        volumes_lm
        volumes_lsq
        em_classify
        grids_pca
        RUNTIME DESTINATION bin )

ENDIF(GSL_FOUND AND HAVE_MINC1)


IF(FFTW3F_FOUND AND HAVE_MINC1)
  ADD_EXECUTABLE( fast_blur
      fast_blur.cpp )

  TARGET_LINK_LIBRARIES( fast_blur
    minc_image_proc ${EZMINC_LIBRARIES})
  
  INSTALL(TARGETS 
        fast_blur
        RUNTIME DESTINATION bin )
ENDIF(FFTW3F_FOUND AND HAVE_MINC1 )


IF(FFTW3F_FOUND AND GSL_FOUND AND HAVE_MINC1)
  ADD_EXECUTABLE( noise_estimate noise_estimate.cpp)

  TARGET_LINK_LIBRARIES(noise_estimate minc_image_proc ${EZMINC_LIBRARIES})

  ADD_EXECUTABLE( sharpness_estimate sharpness_estimate.cpp)

  TARGET_LINK_LIBRARIES(sharpness_estimate minc_image_proc ${EZMINC_LIBRARIES})
  
  INSTALL(TARGETS 
        noise_estimate sharpness_estimate
        RUNTIME DESTINATION bin )
ENDIF(FFTW3F_FOUND AND GSL_FOUND AND HAVE_MINC1)


IF(ITK_FOUND)
  IF(HAVE_MINC4ITK)
    LINK_DIRECTORIES(${MINC4ITK_LIBRARY_DIRS})
    INCLUDE_DIRECTORIES( ${MINC4ITK_INCLUDE_DIRS})  
  ELSE(HAVE_MINC4ITK)
    SET(ITK_LIBRARIES itk4MincHelpers ${ITK_LIBRARIES})
  ENDIF(HAVE_MINC4ITK)
 
  message("tools: ITK_LIBRARIES=${ITK_LIBRARIES}")

  ADD_EXECUTABLE( multiple_volume_similarity multiple_volume_similarity.cpp)
  TARGET_LINK_LIBRARIES( multiple_volume_similarity ${ITK_LIBRARIES}  )

  ADD_EXECUTABLE( itk_distance itk_distance.cpp)
  TARGET_LINK_LIBRARIES(itk_distance ${ITK_LIBRARIES})

  ADD_EXECUTABLE( itk_similarity itk_similarity.cpp )
  TARGET_LINK_LIBRARIES(itk_similarity ${ITK_LIBRARIES})

  ADD_EXECUTABLE( itk_resample itk_resample.cpp )
  TARGET_LINK_LIBRARIES(itk_resample  ${ITK_LIBRARIES})

  ADD_EXECUTABLE( itk_convert itk_convert.cpp )
  TARGET_LINK_LIBRARIES(itk_convert  ${ITK_LIBRARIES})

  ADD_EXECUTABLE( itk_convert_xfm itk_convert_xfm.cpp )
  TARGET_LINK_LIBRARIES(itk_convert_xfm  ${ITK_LIBRARIES})

  ADD_EXECUTABLE( itk_morph
                  itk_morph.cpp )
  TARGET_LINK_LIBRARIES( itk_morph   ${ITK_LIBRARIES})

  ADD_EXECUTABLE( itk_g_morph
                  itk_g_morph.cpp )
  TARGET_LINK_LIBRARIES( itk_g_morph   ${ITK_LIBRARIES} )

  ADD_EXECUTABLE( itk_laplacian_sharpening
                  itk_laplacian_sharpening.cpp )

  TARGET_LINK_LIBRARIES( itk_laplacian_sharpening  ${ITK_LIBRARIES})


  ADD_EXECUTABLE( itk_laplace
                  itk_laplace.cpp )

  TARGET_LINK_LIBRARIES( itk_laplace  ${ITK_LIBRARIES})

  ADD_EXECUTABLE( itk_vesselness
                  itk_vesselness.cpp
                  )

  TARGET_LINK_LIBRARIES( itk_vesselness  ${ITK_LIBRARIES})


  ADD_EXECUTABLE( itk_fill_hole
                  itk_fill_hole.cpp )

  TARGET_LINK_LIBRARIES( itk_fill_hole  ${ITK_LIBRARIES})


  ADD_EXECUTABLE( itk_diffusion
                  itk_diffusion.cpp )

  TARGET_LINK_LIBRARIES( itk_diffusion  ${ITK_LIBRARIES})

  ADD_EXECUTABLE( grid_proc grid_proc.cpp )
  TARGET_LINK_LIBRARIES(grid_proc ${ITK_LIBRARIES})

  INSTALL(TARGETS 
    multiple_volume_similarity 
    itk_distance 
    itk_similarity 
    itk_resample  
    itk_convert 
    itk_convert_xfm
    itk_morph 
    itk_g_morph
    itk_laplace
    itk_laplacian_sharpening
    itk_diffusion
    itk_vesselness
    itk_fill_hole
    grid_proc
  RUNTIME DESTINATION bin )


  IF(HAVE_MINC1)
    ADD_EXECUTABLE( resample_grid
                  resample_grid.cpp )

    TARGET_LINK_LIBRARIES(resample_grid ${ITK_LIBRARIES} ${VOLUME_IO_LIBRARIES} ${LIBMINC_LIBRARIES} )

    ADD_EXECUTABLE( dircos_to_xfm
                    dircos_to_xfm.cpp )
    TARGET_LINK_LIBRARIES( dircos_to_xfm itk4MincHelpers  ${ITK_LIBRARIES} )

    INSTALL(TARGETS dircos_to_xfm resample_grid RUNTIME DESTINATION bin)
  ENDIF(HAVE_MINC1)


ENDIF(ITK_FOUND)

IF(HAVE_MINC1)

  ADD_EXECUTABLE( volume_dwt volume_dwt.cpp)

  TARGET_LINK_LIBRARIES(volume_dwt minc_image_proc ${EZMINC_LIBRARIES})

  ADD_EXECUTABLE( minc_downsample 
                  minc_downsample.cpp )
  TARGET_LINK_LIBRARIES(minc_downsample ${EZMINC_LIBRARIES})

  ADD_EXECUTABLE( minc_nuyl
                  minc_nuyl.cpp )
  TARGET_LINK_LIBRARIES(minc_nuyl minc_image_proc ${EZMINC_LIBRARIES})

  ADD_EXECUTABLE( minc_rank
                  minc_rank.cpp )
  TARGET_LINK_LIBRARIES(minc_rank minc_image_proc ${EZMINC_LIBRARIES})

  ADD_EXECUTABLE( xfm2def xfm2def.c)
  TARGET_LINK_LIBRARIES(xfm2def ${VOLUME_IO_LIBRARIES} ${LIBMINC_LIBRARIES} )
  
  ADD_EXECUTABLE( tag2csv  tag2csv.cpp)
  TARGET_LINK_LIBRARIES(tag2csv ${ITK_LIBRARIES})

  ADD_EXECUTABLE( fuzzy_volume_similarity fuzzy_volume_similarity.cpp)
  TARGET_LINK_LIBRARIES(fuzzy_volume_similarity minc_image_proc ${EZMINC_LIBRARIES})
  
  ADD_EXECUTABLE( volume_gtc_similarity volume_gtc_similarity.cpp)
  TARGET_LINK_LIBRARIES(volume_gtc_similarity minc_image_proc ${EZMINC_LIBRARIES})
  
  ADD_EXECUTABLE( volume_similarity volume_similarity.cpp)
  TARGET_LINK_LIBRARIES(volume_similarity minc_image_proc ${EZMINC_LIBRARIES})
  
  ADD_EXECUTABLE( minc_taylor_reg minc_taylor_reg.c)
  TARGET_LINK_LIBRARIES(minc_taylor_reg ${VOLUME_IO_LIBRARIES} ${LIBMINC_LIBRARIES})


  INSTALL(TARGETS 
            volume_dwt
            minc_downsample
            minc_nuyl
            minc_rank
            xfm2def
            fuzzy_volume_similarity
            volume_gtc_similarity
            volume_similarity
            minc_taylor_reg
          RUNTIME DESTINATION bin )

ENDIF(HAVE_MINC1)
