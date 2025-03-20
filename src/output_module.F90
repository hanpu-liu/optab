
MODULE output_module

  USE ISO_FORTRAN_ENV
  IMPLICIT NONE

CONTAINS  

  SUBROUTINE output_table(np, temp, rho, mean_p, mean_r, temp2, mean_p2)
    USE HDF5
    USE H5LT
    USE MPI
    USE h5LX_module, ONLY : h5LXset_dataset, h5LXwrite_hyperslab_dataset_F64
    USE mpi_module, ONLY : myrk
    REAL(REAL64), INTENT(IN) :: np(:,:)
    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: rho(:)
    REAL(REAL64), INTENT(IN) :: mean_p(:)
    REAL(REAL64), INTENT(IN) :: mean_r(:)
    REAL(REAL64), INTENT(IN) :: temp2
    REAL(REAL64), INTENT(IN) :: mean_p2(:)
    
    INTEGER(HID_T) :: file!, data_temp, data_rho, data_mean_p, data_mean_r, data_mean_p2, data_ndens
    INTEGER(HSIZE_T) :: dims(1), dims2(2)!, dims_global(1), dims2_global(2), start(1), count(1), start2(2), count2(2)
    INTEGER(INT32) :: error, nmax, jmax

    jmax = UBOUND(temp,1)
    
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, mean_p, jmax, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, error)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, mean_r, jmax, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, error)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, mean_p2, jmax, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, error)
    
    IF(myrk == 0) THEN
       CALL h5open_f(error)
       CALL h5Fcreate_f('output/opac_table.h5', H5F_ACC_TRUNC_F, file, error)
       dims = [1]
       CALL h5LTmake_dataset_f(file, 'jmax', 1, dims, H5T_STD_I32LE, jmax, error)
       CALL h5LTmake_dataset_f(file, 'temp2', 1, dims, H5T_IEEE_F64LE, temp2, error)
       
       dims = [jmax]
       CALL h5LTmake_dataset_f(file, 'temp', 1, dims, H5T_IEEE_F64LE, temp, error)
       CALL h5LTmake_dataset_f(file, 'rho', 1, dims, H5T_IEEE_F64LE, rho, error)
       CALL h5LTmake_dataset_f(file, 'mean_p', 1, dims, H5T_IEEE_F64LE, mean_p, error)
       CALL h5LTmake_dataset_f(file, 'mean_r', 1, dims, H5T_IEEE_F64LE, mean_r, error)
       CALL h5LTmake_dataset_f(file, 'mean_p2', 1, dims, H5T_IEEE_F64LE, mean_p2, error)
       nmax = UBOUND(np,1)
       dims2 = [nmax, jmax]
       CALL h5LTmake_dataset_f(file, 'ndens', 2, dims2, H5T_IEEE_F64LE, np(1:nmax,1:jmax), error)
       
       CALL h5Fclose_f(file, error)
       CALL h5close_f(error)
    END IF
    
    RETURN
  END SUBROUTINE output_table

  
  SUBROUTINE output_mono(j, temp, nden, grd, alp, sca, line, temp2, rho, pmean, pmean2, kbin, nbin)
    USE HDF5
    USE H5LT
    USE const_module, ONLY : c2, h, c, pi

    INTEGER(INT32), INTENT(IN) :: j
    REAL(REAL64), INTENT(IN) :: temp
    REAL(REAL64), INTENT(IN) :: nden(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(IN) :: alp(:)
    REAL(REAL64), INTENT(IN) :: sca(:)
    REAL(REAL64), INTENT(IN) :: line(:)
    REAL(REAL64), INTENT(IN) :: temp2
    REAL(REAL64), INTENT(IN) :: rho
    REAL(REAL64), INTENT(IN) :: pmean(:)
    REAL(REAL64), INTENT(IN) :: pmean2(:)
    INTEGER(INT64), INTENT(IN) :: kbin(:)
    INTEGER(INT32), INTENT(IN) :: nbin
    INTEGER(INT32) :: ns, ncurr
    INTEGER :: error
    INTEGER(HID_T) :: file_id!, plist_id, data_grid, data_abs, data_sca
    INTEGER(HSIZE_T), DIMENSION(1) :: dims1
    REAL(REAL64), ALLOCATABLE :: b(:), dbdt(:) ! Planck function
   !  REAL(REAL64) :: nume, deno, weight, pla, pla2, ros, plac, plal, numec, plac2, plal2
    REAL(REAL64) :: nume, deno, weight, numec, numee, numes, numegray, denogray, numecgray, numeegray, numesgray
    REAL(REAL64), ALLOCATABLE :: pla(:), pla2(:), ros(:), plac(:), plac2(:), plal(:), plal2(:), eff(:), scamean(:)
    INTEGER(INT64) :: k, ks, ke, k_total, nbinp1

    ! ------------------------------------------
    nbinp1 = nbin + 1
    ! last element is the gray value
    ALLOCATE(pla(nbinp1), pla2(nbinp1), ros(nbinp1), plac(nbinp1), plac2(nbinp1), plal(nbinp1), plal2(nbinp1), eff(nbinp1), scamean(nbinp1))

    ks = LBOUND(alp,1)    ! 1
    ke = UBOUND(alp,1)    ! 100001

    ! ------------------------------------------
    ! MEAN OPACITIES

    ALLOCATE(b(ks:ke), dbdt(ks:ke))
    
    ! PLANCK FUNCTION
    DO k = ks, ke
       b(k) = grd(k)**3 / (EXP(c2 * grd(k) / temp) - 1d0)
       dbdt(k) = grd(k)**4 / temp**2 * EXP(-c2 * grd(k) / temp) / (EXP(-c2 * grd(k) / temp) - 1d0)**2
    END DO
    b = b * (2d0 * h * c**2)
    dbdt = dbdt * (2d0 * h * c**2 * c2)

    ! PLANCK
    nume = 0d0
    numec= 0d0
    numee= 0d0
    numes= 0d0
    deno = 0d0
    numegray = 0d0
    numecgray= 0d0
    numeegray= 0d0
    numesgray= 0d0
    denogray = 0d0
    ncurr = 1
    DO k = ks, ke
       IF(k == ks) THEN
          weight = (grd(k+1) - grd(k  )) * 0.5d0
       ELSE IF(k == ke) THEN
          weight = (grd(k  ) - grd(k-1)) * 0.5d0
       ELSE
          weight = (grd(k+1) - grd(k-1)) * 0.5d0
       END IF
       ! gray
       numegray = numegray + b(k) * DBLE(alp(k) + line(k)) * weight
       numecgray= numecgray+ b(k) * DBLE(alp(k)) * weight
       numeegray= numeegray+ b(k) * DBLE(SQRT(3d0 * (alp(k) + line(k) + sca(k)) * (alp(k) + line(k)))) * weight
       numesgray= numesgray+ b(k) * DBLE(sca(k)) * weight
       denogray = denogray + b(k)                * weight
       ! multigroup
       IF(ncurr .GT. nbin) CYCLE    ! out of the maximum bin, multigroup finished
       IF(k .GT. kbin(ncurr+1) .OR. k == ke) THEN   ! one bin finished
         pla(ncurr)  = REAL(nume / deno)
         plac(ncurr) = REAL(numec/ deno)
         eff(ncurr) = REAL(numee/ deno)
         scamean(ncurr) = REAL(numes/ deno)
         plal(ncurr) = REAL(pmean(ncurr)/ deno * (2d0 * h * c**2))
         nume = 0d0
         numec= 0d0
         numee= 0d0
         numes= 0d0
         deno = 0d0
         ncurr = ncurr + 1
       ENDIF
       IF(k .GT. kbin(ncurr) .AND. k .LE. kbin(ncurr+1)) THEN
         nume = nume + b(k) * DBLE(alp(k) + line(k)) * weight
         numec= numec+ b(k) * DBLE(alp(k)) * weight
         numee= numee+ b(k) * DBLE(SQRT(3d0 * (alp(k) + line(k) + sca(k)) * (alp(k) + line(k)))) * weight
         numes= numes+ b(k) * DBLE(sca(k)) * weight
         deno = deno + b(k)                * weight
       END IF
    END DO
    pla(nbinp1)  = REAL(numegray / denogray)
    plac(nbinp1) = REAL(numecgray/ denogray)
    eff(nbinp1) = REAL(numeegray/ denogray)
    scamean(nbinp1) = REAL(numesgray/ denogray)
    plal(nbinp1) = REAL(pmean(nbinp1) / denogray * (2d0 * h * c**2))

    ! ROSSELAND
    nume = 0d0
    deno = 0d0
    numegray = 0d0
    denogray = 0d0
    ncurr = 1
    DO k = ks, ke
       IF(alp(k) + line(k) + sca(k) == 0d0) CYCLE
       IF(k == ks) THEN
          weight = (grd(k+1) - grd(k  )) * 0.5d0
       ELSE IF(k == ke) THEN
          weight = (grd(k  ) - grd(k-1)) * 0.5d0
       ELSE
          weight = (grd(k+1) - grd(k-1)) * 0.5d0
       END IF
       ! gray
       numegray = numegray + dbdt(k)                           * weight
       denogray = denogray + dbdt(k) / DBLE((alp(k) + line(k) +  sca(k))) * weight
       ! debug
       ! multigroup
       IF(ncurr .GT. nbin) CYCLE    ! out of the maximum bin, multigroup finished
       IF(k .GT. kbin(ncurr+1) .OR. k == ke) THEN   ! one bin finished
          ros(ncurr) = REAL(nume / deno)
          nume = 0d0
          deno = 0d0
          ncurr = ncurr + 1
       ENDIF
       IF(k .GT. kbin(ncurr) .AND. k .LE. kbin(ncurr+1)) THEN
          nume = nume + dbdt(k)                           * weight
          deno = deno + dbdt(k) / DBLE((alp(k) + line(k) +  sca(k))) * weight
       ENDIF
    END DO
    ros(nbinp1) = REAL(numegray / denogray)


    ! PLANCK FUNCTION with temp2
    DO k = ks, ke
       b(k) = grd(k)**3 / (EXP(c2 * grd(k) / temp2) - 1d0)
    END DO
    b = b * (2d0 * h * c**2)

    ! PLANCK
    nume = 0d0
    numec= 0d0
    deno = 0d0
    numegray = 0d0
    numecgray= 0d0
    denogray = 0d0
    ncurr = 1
    DO k = ks, ke
       IF(k == ks) THEN
          weight = (grd(k+1) - grd(k  )) * 0.5d0
       ELSE IF(k == ke) THEN
          weight = (grd(k  ) - grd(k-1)) * 0.5d0
       ELSE
          weight = (grd(k+1) - grd(k-1)) * 0.5d0
       END IF
       ! gray
       numegray = numegray + b(k) * DBLE(alp(k) + line(k)) * weight
       numecgray= numecgray+ b(k) * DBLE(alp(k)) * weight
       denogray = denogray + b(k)                * weight
       ! multigroup
       IF(ncurr .GT. nbin) CYCLE    ! out of the maximum bin, multigroup finished
       IF(k .GT. kbin(ncurr+1) .OR. k == ke) THEN   ! one bin finished
         pla2(ncurr) = REAL(nume / deno)
         plac2(ncurr)= REAL(numec/ deno)
         plal2(ncurr)= REAL(pmean2(ncurr)/ deno * (2d0 * h * c**2))
         nume = 0d0
         numec= 0d0
         deno = 0d0
         ncurr = ncurr + 1
       ENDIF
       IF(k .GT. kbin(ncurr) .AND. k .LE. kbin(ncurr+1)) THEN
         nume = nume + b(k) * DBLE(alp(k) + line(k)) * weight
         numec= numec+ b(k) * DBLE(alp(k)) * weight
         deno = deno + b(k)                * weight
       END IF
    END DO
    pla2(nbinp1) = REAL(numegray / denogray)
    plac2(nbinp1)= REAL(numecgray/ denogray)
    plal2(nbinp1) = REAL(pmean2(nbinp1) / denogray * (2d0 * h * c**2))

    DEALLOCATE(b, dbdt)
    ! ------------------------------------------

    k_total = ke - ks + 1
    ns = UBOUND(nden,1) - LBOUND(nden,1) + 1

    CALL h5open_f(error)
    CALL h5Fcreate_f('output/mono_'//i2c(j)//'.h5', H5F_ACC_TRUNC_F, file_id, error)
    
    dims1 = [1]
    CALL h5LTmake_dataset_f(file_id, 'temp', 1, dims1, H5T_IEEE_F64LE, temp, error)
    CALL h5LTmake_dataset_f(file_id, 'temp2', 1, dims1, H5T_IEEE_F64LE, temp2, error)
    CALL h5LTmake_dataset_f(file_id, 'rho', 1, dims1, H5T_IEEE_F64LE, rho, error)
    ! mean opacities
    ! grid-based Planck means
    dims1 = [nbinp1]
    CALL h5LTmake_dataset_f(file_id, 'pla', 1, dims1, H5T_IEEE_F64LE, pla, error)
    CALL h5LTmake_dataset_f(file_id, 'pla2', 1, dims1, H5T_IEEE_F64LE, pla2, error)
    ! line-based Planck means and continuum Planck means
    CALL h5LTmake_dataset_f(file_id, 'plac', 1, dims1, H5T_IEEE_F64LE, plac, error)
    CALL h5LTmake_dataset_f(file_id, 'plac2', 1, dims1, H5T_IEEE_F64LE, plac2, error)
    CALL h5LTmake_dataset_f(file_id, 'plal', 1, dims1, H5T_IEEE_F64LE, plal, error)
    CALL h5LTmake_dataset_f(file_id, 'plal2', 1, dims1, H5T_IEEE_F64LE, plal2, error)
    ! Planck-like scattering mean
    CALL h5LTmake_dataset_f(file_id, 'scamean', 1, dims1, H5T_IEEE_F64LE, scamean, error)
    CALL h5LTmake_dataset_f(file_id, 'eff', 1, dims1, H5T_IEEE_F64LE, eff, error)
    ! Rosseland means
    CALL h5LTmake_dataset_f(file_id, 'ros', 1, dims1, H5T_IEEE_F64LE, ros, error)
    dims1 = [ns]
    CALL h5LTmake_dataset_f(file_id, 'nden', 1, dims1, H5T_IEEE_F64LE, nden, error)
    ! monochromatic opacities
    dims1 = [k_total]
    CALL h5LTmake_dataset_f(file_id, 'grd', 1, dims1, H5T_IEEE_F64LE, grd, error)
    CALL h5LTmake_dataset_f(file_id, 'abs', 1, dims1, H5T_IEEE_F64LE, alp, error)
    CALL h5LTmake_dataset_f(file_id, 'sca', 1, dims1, H5T_IEEE_F64LE, sca, error)
    CALL h5LTmake_dataset_f(file_id, 'line', 1, dims1, H5T_IEEE_F64LE, line, error)

    CALL h5Fclose_f(file_id, error)
    CALL h5close_f(error)

    DEALLOCATE(pla, pla2, ros, plac, plac2, plal, plal2, eff, scamean)

    RETURN
  END SUBROUTINE output_mono


  FUNCTION i2c(j)
    IMPLICIT NONE
    CHARACTER*5 :: i2c
    INTEGER(INT32) :: j
    WRITE(i2c,FMT='(I5.5)') j
    RETURN
  END FUNCTION i2c

END MODULE output_module
