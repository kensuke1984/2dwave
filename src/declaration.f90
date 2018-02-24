  double precision :: DELTAX, DELTAY, DELTAT
  integer :: NT

  integer, allocatable, dimension (:,:) :: &
    MASK_CELL, MASK_NODE, MASK1X, MASK1Y, MASK2Y, MASK2X

  double precision, allocatable, dimension(:,:) :: &
    u0x, u0y, u1x, u1y, u2x, u2y, corr_ux, corr_uy

  double precision, allocatable :: C (:,:,:,:)

  double precision,allocatable,  dimension (:,:) :: rho, INV_MASS
  double precision, dimension (3,3) :: gx, gy

  double precision, allocatable :: wavelet(:)

  integer :: NX_source, NY_source

  integer :: itime, irec

  double precision,allocatable,  dimension (:) :: X_rec, Y_rec
  integer, allocatable, dimension (:) :: NX_rec, NY_rec
  
  double precision,allocatable  :: inpo_rec(:,:,:)

  character(len=32) :: snap_file
  character(len=32),allocatable :: outfname_rec(:)

  double precision, allocatable :: buffer_ux(:,:), buffer_uy(:,:)

!----------------valuables for CPML-------------------------------------

  double precision, allocatable, dimension(:) :: &
    a1_lef, a1_rig, a1_btm, b1_lef, b1_rig, b1_btm, &
    a2_lef, a2_rig, a2_btm, b2_lef, b2_rig, b2_btm

  double precision, allocatable, dimension(:,:) :: &
    dux_dx_lef, duy_dx_lef, dux_dx_rig, duy_dx_rig, dux_dy_btm, duy_dy_btm, &
    dTxx_dx_lef, dTyx_dx_lef, dTxx_dx_rig, dTyx_dx_rig, dTxy_dy_btm, dTyy_dy_btm

  double precision, allocatable, dimension(:,:) :: &
    mem1_xx_lef, mem1_xy_lef, mem1_xx_rig, mem1_xy_rig, mem1_yx_btm, mem1_yy_btm, &
    mem2_xx_lef, mem2_xy_lef, mem2_xx_rig, mem2_xy_rig, mem2_yx_btm, mem2_yy_btm

  double precision, allocatable, dimension(:,:) :: &
    corr_dux_dx_lef, corr_duy_dx_lef, corr_dux_dx_rig, corr_duy_dx_rig, &
    corr_dux_dy_btm, corr_duy_dy_btm, corr_dTxx_dx_lef, corr_dTyx_dx_lef, &
    corr_dTxx_dx_rig, corr_dTyx_dx_rig, corr_dTxy_dy_btm, corr_dTyy_dy_btm

  double precision,  allocatable, dimension(:,:):: &
    corr_mem1_xx_lef, corr_mem1_xy_lef, corr_mem1_xx_rig, &
    corr_mem1_xy_rig, corr_mem1_yx_btm, corr_mem1_yy_btm, &
    corr_mem2_xx_lef, corr_mem2_xy_lef, corr_mem2_xx_rig, &
    corr_mem2_xy_rig, corr_mem2_yx_btm, corr_mem2_yy_btm

!-----------END CPML--------------------------------------------------------

 
