  double precision :: DELTAX, DELTAY, DELTAT
  integer :: NT

  integer, allocatable, dimension (:,:) :: &
    MASK_CELL, MASK_NODE, MASK1X, MASK1Y, MASK2Y, MASK2X

  double precision, allocatable, dimension(:,:) :: &
    u0x, u0y, u1x, u1y, u2x, u2y, corr_ux, corr_uy

  double precision, allocatable :: C (:,:,:,:)

  double precision, allocatable,  dimension (:,:) :: rho, INV_MASS
  double precision, dimension (3,3) :: gx, gy

  double precision, allocatable :: wavelet(:)

  integer :: NX_source, NY_source

  integer :: itime, irec

  double precision, allocatable, dimension (:) :: X_rec, Y_rec
  integer, allocatable, dimension (:) :: NX_rec, NY_rec
  
  double precision,allocatable  :: inpo_rec(:,:,:)

  character(len=32) :: snap_file
  character(len=100),allocatable :: outfname_rec(:)

  double precision, allocatable :: buffer_ux(:,:), buffer_uy(:,:)

