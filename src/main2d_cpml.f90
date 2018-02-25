! ============================================================================
! Name        : main2d_cpml.f90
! Author      : Kei Hasegawa
! Version     : 0.0.1
! Copyright   : It is Complicated.
! Description : Wave propagation in 2D space.
! ============================================================================
program main2d

  use parameter_2dwave
  use copt2d_cpml
  implicit none

!---------------------------------------------------------------
  include 'declaration.f90'
  call get_arguments
  call read_parameter
  call allocate_main
  call allocate_cpml(NX, NY, NX_CPML_lef,NX_CPML_rig, NY_CPML_btm)

!###########################################################################
  DELTAX = X_size / dble(NX)
  DELTAY = Y_size / dble(NY)
  DELTAT = CFL * min(DELTAX,DELTAY) / maxvp
  NT = int( tlen / DELTAT)

  print *, 'DELTAX = ', DELTAX
  print *, 'DELTAY = ', DELTAY
  print *, 'DELTAT = ', DELTAT

  print *, 'NX = ', NX
  print *, 'NY = ', NY
  print *, 'NT = ', NT
  call output_parameter


  call set_physical_parameters&
    (NX, NY, rho, C, modelfile, header_model)

  if (need_waveform) then

    call set_receiver_info&
      ( num_rec, X_rec, Y_rec, outfname_rec, recinfo )

    allocate( buffer_ux(num_rec,NT), buffer_uy(num_rec,NT) )
    allocate (NX_rec(num_rec), NY_rec(num_rec))
    allocate(inpo_rec(3,3,num_rec))
    call set_reciever_interpolating_function&
      ( NX,NY,num_rec,DELTAX,DELTAY,X_rec,Y_rec,NX_rec,NY_rec,inpo_rec )
  endif

  call make_mask&
    ( NX,NY,rho,MASK_CELL,MASK_NODE,MASK1X,MASK1Y,MASK2X,MASK2Y )

  call compt_inverse_diagonal_mass&
    ( NX,NY,MASK1Y,MASK2X,rho,DELTAT,INV_MASS )

  allocate( wavelet(NT) )
  call calc_Ricker_wavelet&
    ( NT, DELTAT, tp, ts, wavelet )


  if ( source_type_single ) then

    call set_source_interpolating_function_single&
      ( NX,NY,DELTAX,DELTAY,X_source,Y_source,vx,vy,NX_source,NY_source,gx,gy )

  elseif ( source_type_double_couple ) then

    call set_source_interpolating_function_double_couple&
      ( NX,NY,DELTAX,DELTAY,X_source,Y_source,mxx,mxy,myx,myy,NX_source,NY_source,gx,gy )

  else

    stop 'invalid source type'

  endif

  call set_CPML_parameters&
    ( NX,NY,DELTAT,DELTAX,DELTAY, &
      NX_CPML_lef, NX_CPML_rig, NY_CPML_btm, &
      maxvp, reflection_rate, source_freq, NPOWER, &
      a1_lef, a1_rig, a1_btm, &
      b1_lef, b1_rig, b1_btm, &
      a2_lef, a2_rig, a2_btm, &
      b2_lef, b2_rig, b2_btm )

!#######################compute wavefield##########################################
  do itime = 1,NT

!-----------write snapshot-------------------------------------------
    if ( need_snapshot .and. mod(itime, NT_snap) == 0 ) then

      write(*,*) 'itime = ',itime
      write(snap_file,*) itime

      snap_file = 'field_' // ( trim( adjustl(snap_file) ) ) // '.dat'

      call write_snapshot( NX,NY,MASK_NODE,DELTAX,DELTAY,u1x,u1y,snap_file )

    endif
!-----END------------------------------------------------------------

!---save waveforms to buffer arrays----------------------------------
    if ( need_waveform ) then

      do irec = 1,num_rec

        buffer_ux(irec,itime) = sum( &
          u1x( NX_rec(irec)-2:NX_rec(irec), NY_rec(irec)-2:NY_rec(irec) ) &
          * inpo_rec(:,:,irec) )

        buffer_uy(irec,itime) = sum( &
          u1y( NX_rec(irec)-2:NX_rec(irec), NY_rec(irec)-2:NY_rec(irec) ) &
          * inpo_rec(:,:,irec) )

      enddo

    endif
!-------END--------------------------------------------------------

    call compt_traction_force&
      ( NX,NY,MASK1X,MASK1Y,MASK2X,MASK2Y,u1x,u1y,C,DELTAX,DELTAY, &
        NX_CPML_lef, NX_CPML_rig, NY_CPML_btm, &
        a1_lef, a1_rig, a1_btm, &
        b1_lef, b1_rig, b1_btm, &
        a2_lef, a2_rig, a2_btm, &
        b2_lef, b2_rig, b2_btm, &
        dux_dx_lef, duy_dx_lef, &
        dux_dx_rig, duy_dx_rig, &
        dux_dy_btm, duy_dy_btm, &
        dTxx_dx_lef, dTyx_dx_lef, &
        dTxx_dx_rig, dTyx_dx_rig, &
        dTxy_dy_btm, dTyy_dy_btm, &
        mem1_xx_lef, mem1_xy_lef, &
        mem1_xx_rig, mem1_xy_rig, &
        mem1_yx_btm, mem1_yy_btm, &
        mem2_xx_lef, mem2_xy_lef, &
        mem2_xx_rig, mem2_xy_rig, &
        mem2_yx_btm, mem2_yy_btm, &
        u2x,u2y )

    u2x(NX_source-2:NX_source,NY_source-2:NY_source) &
      = u2x(NX_source-2:NX_source,NY_source-2:NY_source) - gx * wavelet(itime)
    u2y(NX_source-2:NX_source,NY_source-2:NY_source) &
      = u2y(NX_source-2:NX_source,NY_source-2:NY_source) - gy * wavelet(itime)

    u2x = - INV_MASS * u2x
    u2y = - INV_MASS * u2y

!----Dirichlet boundary contion----------------------
    u2x(0,:) = 0.d0
    u2y(0,:) = 0.d0

    u2x(NX,:) = 0.d0
    u2y(NX,:) = 0.d0

    u2x(:,0) = 0.d0
    u2y(:,0) = 0.d0
!-------END------------------------------------------

    call compt_corrector_force&
      ( NX,NY,MASK1X,MASK1Y,MASK2X,MASK2Y,u2x,u2y,C,rho,DELTAT,DELTAX,DELTAY, &
        NX_CPML_lef, NX_CPML_rig, NY_CPML_btm, &
        a1_lef, a1_rig, a1_btm, &
        b1_lef, b1_rig, b1_btm, &
        a2_lef, a2_rig, a2_btm, &
        b2_lef, b2_rig, b2_btm, &
        corr_dux_dx_lef, corr_duy_dx_lef, &
        corr_dux_dx_rig, corr_duy_dx_rig, &
        corr_dux_dy_btm, corr_duy_dy_btm, &
        corr_dTxx_dx_lef, corr_dTyx_dx_lef, &
        corr_dTxx_dx_rig, corr_dTyx_dx_rig, &
        corr_dTxy_dy_btm, corr_dTyy_dy_btm, &
        corr_mem1_xx_lef, corr_mem1_xy_lef, &
        corr_mem1_xx_rig, corr_mem1_xy_rig, &
        corr_mem1_yx_btm, corr_mem1_yy_btm, &
        corr_mem2_xx_lef, corr_mem2_xy_lef, &
        corr_mem2_xx_rig, corr_mem2_xy_rig, &
        corr_mem2_yx_btm, corr_mem2_yy_btm, &
        corr_ux,corr_uy )

    corr_ux = - INV_MASS * corr_ux
    corr_uy = - INV_MASS * corr_uy

    u2x = u2x + corr_ux
    u2y = u2y + corr_uy

!----Dirichlet boundary contion----------------------
    u2x(0,:) = 0.d0
    u2y(0,:) = 0.d0

    u2x(NX,:) = 0.d0
    u2y(NX,:) = 0.d0

    u2x(:,0) = 0.d0
    u2y(:,0) = 0.d0
!------END-------------------------------------------


!---Carjan abrorbing boundary condition (optional)---------
    if ( use_Cerjan ) then

      call compt_absorbing_boundary_Cerjan&
        ( NX,NY,u0x,u0y,u1x,u1y,u2x,u2y,N_Cerjan,dumping_rate_Cerjan )

    endif
!----END------------------------------------------------------

    u2x = u2x + 2.d0 * u1x - u0x
    u2y = u2y + 2.d0 * u1y - u0y

    u0x = u1x
    u0y = u1y
    u1x = u2x
    u1y = u2y

  enddo !end of time loop

!------output waveforms-------------------------------------------

  if ( need_waveform ) then

    do irec = 1,num_rec
      open(17,file=outfname_rec(irec),status='replace')

      write (17,*) 0.d0, 0.d0, 0.d0

      do itime = 1,NT

        write(17,fmt=wave_format) DELTAT*dble(itime), buffer_ux(irec,itime), buffer_uy(irec,itime)

      enddo

      close(17)

    enddo

  endif

!----------------------------------------------------------------

contains

!##################################################################################
  subroutine allocate_main
    allocate(MASK_CELL(NX,NY), MASK_NODE(0:NX,0:NY))
    allocate(MASK1X(NX,0:NY), MASK1Y(0:NX,NY))
    allocate(MASK2X(NX,NY), MASK2Y(NX,NY))
    allocate(u0x(0:NX,0:NY), u0y(0:NX,0:NY))
    allocate(u1x(0:NX,0:NY), u1y(0:NX,0:NY))
    allocate(u2x(0:NX,0:NY), u2y(0:NX,0:NY))
    allocate(corr_ux(0:NX,0:NY), corr_uy(0:NX,0:NY))
    allocate(C(3,3,NX,NY), rho(NX,NY))
    allocate(INV_MASS(0:NX,0:NY))
!    allocate(X_rec(num_rec), Y_rec(num_rec))
!    allocate (NX_rec(num_rec), NY_rec(num_rec))
!    allocate(inpo_rec(3,3,num_rec))
!    allocate(outfname_rec(num_rec))

    u0x = 0.d0
    u0y = 0.d0
    u1x = 0.d0
    u1y = 0.d0

  end subroutine allocate_main




  subroutine write_snapshot( NX,NY,MASK_NODE,DELTAX,DELTAY,ux,uy,fname )

    integer, intent(in) :: &
      NX, NY, MASK_NODE(0:NX,0:NY)

    double precision, intent(in) :: &
      DELTAX, DELTAY, &
      ux(0:NX,0:NY), uy(0:NX,0:NY)

    character(*), intent(in) :: &
      fname

    double precision :: X, Y
    integer :: i, j

    open(17,file=fname,status='replace')

    do i = 0,NX

      X = dble(i) * DELTAX

      do j = 0,NY

        Y = dble(j) * DELTAY

        if ( MASK_NODE(i,j) == 1) then

          write(17,fmt=snap_format) X, Y, ux(i,j), uy(i,j)

        endif

      enddo

      write(17,*) ''

    enddo

    close(17)

  end subroutine

!##################################################################################

  subroutine  calc_Ricker_wavelet( NT,DELTAT,tp,ts,wavelet )

    integer, intent(in) :: &
      NT
    double precision, intent(in) :: &
      DELTAT, tp, ts

    double precision, intent(out) :: &
      wavelet(NT)

    double precision, parameter :: PI = 3.1415926535897932d0
    double precision :: t, b
    integer :: i

    do i = 1,NT

      t = dble(i) * DELTAT

      b = PI * (t-ts) / tp

      wavelet(i) = 0.5d0 * sqrt(PI) * (b * b - 0.5d0) * exp( -b * b )

    enddo

  end subroutine

!##################################################################################

  subroutine set_source_interpolating_function_single&
    ( NX,NY,DELTAX,DELTAY,X_source,Y_source,vx,vy,NX_source,NY_source,gx,gy )

    integer, intent(in) :: &
      NX, NY
    double precision, intent(in) :: &
      DELTAX, DELTAY, X_source, Y_source, vx, vy

    integer, intent(out) :: &
      NX_source, NY_source

    double precision, intent(out) :: &
      gx(3,3), gy(3,3)

    integer :: i,j

    double precision :: a, lx(3), ly(3)

    do i = 1,NX

      if ( dble(i-1) * DELTAX < X_source .and. X_source <= dble(i) * DELTAX ) NX_source = i

    enddo

    a = X_source / DELTAX - dble(NX_source - 1)

    lx(1) = a * (a - 1.d0) / 2.d0
    lx(2) = - (a + 1.d0) * (a - 1.d0)
    lx(3) = (a + 1.d0) * a / 2.d0

    do i = 1,NY

      if ( dble(i-1) * DELTAY < Y_source .and. Y_source <= dble(i) * DELTAY ) NY_source = i

    enddo

    a = Y_source / DELTAY - dble(NY_source - 1)

    ly(1) = a * (a - 1.d0) / 2.d0
    ly(2) = - (a + 1.d0) * (a - 1.d0)
    ly(3) = (a + 1.d0) * a / 2.d0

    do i = 1,3
    do j = 1,3

      gx(i,j) = vx * lx(i) * ly(j) / DELTAX / DELTAY
      gy(i,j) = vy * lx(i) * ly(j) / DELTAX / DELTAY

    enddo
    enddo

  end subroutine

!##################################################################################

  subroutine set_source_interpolating_function_double_couple&
    ( NX,NY,DELTAX,DELTAY,X_source,Y_source,mxx,mxy,myx,myy,NX_source,NY_source,gx,gy )

    integer, intent(in) :: &
      NX, NY
    double precision, intent(in) :: &
      DELTAX, DELTAY, X_source, Y_source
    double precision, intent(in) :: &
      mxx, mxy, myx, myy

    integer, intent(out) :: &
      NX_source, NY_source

    double precision, intent(out) :: &
      gx(3,3), gy(3,3)

    integer :: i,j

    double precision :: a, lx(3), ly(3), dlx(3), dly(3)

    do i = 1,NX

      if ( dble(i-1) * DELTAX < X_source .and. X_source <= dble(i) * DELTAX ) NX_source = i

    enddo

    a = X_source / DELTAX - dble(NX_source - 1)

    lx(1) = a * (a - 1.d0) / 2.d0
    lx(2) = - (a + 1.d0) * (a - 1.d0)
    lx(3) = (a + 1.d0) * a / 2.d0

    dlx(1) = ( 2.0d0 * a - 1.0d0 ) / 2.0d0 / DELTAX
    dlx(2) = - 2.0d0 * a / DELTAX
    dlx(3) = ( 2.0d0 * a + 1.0d0 ) / 2.0d0 / DELTAX

    do i = 1,NY

      if ( dble(i-1) * DELTAY < Y_source .and. Y_source <= dble(i) * DELTAY ) NY_source = i

    enddo

    a = Y_source / DELTAY - dble(NY_source - 1)

    ly(1) = a * (a - 1.d0) / 2.d0
    ly(2) = - (a + 1.d0) * (a - 1.d0)
    ly(3) = (a + 1.d0) * a / 2.d0

    dly(1) = ( 2.0d0 * a - 1.0d0 ) / 2.0d0 / DELTAY
    dly(2) = - 2.0d0 * a / DELTAY
    dly(3) = ( 2.0d0 * a + 1.0d0 ) / 2.0d0 / DELTAY

    do i = 1,3
    do j = 1,3

      gx(i,j) = mxx * dlx(i) * ly(j) + mxy * lx(i) * dly(j)
      gy(i,j) = myx * dlx(i) * ly(j) + myy * lx(i) * dly(j)

    enddo
    enddo

    gx = gx / DELTAX / DELTAY
    gy = gy / DELTAX / DELTAY

  end subroutine

!##################################################################################



  subroutine set_reciever_interpolating_function&
    ( NX,NY,num_rec,DELTAX,DELTAY,X_rec,Y_rec,NX_rec,NY_rec,inpo_rec )

    integer, intent(in) :: &
      NX, NY, num_rec
    double precision, intent(in) :: &
      DELTAX, DELTAY, X_rec(num_rec), Y_rec(num_rec)

    integer, intent(out) :: &
      NX_rec(num_rec), NY_rec(num_rec)

    double precision, intent(out) :: &
      inpo_rec(3,3,num_rec)

    integer :: i, j, irec

    double precision :: a, lx(3), ly(3)

    do irec = 1,num_rec

      do i = 1,NX

        if ( dble(i-1) * DELTAX < X_rec(irec) .and. X_rec(irec) <= dble(i) * DELTAX ) NX_rec(irec) = i

      enddo

      a = X_rec(irec) / DELTAX - dble(NX_rec(irec) - 1)

      lx(1) = a * (a - 1.d0) / 2.d0
      lx(2) = - (a + 1.d0) * (a - 1.d0)
      lx(3) = (a + 1.d0) * a / 2.d0

      do i = 1,NY

        if ( dble(i-1) * DELTAY < Y_rec(irec) .and. Y_rec(irec) <= dble(i) * DELTAY ) NY_rec(irec) = i

      enddo

      a = Y_rec(irec) / DELTAY - dble(NY_rec(irec) - 1)

      ly(1) = a * (a - 1.d0) / 2.d0
      ly(2) = - (a + 1.d0) * (a - 1.d0)
      ly(3) = (a + 1.d0) * a / 2.d0

      do i = 1,3
      do j = 1,3

        inpo_rec(i,j,irec) = lx(i) * ly(j)

      enddo
      enddo

    enddo

  end subroutine set_reciever_interpolating_function
!##################################################################################

  subroutine compt_absorbing_boundary_Cerjan&
    ( NX,NY,u0x,u0y,u1x,u1y,u2x,u2y,N_Cerjan,dumping_rate_Cerjan )

    integer, intent(in) :: &
      NX, NY, N_Cerjan
    double precision, intent(in) :: &
      dumping_rate_Cerjan

    double precision, dimension(0:NX,0:NY), intent(inout) :: &
      u0x, u0y, u1x, u1y, u2x, u2y

    double precision :: r

    integer :: i, j

    do i = 0, N_Cerjan

      r = dble( N_Cerjan - i )
      r = exp( - dumping_rate_Cerjan * dumping_rate_Cerjan * r * r )

      u0x(i,:) = u0x(i,:) * r
      u0y(i,:) = u0y(i,:) * r
      u1x(i,:) = u1x(i,:) * r
      u1y(i,:) = u1y(i,:) * r
      u2x(i,:) = u2x(i,:) * r
      u2y(i,:) = u2y(i,:) * r

      j = NX - j

      u0x(j,:) = u0x(j,:) * r
      u0y(j,:) = u0y(j,:) * r
      u1x(j,:) = u1x(j,:) * r
      u1y(j,:) = u1y(j,:) * r
      u2x(j,:) = u2x(j,:) * r
      u2y(j,:) = u2y(j,:) * r

    enddo

    do i = 0, N_Cerjan

      r = dble( N_Cerjan - i )
      r = exp( - dumping_rate_Cerjan * dumping_rate_Cerjan * r * r )

      u0x(:,i) = u0x(:,i) * r
      u0y(:,i) = u0y(:,i) * r
      u1x(:,i) = u1x(:,i) * r
      u1y(:,i) = u1y(:,i) * r
      u2x(:,i) = u2x(:,i) * r
      u2y(:,i) = u2y(:,i) * r

    enddo

  end subroutine

!##################################################################################


  subroutine make_mask&
    ( NX,NY,rho,MASK_CELL,MASK_NODE,MASK1X,MASK1Y,MASK2X,MASK2Y )

    integer, intent(in) :: &
      NX, NY

    double precision, intent(in) :: &
      rho(NX,NY)

    integer, intent(out) :: &
      MASK_CELL(NX,NY), &
      MASK_NODE(0:NX,0:NY), &
      MASK1X(NX,0:NY), MASK1Y(0:NX,NY), &
      MASK2X(NX,NY), MASK2Y(NX,NY)

    integer :: i, j

! making masks for cells and nodes
    MASK_CELL = 0
    MASK_NODE = 0

    do j = 1,NY
    do i = 1,NX

      if ( rho(i,j) > 0.d0 ) then

        MASK_CELL(i,j) = 1

        MASK_NODE(i-1,j-1) = 1
        MASK_NODE(i-1,j) = 1
        MASK_NODE(i,j-1) = 1
        MASK_NODE(i,j) = 1

      end if

    enddo
    enddo

! making MASK1X
    MASK1X = 0

    do j = 0,NY
    do i = 1,NX

      if ( MASK_NODE(i-1,j) == 1 .and. MASK_NODE(i,j) == 1 ) MASK1X(i,j) = 1
      if ( i == NX ) cycle
      if ( MASK1X(i,j) == 1 .and. MASK_NODE(i+1,j) == 1 ) MASK1X(i,j) = 2

    enddo
    enddo

! making MASK2Y
    MASK2Y = 0

    do j = 1,NY
    do i = 1,NX

      if ( MASK_CELL(i,j) == 1 ) MASK2Y(i,j) = 1
      if ( j == NY ) cycle
      if ( MASK2Y(i,j) == 1 .and. MASK_CELL(i,j+1) == 1 ) MASK2Y(i,j) = 2

    enddo
    enddo

! making MASK1Y
    MASK1Y = 0

    do j = 1,NY
    do i = 0,NX

      if ( MASK_NODE(i,j-1) == 1 .and. MASK_NODE(i,j) == 1 ) MASK1Y(i,j) = 1
      if ( j== NY ) cycle
      if ( MASK1Y(i,j) == 1 .and. MASK_NODE(i,j+1) == 1 ) MASK1Y(i,j) = 2

    enddo
    enddo

! making MASK2X
    MASK2X = 0

    do j = 1,NY
    do i = 1,NX

      if (MASK_CELL(i,j) == 1 ) MASK2X(i,j) = 1
      if ( i == NX ) cycle
      if ( MASK2X(i,j) == 1 .and. MASK_CELL(i+1,j) == 1 ) MASK2X(i,j) = 2

    enddo
    enddo

  end subroutine

!##################################################################################

  subroutine set_physical_parameters&
    ( NX,NY,rho,C,model_csv,header )

    integer, intent(in) :: &
      NX, NY

    character(*), intent(in) :: &
      model_csv

    logical, intent(in) :: &
      header

    double precision, intent(out) :: &
      rho(NX,NY), C(3,3,NX,NY)

    double precision :: &
      density, vp, vs, mu, lambda, elas(3,3)

    integer :: &
      iline, ix, iy, i, j

    rho = - 80704.d0



    open (17, file=model_csv, status='old')

    if ( header ) read (17, '()')  !skip header

    do iline = 1, NX*NY

      read (17, *) ix, iy, density, vp, vs

      mu = density * vs * vs
      lambda = density * ( vp * vp - 2.d0 * vs * vs )

      elas = 0.d0

      do i = 1,2

        elas(i,i) = mu

        do j = 1,2

          elas(i,j) = elas(i,j) + lambda

        enddo

      enddo

      do i = 1,3

        elas(i,i) = elas(i,i) + mu

      enddo

      rho(ix,iy) = density
      C(:,:,ix,iy) = elas(:,:)

    enddo

    close (17)

!--check if csv is correctly inputted----
    do ix = 1,NX
    do iy = 1,NY

        if ( rho(ix,iy) < 0.d0 .or. &
             rho(1,iy) <= 1.d-10 .or. &
             rho(NX,iy) <= 1.d-10 .or. &
             rho(ix,1) <= 1.d-10 ) stop 'this model is invalid'

    enddo
    enddo

!----END-------------------------------

  end subroutine

!##################################################################################
  subroutine set_receiver_info&
   ( num_rec, X_rec, Y_rec, outfname_rec, recinfo )
    implicit none
    character(*), intent(in) :: recinfo

    integer, intent(out) :: num_rec

    double precision,allocatable, intent(out) :: &
      X_rec(:), Y_rec(:)

    character(100),allocatable, intent(out) :: &
      outfname_rec(:)

    character(1024) buffer
    integer :: iline, io, name_length

!    X_rec = - 80704.d0

    open (17, file=recinfo, status='old')
    num_rec =0
!    if ( header ) read (17, '()')  !skip header
    do
      read(17,'(a)',iostat=io) buffer
      buffer = adjustl(buffer)
      name_length = index(buffer,',')
      buffer = buffer(1:1)
      if(buffer =='#' .or. buffer == '!') cycle
      if (io/=0) exit
      if(name_length ==0) stop 'Something wrong.'
      if(name_length >101) stop 'Each length of each output path must be less than 100 characters.'
      num_rec = num_rec + 1
    enddo
    allocate(outfname_rec(num_rec))
    allocate(X_rec(num_rec), Y_rec(num_rec))
    rewind(17)
    iline =1
    do
      read (17,'(a)', iostat = io) buffer
      if (io/=0) exit
      buffer = adjustl(buffer)
      if(buffer(1:1) =='#' .or. buffer(1:1) == '!') cycle
      buffer = trim(buffer)
      read (buffer,*) outfname_rec(iline), X_rec(iline), Y_rec(iline)
      outfname_rec(iline) = trim(adjustl(outfname_rec(iline)))
      if ( X_rec(iline) < 0.d0 .or. Y_rec(iline) <0) stop 'reciever info file is invalid'
      iline = iline +1
    enddo
    close (17)

  end subroutine set_receiver_info

end program
