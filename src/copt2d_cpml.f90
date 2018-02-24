! ============================================================================
! Name        : copt2d_cpml.f90
! Author      : Kei Hasegawa
! Version     : 0.0.1
! Copyright   : It is Complicated.
! Description : Procedure for CPML.
! ============================================================================
module copt2d_cpml
  implicit none

contains

  subroutine compt_inverse_diagonal_mass&
    ( NX,NY,MASK1Y,MASK2X,rho,DELTAT,INV_MASS )

    integer, intent(in) :: &
      NX, NY, &
      MASK1Y(0:NX,NY), MASK2X(NX,NY)

    double precision, intent(in) :: &
      rho(NX,NY), DELTAT

    double precision, intent(out) :: &
      INV_MASS(0:NX,0:NY)

    double precision :: &
      tmp(0:NX,NY), E(3), Ec(2)

    integer :: i, j

    E = (/ 5.d0, 8.d0, -1.d0 /) / 12.d0
    Ec = (/ 1.d0, 1.d0 /) / 2.d0

    tmp = 0.d0
    INV_MASS = 0.d0

    do j = 1,NY
    do i = 1,NX

      if ( MASK2X(i,j) == 2 ) then

        tmp(i-1:i+1,j) = tmp(i-1:i+1,j) + rho(i,j) * E

      elseif ( MASK2X(i,j) == 1 ) then

        tmp(i-1:i,j) = tmp(i-1:i,j) + rho(i,j) * Ec

      endif

    enddo
    enddo

    do j = 1,NY
    do i = 0,NX

      if ( MASK1Y(i,j) == 2 ) then

        INV_MASS(i,j-1:j+1) = INV_MASS(i,j-1:j+1) + tmp(i,j) * E

      elseif ( MASK1Y(i,j) == 1 ) then

        INV_MASS(i,j-1:j) = INV_MASS(i,j-1:j) + tmp(i,j) * Ec

      endif

    enddo
    enddo

    do j = 0,NY
    do i = 0,NX

      if ( INV_MASS(i,j) > 1.d-10 ) &
        INV_MASS(i,j) = DELTAT * DELTAT / INV_MASS(i,j)

    enddo
    enddo

  end subroutine

!##################################################################################

  subroutine compt_corrector_force&
    ( NX,NY,MASK1X,MASK1Y,MASK2X,MASK2Y,ux,uy,C,rho,DELTAT,DELTAX,DELTAY, &
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
      Hux,Huy )

    integer, intent(in) :: &
      NX, NY, &
      MASK1X(NX,0:NY), MASK1Y(0:NX,NY), &
      MASK2X(NX,NY), MASK2Y(NX,NY)

    double precision, intent(in) :: &
      ux(0:NX,0:NY), uy(0:NX,0:NY), &
      C(3,3,NX,NY), rho(NX,NY), &
      DELTAT, DELTAX, DELTAY


!----------------valuables for CPML-------------------------------------
    integer, intent(in) :: &
      NX_CPML_lef, NX_CPML_rig, NY_CPML_btm

    double precision, intent(in) :: &
      a1_lef(NX_CPML_lef), a1_rig(NX - NX_CPML_rig+1:NX), a1_btm(NY_CPML_btm), &
      b1_lef(NX_CPML_lef), b1_rig(NX - NX_CPML_rig+1:NX), b1_btm(NY_CPML_btm), &
      a2_lef(0:NX_CPML_lef), a2_rig(NX - NX_CPML_rig:NX), a2_btm(0:NY_CPML_btm), &
      b2_lef(0:NX_CPML_lef), b2_rig(NX - NX_CPML_rig:NX), b2_btm(0:NY_CPML_btm)

    double precision, intent(inout) :: &
      dux_dx_lef(NX_CPML_lef,0:NY), &
      duy_dx_lef(NX_CPML_lef,0:NY), &
      dux_dx_rig(NX - NX_CPML_rig+1:NX,0:NY), &
      duy_dx_rig(NX - NX_CPML_rig+1:NX,0:NY), &
      dux_dy_btm(0:NX,NY_CPML_btm), &
      duy_dy_btm(0:NX,NY_CPML_btm), &
      dTxx_dx_lef(0:NX_CPML_lef,0:NY), &
      dTyx_dx_lef(0:NX_CPML_lef,0:NY), &
      dTxx_dx_rig(NX - NX_CPML_rig:NX,0:NY), &
      dTyx_dx_rig(NX - NX_CPML_rig:NX,0:NY), &
      dTxy_dy_btm(0:NX,0:NY_CPML_btm), &
      dTyy_dy_btm(0:NX,0:NY_CPML_btm)

    double precision :: &
      dTxx_dx_lef_tmp(0:NX_CPML_lef,0:NY), &
      dTyx_dx_lef_tmp(0:NX_CPML_lef,0:NY), &
      dTxx_dx_rig_tmp(NX - NX_CPML_rig:NX,0:NY), &
      dTyx_dx_rig_tmp(NX - NX_CPML_rig:NX,0:NY), &
      dTxy_dy_btm_tmp(0:NX,0:NY_CPML_btm), &
      dTyy_dy_btm_tmp(0:NX,0:NY_CPML_btm)

    double precision, intent(inout) :: &
      mem1_xx_lef(NX_CPML_lef,0:NY), &
      mem1_xy_lef(NX_CPML_lef,0:NY), &
      mem1_xx_rig(NX - NX_CPML_rig+1:NX,0:NY), &
      mem1_xy_rig(NX - NX_CPML_rig+1:NX,0:NY), &
      mem1_yx_btm(0:NX,NY_CPML_btm), &
      mem1_yy_btm(0:NX,NY_CPML_btm), &
      mem2_xx_lef(0:NX_CPML_lef,0:NY), &
      mem2_xy_lef(0:NX_CPML_lef,0:NY), &
      mem2_xx_rig(NX - NX_CPML_rig:NX,0:NY), &
      mem2_xy_rig(NX - NX_CPML_rig:NX,0:NY), &
      mem2_yx_btm(0:NX,0:NY_CPML_btm), &
      mem2_yy_btm(0:NX,0:NY_CPML_btm)

!-----------END--------------------------------------------------------


    double precision, intent(out) :: &
      Hux(0:NX,0:NY), Huy(0:NX,0:NY)

    double precision :: &
      tmp1_xx(NX,0:NY), tmp1_xy(NX,0:NY), &
      tmp1_yx(0:NX,NY), tmp1_yy(0:NX,NY), &
      tmp2_xx(NX,0:NY), tmp2_xy(NX,0:NY), &
      tmp2_yx(0:NX,NY), tmp2_yy(0:NX,NY), tmp3, &
      Pxx, Pxy, Pyx, Pyy, &
      Qxx, Qxy, Qyx, Qyy, &
      sgm_xx, sgm_xy, sgm_yx, sgm_yy, &
      gmm_xx, gmm_xy, gmm_yx, gmm_yy, &
      dx2_dt2, dy2_dt2

    double precision :: &
      E(2), F(2), Delx(2), Dely(2)

    integer :: i, j

    E = (/ 1.d0, 1.d0 /) / 2.d0
    F = (/ 1.d0, -1.d0 /) / 2.d0
    Delx = (/ -1.d0, 1.d0 /) / DELTAX
    Dely = (/ -1.d0, 1.d0 /) / DELTAY

    Hux = 0.d0
    Huy = 0.d0

    tmp2_xx = 0.d0
    tmp2_xy = 0.d0
    tmp2_yx = 0.d0
    tmp2_yy = 0.d0

    dx2_dt2 = DELTAX * DELTAX / DELTAT / DELTAT
    dy2_dt2 = DELTAY * DELTAY / DELTAT / DELTAT

    do j = 0,NY
    do i = 1,NX

      if ( MASK1X(i,j) /= 0 ) then

        tmp1_xx(i,j) = dot_product( Delx, ux(i-1:i,j) )
        tmp1_xy(i,j) = dot_product( Delx, uy(i-1:i,j) )

      else

        tmp1_xx(i,j) = 0.d0
        tmp1_xy(i,j) = 0.d0

      endif

!-----------implement CPML for x-axis-----------------------------------

      if ( i <= NX_CPML_lef ) then

        mem1_xx_lef(i,j) = a1_lef(i) * mem1_xx_lef(i,j) + b1_lef(i) * ( tmp1_xx(i,j) + dux_dx_lef(i,j) )
        mem1_xy_lef(i,j) = a1_lef(i) * mem1_xy_lef(i,j) + b1_lef(i) * ( tmp1_xy(i,j) + duy_dx_lef(i,j) )

        dux_dx_lef(i,j) = tmp1_xx(i,j)
        duy_dx_lef(i,j) = tmp1_xy(i,j)

        tmp1_xx(i,j) = tmp1_xx(i,j) + mem1_xx_lef(i,j)
        tmp1_xy(i,j) = tmp1_xy(i,j) + mem1_xy_lef(i,j)

      elseif ( i > NX - NX_CPML_rig  ) then

        mem1_xx_rig(i,j) = a1_rig(i) * mem1_xx_rig(i,j) + b1_rig(i) * ( tmp1_xx(i,j) + dux_dx_rig(i,j) )
        mem1_xy_rig(i,j) = a1_rig(i) * mem1_xy_rig(i,j) + b1_rig(i) * ( tmp1_xy(i,j) + duy_dx_rig(i,j) )

        dux_dx_rig(i,j) = tmp1_xx(i,j)
        duy_dx_rig(i,j) = tmp1_xy(i,j)

        tmp1_xx(i,j) = tmp1_xx(i,j) + mem1_xx_rig(i,j)
        tmp1_xy(i,j) = tmp1_xy(i,j) + mem1_xy_rig(i,j)

      endif

!---------------END----------------------------------------------------

    enddo
    enddo

    do j = 1,NY
    do i = 0,NX

      if ( MASK1Y(i,j) /= 0 ) then

        tmp1_yx(i,j) = dot_product( Dely, ux(i,j-1:j) )
        tmp1_yy(i,j) = dot_product( Dely, uy(i,j-1:j) )

      else

        tmp1_yx(i,j) = 0.d0
        tmp1_yy(i,j) = 0.d0

      endif

!------------------implement CPML for y-axis----------------------------
      if ( j <= NY_CPML_btm ) then

        mem1_yx_btm(i,j) = a1_btm(j) * mem1_yx_btm(i,j) + b1_btm(j) * ( tmp1_yx(i,j) + dux_dy_btm(i,j) )
        mem1_yy_btm(i,j) = a1_btm(j) * mem1_yy_btm(i,j) + b1_btm(j) * ( tmp1_yy(i,j) + duy_dy_btm(i,j) )

        dux_dy_btm(i,j) = tmp1_yx(i,j)
        duy_dy_btm(i,j) = tmp1_yy(i,j)

        tmp1_yx(i,j) = tmp1_yx(i,j) + mem1_yx_btm(i,j)
        tmp1_yy(i,j) = tmp1_yy(i,j) + mem1_yy_btm(i,j)

      endif

!---------------END-----------------------------------------------------

    enddo
    enddo

    do j = 1,NY
    do i = 1,NX

      if ( MASK2Y(i,j) /= 0 ) then

        Pxx = dot_product( E, tmp1_xx(i,j-1:j) ) / 12.d0
        Pxy = dot_product( E, tmp1_xy(i,j-1:j) ) / 12.d0
        Qxx = dot_product( F, tmp1_xx(i,j-1:j) ) / 12.d0
        Qxy = dot_product( F, tmp1_xy(i,j-1:j) ) / 12.d0

      else

        Pxx = 0.d0
        Pxy = 0.d0
        Qxx = 0.d0
        Qxy = 0.d0

      endif

      if ( MASK2X(i,j) /= 0 ) then

        Pyx = dot_product( E, tmp1_yx(i-1:i,j) ) / 12.d0
        Pyy = dot_product( E, tmp1_yy(i-1:i,j) ) / 12.d0
        Qyx = dot_product( F, tmp1_yx(i-1:i,j) ) / 12.d0
        Qyy = dot_product( F, tmp1_yy(i-1:i,j) ) / 12.d0

      else

        Pxx = 0.d0
        Pxy = 0.d0
        Qxx = 0.d0
        Qxy = 0.d0

      endif

      tmp3 = Pxy + Pyx
      sgm_xx = C(1,1,i,j) * Pxx + C(1,3,i,j) * tmp3 + C(1,2,i,j) * Pyy
      sgm_yy = C(2,1,i,j) * Pxx + C(2,3,i,j) * tmp3 + C(2,2,i,j) * Pyy
      tmp3   = C(3,1,i,j) * Pxx + C(3,3,i,j) * tmp3 + C(3,2,i,j) * Pyy
      gmm_xx = C(1,1,i,j) * Qxx + C(1,3,i,j) * Qxy
      gmm_xy = C(3,1,i,j) * Qxx + C(3,3,i,j) * Qxy
      gmm_yx = C(3,3,i,j) * Qyx + C(3,2,i,j) * Qyy
      gmm_yy = C(2,3,i,j) * Qyx + C(2,2,i,j) * Qyy

      sgm_xx = sgm_xx - rho(i,j) * dx2_dt2 * Pxx
      sgm_xy = tmp3   - rho(i,j) * dx2_dt2 * Pxy
      sgm_yx = tmp3   - rho(i,j) * dy2_dt2 * Pyx
      sgm_yy = sgm_yy - rho(i,j) * dy2_dt2 * Pyy
      gmm_xx = gmm_xx - rho(i,j) * dx2_dt2 * Qxx
      gmm_xy = gmm_xy - rho(i,j) * dx2_dt2 * Qxy
      gmm_yx = gmm_yx - rho(i,j) * dy2_dt2 * Qyx
      gmm_yy = gmm_yy - rho(i,j) * dy2_dt2 * Qyy

      if ( MASK2Y(i,j) /= 0 ) then

        tmp2_xx(i,j-1:j) = tmp2_xx(i,j-1:j) + sgm_xx * E + gmm_xx * F
        tmp2_xy(i,j-1:j) = tmp2_xy(i,j-1:j) + sgm_xy * E + gmm_xy * F

      endif

      if ( MASK2X(i,j) /= 0 ) then

        tmp2_yx(i-1:i,j) = tmp2_yx(i-1:i,j) + sgm_yx * E + gmm_yx * F
        tmp2_yy(i-1:i,j) = tmp2_yy(i-1:i,j) + sgm_yy * E + gmm_yy * F

      endif

    enddo
    enddo

!-----------set valuables for CPML to zero------------------------------------
    dTxx_dx_lef_tmp = 0.d0
    dTyx_dx_lef_tmp = 0.d0
    dTxx_dx_rig_tmp = 0.d0
    dTyx_dx_rig_tmp = 0.d0
    dTxy_dy_btm_tmp = 0.d0
    dTyy_dy_btm_tmp = 0.d0
!-----------------END--------------------------------------------------------

    do j = 0,NY
    do i = 1,NX

      if ( MASK1X(i,j) /= 0) then

        Hux(i-1:i,j) = Hux(i-1:i,j) + tmp2_xx(i,j) * Delx
        Huy(i-1:i,j) = Huy(i-1:i,j) + tmp2_xy(i,j) * Delx

!-----------implement CPML for x-axis------------------------------------------
        if ( i <= NX_CPML_lef ) then

          dTxx_dx_lef_tmp(i-1:i,j) = dTxx_dx_lef_tmp(i-1:i,j) + tmp2_xx(i,j) * Delx
          dTyx_dx_lef_tmp(i-1:i,j) = dTyx_dx_lef_tmp(i-1:i,j) + tmp2_xy(i,j) * Delx

        elseif ( i > NX - NX_CPML_rig  ) then

          dTxx_dx_rig_tmp(i-1:i,j) = dTxx_dx_rig_tmp(i-1:i,j) + tmp2_xx(i,j) * Delx
          dTyx_dx_rig_tmp(i-1:i,j) = dTyx_dx_rig_tmp(i-1:i,j) + tmp2_xy(i,j) * Delx

        endif
!-----------------END----------------------------------------------------------

      endif

    enddo
    enddo

!--------------implement CPML for x-axis--------------------------

    do i = 0, NX_CPML_lef

      mem2_xx_lef(i,:) = a2_lef(i) * mem2_xx_lef(i,:) &
        + b2_lef(i) * ( dTxx_dx_lef_tmp(i,:) + dTxx_dx_lef(i,:) )

      mem2_xy_lef(i,:) = a2_lef(i) * mem2_xy_lef(i,:) &
        + b2_lef(i) * ( dTyx_dx_lef_tmp(i,:) + dTyx_dx_lef(i,:) )

      Hux(i,:) = Hux(i,:) + mem2_xx_lef(i,:)
      Huy(i,:) = Huy(i,:) + mem2_xy_lef(i,:)

    enddo

    dTxx_dx_lef = dTxx_dx_lef_tmp
    dTyx_dx_lef = dTyx_dx_lef_tmp

    do i = NX - NX_CPML_rig, NX

      mem2_xx_rig(i,:) = a2_rig(i) * mem2_xx_rig(i,:) &
        + b2_rig(i) * ( dTxx_dx_rig_tmp(i,:) + dTxx_dx_rig(i,:) )

      mem2_xy_rig(i,:) = a2_rig(i) * mem2_xy_rig(i,:) &
        + b2_rig(i) * ( dTyx_dx_rig_tmp(i,:) + dTyx_dx_rig(i,:) )

      Hux(i,:) = Hux(i,:) + mem2_xx_rig(i,:)
      Huy(i,:) = Huy(i,:) + mem2_xy_rig(i,:)

    enddo

    dTxx_dx_rig = dTxx_dx_rig_tmp
    dTyx_dx_rig = dTyx_dx_rig_tmp

!-----------END--------------------------------------------------

    do j = 1,NY
    do i = 0,NX

      if ( MASK1Y(i,j) /= 0) then

        Hux(i,j-1:j) = Hux(i,j-1:j) + tmp2_yx(i,j) * Dely
        Huy(i,j-1:j) = Huy(i,j-1:j) + tmp2_yy(i,j) * Dely

!----------------implement CPML for y-axis-----------------------

        if ( j <= NY_CPML_btm ) then

          dTxy_dy_btm_tmp(i,j-1:j) = dTxy_dy_btm_tmp(i,j-1:j) + tmp2_yx(i,j) * Dely
          dTyy_dy_btm_tmp(i,j-1:j) = dTyy_dy_btm_tmp(i,j-1:j) + tmp2_yy(i,j) * Dely

        endif
!-----------END--------------------------------------------------

      endif

    enddo
    enddo

!-----------------implement CPML for y-axis------------------
    do j = 0, NY_CPML_btm

      mem2_yx_btm(:,j) = a2_btm(j) * mem2_yx_btm(:,j) &
        + b2_btm(j) * ( dTxy_dy_btm_tmp(:,j) + dTxy_dy_btm(:,j) )

      mem2_yy_btm(:,j) = a2_btm(j) * mem2_yy_btm(:,j) &
        + b2_btm(j) * ( dTyy_dy_btm_tmp(:,j) + dTyy_dy_btm(:,j) )

      Hux(:,j) = Hux(:,j) + mem2_yx_btm(:,j)
      Huy(:,j) = Huy(:,j) + mem2_yy_btm(:,j)

    enddo

    dTxy_dy_btm = dTxy_dy_btm_tmp
    dTyy_dy_btm = dTyy_dy_btm_tmp


!-----------END----------------------------------------------


  end subroutine


!##################################################################################

  subroutine compt_traction_force&
    ( NX,NY,MASK1X,MASK1Y,MASK2X,MASK2Y,ux,uy,C,DELTAX,DELTAY, &
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
      Hux,Huy )

    integer, intent(in) :: &
      NX, NY, &
      MASK1X(NX,0:NY), MASK1Y(0:NX,NY), &
      MASK2X(NX,NY), MASK2Y(NX,NY)

    double precision, intent(in) :: &
      ux(0:NX,0:NY), uy(0:NX,0:NY), &
      C(3,3,NX,NY), &
      DELTAX, DELTAY

!----------------valuables for CPML-------------------------------------
    integer, intent(in) :: &
      NX_CPML_lef, NX_CPML_rig, NY_CPML_btm

    double precision, intent(in) :: &
      a1_lef(NX_CPML_lef), a1_rig(NX - NX_CPML_rig+1:NX), a1_btm(NY_CPML_btm), &
      b1_lef(NX_CPML_lef), b1_rig(NX - NX_CPML_rig+1:NX), b1_btm(NY_CPML_btm), &
      a2_lef(0:NX_CPML_lef), a2_rig(NX - NX_CPML_rig:NX), a2_btm(0:NY_CPML_btm), &
      b2_lef(0:NX_CPML_lef), b2_rig(NX - NX_CPML_rig:NX), b2_btm(0:NY_CPML_btm)

    double precision, intent(inout) :: &
      dux_dx_lef(NX_CPML_lef,0:NY), &
      duy_dx_lef(NX_CPML_lef,0:NY), &
      dux_dx_rig(NX - NX_CPML_rig+1:NX,0:NY), &
      duy_dx_rig(NX - NX_CPML_rig+1:NX,0:NY), &
      dux_dy_btm(0:NX,NY_CPML_btm), &
      duy_dy_btm(0:NX,NY_CPML_btm), &
      dTxx_dx_lef(0:NX_CPML_lef,0:NY), &
      dTyx_dx_lef(0:NX_CPML_lef,0:NY), &
      dTxx_dx_rig(NX - NX_CPML_rig:NX,0:NY), &
      dTyx_dx_rig(NX - NX_CPML_rig:NX,0:NY), &
      dTxy_dy_btm(0:NX,0:NY_CPML_btm), &
      dTyy_dy_btm(0:NX,0:NY_CPML_btm)

    double precision :: &
      dTxx_dx_lef_tmp(0:NX_CPML_lef,0:NY), &
      dTyx_dx_lef_tmp(0:NX_CPML_lef,0:NY), &
      dTxx_dx_rig_tmp(NX - NX_CPML_rig:NX,0:NY), &
      dTyx_dx_rig_tmp(NX - NX_CPML_rig:NX,0:NY), &
      dTxy_dy_btm_tmp(0:NX,0:NY_CPML_btm), &
      dTyy_dy_btm_tmp(0:NX,0:NY_CPML_btm)


    double precision, intent(inout) :: &
      mem1_xx_lef(NX_CPML_lef,0:NY), &
      mem1_xy_lef(NX_CPML_lef,0:NY), &
      mem1_xx_rig(NX - NX_CPML_rig+1:NX,0:NY), &
      mem1_xy_rig(NX - NX_CPML_rig+1:NX,0:NY), &
      mem1_yx_btm(0:NX,NY_CPML_btm), &
      mem1_yy_btm(0:NX,NY_CPML_btm), &
      mem2_xx_lef(0:NX_CPML_lef,0:NY), &
      mem2_xy_lef(0:NX_CPML_lef,0:NY), &
      mem2_xx_rig(NX - NX_CPML_rig:NX,0:NY), &
      mem2_xy_rig(NX - NX_CPML_rig:NX,0:NY), &
      mem2_yx_btm(0:NX,0:NY_CPML_btm), &
      mem2_yy_btm(0:NX,0:NY_CPML_btm)

!-----------END--------------------------------------------------------

    double precision, intent(out) :: &
      Hux(0:NX,0:NY), Huy(0:NX,0:NY)

    double precision :: &
      tmp1_xx(NX,0:NY), tmp1_xy(NX,0:NY), &
      tmp1_yx(0:NX,NY), tmp1_yy(0:NX,NY), &
      tmp2_xx(NX,0:NY), tmp2_xy(NX,0:NY), &
      tmp2_yx(0:NX,NY), tmp2_yy(0:NX,NY), &
      Pxx, Pxy, Pyx, Pyy, &
      Qxx, Qxy, Qyx, Qyy, &
      sgm_xx, sgm_xy, sgm_yy, &
      gmm_xx, gmm_xy, gmm_yx, gmm_yy

    double precision :: &
      E(3), Ec(2), F(3), Fc(2), Delx(2), Dely(2)

    integer :: i, j

    E = (/ 5.d0, 8.d0, -1.d0 /) / 12.d0
    Ec = (/ 1.d0, 1.d0 /) / 2.d0
    F = (/ 1.d0, -2.d0, 1.d0 /) * sqrt(5.d0) / 12.d0
    Fc = (/ 1.d0, -1.d0 /) / 2.d0
    Delx = (/ -1.d0, 1.d0 /) / DELTAX
    Dely = (/ -1.d0, 1.d0 /) / DELTAY

    tmp2_xx = 0.d0
    tmp2_xy = 0.d0
    tmp2_yx = 0.d0
    tmp2_yy = 0.d0

    Hux = 0.d0
    Huy = 0.d0


    do j = 0,NY
    do i = 1,NX

      if ( MASK1X(i,j) /= 0 ) then

        tmp1_xx(i,j) = dot_product( Delx, ux(i-1:i,j) )
        tmp1_xy(i,j) = dot_product( Delx, uy(i-1:i,j) )

      else

        tmp1_xx(i,j) = 0.d0
        tmp1_xy(i,j) = 0.d0

      endif


!-----------implement CPML for x-axis-----------------------------------
      if ( i <= NX_CPML_lef ) then

        mem1_xx_lef(i,j) = a1_lef(i) * mem1_xx_lef(i,j) + b1_lef(i) * ( tmp1_xx(i,j) + dux_dx_lef(i,j) )
        mem1_xy_lef(i,j) = a1_lef(i) * mem1_xy_lef(i,j) + b1_lef(i) * ( tmp1_xy(i,j) + duy_dx_lef(i,j) )

        dux_dx_lef(i,j) = tmp1_xx(i,j)
        duy_dx_lef(i,j) = tmp1_xy(i,j)

        tmp1_xx(i,j) = tmp1_xx(i,j) + mem1_xx_lef(i,j)
        tmp1_xy(i,j) = tmp1_xy(i,j) + mem1_xy_lef(i,j)

      elseif ( i > NX - NX_CPML_rig  ) then

        mem1_xx_rig(i,j) = a1_rig(i) * mem1_xx_rig(i,j) + b1_rig(i) * ( tmp1_xx(i,j) + dux_dx_rig(i,j) )
        mem1_xy_rig(i,j) = a1_rig(i) * mem1_xy_rig(i,j) + b1_rig(i) * ( tmp1_xy(i,j) + duy_dx_rig(i,j) )

        dux_dx_rig(i,j) = tmp1_xx(i,j)
        duy_dx_rig(i,j) = tmp1_xy(i,j)

        tmp1_xx(i,j) = tmp1_xx(i,j) + mem1_xx_rig(i,j)
        tmp1_xy(i,j) = tmp1_xy(i,j) + mem1_xy_rig(i,j)

      endif

!---------------END----------------------------------------------------

    enddo
    enddo

    do j = 1,NY
    do i = 0,NX

      if ( MASK1Y(i,j) /= 0 ) then

        tmp1_yx(i,j) = dot_product( Dely, ux(i,j-1:j) )
        tmp1_yy(i,j) = dot_product( Dely, uy(i,j-1:j) )

      else

        tmp1_yx(i,j) = 0.d0
        tmp1_yy(i,j) = 0.d0

      endif

!------------------implement CPML for y-axis----------------------------
      if ( j <= NY_CPML_btm ) then

        mem1_yx_btm(i,j) = a1_btm(j) * mem1_yx_btm(i,j) + b1_btm(j) * ( tmp1_yx(i,j) + dux_dy_btm(i,j) )
        mem1_yy_btm(i,j) = a1_btm(j) * mem1_yy_btm(i,j) + b1_btm(j) * ( tmp1_yy(i,j) + duy_dy_btm(i,j) )

        dux_dy_btm(i,j) = tmp1_yx(i,j)
        duy_dy_btm(i,j) = tmp1_yy(i,j)

        tmp1_yx(i,j) = tmp1_yx(i,j) + mem1_yx_btm(i,j)
        tmp1_yy(i,j) = tmp1_yy(i,j) + mem1_yy_btm(i,j)

      endif

!---------------END-----------------------------------------------------

    enddo
    enddo


    do j = 1,NY
    do i = 1,NX

      if ( MASK2Y(i,j) == 2 ) then

        Pxx = dot_product( E, tmp1_xx(i,j-1:j+1) )
        Pxy = dot_product( E, tmp1_xy(i,j-1:j+1) )
        Qxx = dot_product( F, tmp1_xx(i,j-1:j+1) )
        Qxy = dot_product( F, tmp1_xy(i,j-1:j+1) )

      elseif ( MASK2Y(i,j) == 1 ) then
        Pxx = dot_product( Ec, tmp1_xx(i,j-1:j) )
        Pxy = dot_product( Ec, tmp1_xy(i,j-1:j) )
        Qxx = dot_product( Fc, tmp1_xx(i,j-1:j) )
        Qxy = dot_product( Fc, tmp1_xy(i,j-1:j) )

      else

        Pxx = 0.d0
        Pxy = 0.d0
        Qxx = 0.d0
        Qxy = 0.d0

      endif

      if ( MASK2X(i,j) == 2 ) then

        Pyx = dot_product( E, tmp1_yx(i-1:i+1,j) )
        Pyy = dot_product( E, tmp1_yy(i-1:i+1,j) )
        Qyx = dot_product( F, tmp1_yx(i-1:i+1,j) )
        Qyy = dot_product( F, tmp1_yy(i-1:i+1,j) )

      elseif ( MASK2X(i,j) == 1 ) then

        Pyx = dot_product( Ec, tmp1_yx(i-1:i,j) )
        Pyy = dot_product( Ec, tmp1_yy(i-1:i,j) )
        Qyx = dot_product( Fc, tmp1_yx(i-1:i,j) )
        Qyy = dot_product( Fc, tmp1_yy(i-1:i,j) )

      else

        Pxx = 0.d0
        Pxy = 0.d0
        Qxx = 0.d0
        Qxy = 0.d0

      endif


      Pxy = Pxy + Pyx
      sgm_xx = C(1,1,i,j) * Pxx + C(1,3,i,j) * Pxy + C(1,2,i,j) * Pyy
      sgm_xy = C(3,1,i,j) * Pxx + C(3,3,i,j) * Pxy + C(3,2,i,j) * Pyy
      sgm_yy = C(2,1,i,j) * Pxx + C(2,3,i,j) * Pxy + C(2,2,i,j) * Pyy

      gmm_xx = C(1,1,i,j) * Qxx + C(1,3,i,j) * Qxy
      gmm_xy = C(3,1,i,j) * Qxx + C(3,3,i,j) * Qxy
      gmm_yx = C(3,3,i,j) * Qyx + C(3,2,i,j) * Qyy
      gmm_yy = C(2,3,i,j) * Qyx + C(2,2,i,j) * Qyy


      if ( MASK2Y(i,j) == 2 ) then

        tmp2_xx(i,j-1:j+1) = tmp2_xx(i,j-1:j+1) + sgm_xx * E + gmm_xx * F
        tmp2_xy(i,j-1:j+1) = tmp2_xy(i,j-1:j+1) + sgm_xy * E + gmm_xy * F

      elseif ( MASK2Y(i,j) == 1 ) then

        tmp2_xx(i,j-1:j) = tmp2_xx(i,j-1:j) + sgm_xx * Ec + gmm_xx * Fc
        tmp2_xy(i,j-1:j) = tmp2_xy(i,j-1:j) + sgm_xy * Ec + gmm_xy * Fc

      endif

      if ( MASK2X(i,j) == 2 ) then

        tmp2_yx(i-1:i+1,j) = tmp2_yx(i-1:i+1,j) + sgm_xy * E + gmm_yx * F
        tmp2_yy(i-1:i+1,j) = tmp2_yy(i-1:i+1,j) + sgm_yy * E + gmm_yy * F

      elseif ( MASK2X(i,j) == 1 ) then

        tmp2_yx(i-1:i,j) = tmp2_yx(i-1:i,j) + sgm_xy * Ec + gmm_yx * Fc
        tmp2_yy(i-1:i,j) = tmp2_yy(i-1:i,j) + sgm_yy * Ec + gmm_yy * Fc

      endif

    enddo
    enddo


!-----------set valuables for CPML to zero------------------------------------
    dTxx_dx_lef_tmp = 0.d0
    dTyx_dx_lef_tmp = 0.d0
    dTxx_dx_rig_tmp = 0.d0
    dTyx_dx_rig_tmp = 0.d0
    dTxy_dy_btm_tmp = 0.d0
    dTyy_dy_btm_tmp = 0.d0
!-----------------END--------------------------------------------------------

    do j = 0,NY
    do i = 1,NX

      if ( MASK1X(i,j) /= 0) then

        Hux(i-1:i,j) = Hux(i-1:i,j) + tmp2_xx(i,j) * Delx
        Huy(i-1:i,j) = Huy(i-1:i,j) + tmp2_xy(i,j) * Delx


!-----------implement CPML for x-axis------------------------------------------
          if ( i <= NX_CPML_lef ) then

            dTxx_dx_lef_tmp(i-1:i,j) = dTxx_dx_lef_tmp(i-1:i,j) + tmp2_xx(i,j) * Delx
            dTyx_dx_lef_tmp(i-1:i,j) = dTyx_dx_lef_tmp(i-1:i,j) + tmp2_xy(i,j) * Delx

          elseif ( i > NX - NX_CPML_rig  ) then

            dTxx_dx_rig_tmp(i-1:i,j) = dTxx_dx_rig_tmp(i-1:i,j) + tmp2_xx(i,j) * Delx
            dTyx_dx_rig_tmp(i-1:i,j) = dTyx_dx_rig_tmp(i-1:i,j) + tmp2_xy(i,j) * Delx

          endif
!-----------------END--------------------------------------------------------

      endif

    enddo
    enddo

!--------------implement CPML for x-axis--------------------------

    do i = 0, NX_CPML_lef

      mem2_xx_lef(i,:) = a2_lef(i) * mem2_xx_lef(i,:) &
        + b2_lef(i) * ( dTxx_dx_lef_tmp(i,:) + dTxx_dx_lef(i,:) )

      mem2_xy_lef(i,:) = a2_lef(i) * mem2_xy_lef(i,:) &
        + b2_lef(i) * ( dTyx_dx_lef_tmp(i,:) + dTyx_dx_lef(i,:) )

      Hux(i,:) = Hux(i,:) + mem2_xx_lef(i,:)
      Huy(i,:) = Huy(i,:) + mem2_xy_lef(i,:)

    enddo

    dTxx_dx_lef = dTxx_dx_lef_tmp
    dTyx_dx_lef = dTyx_dx_lef_tmp

    do i = NX - NX_CPML_rig, NX

      mem2_xx_rig(i,:) = a2_rig(i) * mem2_xx_rig(i,:) &
        + b2_rig(i) * ( dTxx_dx_rig_tmp(i,:) + dTxx_dx_rig(i,:) )

      mem2_xy_rig(i,:) = a2_rig(i) * mem2_xy_rig(i,:) &
        + b2_rig(i) * ( dTyx_dx_rig_tmp(i,:) + dTyx_dx_rig(i,:) )

      Hux(i,:) = Hux(i,:) + mem2_xx_rig(i,:)
      Huy(i,:) = Huy(i,:) + mem2_xy_rig(i,:)

    enddo

    dTxx_dx_rig = dTxx_dx_rig_tmp
    dTyx_dx_rig = dTyx_dx_rig_tmp

!-----------END--------------------------------------------------


    do j = 1,NY
    do i = 0,NX

      if ( MASK1Y(i,j) /= 0) then

        Hux(i,j-1:j) = Hux(i,j-1:j) + tmp2_yx(i,j) * Dely
        Huy(i,j-1:j) = Huy(i,j-1:j) + tmp2_yy(i,j) * Dely

!----------------implement CPML for y-axis-----------------------

        if ( j <= NY_CPML_btm ) then

          dTxy_dy_btm_tmp(i,j-1:j) = dTxy_dy_btm_tmp(i,j-1:j) + tmp2_yx(i,j) * Dely
          dTyy_dy_btm_tmp(i,j-1:j) = dTyy_dy_btm_tmp(i,j-1:j) + tmp2_yy(i,j) * Dely

        endif
!-----------END--------------------------------------------------

      endif

    enddo
    enddo

!-----------------implement CPML for y-axis------------------
    do j = 0, NY_CPML_btm

      mem2_yx_btm(:,j) = a2_btm(j) * mem2_yx_btm(:,j) &
        + b2_btm(j) * ( dTxy_dy_btm_tmp(:,j) + dTxy_dy_btm(:,j) )

      mem2_yy_btm(:,j) = a2_btm(j) * mem2_yy_btm(:,j) &
        + b2_btm(j) * ( dTyy_dy_btm_tmp(:,j) + dTyy_dy_btm(:,j) )

      Hux(:,j) = Hux(:,j) + mem2_yx_btm(:,j)
      Huy(:,j) = Huy(:,j) + mem2_yy_btm(:,j)

    enddo

    dTxy_dy_btm = dTxy_dy_btm_tmp
    dTyy_dy_btm = dTyy_dy_btm_tmp
!-----------END----------------------------------------------

  end subroutine
!##################################################################################

  subroutine set_CPML_parameters&
    ( NX,NY,DELTAT,DELTAX,DELTAY, &
      NX_CPML_lef, NX_CPML_rig, NY_CPML_btm, &
      vp, reflection_rate, source_freq, NPOWER, &
      a1_lef, a1_rig, a1_btm, &
      b1_lef, b1_rig, b1_btm, &
      a2_lef, a2_rig, a2_btm, &
      b2_lef, b2_rig, b2_btm )

    integer, intent(in) :: &
      NX, NY, &
      NX_CPML_lef, NX_CPML_rig, NY_CPML_btm, &
      NPOWER

    double precision, intent(in) :: &
      DELTAT, DELTAX, DELTAY

    double precision, intent(in) :: &
      vp, reflection_rate, source_freq

    double precision, intent(out) :: &
      a1_lef(NX_CPML_lef), a1_rig(NX - NX_CPML_rig+1:NX), a1_btm(NY_CPML_btm), &
      b1_lef(NX_CPML_lef), b1_rig(NX - NX_CPML_rig+1:NX), b1_btm(NY_CPML_btm), &
      a2_lef(0:NX_CPML_lef), a2_rig(NX - NX_CPML_rig:NX), a2_btm(0:NY_CPML_btm), &
      b2_lef(0:NX_CPML_lef), b2_rig(NX - NX_CPML_rig:NX), b2_btm(0:NY_CPML_btm)

    double precision :: &
      L_lef, L_rig, L_btm, &
      d0_lef, d0_rig, d0_btm, alph0, &
      d, alph

    integer :: i, j

    double precision, parameter :: PI = 3.1415926535897932d0

    L_lef = dble(NX_CPML_lef) * DELTAX
    L_rig = dble(NX_CPML_rig) * DELTAX
    L_btm = dble(NY_CPML_btm) * DELTAY


    d0_lef = - dble(NPOWER+1) * vp * log(reflection_rate) / 2.d0 / L_lef
    d0_rig = - dble(NPOWER+1) * vp * log(reflection_rate) / 2.d0 / L_rig
    d0_btm = - dble(NPOWER+1) * vp * log(reflection_rate) / 2.d0 / L_btm

    alph0 = PI * source_freq
!------------------the left side of x-axis--------------------------
    do i = 1,NX_CPML_lef

      j = NX_CPML_lef - i
      d = d0_lef * ( ( dble(j) + 0.5d0 ) / dble(NX_CPML_lef) ) ** NPOWER
      alph = alph0 * ( 1.d0 - ( dble(j) + 0.5d0 ) / dble(NX_CPML_lef) )

      a1_lef(i) = (2.d0 - (d + alph) * DELTAT) / (2.d0 + (d + alph) * DELTAT)

      b1_lef(i) = - d * DELTAT / (2.d0 + (d + alph) * DELTAT)

    enddo

    do i = 0,NX_CPML_lef

      j = NX_CPML_lef - i

      d = d0_lef * ( dble(j)/ dble(NX_CPML_lef) ) ** NPOWER
      alph = alph0 * ( 1.d0 - dble(j) / dble(NX_CPML_lef) )

      a2_lef(i) = (2.d0 - (d + alph) * DELTAT) / (2.d0 + (d + alph) * DELTAT)
      b2_lef(i) = - d * DELTAT / (2.d0 + (d + alph) * DELTAT)

    enddo
!------------------the right side of x-axis--------------------------
    do i = NX - NX_CPML_rig+1, NX

      j = i - NX + NX_CPML_rig

      d = d0_rig * ( ( dble(j) - 0.5d0 ) / dble(NX_CPML_rig) ) ** NPOWER
      alph = alph0 * ( 1.d0 - ( dble(j) - 0.5d0 ) / dble(NX_CPML_rig) )

      a1_rig(i) = (2.d0 - (d + alph) * DELTAT) / (2.d0 + (d + alph) * DELTAT)
      b1_rig(i) = - d * DELTAT / (2.d0 + (d + alph) * DELTAT)

    enddo

    do i = NX - NX_CPML_rig, NX

      j = i - NX + NX_CPML_rig

      d = d0_rig * ( dble(j)/ dble(NX_CPML_rig) ) ** NPOWER
      alph = alph0 * ( 1.d0 - dble(j) / dble(NX_CPML_rig) )

      a2_rig(i) = (2.d0 - (d + alph) * DELTAT) / (2.d0 + (d + alph) * DELTAT)
      b2_rig(i) = - d * DELTAT / (2.d0 + (d + alph) * DELTAT)

    enddo
!------------------the bottom of y-axis------------------------------
    do i = 1,NY_CPML_btm

      j = NY_CPML_btm - i

      d = d0_btm * ( ( dble(j) + 0.5d0 ) / dble(NY_CPML_btm) ) ** NPOWER
      alph = alph0 * ( 1.d0 - ( dble(j) + 0.5d0 ) / dble(NY_CPML_btm) )

      a1_btm(i) = (2.d0 - (d + alph) * DELTAT) / (2.d0 + (d + alph) * DELTAT)
      b1_btm(i) = - d * DELTAT / (2.d0 + (d + alph) * DELTAT)

    enddo

    do i = 0,NY_CPML_btm

      j = NY_CPML_btm - i

      d = d0_btm * ( dble(j)/ dble(NY_CPML_btm) ) ** NPOWER
      alph = alph0 * ( 1.d0 - dble(j) / dble(NY_CPML_btm) )

      a2_btm(i) = (2.d0 - (d + alph) * DELTAT) / (2.d0 + (d + alph) * DELTAT)
      b2_btm(i) = - d * DELTAT / (2.d0 + (d + alph) * DELTAT)

    enddo

  end subroutine

end module







