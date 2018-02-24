! ============================================================================
! Name        : parameter_copt_run.f90
! Author      : Kei Hasegawa, Kensuke Konishi
! Version     : 0.0.1
! Copyright   : It is Complicated.
! Description : Parameter set for 2dwave
! ============================================================================

module parameter_2dwave
  
  implicit none 
! model
  character(len=32) :: modelfile = 'kmodel.csv'
  logical :: header_model = .true.
  namelist /structure/modelfile, header_model

! size of region
  double precision :: X_size = 8.d0
  double precision :: Y_size = 4.d0
! number of boxes in the region
  integer :: NX = 600
  integer :: NY = 300
  namelist/geometry/X_size, Y_size, NX, NY

! Information for source
  !location
  double precision :: X_source = 1.0d0
  double precision :: Y_source = 3.9d0

  logical :: source_type_single = .false.
  double precision :: vx = 1.d0
  double precision :: vy = 0.d0

  logical :: source_type_double_couple = .true.
  double precision :: mxx = 1.d0
  double precision :: mxy = 0.d0
  double precision :: myx = 0.d0
  double precision :: myy = 1.d0
  namelist/source/X_source, Y_source, source_type_single, source_type_double_couple, &
    vx,vy, mxx, mxy, myx, myy

!time length [s]
  double precision :: tlen = 10.d0
!central period of ricker function [s]
  double precision :: tp = 1.d0/ 6.5d0
!shift time [s]
!  double precision :: ts = 2.d0 * tp
  double precision :: ts = 2.d0 * 1 /6.5d0
  double precision :: CFL = 0.3d0
  double precision :: maxvp = 2.8d0
  namelist/computation/ tlen, tp, ts, CFL, maxvp

  logical :: use_Cerjan = .false.
  integer :: N_Cerjan = 5
  double precision :: dumping_rate_Cerjan = 0.0001d0
  namelist/cerjan/use_Cerjan, N_Cerjan, dumping_rate_Cerjan

  character(:), allocatable :: arg

  logical :: need_snapshot = .true.
  integer :: NT_snap = 100
  logical :: need_waveform = .true.
!number of recievers
  integer :: num_rec = 59
!recievers names and the locations (X,Y) are included  station
  character(len=32) :: recinfo = 'recinfo.csv'
  logical :: header_recinfo = .false.
  namelist/output/need_snapshot, NT_snap, need_waveform, num_rec, recinfo, header_recinfo

!----------------CPML settings-------------------------------------
  integer :: NX_CPML_lef, NX_CPML_rig, NY_CPML_btm
  double precision :: reflection_rate = 0.001d0
!  double precision :: source_freq = 1.1283791671d0 / tp
  double precision :: source_freq = 1.1283791671d0 *6.5
  integer :: NPOWER = 2
  namelist /cpml/NX_CPML_lef,NX_CPML_rig, NY_CPML_btm,reflection_rate,source_freq, NPOWER
!-----------CPML END--------------------------------------------------------
  contains

    subroutine read_parameter
      implicit none
      print *, 'Reading the input file'
      open (1, file=arg, status='old')
      read (1, structure)
      read (1, geometry)
      read (1, source)
      read (1, computation)
      read (1, cerjan)
      read (1, output)
      read (1, cpml)
      close (1)
    end subroutine read_parameter

    subroutine output_parameter
       write(*,structure)
       write(*,geometry)
       write(*,source)
       write(*,computation)
       write(*,cerjan)
       write(*,output)
       write(*,cpml)
    end subroutine output_parameter

    subroutine get_arguments
      integer :: i, length, status, nargs
      intrinsic :: command_argument_count, get_command_argument
      if (command_argument_count() > 1) stop 'Only one argument (input file name) is available'
      if (command_argument_count() == 0) stop 'Input file name is required'

      call get_command_argument(1, length = length, status = status)
      if (status == 0) then
        allocate (character(length) :: arg)
        call get_command_argument(1, arg, status = status)
        if (status == 0) print *, 'Input file is "', arg, '"'
!        deallocate (arg)
      end if
      if (status /= 0) stop 'Error on argument'
    end subroutine get_arguments


end module parameter_2dwave

