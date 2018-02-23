module parameter_copt_run

  implicit none
   
! model
  character(len=32), parameter :: model_csv = 'kmodel.csv'
  logical, parameter :: header_model = .true.

! size of region
  double precision, parameter :: X_size = 11118.d0
  double precision, parameter :: Y_size = 4.d0

! source location
  double precision, parameter :: X_source = 1.0d0
  double precision, parameter :: Y_source = 3.9d0

!source
  logical, parameter :: source_type_single = .false.
  double precision, parameter :: vx = 1.d0
  double precision, parameter :: vy = 0.d0

  logical, parameter :: source_type_double_couple = .true.
  double precision, parameter :: mxx = 1.d0
  double precision, parameter :: mxy = 0.d0
  double precision, parameter :: myx = 0.d0
  double precision, parameter :: myy = 1.d0

!time length [s]
  double precision, parameter :: tlen = 10.d0
!central period of ricker function [s]
  double precision, parameter :: tp = 1.d0/ 6.5d0
!shift time [s]
  double precision, parameter :: ts = 2.d0 * tp

!number of boxes in the region
  integer, parameter :: NX = 600
  integer, parameter :: NY = 300

  logical, parameter :: use_Cerjan = .false.
  integer, parameter :: N_Cerjan = 5
  double precision, parameter :: dumping_rate_Cerjan = 0.0001d0

  double precision, parameter :: CFL = 0.3d0
  double precision, parameter :: maxvp = 2.8d0

  logical, parameter :: need_snapshot = .true.
  integer, parameter :: NT_snap = 100


  logical, parameter :: need_waveform = .true.
!number of recievers
  integer, parameter :: num_rec = 59
!recievers names and the locations (X,Y) are included  station
  character(len=32), parameter :: recinfo = 'recinfo.csv'
  logical, parameter :: header_recinfo = .false.

!----------------CPML settings-------------------------------------
  integer, parameter :: NX_CPML_lef = 10
  integer, parameter :: NX_CPML_rig = 10
  integer, parameter :: NY_CPML_btm = 10

  double precision, parameter :: reflection_rate = 0.001d0

  double precision, parameter :: source_freq = 1.1283791671d0 / tp

  integer, parameter :: NPOWER = 2

!-----------END--------------------------------------------------------

end module
