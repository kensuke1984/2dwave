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
    namelist /structure/modelfile
    double precision, allocatable :: &
        rho(:,:), C(:,:,:,:)

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

    !----------------parameters for output
    logical :: need_snapshot = .true.
    integer :: NT_snap = 100
    !should be four outputs (x y ux uy)
    character(12) :: snap_format = '(4G15.5)'
    logical :: need_waveform = .true.
    !should be three outputs (t ux uy)
    character(12) :: wave_format = '(3G15.5)'
    !number of recievers
    integer :: num_rec
    !recievers names and the locations (X,Y) are included  station
    character(len=32) :: recinfo = 'recinfo.csv'
    namelist/output/need_snapshot, NT_snap, snap_format, &
        need_waveform, wave_format, recinfo
     !output files in order
    character(len=100), allocatable :: outfname_rec(:)
    !position for receivers in order
    double precision, allocatable, dimension (:) :: X_rec, Y_rec
    !---------------------------------------------------------------------

    !----------------CPML settings-------------------------------------
    integer :: NX_CPML_lef, NX_CPML_rig, NY_CPML_btm
    double precision :: reflection_rate = 0.001d0
    !  double precision :: source_freq = 1.1283791671d0 / tp
    double precision :: source_freq = 1.1283791671d0 *6.5
    integer :: NPOWER = 2
    namelist /cpml/NX_CPML_lef,NX_CPML_rig, NY_CPML_btm,reflection_rate,source_freq, NPOWER
    !-----------CPML END--------------------------------------------------------

    ! arguments for the process (means parameter files).
    !  Each file name must be under 100 characters.
    !    character(:), allocatable :: arg
    character (100), allocatable :: args(:)
contains

    !#####################
    subroutine read_structure_csv
               !##########################
        double precision, allocatable, dimension (:,:) :: vp, vs

        double precision :: read_rho, read_vp, read_vs

        double precision :: mu, lambda, elas(3,3)

        integer :: i, j, ix, iy, io, read_count

        character (100) :: buffer

        allocate(rho(NX,NY), vp(NX,NY), vs(NX,NY))
        allocate(C(3,3,NX,NY))
        rho = - 80704.d0
        open (17, file=modelfile, status='old')
        read_count = 0
        do
            read (17,'(a)' , iostat=io) buffer
            buffer = adjustl(buffer)
            if(buffer(1:1) =='#' .or. buffer(1:1) == '!') cycle
            if (io/=0) exit
            read (buffer, *) ix, iy, read_rho, read_vp, read_vs
            rho(ix,iy)=read_rho
            vp(ix,iy)=read_vp
            vs(ix,iy)=read_vs
            read_count = read_count+1
        enddo
        close(17)

        if (read_count /= NX*NY) stop 'this model file (csv) is invalid.'

        do ix = 1, NX
            do iy = 1, NY
                !--check if csv is correctly inputted----
                if ( rho(ix,iy) < 0.d0 .or. &
                    rho(1,iy) <= 1.d-10 .or. &
                    rho(NX,iy) <= 1.d-10 .or. &
                    rho(ix,1) <= 1.d-10 ) stop 'this model is invalid. (containing impossible values)'

                mu = rho(ix,iy) * vs(ix,iy) * vs(ix,iy)
                lambda = rho(ix,iy) * ( vp(ix,iy) * vp(ix,iy) - 2.d0 * vs(ix,iy) * vs(ix,iy) )

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

                C(:,:,ix,iy) = elas(:,:)

            enddo
        enddo

    end subroutine read_structure_csv

    subroutine read_parameter(i)
        integer, intent (in) :: i
        print *, 'Reading the input file: ',args(i)
        open (1, file=args(i), status='old')
        read (1, structure)
        read (1, geometry)
        read (1, source)
        read (1, computation)
        read (1, cerjan)
        read (1, output)
        read (1, cpml)
        close (1)
        call read_structure_csv
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
        integer :: sta, i
        intrinsic :: command_argument_count, get_command_argument
        !        if (command_argument_count() > 1) stop 'Only one argument (input file name) is available'
        if (command_argument_count() == 0) stop 'Input file name is required'
        allocate(args(command_argument_count()))

        sta = 0
        do i = 1, size(args)
            !        call get_command_argument(i, length = length, status = status)
            call get_command_argument(i, args(i), status = sta)
            if (sta /= 0) stop 'Error on argument'
        enddo
    !        write(*,*) 'Input file: ',(args(i),', ',i=1,size(args))
    end subroutine get_arguments

    subroutine set_receiver_info
        character(1024) buffer
        integer :: iline, io, name_length

        !    X_rec = - 80704.d0

        open (17, file=recinfo, status='old')
        num_rec =0
        !    if ( header ) read (17, '()')  !skip header
        do
            read(17, '(a)', iostat=io) buffer
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
            name_length = index(buffer, ",")
            buffer = trim(buffer)
            outfname_rec(iline) = buffer(1:name_length-1)
            read (buffer(name_length+1:),*) X_rec(iline), Y_rec(iline)
            outfname_rec(iline) = trim(adjustl(outfname_rec(iline)))
            if ( X_rec(iline) < 0.d0 .or. Y_rec(iline) <0) stop 'reciever info file is invalid'
            iline = iline +1
        enddo
        close (17)
    end subroutine set_receiver_info
end module parameter_2dwave

