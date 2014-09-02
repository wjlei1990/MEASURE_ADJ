!==============================================================================
!Program: Global_Tomography_Data_Processing
!Developer: Princeton Global Tomography Group(PGTG)
!Group Member: Wenjie Lei(lei@princeton.edu), Ebru Bozdag(bozdag@princeton.edu),
!James A. Smith(jas11@princeton.edu)
!===================
!Functions:
!1) Read ADIOS
!2) Processing: rtrend, rmean, taper, inst_remove(filter)
!3) data_quality: generate the useful data list
!4) flexwin: window selection
!5) Measure_adj: make measurements and adjoint sources
!5) write out the new ADIOS file(filtered, selected and window selected)
!===================
!Bug Report: lei@princeton.edu
!===============================================================================

program main

  use asdf_data
  use flexwin_struct
  use ma_struct

  use win_io_subs
  use asdf_read_subs
  use asdf_write_subs

  use ma_interface
  use rotate_subs
  use mpi_weighting_subs

  use var_main
  use main_subs

  implicit none
  include 'mpif.h'

  type(asdf_event)        :: synt_all, obsd_all, synt_phydisp_all
  type(asdf_event)        :: adj_all, adj_all_rotate
  type(win_info),allocatable      :: win_all(:)
  type(win_chi_info), allocatable :: win_chi_all(:)
  type(ma_par_struct_all)         :: measure_adj_par_all
  type(ma_weighting_par_struct)   :: ma_weighting_par

  character(len=200) :: ADJ_FILE

  integer :: nrecords
  character(len=20) :: station(MAXDATA_PER_PROC), network(MAXDATA_PER_PROC)
  character(len=20) :: component(MAXDATA_PER_PROC), receiver_id(MAXDATA_PER_PROC)
  character(len=150) :: ma_outdir

  double precision, allocatable :: adj_source(:)
  !mpi_var
  integer                 :: nproc,comm,rank
  integer                 :: ierr, adios_err
  !adios_var
  integer(kind=8)         :: adios_group

  integer                 :: i

  double precision :: t1, t2

  character(len=20) :: p1_string, p2_string

  t1=MPI_WTIME()
  !----------.
  !init mpi  !
  !----------'
  call mpi_init(ierr)
  call mpi_comm_dup(mpi_comm_world,comm,ierr)
  call mpi_comm_rank(comm,rank,ierr)
  call mpi_comm_size(comm,nproc,ierr)
  if(rank.eq.0) then
    print *, "Start Measure_Adj..."
    print *, "NPROC:", nproc
  endif

  !init adios
  call adios_init_noxml(comm, adios_err)
  call adios_allocate_buffer(600, adios_err)
  call adios_declare_group(adios_group, "EVENTS", "iter", 1, adios_err)
  call adios_select_method(adios_group, "MPI", "", "", adios_err)

  !--------------------------.
  !read main parfile         !
  !--------------------------'
  if(rank.eq.0)then
    print *, "-----------------"
    print *,"Read in main Parfile..."
    print *, "-----------------"
  endif
  call read_main_parfile_mpi(rank,comm,ierr)
  !stop

  !--------------------------.
  !read in asdf data         !
  !--------------------------'
  if(rank.eq.0) then
    print *, "-----------------"
    print *, "Read in file"
    print *, "OBSD_FILE: ",trim(OBSD_FILE)
    print *, "SYNT_FILE: ",trim(SYNT_FILE)
    print *, "-----------------"
  endif
  call read_asdf_file(OBSD_FILE, obsd_all, nrecords, &
    station, network, component, receiver_id, 0, &
    rank, nproc, comm, ierr)
  print *, "read obsd finished!"
  call read_asdf_file(SYNT_FILE, synt_all, nrecords, &
    station, network, component, receiver_id, 1, &
    rank, nproc, comm, ierr)
  print *, "read synt finished!"
  if(USE_PHYDISP)then
    !if use-phydisp, then read in phydisp file
    print *, "read phydisp!"
    call read_asdf_file(SYNT_PHYDISP_FILE, synt_phydisp_all, nrecords, &
    station, network, component, receiver_id, 1, &
    rank, nproc, comm, ierr)
  endif
  if(rank.eq.0) then
    print *, "/event:", trim(obsd_all%event)
    print *, "/nrecords:", obsd_all%nrecords
  endif

  !should be removed in the future
  obsd_all%min_period = MIN_PERIOD
  obsd_all%max_period = MAX_PERIOD

  allocate(win_all(obsd_all%nrecords))
  allocate(win_chi_all(obsd_all%nrecords))

  !----------------------------------------.
  !read in the window information           
  !----------------------------------------'
  if(rank.eq.0)then
    print *,"----------------------------------"
    print *,"READ WIN FILE                     "
    print *,"----------------------------------"
  endif
  call win_read(WIN_DIR, obsd_all%event, obsd_all%min_period,&
                obsd_all%max_period, win_all,&
                obsd_all%nrecords, rank, ierr)
 
  !----------------------------------------------------------------
  !measure_adj interface for fortran program
  !----------------------------------------------------------------
  call measure_adj_interface(obsd_all, synt_all, synt_phydisp_all, &
          win_all, adj_all, adios_group, rank, nproc, comm, ierr)

  !--------------------------.
  !finalize mpi              !
  !--------------------------'
  call MPI_Barrier(comm,ierr)
  call adios_finalize(rank, ierr)
  call mpi_finalize(ierr)

end program main
