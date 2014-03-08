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
!5) write out the new ADIOS file(filtered, selected and window selected)
!===================
!Bug Report: lei@princeton.edu
!===============================================================================

!-------------------------------------------------------------
!main: test the measure_adj asdf version
program main

  use asdf_data
  use flexwin_struct
  use ma_struct
  use win_io
  use asdf_subs
  use measure_adj_subs
  use adios_write_mod
  use adios_read_mod
  !use asdf_subs
	!use flexwin_subs
  implicit none
  include 'mpif.h'

  !this is the input file name
  !character(len=120) :: input_fn='input_cmt_LHZ_17_60'
  character(len=150)      :: OBSD_FILE, SYNT_FILE, SYNT_PHYDISP_FILE, ADJ_FILE
  character(len=150)      :: WIN_FILE
  !this is the output directory
  character(len=150)      :: OUTDIR
  !just for test use
  character(len=150)      :: OUTFN
  !variables
  type(asdf_event)        :: synt_all,synt_phydisp_all,obsd_all,adj_all
  type(win_info),allocatable      :: win_all(:)
  type(win_chi_info), allocatable :: win_chi_all(:)
  type(ma_par_struct_all) :: measure_adj_par_all
  type(ma_weighting_par_struct) :: ma_weighting_par
  !type(ma_)
  double precision        :: all_chi
  real, allocatable :: adj_source(:)
  !mpi_var
  integer                 :: nproc,comm,rank
  integer                 :: ierr, adios_err
  integer                 :: i, irecord
  !adios_var
  integer(kind=8)         :: adios_groupsize, adios_totalsize, varid
  integer(kind=8)         :: adios_handle, adios_group

  !init mpi
  call mpi_init(ierr)
  call mpi_comm_dup(mpi_comm_world,comm,ierr)
  call mpi_comm_rank(comm,rank,ierr)
  call mpi_comm_size(comm,nproc,ierr)
  !split the job for every processor
  !call adios_init_noxml(comm)
  !call adios_allocate_buffer(100, adios_err)
  !call adios_declare_group(adios_group,"EVENTS","iter",1,adios_err)
  !call adios_select_method(adios_group,"MPI","","",adios_err)

  OBSD_FILE='./DATA/test/200801151752A_obsd.bp'
  SYNT_FILE='./DATA/test/200801151752A_synt.bp'
  SYNT_PHYDISP_FILE=''
  ADJ_FILE ='./output/200801151752A_adj.bp'
  WIN_FILE='./200801151752A_1.win.mat'
  OUTDIR='./output'

  print *,"Read ma parfile Begin!"
  call read_ma_parfile(measure_adj_par_all, 17.0, 60.0)
  !==============================================
  !read in the obsd and synt
  !call ADIOS_read(obsd_all,synt_all,input_fn)
  call read_asdf_file(OBSD_FILE,obsd_all,comm)
  print *, "read obsd finished!"
  !print *,"Reading OBSD finished"
  call read_asdf_file(SYNT_FILE,synt_all,comm)
  print *, "read synt finished!"

  !stop
  !==============================================

  print *, "number of records: ", obsd_all%nrecords
  allocate(win_all(obsd_all%nrecords))
  allocate(win_chi_all(obsd_all%nrecords))

  !==============================================
  !FLEXWIN
  !call flexwin(obsd_all,synt_all,win,myid,numprocs,comm)
  !==============================================

  !==============================================
  !read in window
  call win_read(WIN_FILE, win_all, obsd_all%nrecords)
  !==============================================

  !OUTFN='check.win'
  !call win_write_demo(OUTFN,obsd_all%nrecords,win_all)

  !==============================================
  !start flexwin part
  !call read_flexwin_parfile(flexwin_par_all, obsd_all%min_period, &
  !        obsd_all%max_period)

  !do i=1, obsd_all%nrecords
    !call flexwin subroutine
    !call flexwin(obsd_all%records(i)%record,obsd_all%npoints(i),obsd_all%sample_rate(i),obsd_all%begin_value(i),&
    !synt_all%records(i)%record, synt_all%npoints(i),synt_all%sample_rate(i),synt_all%begin_value(i),&
    !obsd_all%event_lat, obsd_all%event_lo, obsd_all%event_dpt,obsd_all%receiver_lat(i),obsd_all%receiver_lo(i),&
    !obsd_all%receiver_name,obsd_all%network,obsd_all%component,obsd_all%P_pick(i),obsd_all%S_pick(i),&
    !flexwin_par_all,win_all(i))
  !enddo

  !==============================================
  !start measure_adj part
  !call read_ma_parfile(measure_adj_par_all,obsd_all%min_period,&
  !                          obsd_all%max_period)
  print *,"Read ma parfile finished!"
  print *,"Weighting Begin!"
  call setup_measure_adj_weighting_asdf(win_all,obsd_all%nrecords, &
          obsd_all%great_circle_dist, obsd_all%component_array, &
          ma_weighting_par)
  print *, "Weighting finished!"

  allocate(adj_source(measure_adj_par_all%Z%nn))

  print *, "GO into the measure_adj"
  do i=1, obsd_all%nrecords
    !call measure_adj subroutine
    call measure_adj(obsd_all%records(i)%record,obsd_all%npoints(i),obsd_all%begin_value(i),obsd_all%sample_rate(i),&
      synt_all%records(i)%record,synt_all%npoints(i),synt_all%begin_value(i),synt_all%sample_rate(i),&
      obsd_all%great_circle_dist(i),obsd_all%receiver_name_array(i),obsd_all%network_array(i),obsd_all%component_array(i),&
      win_all(i),measure_adj_par_all, ma_weighting_par, &
      win_chi_all(i), adj_source)

      !adj_all%records(i)%record=adj_source(1:measure_adj_par%Z%nn)
  end do
  !==============================================


  !==============================================
  !write out
  !>begin write out the adj_all
  !adios_groupsize = 0.
  !call define_asdf_data(adios_group, adios_groupsize,&
  !                      adj_all%nrecords, adj_all%npoints)

  !write out the adj_all to "****_adj.bp" file
  !call adios_open(adios_handle,"EVENTS",ADJ_FILE,"w",comm,adios_err)
  
  !call adios_group_size (adios_handle, adios_groupsize, adios_totalsize,adios_err)
  !call write_asdf_file (ADJ_FILE,adj_all,adios_handle,adios_group,&
  !                        adios_groupsize,comm)
  !call adios_close(adios_handle, adios_err)
  !==============================================


  !finalize mpi
  call MPI_Barrier(comm,ierr)
  call adios_finalize(rank,adios_err)
  call mpi_finalize(ierr)

end program main
