module main_subs

  implicit none

contains

subroutine read_main_parfile_mpi(rank, comm, ierr)

  use var_main

  include 'mpif.h'
  integer :: rank, comm, ierr
  
  if(rank.eq.0)then
    print *, "Read in master node:"
    call read_main_parfile(ierr)
  endif

  print *,"Bcast the par..."
  call MPI_Bcast(WRITE_ADJ_ASDF,1,MPI_LOGICAL,0,comm,ierr)
  call MPI_Bcast(ROTATE_COMP,1,MPI_LOGICAL,0,comm,ierr)
  call MPI_Bcast(WRITE_NORMAL_OUTPUT,1,MPI_LOGICAL,0,comm,ierr)

  call MPI_Bcast(OBSD_FILE,150,MPI_CHARACTER,0,comm,ierr)
  call MPI_Bcast(SYNT_FILE,150,MPI_CHARACTER,0,comm,ierr)
  call MPI_Bcast(FLEXWIN_OUTDIR,150,MPI_CHARACTER,0,comm,ierr)
  call MPI_Bcast(MEASURE_ADJ_OUTDIR,150,MPI_CHARACTER,0,comm,ierr)

  call MPI_Bcast(WEIGHTING_OPTION,1,MPI_INTEGER,0,comm,ierr)
	!if(rank.eq.1) then
	!	print *, "MPI_staff"
	!	print *, RUN_FLEXWIN, RUN_MEASURE_ADJ, WRITE_ADJ_ASDF,&
	!			ROTATE_COMP, WRITE_NORMAL_OUTPUT
	!		print *, trim(OBSD_FILE), 	
	!   PRINT *, trim(MEASURE_ADJ_OUTDIR)
	!endif

end subroutine read_main_parfile_mpi


subroutine read_main_parfile(ierr)

  !read the parfile for the main(some flags)
  use var_main

  integer :: dummy_row
  integer :: ierr
  integer :: IIN=21
  integer :: i

  character(len=30) :: dummy_string

  !print *,"Read main par"
  dummy_row = 8
  
  open(UNIT=IIN,FILE="PAR_FILE_MAIN",iostat=ierr)
  if(ierr.ne.0)then
    print *,"Can't find PAR_FILE_MAIN. Stop! "
    stop
  endif

  do i=1,dummy_row
    read(IIN,*)
  enddo

  !print *,"HERE"

  read(IIN,3) dummy_string, WRITE_ADJ_ASDF
  print *,"WRITE_ADJ_ASDF: ", WRITE_ADJ_ASDF
  read(IIN,3) dummy_string, ROTATE_COMP 
  print *,"ROTATE_COMP: ", ROTATE_COMP
  read(IIN,3) dummy_string, WRITE_NORMAL_OUTPUT 
  print *,"WRITE_NORMAL_OUTPUT",WRITE_NORMAL_OUTPUT

	do i=1,2
		read(IIN,*)
	enddo

	read(IIN,2) dummy_string, OBSD_FILE
	print *, "OBSD_FILE: ", trim(OBSD_FILE)
	read(IIN,2) dummy_string, SYNT_FILE
	print *, "SYNT_FILE: ", trim(SYNT_FILE)
	read(IIN,2) dummy_string, FLEXWIN_OUTDIR
	print *, "FLEXWIN_OUTDIR: ", trim(FLEXWIN_OUTDIR)
	read(IIN,2) dummy_string, MEASURE_ADJ_OUTDIR
	print *, "MEASURE_ADJ_OUTDIR: ", trim(MEASURE_ADJ_OUTDIR)

  do i=1,2
    read(IIN,*)
  enddo

  read(IIN,4) dummy_string, weighting_option
  print *, "weighting_option: ", weighting_option

2 format(a,a)
3 format(a,l20)
4 format(a,i)

  close(IIN)
	!stop

end subroutine read_main_parfile

subroutine read_ma_parfile_mpi(ma_par_all,min_period,&
              max_period,event_dpt, rank, comm, ierr)

	use ma_struct
  use measure_adj_subs

  include 'mpif.h'

	type(ma_par_struct_all) :: ma_par_all
	real :: min_period, max_period, event_dpt
	integer :: rank, comm, ierr

  integer, parameter :: NDIM_PAR=3
  integer, parameter :: BLOCK_PER_DIM=4
  type(ma_par_struct) :: ma_par_temp(NDIM_PAR)

	integer :: oldtype(BLOCK_PER_DIM), newtype
  integer :: offset(BLOCK_PER_DIM), blockcount(BLOCK_PER_DIM)
	integer :: extent

	integer :: tag=1, i, source, loc
	integer :: stat(MPI_STATUS_SIZE)

	!print *,"SET UP"
	!setup description of the flexwin_par
	!call read_flexwin_parfile_mpi(flexwin_par, fstart, fend, rank, nproc, comm)

  blockcount = (/5,5,19,10/)
  !LOGICAL
  offset(1) = 0
	oldtype(1) = MPI_LOGICAL

  !INTEGER
	call MPI_TYPE_EXTENT(MPI_LOGICAL, extent, ierr)
	offset(2) = offset(1)+blockcount(1)*extent
	oldtype(2) = MPI_INTEGER

  !DOUBLE_PRECISION
	call MPI_TYPE_EXTENT(MPI_INTEGER,extent, ierr)
	offset(3) = offset(2)+blockcount(2)*extent
	oldtype(3) = MPI_DOUBLE_PRECISION

  !CHARACTER
	call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent, ierr)
	offset(4) = offset(3)+blockcount(3)*extent
	oldtype(4) = MPI_CHARACTER
	!if(rank.eq.0) then
	!	print *,"blockcount",loc,blockcount(loc)
	!endif
	!if(rank.eq.1) then
	!	print *,"blockcount:",blockcount(loc)
	!endif

	!print *,i, blockcount(1), blockcount(2), blockcount(3), blockcount(4),&
	!					blockcount(5), blockcount(6)

	!now define and commit
	call MPI_TYPE_STRUCT(BLOCK_PER_DIM, blockcount, offset, oldtype, newtype, ierr)
	call MPI_TYPE_COMMIT(newtype, ierr)

	print *, "SET UP finished!"
	if(rank.eq.0) then
		call read_ma_parfile(ma_par_all,min_period,max_period,event_dpt)
    ma_par_temp(1)=ma_par_all%R
    ma_par_temp(2)=ma_par_all%T
    ma_par_temp(3)=ma_par_all%Z
	endif

  call MPI_Bcast(ma_par_temp, 3, newtype, 0, comm, ierr)

	!if(rank==1) then
	!	print *,"HERE check"
	!	print *, ma_par_temp(:)%TLONG, ma_par_temp(:)%TSHORT
	!endif
	!call MPI_Barrier(comm, ierr)

  ma_par_all%R=ma_par_temp(1)
  ma_par_all%T=ma_par_temp(2)
  ma_par_all%Z=ma_par_temp(3)

	!if(rank==1) then
		!print *, "CHECK"
		!print *, ma_par_all%Z%WTR
		!print *, ma_par_all%Z%TLONG, ma_par_all%Z%TSHORT
		!print *, trim(ma_par_all%Z%chan)
		!print *, ma_par_all%Z%RUN_BANDPASS
	!endif
  print *, "finalize"

	!stop

end subroutine read_ma_parfile_mpi

subroutine write_win_chi(MEASURE_ADJ_OUTDIR, nrecords, &
              event_name, p1, p2, sta, net,&
              chan_syn, win_chi_all, win_all,&
              rank, ierr)

  use ma_struct
  use flexwin_struct
  implicit none
  
  type(win_chi_info),dimension(:),intent(in) :: win_chi_all
  type(win_info),dimension(:), intent(in) :: win_all
  character(len=*) :: event_name
  real :: p1, p2
  character(len=*) :: sta(:),net(:),chan_syn(:)
  integer :: nrecords
  character(len=*) :: MEASURE_ADJ_OUTDIR
  integer :: rank, ierr

  integer :: IIN=110
  integer :: i,j,k

  character(len=32) :: file_prefix 
  character(len=250) :: fn, p1_string, p2_string, myid_string

  write(myid_string, '(I8)') rank
  write(p1_string, '(I8)') int(p1)
  write(p2_string, '(I8)') int(p2)
  myid_string=adjustl(myid_string)
  p1_string=adjustl(p1_string)
  p2_string=adjustl(p2_string)
  
  call system('mkdir -p '//trim(MEASURE_ADJ_OUTDIR)//'/'//trim(event_name)//'')
  fn=trim(MEASURE_ADJ_OUTDIR)//'/'//trim(event_name)//'/'//&
      trim(event_name)//'_'//trim(p1_string)//'_'//&
      trim(p2_string)//'.'//trim(myid_string)//'.winchi'
  open(UNIT=IIN,file=fn)

  do i=1,nrecords
    print *, "nrecords:", nrecords
    !print *,"i",i
    print *, win_all(i)%num_win
    do j=1,win_all(i)%num_win
      !print *, "irecords, num of window:", i,j
      !print *, "sta:", trim(sta(i))
      !print *, "net:", trim(net(i))
      !print *, "comp:", trim(chan_syn(i))
      file_prefix=trim(sta(i))//"."//trim(net(i))//"."//trim(chan_syn(i))
      print *,"file_prefix: ", file_prefix
      write(IIN,'(a14,a8,a3,a5,i4,i4,2e14.6,20e14.6,2e14.6,2f14.6)') &
           trim(file_prefix),trim(sta(i)),trim(net(i)),trim(chan_syn(i)),j,&
           win_chi_all(i)%imeas(j),&
           win_all(i)%t_start(j),win_all(i)%t_end(j),&
           (win_chi_all(i)%chi(j,k),k=1,20),&
           win_chi_all(i)%tr_chi(j),win_chi_all(i)%am_chi(j),&
           win_chi_all(i)%T_pmax_dat(j),win_chi_all(i)%T_pmax_syn(j)
      print *, '   tr_chi = ', sngl(win_chi_all(i)%tr_chi(j)),&
               '   am_chi = ', sngl(win_chi_all(i)%am_chi(j))
    enddo
  enddo

end subroutine write_win_chi


subroutine write_ascii_output(my_asdf, outdir)

  use asdf_data
  use ascii_rw

  type(asdf_event) :: my_asdf
  character(len=150) :: outdir

  double precision, allocatable :: data(:)
  double precision :: b, dt
  integer :: npt

  integer :: i, j

  character(len=300) :: fn, file_prefix

	!do a channel name modify here
	do i=1,my_asdf%nrecords
		my_asdf%component_array(i)(1:2)="LH"
	enddo

  do i=1,my_asdf%nrecords
    file_prefix=trim(my_asdf%receiver_name_array(i))//"."//&
          trim(my_asdf%network_array(i))//"."//&
          trim(my_asdf%component_array(i))
    fn=trim(outdir)//"/"//trim(file_prefix)//".adj"
    print *, "fn:", trim(fn)

    allocate(data(my_asdf%npoints(i)))
    data(:)=dble(my_asdf%records(i)%record)
    b=dble(my_asdf%begin_value(i))
    dt=dble(my_asdf%sample_rate(i))
    npt=my_asdf%npoints(i)
    call dwascii(fn, data, npt, b, dt)
    deallocate(data)
  enddo

end subroutine write_ascii_output

subroutine copy_general_info_to_adj(obsd, adj)
  
  use asdf_data
  implicit none

  type(asdf_event) :: obsd, adj
  integer :: i

  adj%event_lat=obsd%event_lat
  adj%event_lo=obsd%event_lo
  adj%event_dpt=obsd%event_dpt
  adj%event=obsd%event

  do i=1,obsd%nrecords
    adj%receiver_lat(:)=obsd%receiver_lat(:)
    adj%receiver_lo(:)=obsd%receiver_lo(:)
    adj%scale_factor(:)=obsd%scale_factor(:)
    adj%receiver_name_array(:)=obsd%receiver_name_array(:)
    adj%network_array(:)=obsd%network_array(:)
    adj%component_array(:)=obsd%component_array(:)
    !adj%=obsd%
    !adj%=obsd%
    !adj%=obsd%
  enddo

end subroutine copy_general_info_to_adj

end module main_subs
