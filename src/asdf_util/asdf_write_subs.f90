module asdf_write_subs

!> The module asdf_subs contains Parallel ASDF I/O API:
!! 1)define_asdf_data
!! 3)write_asdf_file
  implicit none

contains

!! \param nreceivers The number of receivers
!! \param adios_group adios group
!! \param my_group_size Stores the adios group size

subroutine define_asdf_data (adios_group, my_group_size, my_asdf, &
								rank, nproc, comm, ierr)

  use adios_write_mod
	use adios_helpers_mod
	use asdf_data
  !use adios_read_mod
  implicit none

  integer(kind=8), intent(in) :: adios_group
	integer(kind=8) :: my_group_size
	type(asdf_event), intent(in) :: my_asdf
	integer, intent(in) :: rank, nproc, comm
	integer :: ierr

  integer :: i, nerr, string_total_length
  integer, parameter :: STRING_COMMON_LENGTH = 20
  integer :: adios_err, stat

  integer(kind=8) :: varid

  integer :: nrecords

  !character                    :: data_type, dummy_blank
  character(len=2)             :: data_type
  character(len=32)            :: header, record
  character(len=6)             :: npts_string
  character(len=10)            :: i_string
  character(len=200)           :: command, dummy, record_path

  integer :: dum_int, int_array(10)
  real    :: dum_real, real_array(10)
  character(len=10) :: dum_string

	integer :: nrecords_total, offset
	!gather info. Here, we only need nrecords_total
	nrecords=my_asdf%nrecords
	call gather_offset_info(nrecords,nrecords_total,offset,&
					rank, nproc, comm, ierr)

  call define_adios_local_string_1d_array (adios_group, my_group_size, &
												13, "", "event", dummy)
	!print *,"here", nrecords

  !nrecords info
  call define_adios_scalar (adios_group, my_group_size, "", "nreceivers",&
                        dum_int)
  call define_adios_scalar (adios_group, my_group_size, "", "nrecords",&
                        dum_int)
  !frequency(period) info
  call define_adios_scalar (adios_group, my_group_size, "", "min_period", &
                        dum_real)
  call define_adios_scalar (adios_group, my_group_size, "", "max_period", &
                        dum_real)

  !string info
  call define_adios_scalar (adios_group, my_group_size, "", "receiver_name_len", &
                        dum_int)
  call define_adios_scalar (adios_group, my_group_size, "", "network_len", &
                        dum_int)
  call define_adios_scalar (adios_group, my_group_size, "", "receiver_id_len", &
                        dum_int)
  call define_adios_scalar (adios_group, my_group_size, "", "component_len", &
                        dum_int)

  !print *, "TAG"
  !HEADER info
  open(5, file="./src/asdf_util/ASDF_HEADERS", iostat=stat, status='old')
  if(stat.ne.0)then
    print *,"Can not find ASDF_HEADERS"
    print *,"The default path: ./src/asdf_util/ASDF_HEADERS"
    print *,"Quit!"
    stop
  endif

  do
    read (5, *, iostat=stat) data_type, header
    if (stat /= 0) exit
    !print *,"=="
    !print *,trim(data_type), len_trim(data_type), trim(header), len_trim(header)
    select case (data_type(1:1))
      case ("i")
        call define_adios_global_integer_1d_array (adios_group, my_group_size,&
                    nrecords, "", trim(header), int_array)
      case ("r")
        call define_adios_global_real_1d_array (adios_group, my_group_size, &
                    nrecords, "", trim(header), real_array)
      case ("s")
        !Needs to pay attention in the future...here
        !potential bugs
        string_total_length = STRING_COMMON_LENGTH * nrecords_total
        call define_adios_local_string_1d_array (adios_group, my_group_size,&
                    string_total_length, "", trim(header), dum_string)
    end select
  enddo
  close(5)

  !DISPLACEMENT
  do i = 1, nrecords
	  !print *,"begin define:",i, my_asdf%npoints(i)
    write(i_string, '(I10)' ) i+offset
		record=trim(my_asdf%receiver_name_array(i))//"."//&
						trim(my_asdf%network_array(i))//"."//&
						trim(my_asdf%component_array(i))//"."//&
						trim(my_asdf%receiver_id_array(i))
    call define_adios_global_real_1d_array (adios_group, my_group_size,&
          my_asdf%npoints(i), "", trim(record),&
          real_array)
  enddo

  !define attribute
  call adios_define_attribute ( adios_group , "nreceivers", "desc", &
        adios_string, "Number of receivers ", "" , adios_err )
  call adios_define_attribute ( adios_group , "nrecords", "desc", &
        adios_string, "Number of records ", "" , adios_err ) 
  call adios_define_attribute ( adios_group , "min_period", "desc", &
        adios_string, "Low pass filter in Hz (0 if none applied)  ", "" , adios_err )
  call adios_define_attribute ( adios_group , "max_period", "desc", &
        adios_string, "High pass filter in Hz (0 if none applied)  ", "" , adios_err )
  call adios_define_attribute ( adios_group , "event_lat", "desc", adios_string, &
        "Event CMT latitude (degrees, north positive) ", "", adios_err )
  call adios_define_attribute ( adios_group , "event_lo", "desc", adios_string, &
        "Event CMT longitude (degrees, east positive) ", "", adios_err )
  call adios_define_attribute ( adios_group , "event_dpt", "desc", adios_string, &
        "Event CMT depth (km) ", "" , adios_err )
  call adios_define_attribute ( adios_group , "event_dpt", "desc", adios_string, &
        "Event CMT depth (km) ", "" , adios_err )
  call adios_define_attribute ( adios_group , "component", "desc", adios_string, &
        "Record component ", "" , adios_err)
  call adios_define_attribute ( adios_group, "gmt_year", "desc", adios_string, &
        "GMT year corresponding to reference (zero) time in file. ", "" , adios_err)
  call adios_define_attribute ( adios_group, "gmt_day", "desc", adios_string, &
        "GMT julian day corresponding to reference (zero) time in file. ", "" , adios_err)
  call adios_define_attribute ( adios_group, "gmt_hour", "desc", adios_string, &
        "GMT hour corresponding to reference (zero) time in file. ", "" , adios_err)
  call adios_define_attribute ( adios_group, "gmt_min", "desc", adios_string, &
        "GMT minute corresponding to reference (zero) time in file. ", "" , adios_err)
  call adios_define_attribute ( adios_group, "gmt_sec", "desc", adios_string, &
        "GMT second corresponding to reference (zero) time in file. ", "" , adios_err)
  call adios_define_attribute ( adios_group, "gmt_msec", "desc", adios_string, &
        "GMT millisecond corresponding to reference (zero) time in file. ", "" , adios_err)
  call adios_define_attribute ( adios_group , "receiver_lat", "desc", adios_string, &
        "Receiver latitude (degrees, north positive)  ", "" , adios_err )
  call adios_define_attribute ( adios_group , "receiver_lo", "desc", adios_string, &
        "Receiver longitude (degrees, east positive) ", "" , adios_err )
  call adios_define_attribute ( adios_group , "receiver_dpt", "desc", adios_string, &
        "Receiver depth below surface (meters) ", "" , adios_err )
  call adios_define_attribute ( adios_group , "receiver_el", "desc", adios_string, &
        "Receiver elevation (meters) ", "" , adios_err )
  call adios_define_attribute ( adios_group , "begin_value", "desc", adios_string, &
        "Beginning value of time array ", "" , adios_err )
  call adios_define_attribute ( adios_group , "end_value", "desc", adios_string, &
        "End value of time array ", "" , adios_err )
  call adios_define_attribute ( adios_group , "cmp_azimuth", "desc", adios_string, &
        "Component azimuth (degrees clockwise from north) ", "", adios_err )
  call adios_define_attribute ( adios_group , "cmp_incident_ang", "desc", adios_string,&
        "Component incident angle (degrees from vertical) ", "", adios_err )
  call adios_define_attribute ( adios_group , "sample_rate", "desc", adios_string, &
        "Sampling rate (s) ", "" , adios_err )
  call adios_define_attribute ( adios_group , "scale_factor", "desc", adios_string, &
        "Scale factor to convert the unit of synthetics from meters to nanometer ", &
        "" , adios_err )
  call adios_define_attribute ( adios_group , "ev_to_sta_AZ", "desc", adios_string, &
        "Event to station azimuth (degrees) ", "" , adios_err )
  call adios_define_attribute ( adios_group , "sta_to_ev_AZ", "desc", adios_string, &
        "Station to event azimuth (backazimuth, degrees) ", "", adios_err )
  call adios_define_attribute ( adios_group , "great_circle_dist", "desc", adios_string, &
        "Great circle distance between event and station (degrees) ", "", adios_err )
  call adios_define_attribute ( adios_group , "receiver_name", "desc", adios_string, &
        "Receiver name ", "" , adios_err )
  call adios_define_attribute( adios_group , "network", "desc", adios_string, &
        "Receiver network name ", "" , adios_err )
  call adios_define_attribute( adios_group , "receiver_id", "desc", adios_string, &
        "Receiver number ", "" , adios_err )
  call adios_define_attribute ( adios_group , "component", "desc", adios_string,&
        "Receiver component name ", "" , adios_err )

end subroutine define_asdf_data

!> Writes sac data to an asdf data file
!! \param file_name The file will be saved as file_name.
!! \param comm Size of the group associated with the MPI communicator

subroutine write_asdf_file(asdf_fn, my_asdf, adios_group, rank, nproc, comm, ierr)

  use asdf_data
  use adios_write_mod

  character(len=*) :: asdf_fn 
  type(asdf_event) :: my_asdf
  integer :: rank, nproc, comm, ierr

  integer        :: adios_err
  integer(kind=8)         :: adios_groupsize, adios_totalsize, varid
  integer(kind=8)         :: adios_handle, adios_group

	!print *,"Write out file: ", trim(asdf_fn)
  !print *,"comm:", comm
  !adios write init
  !call adios_init_noxml (comm, adios_err)
	!print *,"Write out file: ", trim(asdf_fn)
  !call adios_allocate_buffer (600, adios_err)
	!print *,"Write out file: ", trim(asdf_fn)
  !call adios_declare_group (adios_group, "EVENTS", "iter", 1, adios_err)
	!print *,"Write out file: ", trim(asdf_fn)
  !call adios_select_method (adios_group, "MPI", "", "", adios_err)

  !calculate size
  adios_groupsize = 0
	print *,"Write out file: ", trim(asdf_fn)
  print *, "Define adios data structure..."
  call define_asdf_data (adios_group, adios_groupsize, my_asdf,&
						rank, nproc, comm, ierr)
  !print *, "define finished!"
  call adios_open (adios_handle, "EVENTS", asdf_fn, "w", comm, adios_err)
  call adios_group_size (adios_handle, adios_groupsize, adios_totalsize, adios_err)

  !call the write sub
  call write_asdf_file_sub (my_asdf, adios_handle, adios_group,&
						adios_groupsize, rank, nproc, comm, ierr)

  !adios close
  call adios_close(adios_handle, adios_err)
  !print *, "adios_err", adios_err
  if(adios_err.eq.0) then
    print *, "Finish writing file:", trim(asdf_fn)
  endif
  !call adios_finalize (rank, adios_err)

end subroutine write_asdf_file


subroutine write_asdf_file_sub (my_asdf, adios_handle, my_adios_group, adios_groupsize, rank, nproc, comm, ierr)

  use adios_write_mod
  use asdf_data
	use adios_helpers_writers_mod
  !use seismo_variables

  implicit none
  integer                       :: adios_err, i
  integer(kind=8),intent(in)    :: my_adios_group, adios_groupsize
  integer(kind=8),intent(in)    :: adios_handle
  integer,intent(in)            :: rank, nproc, comm, ierr
	integer :: nrecords_total, offset
  integer :: receiver_name_len, network_len, component_len, receiver_id_len
	integer :: rn_len_total, nw_len_total, rid_len_total, comp_len_total
	integer :: rn_offset, nw_offset, rid_offset, comp_offset
  character(len=32)              :: loc_string

	character(len=:), allocatable :: receiver_name, network, component, receiver_id
	character(len=:), allocatable :: receiver_name_total, network_total, &
                                  component_total, receiver_id_total

  type(asdf_event), intent(inout) :: my_asdf

  !gather array offset info
	call gather_offset_info(my_asdf%nrecords,nrecords_total,offset,&
					rank, nproc, comm, ierr)

  !ensemble the string for receiver_name, network, componen and receiver_id
  allocate(character(len=6*my_asdf%nrecords) :: receiver_name)
  allocate(character(len=6*my_asdf%nrecords) :: network)
  allocate(character(len=6*my_asdf%nrecords) :: component)
  allocate(character(len=6*my_asdf%nrecords) :: receiver_id)
  receiver_name=''
  network=''
  component=''
  receiver_id=''
	do i=1, my_asdf%nrecords
		receiver_name=trim(receiver_name)//trim(my_asdf%receiver_name_array(i))//'.'
		network=trim(network)//trim(my_asdf%network_array(i))//'.'
		component=trim(component)//trim(my_asdf%component_array(i))//'.'
		receiver_id=trim(receiver_id)//trim(my_asdf%receiver_id_array(i))//'.'
	enddo
	receiver_name_len = len_trim(receiver_name)
	network_len = len_trim(network)
	component_len = len_trim(component)
	receiver_id_len = len_trim(receiver_id)

	call gather_string_offset_info(receiver_name_len, rn_len_total, rn_offset, &
					receiver_name, receiver_name_total,&
					rank, nproc, comm, ierr)
	call gather_string_offset_info(network_len, nw_len_total, nw_offset, &
					network, network_total,&
					rank, nproc, comm, ierr)
	call gather_string_offset_info(receiver_id_len, rid_len_total, rid_offset, &
					receiver_id, receiver_id_total,&
					rank, nproc, comm, ierr)
	call gather_string_offset_info(component_len, comp_len_total, comp_offset, &
					component, component_total,&
					rank, nproc, comm, ierr)

  !===========================
  !write out the string info
  print *,"write string"
	if(rank.eq.0)then
		!print *,"string_gathered:", trim(receiver_name)
  	call adios_write(adios_handle, "receiver_name", trim(receiver_name), adios_err)
  	call adios_write(adios_handle, "network", trim(network), adios_err)
  	call adios_write(adios_handle, "component", trim(component), adios_err)
  	call adios_write(adios_handle, "receiver_id", trim(receiver_id), adios_err) 
	endif

  !===========================
  print *,"Write seismic record"
  do i = 1, my_asdf%nrecords
  	write( loc_string, '(I10)' ) i+offset
		loc_string=trim(my_asdf%receiver_name_array(i))//"."//&
						trim(my_asdf%network_array(i))//"."//&
						trim(my_asdf%component_array(i))//"."//&
						trim(my_asdf%receiver_id_array(i))
    call write_adios_global_real_1d_array(adios_handle, rank, nproc, &
	 				my_asdf%npoints(i), my_asdf%npoints(i), 0, &
					loc_string, my_asdf%records(i)%record)
  enddo

  !===========================
  !scalar
  print *,"write scalar"
	if(rank.eq.0)then
  	call adios_write(adios_handle, "nrecords", nrecords_total, adios_err)

    call adios_write(adios_handle, "receiver_name_len", rn_len_total, adios_err)
  	call adios_write(adios_handle, "network_len", nw_len_total, adios_err)
  	call adios_write(adios_handle, "component_len", comp_len_total, adios_err)
  	call adios_write(adios_handle, "receiver_id_len", rid_len_total, adios_err)

  	call adios_write(adios_handle, "nreceivers", my_asdf%nreceivers, adios_err)

  	call adios_write(adios_handle, "min_period", my_asdf%min_period, adios_err) 
  	call adios_write(adios_handle, "max_period", my_asdf%max_period, adios_err) 
	
		call adios_write(adios_handle, "event", my_asdf%event, adios_err)
		!print *, "tag:",trim(my_asdf%event)
	endif

  !===========================
  !write out the array using the offset info
  print *,"write array"
  call write_adios_global_integer_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "npoints", my_asdf%npoints)

  call write_adios_global_integer_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "gmt_year", my_asdf%gmt_year)
  call write_adios_global_integer_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "gmt_day", my_asdf%gmt_day)
  call write_adios_global_integer_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "gmt_hour", my_asdf%gmt_hour)
  call write_adios_global_integer_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "gmt_min", my_asdf%gmt_min)
  call write_adios_global_integer_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "gmt_sec", my_asdf%gmt_sec)
  call write_adios_global_integer_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "gmt_msec", my_asdf%gmt_msec)

  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords, &
        nrecords_total, offset, "event_lat", my_asdf%event_lat)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords, &
        nrecords_total, offset, "event_lo", my_asdf%event_lo)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords, &
        nrecords_total, offset, "event_dpt", my_asdf%event_dpt)

  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords, &
        nrecords_total, offset, "receiver_lat", my_asdf%receiver_lat)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,& 
				nrecords_total, offset, "receiver_lo", my_asdf%receiver_lo)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
				nrecords_total, offset, "receiver_el", my_asdf%receiver_el)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
				nrecords_total, offset, "receiver_dpt", my_asdf%receiver_dpt)

  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,& 
				nrecords_total, offset, "begin_value", my_asdf%begin_value)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,& 
				nrecords_total, offset, "end_value", my_asdf%end_value)

  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,& 
				nrecords_total, offset, "cmp_azimuth", my_asdf%cmp_azimuth)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,& 
				nrecords_total, offset, "cmp_incident_ang", my_asdf%cmp_incident_ang)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,& 
				nrecords_total, offset, "sample_rate", my_asdf%sample_rate)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
				nrecords_total, offset, "scale_factor", my_asdf%scale_factor)

  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
				nrecords_total, offset, "ev_to_sta_AZ", my_asdf%ev_to_sta_AZ)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,& 
				nrecords_total, offset, "sta_to_ev_AZ", my_asdf%sta_to_ev_AZ)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,& 
				nrecords_total, offset, "great_circle_arc", my_asdf%great_circle_arc)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,& 
				nrecords_total, offset, "dist", my_asdf%dist)

  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,& 
				nrecords_total, offset, "P_pick", my_asdf%P_pick)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,& 
				nrecords_total, offset, "S_pick", my_asdf%S_pick)

  !end of write_asdf_file_sub

end subroutine write_asdf_file_sub

subroutine gather_offset_info(local_dim, global_dim, offset,&
						rank, nproc, comm, ierr)

	use mpi
	implicit none

	integer :: local_dim, global_dim, offset
	integer :: rank, nproc, comm, ierr

	integer, allocatable :: local_dim_all_proc(:)
	integer, allocatable :: offset_all_proc(:)
	integer :: i

	!if(rank.eq.0)then
		allocate(local_dim_all_proc(nproc))
		allocate(offset_all_proc(nproc))
	!endif
	
	call MPI_Barrier(comm, ierr)

	call MPI_Gather(local_dim, 1, MPI_INTEGER, local_dim_all_proc, 1, &
					MPI_INTEGER, 0, comm, ierr)

	if(rank.eq.0)then
		offset_all_proc(1)=0
		do i=2, nproc
			offset_all_proc(i)=sum(local_dim_all_proc(1:(i-1)))
		enddo
		global_dim=sum(local_dim_all_proc(1:nproc))
		!print *, "offset_all_proc:", offset_all_proc(:)
	endif

	call MPI_Scatter(offset_all_proc, 1, MPI_INTEGER, offset, &
					1, MPI_INTEGER, 0, comm, ierr)
	call MPI_Bcast(global_dim, 1, MPI_INTEGER, 0, comm, ierr)

	!print *,"rank, local dim, global_dim,offset:", rank, local_dim, &
	!						global_dim, offset

end subroutine gather_offset_info


subroutine gather_string_offset_info(local_dim, global_dim, offset,&
					 string_piece, string_total,&
						rank, nproc, comm, ierr)

	use mpi
	implicit none

	integer :: local_dim, global_dim, offset
	character(len=*) :: string_piece
	character(len=:), allocatable :: string_total
	character(len=10000) :: buffer_string
	!character(len=:), allocatable :: buffer_string
	integer :: rank, nproc, comm, ierr

	integer, allocatable :: local_dim_all_proc(:)
	integer, allocatable :: offset_all_proc(:)
	integer :: i, tag, mpi_status(MPI_STATUS_SIZE)

	!if(rank.eq.0)then
		allocate(local_dim_all_proc(nproc))
		allocate(offset_all_proc(nproc))
	!endif
	
	call MPI_Barrier(comm, ierr)

	call MPI_Gather(local_dim, 1, MPI_INTEGER, local_dim_all_proc, 1, &
					MPI_INTEGER, 0, comm, ierr)

	!allocate(character(len=10000) :: buffer_string )

	if(rank.eq.0)then
		offset_all_proc(1)=0
		do i=2, nproc
			offset_all_proc(i)=sum(local_dim_all_proc(1:(i-1)))
		enddo
		global_dim=sum(local_dim_all_proc(1:nproc))
		!print *, "offset_all_proc:", offset_all_proc(:)
		allocate(character(len=global_dim) :: string_total)
		!allocate(character(len=global_dim) :: buffer_string)
		string_total=""
		buffer_string=""
		string_total=trim(string_total)//trim(string_piece(1:local_dim))
	endif
	
	!print *,"TAG1"
	!if(rank.eq.0) then
!		print *,"global_dim",global_dim
!	endif

	if(rank.eq.0)then
		do i=1,nproc-1
			!print *, "buffer_before:",trim(buffer_string)
			!print *, "local_dim_all_proc:",local_dim_all_proc(i+1)
			call MPI_Recv(buffer_string, local_dim_all_proc(i+1), MPI_CHARACTER,&
							i, 1, comm, mpi_status, ierr)
			!print *,"buffer_string:", trim(buffer_string)
			string_total=trim(string_total)//buffer_string(1:local_dim_all_proc(i+1))
		enddo
	else
		!print *, "local_dim:", local_dim
		!print *,"string_piece:", trim(string_piece)
		call MPI_Send(string_piece, local_dim, MPI_CHARACTER,&
							0, 1, comm, ierr)
	endif
	!print *,"TAG", rank

	call MPI_Scatter(offset_all_proc, 1, MPI_INTEGER, offset, &
					1, MPI_INTEGER, 0, comm, ierr)
	call MPI_Bcast(global_dim, 1, MPI_INTEGER, 0, comm, ierr)

	!print *,"rank, local dim, global_dim,offset:", rank, local_dim, &
!							global_dim, offset

end subroutine gather_string_offset_info


end module asdf_write_subs
