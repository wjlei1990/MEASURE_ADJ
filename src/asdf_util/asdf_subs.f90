module asdf_subs

!> The module asdf_subs contains Parallel ASDF I/O API:
!! 1)define_asdf_data
!! 2)read_asdf_file
!! 3)write_asdf_file

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
  !use seismo_variables
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
	!print *,"herere", nrecords

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
  !print *, "TAG"
  !time info
  !call define_adios_scalar (adios_group, my_group_size, "", "gmt_year", &
  !                      dum_int)
  !call define_adios_scalar (adios_group, my_group_size, "", "gmt_month", &
  !                      dum_int)
  !call define_adios_scalar (adios_group, my_group_size, "", "gmt_day", &
  !                      dum_int)
  !call define_adios_scalar (adios_group, my_group_size, "", "gmt_hour", &
  !                      dum_int)
  !call define_adios_scalar (adios_group, my_group_size, "", "gmt_min", &
  !                      dum_int)
  !call define_adios_scalar (adios_group, my_group_size, "", "gmt_sec", &
  !                      dum_int)
  !call define_adios_scalar (adios_group, my_group_size, "", "gmt_msec", &
  !                      dum_int)
  !location info
  !call define_adios_scalar (adios_group, my_group_size, "", "event_lat", &
  !                      dum_real)
  !call define_adios_scalar (adios_group, my_group_size, "", "event_lo", &
  !                      dum_real)
  !call define_adios_scalar (adios_group, my_group_size, "", "event_dpt", &
  !                      dum_real)

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
    !print *,"=="
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
  !stop

  !DISPLACEMENT
  do i = 1, nrecords
	  !print *,"begin define:",i, my_asdf%npoints(i)
    write (i_string, '(I10)' ) i+offset
		record=trim(my_asdf%receiver_name_array(i))//"."//&
						trim(my_asdf%network_array(i))//"."//&
						trim(my_asdf%component_array(i))//"."//&
						trim(my_asdf%receiver_id_array(i))
		!print *, trim(record)
    !write( npts_string, '(I6)' ) npts(i)
    !call define_adios_global_real_1d_array (adios_group, my_group_size,&
    !      my_asdf%npoints(i), "", trim('DISPLACEMENT#'//trim(adjustl(i_string))),&
    !      real_array)
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

subroutine write_asdf_file(asdf_fn, my_asdf, rank, nproc, comm, ierr)

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
  call adios_init_noxml (comm, adios_err)
	!print *,"Write out file: ", trim(asdf_fn)
  call adios_allocate_buffer (200, adios_err)
	!print *,"Write out file: ", trim(asdf_fn)
  call adios_declare_group (adios_group, "EVENTS", "iter", 1, adios_err)
	!print *,"Write out file: ", trim(asdf_fn)
  call adios_select_method (adios_group, "MPI", "", "", adios_err)

  !calculate size
  adios_groupsize = 0
	print *,"Write out file: ", trim(asdf_fn)
  call define_asdf_data (adios_group, adios_groupsize, my_asdf,&
						rank, nproc, comm, ierr)
  print *, "define finished!"
  call adios_open (adios_handle, "EVENTS", asdf_fn, "w", comm, adios_err)
  call adios_group_size (adios_handle, adios_groupsize, adios_totalsize, adios_err)

  !call the write sub
  call write_asdf_file_sub (my_asdf, adios_handle, adios_group,&
						adios_groupsize, rank, nproc, comm, ierr)

  !adios close
  call adios_close(adios_handle, adios_err)
  call adios_finalize (rank, adios_err)

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
	integer :: rn_len_total, nw_len_total, rid_len_total, comp_len_total
	integer :: rn_offset, nw_offset, rid_offset, comp_offset
  character(len=32)              :: loc_string

	character(len=:), allocatable :: receiver_name, network, component, receiver_id

  type(asdf_event), intent(inout) :: my_asdf

	call gather_offset_info(my_asdf%nrecords,nrecords_total,offset,&
					rank, nproc, comm, ierr)
	
	call gather_string_offset_info(my_asdf%receiver_name_len, rn_len_total, rn_offset, &
					my_asdf%receiver_name, receiver_name,&
					rank, nproc, comm, ierr)
	call gather_string_offset_info(my_asdf%network_len, nw_len_total, nw_offset, &
					my_asdf%network, network,&
					rank, nproc, comm, ierr)
	call gather_string_offset_info(my_asdf%receiver_id_len, rid_len_total, rid_offset, &
					my_asdf%receiver_id, receiver_id,&
					rank, nproc, comm, ierr)
	call gather_string_offset_info(my_asdf%component_len, comp_len_total, comp_offset, &
					my_asdf%component, component,&
					rank, nproc, comm, ierr)

  !string
  print *,"write string"
	if(rank.eq.0)then
		!print *,"string_gathered:", trim(receiver_name)
  	call adios_write(adios_handle, "receiver_name", trim(receiver_name), adios_err)
  	call adios_write(adios_handle, "network", trim(network), adios_err)
  	call adios_write(adios_handle, "component", trim(component), adios_err)
  	call adios_write(adios_handle, "receiver_id", trim(receiver_id), adios_err) 
	endif

  print *,"Write seismic record"
  !call adios_write(adios_handle, "event", my_asdf%event, adios_err)
  do i = 1, my_asdf%nrecords
  	write( loc_string, '(I10)' ) i+offset
	 	!loc_string=trim("DISPLACEMENT#")//trim(adjustl(loc_string))
		loc_string=trim(my_asdf%receiver_name_array(i))//"."//&
						trim(my_asdf%network_array(i))//"."//&
						trim(my_asdf%component_array(i))//"."//&
						trim(my_asdf%receiver_id_array(i))
	  !print *, trim(loc_string)
    call write_adios_global_real_1d_array(adios_handle, rank, nproc, &
	 				my_asdf%npoints(i), my_asdf%npoints(i), 0, &
					loc_string, my_asdf%records(i)%record)
  enddo
  !scalar
  print *,"write scalar"
	if(rank.eq.0)then
  	call adios_write(adios_handle, "nrecords", nrecords_total, adios_err)

    call adios_write(adios_handle, "receiver_name_len", rn_len_total, adios_err)
  	call adios_write(adios_handle, "network_len", nw_len_total, adios_err)
  	call adios_write(adios_handle, "component_len", comp_len_total, adios_err)
  	call adios_write(adios_handle, "receiver_id_len", rid_len_total, adios_err)

  	call adios_write(adios_handle, "nreceivers", my_asdf%nreceivers, adios_err)

  	!call adios_write(adios_handle, "gmt_year", my_asdf%gmt_year, adios_err)
		!call adios_write(adios_handle, "gmt_month", my_asdf%gmt_month, adios_err)
 	 	!call adios_write(adios_handle, "gmt_day", my_asdf%gmt_day, adios_err)
  	!call adios_write(adios_handle, "gmt_hour", my_asdf%gmt_hour, adios_err)
  	!call adios_write(adios_handle, "gmt_min", my_asdf%gmt_min, adios_err)
  	!call adios_write(adios_handle, "gmt_sec", my_asdf%gmt_sec, adios_err)
  	!call adios_write(adios_handle, "gmt_msec", my_asdf%gmt_msec, adios_err) 

  	!call adios_write(adios_handle, "event_lat", my_asdf%event_lat, adios_err) 
  	!call adios_write(adios_handle, "event_lo", my_asdf%event_lo, adios_err) 
  	!call adios_write(adios_handle, "event_dpt", my_asdf%event_dpt, adios_err) 
		!string needs special attention in adios: event

  	call adios_write(adios_handle, "min_period", my_asdf%min_period, adios_err) 
  	call adios_write(adios_handle, "max_period", my_asdf%max_period, adios_err) 
	
		!call write_adios_global_string_1d_array(adios_handle, rank, nproc, &
		!		 	13, 13, 0, "event", my_asdf%event)
		call adios_write(adios_handle, "event", my_asdf%event, adios_err)
		!print *, "tag:",trim(my_asdf%event)
	endif


  !array
	!call write_adios_global_string_1d_array(adios_handle, rank, nproc,&
	!			my_asdf%receiver_name_len, rn_len_total, rn_offset,&
	!			"receiver_name", my_asdf%receiver_name)
	!call write_adios_global_string_1d_array(adios_handle, rank, nproc,&
	!			my_asdf%network_len, nw_len_total, nw_offset,&
	!			"network", my_asdf%network)
	!call write_adios_global_string_1d_array(adios_handle, rank, nproc,&
	!			my_asdf%receiver_id_len, rid_len_total, rid_offset,&
	!			"receiver_id", my_asdf%receiver_id)
	!call write_adios_global_string_1d_array(adios_handle, rank, nproc,&
!				my_asdf%component_len, comp_len_total, comp_offset,&
!				"component", my_asdf%component)

  print *,"write_array"
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
			print *, "buffer_before:",trim(buffer_string)
			print *, "local_dim_all_proc:",local_dim_all_proc(i+1)
			call MPI_Recv(buffer_string, local_dim_all_proc(i+1), MPI_CHARACTER,&
							i, 1, comm, mpi_status, ierr)
			print *,"buffer_string:", trim(buffer_string)
			string_total=trim(string_total)//buffer_string(1:local_dim_all_proc(i+1))
		enddo
	else
		print *, "local_dim:", local_dim
		print *,"string_piece:", trim(string_piece)
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


!! \param file_name The name of the asdf file to read
!! \param my_asdf The asdf object that will be read
!! \param comm Size of the group associated with the MPI communicator

subroutine read_asdf_file (file_name, my_asdf, rank, nproc, comm, ierr)

  use adios_read_mod
  use asdf_data
  implicit none

  character(len=100),intent(in) :: file_name
  type(asdf_event), intent(inout) :: my_asdf
  integer,intent(in) :: rank, nproc, comm
	integer :: ierr
  !integer, dimension(:), allocatable :: loc_begin, loc_end

  integer                 :: i, j, scalarnum, local_index
  integer                 :: vcnt, acnt, tfirst, tlast
  integer                 :: vrank, vtype, vsteps, adios_err

  integer(kind=8)         :: fh, sel=0
  integer(kind=8),dimension(1)   :: start, count
  integer(kind=8),dimension(10)  :: readsize, dims
  integer :: dim_array
  integer :: nrecords_total, nrecords_local, nrecords
  integer :: loc_begin, loc_end

  character(len=128),dimension(:),allocatable :: vnamelist

  character(len=20), allocatable :: receiver_name_array_temp(:)
  character(len=20), allocatable :: network_array_temp(:)
  character(len=20), allocatable :: component_array_temp(:)
  character(len=20), allocatable :: receiver_id_array_temp(:)

  character(len=30) :: loc_string
  integer :: receiver_name_len, network_len, component_len, receiver_id_len
  character(len=:), allocatable :: receiver_name, network
  character(len=:), allocatable :: component, receiver_id 


  real,dimension(20000)::temp_record

  !print *,"checkpoint 6"

  !>Initialization,Get Varname and Varnumber
  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, "verbose=2", ierr)
  !print *,"checkpoint 6"
  !print *,"filename",trim(adjustl(file_name))
  call adios_read_open_file (fh, trim(adjustl(file_name)), 0, comm, ierr)
  !print *,"checkpoint 6"
  call adios_inq_file (fh, vcnt, acnt, tfirst, tlast, ierr)

  !print *,"checkpoint 5"

  allocate (vnamelist(vcnt))
  call adios_inq_varnames (fh, vnamelist, ierr)
	if(rank==0) then
  	write (*, '("Number of variables : ", i0)') vcnt
  	do i = 1, vcnt
    	write (*, "(i5, a, a)") i,") ", trim(vnamelist(i))
  	enddo
	endif

  !print *,"checkpoint 2"

  !if(rank==0)then
  !>get the number of records and npts
    call adios_get_scalar (fh,"/nrecords",nrecords_total,ierr)
		if(rank==0)then
    	print *,"nrecords_total:", nrecords_total
		endif

    allocate(receiver_name_array_temp(nrecords_total))
    allocate(network_array_temp(nrecords_total))
    allocate(component_array_temp(nrecords_total))
    allocate(receiver_id_array_temp(nrecords_total))
    !print *,"here"
    !call split_job_mpi_simple(nrecords_total, nrecords_local, loc_begin, loc_end, &
    !                    rank,nproc)
    
    !introduce a complex split_job_mpi subroutine
    !which requires the three components stay together
    call adios_get_scalar (fh,"/receiver_name_len",receiver_name_len,ierr)
    call adios_get_scalar (fh,"/network_len",network_len,ierr)
    call adios_get_scalar (fh,"/component_len",component_len,ierr)
    call adios_get_scalar (fh,"/receiver_id_len",receiver_id_len,ierr)
    allocate(character(len=6*nrecords_total) :: receiver_name)
    allocate(character(len=6*nrecords_total) :: network)
    allocate(character(len=6*nrecords_total) :: component)
    allocate(character(len=6*nrecords_total) :: receiver_id)
    call adios_get_scalar (fh, "/receiver_name", receiver_name, ierr)
    call adios_get_scalar (fh, "/network", network, ierr)
    call adios_get_scalar (fh, "/component", component, ierr)
    call adios_get_scalar (fh, "/receiver_id", receiver_id, ierr)
    !print *, "receiver_name: ", trim(receiver_name)
    !stop
    call split_job_mpi_complex(nrecords_total,receiver_name, receiver_name_len, &
              nrecords_local, loc_begin, loc_end, rank, nproc)
    print *,"rank, nproc, loc_begin,loc_end,nrecords_local:", &
									rank, nproc, loc_begin, loc_end, nrecords_local
    nrecords=nrecords_local
		!stop
  !allocate (my_asdf%npoints(my_asdf%nrecords))
  !call adios_get_scalar (fh, "/npoints",my_asdf%npoints,ierr)
  !call adios_get_scalar (fh,"/nrecords",my_asdf%nrecords,ierr)
  !print *, "Number of Records:",nrecords
  !scalarnum = 2 + nrecords 

  !call adios_get_scalar (fh, vnamelist(scalarnum+1), my_asdf%nrecords, ierr)
  !print *,"Number of Records:", my_asdf%nrecords
  !>allocate variables
  call init_asdf_data(my_asdf, nrecords)

	!print *, "receiver_name:", trim(receiver_name)
	!print *, "receiver_name_len", receiver_name_len
  !split the long string
  !print *, "my_asdf%receiver_name", trim(my_asdf%receiver_name)
  call split_string(receiver_name,receiver_name_len, &
                            receiver_name_array_temp,dim_array,'.')
  call split_string(network,network_len, &
                            network_array_temp,dim_array,'.')
  call split_string(component,component_len,&
                            component_array_temp,dim_array,'.')
  call split_string(receiver_id,receiver_id_len,&
                            receiver_id_array_temp,dim_array,'.')

	!print *,trim(receiver_name_array_temp(1))

  my_asdf%receiver_name_array(1:nrecords_local) = &
                          receiver_name_array_temp(loc_begin:loc_end)
  my_asdf%network_array(1:nrecords_local) = &
                          network_array_temp(loc_begin:loc_end)
  my_asdf%component_array(1:nrecords_local) = &
                          component_array_temp(loc_begin:loc_end)
  my_asdf%receiver_id_array(1:nrecords_local) = &
                          receiver_id_array_temp(loc_begin:loc_end)

!-----------------------------------------------
!>read all the records 
  do i=loc_begin, loc_end
		local_index=i-loc_begin+1
		!print *, local_index
		!print *,trim(my_asdf%receiver_name_array(local_index))
    !print *, "i,vnamelist:", i, trim(vnamelist(i+1))
    !write(loc_string,'(I6)') i 
    !loc_string="/DISPLACEMENT#"//trim(adjustl(loc_string))
    !print *, "loc_string:", trim(loc_string)
		loc_string=trim(my_asdf%receiver_name_array(local_index))//"."//&
						trim(my_asdf%network_array(local_index))//"."//&
						trim(my_asdf%component_array(local_index))//"."//&
						trim(my_asdf%receiver_id_array(local_index))
    !call adios_inq_var(fh,loc_string,vtype,vsteps,vrank,dims,ierr)
		call adios_get_scalar(fh, trim(loc_string)//"/global_dim",dims(1), ierr)
		!print *, trim(loc_string), dims(1)
		!stop
    allocate (my_asdf%records(local_index)%record(dims(1)))
    start(1) = 0
    count(1) = dims(1)
    call adios_selection_boundingbox (sel, 1 , start , count )
    !call adios_schedule_read (fh, sel, vnamelist(i+1), 0, 1, my_asdf%record(i,:), ierr)
    !call adios_schedule_read (fh, sel, vnamelist(i+1), 0, 1, my_asdf%records(i)%record, ierr)
    call adios_schedule_read (fh, sel, trim(loc_string)//"/array", 0, 1, &
						my_asdf%records(local_index)%record, ierr)
    !my_asdf%record(i,:)=temp_record(:)
  enddo

  !print *," checkpoint 1"

 ! stop 

  !print *,"reading records finished!"

!  j = 2
!  do i = 1, my_asdf%nreceivers
!    call adios_inq_var (fh, vnamelist(j), vtype, vsteps, vrank, dims, ierr)
!    !allocate (my_asdf%records(i)%record1(dims(1)))
!    start(1) = 0
!    count(1) = dims(1)
!    call adios_selection_boundingbox (sel, 1 , start , count )
!
!    call adios_inq_var (fh, vnamelist(j+1), vtype, vsteps, vrank, dims, ierr)
!    allocate (my_asdf%records(i)%record2(dims(1)))
!
!    call adios_inq_var (fh, vnamelist(j+2), vtype, vsteps, vrank, dims, ierr)
!    allocate (my_asdf%records(i)%record3(dims(1)))
!   
!    call adios_schedule_read (fh, sel, vnamelist(j), 0, 1, my_asdf%records(i)%record1, ierr)
!    call adios_schedule_read (fh, sel, vnamelist(j+1), 0, 1, my_asdf%records(i)%record2, ierr)
!    call adios_schedule_read (fh, sel, vnamelist(j+2), 0, 1, my_asdf%records(i)%record3, ierr)
!    j = j + 3
!  enddo

!--------------------------------------
!>read in earthquake information
  !call adios_get_scalar (fh, vnamelist(1), my_asdf%event, ierr)  
  !call adios_get_scalar (fh, vnamelist(scalarnum),   my_asdf%nreceivers, ierr)
  call adios_get_scalar (fh, "/event", my_asdf%event, ierr)  
  !print *,"/event",my_asdf%event
  !stop
  call adios_get_scalar (fh, "/nreceivers", my_asdf%nreceivers, ierr)
  !print *,"nreceiver:",my_asdf%nreceivers

  !call adios_get_scalar (fh,"/receiver_name_len",my_asdf%receiver_name_len,ierr)
  !call adios_get_scalar (fh,"/network_len",my_asdf%network_len,ierr)
  !call adios_get_scalar (fh,"/component_len",my_asdf%component_len,ierr)
  !call adios_get_scalar (fh,"/receiver_id_len",my_asdf%receiver_id_len,ierr)
  !print *,"checkpoint 3"
  !allocate(character(my_asdf%receiver_name_len) :: my_asdf%receiver_name)
  !allocate(character(my_asdf%network_len) :: my_asdf%network)
  !allocate(character(my_asdf%component_len) :: my_asdf%component)
  !allocate(character(my_asdf%receiver_id_len) :: my_asdf%receiver_id)

  !call adios_get_scalar (fh, "/receiver_name", my_asdf%receiver_name, ierr)
  !call adios_get_scalar (fh, "/network", my_asdf%network, ierr)
  !call adios_get_scalar (fh, "/component", my_asdf%component, ierr)
  !call adios_get_scalar (fh, "/receiver_id", my_asdf%receiver_id, ierr)
  !print *,"event:", trim(my_asdf%event)
  !print *, "recevier_name:", trim(my_asdf%receiver_name)

  !print *,"len:", my_asdf%receiver_name_len, my_asdf%
  !stop

  print *, "my_asdf%event:",my_asdf%event

  !stop

  !call adios_get_scalar (fh, "/gmt_year", my_asdf%gmt_year, ierr)
	!call adios_get_scalar (fh, "/gmt_month", my_asdf%gmt_month, ierr)
  !call adios_get_scalar (fh, "/gmt_day", my_asdf%gmt_day, ierr)
  !call adios_get_scalar (fh, "/gmt_hour", my_asdf%gmt_hour, ierr)
  !print *,"gmt_year,day,hour",my_asdf%gmt_year,my_asdf%gmt_day,my_asdf%gmt_hour
  !print *, "my_asdf%event:",my_asdf%event
  !call adios_get_scalar (fh, "/gmt_min", my_asdf%gmt_min, ierr)
  !call adios_get_scalar (fh, "/gmt_sec", my_asdf%gmt_sec, ierr)
  !call adios_get_scalar (fh, "/gmt_msec", my_asdf%gmt_msec, ierr)
  !print *, "my_asdf%event:",my_asdf%event
  !call adios_get_scalar (fh, "/event_lat", my_asdf%event_lat, ierr)
  !print *, "my_asdf%event:",my_asdf%event
  !call adios_get_scalar (fh, "/event_lo", my_asdf%event_lo, ierr)
  !print *, "my_asdf%event:",my_asdf%event
  !call adios_get_scalar (fh, "/event_dpt", my_asdf%event_dpt, ierr)

	call adios_get_scalar (fh, "/max_period", my_asdf%max_period, ierr)
	call adios_get_scalar (fh, "/min_period", my_asdf%min_period, ierr)
  !stop
  !print *, "my_asdf%event:",my_asdf%event
  !print *,"here, error?"


  start(1) = loc_begin-1 
  count(1) = my_asdf%nrecords
  call adios_selection_boundingbox (sel, 1 , start , count )
  call adios_schedule_read (fh, sel, "npoints/array", 0, 1, my_asdf%npoints, ierr)

  call adios_schedule_read (fh, sel, "gmt_year/array", 0, 1, my_asdf%gmt_year, ierr)
  call adios_schedule_read (fh, sel, "gmt_day/array", 0, 1, my_asdf%gmt_day, ierr)
  call adios_schedule_read (fh, sel, "gmt_hour/array", 0, 1, my_asdf%gmt_hour, ierr)
  call adios_schedule_read (fh, sel, "gmt_min/array", 0, 1, my_asdf%gmt_min, ierr)
  call adios_schedule_read (fh, sel, "gmt_sec/array", 0, 1, my_asdf%gmt_sec, ierr)
  call adios_schedule_read (fh, sel, "gmt_msec/array", 0, 1, my_asdf%gmt_msec, ierr)

  call adios_schedule_read (fh, sel, "event_lat/array", 0, 1, my_asdf%event_lat, ierr)
  call adios_schedule_read (fh, sel, "event_lo/array", 0, 1, my_asdf%event_lo, ierr)
  call adios_schedule_read (fh, sel, "event_dpt/array", 0, 1, my_asdf%event_dpt, ierr)

  call adios_schedule_read (fh, sel, "receiver_lat/array", 0, 1, my_asdf%receiver_lat, ierr)
  call adios_schedule_read (fh, sel, "receiver_lo/array", 0, 1, my_asdf%receiver_lo, ierr)
  call adios_schedule_read (fh, sel, "receiver_el/array", 0, 1, my_asdf%receiver_el, ierr)
  call adios_schedule_read (fh, sel, "receiver_dpt/array", 0, 1, my_asdf%receiver_dpt, ierr)
  !call adios_schedule_read (fh, sel, "/receiver_name", 0, 1, my_asdf%receiver_dpt, ierr)
  !call adios_schedule_read (fh, sel, "/network", 0, 1, my_asdf%network, ierr)
  !call adios_schedule_read (fh, sel, "/component", 0, 1, my_asdf%component, ierr)
  call adios_schedule_read (fh, sel, "begin_value/array", 0, 1, my_asdf%begin_value, ierr)
  call adios_schedule_read (fh, sel, "end_value/array", 0, 1, my_asdf%end_value, ierr)
  call adios_schedule_read (fh, sel, "cmp_azimuth/array", 0, 1, my_asdf%cmp_azimuth, ierr)
  call adios_schedule_read (fh, sel, "cmp_incident_ang/array", 0, 1, my_asdf%cmp_incident_ang, ierr)
  call adios_schedule_read (fh, sel, "sample_rate/array", 0, 1, my_asdf%sample_rate, ierr)
  call adios_schedule_read (fh, sel, "scale_factor/array", 0, 1, my_asdf%scale_factor, ierr)
  call adios_schedule_read (fh, sel, "ev_to_sta_AZ/array", 0, 1, my_asdf%ev_to_sta_AZ, ierr)
  call adios_schedule_read (fh, sel, "sta_to_ev_AZ/array", 0, 1, my_asdf%sta_to_ev_AZ, ierr)
  call adios_schedule_read (fh, sel, "great_circle_arc/array", 0, 1, my_asdf%great_circle_arc, ierr)
  call adios_schedule_read (fh, sel, "dist/array", 0, 1, my_asdf%dist, ierr)
  call adios_schedule_read (fh, sel, "P_pick/array", 0, 1, my_asdf%P_pick, ierr)
  call adios_schedule_read (fh, sel, "S_pick/array", 0, 1, my_asdf%S_pick, ierr)

  !print *, "my_asdf%event:",my_asdf%event
  !print *, "XXXX"

  call adios_perform_reads (fh, ierr)
  !print *, "XXXX"
  call adios_read_close(fh, ierr)
  !print *, "XXXX"
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, adios_err)

  !print *, "XXXX"

  !do i=1,my_asdf%nrecords
  !  print *,trim(vnamelist(i+1)),my_asdf%npoints(i)
  !enddo
 
  deallocate(vnamelist)

  !split the long string
  !print *, "my_asdf%receiver_name", trim(my_asdf%receiver_name)
  !call split_string(receiver_name,my_asdf%receiver_name_len, &
  !                          receiver_name_array_temp,dim_array,'.')
  !call split_string(network,my_asdf%network_len, &
  !                          network_array_temp,dim_array,'.')
  !call split_string(component,my_asdf%component_len,&
  !                          component_array_temp,dim_array,'.')
  !call split_string(receiver_id,my_asdf%receiver_id_len,&
  !                          receiver_id_array_temp,dim_array,'.')
!
  !my_asdf%receiver_name_array(1:nrecords_local) = &
  !                        receiver_name_array_temp(loc_begin:loc_end)
  !my_asdf%network_array(1:nrecords_local) = &
  !                        network_array_temp(loc_begin:loc_end)
  !my_asdf%component_array(1:nrecords_local) = &
  !                        component_array_temp(loc_begin:loc_end)
  !my_asdf%receiver_id_array(1:nrecords_local) = &
  !                        receiver_id_array_temp(loc_begin:loc_end)

	my_asdf%receiver_name = ''
	my_asdf%network = ''
	my_asdf%component = ''
	my_asdf%receiver_id = ''
	!Ensemble strings...
	do i=1, my_asdf%nrecords
		my_asdf%receiver_name=trim(my_asdf%receiver_name)//trim(my_asdf%receiver_name_array(i))//'.'
		my_asdf%network=trim(my_asdf%network)//trim(my_asdf%network_array(i))//'.'
		my_asdf%component=trim(my_asdf%component)//trim(my_asdf%component_array(i))//'.'
		my_asdf%receiver_id=trim(my_asdf%receiver_id)//trim(my_asdf%receiver_id_array(i))//'.'
	enddo
	my_asdf%receiver_name_len = len_trim(my_asdf%receiver_name)
	my_asdf%network_len = len_trim(my_asdf%network)
	my_asdf%component_len = len_trim(my_asdf%component)
	my_asdf%receiver_id_len = len_trim(my_asdf%receiver_id)

	print *,"receiver_name_len:",my_asdf%receiver_name_len
	print *,trim(my_asdf%receiver_name)
	print *,"network_len:",my_asdf%network_len
	print *,trim(my_asdf%network)
	print *,"receiver_id_len:",my_asdf%receiver_id_len
	print *,trim(my_asdf%receiver_id)
	print *,"component:", my_asdf%component_len
	print *,trim(my_asdf%component)
		
  !print *, "my_asdf%event:",my_asdf%event
	!call MPI_Barrier(comm, ierr)
	!do j=1, 10000000000*rank
	!enddo
	!print *, "rank, nrecords:", rank, my_asdf%nrecords
	!do i=1, my_asdf%nrecords
	!	print *, trim(my_asdf%receiver_name_array(i))
	!enddo

	!call MPI_Barrier(comm, ierr)
  

end subroutine read_asdf_file

subroutine split_string(string, string_len, string_array, dim_array, delimiter)

  implicit none

  character(len=*) :: string
  character(len=*) :: delimiter
  character(len=*) :: string_array(:)
  integer :: dim_array
  integer :: string_len
  character(len=:), allocatable :: string_temp

  integer :: i1,i2,i3,i

  !print *,"split begins"

  allocate(character(len=len(string))::string_temp)
  
  string_array(:)=''
  string_temp=''
  dim_array=0


  i1=1
  i3=string_len
  !print *, " string:" ,string
  !print *, "len:", i3
  !print *,"len_trim:",len_trim(string)

  do while (i1<=i3)
    i2=index(string(i1:i3),delimiter)
    dim_array=dim_array+1
    !print *,"dim_array:",dim_array
    string_array(dim_array)(1:(i2-1))=string(i1:(i1+i2-2))
    !print *, string_array(dim_array)
    !print *,"i1,i2,i3,dim,string_array():"
    !print *,i1,i2,i3,dim_array,trim(string_array(dim_array))
    i1=i1+i2
  enddo
  !print *, dim_array

end subroutine split_string

subroutine init_asdf_data(my_asdf,nrecords)
!receiver_name, network, component and receiver_id are not allocated

  use asdf_data
  type(asdf_event) :: my_asdf
  integer :: nrecords
  integer :: len_temp

  my_asdf%event = ""

  my_asdf%nrecords=nrecords

  !call adios_get_scalar (fh, vnamelist(scalarnum+1), my_asdf%nrecords, ierr)
  !print *,"Number of Records:", my_asdf%nrecords
  !>allocate variables
  allocate (my_asdf%npoints(my_asdf%nrecords))

  allocate (my_asdf%gmt_year(my_asdf%nrecords))
  allocate (my_asdf%gmt_hour(my_asdf%nrecords))
  allocate (my_asdf%gmt_day(my_asdf%nrecords))
  allocate (my_asdf%gmt_min(my_asdf%nrecords))
  allocate (my_asdf%gmt_sec(my_asdf%nrecords))
  allocate (my_asdf%gmt_msec(my_asdf%nrecords))

  allocate (my_asdf%event_lat(my_asdf%nrecords))
  allocate (my_asdf%event_lo(my_asdf%nrecords))
  allocate (my_asdf%event_dpt(my_asdf%nrecords))

  allocate (my_asdf%receiver_lat(my_asdf%nrecords))
  allocate (my_asdf%receiver_lo(my_asdf%nrecords))
  allocate (my_asdf%receiver_el(my_asdf%nrecords))
  allocate (my_asdf%receiver_dpt(my_asdf%nrecords))
  !allocate (my_asdf%receiver_name(my_asdf%nrecords))
  !allocate (my_asdf%network(my_asdf%nrecords))
  !allocate (my_asdf%component(my_asdf%nrecords))
  allocate (my_asdf%begin_value(my_asdf%nrecords))
  allocate (my_asdf%end_value(my_asdf%nrecords))
  allocate (my_asdf%cmp_azimuth(my_asdf%nrecords))
  allocate (my_asdf%cmp_incident_ang(my_asdf%nrecords))
  allocate (my_asdf%sample_rate(my_asdf%nrecords))
  allocate (my_asdf%scale_factor(my_asdf%nrecords))
  allocate (my_asdf%ev_to_sta_AZ(my_asdf%nrecords))
  allocate (my_asdf%sta_to_ev_AZ(my_asdf%nrecords))
  allocate (my_asdf%great_circle_arc(my_asdf%nrecords))
  allocate (my_asdf%dist(my_asdf%nrecords))
  allocate (my_asdf%P_pick(my_asdf%nrecords))
  allocate (my_asdf%S_pick(my_asdf%nrecords))
  !>the kernel part: allocate the record
  allocate (my_asdf%records(my_asdf%nrecords))

  allocate (my_asdf%receiver_name_array(my_asdf%nrecords))
  allocate (my_asdf%network_array(my_asdf%nrecords))
  allocate (my_asdf%component_array(my_asdf%nrecords))
  allocate (my_asdf%receiver_id_array(my_asdf%nrecords))

  len_temp=6*nrecords
  allocate (character(len=len_temp) :: my_asdf%receiver_name)
  allocate (character(len=len_temp) :: my_asdf%network)
  allocate (character(len=len_temp) :: my_asdf%component)
  allocate (character(len=len_temp) :: my_asdf%receiver_id)
!-----------------------------------------------
end subroutine init_asdf_data

subroutine split_job_mpi_simple(nrecords_total, nrecords_local, loc_begin, loc_end, &
              rank,nproc)

  !simplified version the mpi_job_split
  !future version should be based on three components data
  implicit none

  integer :: nrecords_total, nrecords_local, loc_begin, loc_end, rank, nproc

  integer :: nrecords_per_proc

  nrecords_per_proc=int(nrecords_total/nproc)

  loc_begin = nrecords_per_proc*rank+1
  loc_end = nrecords_per_proc*(rank+1)

  if(rank.eq.(nproc-1)) then !last node
    loc_end=nrecords_total
  endif

  nrecords_local=loc_end-loc_begin+1
  
  print *, "rank, loc_begin, loc_end:", rank, loc_begin, loc_end, &
              nrecords_per_proc, nproc, rank

  !stop
end subroutine split_job_mpi_simple


subroutine split_job_mpi_complex(nrecords_total,receiver_name, receiver_name_len, &
              nrecords_local, loc_begin, loc_end, rank, nproc)

  integer :: nrecords_total
  character(len=*) :: receiver_name
  integer :: receiver_name_len
  
  integer :: nrecords_local
  integer :: loc_begin, loc_end
  integer :: rank, nproc

  integer :: nrecords_per_proc
  character(len=20), allocatable :: receiver_name_array(:)
  integer :: dim_info

  integer :: loc_b(nproc), loc_e(nproc)
  integer :: i

  allocate(receiver_name_array(nrecords_total))

  call split_string(receiver_name, receiver_name_len, receiver_name_array, &
                        dim_info, '.')

  nrecords_per_proc=int(nrecords_total/nproc)

  do i=1,nproc
    loc_b(i) = nrecords_per_proc*(i-1)+1
    loc_e(i) = nrecords_per_proc*(i)
  enddo

  loc_e(nproc)=nrecords_total
  !if(rank.eq.(nproc-1)) then !last node
  !  loc_e()=nrecords_total
  !endif

  !further adjustment
  !if(rank.ne.0) then
  !  do while(trim(receiver_name_array(loc_begin)).eq.&
  !                      trim(receiver_name_array(loc_begin+1)))
  !    loc_begin=loc_begin+1
  !  enddo
  !endif

  do i=1, nproc-1
    do while(trim(receiver_name_array(loc_e(i))).eq.&
                        trim(receiver_name_array(loc_e(i)+1)))
      loc_e(i)=loc_e(i)+1
    enddo
  enddo

  do i=2, nproc
    loc_b(i)=loc_e(i-1)+1
  enddo

  loc_begin=loc_b(rank+1)
  loc_end=loc_e(rank+1)
  !print *,"loc_begin, loc_end:", loc_begin, loc_end
  nrecords_local=loc_end-loc_begin+1

end subroutine split_job_mpi_complex

end module asdf_subs
