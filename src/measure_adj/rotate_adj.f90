module rotate_subs

implicit none
  double precision, parameter :: PI=3.1415926

contains

subroutine rotate_adj(my_asdf, my_asdf_rotate)
      
    use asdf_data
    use asdf_read_subs
    implicit none

    type(asdf_event) :: my_asdf, my_asdf_rotate

    integer :: i, j, k
    integer :: num_stations, sta_index

    integer :: nrecords

    logical :: z_exist, t_exist, r_exist
    integer :: loc_z, loc_r, loc_t, loc_e, loc_n
    integer :: nn, nn_r, nn_z, nn_t

    real :: azm,bzm,ddg,dkm
    real :: receiver_lat, receiver_lo, receiver_dpt
    real :: costh, sinth
    real(kind=8) :: cmt_lat2, stlat2, rlo, delta, azep, azst

    real, allocatable :: zdata(:), rdata(:), tdata(:), ndata(:), edata(:)

    character(len=32), allocatable :: sta_list(:), nw_list(:)
    character(len=32) :: sta, nw

    real :: rad, fl, ecc
    real :: tt=0, dtt=0
    real :: baz
    integer :: icomp

    double precision :: rotang, rotaz
    integer :: nerr

    real :: TORAD

    TORAD = 3.1415926/180.

    nrecords=my_asdf%nrecords

    print *, "******************"
    print *, "begin rotate..."
    allocate(sta_list(my_asdf%nrecords))
    allocate(nw_list(my_asdf%nrecords))

    !********************************
    !first see how many stations in the my_asdf. then allocate the rotated data
    !structure by 3*num_stations
    sta_index=0
    sta_list(1)=my_asdf%receiver_name_array(1)
    nw_list(1)=my_asdf%network_array(1)
    do i=1, my_asdf%nrecords
      sta=my_asdf%receiver_name_array(i)
      nw=my_asdf%network_array(i)
      if(sta_exist(sta, nw, sta_list, nw_list, sta_index))then
        cycle
      else
       sta_index=sta_index+1
       sta_list(sta_index)=my_asdf%receiver_name_array(i)
       nw_list(sta_index)=my_asdf%network_array(i)
      endif
    enddo
    num_stations=sta_index

    call init_asdf_data(my_asdf_rotate,3*num_stations, .false.)

    print *, "Num_stations:", num_stations

    !*********************************
    !recalculate the az and baz
    rad=PI/180.0
    fl=0.00335293
    ecc=(1.-fl)**2.
    !cmt_lat2=atan(ecc*tan(cmt_lat*rad))*(180/PI)
    do i=1, my_asdf%nrecords
      print *, "ii:",i
      !stlat2=atan(ecc*tan(my_asdf%receiver_lat(i)*rad))*(180/PI)
      call distaz(sngl(my_asdf%event_lat(i)), sngl(my_asdf%event_lo(i)), &
            sngl(my_asdf%receiver_lat(i)), sngl(my_asdf%receiver_lo(i)),&
            azm, bzm, ddg, dkm)
      print *, "info:", my_asdf%event_lat(i), my_asdf%event_lo(i), &
        my_asdf%receiver_lat(i), my_asdf%receiver_lo(i), azm, bzm, ddg, dkm
    !  call SPH_AZI(cmt_lat2,cmt_lon,stlat2,rlo,delta,azep,azst)
      my_asdf%ev_to_sta_AZ(i)=azm
      my_asdf%sta_to_ev_AZ(i)=bzm
    enddo

    print *, "Calculate azi done"

    !**********************************
    !loop over sta list to find the specific stations
    do i=1, num_stations
      print *, "Station and network:", trim(sta_list(i)), trim(nw_list(i))
      !locate three components
      call locate_record(sta_list(i), nw_list(i), "Z", my_asdf%receiver_name_array, &
        my_asdf%network_array, my_asdf%component_array, nrecords, loc_z)
      call locate_record(sta_list(i), nw_list(i), "R", my_asdf%receiver_name_array, &
        my_asdf%network_array, my_asdf%component_array, nrecords, loc_r)
      call locate_record(sta_list(i), nw_list(i), "T", my_asdf%receiver_name_array, &
        my_asdf%network_array, my_asdf%component_array, nrecords, loc_t)
      print *, "loc_r, t, z:", loc_r, loc_t, loc_z
      
      !if found, set the npionts; otherwise, set it to zero
      if(loc_r>0)then
        nn_r=my_asdf%npoints(loc_r)
      else
        nn_r=0
      endif
      if(loc_t>0)then
        nn_t=my_asdf%npoints(loc_t)
      else
        nn_t=0
      endif
      if(loc_z>0)then
        nn_z=my_asdf%npoints(loc_z)
      else
        nn_z=0
      endif

      !set bzm, tt, dtt for rotate
      if(loc_z.gt.0)then
        receiver_lat=my_asdf%receiver_lat(loc_z)
        receiver_lo=my_asdf%receiver_lo(loc_z)
        receiver_dpt=my_asdf%receiver_dpt(loc_z)

        azm=my_asdf%ev_to_sta_AZ(loc_z)
        bzm=my_asdf%sta_to_ev_AZ(loc_z)
        tt=my_asdf%begin_value(loc_z)
        dtt=my_asdf%sample_rate(loc_z)
      elseif(loc_r.gt.0)then
        receiver_lat=my_asdf%receiver_lat(loc_r)
        receiver_lo=my_asdf%receiver_lo(loc_r)
        receiver_dpt=my_asdf%receiver_dpt(loc_r)

        azm=my_asdf%ev_to_sta_AZ(loc_r)
        bzm=my_asdf%sta_to_ev_AZ(loc_r)
        tt=my_asdf%begin_value(loc_r)
        dtt=my_asdf%sample_rate(loc_r)
      elseif(loc_t.gt.0)then
        receiver_lat=my_asdf%receiver_lat(loc_t)
        receiver_lo=my_asdf%receiver_lo(loc_t)
        receiver_dpt=my_asdf%receiver_dpt(loc_t)

        azm=my_asdf%ev_to_sta_AZ(loc_t)
        bzm=my_asdf%sta_to_ev_AZ(loc_t)
        tt=my_asdf%begin_value(loc_t)
        dtt=my_asdf%sample_rate(loc_t)
      endif

      !============================================
      !deal with North and East component first
      print *, "Rotate East and North"
      if((loc_r.eq.0).and.(loc_t.eq.0))then
        !both north and east component is missing. Set them to zero and their
        !length same to Z component
        print *, "Both component missing"
        nn=nn_z
        allocate(my_asdf_rotate%records(3*i-2)%record(nn))
        allocate(my_asdf_rotate%records(3*i-1)%record(nn))
        my_asdf_rotate%records(3*i-2)%record(:)=0.
        my_asdf_rotate%records(3*i-1)%record(:)=0.
      elseif((loc_r.eq.0).or.(loc_t.eq.0))then
        !if one of the components is missing, set them to zero and their length
        !same as the other one
        print *, "One of the components is missing"
        nn=max(nn_r, nn_t)
        allocate(my_asdf_rotate%records(3*i-2)%record(nn))
        allocate(my_asdf_rotate%records(3*i-1)%record(nn))
        my_asdf_rotate%records(3*i-2)%record(:)=0.
        my_asdf_rotate%records(3*i-1)%record(:)=0.
      else
        !if both of them exist, then rotate them...
        !It is the only situation we need to rotate
        print *, "Both components exist"
        nn=min(nn_r, nn_t)
        allocate(rdata(nn))
        allocate(tdata(nn))
        allocate(edata(nn))
        allocate(ndata(nn))
        allocate(my_asdf_rotate%records(3*i-2)%record(nn))
        allocate(my_asdf_rotate%records(3*i-1)%record(nn))
        rdata(1:nn)=my_asdf%records(loc_r)%record(1:nn)
        tdata(1:nn)=my_asdf%records(loc_t)%record(1:nn)
        !bzm=326.2812
        rotaz=bzm
        print *, "bzm, rotaz:", bzm, rotaz
        !*********************
        !my own rotate subroutine
        costh = cos(TORAD*rotaz)
        sinth = sin(TORAD*rotaz)
        print *, "costh, sinth:", costh, sinth
        edata(1:nn) = -costh * tdata(1:nn) - sinth * rdata(1:nn)
        ndata(1:nn) = sinth * tdata(1:nn) - costh * rdata(1:nn)
        !call rotate(ndata, edata, nn, rotaz, .true., .true., rdata, tdata)
        my_asdf_rotate%records(3*i-2)%record(1:nn)= edata(1:nn)
        my_asdf_rotate%records(3*i-1)%record(1:nn)= ndata(1:nn)
        deallocate(rdata)
        deallocate(tdata)
        deallocate(edata)
        deallocate(ndata)
      endif

      !==========================================
      print *, "Copy Z component"
      if(loc_z>0)then
        !if Z exists, then just copy the data
        allocate(my_asdf_rotate%records(3*i)%record(nn_z))
        my_asdf_rotate%records(3*i)%record(1:nn_z) = &
                     my_asdf%records(loc_z)%record(1:nn_z)
      else
        !if Z does not exist, then set the array length same as N and E, and
        !record to zero
        nn_z=nn
        allocate(my_asdf_rotate%records(3*i)%record(nn_z))
        my_asdf_rotate%records(3*i)%record(1:nn_z) = 0.0
      endif


      print *, "nn", nn, bzm, tt, dtt

      !************************************
      print *, "fill other info in rotate asdf structure"
      !fill other info
      my_asdf_rotate%component_array(3*i-2) = "LHE"
      my_asdf_rotate%component_array(3*i-1) = "LHN"
      my_asdf_rotate%component_array(3*i) = "LHZ"
      print *, "locations: ", loc_e, loc_n, loc_z
      do icomp=0,2
        my_asdf_rotate%npoints(3*i-icomp) = nn
        my_asdf_rotate%sample_rate(3*i-icomp) = dtt
        my_asdf_rotate%begin_value(3*i-icomp) = tt
        my_asdf_rotate%receiver_name_array(3*i-icomp) = sta_list(i)
        my_asdf_rotate%network_array(3*i-icomp) = nw_list(i)

        !set time to cmt time
        !my_asdf_rotate%gmt_year(3*i-icomp) = gmt_year
        !my_asdf_rotate%gmt_day(3*i-icomp) = gmt_day
        !my_asdf_rotate%gmt_hour(3*i-icomp) = gmt_hour
        !my_asdf_rotate%gmt_min(3*i-icomp) = gmt_min
        !my_asdf_rotate%gmt_sec(3*i-icomp) = gmt_sec
        !my_asdf_rotate%gmt_msec(3*i-icomp) = gmt_msec

        !set event to cmt
        !my_asdf_rotate%event_lat(3*i-icomp) = cmt_lat
        !my_asdf_rotate%event_lo(3*i-icomp) = cmt_lon
        !my_asdf_rotate%event_dpt(3*i-icomp) = cmt_depth

        !set receiver location info
        my_asdf_rotate%receiver_lat(3*i-icomp) = receiver_lat
        my_asdf_rotate%receiver_lo(3*i-icomp) = receiver_lo
        my_asdf_rotate%receiver_dpt(3*i-icomp) =receiver_dpt

        !set azimuth and back azimuth
        my_asdf_rotate%ev_to_sta_AZ(3*i-icomp) = azm
        my_asdf_rotate%sta_to_ev_AZ(3*i-icomp) = bzm

        !set receiver id to null
        my_asdf_rotate%receiver_id_array(3*i-icomp) = ""
      enddo
      !if (loc_e>0) then
      !  my_asdf_rotate%receiver_id_array(3*i-2) = &
      !          my_asdf%receiver_id_array(loc_e)
      !else
      !  my_asdf_rotate%receiver_id_array(3*i-2) = ""
                
      !endif
      !if (loc_n>0) then
      !  my_asdf_rotate%receiver_id_array(3*i-1) = &
      !          my_asdf%receiver_id_array(loc_n)
      !else
      !  my_asdf_rotate%receiver_id_array(3*i-1) =""
      !endif
      !if (loc_z>0) then
      !  my_asdf_rotate%receiver_id_array(3*i) = &
      !          my_asdf%receiver_id_array(loc_z)
      !else
      !  my_asdf_rotate%receiver_id_array(3*i)=""
      !endif

     !print *, "receiver_id:", &
     !   trim(my_asdf_rotate%receiver_id_array(1)), &
     !   trim(my_asdf_rotate%receiver_id_array(2)), &
     !   trim(my_asdf_rotate%receiver_id_array(3))

     !call wsac1("E.test", edata, nn, 0, 1.0, nerr)
     !call wsac1("N.test", ndata, nn, 0, 1.0, nerr)
     !call wsac1("R.test", rdata, nn, 0, 1.0, nerr)
     !call wsac1("T.test", tdata, nn, 0, 1.0, nerr)
     !call wsac1("Z.test", zdata, nn, 0, 1.0, nerr)
    enddo

    !fill common header infor
    my_asdf_rotate%event=my_asdf%event

    !********************************************

end subroutine rotate_adj

!======================================================================
!> Interpolates a seismogram to a given sample rate
!! \param syn The seismogram to interpolate
!! \param t1 The start time of the seismogram
!! \param dt1 The sample rate of the seismogram
!! \param npt1 The number of points sampled in the seismogram
        

  logical function sta_exist(receiver, network, sta, ntw, sta_index)

    character(len=*) :: receiver, network
    character(len=*) :: sta(:), ntw(:)
    integer :: sta_index
    integer :: i

    if(sta_index.eq.0)then
      sta_exist=.false.
      return
    endif

    do i=1, sta_index
      if( (trim(sta(i)).eq.trim(receiver)) .and. &
                            (trim(ntw(i)).eq.trim(network)) ) then
        sta_exist=.true.
        return
      endif
    enddo

    sta_exist=.false.
    return

  end function sta_exist

  !subroutine that found the specific record
  !return the loc
  subroutine  locate_record(receiver, network, component, sta_array, &
        nw_array, comp_array, dim_array, loc) 

    character(len=*) :: receiver, network, component
    character(len=*) :: sta_array(:), nw_array(:), comp_array(:)
    integer :: dim_array, loc

    integer :: i

    if(dim_array.eq.0)then
      loc=0
      return
    endif

    loc=0
    do i=1, dim_array
      if( (trim(sta_array(i)).eq.trim(receiver)) &
          .and. (trim(nw_array(i)).eq.trim(network)) )then
          if(trim(component).eq.trim(comp_array(i)(3:3)))then
            loc=i
          endif
      endif
    enddo

  end subroutine locate_record

end module rotate_subs
