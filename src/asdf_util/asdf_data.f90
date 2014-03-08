module asdf_data

  type asdf_record
    real, allocatable :: record(:)
  end type asdf_record

  type asdf_event
    ! scalars
    character(len=13)     :: event
    real, allocatable     :: event_lat(:), event_lo(:), event_dpt(:)

    !size info
    integer           :: nrecords
    integer           :: nreceivers

    !Processing info
    real              :: min_period, max_period

    !time info
    integer, allocatable    :: gmt_year(:), gmt_day(:), gmt_hour(:)
    integer, allocatable    :: gmt_min(:), gmt_sec(:), gmt_msec(:)

    !seismic record info
    integer, allocatable    :: npoints(:)
    real, allocatable       :: receiver_lat(:), receiver_lo(:)
    real, allocatable       :: receiver_el(:),  receiver_dpt(:)
    real, allocatable       :: begin_value(:),  end_value(:)
    real, allocatable       :: cmp_azimuth(:),  cmp_incident_ang(:)
    real, allocatable       :: sample_rate(:),  scale_factor(:)

    real, allocatable       :: ev_to_sta_AZ(:), sta_to_ev_AZ(:)
    real, allocatable       :: great_circle_arc(:) 
    real, allocatable       :: dist(:)
    real, allocatable       :: P_pick(:), S_pick(:)

    !string, name, etc
    !integer    :: receiver_name_len, network_len
    !integer    :: component_len,receiver_id_len
    !character(len=:),allocatable :: receiver_name
    !character(len=:),allocatable :: network
    !character(len=:),allocatable :: component
    !character(len=:),allocatable :: receiver_id

    character(len=20),allocatable :: receiver_name_array(:), network_array(:)
    character(len=20),allocatable :: component_array(:), receiver_id_array(:)

    !seismograms
    type (asdf_record), allocatable :: records(:)

  end type asdf_event

end module asdf_data
