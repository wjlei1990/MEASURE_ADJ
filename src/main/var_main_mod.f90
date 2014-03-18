module var_main

  logical :: WRITE_ADJ_ASDF
  logical :: ROTATE_COMP, WRITE_NORMAL_OUTPUT

  character(len=150) :: OBSD_FILE, SYNT_FILE, SYNT_PHYDISP_FILE
  character(len=150) :: WIN_DIR

  character(len=150) :: MEASURE_ADJ_OUTDIR
  
  real :: MIN_PERIOD, MAX_PERIOD

  integer :: weighting_option

  logical :: USE_PHYDISP

end module var_main
