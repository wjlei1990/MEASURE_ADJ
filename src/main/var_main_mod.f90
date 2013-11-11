module var_main

  logical :: WRITE_ADJ_ASDF
  logical :: ROTATE_COMP, WRITE_NORMAL_OUTPUT

  character(len=150) :: OBSD_FILE, SYNT_FILE
  character(len=150) :: WIN_DIR

  character(len=150) :: MEASURE_ADJ_OUTDIR

  integer :: weighting_option

end module var_main
