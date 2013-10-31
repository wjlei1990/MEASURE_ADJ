module var_main

  logical :: RUN_FLEXWIN, RUN_MEASURE_ADJ, WRITE_ADJ_ASDF
  logical :: ROTATE_COMP, WRITE_NORMAL_OUTPUT

  character(len=150) :: OBSD_FILE, SYNT_FILE, ADJ_FILE
  character(len=150) :: WIN_FILE

  character(len=150) :: FLEXWIN_OUTDIR, MEASURE_ADJ_OUTDIR

  integer :: weighting_option

end module var_main
