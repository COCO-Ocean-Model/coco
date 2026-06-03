#define READ_NAMELIST(NML)                                       \
  call rewnml(ifpar, jfpar);                                     \
  read(ifpar, NML, iostat=istat);                                \
  call cstnml(jfpar, __FILE__, __LINE__, istat);                 \
  write(jfpar, NML)
