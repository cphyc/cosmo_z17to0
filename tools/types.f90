module types
  use, intrinsic :: iso_fortran_env
  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(15, 307)

  integer, parameter :: isp = 4
  integer, parameter :: idp = 8
end module types
