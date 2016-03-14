module convert
  implicit none

  private

  public :: itos, stoi

contains
  FUNCTION stoi(str) RESULT(num)
    CHARACTER(len=*), INTENT(in) :: str
    INTEGER                      :: num
    integer :: tmp

    write(tmp, '(i8)') str
    num = tmp

  END FUNCTION stoi

  FUNCTION itos(num) RESULT(strout)
    INTEGER, INTENT(in)          :: num
    CHARACTER(len=10)            :: strout, tmp2
    integer :: tmp
    tmp = num

    read(tmp2, '(i8)') tmp
    strout = tmp2
  END FUNCTION itos
end module convert
