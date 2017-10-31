! one uniform deviate
function dranf() result(u)
  implicit none
  real(kind=8) :: u
  call random_number(u)
  return
end function dranf
