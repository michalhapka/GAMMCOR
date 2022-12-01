subroutine pivoted_cholesky( A, rank, tol, ndim, U)
!
! A = U**T * U
!
! matrix A is destroyed inside this subroutine
! Cholesky vectors are stored in U
! dimension of U: U(1:rank, 1:n)
! U is allocated inside this subroutine
! rank is the number of Cholesky vectors depending on tol
!
integer :: ndim
integer, intent(inout)                                        :: rank
double precision, dimension(ndim, ndim), intent(inout)        :: A
double precision, dimension(ndim, rank), intent(out)          :: U
double precision, intent(in)                                  :: tol

integer, dimension(:), allocatable          :: piv
double precision, dimension(:), allocatable :: work
character, parameter :: uplo = "U"
integer :: N, LDA
integer :: info
integer :: k, l, rank0
external :: dpstrf

rank0 = rank
N = size(A, dim=1)
LDA = N
allocate(piv(N))
allocate(work(2*N))
call dpstrf(uplo, N, A, LDA, piv, rank, tol, work, info)

if (rank > rank0) then
  print *, 'Bug: rank > rank0 in pivoted cholesky. Increase rank before calling'
  stop
end if

do k = 1, N
  A(k+1:, k) = 0.00D+0
end do
! TODO: It should be possible to use only one vector of size (1:rank) as a buffer
! to do the swapping in-place
U = 0.00D+0
do k = 1, N
  l = piv(k)
  U(l, :) = A(1:rank, k)
end do

end subroutine pivoted_cholesky

