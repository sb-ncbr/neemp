
!! Test programm for Subroutine minimize_with_newuoa.
!! To compile and run the program with gfortran do on the console:
!! - 1) gfortran -c newuoa_module.f95
!! - 2) gfortran -c test_newuoa_module.f95
!! - 3) gfortran newuoa_module.o test_newuoa_module.o
!!
!! Or in one step:
!! - gfortran newuoa_module.f95 test_newuoa_module.f95
!!
!! This produces an executable a.out (linux) or a.exe (windows).


PROGRAM test_newuoa_module

   use newuoa_module, only : minimize_with_newuoa
   
   integer :: I, N, NPT
   real(8) :: RHOEND, RHOBEG
   real(8), dimension(:), allocatable :: X
   
      RHOEND=1.0D-6
      DO N=2,8,2
         allocate(X(N))
         NPT=2*N+1
         forall(I=1:N) X(I)=DFLOAT(I)/DFLOAT(N+1)
         RHOBEG=0.2D0*X(1)
         write(*,fmt='(/,a,I3,a,I5)') 'Results with N =', N, &
               ' and NPT =', NPT
         call minimize_with_newuoa(CALFUN, X, RHOBEG, RHOEND, &
                              i_opt_NPT=NPT, i_opt_IPRINT=2, i_opt_MAXFUN=50000)
         deallocate(X)
      END DO
   

 CONTAINS


      SUBROUTINE CALFUN(X,F)
      IMPLICIT REAL*8 (A-H,O-Z)
      real(kind=8), intent(in), dimension(:) :: X
      real(kind=8), intent(out) :: F

      real(kind=8), dimension(size(X)+1,size(X)) :: Y
      N = size(X)
      DO 10 J=1,N
      Y(1,J)=1.0D0
   10 Y(2,J)=2.0D0*X(J)-1.0D0
      DO 20 I=2,N
      DO 20 J=1,N
   20 Y(I+1,J)=2.0D0*Y(2,J)*Y(I,J)-Y(I-1,J)
      F=0.0D0
      NP=N+1
      IW=1
      DO 40 I=1,NP
      SUM=0.0D0
      DO 30 J=1,N
   30 SUM=SUM+Y(I,J)
      SUM=SUM/DFLOAT(N)
      IF (IW .GT. 0) SUM=SUM+1.0D0/DFLOAT(I*I-2*I)
      IW=-IW
   40 F=F+SUM*SUM
      RETURN
      END SUBROUTINE



END PROGRAM test_newuoa_module