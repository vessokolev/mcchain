module mod_rot

   use mod_t
   use mod_distance

   implicit none

contains

   subroutine rotateGroupOfAtoms(atoms,pdihs,distM1D,atomIDArray,pdih,angle)
   type(atom_t), dimension(:), intent(inout) :: atoms
   type(pdih_t), dimension(:), intent(in) :: pdihs
   real(kind=real64), dimension(:), intent(inout) :: distM1D
   integer(kind=int64), dimension(:), intent(in) :: atomIDArray
   type(pdih_t), intent(in) :: pdih
   real(kind=real64), intent(in) :: angle
   integer(kind=int64) :: i
   real(kind=real64), dimension(3,3) :: rotMatrix
   real(kind=real64), dimension(3) :: tmp,p1,p2

   p1(1)=atoms(pdih%member_id(2))%x
   p1(2)=atoms(pdih%member_id(2))%y
   p1(3)=atoms(pdih%member_id(2))%z

   p2(1)=atoms(pdih%member_id(3))%x
   p2(2)=atoms(pdih%member_id(3))%y
   p2(3)=atoms(pdih%member_id(3))%z

   rotMatrix=getRotMatrix(p1,p2,angle)

   call translateAtomCoords(atoms,p1)

   do i=1,size(atomIDArray,1)
      tmp=(/atoms(atomIDArray(i))%x,atoms(atomIDArray(i))%y,atoms(atomIDArray(i))%z/)
      tmp=rotatePoint(tmp,rotMatrix)
      atoms(atomIDArray(i))%x=tmp(1)
      atoms(atomIDArray(i))%y=tmp(2)
      atoms(atomIDArray(i))%z=tmp(3)
      call updateDistanceMatrix1Atom(atoms,distM1D,i)
   end do

   call translateAtomCoords(atoms,-p1)

   end subroutine rotateGroupOfAtoms


   subroutine translateAtomCoords(atoms,transVector)
   type(atom_t), dimension(:), intent(inout) :: atoms
   real(kind=real64), dimension(3), intent(in) :: transVector
   integer(kind=int64) :: i
   
   do i=1,size(atoms,1)
      atoms(i)%x=atoms(i)%x-transVector(1)
      atoms(i)%y=atoms(i)%y-transVector(2)
      atoms(i)%z=atoms(i)%z-transVector(3)
   end do

   end subroutine translateAtomCoords


   function getRotMatrix(point1,point2,angle) result(matrix)
   real(kind=real64), dimension(3), intent(in) :: point1,point2
   real(kind=real64), intent(in) :: angle
   real(kind=real64), dimension(3,3) :: matrix
   real(kind=real64), dimension(3) :: norm
   real(kind=real64) :: cos_,sin_,omcos_,n12,n13,n23

   norm=(/point2(1)-point1(1),point2(2)-point1(2),point2(3)-point1(3)/)
   norm(:)=norm(:)/norm2(norm)

   cos_=cos(angle)
   sin_=sin(angle)

   omcos_=1-cos_
   n12=norm(1)*norm(2)
   n13=norm(1)*norm(3)
   n23=norm(2)*norm(3)

   matrix(1,1)=cos_+norm(1)*norm(1)*omcos_
   matrix(1,2)=n12*omcos_-norm(3)*sin_
   matrix(1,3)=n13*omcos_+norm(2)*sin_

   matrix(2,1)=n12*omcos_+norm(3)*sin_
   matrix(2,2)=cos_+norm(2)*norm(2)*omcos_
   matrix(2,3)=n23*omcos_-norm(1)*sin_

   matrix(3,1)=n13*omcos_-norm(2)*sin_
   matrix(3,2)=n23*omcos_+norm(1)*sin_
   matrix(3,3)=cos_+norm(3)*norm(3)*omcos_

   end function getRotMatrix


   function rotatePoint(point,rotMatrix) result(rotated)
   real(kind=real64), dimension(3), intent(in) :: point
   real(kind=real64), dimension(3,3), intent(in) :: rotMatrix
   real(kind=real64), dimension(3) :: rotated
   integer :: i,j

   do i=1,3
      rotated(i)=0.0d0
      do j=1,3
         rotated(i)=rotated(i)+rotMatrix(i,j)*point(j)
      end do
   end do

   end function rotatePoint


   function translatePoint(point,bench) result(translated)
   real(kind=real64), dimension(3), intent(in) :: point,bench
   real(kind=real64), dimension(3) :: translated

   translated(:)=point(:)-bench(:)

   end function


   function inv(A) result(Ainv)
   real(kind=real64), dimension(:,:), intent(in) :: A
   real(kind=real64), dimension(size(A,1),size(A,2)) :: Ainv

   real(kind=real64), dimension(size(A,1)) :: work  ! work array for LAPACK
   integer, dimension(size(A,1)) :: ipiv   ! pivot indices
   integer(kind=int8) :: n, info

   ! External procedures defined in LAPACK
   external DGETRF
   external DGETRI

   ! Store A in Ainv to prevent it from being overwritten by LAPACK
   Ainv = A
   n = size(A,1)

   ! DGETRF computes an LU factorization of a general M-by-N matrix A
   ! using partial pivoting with row interchanges.
   call DGETRF(n, n, Ainv, n, ipiv, info)

   if (info /= 0) then
      stop 'Matrix is numerically singular!'
   end if

   ! DGETRI computes the inverse of a matrix using the LU factorization
   ! computed by DGETRF.
   call DGETRI(n, Ainv, n, ipiv, work, n, info)

   if (info /= 0) then
      stop 'Matrix inversion failed!'
   end if
   end function inv


end module mod_rot
