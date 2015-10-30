module mod_rot
   ! Created by Vesselin Kolev <vesso.kolev@gmail.com>
   ! 20151008231902
   use mod_t
   use mod_array
   use mod_distance

   implicit none

contains

   subroutine rotateGroupOfAtoms(atoms,pdihs,distM1D,atomIDArray,pdih,angle)
   !
   ! This subroutine rotates a group of atoms with respect to a given axis
   ! (the axis is defined by using the j and k atoms of a certain dihedral
   ! angle) by implementing Rodrigues' formula.
   !
   ! Here atomIDArray supplies the indexes of the atoms that have to be rotated.
   ! The angle (angle) is given in radians (not in degrees!!!). After each
   ! rotation the subroutine updates the distance matrix distM1D.
   !
   type(atom_t), dimension(:), intent(inout) :: atoms
   type(pdih_t), dimension(:), intent(in) :: pdihs
   real(kind=real32), dimension(:), intent(inout) :: distM1D
   integer(kind=int64), dimension(:), intent(in) :: atomIDArray
   type(pdih_t), intent(in) :: pdih
   real(kind=real32), intent(in) :: angle
   integer(kind=int64) :: i
   real(kind=real32), dimension(3,3) :: rotMatrix
   real(kind=real32), dimension(3) :: tmp,p1,p2

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
   !
   ! This subroutine translates the atomic coordinates with respect to
   ! a point given by transVector. Note that the subroutine translates
   ! ALL atomic coordinates (the coordinates of all atoms in the set).
   !
   type(atom_t), dimension(:), intent(inout) :: atoms
   real(kind=real32), dimension(3), intent(in) :: transVector
   integer(kind=int64) :: i
   
   do i=1,size(atoms,1)
      atoms(i)%x=atoms(i)%x-transVector(1)
      atoms(i)%y=atoms(i)%y-transVector(2)
      atoms(i)%z=atoms(i)%z-transVector(3)
   end do

   end subroutine translateAtomCoords


   function getRotMatrix(point1,point2,angle) result(matrix)
   !
   ! Calculates the rotation matrix elements that allows the application of the
   ! Rodrigues' formula. For more details see:
   ! https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
   !
   real(kind=real32), dimension(3), intent(in) :: point1,point2
   real(kind=real32), intent(in) :: angle
   real(kind=real32), dimension(3,3) :: matrix
   real(kind=real32), dimension(3) :: norm
   real(kind=real32) :: cos_,sin_,omcos_,n12,n13,n23

   norm=(/point2(1)-point1(1),point2(2)-point1(2),point2(3)-point1(3)/)
   norm(:)=norm(:)/norm02(norm)

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
   !
   ! Rotates the coordinates of a point in 3D space by implementing the
   ! Rodrugues' formula. The rotation matrix is supplied by the subroutine
   ! getRotMatrix.
   !
   real(kind=real32), dimension(3), intent(in) :: point
   real(kind=real32), dimension(3,3), intent(in) :: rotMatrix
   real(kind=real32), dimension(3) :: rotated
   integer :: i,j

   do i=1,3
      rotated(i)=0.0d0
      do j=1,3
         rotated(i)=rotated(i)+rotMatrix(i,j)*point(j)
      end do
   end do

   end function rotatePoint

end module mod_rot
