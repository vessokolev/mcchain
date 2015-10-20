module mod_dihedral
   !
   ! Created by Veselin Kolev <vesso.kolev@gmail.com>
   ! 20151012043512
   !
   use mod_t
   use mod_array

   implicit none

contains


   subroutine calculateDihedralAngles(atoms,pdihs)
   !
   ! Calculates the dihedral (torsion) angles corresponding to each and every
   ! dihedral definition given in the input files.
   !
   type(atom_t), dimension(:), intent(in) :: atoms
   type(pdih_t), dimension(:), intent(inout) :: pdihs
   integer(kind=int64) :: i

   do i=1,pdihs_s
      call calculateDihedralAngle(atoms,pdihs,i)
   end do

   end subroutine calculateDihedralAngles


   subroutine calculateDihedralAngle(atoms,pdihs,pdihID)
   !
   ! Calculates the dihedral angle. If i,j,k,l is a sequence of vectors
   ! each one representing the position of an atom, the ternary i,j,k
   ! lies on a same plane - reference plane. The dihedral angle is defined
   ! then as the angle between that reference plane and the plane fixed by
   ! the ternary j,k,l. Note than the provided here routine calculates first
   ! the cosine of that angle and then by using the inverse function obtains
   ! the actual dihedral angle.
   ! This routine is supposed to be helpful only to handle the input data.
   ! During the runtime one should recalculate the angle as sum of the previous
   ! value and its change.
   !
   type(atom_t), dimension(:), intent(in) :: atoms
   type(pdih_t), dimension(:), intent(inout) :: pdihs
   integer(kind=int64), intent(in) :: pdihID
   real(kind=real64), dimension(3) :: m,n,ri,rj,rk,rl
   real(kind=real64), dimension(3) :: rji,rjk,rkl
   real(kind=real64) :: dummy,cos_
   integer(kind=int64) :: di,dj,dk,dl

   di=pdihs(pdihID)%member_id(1)
   dj=pdihs(pdihID)%member_id(2)
   dk=pdihs(pdihID)%member_id(3)
   dl=pdihs(pdihID)%member_id(4)

   ri=(/atoms(di)%x,atoms(di)%y,atoms(di)%z/)
   rj=(/atoms(dj)%x,atoms(dj)%y,atoms(dj)%z/)
   rk=(/atoms(dk)%x,atoms(dk)%y,atoms(dk)%z/)
   rl=(/atoms(dl)%x,atoms(dl)%y,atoms(dl)%z/)

   rji=rj-ri
   rjk=rj-rk
   rkl=rk-rl

   m=vectorCrossProduct(rji,rjk)
   n=vectorCrossProduct(-rjk,rkl)
   dummy=norm2(m)*norm2(n)

   cos_=dot_product(m,n)/dummy
   !
   ! Note that the dihedral angle may be calculated also by using the inverse
   ! of tan. In that case one should define the sine of the dihedral angle as
   ! sin_=norm2(vectorCrossProduct(m,n))/dummy
   !
   pdihs(pdihID)%currentAngle=acos(cos_)

   end subroutine calculateDihedralAngle


end module mod_dihedral
