module mod_dihedral
   !
   ! Created by Veselin Kolev <vesso.kolev@gmail.com>
   ! 20151012043512
   !
   use mod_t
   use mod_array

   implicit none

contains


   subroutine updateTorsionEnergy(pdihs,ETTot,pdihID)
   !
   ! After the rotation to the axis defined by the atoms j and k of the dihedral
   ! angle pdihs(pdihID)%currentAngle should be changed in the data set located
   ! handled into the memory. This routine should be called afterwards to
   ! calculate the torsion energy of the dihedral angle and update EETor.
   !
   type(pdih_t), dimension(:), intent(inout) :: pdihs
   real(kind=real64), intent(inout) :: ETTot
   integer(kind=int64), intent(in) :: pdihID

   ETTot=ETTot-pdihs(pdihID)%energy

   call calculateTorsionEnergy(pdihs,pdihID)

   ETTot=ETTot+pdihs(pdihID)%energy

   end subroutine updateTorsionEnergy


   subroutine calculateTotalTorsionEnergy(pdihs,ETTot)
   !
   ! Calculates the sum of all torsion energies that are already calculated
   ! by the subroutne calculateTorsionEnergies. This routine should be called
   ! only once - before starting the real simultion. Later the subroutine
   ! updateTorsionEnergy should be used in order to reduce the total
   ! computational time.
   !
   type(pdih_t), dimension(:), intent(in) :: pdihs
   real(kind=real64), intent(out) :: ETTot
   integer(kind=int64) :: i

   ETTot=0.0_real64

   do i=1,size(pdihs,1)
      ETTot=ETTot+pdihs(i)%energy
   end do

   end subroutine calculateTotalTorsionEnergy


   subroutine calculateDihedralAngles(atoms,pdihs)
   !
   ! Calculates the dihedral (torsion) angles corresponding to each and every
   ! dihedral definition given in the input files.
   !
   type(atom_t), dimension(:), intent(in) :: atoms
   type(pdih_t), dimension(:), intent(inout) :: pdihs
   integer(kind=int64) :: i

   do i=1,size(pdihs,1)
      call calculateDihedralAngle(atoms,pdihs,i)
   end do

   end subroutine calculateDihedralAngles


   subroutine calculateTorsionEnergies(pdihs)
   !
   ! Calculates the torsion energy corresponding to each and every of the
   ! dihedral angles defined in the input files. All dihedral angles have to
   ! be calculated!!! See subroutine calculateDihedralAngles.
   !
   type(pdih_t), dimension(:), intent(inout) :: pdihs
   integer(kind=int64) :: i

   do i=1,size(pdihs,1)
      call calculateTorsionEnergy(pdihs,i)
   end do

   end subroutine calculateTorsionEnergies


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


   subroutine calculateTorsionEnergy(pdihs,pdihID)
   !
   ! Calculates the torsion energy corresponding to a certain dihedral angle
   ! (selected by using its ID - pdihID). Note that the torsion energy is
   ! defined as sum(i=1,i=N) kd*[1+cos(mult*phi-phi0)]. Here kd, mult, and
   ! phi0 are supplied by the force field while phi is supplied by the
   ! subroutine calculateTorsionEnergy.
   !
   type(pdih_t), dimension(:), intent(inout) :: pdihs
   integer(kind=int64), intent(in) :: pdihID
   integer(kind=int64) :: i

   pdihs(pdihID)%energy=0.0_real64

   do i=1,size(pdihs(pdihID)%params,1)
      pdihs(pdihID)%energy=pdihs(pdihID)%energy+&
      pdihs(pdihID)%params(i)%kd*(1.0_real64+&
      cos(pdihs(pdihID)%params(i)%mult*pdihs(pdihID)%currentAngle-&
      pdihs(pdihID)%params(i)%phi0))
   end do

   end subroutine calculateTorsionEnergy

end module mod_dihedral
