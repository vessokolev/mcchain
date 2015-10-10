module mod_input

   use mod_t
   use mod_array

contains


   subroutine readInput(files,atoms,bonds,pdihs,eps_M,sigma_M)
   type(file_t), intent(in) :: files
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   type(pdih_t), dimension(:), allocatable, intent(out) :: pdihs
   real(kind=real64), dimension(:,:), allocatable :: eps_M
   real(kind=real64), dimension(:,:), allocatable :: sigma_M

   call readAtoms(files,atoms)
   call readBonds(files,bonds)
   call readPropDihedrals(files,pdihs)
   call updateAtomAddBondsArray(atoms,bonds)

   call calculateDihedralAnglesEnergies(atoms,pdihs)
   call putAtomsInEpsGroups(atoms,eps_M)
   call putAtomsInSigmaGroups(atoms,sigma_M)

   end subroutine readInput



   subroutine readAtoms(files,atoms)
   type(file_t), intent(in) :: files
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   integer(kind=int64) :: i,numAtoms

   open(unit=1,file=files%atoms,form="FORMATTED")

   read(1,fmt="(i10)") numAtoms

   allocate(atoms(numAtoms))

   do i=1,numAtoms
      read(1,fmt="(i10,1x,a5,1x,i10,1x,a5,1x,a5,1x,i10,1x,f12.6,1x,&
                   f12.6,1x,f12.6,1x,f12.6,1x,f12.6,1x,f12.6,1x,f12.6)")&
      atoms(i)%anum,atoms(i)%atype,atoms(i)%resnum,atoms(i)%resname,&
      atoms(i)%aname,atoms(i)%cgrp,atoms(i)%charge,atoms(i)%mass,&
      atoms(i)%sigma,atoms(i)%eps,atoms(i)%x,atoms(i)%y,atoms(i)%z
   end do

   close(1)

   end subroutine readAtoms


   subroutine readBonds(files,bonds)
   type(file_t), intent(in) :: files
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   integer(kind=int64) :: i,numBonds

   open(unit=1,file=files%bonds,form="FORMATTED")

   read(1,fmt="(i15)") numBonds

   allocate(bonds(numBonds))

   do i=1,numBonds
      read(1,fmt="(i10,1x,i10,1x,a5,1x,a5)")&
      bonds(i)%left_a,bonds(i)%right_a,&
      bonds(i)%left_a_t,bonds(i)%right_a_t
   end do

   close(1)

   end subroutine readBonds


   subroutine readPropDihedrals(files,pdihs)
   type(file_t), intent(in) :: files
   type(pdih_t), dimension(:), allocatable, intent(out) :: pdihs
   integer(kind=int64) :: i,j,numPDihrs,dummy
   integer(kind=int16) :: numDihs


   open(unit=1,file=files%pdihsDefs,form="FORMATTED")
   open(unit=2,file=files%pdihsParams,form="FORMATTED")

   read(1,fmt="(i10)") numPDihrs

   allocate(pdihs(numPDihrs))

   do i=1,numPDihrs
      read(1,fmt="(i10,1x,i10,1x,i10,1x,i10,1x,i10,1x,i2)") &
      pdihs(i)%id,&
      pdihs(i)%member_id(1),pdihs(i)%member_id(2),pdihs(i)%member_id(3),&
      pdihs(i)%member_id(4),numDihs
      allocate(pdihs(i)%params(numDihs))
      do j=1,numDihs
         read(2,fmt="(i10,1x,f12.6,1x,f12.6,1x,i2)") dummy,&
         pdihs(i)%params(j)%phi0,pdihs(i)%params(j)%kd,&
         pdihs(i)%params(j)%mult
      end do
   end do

   close(2)
   close(1)

   end subroutine readPropDihedrals


   subroutine updateAtomAddBondsArray(atoms,bonds)
   type(atom_t), dimension(:), intent(inout) :: atoms
   type(bond_t), dimension(:), intent(in) :: bonds
   integer(kind=int64) :: i
   logical, dimension(size(atoms,1)) :: flag

   flag(:)=.True.

   do i=1,size(atoms,1)
      allocate(atoms(i)%bonds(1))
   end do

   do i=1,size(bonds,1)
         if (flag(bonds(i)%left_a)) then
            flag(bonds(i)%left_a)=.False.
            atoms(bonds(i)%left_a)%bonds(1)=bonds(i)%right_a
         else
            call extend1DIntArray(atoms(bonds(i)%left_a)%bonds,&
                                  bonds(i)%right_a)
         end if
         if (flag(bonds(i)%right_a)) then
            flag(bonds(i)%right_a)=.False.
            atoms(bonds(i)%right_a)%bonds(1)=bonds(i)%left_a
         else
            call extend1DIntArray(atoms(bonds(i)%right_a)%bonds,&
                                  bonds(i)%left_a)
         end if
   end do

   end subroutine updateAtomAddBondsArray


!   subroutine extend1DIntArray(array,newElement)
!   integer(kind=int64), dimension(:), allocatable, intent(inout) :: array
!   integer(kind=int64), intent(in) :: newElement
!   integer(kind=int64), dimension(:), allocatable :: temp
!   integer(kind=int64), dimension(1) :: array_s
!   integer(kind=int64) :: i

!   array_s=shape(array)

!   call move_alloc(array,temp)

!   array_s(1)=array_s(1)+1

!   allocate(array(array_s(1)))

!   array(1:array_s(1)-1)=temp(:)

!   deallocate(temp)

!   array(array_s(1))=newElement

!   end subroutine extend1DIntArray


   function calculateDihedralAngle(indx,atoms,pdihs) result(angle)
   integer(kind=int64), intent(in) :: indx
   type(atom_t), dimension(:), intent(in) :: atoms
   type(pdih_t), dimension(:), intent(in) :: pdihs
   real(kind=real64), dimension(3) :: m,n,ri,rj,rk,rl
   real(kind=real64), dimension(3) :: rji,rjk,rkj,rkl
   real(kind=real64) :: angle,dummy,sin_,cos_
   integer(kind=int64) :: di,dj,dk,dl

   di=pdihs(indx)%member_id(1)
   dj=pdihs(indx)%member_id(2)
   dk=pdihs(indx)%member_id(3)
   dl=pdihs(indx)%member_id(4)

   ri=(/atoms(di)%x,atoms(di)%y,atoms(di)%z/)
   rj=(/atoms(dj)%x,atoms(dj)%y,atoms(dj)%z/)
   rk=(/atoms(dk)%x,atoms(dk)%y,atoms(dk)%z/)
   rl=(/atoms(dl)%x,atoms(dl)%y,atoms(dl)%z/)

   rji=rj-ri
   rjk=rj-rk
   rkj=rk-rj
   rkl=rk-rl

   m=vectorCrossProduct(rji,rjk)
   n=vectorCrossProduct(rkj,rkl)
   dummy=norm2(m)*norm2(n)

   cos_=dot_product(m,n)/dummy
   sin_=norm2(vectorCrossProduct(m,n))/dummy

   angle=acos(cos_)

   end function calculateDihedralAngle


!   function vectorCrossProduct(vect1,vect2) result(crossprod)
!   real(kind=real64), dimension(3), intent(in) :: vect1,vect2
!   real(kind=real64), dimension(3) :: crossprod

!   crossprod(1)=vect1(2)*vect2(3)-vect1(3)*vect2(2)
!   crossprod(2)=vect1(3)*vect2(1)-vect1(1)*vect2(3)
!   crossprod(3)=vect1(1)*vect2(2)-vect1(2)*vect2(1)

!   end function vectorCrossProduct


   subroutine calculateDihedralAnglesEnergies(atoms,pdihs)
   type(atom_t), dimension(:), intent(in) :: atoms
   type(pdih_t), dimension(:), intent(inout) :: pdihs
   integer(kind=int64) :: i,j
   real(kind=real64) :: energy

   do i=1,size(pdihs,1)
      pdihs(i)%currentAngle=calculateDihedralAngle(i,atoms,pdihs)
      energy=0.0_real64
      do j=1,size(pdihs(i)%params,1)
         energy=energy+pdihs(i)%params(j)%kd*(1.0_real64+&
                cos(pdihs(i)%params(j)%mult*pdihs(i)%currentAngle-&
                    pdihs(i)%params(j)%phi0))
      end do
      pdihs(i)%energy=energy
   end do

   end subroutine calculateDihedralAnglesEnergies


   subroutine putAtomsInEpsGroups(atoms,eps_M)
   type(atom_t), dimension(:), intent(inout) :: atoms
   real(kind=real64), dimension(:,:), allocatable, intent(out) :: eps_M
   real(kind=real64), dimension(:), allocatable :: u_eps
   integer(kind=int64), dimension(1) :: atoms_s
   integer(kind=int64) :: i

   atoms_s=shape(atoms)

   allocate(u_eps(atoms_s(1)))

   do i=1,atoms_s(1)
      u_eps(i)=atoms(i)%eps
   end do

   call findUniqueElements1Darray(u_eps)

   do i=1,atoms_s(1)
      atoms(i)%eps_g=searchIn1DFloatArray(u_eps,atoms(i)%eps)
   end do

   call constructEpsMatrix(u_eps,eps_M)

   end subroutine putAtomsInEpsGroups


   subroutine putAtomsInSigmaGroups(atoms,sigma_M)
   type(atom_t), dimension(:), intent(inout) :: atoms
   real(kind=real64), dimension(:,:), allocatable, intent(out) :: &
   sigma_M
   real(kind=real64), dimension(:), allocatable :: u_sgm
   integer(kind=int64), dimension(1) :: atoms_s
   integer(kind=int64) :: i

   atoms_s=shape(atoms)

   allocate(u_sgm(atoms_s(1)))

   do i=1,atoms_s(1)
      u_sgm(i)=atoms(i)%sigma
   end do

   call findUniqueElements1Darray(u_sgm)

   do i=1,atoms_s(1)
      atoms(i)%sigma_g=searchIn1DFloatArray(u_sgm,atoms(i)%sigma)
   end do

   call constructSigmaMatrix(u_sgm,sigma_M)

   end subroutine putAtomsInSigmaGroups


   subroutine findUniqueElements1Darray(array)
   real(kind=real64), dimension(:), allocatable, intent(inout) :: array
   real(kind=real64), dimension(:), allocatable :: tmp
   integer(kind=int64),dimension(1) :: array_s
   integer(kind=int64) :: i,counter,dummy_i
   real(kind=real64) :: dummy

   array_s=shape(array)
   array_s=array_s-1

   dummy_i=1
   do while (dummy_i .gt. 0)
      dummy_i=0
      do i=1,array_s(1)
         if (array(i)>array(i+1)) then
            dummy=array(i+1)
            array(i+1)=array(i)
            array(i)=dummy
            dummy_i=dummy_i+1
         end if
      end do
   end do

   allocate(tmp(1))
   tmp(1)=array(1)

   do i=1,array_s(1)
      if (array(i+1) .ne. array(i)) then
         call extend1DFloatArray(tmp,array(i+1))
      end if

   end do

   deallocate(array)

   call move_alloc(tmp,array)

   end subroutine findUniqueElements1Darray


   subroutine extend1DFloatArray(array,newElement)
   real(kind=real64), dimension(:), allocatable, intent(inout) :: array
   real(kind=real64), intent(in) :: newElement
   real(kind=real64), dimension(:), allocatable :: temp
   integer(kind=int64), dimension(1) :: array_s
   integer(kind=int64) :: i

   array_s=shape(array)

   call move_alloc(array,temp)

   array_s(1)=array_s(1)+1

   allocate(array(array_s(1)))

   array(1:array_s(1)-1)=temp(:)

   deallocate(temp)

   array(array_s(1))=newElement

   end subroutine extend1DFloatArray


   function searchIn1DFloatArray(array,numToSearchFor) result(found)
   real(kind=real64), dimension(:), intent(in) :: array
   real(kind=real64), intent(in) :: numToSearchFor
   logical :: flag
   integer(kind=int64) :: found,dummy
   integer(kind=int64), dimension(1) :: array_s

   array_s=shape(array)
   array_s(1)=array_s(1)+1

   flag=.True.

   dummy=1

   do while (flag)
      if (array(dummy) .eq. numToSearchFor) then
         flag=.False.
         found=dummy
      end if
      dummy=dummy+1
      if (dummy .eq. array_s(1)) then
         flag=.False.
      end if
   end do

   end function searchIn1DFloatArray


   subroutine constructEpsMatrix(u_eps,eps_M)
   real(kind=real64), dimension(:), intent(in) :: u_eps
   real(kind=real64), dimension(:,:), allocatable, intent(out) :: eps_M
   integer(kind=int64), dimension(1) :: u_eps_s
   integer(kind=int64) :: i,j

   u_eps_s=shape(u_eps)

   allocate(eps_M(u_eps_s(1),u_eps_s(1)))

   ! Remarks: Do multiply the epsilon matrix elements by 4.
   ! It is needed by the formula used for computing VdW potential:
   ! V(i,j)=4.0*epsilon(i,j)*((sigma(i,j)/r(i,j))**12-(sigma(i,j)/r(i,j))**6)

   do i=1,u_eps_s(1)
      eps_M(i,i)=4.0_real64*u_eps(i)
      do j=i+1,u_eps_s(1)
         eps_M(i,j)=4.0_real64*sqrt(u_eps(i)*u_eps(j))
         eps_M(j,i)=eps_M(i,j)
      end do
   end do

   end subroutine constructEpsMatrix


   subroutine constructSigmaMatrix(u_sgm,sigma_M)
   real(kind=real64), dimension(:), intent(in) :: u_sgm
   real(kind=real64), dimension(:,:), allocatable, intent(out) :: &
   sigma_M
   integer(kind=int64), dimension(1) :: u_sgm_s
   integer(kind=int64) :: i,j

   u_sgm_s=shape(u_sgm)

   allocate(sigma_M(u_sgm_s(1),u_sgm_s(1)))

   do i=1,u_sgm_s(1)
      sigma_M(i,i)=u_sgm(i)
      do j=i+1,u_sgm_s(1)
         sigma_M(i,j)=0.5_real64*(u_sgm(i)+u_sgm(j))
         sigma_M(j,i)=sigma_M(i,j)
      end do
   end do

   end subroutine constructSigmaMatrix


end module mod_input
