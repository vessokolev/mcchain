module mod_input
   !
   ! Create by Veselin Kolev <vesso.kolev@gmail.com>
   ! 20151008233108
   !
   use mod_t
   use mod_array
   use mod_dihedral

contains


   subroutine readInput(files,atoms,bonds,pdihs,eps_M,sigma_M)
   type(file_t), intent(in) :: files
   type(atom_t), dimension(:), allocatable, intent(out) :: atoms
   type(bond_t), dimension(:), allocatable, intent(out) :: bonds
   type(pdih_t), dimension(:), allocatable, intent(out) :: pdihs
   real(kind=real64), dimension(:,:), allocatable :: eps_M
   real(kind=real64), dimension(:,:), allocatable :: sigma_M

   call readAtoms(files,atoms)

   atoms_s=size(atoms,1)

   call readBonds(files,bonds)

   bonds_s=size(bonds,1)

   call readPropDihedrals(files,pdihs)

   pdihs_s=size(pdihs,1)

   call updateAtomAddBondsArray(atoms,bonds)

   call calculateDihedralAngles(atoms,pdihs)
   call calculateTorsionEnergies(pdihs)

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


   subroutine createDistanceMatrix(atoms,distM)
   !
   ! VERY IMPORTANT! To save memory we don't use two dimensional distance
   ! matrix. Instead an one-dimensional representation of the upper triangular
   ! matrix (excluding its diagonal elements) is implemented. The connection
   ! between the index variable of the 1D-array (q) and the indexes of the
   ! matrix
   ! elements (i and j) is:
   !
   ! q=N*(i-1)-i*(i+1)/2+j
   !
   ! where N is the size of the distance matrix (the number of its rows, or
   ! columns, or diagonal elements) - equal to the total number of the atoms.
   !
   type(atom_t), dimension(:), intent(in) :: atoms
   real(kind=real64), dimension(:), allocatable, intent(out) :: distM
   integer(kind=int64) :: i,j,counter

   allocate(distM(atoms_s*(atoms_s-1)/2))

   counter=1

   do i=1,atoms_s-1
      do j=i+1,atoms_s
            distM(counter)=norm2((/atoms(i)%x,atoms(i)%y,atoms(i)%z/)-&
                                 (/atoms(j)%x,atoms(j)%y,atoms(j)%z/))
         counter=counter+1
      end do
   end do

   end subroutine createDistanceMatrix


end module mod_input
