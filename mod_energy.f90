module mod_energy
   !
   ! Created by Vesselin Kolev <vesso.kolev@gmail.com>
   ! 20151013015825
   !
   use mod_t
   use mod_array

   implicit none

contains

   subroutine selectAtomConnections(atoms,bonds,connections,numconns)
   !
   ! When calculating vdW and Coulomb potential we need not to take into account
   ! the bonded atoms. I.e. the non-bonded interactions between bonded atoms
   ! have to be skipped! The routine defines a matrix and stores the connections
   ! of each atom in a row. The index of the row corresponds to the number of
   ! the atom in the storage array as defined by mod_input.
   !
   type(atom_t), dimension(:), intent(in) :: atoms
   type(bond_t), dimension(:), intent(in) :: bonds
   integer(kind=int64), dimension(:,:), allocatable, intent(out) :: connections
   integer(kind=int64), dimension(:), allocatable, intent(out) :: numconns
   integer(kind=int64) :: i,j,numconns_s
   integer(kind=int64), dimension(2) :: connections_s
   logical :: flag

   allocate(connections(1,1))
   allocate(numconns(1))

   flag=.True.
   numconns_s=1
   connections_s=(/1,1/)

   do i=1,size(atoms,1)

      if (.not. flag) then
         call extendArrayAddRow(connections)
         call extend1DIntArray(numconns,0)
         numconns_s=numconns_s+1
         numconns(numconns_s)=0
         connections_s(1)=connections_s(1)+1
      end if

      do j=1,size(bonds,1)
         if (bonds(j)%left_a .eq. i) then
            call update(bonds(j)%right_a)
         end if
         if (bonds(j)%right_a .eq. i) then
            call update(bonds(j)%left_a)
         end if
      end do

   end do

   contains

      subroutine update(atomIDtoAdd)
      integer(kind=int64), intent(in) :: atomIDtoAdd

      if (flag) then
         connections(1,1)=atomIDtoAdd
         numconns(1)=1
         flag=.False.
      else
         numconns(i)=numconns(i)+1
         if (numconns(i) .gt. connections_s(2)) then
            connections_s(2)=connections_s(2)+1
            call extendArrayAddColumn(connections)
         end if
         connections(i,numconns(i))=atomIDtoAdd
      end if

      end subroutine update


   end subroutine selectAtomConnections


   subroutine getTotvdWContribToEnergy(atoms,connections,numconns,eps_M,sigma_M,energy)
   type(atom_t), dimension(:), intent(in) :: atoms
   integer(kind=int64), dimension(:,:), intent(in) :: connections
   integer(kind=int64), dimension(:), intent(in) :: numconns
   real(kind=real64), dimension(:,:), intent(in) :: eps_M
   real(kind=real64), dimension(:,:), intent(in) :: sigma_M
   real(kind=real64), intent(out) :: energy
   integer(kind=int64) :: i,j,atoms_s

   atoms_s=size(atoms,1)

   energy=0.0_real64

   do i=1,atoms_s-1
      do j=i+1,atoms_s
         energy=energy+getvdWPairEnergy(atoms,i,j,eps_M,sigma_M)
      end do
   end do

   end subroutine getTotvdWContribToEnergy


   function getvdWPairEnergy(atoms,i,j,eps_M,sigma_M) result(pairE)
   type(atom_t), dimension(:), intent(in) :: atoms
   integer(kind=int64), intent(in) :: i
   integer(kind=int64), intent(in) :: j
   real(kind=real64), dimension(:,:), intent(in) :: eps_M
   real(kind=real64), dimension(:,:), intent(in) :: sigma_M
   real(kind=real64) :: pairE

   pairE=norm2((/atoms(i)%x,atoms(i)%y,atoms(i)%z/)-&
               (/atoms(j)%x,atoms(j)%y,atoms(j)%z/))

   pairE=(sigma_M(atoms(i)%sigma_g,atoms(j)%sigma_g)/pairE)**6.0_real64
   pairE=eps_M(atoms(i)%sigma_g,atoms(j)%sigma_g)*pairE*(pairE-1)

   end function getvdWPairEnergy


end module mod_energy
