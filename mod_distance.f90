module mod_distance

   use mod_t

   implicit none
   !
   ! Note that atoms_s (the number of the atoms in the system) is defined in
   ! mod_t.f90 and its value is assigned in mod_input.f90 (see subroutine
   ! readInput there).
   !


contains

   subroutine updateDistanceMatrix(atoms,distM,atomIDs)
   !
   ! If the position of one atom got changed one should update the row in the
   ! distance matrix that corresonds to the neighbours of the atom. Because here
   ! the diagonalization matrix is represented as 1D array, (see the subroutine
   ! createDistanceMatrix) the corresponding array elements should be found and
   ! changed.
   !
   type(atom_t), dimension(:), intent(in) :: atoms
   real(kind=real64), dimension(:), intent(inout) :: distM
   integer(kind=int64), dimension(:), intent(in) :: atomIDs
   integer(kind=int64) :: i

   do i=1,size(atomIDs,1)
      call updateDistanceMatrix1(atoms,distM,atomIDs(i))
   end do

   end subroutine updateDistanceMatrix


   subroutine updateDistanceMatrix1(atoms,distM,atomID)
   !
   ! At each step of the execution subroutine updateDistanceMatrix calls this
   ! subroutine to change the distance matrix 1D-array representation.
   !
   type(atom_t), dimension(:), intent(in) :: atoms
   real(kind=real64), dimension(:), intent(inout) :: distM
   integer(kind=int64), intent(in) :: atomID
   integer(kind=int64) :: i

   do i=atomID+1,atoms_s
      distM(get1DdimMIndex(atomID,i))=norm2((/atoms(i)%x,atoms(i)%y,atoms(i)%z/)-&
      (/atoms(atomID)%x,atoms(atomID)%y,atoms(atomID)%z/))
   end do

   end subroutine updateDistanceMatrix1


   function get1DdimMIndex(i,j) result(indx)
   integer(kind=int64), intent(in) :: i,j
   integer(kind=int64) :: indx

   indx=atoms_s*(i-1)-i*(i+1)/2+j

   end function get1DdimMIndex


end module mod_distance
