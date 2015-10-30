module mod_distance
   ! Created by Veselin Kolev <vesso.kolev@gmail.com>
   ! 20151020031723
   use mod_t
   use mod_array

   implicit none
   !
   ! Note that atoms_s (the number of the atoms in the system) is defined as a
   ! global constant in the body of module mod_t.f90 and its value is being 
   ! assigned in mod_input.f90 (see subroutine readInput there).
   !


contains

   subroutine updateDistanceMatrix(atoms,distM,atomIDs)
   !
   ! This subroutine calls updateDistanceMatrix1Atom for each of the atom IDs
   ! supplied by the array atomIDs. Its role is to update the distance matrix at
   ! complex atomic coordinate changes (like rotation over axis).
   !
   type(atom_t), dimension(:), intent(in) :: atoms
   real(kind=real32), dimension(:), intent(inout) :: distM
   integer(kind=int64), dimension(:), intent(in) :: atomIDs
   integer(kind=int64) :: i

   do i=1,size(atomIDs,1)
      call updateDistanceMatrix1Atom(atoms,distM,atomIDs(i))
   end do

   end subroutine updateDistanceMatrix


   subroutine updateDistanceMatrix1Atom(atoms,distM1D,atomID)
   !
   ! Updates the distance matrix with respect to the changes of the coordinates
   ! of the atom with atomID. This is an operation that reqiures changing of
   ! some of the elements of the distance matrix. Note that here we use 1D
   ! representation of the distance matrix.
   !
   type(atom_t), dimension(:), intent(in) :: atoms
   real(kind=real32), dimension(:), intent(inout) :: distM1D
   integer(kind=int64), intent(in) :: atomID
   integer(kind=int64) :: j,k

   do j=1,atomID
      if (j .eq. atomID) then
         do k=atomID+1,atoms_s
            distM1D(get1DdimMIndex(atomID,k))=&
            norm02((/atoms(atomID)%x,atoms(atomID)%y,atoms(atomID)%z/)-&
                  (/atoms(k)%x,atoms(k)%y,atoms(k)%z/))
         end do
      else
         distM1D(get1DdimMIndex(j,atomID))=&
         norm02((/atoms(j)%x,atoms(j)%y,atoms(j)%z/)-&
               (/atoms(atomID)%x,atoms(atomID)%y,atoms(atomID)%z/))
      end if
   end do

   end subroutine updateDistanceMatrix1Atom


   function get1DdimMIndex(i,j) result(indx)
   !
   ! This function supplies the index of a one-dimensional array that
   ! corresponds to the element of a square matrix located above its diagonal.
   !
   integer(kind=int64), intent(in) :: i,j
   integer(kind=int64) :: indx

   indx=atoms_s*(i-1)-i*(i+1)/2+j

   end function get1DdimMIndex


   subroutine create1DdistM(atoms,distM1D)
   !
   ! This subroutine creates a distance matrix as an 1-dimensional array to save
   ! a lot of memory. Note that instead of calling function get1DdimMIndex(i,j)
   ! it is possible only during the creation of the matrix and assigning values
   ! to the elements, to use increment of integer variable (see how the variable
   ! counter is used bellow).
   !
   type(atom_t), dimension(:), intent(in) :: atoms
   real(kind=real32), dimension(:), allocatable, intent(inout) :: distM1D
   integer(kind=int64) :: i,j,counter

   allocate(distM1D(atoms_s*(atoms_s-1)/2))

   counter=1

   do i=1,atoms_s-1
      do j=i+1,atoms_s
         distM1D(counter)=&
         norm02((/atoms(i)%x,atoms(i)%y,atoms(i)%z/)-&
               (/atoms(j)%x,atoms(j)%y,atoms(j)%z/))
         counter=counter+1
      end do
   end do

   end subroutine create1DDistM

end module mod_distance
