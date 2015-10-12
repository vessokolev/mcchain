module mod_select

   ! Created by Veselin Kolev <vesso.kolev@gmail.com>
   ! 20151011034812

   use mod_t
   use mod_array

   implicit none

contains


   subroutine createAtomsToRotateStruct(atoms,bonds,pdihs,dihsOfInter,struct,&
                                        size_s)
   !
   ! This subroutine creates the matrix struct and the vector size_s.
   ! Each line of the matrix contains the atomIDs that should be rotated
   ! with respect to the axis of dihsOfInter(i) proper dihedral.
   ! Because the matrix may contain zero-elements due to the different sizes
   ! of the set of rotated atoms, the size_s vector contains an information
   ! how many of the matrix elements in the respective row are non-zero
   ! sequence. I.e. size_s=struct(i,1:number_of_atomIDs_to_rotate).
   !
   type(atom_t), dimension(:), intent(in) :: atoms
   type(bond_t), dimension(:), intent(in) :: bonds
   type(pdih_t), dimension(:), intent(in) :: pdihs
   integer(kind=int64), dimension(:), intent(in) :: dihsOfInter
   integer(kind=int64), dimension(:,:), allocatable, intent(out) :: struct
   integer(kind=int64), dimension(:), allocatable, intent(out) :: size_s
   integer(kind=int64), dimension(:), allocatable :: atomIDsToRotate
   integer(kind=int64) :: i,j,dihsOfInter_s,atomIDsToRotate_s
   integer(kind=int64), dimension(2) :: struct_s
   logical :: flag

   dihsOfInter_s=size(dihsOfInter,1)

   flag=.True.

   do i=1,dihsOfInter_s
      call optimizeAtomIDsToRotateLen(atoms,bonds,pdihs,dihsOfInter(i),atomIDsToRotate)
      atomIDsToRotate_s=size(atomIDsToRotate,1)
      if (flag) then
         allocate(struct(1,atomIDsToRotate_s))
         allocate(size_s(1))
         struct(1,:)=atomIDsToRotate(:)
         size_s(1)=atomIDsToRotate_s
         flag=.False.
      else
         struct_s=shape(struct)
         if (atomIDsToRotate_s .gt. struct_s(2)) then
            do j=1,atomIDsToRotate_s-struct_s(2)
               call extendArrayAddColumn(struct)
            end do
         end if
         call extendArrayAddRow(struct)
         struct(i,1:atomIDsToRotate_s)=atomIDsToRotate(:)
         call extend1DIntArray(size_s,atomIDsToRotate_s)
      end if
   end do

   end subroutine createAtomsToRotateStruct


   subroutine optimizeAtomIDsToRotateLen(atoms,bonds,pdihs,pdihID,atomIDsToRotate)
   !
   ! This is a helper routine. It helps createAtomsToRotateStruct to minimize
   ! the number of atom which coordinates are rotated.
   !
   type(atom_t), dimension(:), intent(in) :: atoms
   type(bond_t), dimension(:), intent(in) :: bonds
   type(pdih_t), dimension(:), intent(in) :: pdihs
   integer(kind=int64), intent(in) :: pdihID
   integer(kind=int64), dimension(:), allocatable, intent(out) :: atomIDsToRotate
   integer(kind=int64), dimension(:), allocatable :: tmp
   integer(kind=int64) :: i,atoms_s
   logical :: flag
   
   atoms_s=size(atoms,1)

   call getAtomsToRotate(bonds,pdihs,pdihID,atomIDsToRotate)

   if (size(atomIDsToRotate,1) .gt. atoms_s/2) then
      allocate(tmp(1))
      flag=.True.
      do i=1,atoms_s
         if (.not. checkIfElelemntIsInArray(atomIDsToRotate,i)) then
            if (flag) then
               tmp(1)=i
               flag=.False.
            else
               call extend1DIntArray(tmp,i)
            end if
         end if
      end do
      deallocate(atomIDsToRotate)
      allocate(atomIDsToRotate(size(tmp,1)))
      call move_alloc(tmp,atomIDsToRotate)
   end if

   end subroutine optimizeAtomIDsToRotateLen


   subroutine getAtomsToRotate(bonds,pdihs,dihID,atomIDsToRotate)
   !
   ! It is an enhancement of getAtomIDsToRotate subroutine.
   !
   type(bond_t), dimension(:), intent(in) :: bonds
   type(pdih_t), dimension(:), intent(in) :: pdihs
   integer(kind=int64) :: dihID
   integer(kind=int64), dimension(:), allocatable, intent(out) :: &
   atomIDsToRotate
   integer(kind=int64), dimension(:), allocatable :: excludes

   allocate(excludes(1))

   excludes(1)=pdihs(dihID)%member_id(2)

   call getAtomIDsToRotate(bonds,pdihs(dihID)%member_id(3),&
   excludes,atomIDsToRotate)

   deallocate(excludes)

   end subroutine getAtomsToRotate


   subroutine getAtomIDsToRotate(bonds,atomID,excludes,atomIDsToRotate)
   type(bond_t), dimension(:), intent(in) :: bonds
   integer(kind=int64), intent(in) :: atomID
   integer(kind=int64), dimension(:), allocatable, intent(inout) :: excludes
   integer(kind=int64), dimension(:), allocatable, intent(out) :: &
   atomIDsToRotate
   integer(kind=int64), dimension(:,:), allocatable :: atomsInPath

   call getTreePaths(bonds,atomID,excludes,atomsInPath)
   call unique1DArray(atomsInPath,atomIDsToRotate)

   end subroutine getAtomIDsToRotate


   subroutine getTreePaths(bonds,atomID,excludes,atomsInPath)
   type(bond_t), dimension(:), intent(in) :: bonds
   integer(kind=int64), intent(in) :: atomID
   integer(kind=int64), dimension(:), allocatable, intent(inout) :: excludes
   integer(kind=int64), dimension(:,:), allocatable, intent(out) :: atomsInPath
   integer(kind=int64), dimension(:,:), allocatable :: tmp
   integer(kind=int64), dimension(:), allocatable :: neighs,excl
   integer(kind=int64) :: i,excl_s,excludes_s
   integer(kind=int8) :: j
   integer(kind=int64), dimension(2) :: atomsInPath_s,tmp_s
   logical :: flag,flag1

   allocate(atomsInPath(1,1))
   atomsInPath(1,1)=atomID

   call extend1DIntArray(excludes,atomID)

   excludes_s=size(excludes,1)

   atomsInPath_s=shape(atomsInPath)

   flag=.True.

   if (ALL(atomsInPath(:,atomsInPath_s(2)) .gt. 0)) then

      do while (flag)

         tmp_s=shape(atomsInPath)

         call move_alloc(atomsInPath,tmp)

         allocate(atomsInPath(1,tmp_s(2)+1))

         atomsInPath_s=shape(atomsInPath)

         flag1=.True.

         do i=1,tmp_s(1)

            excl_s=tmp_s(2)+excludes_s

            allocate(excl(excl_s))

            excl(1:tmp_s(2))=tmp(i,:)
            excl(tmp_s(2)+1:excl_s)=excludes(:)

            call getNeighs(bonds,tmp(i,tmp_s(2)),excl,neighs)

            if (sum(neighs) .gt. 0) then

               do j=1,size(neighs,1)
                  if (flag1) then
                     atomsInPath(atomsInPath_s(1),1:atomsInPath_s(2)-1)=tmp(i,:)
                     atomsInPath(atomsInPath_s(1),atomsInPath_s(2))=neighs(j)
                     flag1=.False.
                  else
                     call extendArrayAddRow(atomsInPath)
                     atomsInPath_s=shape(atomsInPath)
                     atomsInPath(atomsInPath_s(1),1:atomsInPath_s(2)-1)=tmp(i,:)
                     atomsInPath(atomsInPath_s(1),atomsInPath_s(2))=neighs(j)
                  end if
               end do
            else
               if (flag1) then
                  atomsInPath(atomsInPath_s(1),1:atomsInPath_s(2)-1)=tmp(i,:)
                  atomsInPath(atomsInPath_s(1),atomsInPath_s(2))=0
                  flag1=.False.
               else
                  call extendArrayAddRow(atomsInPath)
                  atomsInPath_s=shape(atomsInPath)
                  atomsInPath(atomsInPath_s(1),1:atomsInPath_s(2)-1)=tmp(i,:)
                  atomsInPath(atomsInPath_s(1),atomsInPath_s(2))=0
               end if
            end if
            deallocate(excl)
         end do
         flag=(sum(atomsInPath(:,atomsInPath_s(2))) .ne. 0)
      end do
   end if

   end subroutine getTreePaths


   subroutine getNeighs(bonds,atomID,excludes,neighs)
   type(bond_t), dimension(:), intent(in) :: bonds
   integer(kind=int64), intent(in) :: atomID
   integer(kind=int64), dimension(:), intent(in) :: excludes
   integer(kind=int64), dimension(:), allocatable, intent(out) :: neighs
   integer(kind=int64) :: i
   logical :: flag

   allocate(neighs(1))

   neighs(1)=0

   flag=.True.

   do i=1,size(bonds,1)
       if ((bonds(i)%left_a .eq. atomID)) then
          if (.not. ANY(excludes(:) .eq. bonds(i)%right_a)) then
             if (flag) then
                neighs(1)=bonds(i)%right_a
                flag=.False.
             else
                call extend1DIntArray(neighs,bonds(i)%right_a)
             end if
          end if
       end if

       if ((bonds(i)%right_a .eq. atomID)) then
          if (.not. ANY(excludes(:) .eq. bonds(i)%left_a)) then
             if (flag) then
                neighs(1)=bonds(i)%left_a
                flag=.False.
             else
                call extend1DIntArray(neighs,bonds(i)%left_a)
             end if
          end if
       end if
   end do

   end subroutine getNeighs



end module mod_select
