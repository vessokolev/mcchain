module mod_select

   use mod_t
   use mod_array

   implicit none

contains

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
