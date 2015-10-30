module mod_tmp

   use mod_t
   use mod_array

   implicit none


contains

   subroutine getAtomIDsToRotate(bonds,pdihs,dihID,atomsIDsToRotate)
   type(bond_t), dimension(:), intent(in) :: bonds
   type(pdih_t), dimension(:), intent(in) :: pdihs
   integer(kind=int64), intent(in) :: dihID
   integer(kind=int64), dimension(:), allocatable, intent(out) :: atomsIDsToRotate
   integer(kind=int64), dimension(:,:), allocatable :: selection
   integer(kind=int64), dimension(2) :: selection_s
   integer(kind=int64), dimension(:), allocatable :: excludedAtoms
   logical :: flag

   call initSelection(bonds,pdihs(dihID)%member_id(2),pdihs(dihID)%member_id(3),selection)

   selection_s=shape(selection)

   if (all(selection_s(:) .ne. 0)) then

      flag=.True.

      do while (flag)
         call updateSelection(bonds,selection,pdihs(dihID)%member_id(2))
         selection_s=shape(selection)
         flag=any(selection(:,selection_s(2)) .gt. 0)
      end do

      allocate(excludedAtoms(5))

      excludedAtoms(1)=-2
      excludedAtoms(2)=-1
      excludedAtoms(3)=0
      excludedAtoms(4)=pdihs(dihID)%member_id(2)
      excludedAtoms(5)=pdihs(dihID)%member_id(3)

      call getUniqueElements(selection,excludedAtoms,atomsIDsToRotate)

   else
      allocate(atomsIDsToRotate(1))
      atomsIDsToRotate(1)=0
   end if

   end subroutine getAtomIDsToRotate


   subroutine initSelection(bonds,atomID,exlcAtomID,selection)
   type(bond_t), dimension(:), intent(in) :: bonds
   integer(kind=int64), intent(in) :: atomID, exlcAtomID
   integer(kind=int64), dimension(:,:), allocatable, intent(out) :: selection
   integer(kind=int64), dimension(:), allocatable :: neighs
   integer(kind=int64), dimension(1) :: neighs_s
   integer(kind=int64) :: i

   call getNeighbours(bonds,atomID,exlcAtomID,neighs)

   neighs_s=shape(neighs)

   allocate(selection(neighs_s(1),2))

   do i=1,neighs_s(1)
      selection(i,1)=atomID
      selection(i,2)=neighs(i)
   end do

   deallocate(neighs)

   end subroutine initSelection


   subroutine updateSelection(bonds,selection,exclAtomID)
   type(bond_t), dimension(:), intent(in) :: bonds
   integer(kind=int64), dimension(:,:), allocatable, intent(inout) :: selection
   integer(kind=int64), intent(in) :: exclAtomID
   integer(kind=int64), dimension(1) :: neighs_s
   integer(kind=int64), dimension(2) :: selection_s,subsel_s
   integer(kind=int64), dimension(2) :: newSelection_s
   integer(kind=int64) :: i,j,k
   integer(kind=int64) :: counter,counterNewSel
   integer(kind=int64) :: actElem
   integer(kind=int64), dimension(:,:), allocatable :: subsel,newSelection
   integer(kind=int64), dimension(:), allocatable :: neighs
   logical :: flag,flagNewSel

   selection_s=shape(selection)

   flag=.True.
   counter=1

   do i=1,selection_s(1)
      if (selection(i,selection_s(2)) .gt. 0) then
         if (flag) then
            allocate(subsel(counter,2))
            subsel(counter,1)=selection(i,selection_s(2))
            subsel(counter,2)=i
            flag=.False.
            counter=counter+1
         else
            call extendArrayAddRow(subsel)
            subsel(counter,1)=selection(i,selection_s(2))
            subsel(counter,2)=i
            counter=counter+1
         end if
      end if
   end do

   subsel_s=shape(subsel)

   if (.not. flag) then

      subsel_s=shape(subsel)
      flagNewSel=.True.
      counter=1
      counterNewSel=1

      allocate(newSelection(1,selection_s(2)+1))

      newSelection_s=shape(newSelection)

      do i=1,selection_s(1)
         if (all(subsel(:,2) .ne. i)) then
            if (flagNewSel) then
               newSelection(counter,1:selection_s(2))=selection(i,:)
               newSelection(counter,selection_s(2)+1)=0
               counter=counter+1
               flagNewSel=.False.
            else
               call extendArrayAddRow(newSelection)
               newSelection(counter,1:selection_s(2))=selection(i,:)
               newSelection(counter,selection_s(2)+1)=0
               counter=counter+1
            end if
         else
            call getNeighbours(bonds,subsel(counterNewSel,1),exclAtomID,neighs)
            neighs_s=shape(neighs)
            counterNewSel=counterNewSel+1
            do j=1,neighs_s(1)
               if (neighs(j) .ne. selection(i,selection_s(2)-1)) then
                  if (flagNewSel) then
                     flagNewSel=.False.
                  else
                     call extendArrayAddRow(newSelection)
                  end if
                  newSelection(counter,1:selection_s(2))=selection(i,:)
                  if (neighs(j) .eq. exclAtomID) then
                     newSelection(counter,selection_s(2)+1)=-2
                  else
                     if (neighs(j) .eq. selection(i,1)) then
                        newSelection(counter,selection_s(2)+1)=-1
                     else
                        newSelection(counter,selection_s(2)+1)=neighs(j)
                     end if
                  end if
                  counter=counter+1
               else
                  if (neighs_s(1) .eq. 1) then
                     if (flagNewSel) then
                        flagNewSel=.False.
                     else
                        call extendArrayAddRow(newSelection)
                     end if
                     newSelection(counter,1:selection_s(2))=selection(i,:)
                     newSelection(counter,selection_s(2)+1)=0
                     counter=counter+1
                  end if
               end if
            end do
         end if
      end do
   end if

   newSelection_s=shape(newSelection)

   deallocate(selection)

   call move_alloc(newSelection,selection)

   deallocate(subsel)
   deallocate(neighs)

   end subroutine updateSelection


   subroutine extend1DArray(array,newElement)
   integer(kind=int64), dimension(:), allocatable, intent(inout) :: array
   integer(kind=int64), intent(in) :: newElement
   integer(kind=int64), dimension(:), allocatable :: temp
   integer(kind=int64), dimension(1) :: array_s
   integer(kind=int64) :: i

   array_s=shape(array)

   call move_alloc(array,temp)

   array_s(1)=array_s(1)+1

   allocate(array(array_s(1)))

   array(1:array_s(1)-1)=temp(:)

   deallocate(temp)

   array(array_s(1))=newElement

   end subroutine extend1DArray


   subroutine getNeighbours(bonds,inclAtomID,exclAtomID,selection)
   type(bond_t), dimension(:), intent(in) :: bonds
   integer(kind=int64), intent(in) :: inclAtomID,exclAtomID
   integer(kind=int64), dimension(:), allocatable, intent(out) :: selection
   integer(kind=int64) :: i
   logical :: flag

   flag=.True.

   allocate(selection(1))

   selection(1)=0

   do i=1,size(bonds,1)
      if ((bonds(i)%left_a .ne. exclAtomID) .and. (bonds(i)%right_a .ne. &
           exclAtomID)) then
         if (bonds(i)%left_a .eq. inclAtomID) then
            if (flag) then
               selection(1)=bonds(i)%right_a
               flag=.False.
            else
               call extend1DIntArray(selection,bonds(i)%right_a)
            end if
         end if
         if (bonds(i)%right_a .eq. inclAtomID) then
            if (flag) then
               selection(1)=bonds(i)%left_a
               flag=.False.
            else
               call extend1DIntArray(selection,bonds(i)%left_a)
            end if
         end if
      end if
   end do

   end subroutine getNeighbours


   subroutine getUniqueElements(selection,excludedAtoms,uniqueArray)
   integer(kind=int64), dimension(:,:), intent(in) :: selection
   integer(kind=int64), dimension(:), intent(in) :: excludedAtoms
   integer(kind=int64), dimension(:), allocatable, intent(out) :: uniqueArray
   integer(kind=int64), dimension(2) :: selection_s
   integer(kind=int64) :: i,j
   logical :: flag

   selection_s=shape(selection)

   allocate(uniqueArray(1))
   uniqueArray(1)=0

   flag=.True.

   do i=1,selection_s(1)
      do j=1,selection_s(2)
         if (all(uniqueArray(:) .ne. selection(i,j)) .and. &
            (all(excludedAtoms(:) .ne. selection(i,j)))) then
            if (flag) then
               uniqueArray(1)=selection(i,j)
               flag=.False.
            else
               call extend1DArray(uniqueArray,selection(i,j))
            end if
         end if
      end do
   end do

   end subroutine getUniqueElements


!   subroutine extendArrayAddRow(array)
!   integer(kind=int64), dimension(:,:), allocatable, intent(inout) :: array
!   integer(kind=int64), dimension(:,:), allocatable :: temp
!   integer(kind=int64), dimension(2) :: array_s
!   integer(kind=int64) :: i

!   array_s=shape(array)
!   call move_alloc(array,temp)
!   array_s(1)=array_s(1)+1
!   allocate(array(array_s(1),array_s(2)))

!   do i=1,array_s(1)-1
!      array(i,:)=temp(i,:)
!   end do

!   array(array_s(1),:)=0
!   deallocate(temp)

!   end subroutine

end module mod_tmp
