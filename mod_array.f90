module mod_array

   use mod_t

   implicit none

contains

   function checkIfElelemntIsInArray(array,element) result(check)
   integer(kind=int64), dimension(:), intent(in) :: array
   integer(kind=int64), intent(in) :: element
   integer(kind=int64) :: i,array_s
   logical :: check,flag

   check=.False.

   i=1
   array_s=size(array,1)

   check=.False.
   flag=.True.

   do while (flag)
      if (array(i) .eq. element) then
         check=.True.
         flag=.False.
      end if
      i=i+1
      if (i .gt. array_s) then
         flag=.False.
      end if
   end do

   end function checkIfElelemntIsInArray


   subroutine extend1DIntArray(array,newElement)
   integer(kind=int64), dimension(:), allocatable, intent(inout) :: array
   integer(kind=int64), intent(in) :: newElement
   integer(kind=int64), dimension(:), allocatable :: temp
   integer(kind=int64) :: array_s
   integer(kind=int64) :: i

   array_s=size(array,1)

   call move_alloc(array,temp)

   array_s=array_s+1

   allocate(array(array_s))

   array(1:array_s-1)=temp(:)

   deallocate(temp)

   array(array_s)=newElement

   end subroutine extend1DIntArray


   function vectorCrossProduct(vect1,vect2) result(crossprod)
   real(kind=real64), dimension(3), intent(in) :: vect1,vect2
   real(kind=real64), dimension(3) :: crossprod

   crossprod(1)=vect1(2)*vect2(3)-vect1(3)*vect2(2)
   crossprod(2)=vect1(3)*vect2(1)-vect1(1)*vect2(3)
   crossprod(3)=vect1(1)*vect2(2)-vect1(2)*vect2(1)

   end function vectorCrossProduct


   function searchIn1DIntArray(array,numToSearchFor) result(found)
   integer(kind=int64), dimension(:), intent(in) :: array
   integer(kind=int64), intent(in) :: numToSearchFor
   logical :: found

   found=ANY(array(:)-numToSearchFor .eq. 0)

   end function searchIn1DIntArray


   subroutine extendArrayAddRow(array)
   integer(kind=int64), dimension(:,:), allocatable, intent(inout) :: array
   integer(kind=int64), dimension(:,:), allocatable :: temp
   integer(kind=int64), dimension(2) :: array_s
   integer(kind=int64) :: i

   array_s=shape(array)
   call move_alloc(array,temp)
   array_s(1)=array_s(1)+1
   allocate(array(array_s(1),array_s(2)))

   do i=1,array_s(1)-1
      array(i,:)=temp(i,:)
   end do

   array(array_s(1),:)=0
   deallocate(temp)

   end subroutine extendArrayAddRow


   subroutine extendArrayAddColumn(array)
   integer(kind=int64), dimension(:,:), allocatable, intent(inout) :: array
   integer(kind=int64), dimension(:,:), allocatable :: temp
   integer(kind=int64), dimension(2) :: array_s
   integer(kind=int64) :: i,dummy

   array_s=shape(array)
   call move_alloc(array,temp)
   array_s(2)=array_s(2)+1
   allocate(array(array_s(1),array_s(2)))

   dummy=array_s(2)-1

   do i=1,array_s(1)
      array(i,1:dummy)=temp(i,:)
   end do

   array(:,array_s(2))=0
   deallocate(temp)

   end subroutine extendArrayAddColumn


   subroutine unique1DArray(array2D,array1D)
   ! Transforms array2D into array1D by removing the
   ! repetitions and 0-elements.
   integer(kind=int64), dimension(:,:), intent(in) :: array2D
   integer(kind=int64), dimension(:), allocatable, intent(out) :: array1D
   logical :: flag
   integer(kind=int64) :: i,j,counter

   allocate(array1D(1))
   flag=.True.
   counter=1

   do i=1,size(array2D,1)
      do j=1,size(array2D,2)
         if (array2D(i,j) .gt. 0) then
            if (flag) then
               array1D(counter)=array2D(i,j)
               counter=counter+1
               flag=.False.
            else
               if (.not. checkIfElelemntIsInArray(array1D,array2D(i,j))) then
                  call extend1DIntArray(array1D,array2D(i,j))
                  array1D(counter)=array2D(i,j)
                  counter=counter+1
               end if
            end if
         end if
      end do
   end do

   end subroutine unique1DArray


   function searchIn1DFloatArray(array,numToSearchFor) result(found)
   real(kind=real64), dimension(:), intent(in) :: array
   real(kind=real64), intent(in) :: numToSearchFor
   logical :: flag
   integer(kind=int64) :: found,dummy,array_s

   array_s=size(array,1)
   array_s=array_s+1

   flag=.True.

   dummy=1

   do while (flag)
      if (array(dummy) .eq. numToSearchFor) then
         flag=.False.
         found=dummy
      end if
      dummy=dummy+1
      if (dummy .eq. array_s) then
         flag=.False.
      end if
   end do

   end function searchIn1DFloatArray


   subroutine extend1DFloatArray(array,newElement)
   real(kind=real64), dimension(:), allocatable, intent(inout) :: array
   real(kind=real64), intent(in) :: newElement
   real(kind=real64), dimension(:), allocatable :: temp
   integer(kind=int64) :: i,array_s

   array_s=size(array,1)

   call move_alloc(array,temp)

   array_s=array_s+1

   allocate(array(array_s))

   array(1:array_s-1)=temp(:)

   deallocate(temp)

   array(array_s)=newElement

   end subroutine extend1DFloatArray


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


end module mod_array
