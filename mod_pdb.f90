module mod_pdb

   use mod_t

   character(len=26) :: cryst1fmt="(A6,3F9.3,3F7.2,1X,A11,I4)"
   character(len=56) :: atomfmt="(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,6X,A4,2A2)"

contains

   subroutine saveAtomsOnPDB(atoms,inputfile,outputfile)
   type(atom_t), dimension(:), intent(in) :: atoms
   character(len=*), intent(in) :: inputfile
   character(len=*), intent(in) :: outputfile
   type(pdb_cryst1_t) :: pdbcryst1
   type(pdb_atom_t), dimension(:), allocatable :: pdbatoms
   integer(kind=int64) :: i

   call pdbread(inputfile,pdbcryst1,pdbatoms)

   do i=1,size(atoms,1)
      pdbatoms(i)%x=atoms(i)%x
      pdbatoms(i)%y=atoms(i)%y
      pdbatoms(i)%z=atoms(i)%z
   end do

   pdbatoms(:)%x=pdbatoms(:)%x*10.0_real64
   pdbatoms(:)%y=pdbatoms(:)%y*10.0_real64
   pdbatoms(:)%z=pdbatoms(:)%z*10.0_real64

   call pdbwrite(outputfile,pdbcryst1,pdbatoms)

   end subroutine saveAtomsOnPDB


   subroutine pdbread(inputfile,pdbcryst1,pdbatoms)
   character(len=*), intent(in) :: inputfile
   type(pdb_cryst1_t), intent(out) :: pdbcryst1
   type(pdb_atom_t), dimension(:), allocatable, intent(out) :: pdbatoms
   logical :: flag
   integer(kind=int64) :: numatoms
   type(pdb_atom_t) :: tmp

   open(unit=1,file=inputfile,form="FORMATTED")

   allocate(pdbatoms(1))

   flag=.True.
   numatoms=0

   do while (flag)
      read(1,fmt="(A6)") tmp%ident

      if (tmp%ident .eq. "CRYST1") then
         backspace(1)
         read(1,fmt=cryst1fmt) pdbcryst1%rname,pdbcryst1%a,pdbcryst1%b,&
         pdbcryst1%c,pdbcryst1%alpha,pdbcryst1%beta,pdbcryst1%gamma,&
         pdbcryst1%sgroup,pdbcryst1%zval
      end if

      if ((tmp%ident .eq. "ATOM") .or. (tmp%ident .eq. "HETATM")) then
         backspace(1)
         read(1,fmt=atomfmt) &
         tmp%ident,tmp%anum,tmp%aname,tmp%altloc,tmp%resname,tmp%chainid,&
         tmp%resnum,tmp%inscode,tmp%x,tmp%y,tmp%z,tmp%occup,tmp%tempf,&
         tmp%segmentid,tmp%elsymb,tmp%charge
         numatoms=numatoms+1
         if (numatoms .eq. 1) then
            pdbatoms(1)=tmp
         else
            call extendpdbatoms(pdbatoms,tmp)
         end if
      else
         if ((tmp%ident .eq. "END") .or. (tmp%ident .eq. "ENDMDL")) then
            flag=.False.
         end if
      end if
   end do

   close(1)

   contains

      subroutine extendpdbatoms(pdbatoms,newElement)
      type(pdb_atom_t) :: newElement
      type(pdb_atom_t), dimension(:), allocatable, intent(inout) :: pdbatoms
      integer(kind=int64), dimension(1) :: array_s
      type(pdb_atom_t), dimension(:), allocatable :: temp

      array_s=shape(pdbatoms)

      call move_alloc(pdbatoms,temp)

      array_s(1)=array_s(1)+1

      allocate(pdbatoms(array_s(1)))

      pdbatoms(1:array_s(1)-1)=temp(:)

      deallocate(temp)

      pdbatoms(array_s(1))=newElement

      end subroutine extendpdbatoms

   end subroutine pdbread


   subroutine pdbwrite(outputfile,pdbcryst1,pdbatoms)
   character(len=*), intent(in) :: outputfile
   type(pdb_cryst1_t), intent(in) :: pdbcryst1
   type(pdb_atom_t), dimension(:), intent(in) :: pdbatoms
   logical :: flag
   integer(kind=int64) :: i
   type(pdb_atom_t) :: tmp

   open(unit=1,file=outputfile,form="FORMATTED")

   write(1,fmt=cryst1fmt) pdbcryst1

   do i=1,size(pdbatoms,1)
      tmp=pdbatoms(i)
      write(1,fmt=atomfmt) &
      tmp%ident,tmp%anum,tmp%aname,tmp%altloc,tmp%resname,tmp%chainid,&
      tmp%resnum,tmp%inscode,tmp%x,tmp%y,tmp%z,tmp%occup,tmp%tempf,&
      tmp%segmentid,adjustr(tmp%elsymb),tmp%charge
   end do

   write(1,fmt="(A3)") "END"

   close(1)

   end subroutine pdbwrite

end module mod_pdb
