program dev_do_select

   use mod_t ! All types are defined here
   use mod_pdb ! Reading and writing PDB files
   use mod_rot ! Rotating atomic coordinates
   use mod_input ! Reading the input files (not PDB)
   use mod_select ! Selects the atoms to rotate
   use mod_array

   implicit none

   type(atom_t), dimension(:), allocatable :: atoms
   type(bond_t), dimension(:), allocatable :: bonds
   type(pdih_t), dimension(:), allocatable :: pdihs
   real(kind=real64), dimension(:,:), allocatable :: eps_M
   real(kind=real64), dimension(:,:), allocatable :: sigma_M

   real(kind=real64), dimension(3) :: temp,p1,p2,new
   real(kind=real64), dimension(3,3) :: rotMatrix
   real(kind=real64) :: phi,r(5)
   integer(kind=int64), dimension(:,:), allocatable :: bonds_f
   integer(kind=int64), dimension(:), allocatable :: atomIDsToRotate
   integer(kind=int64), dimension(:,:), allocatable :: selection
   integer(kind=int64), dimension(:), allocatable :: pdihssel
   integer(kind=int64), dimension(2) :: axis
   real(kind=real64), dimension(3) :: tmp
   integer(kind=int64) :: i,selectedDihedralID
   type(pdih_t) :: selectedDihedral
   type(file_t) :: files


   files%atoms="/home/vesso/Code/MC_dih/current/files/atoms.txt"
   files%bonds="/home/vesso/Code/MC_dih/current/files/bonds.txt"
   files%pdihsDefs="/home/vesso/Code/MC_dih/current/files/pdhdr.txt"
   files%pdihsParams="/home/vesso/Code/MC_dih/current/files/pdhdr-params.txt"
   files%pdbinp="/home/vesso/Code/MC_dih/current/files/opt.pdb"
   files%pdbout="/home/vesso/Code/MC_dih/current/files/tmp.pdb"

   call readInput(files,atoms,bonds,pdihs,eps_M,sigma_M)

   allocate(pdihssel(10))

   pdihssel(1)=1
   pdihssel(2)=134
   pdihssel(3)=95

   pdihssel(4)=267
   pdihssel(5)=243
   pdihssel(6)=317

   pdihssel(7)=322
   pdihssel(8)=404
   pdihssel(9)=382

   pdihssel(10)=421

   do i=1,size(pdihssel,1)

      call random_number(r)

      selectedDihedralID=pdihssel(i)

      call getAtomsToRotate(bonds,pdihs,selectedDihedralID,atomIDsToRotate)

      call rotateGroupOfAtoms(atoms,pdihs,atomIDsToRotate,pdihs(selectedDihedralID),pix2*r(1))

   end do

   call saveAtomsOnPDB(atoms,files%pdbinp,files%pdbout)


contains

   subroutine getAtomsToRotate(bonds,pdihs,dihID,atomIDsToRotate)
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

end program dev_do_select
