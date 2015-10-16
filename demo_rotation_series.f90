program demo_rotation_series

   use mod_t       ! Defines the base types
   use mod_pdb     ! Reads and write PDB files
   use mod_rot     ! Rotates atomic coordinates
   use mod_input   ! Reads the input files and creates 
                   ! the data structures into the memory.
   use mod_select  ! Selects the atomic coordinates that
                   ! need to be rotated.
   use mod_array   ! Helper module that assists the operations
                   ! with arrays

   implicit none

   type(atom_t), dimension(:), allocatable :: atoms
   type(bond_t), dimension(:), allocatable :: bonds
   type(pdih_t), dimension(:), allocatable :: pdihs
   real(kind=real64), dimension(:,:), allocatable :: eps_M
   real(kind=real64), dimension(:,:), allocatable :: sigma_M

   integer(kind=int64), dimension(:), allocatable :: size_s
   integer(kind=int64), dimension(:,:), allocatable :: struct
   integer(kind=int64), dimension(:), allocatable ::  dihsOfInter
   integer(kind=int16), dimension(1) :: seed=(/3232/)
   integer(kind=int64) :: i
   type(file_t) :: files
   real(kind=real64) :: r, start, finish

   !
   ! The files bellow are created by running create_input.py
   ! 
   files%atoms="/home/vesso/Code/MC_dih/current/files/atoms.txt"
   files%bonds="/home/vesso/Code/MC_dih/current/files/bonds.txt"
   files%pdihsDefs="/home/vesso/Code/MC_dih/current/files/pdhdr.txt"
   files%pdihsParams="/home/vesso/Code/MC_dih/current/files/pdhdr-params.txt"
   !
   ! The file opt.pdb contains the initial coordinates of the atoms.
   ! Note that it is the same file that create_input.py reads in order
   ! to get the atomic coordinates.
   ! The file tmp.pdb conatins the result of the torsion angles changes.
   !
   files%pdbinp="/home/vesso/Code/MC_dih/current/files/opt.pdb"
   files%pdbout="/home/vesso/Code/MC_dih/current/files/tmp.pdb"

   !
   ! defined in mod_input.f90
   !
   ! This routine reads the set of input files supplied by create_input.py.
   ! It creates into the memory arrays of types for the atoms, bonds, and
   ! the proper dihedrals. It also calculates the matrix of eps and sigma
   ! vdW parameters that later are used to calculate the energy.
   !
   call readInput(files,atoms,bonds,pdihs,eps_M,sigma_M)

   allocate(dihsOfInter(19))

   !
   ! This array conatins the indexes of the dihedral angles, which j and k atoms
   ! define the axis of rotation. The sequence bellow contains only the proper
   ! dihedrals which axis of rotation includes C-alpha atoms.
   !
   dihsOfInter=(/21,44,69,94,119,144,169,194,219,244,269,294,319,344,369,394,419,444,469/)
   
   call random_seed (PUT=seed)

   !
   ! defined in: mod_select.f90
   !
   ! This routine creates a matrix containing the indexes of atoms that should
   ! be rotated with over the defines axis. Each row of the matrix contains the
   ! atom IDs corresponding to the roration over one dihedral axis. See the
   ! comment in mod_select.f90.
   !
   call createAtomsToRotateStruct(atoms,bonds,pdihs,dihsOfInter,struct,size_s)

   do i=1,size(dihsOfInter,1)

      call cpu_time(start)

      call random_number(r) 

      !
      ! defined in: mod_rot.f90
      !
      ! Rotates a series of atomic coordinates over a given axis.
      !
      call rotateGroupOfAtoms(atoms,pdihs,struct(i,1:size_s(i)),pdihs(dihsOfInter(i)),pix2*r)

      call cpu_time(finish)

      print *,"CPU time per rotation of",size_s(i),"atoms:",finish-start,"s"

   end do

   !
   ! defined in: mod_pdb.f90
   !
   ! Writes the transformed coordinates in PDB format.
   !
   call saveAtomsOnPDB(atoms,files%pdbinp,files%pdbout)

end program demo_rotation_series
