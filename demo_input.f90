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

   use mod_dihedral
   use mod_energy

   implicit none

   type(atom_t), dimension(:), allocatable :: atoms
   type(bond_t), dimension(:), allocatable :: bonds
   type(pdih_t), dimension(:), allocatable :: pdihs,pdihs_tmp
   real(kind=real64), dimension(:,:), allocatable :: eps_M
   real(kind=real64), dimension(:,:), allocatable :: sigma_M

   integer(kind=int64), dimension(:), allocatable :: size_s,excludes
   integer(kind=int64), dimension(:,:), allocatable :: struct,connections,&
   atomsInPath,nonBondedExclusions
   integer(kind=int64), dimension(:), allocatable ::  dihsOfInter,numconns,&
   nonBondedExcl
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

   eps_M(:,:)=eps_M(:,:)/4.0_real64

   call createAtomsToRotateStruct(atoms,bonds,pdihs,dihsOfInter,struct,size_s)

   allocate(dihsOfInter(19))

   dihsOfInter=(/21,44,69,94,119,144,169,194,219,244,269,294,319,344,369,394,419,444,469/)

   call getNonBondedExcl(atoms,bonds,nonBondedExclusions)


   deallocate(size_s)
   deallocate(struct)
   deallocate(dihsOfInter)
   deallocate(nonBondedExclusions)
   deallocate(eps_M)
   deallocate(sigma_M)
   deallocate(atoms)
   deallocate(bonds)
   deallocate(pdihs)

end program demo_rotation_series
