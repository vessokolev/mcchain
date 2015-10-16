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
   use mod_distance

   implicit none

   real(kind=real64),parameter :: eps_r=1.0_real64
   real(kind=real64),parameter :: f=138.935485_real64/eps_r

   ! Note that the parameter used to handle the electrostatic
   ! interactions - f can be calculated as:
   !
   ! f=NA*c*c/4/pi/eps0
   !
   ! where NA is the Avogadro's constant:
   ! NA=6.0221409e+23 (in mol^-1),
   ! c is the charge of the electron:
   ! c=1.60217662e-19 (in C)
   ! and eps0 is the vacuum permitivity:
   ! eps0=8.854187817e-12 (in C^2.N^-1.m^-2).
   !
   ! In order to work with interatomic distances in nm:
   !
   ! f=NA*c*c/4/pi/eps0/1e-9
   !
   ! which gives f=138.935485 kJ.mol^-1.nm.e^-2
   !
   ! By involving f in the computation of the electrostatic potential the number
   ! of calculations can be reduced significantly by using the formula:
   !
   ! V(i,j) = f*z(i)*z(j)/r(i,j)
   !
   ! where z(i) and z(j) are the relative charges (in e) of the atoms and r(i,j) is the
   ! distance between the atoms.
   !

   type(atom_t), dimension(:), allocatable :: atoms
   type(bond_t), dimension(:), allocatable :: bonds
   type(pdih_t), dimension(:), allocatable :: pdihs,pdihs_tmp
   real(kind=real64), dimension(:,:), allocatable :: eps_M
   real(kind=real64), dimension(:,:), allocatable :: sigma_M
   real(kind=real64), dimension(:), allocatable :: distM

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

   print *,"Processing the input data..."
   call readInput(files,atoms,bonds,pdihs,eps_M,sigma_M)
   call createDistanceMatrix(atoms,distM)
   call getNonBondedExcl(atoms,bonds,nonBondedExclusions)

   allocate(dihsOfInter(19))

   dihsOfInter=(/21,44,69,94,119,144,169,194,219,244,269,294,319,344,369,394,419,444,469/)

   call createAtomsToRotateStruct(atoms,bonds,pdihs,dihsOfInter,struct,size_s)

   print *,"Starting the computations..."

   print *,"r=",norm2((/atoms(1)%x,atoms(1)%y,atoms(1)%z/)-&
                 (/atoms(3)%x,atoms(3)%y,atoms(3)%z/))

   print *,"sigma calc",0.5_real64*(atoms(1)%sigma+atoms(3)%sigma),&
           sigma_M(atoms(1)%sigma_g,atoms(3)%sigma_g)
   print *,"eps   calc",sqrt(atoms(1)%eps*atoms(3)%eps)*4.0_real64,&
           eps_M(atoms(1)%eps_g,atoms(3)%eps_g)

   print *,"charges :",atoms(1)%charge,atoms(3)%charge

   print *,"vdW Energy :",getvdWPairEnergy(atoms,eps_M,sigma_M,distM,1,3)
   print *,"Electrost E:",getElPairEnergy(atoms,distM,1,3)*f

!   call createDistanceMatrix(atoms,distM)

!   call updateDistanceMatrix(atoms,distM,(/2/))

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
