program demo_rotation_series

   use mod_t
   use mod_pdb
   use mod_rot
   use mod_input
   use mod_select
   use mod_array

   implicit none

   type(atom_t), dimension(:), allocatable :: atoms
   type(bond_t), dimension(:), allocatable :: bonds
   type(pdih_t), dimension(:), allocatable :: pdihs
   real(kind=real64), dimension(:,:), allocatable :: eps_M
   real(kind=real64), dimension(:,:), allocatable :: sigma_M

   real(kind=real64) :: phi,r
   integer(kind=int64), dimension(:,:), allocatable :: bonds_f
   integer(kind=int64), dimension(:), allocatable :: atomIDsToRotate,size_s
   integer(kind=int64), dimension(:,:), allocatable :: selection,struct
   integer(kind=int64), dimension(:), allocatable :: pdihssel
   integer(kind=int64), dimension(:), allocatable ::  dihsOfInter
   integer(kind=int16), dimension(1) :: seed=(/3232/)
   real(kind=real64), dimension(3) :: tmp
   integer(kind=int64) :: i,selectedDihedralID
   type(pdih_t) :: selectedDihedral
   type(file_t) :: files
   real(kind=real64) :: start, finish


   files%atoms="/home/vesso/Code/MC_dih/current/files/atoms.txt"
   files%bonds="/home/vesso/Code/MC_dih/current/files/bonds.txt"
   files%pdihsDefs="/home/vesso/Code/MC_dih/current/files/pdhdr.txt"
   files%pdihsParams="/home/vesso/Code/MC_dih/current/files/pdhdr-params.txt"
   files%pdbinp="/home/vesso/Code/MC_dih/current/files/opt.pdb"
   files%pdbout="/home/vesso/Code/MC_dih/current/files/tmp.pdb"

   call readInput(files,atoms,bonds,pdihs,eps_M,sigma_M)

   allocate(dihsOfInter(19))

   dihsOfInter=(/21,44,69,94,119,144,169,194,219,244,269,294,319,344,369,394,419,444,469/)
   
   call random_seed (PUT=seed)

   call createAtomsToRotateStruct(bonds,pdihs,dihsOfInter,struct,size_s)

   print *,shape(struct)

   do i=1,size(dihsOfInter,1)

      call cpu_time(start)

      call random_number(r) 

      call rotateGroupOfAtoms(atoms,pdihs,struct(i,1:size_s(i)),pdihs(dihsOfInter(i)),pix2*r)

      call cpu_time(finish)

      print *,"CPU time per segment rotation:",finish-start,"s"

   end do

   call saveAtomsOnPDB(atoms,files%pdbinp,files%pdbout)


contains

end program demo_rotation_series
