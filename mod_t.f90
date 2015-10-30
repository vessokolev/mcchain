module mod_t
   !
   ! Created by Vesselin Kolev <vesso.kolev@gmail.com>
   ! 20151005123114
   !
   implicit none
   !
   ! Define the kind of the floating point and integer parameters
   ! Note: if one uses gfortran of Intel Fortran compilers it is possible
   ! to get these numbers definied in the module iso_fortran_env. It might
   ! be called as:
   !
   ! use iso_fortran_env, only: real32, real64, int8, int16, int32, int64
   !
   ! But PGI Fortran compiler does not provided iso_fortran_env. That's why
   ! for the sake of compatibility of the code the parameters are defined in
   ! the code bellow.
   !

   integer, parameter :: real32=4
   integer, parameter :: real64=8

   integer, parameter :: int8=1
   integer, parameter :: int16=2
   integer, parameter :: int32=4
   integer, parameter :: int64=8


   type atom_t
      integer(kind=int64) :: anum
      character(5) :: atype
      integer(kind=int16) :: resnum
      character(4) :: resname
      character(4) :: aname
      integer(kind=int64) :: cgrp
      real(kind=real32) :: charge
      real(kind=real32) :: mass
      real(kind=real32) :: sigma
      real(kind=real32) :: eps
      real(kind=real32) :: x
      real(kind=real32) :: y
      real(kind=real32) :: z
      integer(kind=int16) :: sigma_g
      integer(kind=int16) :: eps_g
      integer(kind=int64), dimension(:), allocatable :: bonds
   end type atom_t


   type bond_t
      integer(kind=int64) :: left_a
      integer(kind=int64) :: right_a
      character(5) :: left_a_t
      character(5) :: right_a_t
   end type bond_t


   type file_t
      character(len=4096) :: atoms
      character(len=4096) :: bonds
      character(len=4096) :: pdihsDefs
      character(len=4096) :: pdihsParams
      character(len=4096) :: pdbcryst1
      character(len=4096) :: pdbinp
      character(len=4096) :: pdbout
      character(len=4096) :: traj
      character(len=4096) :: trajout
   end type file_t


   type pdih_p_t
      real(kind=real32) :: phi0
      real(kind=real32) :: kd
      integer(kind=int8) :: mult
   end type pdih_p_t


   type pdih_t
      integer(kind=int64) :: id
      integer(kind=int64), dimension(4) :: member_id
      real(kind=real32) :: currentAngle
      real(kind=real32) :: energy
      type(pdih_p_t), dimension(:), allocatable :: params
   end type pdih_t


   type pdb_cryst1_t
      character(len=6) :: rname
      real(kind=real32) :: a
      real(kind=real32) :: b
      real(kind=real32) :: c
      real(kind=real32) :: alpha
      real(kind=real32) :: beta
      real(kind=real32) :: gamma
      character(len=11) :: sgroup
      integer(kind=int16) :: zval
   end type pdb_cryst1_t


   type pdb_atom_t
      character(len=6) :: ident
      integer(kind=int16) :: anum
      character(len=4) :: aname
      character(len=1) :: altloc
      character(len=3) :: resname
      character(len=2) :: chainid
      integer(kind=int16) :: resnum
      character(len=1) :: inscode
      real(kind=real32) :: x
      real(kind=real32) :: y
      real(kind=real32) :: z
      real(kind=real32) :: occup
      real(kind=real32) :: tempf
      character(len=4) :: segmentid
      character(len=2) :: elsymb
      character(len=2) :: charge
   end type pdb_atom_t

   real(kind=real32), parameter :: &
   pi=3.141592653589793116_real32
   real(kind=real32), parameter :: &
   pix2=6.283185307179586232_real32
   real(kind=real32), parameter :: &
   pix3div2=4.712388980384689674_real32

   !
   ! The following parameters should be defined:
   !
   ! - in the module mod_input.f90 (see subroutine readInput there):
   ! atoms_s=size(atoms,1)
   ! bonds_s=size(bonds,1)
   ! pdihs_s=size(pdihs,1)
   !
   integer(kind=int64) :: atoms_s,bonds_s,pdihs_s


contains

end module mod_t
