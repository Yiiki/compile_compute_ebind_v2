program nonscf2occ
implicit none
integer*4 :: new_unit,ios,ntot,nhig,pos,i
character*256 :: str_key, str
real*8,allocatable :: eigs(:)
real*8 :: occ
str_key="eigen energies, in ev"
call getarg(1,str)
read(str,*) ntot
allocate(eigs(ntot))
call getarg(2,str)
read(str,*) nhig  ! highest occ level index

open(newunit=new_unit,file="out")
rewind(new_unit)
do 
read(new_unit,'(A)',iostat=ios) str
if(ios.ne.0) exit

pos=index(str,trim(str_key))
if(pos.gt.0) then
read(new_unit,*,iostat=ios) (eigs(i), i=1, ntot)
if(ios.ne.0) then
  write(6,*) "something wrong when reading eigs, stop", ios
  stop
end if
end if

end do
close(new_unit)

open(newunit=new_unit,file="OUT.OCC.NONSCF")
rewind(new_unit)
write(new_unit,*) "KPOINTS      1:    0.0000    0.0000    0.0000"
write(new_unit,*) " NO.   ENERGY(eV) OCCUPATION"
do i=1,ntot
occ=2.d0
if(i.gt.nhig) occ=0.d0
write(new_unit,141) i, eigs(i), occ
141 format(I7,F12.4,F10.5)
end do
close(new_unit)

stop
end program nonscf2occ
