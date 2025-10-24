program compute_dipole
implicit none
character*256 :: frho,fcfg,fout

integer*4 :: num_ele
integer*4,allocatable :: type_list(:)
real*8,allocatable :: zval_list(:)

type config
  integer*4 :: natom
  integer*4,allocatable :: tatom(:),matom(:,:)
  real*8 :: AL(3,3)
  real*8,allocatable :: xatom(:,:)
end type config

type(config) :: ac

integer*4 :: ipln, num_pln, num_lat, iatom, n1,n2,n3
real*8 :: dipole,xtmp,xcut,dvol,a_au
real*8,allocatable :: dipole_layer(:), zval(:), &
                      dipole_rho(:),rho(:,:,:)

real*8 :: tot_zval
real*8 :: AL(3,3), AL_out(3,3)
real*8,parameter :: one_third=1.d0/3.d0
integer*4 :: new_unit,i,j,k,ierr,nnodest
real*8,allocatable :: rho_out(:,:,:)
real*8 :: V_corner,V_bulk,E_6,E6_1,E6_2,E6_3,Evbm,E_bind,E_bind_2,m33
real*8,allocatable :: occ(:)
integer*4 :: num_band
real*8,parameter :: E_AU=27.211386d0
logical :: ist

call read_rho("~/work/single_doped.B/scf.PBE.DIJ.pure/OUT.VR",rho,AL,n1,n2,n3,nnodest)
V_bulk=rho(1,1,1)*E_AU

inquire(file="OUT.VR",exist=ist)
if(ist) then
call read_rho("OUT.VR",rho,AL,n1,n2,n3,nnodest)
V_corner=rho(1,1,1)*E_AU
else
write(6,*) "OUT.VR not found, use V_bulk by default"
V_corner=V_bulk
end if

call read_occ("~/work/single_doped.B/scf.PBE.DIJ.pure/OUT.OCC",occ)
E6_1=sum(occ(1022:1024))/3.d0
E6_2=sum(occ(1019:1021))/3.d0
E6_3=sum(occ(1018:1018))/1.d0

! E_6=min(E6_1,E6_2)
E_6 =E6_1


call read_occ("OUT.OCC",occ)
E6_1=sum(occ(1022:1024))/3.d0

Evbm=E6_1

E_bind_2=E_6+V_corner-V_bulk-Evbm   ! but I think this is truth
E_bind=E_6-V_corner+V_bulk-Evbm ! suggested by LWW

call read_parameter("parameter.input",m33)

call write_output()

contains
subroutine read_occ(file,occ)
implicit none
integer*4 :: nline,iE,iline,ios
real*8 :: f
character(len=*) :: file
real*8,allocatable :: occ(:)
print *, "we only address single K-point ! no spin !"
nline=0
  open(newunit=new_unit,file=trim(file))
  rewind(new_unit)
loop_line : do
  read(new_unit,*,iostat=ios) 
  if(ios.ne.0) exit loop_line
  nline=nline+1
end do loop_line
  close(new_unit)
print *, "total line", nline
if(allocated(occ)) deallocate(occ)
allocate(occ(nline-2))
  open(newunit=new_unit,file=trim(file))
  rewind(new_unit)
  read(new_unit,*,iostat=ios) 
  read(new_unit,*,iostat=ios) 
do iline=1,nline-2
  read(new_unit,*,iostat=ios) iE, occ(iline), f 
end do
  close(new_unit)
end subroutine read_occ

subroutine read_parameter(file,m33)
implicit none
character(len=*) :: file
real*8 :: m33
open(newunit=new_unit,file=trim(file))
rewind(new_unit)
read(new_unit,*) m33
close(new_unit)
end subroutine read_parameter

subroutine write_output()
implicit none
open(newunit=new_unit,file="data.out")
rewind(new_unit)
write(new_unit,190) "mii", "E_bind(eV)=", "(E_6-Evbm)-", "(V_bulk-V_corner)", "E_bind_LWW (eV)"
write(new_unit,191) m33, E_bind, E_6-Evbm, V_bulk-V_corner, E_bind_2
191 format(5(1x,E18.10))
190 format(5(1x,A18   ))
close(new_unit)
end subroutine write_output

subroutine read_rho(file,rho,AL,n1,n2,n3,nnodest)
implicit none
character(len=*),intent(in) :: file
real*8,allocatable,intent(out) :: rho(:,:,:)
integer*4,intent(out) :: n1,n2,n3
real*8,intent(out) :: AL(3,3)
integer*4,intent(out) :: nnodest
integer*4 :: nr, nr_n, iread, ii, jj, i, j, k
real*8,allocatable :: vr_tmp(:)
open(13,file=trim(file),form="unformatted",&
    status='old',action='read',iostat=ierr)
rewind(13)
read(13) n1 ,n2 ,n3 ,nnodest
read(13) AL
if(allocated(rho)) deallocate(rho)
allocate(rho(n1,n2,n3))
nr=n1*n2*n3
nr_n=(n1*n2*n3)/nnodest
allocate(vr_tmp(nr_n))
do iread=1,nnodest
   read(13) vr_tmp
   do ii=1,nr_n
      jj=ii+(iread-1)*nr_n
      i=(jj-1)/(n2*n3)+1
      j=(jj-1-(i-1)*n2*n3)/n3+1
      k=jj-(i-1)*n2*n3-(j-1)*n3
      rho(i,j,k)=vr_tmp(ii)
   enddo
enddo
close(13)
deallocate(vr_tmp)
end subroutine read_rho

subroutine outp_rho(file,rho,AL,n1,n2,n3,nnodest)
implicit none
character(len=*),intent(in) :: file
real*8,allocatable,intent(in) :: rho(:,:,:)
integer*4,intent(in) :: n1,n2,n3
real*8,intent(in) :: AL(3,3)
integer*4,intent(in) :: nnodest
integer*4 :: nr, nr_n, iread, ii, jj, i, j, k, nnodes_tmp
real*8,allocatable :: vr_tmp(:)

nnodes_tmp=nnodest
open(13,file=trim(file),form="unformatted")
rewind(13)
write(13) n1,n2,n3,nnodes_tmp
write(13) AL
nr_n=(n1*n2*n3)/nnodes_tmp
allocate(vr_tmp(nr_n))
do iread=1,nnodes_tmp
   do ii=1,nr_n
      jj=ii+(iread-1)*nr_n
      i=(jj-1)/(n2*n3)+1
      j=(jj-1-(i-1)*n2*n3)/n3+1
      k=jj-(i-1)*n2*n3-(j-1)*n3
      vr_tmp(ii)=rho(i,j,k)
   enddo
   write(13) vr_tmp
enddo
close(13)
deallocate(vr_tmp)

end subroutine outp_rho

subroutine add_rho()
implicit none
do k=1,n3
do j=1,n2
do i=1,n1
! rho(i,j,k)=rho(i,j,k)-2.392108858132828d-05*(i-1.d0-192)
! rho(i,j,k)=rho(i,j,k)-0.035857327586521d-3*(i-1.d0-192)
! rho(i,j,k)=rho(i,j,k)-3*0.035857327586521d-3*(i-1.d0-192)
! rho(i,j,k)=rho(i,j,k)-2*0.035857327586521d-3*(i-1.d0-192)
! rho(i,j,k)=rho(i,j,k)-1.5d0*0.035857327586521d-3*(i-1.d0-192)
! rho(i,j,k)=rho(i,j,k)-1.25d0*0.035857327586521d-3*(i-1.d0-192)
rho_out(i,j,k)=rho(i,j,k)
if(i.gt.1) rho_out(2*n1+2-i,j,k)=rho(i,j,k)
end do
end do
end do
end subroutine add_rho

subroutine read_config(file,ac)
implicit none
character(len=*),intent(in) :: file
type(config),intent(out) :: ac
integer*4 :: new_unit,iatom
open(newunit=new_unit,file=trim(file))
rewind(new_unit)
read(new_unit,*) ac%natom
read(new_unit,*) ! LATTICE
read(new_unit,*) ac%AL(:,1)
read(new_unit,*) ac%AL(:,2)
read(new_unit,*) ac%AL(:,3)
read(new_unit,*) ! POSITION
allocate(ac%tatom(ac%natom),&
ac%xatom(3,ac%natom),&
ac%matom(3,ac%natom))
do iatom=1,ac%natom
read(new_unit,*) ac%tatom(iatom),ac%xatom(1:3,iatom),ac%matom(1:3,iatom) 
end do
close(new_unit)
end subroutine read_config

end program compute_dipole
