program compute_ebind 
implicit none
! use the outermost box surface mesh grid average of the Vbulk and Vcorner
! rather than the single corner point

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

integer*4 :: ipln, num_pln, num_lat, iatom, n1,n2,n3, m1, m2, m3
real*8 :: dipole,xtmp,xcut,dvol,a_au
real*8,allocatable :: dipole_layer(:), zval(:), &
                      dipole_rho(:),rho(:,:,:), rho_bulk(:,:,:), rho_dope(:,:,:)

real*8 :: tot_zval
real*8 :: AL(3,3), AL_out(3,3)
real*8,parameter :: one_third=1.d0/3.d0
integer*4 :: new_unit,i,j,k,ierr,nnodest
real*8,allocatable :: rho_out(:,:,:)
real*8 :: V_corner,V_bulk,E_6,E6_1,E6_2,E6_3,Evbm,E_bind,E_bind2,E_d, E_7, Ecbm, E_e
real*8,allocatable :: occ(:,:)
integer*4 :: num_band, ith_eigen, ith_bedge
real*8,parameter :: E_AU=27.211386d0
logical :: ist
character*256 :: band_type
real*8 :: sign_e
logical :: iqst
integer*4 :: n1_ref, n2_ref, n3_ref
integer*4 :: n1_tmp, n2_tmp, n3_tmp
character*256 :: dir_ref,fref,dir_tmp,ftmp,m33
integer*4 :: ikpt_vbm, ikpt_cbm, ikpt_acc, ikpt_uoc
real*8 :: xs(3),xd(3)
integer*4 :: ith_cv
logical :: short_occ  ! missing Ec (for acceptor case) or Ev (for donor case) levels in occ
call read_parameter("parameter.input",m33)

select case (trim(band_type))
  case ('VBM')
    sign_e=1.d0
  case ('CBM')
    sign_e=-1.d0
  case DEFAULT
    sign_e=1.d0
    write(6,*) "band_type = VBM or CBM, please check", trim(band_type)
    stop
end select


write(fref,'(A,A)') trim(dir_ref),"/OUT.VR"

inquire(file=trim(fref),exist=iqst)

if(.not.iqst) then
  write(6,*) trim(fref), "not exist, stop"
  stop
end if

call read_rho(fref,rho_bulk,AL,n1,n2,n3,nnodest)


write(ftmp,'(A,A)') trim(dir_tmp),"/OUT.VR"

inquire(file=trim(ftmp),exist=ist)
if(ist) then
call read_rho(ftmp,rho_dope,AL,m1,m2,m3,nnodest)
! V_corner=rho_dope(n1_tmp, n2_tmp, n3_tmp)*E_AU
xd=[1.d0, 1.d0, 1.d0]/2.d0
xs=[1.d0, 1.d0, 1.d0]/2.d0
V_corner = bsa_bulk(m1,m2,m3,m1,m2,m3,xd,xs,rho_dope) * E_AU
print *, "V_corner", V_corner
else
! write(6,*) "OUT.VR not found, use V_bulk by default"
! V_corner=V_bulk
write(6,*) "OUT.VR not found, stop"
stop
end if

! single point scheme
! V_bulk=rho_bulk(n1_ref, n2_ref, n3_ref )*E_AU

! box shell average scheme
xd=[0.d0, 0.d0, 0.d0]
xs=[1.d0, 1.d0, 1.d0]/2.d0
V_bulk = bsa_bulk(n1,n2,n3,m1,m2,m3,xd,xs,rho_bulk) * E_AU
print *, "V_bulk", V_bulk  

write(fref,'(A,A)') trim(dir_ref),"/OUT.OCC"

inquire(file=trim(fref),exist=iqst)

if(.not.iqst) then
  write(6,*) trim(fref), "not exist, stop"
  stop
end if

call read_occ(fref,occ)
! ith_bedge=1024
E_e=occ(ith_bedge, ikpt_vbm)
E_d=occ(ith_bedge+nint(sign_e), ikpt_cbm)

write(ftmp,'(A,A)') trim(dir_tmp),"/OUT.OCC"

inquire(file=trim(ftmp),exist=iqst)

if(.not.iqst) then
  write(6,*) trim(ftmp), "not exist, stop"
  stop
end if


call read_occ(ftmp,occ)
! ith_eigen=1024
E_6=occ(ith_eigen,ikpt_acc)
ith_cv=ith_eigen+nint(sign_e)
if(ith_cv.lt.1.or.ith_cv.gt.size(occ(:,ikpt_uoc))) then
  short_occ=.true.
else
  short_occ=.false.
  E_7=occ(ith_cv, ikpt_uoc)
end if


E_bind=sign_e*(E_6-V_corner+V_bulk-E_e) ! suggested by LWW

if(.not.short_occ) then
E_bind2 = sign_e*(E_6 - (E_7 - (E_d - E_e) )) ! use the unocc orbit as reference
else
E_bind2 = -999.d0
end if

call write_output()

contains

function bsa_bulk(n1,n2,n3, m1,m2,m3, xd, xs, rho) result(avg)
implicit none
integer*4 :: n1,n2,n3
integer*4 :: m1,m2,m3
real*8 :: xd(3), xs(3), avg
real*8 :: rho(n1,n2,n3)
integer*4 :: nx, ny, nz, rnx, rny, rnz, nc1, nc2, nc3, cnt
logical :: shell_x, shell_y, shell_z
integer*4 :: k,j,i
integer*4 :: w,v,u
nx=m1/n1
ny=m2/n2
nz=m3/n3
if(&
m1.ne.nx*n1.or.&
m2.ne.ny*n2.or.&
m3.ne.nz*n3) then
  write(6,*) "not evenly divided m  /  n"
  write(6,*) "m123=", m1,m2,m3
  write(6,*) "n123=", n1,n2,n3
  stop
end if
if(&
mod(n1,2).ne.0.or.&
mod(n2,2).ne.0.or.&
mod(n3,2).ne.0) then
  write(6,*) "not evenly divided n / 2"
  write(6,*) "n123=", n1,n2,n3
  stop
end if
rnx=m1/2
rny=m2/2
rnz=m3/2
nc1=nint(m1*xs(1))
nc2=nint(m2*xs(2))
nc3=nint(m3*xs(3))
cnt=0
do k=1, m3+1
  shell_z=nint(dabs(dist_min((k-1-nc3)*1.d0/m3)*m3)).eq.rnz
  w=mod(k-1-nc3+100*n3,n3)+ nint(xd(3)*n3)
  w=mod(w+100*n3,n3)+1
  do j=1, m2+1
    shell_y=nint(dabs(dist_min((j-1-nc2)*1.d0/m2)*m2)).eq.rny
    v=mod(j-1-nc2+100*n2,n2)+ nint(xd(2)*n2)
    v=mod(v+100*n2,n2)+1
    do i=1, m1+1
      shell_x=nint(dabs(dist_min((i-1-nc1)*1.d0/m1)*m1)).eq.rnz
      u=mod(i-1-nc1+100*n1,n1)+ nint(xd(1)*n1)
      u=mod(u+100*n1,n1)+1
      if(shell_z.or.shell_y.or.shell_x) then
        avg=avg+rho(u,v,w)
        cnt=cnt+1
      end if
    end do
  end do
end do
avg=avg/cnt

end function bsa_bulk 

function dist_min(x) result(y)
        implicit none
        real*8,intent(in) :: x
        real*8 :: y
        y=mod(mod(x,1.d0)+1.5d0,1.d0)-0.5d0
end function dist_min


subroutine read_occ(file,occ)
implicit none
integer*4 :: nline,iE,iline,ios
real*8 :: f
character(len=*) :: file
real*8,allocatable :: occ(:,:)
character*256 :: lines
integer*4 :: ikpt, nkpt,neig,pos
! print *, "we only address single K-point ! no spin !"
print *, "we only address no spin ! multiple K-point allowed now"
nline=0
nkpt=0
  open(newunit=new_unit,file=trim(file))
  rewind(new_unit)
loop_line : do
  read(new_unit,*,iostat=ios) lines
  pos=index(lines,"KPOINTS")
  if(ios.ne.0) exit loop_line
  if(pos>0) nkpt=nkpt+1
  nline=nline+1
end do loop_line
  close(new_unit)
print *, "total line", nline
print *, "num kpoint", nkpt 
neig=nline/nkpt-2
if(nkpt*(neig+2).ne.nline) then
  write(6,*) "wrong line count, quit"
  stop
end if 
if(allocated(occ)) deallocate(occ)
allocate(occ(neig,nkpt))
  open(newunit=new_unit,file=trim(file))
  rewind(new_unit)
do ikpt=1,nkpt
  read(new_unit,*,iostat=ios) 
  read(new_unit,*,iostat=ios) 
do iline=1,neig
  read(new_unit,*,iostat=ios) iE, occ(iline,ikpt), f 
end do
end do
  close(new_unit)
end subroutine read_occ

subroutine read_parameter(file,m33)
implicit none
character(len=*) :: file
!real*8 :: m33
character(len=*) :: m33
open(newunit=new_unit,file=trim(file))
rewind(new_unit)
read(new_unit,'(a)') dir_ref
read(new_unit,*) n1_ref, n2_ref, n3_ref
read(new_unit,*) ith_bedge, ikpt_vbm, ikpt_cbm  ! what type of this ith bedge is
read(new_unit,*) band_type  ! what type of this ith bedge is
read(new_unit,'(a)') dir_tmp
read(new_unit,*) n1_tmp, n2_tmp, n3_tmp
read(new_unit,*) ith_eigen, ikpt_acc, ikpt_uoc
read(new_unit,'(A)') m33
close(new_unit)
write(6,*) "band_type = ", trim(band_type)
if(trim(band_type)=="VBM") then
write(6,*) "band_edge, ikpt(VBM), ikpt(CBM) = ", ith_bedge, ikpt_vbm, ikpt_cbm
write(6,*) "defect ei, ikpt(acc), ikpt(CBM) = ", ith_eigen, ikpt_acc, ikpt_uoc
end if
if(trim(band_type)=="CBM") then
write(6,*) "band_edge, ikpt(CBM), ikpt(VBM) = ", ith_bedge, ikpt_vbm, ikpt_cbm
write(6,*) "defect ei, ikpt(dnr), ikpt(VBM) = ", ith_eigen, ikpt_acc, ikpt_uoc
end if

write(6,*) "mii= ", trim(m33)
end subroutine read_parameter

subroutine write_output()
implicit none
open(newunit=new_unit,file="data.out")
rewind(new_unit)
write(new_unit,190) "mii", "E_bind2 (eV)", "E_bind(eV)=", "[(E_6- E_e) +", "(V_bulk-V_corner)]*", "sign_e"
write(new_unit,191) trim(m33), E_bind2, E_bind, E_6-E_e, V_bulk-V_corner, sign_e   
191 format(1x,A19,5(1x,E19.10))
190 format(6(1x,A19   ))
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

end program compute_ebind 
