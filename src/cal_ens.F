program cal_ens

use netcdf
use mpi
use, intrinsic :: iso_c_binding

implicit none

include 'fftw3.f03'

integer, parameter :: nx=5120,ny=1,nz=600,nt=100
real, parameter :: dx=10., dz=10. 
integer :: i,j,k,m,n,o,t
integer :: ierr, ncid1, varid1
integer :: myid,nproc,ts,te,tt
integer :: idum1,idum2
integer :: is,ie,js,je,ks,ke
integer, dimension(:), allocatable :: mpi_s,mpi_n
real :: dum1,dum2,dum3,dum4,dum5
real :: start, finish
real, dimension(nz+1) :: zz,zc,rho,rhoz,thbar,thvbar,pbar,pibar,q1,q2
real, dimension(nx,ny,nz) :: ens
real, dimension(nx,ny,nz,3) :: vor,vorm
! fftw
real(c_double), dimension(nx,nz+200) :: tmp, filter
complex(c_double_complex), dimension(nx/2+1,nz+200) :: c1, c2
type(c_ptr) :: plan_forward, plan_backward
! output
real*4, dimension(:,:,:,:), allocatable :: ens_o
character(200) :: path,run,filename,workpath

call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world,nproc,ierr)
call mpi_comm_rank(mpi_comm_world,myid,ierr)

allocate(mpi_s(nproc),mpi_n(nproc))

path="/data/der0318/tur_ens/"
workpath="/data/der0318/tur_ens/"
run="les_bubble"

idum2=nt/nproc
ts=idum2*myid+1
if (ts==1) then
  ts=ts-1
  idum2=idum2+1
endif
if (myid<mod(nt,nproc)) then
  ts=ts+myid
  idum2=idum2+1
else
  ts=ts+mod(nt,nproc)
endif
te=ts+idum2-1

write(*,*) myid, ts, te

call mpi_allgather(idum2,1,mpi_int,mpi_n,1,mpi_int,mpi_comm_world,ierr)
call mpi_allgather(ts,1,mpi_int,mpi_s,1,mpi_int,mpi_comm_world,ierr)

if (myid==0) then
  allocate(ens_o(nx,nz,10,nt+1))
else
  allocate(ens_o(nx,nz,10,mpi_n(myid+1)))
endif

WRITE(filename,111) trim(path),"/",trim(run),"/fort.98"
111 FORMAT(A,A,A,A)

! read mean profiles
OPEN(10,FILE=filename)
DO i=1,188
  READ(10,*)
ENDDO
DO i=1,nz+1
  READ(10,*) j, zz(i), zc(i)
ENDDO
DO i=1,3
  READ(10,*)
ENDDO
DO i=1,nz+1
  READ(10,*) j, rho(i), thbar(i), pbar(i), pibar(i)
ENDDO
DO i=1,3
  READ(10,*)
ENDDO
DO i=1,nz+1
  READ(10,*) j, rhoz(i), rhoz(i), q1(i), q2(i)
ENDDO
DO i=1,3
  READ(10,*)
ENDDO
DO i=2,nz+1
  READ(10,*) j, rhoz(i)
ENDDO
CLOSE(10)

call mpi_barrier(mpi_comm_world,ierr)


filter=0.
do k=1,(nz+200)/2
do i=1,nx/2
  dum1=(i-1)*dx+dx/2.
  dum2=(k-1)*dz+dz/2.
  !dum3=sqrt(dum1**2.+dum2**2.)
  !if (dum3<=500.) filter(i,k)=1.
  if (dum1<=250. .and. dum2<=250.) filter(i,k)=1.
enddo
enddo
filter(nx/2+1:nx,1:(nz+200)/2)=filter(nx/2:1:-1,1:(nz+200)/2)
filter(:,(nz+200)/2+1:nz+200)=filter(:,(nz+200)/2:1:-1)
filter=filter/sum(sum(filter,2),1)


ts=100
te=100

do t=ts,te

call cpu_time(start)

tt=t-ts+1

! read 3-D data field
WRITE(filename,112) trim(path), "/", trim(run),"/archive/",&
trim(run),".L.Dynamic-",t,".nc"
112 FORMAT(6A,I6.6,A)

ierr = nf90_open(trim(filename),NF90_NOWRITE,ncid1)
IF (ierr/=nf90_noerr) WRITE(*,*) "open fail"
ierr = nf90_inq_varid(ncid1,"xi",varid1)
IF (ierr/=nf90_noerr) WRITE(*,*) "var_inq fail"
ierr = nf90_get_var(ncid1,varid1,vor(:,:,:,1),start=(/ 1,1,1,1 /),count=(/ nx,ny,nz,1 /))
IF (ierr/=nf90_noerr) WRITE(*,*) "read fail"
ierr = nf90_inq_varid(ncid1,"eta",varid1)
IF (ierr/=nf90_noerr) WRITE(*,*) "var_inq fail"
ierr = nf90_get_var(ncid1,varid1,vor(:,:,:,2),start=(/ 1,1,1,1 /),count=(/ nx,ny,nz,1 /))
IF (ierr/=nf90_noerr) WRITE(*,*) "read fail"
ierr = nf90_inq_varid(ncid1,"zeta",varid1)
IF (ierr/=nf90_noerr) WRITE(*,*) "var_inq fail"
ierr = nf90_get_var(ncid1,varid1,vor(:,:,:,3),start=(/ 1,1,1,1 /),count=(/ nx,ny,nz,1 /))
IF (ierr/=nf90_noerr) WRITE(*,*) "read fail"
ierr = nf90_close(ncid1)

! padding for z direction
do n=1,3
  tmp=0.
  tmp(:,101:100+nz)=vor(:,1,:,n)

  call dfftw_plan_dft_r2c_2d(plan_forward,nx,nz+200,tmp,c1,fftw_measure)
  call dfftw_execute_dft_r2c(plan_forward,tmp,c1)
  call dfftw_destroy_plan(plan_forward)
  
  call dfftw_plan_dft_r2c_2d(plan_forward,nx,nz+200,filter,c2,fftw_measure)
  call dfftw_execute_dft_r2c(plan_forward,filter,c2)
  call dfftw_destroy_plan(plan_forward)
  
  c1=c1*c2
  call dfftw_plan_dft_c2r_2d(plan_backward,nx,nz+200,c1,tmp,fftw_measure)
  call dfftw_execute_dft_c2r(plan_backward,c1,tmp)
  call dfftw_destroy_plan(plan_backward)

  vorm(:,1,:,n)=tmp(:,101:100+nz)/real(nx*(nz+200))
enddo



call cpu_time(finish)
write(*,*) t, finish-start
enddo 
    
ens_o(:,:,1:3,ts)=vor(:,1,:,:)
ens_o(:,:,4:6,ts)=vorm(:,1,:,:)
ens_o(:,:,7,ts)=filter(:,1:(nz+200)/2)
ens_o(:,:,8,ts)=filter(:,(nz+200)/2+1:nz+200)
ens_o(:,:,9:10,ts)=0.

open(10,file=trim(workpath)//"/data/ens_data.dat",access="direct",recl=nx*nz*10)
n=1
do t=ts,te
  write(10,rec=n) ens_o(:,:,:,t)
  n=n+1
enddo
close(10)

end program cal_ens    
