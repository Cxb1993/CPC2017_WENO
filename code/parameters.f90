 module OCFD_precision
 implicit none
 include "mpif.h"
      integer,parameter::OCFD_REAL_KIND=8,  OCFD_DATA_TYPE=MPI_DOUBLE_PRECISION   ! double precison computing
 end module OCFD_precision
 
 
 
 !------parameters used in OpenCFD---------------------------------------
 module OCFD_constants
   Use OCFD_precision
   implicit none
!---------The common (global) variables in original opencfd.h ----------------------
!--------MPI-------------------
    integer,save:: my_id,npx,npy,npz,npx0,npy0,npz0, ID_XP1,ID_XM1,ID_YP1,ID_YM1,ID_ZP1,ID_ZM1,       &
                  MSG_BLOCK_SIZE,BLOCK_COMPACT_Z,BLOCK_COMPACT_Y, MPI_COMM_X,MPI_COMM_Y,MPI_COMM_Z,  &
                  MPI_COMM_XY,MPI_COMM_XZ,MPI_COMM_YZ, &
                  TYPE_LAPX1,TYPE_LAPY1,TYPE_LAPZ1,TYPE_LAPX2,TYPE_LAPY2,TYPE_LAPZ2,np_size
   integer,save::  i_offset(0:1024),j_offset(0:1024),k_offset(0:1024),i_nn(0:1024),j_nn(0:1024),k_nn(0:1024)
!--------------------------------
!-------------------------------------------
 end module OCFD_constants


 !--------------------------------------------------------------------------------------------------------
 !--------------------------------------------------------------------------------------------------------
  module  OCFD_flow_parameters
  use OCFD_constants
  implicit none
 !-------- integer parameters for MPI -----------------------------
        
!------------------
! integer parameters for goemtry
   integer,save::    nx_global,ny_global,nz_global,nx,ny,nz,LAP
! parameters related with   Flow
   integer,save::    IFLAG_G,NUM_TIME_Advance, &
        Iflag_grid(3),Iperiodic(3),Kstep_save,Ksave_style,Kread_style,    &
        Kstep_show,Nvars,Istep,KRK
   real(kind=OCFD_REAL_KIND),save:: dt,end_time,tt

!  Coordinate parameters 
   real(kind=OCFD_REAL_KIND),save:: SLx,SLy,SLz,hx,hy,hz
   real(kind=OCFD_REAL_KIND),allocatable,save,dimension(:) :: xx,yy,zz,xx0,yy0,zz0,sx,sy,sz,sx0,sy0,sz0

!--------------------------------------------------
  end module  OCFD_flow_parameters

 ! initial of data , OpenCFD Ver 1.8  by Li Xinliang , Institute of Mechanics, CAS 
  subroutine init(u,u1)
  Use OCFD_flow_parameters
  implicit none
  integer i,j,k,i1,k1,ierr
  real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP):: u
  real(kind=OCFD_REAL_KIND),dimension(nx,ny,nz):: u1
  real(kind=OCFD_REAL_KIND) nn1
  integer x_lxj,y_lxj,z_lxj
  integer incx,incy,incz
  nn1=nx_global*ny_global*nz_global*1.0d0
   u=0.0001d0
!   u1=0.0001d0
incx=mod(nx_global,npx0)
incy=mod(ny_global,npy0)
incz=mod(nz_global,npz0)
do k=1,nz

    if(npz<incz)then
       z_lxj=npz*nz+k
    else 
       z_lxj=npz*nz+k+incz
    endif
   do j=1,ny

       if(npy<incy)then
           y_lxj=npy*ny+j
       else 
           y_lxj=npy*ny+j+incy
       endif

        do i=1,nx
         if(npx<incx)then
             x_lxj=npx*nx+i
         else 
             x_lxj=npx*nx+i+incx
         endif

           u(i,j,k)=x_lxj*y_lxj*z_lxj/nn1-0.0001d0
           u1(i,j,k)=u(i,j,k)
       enddo
    enddo
enddo


   allocate ( xx(1-LAP:nx+LAP),yy(1-LAP:ny+LAP),zz(1-LAP:nz+LAP),  &
       xx0(nx_global),yy0(ny_global),zz0(nz_global),        &
       sx(1-LAP:nx+LAP),sy(1-LAP:ny+LAP),sz(1-LAP:nz+LAP),  &
       sx0(nx_global),sy0(ny_global),sz0(nz_global) )
!------define hx, hy, hz -----------------------
         if(Iperiodic(1).eq.0) then
          hx=SLx/(nx_global-1.d0)
         else 
          hx=SLX/nx_global
         endif  

       do k=1-LAP,nx+LAP
        k1=i_offset(npx)+k-1
        if(k1.ge.1 .and. k1 .le. nx_global) then
          if(Iflag_grid(1).eq.1) then
           xx(k)=xx0(k1)
           sx(k)=sx0(k1)
              else
           xx(k)=(k1-1.d0)*hx
           sx(k)=1.d0
            endif
        endif
       enddo
!--------------------------------------
         if(Iperiodic(2).eq.0) then
           hy=SLy/(ny_global-1.d0)
         else 
           hy=SLy/ny_global
         endif

       do k=1-LAP,ny+LAP
       k1=j_offset(npy)+k-1
       if(k1.ge.1 .and. k1 .le. ny_global) then
        if(Iflag_grid(2).eq.1) then
          yy(k)=yy0(k1)
          sy(k)=sy0(k1)
         else
           yy(k)=(k1-1.d0)*hy
          sy(k)=1.d0
         endif
       endif
       enddo
!--------------------------------------
        if(Iperiodic(3).eq.0) then
           hz=SLz/(nz_global-1.d0)
        else 
           hz=SLz/nz_global       
        endif


       do k=1-LAP,nz+LAP
       k1=k_offset(npz)+k-1
       if(k1.ge.1 .and. k1 .le. nz_global) then
        if(Iflag_grid(3).eq.1) then
           zz(k)=zz0(k1)
          sz(k)=sz0(k1)
         else
          zz(k)=(k1-1.d0)*hz
          sz(k)=1.d0
         endif
       endif
       enddo

 end
