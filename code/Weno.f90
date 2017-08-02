subroutine weno
   Use OCFD_flow_parameters
   implicit none
   integer i,j,k,m, ierr,Negative_T
   real*8 wstar,wend,wstar0
   real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:) :: u,u1
! allocate 
    allocate(u(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP),stat=ierr)
    allocate(u1(nx,ny,nz),stat=ierr)
!----------------initial---------------------------------------------------------
    call init(u,u1) 
    call exchange_boundary_x_standard(u,Iperiodic(1))
    call exchange_boundary_y_standard(u,Iperiodic(2))
    call exchange_boundary_z_standard(u,Iperiodic(3))
!-----------------------------------------------------------------------
    wstar=MPI_wtime() 
    wstar0=wstar              ! Initial Wall time
!-----------------------------------------------------------------------        
100      continue       
    call weno7(u,u1,nx,ny,nz,hx,LAP,dt)
do k=1,nz
do j=1,ny
do i=1,nx
  u(i,j,k)=u1(i,j,k)
enddo
enddo
enddo
    call exchange_boundary_x_standard(u,Iperiodic(1))
    call exchange_boundary_y_standard(u,Iperiodic(2))
    call exchange_boundary_z_standard(u,Iperiodic(3))
    Istep=Istep+1
        tt=tt+dt
!c ------print CPU time  -------------------------------------------
    if(mod(Istep,Kstep_show).eq.0 ) then     
    if(my_id.eq.0)  wend=MPI_wtime()
    if(my_id.eq.0)  then
          print*, '-------Istep=',Istep
          print*, 'CPU (wall) time for this step is',wend-wstar
	  print*, 'Total CPU (wall) time is', wend-wstar0
    endif
    if(my_id .eq. 0) wstar=wend
    endif
call MPI_barrier(MPI_COMM_WORLD,ierr)
!-----------save data---------------------------------------------
    if(mod(Istep,Kstep_save).eq.0) then
         call OCFD_save(Istep,tt,Ksave_style,IFLAG_G,   &
          u,nx,ny,nz,LAP,nx_global,ny_global,nz_global,0,Istep,dt,end_time,Kstep_save)
      if(end_time .le. 0.d0) goto 199      ! end_time .le. 0  means that stop computation just after saving files 
!-------------------------------------------------------------------------
    endif
!-----------------------------------------
      if( tt .lt. end_time  .or. end_time .le. 0.d0) goto 100
!--------------------------------------------------------------------------------------
199  continue
    if(my_id.eq.0) then 
        print*, 'OK! Weno finished'
    endif
  end
     subroutine weno7(f,fx,nx,ny,nz,hx,LAP,dt)
       Use OCFD_constants
       implicit none
       integer nx,ny,nz,LAP,i,j,k
           real(kind=OCFD_REAL_KIND)::f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), fx(nx,ny,nz)
       real*8 dt
           real(kind=OCFD_REAL_KIND)::hx,hx_1,hx_6,hx_12,hx_60,hx_420,hx_840,hx_2520
     


     integer i_begin,i_end
      real(kind=OCFD_REAL_KIND):: hj(0:nx), S0,S1,S2,S3,  &
             s10,s11,s12,s13,s20,s21,s22,s23,s30,s31,s32,s33,  &
             a0,a1,a2,a3,am,q0,q1,q2,q3,ep
      real(kind=OCFD_REAL_KIND),parameter::&   
         C0=1.d0/35.d0, C1=12.d0/35.d0, C2=18.d0/35.d0,  C3=4.d0/35.d0,&
         a11=-2.d0/6.d0,a12=9.d0/6.d0,a13=-18.d0/6.d0,a14=11.d0/6.d0, &
         a21=1.d0/6.d0,               a23=3.d0/6.d0,a24=2.d0/6.d0, &
         a31=-2.d0/6.d0,a32=-3.d0/6.d0,            a34=-1.d0/6.d0,   &   
         a41=-11.d0/6.d0,a42=18.d0/6.d0,a43=-9.d0/6.d0,a44=2.d0/6.d0,  &
         b12=4.d0,b13=-5.d0,b14=2.d0, &    
         b22= -2.d0,      &   
         b41=2.d0,b42=-5.d0,b43=4.d0,  &    
         c12=3.d0, &
        d12=13.d0/12.d0,d13=1043.d0/960.d0,d14=1.d0/12.d0
      real(kind=OCFD_REAL_KIND) e11,e12,e13,e14,e21,e22,e23,e24,e31,e32,e33,e34,e41,e42,e43,e44



!----constant-------------------------------------
           hx_1=1.d0/hx
           hx_6=1.d0/(6.d0*hx)
           hx_12=1.d0/(12.d0*hx)
           hx_60=1.d0/(60.d0*hx)
           hx_420=1.d0/(420.d0*hx)
           hx_840=1.d0/(840.d0*hx)
           hx_2520=1.d0/(2520.d0*hx)

! -------- inner scheme -------------------------          


         e11=-3.d0/12.d0*hx_1;e12=13.d0/12.d0*hx_1;e13=-23.d0/12.d0*hx_1;e14=25.d0/12.d0*hx_1
         e21=1.d0/12.d0*hx_1;e22=-5.d0/12.d0*hx_1;e23=13.d0/12.d0*hx_1;e24=3.d0/12.d0*hx_1
         e31=-1.d0/12.d0*hx_1;e32=7.d0/12.d0*hx_1;e33=7.d0/12.d0*hx_1;e34=-1.d0/12.d0*hx_1
         e41=3.d0/12.d0*hx_1;e42=13.d0/12.d0*hx_1;e43=-5.d0/12.d0*hx_1;e44=1.d0/12.d0*hx_1

         ep=1.d-16    !! WENO-JS
      i_begin=1
      i_end=nx
      do k=1,nz
      do j=1,ny
      do i=i_begin,i_end+1
         S10=a11*f(i-3,j,k)+a12*f(i-2,j,k)+a13*f(i-1,j,k) +a14*f(i,j,k)
         S11=a21*f(i-2,j,k) -   f(i-1,j,k)+a23*f(i,j,k)+a24*f(i+1,j,k)
         S12=a31*f(i-1,j,k)+a32*f(i,j,k)  +    f(i+1,j,k)+a34*f(i+2,j,k)
         S13=a41*f(i,j,k)  +a42*f(i+1,j,k)+a43*f(i+2,j,k)+a44*f(i+3,j,k)
         S20=-f(i-3,j,k)+b12*f(i-2,j,k)+b13*f(i-1,j,k)+b14*f(i,j,k)
         S21=             f(i-1,j,k)+b22*f(i,j,k)  +f(i+1,j,k)
         S22=             f(i,j,k)  +b22*f(i+1,j,k)+f(i+2,j,k)
         S23=b41*f(i,j,k)+b42*f(i+1,j,k)+b43*f(i+2,j,k)-f(i+3,j,k)
         S30=-f(i-3,j,k)+c12*(f(i-2,j,k)-f(i-1,j,k)) +f(i,j,k)
         S31=-f(i-2,j,k)+c12*(f(i-1,j,k)-f(i,j,k))   +f(i+1,j,k)
         S32=-f(i-1,j,k)+c12*(f(i,j,k)-f(i+1,j,k))   +f(i+2,j,k)
         S33=-f(i,j,k)  +c12*(f(i+1,j,k)-f(i+2,j,k)) +f(i+3,j,k)

       S0=S10*S10+d12*S20*S20  +d13*S30*S30 +d14*S10*S30
       S1=S11*S11+d12*S21*S21  +d13*S31*S31 +d14*S11*S31
       S2=S12*S12+d12*S22*S22  +d13*S32*S32 +d14*S12*S32
       S3=S13*S13+d12*S23*S23  +d13*S33*S33 +d14*S13*S33
       a0=C0/((ep+S0)*(ep+S0))
       a1=C1/((ep+S1)*(ep+S1))
       a2=C2/((ep+S2)*(ep+S2))
       a3=C3/((ep+S3)*(ep+S3))
     am=a0+a1+a2+a3
     q0=e11*f(i-3,j,k)+e12*f(i-2,j,k)+e13*f(i-1,j,k) +e14*f(i,j,k)
     q1=e21*f(i-2,j,k)+e22*f(i-1,j,k)+e23*f(i,j,k)   +e24*f(i+1,j,k)
     q2=e31*f(i-1,j,k)+e32*f(i,j,k)  +e33*f(i+1,j,k) +e34*f(i+2,j,k)
     q3=e41*f(i,j,k)  +e42*f(i+1,j,k)+e43*f(i+2,j,k) +e44*f(i+3,j,k)
     hj(i-1)=(a0*q0+a1*q1+a2*q2+a3*q3)/am

    enddo
       do i=i_begin,i_end
       fx(i,j,k)=fx(i,j,k)+(hj(i)-hj(i-1))*dt
       enddo
     enddo
     enddo
     end  subroutine weno7
