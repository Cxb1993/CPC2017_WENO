!OpenCFD ver 1.4, CopyRight by Li Xinliang, LNM, Institute of Mechanics, CAS, Beijing, Email: lixl@lnm.imech.ac.cn
!MPI Subroutines, such as computational domain partation, MPI message send and recv   

     subroutine part
! Domain partation----------------------------------------------------------------------------
     Use OCFD_flow_parameters
     implicit none
     integer k,ka,ierr
     integer npx1,npy1,npz1,npx2,npy2,npz2,my_mod1
!---------------------------------------------------------------------------------------------
        if(np_size .ne. npx0*npy0*npz0) then
         if(my_id.eq.0) print*, 'The Number of total Processes is not equal to npx0*npy0*npz0 !'
         call mpi_finalize(ierr)
         stop
        endif
         
        npx=mod(my_id,npx0)
        npy=mod(my_id,npx0*npy0)/npx0
        npz=my_id/(npx0*npy0)
!------commonicators-----------------------------------------------------------------------------
      
       CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,  npz*npx0*npy0+npy*npx0, npx,MPI_COMM_X,ierr)   ! 1-D 
       CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,  npz*npx0*npy0+npx, npy,MPI_COMM_Y,ierr)
       CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,  npy*npx0+npx, npz,MPI_COMM_Z,ierr)
       CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,  npz, npy*npx0+npx,MPI_COMM_XY,ierr)            ! 2-D 
       CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,  npy, npz*npx0+npx,MPI_COMM_XZ,ierr)
       CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,  npx, npz*npy0+npy,MPI_COMM_YZ,ierr)
      

!------------------------------------------------------------------------------------------------
!      均匀分配网格， 如果nx_global不能被npx0整除，将余下的网格点分到靠前的节点
!------------------------------------------------------------------------------------------------
      nx=nx_global/npx0
      ny=ny_global/npy0
      nz=nz_global/npz0
      if(npx .lt. mod(nx_global,npx0)) nx=nx+1
      if(npy .lt. mod(ny_global,npy0)) ny=ny+1
      if(npz .lt. mod(nz_global,npz0)) nz=nz+1

!------npx=k的节点上x方向网格点的个数，起始位置
!--------------------------------------------------------------------   
     do k=0,npx0-1
        ka=min(k,mod(nx_global,npx0))
        i_offset(k)=int(nx_global/npx0)*k+ka+1
        i_nn(k)=nx_global/npx0
        if(k .lt. mod(nx_global,npx0)) i_nn(k)=i_nn(k)+1
     enddo

     do k=0,npy0-1
        ka=min(k,mod(ny_global,npy0))
        j_offset(k)=int(ny_global/npy0)*k+ka+1
        j_nn(k)=ny_global/npy0
        if(k .lt. mod(ny_global,npy0)) j_nn(k)=j_nn(k)+1
     enddo

     do k=0,npz0-1
        ka=min(k,mod(nz_global,npz0))
        k_offset(k)=int(nz_global/npz0)*k+ka+1
        k_nn(k)=nz_global/npz0
        if(k .lt. mod(nz_global,npz0)) k_nn(k)=k_nn(k)+1
     enddo
!--------------------------------------------------------------------------------
!-------New Data TYPE------------------------------------------------------------
       call New_MPI_datatype

!--------define proc id:  the right, lift, up, bottom, frint and backward  procs
       npx1=my_mod1(npx-1,npx0)
       npx2=my_mod1(npx+1,npx0)
       ID_XM1=npz*(npx0*npy0)+npy*npx0+npx1    ! -1 proc in x-direction
       ID_XP1=npz*(npx0*npy0)+npy*npx0+npx2    ! +1 proc in x-direction
       if(Iperiodic(1) .eq.0 .and. npx .eq. 0) ID_XM1=MPI_PROC_NULL     ! if not periodic, 0 node donot send mesg to npx0-1 node
       if(Iperiodic(1) .eq.0 .and. npx .eq. npx0-1) ID_XP1=MPI_PROC_NULL
      
       npy1=my_mod1(npy-1,npy0)
       npy2=my_mod1(npy+1,npy0)
       ID_YM1=npz*(npx0*npy0)+npy1*npx0+npx
       ID_YP1=npz*(npx0*npy0)+npy2*npx0+npx
       if(Iperiodic(2).eq.0 .and. npy .eq. 0) ID_YM1=MPI_PROC_NULL     ! if not periodic, 0 node donot send mesg to npy0-1 node
       if(Iperiodic(2).eq.0 .and. npy .eq. npy0-1) ID_YP1=MPI_PROC_NULL

       npz1=my_mod1(npz-1,npz0)
       npz2=my_mod1(npz+1,npz0)
       ID_ZM1=npz1*(npx0*npy0)+npy*npx0+npx
       ID_ZP1=npz2*(npx0*npy0)+npy*npx0+npx
       if(Iperiodic(3).eq.0 .and. npz .eq. 0) ID_ZM1=MPI_PROC_NULL     ! if not periodic, 0 node donot send mesg to npz0-1 node
       if(Iperiodic(3).eq.0 .and. npz .eq. npz0-1) ID_ZP1=MPI_PROC_NULL


!--------------------------------------------------------------
       call MPI_barrier(MPI_COMM_WORLD,ierr)

      end            

!--------------------------------------------------------------------------------
         function my_mod1(i,n)
         implicit none
         integer my_mod1,i,n
           if(i.lt.0) then
             my_mod1=i+n
           else if (i.gt.n-1) then
             my_mod1=i-n
           else
             my_mod1=i
           endif
         end
!-----------------------------------------------------------------------------------------------
! Send Recv non-continuous data using derivative data type
   subroutine New_MPI_datatype
   Use OCFD_flow_parameters
   implicit none
   integer:: ierr,TYPE_tmp
  
   call MPI_TYPE_Vector(ny,LAP,nx+2*LAP,OCFD_DATA_TYPE,TYPE_LAPX1,ierr)
   call MPI_TYPE_Vector(LAP,nx,nx+2*LAP,OCFD_DATA_TYPE,TYPE_LAPY1,ierr)
   call MPI_TYPE_Vector(LAP,nx,(nx+2*LAP)*(ny+2*LAP),OCFD_DATA_TYPE,TYPE_LAPZ1,ierr)

   call MPI_TYPE_Vector(ny,nx,nx+2*LAP,OCFD_DATA_TYPE,TYPE_tmp,ierr)

   call MPI_TYPE_HVector(nz,1,(nx+2*LAP)*(ny+2*LAP)*OCFD_REAL_KIND,TYPE_LAPX1,TYPE_LAPX2,ierr)
   call MPI_TYPE_HVector(nz,1,(nx+2*LAP)*(ny+2*LAP)*OCFD_REAL_KIND,TYPE_LAPY1,TYPE_LAPY2,ierr)
   call MPI_TYPE_HVector(LAP,1,(nx+2*LAP)*(ny+2*LAP)*OCFD_REAL_KIND,TYPE_tmp,TYPE_LAPZ2,ierr)

   call MPI_TYPE_COMMIT(TYPE_LAPX1,ierr)
   call MPI_TYPE_COMMIT(TYPE_LAPY1,ierr)
   call MPI_TYPE_COMMIT(TYPE_LAPZ1,ierr)

   call MPI_TYPE_COMMIT(TYPE_LAPX2,ierr)
   call MPI_TYPE_COMMIT(TYPE_LAPY2,ierr)
   call MPI_TYPE_COMMIT(TYPE_LAPZ2,ierr)

   call MPI_barrier(MPI_COMM_WORLD,ierr)
   end
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
       subroutine get_k_node(k_global,node_k,k_local)
       Use OCFD_constants
       implicit none
       integer k_global,node_k,k_local,ka

         node_k=npz0-1
       do ka=0,npz0-2
        if(k_global .ge. k_offset(ka) .and. k_global .lt. k_offset(ka+1) ) node_k=ka
       enddo
         k_local=k_global-k_offset(node_k)+1
       end

!------------------------------------------------------------------------------------
       function get_id(npx1,npy1,npz1)
       Use OCFD_constants
       implicit none
       integer get_id,npx1,npy1,npz1
       get_id=npz1*(npx0*npy0)+npy1*npx0+npx1
       return
       end
!=========================================================================================================
!  Boundary message communication (exchange message)  
!=========================================================================================================
!  Standard (most used)
       subroutine exchange_boundary_x_standard(f,Iperiodic1)
       Use OCFD_flow_parameters
       implicit none
       integer:: Iperiodic1
!      send and recv mesg, to exchange_boundary array in x direction.
!      To avoid msg block, cutting long msg to short msgs   
       integer Status(MPI_status_Size),ierr
       integer i,j,k,k1,nsize
       real(kind=OCFD_REAL_KIND):: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP),  &
            tmp_send1(LAP*ny*nz),tmp_send2(LAP*ny*nz),  &
            tmp_recv1(LAP*ny*nz),tmp_recv2(LAP*ny*nz)
            nsize=LAP*ny*nz
       do k=1,nz
       do j=1,ny
       do i=1,LAP
       k1=(k-1)*ny*LAP+(j-1)*LAP+i
       tmp_send1(k1)=f(i,j,k)
       tmp_send2(k1)=f(nx-LAP+i,j,k)
       enddo
       enddo
       enddo
   call MPI_Sendrecv(tmp_send1(1),nsize,  OCFD_DATA_TYPE, ID_XM1, 9000,  &
         tmp_recv2(1),nsize, OCFD_DATA_TYPE, ID_XP1, 9000,MPI_COMM_WORLD,Status,ierr)

   call MPI_Sendrecv(tmp_send2(1),nsize, OCFD_DATA_TYPE, ID_XP1, 8000, &
        tmp_recv1(1),nsize,  OCFD_DATA_TYPE,  ID_XM1, 8000,MPI_COMM_WORLD,Status,ierr)
     
! if not periodic, node npx=0 Do Not need f(i-LAP,j,k) 
 if(npx.ne.0 .or. Iperiodic1 .eq. 1) then    
      do k=1,nz
      do j=1,ny
      do i=1,LAP
          k1=(k-1)*ny*LAP+(j-1)*LAP+i
          f(i-LAP,j,k)=tmp_recv1(k1)
      enddo
      enddo
      enddo
endif

  if(npx.ne.npx0-1 .or. Iperiodic1 .eq. 1) then    
      do k=1,nz
      do j=1,ny
      do i=1,LAP
          k1=(k-1)*ny*LAP+(j-1)*LAP+i
          f(nx+i,j,k)=tmp_recv2(k1)
      enddo
      enddo
      enddo
  endif
        
  end
!------------------------------------------------------
       subroutine exchange_boundary_y_standard(f,Iperiodic1)
       Use OCFD_flow_parameters
       implicit none
       integer:: Iperiodic1
       integer Status(MPI_status_Size),ierr
       integer i,j,k,k1,nsize
       real(kind=OCFD_REAL_KIND):: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP),   &
          tmp_send1(LAP*nx*nz),tmp_send2(LAP*nx*nz),  &
          tmp_recv1(LAP*nx*nz),tmp_recv2(LAP*nx*nz)

     
       nsize=LAP*nx*nz
       do k=1,nz
       do j=1,LAP
       do i=1,nx
       k1=(k-1)*(LAP*nx)+(j-1)*nx+i
       tmp_send1(k1)=f(i,j,k)
       tmp_send2(k1)=f(i,ny+j-LAP,k)
       enddo
       enddo
       enddo
    call MPI_Sendrecv(tmp_send1(1),nsize,  OCFD_DATA_TYPE, ID_YM1, 9000,  &
          tmp_recv2(1),nsize,   OCFD_DATA_TYPE,  ID_YP1,   9000,MPI_COMM_WORLD,Status,ierr)
    call MPI_Sendrecv(tmp_send2(1),nsize, OCFD_DATA_TYPE,  ID_YP1,  8000,  &
          tmp_recv1(1),nsize, OCFD_DATA_TYPE,  ID_YM1, 8000,MPI_COMM_WORLD,Status,ierr)

  if(npy.ne.0 .or. Iperiodic1 .eq. 1) then    
          do k=1,nz
          do j=1,LAP
          do i=1,nx
          k1=(k-1)*(LAP*nx)+(j-1)*nx+i
          f(i,j-LAP,k)=tmp_recv1(k1)
          enddo
          enddo
          enddo
  endif

  if(npy.ne.npy0-1 .or. Iperiodic1 .eq. 1) then    
          do k=1,nz
          do j=1,LAP
          do i=1,nx
          k1=(k-1)*(LAP*nx)+(j-1)*nx+i
          f(i,ny+j,k)=tmp_recv2(k1)
          enddo
          enddo
          enddo
   endif

  end
!------------------------------------------------------------
       subroutine exchange_boundary_z_standard(f, Iperiodic1)
       Use OCFD_flow_parameters
       implicit none
       integer:: Iperiodic1
       integer Status(MPI_status_Size),ierr
       integer i,j,k,k1,nsize
       real(kind=OCFD_REAL_KIND):: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP),   &
           tmp_send1(LAP*nx*ny),tmp_send2(LAP*nx*ny),   &
           tmp_recv1(LAP*nx*ny),tmp_recv2(LAP*nx*ny)

!-------------------------------------------------------------------------------     
     
       nsize=LAP*nx*ny
       do k=1,LAP
       do j=1,ny
       do i=1,nx
       k1=(k-1)*nx*ny+(j-1)*nx+i
       tmp_send1(k1)=f(i,j,k)
       tmp_send2(k1)=f(i,j,nz-LAP+k)
       enddo
       enddo
       enddo


   call MPI_Sendrecv(tmp_send1(1),nsize,   OCFD_DATA_TYPE,   ID_ZM1,  9000,   &
        tmp_recv2(1),nsize,   OCFD_DATA_TYPE,  ID_ZP1,  9000,MPI_COMM_WORLD,Status,ierr)
   call MPI_Sendrecv(tmp_send2(1),nsize,   OCFD_DATA_TYPE,  ID_ZP1, 8000,  &
        tmp_recv1(1),nsize,  OCFD_DATA_TYPE,  ID_ZM1,   8000,MPI_COMM_WORLD,Status,ierr)
 
   if(npz.ne. 0 .or. Iperiodic1 .eq. 1) then    
          do k=1,LAP
          do j=1,ny
          do i=1,nx
          k1=(k-1)*nx*ny+(j-1)*nx+i
          f(i,j,k-LAP)=tmp_recv1(k1)
          enddo
          enddo
          enddo
    endif

   if(npz.ne.npz0-1 .or. Iperiodic1 .eq. 1) then    
          do k=1,LAP
          do j=1,ny
          do i=1,nx
          k1=(k-1)*nx*ny+(j-1)*nx+i
          f(i,j,nz+k)=tmp_recv2(k1)
          enddo
          enddo
          enddo
   endif

  end

    subroutine OCFD_save(Istep,tt,Ksave_style,IFLAG_G,   &
          u,nx,ny,nz,LAP,nx_global,ny_global,nz_global,Iflag_av,Istep_name,dt,end_time,Kstep_save)     
! Iflag_av==0, write opencfd file; ==1, write averaged data file
      use OCFD_constants
      implicit none
      integer Istep,nx,ny,nz,LAP,nx_global,ny_global,nz_global,IFLAG_G
      integer i,j,k,ierr,Ksave_style,Iflag_av,Istep_name
      real(kind=OCFD_REAL_KIND):: tt
      real(kind=OCFD_REAL_KIND)::dt,end_time
      integer Kstep_save
      real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)::u

        character(len=120) filename1
        character(len=7) Cnum
!-------------------------------------------      
      write(Cnum,'(I7.7)') Istep_name
     
      filename1='OCFD'//Cnum//'.dat'

      if(my_id.eq.0) then
        open(155,file=filename1,form='unformatted')
        write(155) nx_global,ny_global,nz_global
        write(155) dt,end_time,Kstep_save
        write(155) Istep,tt
      endif
      call write_3d1(155,u,nx,ny,nz,LAP,nx_global,ny_global,nz_global)
       
       if(my_id.eq.0 )  then
         close(155)
       endif
     end


!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
       subroutine write_3d1(file_no,U,nx,ny,nz,LAP,nx_global,ny_global,nz_global)
        use OCFD_constants
        implicit none
           integer nx,ny,nz,nx_global,ny_global,nz_global,LAP,file_no,i,j,k,ierr
       real(kind=OCFD_REAL_KIND):: U(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
           real(kind=OCFD_REAL_KIND),allocatable:: U1(:,:,:)
           allocate(U1(nx,ny,nz))
       do k=1,nz
       do j=1,ny
       do i=1,nx
       U1(i,j,k)=U(i,j,k)
       enddo
       enddo
       enddo
       call write_3d(file_no,U1,nx,ny,nz,nx_global,ny_global,nz_global)
           deallocate(U1)
       end

!-----------------------------------------------------------
     subroutine write_3d(file_no,U,nx,ny,nz,nx_global,ny_global,nz_global)
      use OCFD_constants
      implicit none
       integer nx,ny,nz,nx_global,ny_global,nz_global,file_no,ierr
       real(kind=OCFD_REAL_KIND):: U(nx,ny,nz)
 
       call write_3d_message_gather(file_no,U,nx,ny,nz,nx_global,ny_global,nz_global)  

     end  

!------------------------------------------------------------------------------------------
     subroutine write_3d_message_gather(file_no,U,nx,ny,nz,nx0,ny0,nz0)
      use OCFD_constants
      implicit none
         integer nx,ny,nz,nx0,ny0,nz0
         integer Status(MPI_status_Size),ierr,i1,j1,k1,nk,i0,j0,ia,i,j,k,kk
         integer file_no,npx1,npy1,npz1,recvcount
     real(kind=OCFD_REAL_KIND):: U(nx,ny,nz) 
     real(kind=OCFD_REAL_KIND),allocatable:: buff2d(:,:),buff1(:,:),buff2(:),buff_send(:)
     integer,allocatable:: recvcounts1(:), displs1(:),recvcounts2(:), displs2(:) 
!---------------------------------------------------------------
       allocate(buff2d(nx0,ny0),buff1(nx0,ny),buff2(nx0*ny),buff_send(nx*ny))
       allocate(recvcounts1(npy0),displs1(npy0),recvcounts2(npx0),displs2(npx0))

    if(my_id.eq.0)    print*, 'write 3d data ...'

           do npy1=0,npy0-1
              recvcounts1(npy1+1)=nx0*j_nn(npy1)
              if(npy1 .eq. 0) then
                   displs1(npy1+1)=0
              else
                   displs1(npy1+1)=displs1(npy1)+recvcounts1(npy1)
              endif
            enddo


            do npx1=0, npx0-1
              recvcounts2(npx1+1)=i_nn(npx1)*ny
              if(npx1 .eq. 0) then
                   displs2(npx1+1)=0
              else
                   displs2(npx1+1)=displs2(npx1)+recvcounts2(npx1)
              endif
           enddo

    do kk=1,nz0
      call get_k_node(kk,nk,k1)
 
      if(npz .eq. nk ) then
             do j1=1,ny
             do i1=1,nx
               buff_send(i1+(j1-1)*nx)=U(i1,j1,k1)
             enddo
             enddo
               call MPI_gatherv(buff_send,nx*ny,OCFD_DATA_TYPE, buff2, recvcounts2,displs2,OCFD_DATA_TYPE, 0,MPI_COMM_X,ierr)

                 ia=0
                do npx1=0,npx0-1
                 do i=1,i_nn(npx1)*ny
                  i0=i_offset(npx1)+mod(i-1,i_nn(npx1))     ! i_offset(npx1)+(mod(i-1,i_nn(npx1))+1) -1
                  j0=int((i-1)/i_nn(npx1))+1
                  buff1(i0,j0)=buff2(ia+i)
                enddo
                 ia=ia+i_nn(npx1)*ny
                enddo
           
                 if(npx .eq. 0) then
                      call MPI_gatherv(buff1,nx0*ny,OCFD_DATA_TYPE,buff2d,recvcounts1,displs1,OCFD_DATA_TYPE,0,MPI_COMM_Y,ierr)
                 endif
         endif

            if(nk .ne. 0 ) then
               if(npx.eq.0 .and. npy .eq. 0 .and. npz .eq. nk) call   MPI_send(buff2d, nx0*ny0,  OCFD_DATA_TYPE, 0, 6666, MPI_COMM_WORLD,ierr)
               if(my_id .eq. 0)  call   MPI_recv(buff2d, nx0*ny0,  OCFD_DATA_TYPE, nk*(npx0*npy0), 6666, MPI_COMM_WORLD,status,ierr)
            endif

            
             if(my_id.eq.0) then
                       write(file_no) buff2d
             endif
      enddo

           deallocate(buff2d,buff1,buff2,buff_send)
       deallocate(recvcounts1,displs1,recvcounts2,displs2)
   end
      subroutine read_parameters
      Use OCFD_flow_parameters
      implicit none
      integer i,j,k,nk,nr,ierr, ntmp(100)
      real(kind=OCFD_REAL_KIND)::   rtmp(100)
!-------------------------------------------
        ntmp=0
      if(my_id.eq.0) then
      open(30,file='weno.in')
      read(30,*)
      read(30,*)
      read(30,*) nx_global,ny_global,nz_global
      read(30,*)
      read(30,*) npx0,npy0,npz0,LAP
      read(30,*)
      read(30,*) Iflag_grid(1),Iflag_grid(2),Iflag_grid(3),Iperiodic(1),Iperiodic(2),Iperiodic(3)
      read(30,*)
      read(30,*) SLx,SLy,SLz
      read(30,*)
      read(30,*) dt,end_time,NUM_TIME_Advance,Kstep_show, Kstep_save,Kread_style,Ksave_style
      close(30)
      endif
   
!     Boardcast integer and real parameters to all proc      
      ntmp(1)=nx_global
      ntmp(2)=ny_global
      ntmp(3)=nz_global
      ntmp(4)=npx0
      ntmp(5)=npy0
      ntmp(6)=npz0
      ntmp(7)=LAP 
      ntmp(8)=Kstep_show
      ntmp(9)=Kstep_save
      ntmp(10)=Ksave_style
      ntmp(11)=Kread_style

     call MPI_bcast(ntmp(1),100,MPI_INTEGER,0, MPI_COMM_WORLD,ierr)

     nx_global= ntmp(1)
     ny_global= ntmp(2)
     nz_global= ntmp(3)
     npx0 =ntmp(4)
     npy0 =ntmp(5)
     npz0 =ntmp(6)
     LAP =ntmp(7)    
     Kstep_show=ntmp(8)
     Kstep_save = ntmp(9)
     Ksave_style= ntmp(10)
     Kread_style= ntmp(11)
        IFLAG_G=0

!c----------------------------------------------------
      rtmp(1)=SLx
      rtmp(2)=SLy
      rtmp(3)=SLz
      rtmp(4)=dt
      rtmp(5)=end_time

     call MPI_bcast(rtmp(1),100,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)
      SLx=rtmp(1)
      SLy=rtmp(2)
      SLz=rtmp(3)
      dt=rtmp(4)
      end_time=rtmp(5)
! ----------------------------------
!c----------------------------------------------------------------------------

    call MPI_bcast(Iflag_grid(1),3,MPI_INTEGER,0,  MPI_COMM_WORLD,ierr)
    call MPI_bcast(Iperiodic(1),3,MPI_INTEGER,0,  MPI_COMM_WORLD,ierr)

!---------------------------------------------------------------------------------
  end
