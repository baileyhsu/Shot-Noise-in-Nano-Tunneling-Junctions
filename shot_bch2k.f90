PROGRAM SHOTMPI
USE glob_var

!===========================================================
!THIS CODE CALCULATES THE SHOTNOISE IN A MOLECULAR JUNCTION!
!===========================================================

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       include 'mpif.h'
       REAL*8 PI,K2MXB2,LL2(2),AAA
       INTEGER, DIMENSION(:), allocatable :: NFILE,NFILE1
       REAL*8, DIMENSION(:), allocatable ::EE,WGTE,NLOC,ZV
       REAL*8, DIMENSION(:), allocatable ::PHI,WGTPHI,KAPPAX,WGTKPX,XX
       REAL*8, DIMENSION(:,:), allocatable ::WGTK
       REAL*8, DIMENSION(:,:,:), allocatable ::RJ
       REAL*8 DXL
       REAL*8 EPAN(50,3),NORDPN(50,3),NPAN(3)
       DIMENSION NIE(0:NEMAX),IELIST(0:NEMAX-1,NEMAX)
       COMPLEX*16, DIMENSION(:,:,:,:,:,:), allocatable :: PSI0_R,PSI_R
       COMPLEX*16, DIMENSION(:,:,:,:,:,:), allocatable :: PSI0_L,PSI_L
       INTEGER NZV,N_X,N_Y,N_Z,NE,NPHI,NKAPAX,Nup,Ndn,NUNIT,FFT_MAX_Z,FFT_Z
       INTEGER FFT_MAX_XY
       REAL*8  CL,lpar1,lpar2,L,Lav,Ef_R,Ef_L,Eb_R,Eb_L,Eup,Edn    
       CHARACTER(len=50) path,path1,path2,path3,path4,path5
       NAMELIST/PARAM/NZV,lpar1,lpar2,L,ZLEFT,Edn,Ndn,Eup,Nup,&
                      N_X,N_Y,N_Z,NPHI,NKAPAX,NE1,NE2,NE3,Ef_L,&
                      Ef_R,Eb_L,Eb_R,FFT_MAX_XY,FFT_MAX_Z
       NAMELIST/ATOM/RJ
       NAMELIST/PANELS/NPAN,NORDPN,EPAN
       INTEGER NE3
       REAL*8, allocatable, DIMENSION(:,:) :: ZVAL,WGTZ
       REAL*8, allocatable, DIMENSION(:) :: ZTEMP,WZTEMP,LVX,LVY,LVZ,W
       COMPLEX*16, allocatable, DIMENSION(:) ::CIKZ  
       COMPLEX*16,ALLOCATABLE,DIMENSION(:,:,:) :: CTMP
       INTEGER LX,LY,LZ,ISZ
!====================================================

      REAL*8 pGIZT,pGIZ0,pGIZD,GIZT,GIZ0,GIZD
      REAL*8 pGSHTZT,pGSHTZ0,pGSHTZD
      REAL*8 GSHTZT,GSHTZ0,GSHTZD
      REAL*8 pGIZET, pGIZE0, pGIZED
      REAL*8 TMP1,TMP2
      COMPLEX*16, DIMENSION(:,:,:,:,:), allocatable:: SILRT,SILR0,SILRD,SIRRT,SIRR0,SIRRD
      REAL*8,DIMENSION(:),allocatable :: pcnt_t,pcnt_0,pcnt_d,cnt_t,cnt_0,cnt_d,cntx_t,cntx_0,cntx_d
      REAL*8,DIMENSION(:),allocatable ::psht_t,psht_0,psht_d,sht_t,sht_0,sht_d,shtx_t,shtx_0,shtx_d
!END SHOT
!==========================================================
      INTEGER status(MPI_STATUS_SIZE)
allocate(pcnt_t(1024),pcnt_0(1024),pcnt_d(1024),cnt_t(1024),cnt_0(1024),cnt_d(1024),cntx_t(3072),cntx_0(3072),cntx_d(3072),&
psht_t(1024),psht_0(1024),psht_d(1024),sht_t(1024),sht_0(1024),sht_d(1024),shtx_t(3072),shtx_0(3072),shtx_d(3072),STAT=ialloc)
if(ialloc.ne.0)stop 'Unable to allocate memory for SHOT3'

!=======================================================
      PI=ACOS(-1.D0)       
      call MPI_INIT(IERR)
      call MPI_Comm_size(MPI_COMM_WORLD,nNodes,IERR)
      call MPI_Comm_rank(MPI_COMM_WORLD,myNode,IERR)
!=====READING OLD PARAM.TXT==========================
      path4=trim('./wave/')//trim('param.txt')
      IF(mynode.eq.0) THEN
      open(993,file='SZ_bch2K.out')
      open(16,file='print_bch2K')
      END IF
      OPEN(100,file=path4)
      READ(100,PARAM)
      CL=4.D0*lpar1*lpar2
      NE=NE1+NE2+NE3
      allocate(EE(NE),WGTE(NE),NLOC(NZV),RJ(3,20,NZV),ZV(NZV),&
      PHI(NPHI),WGTPHI(NPHI),KAPPAX(NKAPAX),WGTKPX(NKAPAX),STAT=ialloc)
      if(ialloc.ne.0)stop 'Unable to allocate memory for COORD'
      READ(100,ATOM)
      READ(100,PANELS)
       DO I=1,NE
         READ(100,*)EE(I),WGTE(I)
       END DO
      CLOSE(100)

!============================
! WGTZ weight cal
!============================
     FFT_Z=FFT_MAX_Z
  
     allocate(ZVAL(NPTZ,2),WGTZ(NPTZ,2),ZTEMP(NPTZ),WZTEMP(NPTZ),XX(FFT_Z),  &
              STAT=ialloc)
     if(ialloc.ne.0)stop 'unable to allocate memory for dwave'

     allocate(LVX(-N_X:N_X),LVY(-N_Y:N_Y),LVZ(-N_Z:N_Z),CIKZ(-N_Z:N_Z),  &
              STAT=ialloc)
     if(ialloc.ne.0)stop 'unable to allocate memory for dwave'


         Z1=-L
         Z2=L
         Z3=-L+L*(DINT(N_Z/2.d0)+1.d0)/dfloat(N_Z+1)
         Z4=L-L*(DINT(N_Z/2.d0)+1.d0)/dfloat(N_Z+1)
           LL2(1)=Z2-Z1
           LL2(2)=Z4-Z3

         CALL GAUSSLEG(Z1,Z2,NPTZ,ZTEMP,WZTEMP)
            DO ISZ=1,NPTZ
              ZVAL(ISZ,1)=ZTEMP(ISZ)
              WGTZ(ISZ,1)=WZTEMP(ISZ)
            END DO
!
         CALL GAUSSLEG(Z3,Z4,NPTZ,ZTEMP,WZTEMP)
            DO ISZ=1,NPTZ
              ZVAL(ISZ,2)=ZTEMP(ISZ)
              WGTZ(ISZ,2)=WZTEMP(ISZ)
            End Do

      DXL=2.D0*L/DFLOAT(FFT_Z-1)
      DO I=1,FFT_Z
         XX(I)=-L+DFLOAT(I-1)*DXL
      End Do

!======================================================
!LACZOS CONVERGENCE FACTOR
!======================================================

         DO LX=-N_X,N_X
            CVX=DFLOAT(LX)*PI/DFLOAT(N_X+1)
            IF(LX.EQ.0) THEN
               LVX(LX)=1.d0
            ELSE
               LVX(LX)=DSIN(CVX)/CVX
            END IF
             END DO
!
         DO LY=-N_Y,N_Y
            CVY=DFLOAT(LY)*PI/DFLOAT(N_Y+1)
            IF(LY.EQ.0) THEN
               LVY(LY)=1.d0
            ELSE
            LVY(LY)=DSIN(CVY)/CVY
            END IF
         END DO
!
          DO LZ=-N_Z,N_Z
            CVZ=DFLOAT(LZ)*PI/DFLOAT(N_Z+1)
              CIKZ(LZ)=(0.D0,1.D0)*DFLOAT(LZ)*PI/L
            IF(LZ.EQ.0) THEN
               LVZ(LZ)=1.d0
            ELSE
              LVZ(LZ)=DSIN(CVZ)/CVZ
            END IF
        END DO

!================================================
      NES=NE-NE3+1-Ndn
      
      allocate(PSI0_R(-N_X:N_X,-N_Y:N_Y,-N_Z:N_Z,1:NPHI,1:NKAPAX,NES:NE),STAT=ialloc)
      allocate(PSI_R(-N_X:N_X,-N_Y:N_Y,-N_Z:N_Z,1:NPHI,1:NKAPAX,NES:NE),STAT=ialloc)     
     allocate(PSI0_L(-N_X:N_X,-N_Y:N_Y,-N_Z:N_Z,1:NPHI,1:NKAPAX,NES:NE))
      allocate(PSI_L(-N_X:N_X,-N_Y:N_Y,-N_Z:N_Z,1:NPHI,1:NKAPAX,NES:NE))


      allocate(NFILE1(NES:NE),STAT=ialloc)
      allocate(NFILE(NES:NE),STAT=ialloc)
      if(ialloc.ne.0)stop 'Unable to allocate memory for wavef unctions'

         DO IE=NES,NE 
          NFILE(IE)=700+IE
          NFILE1(IE)=300+IE
          IECOUNT=INT(DABS(EE(IE))*1.d5)
          write(path,*) IECOUNT
          path1=ADJUSTL(path)
          path3=trim('./wave/')//trim(path1)//trim('r.bin')
          path5=trim('./wave/')//trim(path1)//trim('l.bin')
          

          open(unit=NFILE(IE),file=path3,form="unformatted",  &
               access="sequential",status="old")
          open(unit=NFILE1(IE),file=path5,form="unformatted",  &
               access="sequential",status="old")

           DO IKAPAX=1,NKAPAX
            DO IPHI=1,NPHI
             DO LZ=-N_Z,N_Z
              DO LY=-N_Y,N_Y
                DO LX=-N_X,N_X
                 READ(unit=NFILE(IE))PSI0_R(LX,LY,LZ,IPHI,IKAPAX,IE), & 
                       PSI_R(LX,LY,LZ,IPHI,IKAPAX,IE) 
                 READ(unit=NFILE1(IE))PSI0_L(LX,LY,LZ,IPHI,IKAPAX,IE), &
                       PSI_L(LX,LY,LZ,IPHI,IKAPAX,IE)
                End Do
               End DO
             End DO
            End Do
           END DO
          END DO 


!============

     pGIZT=0.d0
     pGIZ0=0.d0
     pGIZD=0.d0

     pGIZET=0.d0
     pGIZE0=0.d0
     pGIZED=0.d0
     pGSHTZT=0.d0
     pGSHTZ0=0.d0
     pGSHTZD=0.d0
     
       pcnt_0=0.d0
       pcnt_d=0.d0
       pcnt_t=0.d0


!=========SHOT NOISE CALCULATION=====================
         allocate(SILRT(1:NPHI,1:NKAPAX,1:NPHI,1:NKAPAX,NES:NE))
         allocate(SILR0(1:NPHI,1:NKAPAX,1:NPHI,1:NKAPAX,NES:NE))
         allocate(SILRD(1:NPHI,1:NKAPAX,1:NPHI,1:NKAPAX,NES:NE))
         allocate(SIRRT(1:NPHI,1:NKAPAX,1:NPHI,1:NKAPAX,NES:NE))
         allocate(SIRR0(1:NPHI,1:NKAPAX,1:NPHI,1:NKAPAX,NES:NE))
         allocate(SIRRD(1:NPHI,1:NKAPAX,1:NPHI,1:NKAPAX,NES:NE))
         allocate(CTMP(-N_Z:N_Z,-N_Z:N_Z,1:NPTZ))
         ALLOCATE(WGTK(1:NPHI,1:NKAPAX))

call MPI_Barrier(MPI_COMM_WORLD,IERR)
!===============================================
!DISTRIBUTE OVER ENERGY PROCESSORS
!===============================================
      
      IELIST = 0
      NIE = 0

      If (MyNode == 0) then
         NIE = 0
         IELIST = 0
         Do i=NE,NES,-1
            Call MPI_RECV(j, 1, MPI_INTEGER, MPI_ANY_SOURCE, 100, &
                 MPI_COMM_WORLD, status, ierr)

            j = status(MPI_SOURCE)

            NIE(j) = NIE(j) + 1
            IELIST(j,NIE(j)) = i

            Write(16,*) 'Assigning ',i,' to ',j
            Call Flush(16)
            Call MPI_Send(i, 1, MPI_INTEGER, j, 100, MPI_COMM_WORLD, ierr)
         End Do

         !*** Now wait for everybody to finish ***
         Do i=nNodes-2,0,-1
            Call MPI_RECV(j, 1, MPI_INTEGER, MPI_ANY_SOURCE, 100, &
                 MPI_COMM_WORLD, status, ierr)

            j = status(MPI_SOURCE)

            Write(16,*) 'Stopping ',j, ' * left:',i
            Call FLush(16)
            k = -1
            Call MPI_Send(k, 1, MPI_INTEGER, j, 100, MPI_COMM_WORLD, ierr)
         End Do

         IE = -1
      else
         i = 1
         Call MPI_Send(i, 1, MPI_INTEGER, 0, 100, MPI_COMM_WORLD, ierr)
         Call MPI_RECV(IE, 1, MPI_INTEGER, 0, 100, MPI_COMM_WORLD, status, ierr)
      End If


      JJJ =0
      DO WHILE (IE>0)
        JJJ= JJJ + 1
         E=EE(IE)

!============================ 
       CALL GAUSSLEG(0.d0,2.d0*PI,NPHI,PHI,WGTPHI)
       K2MXB2 = E-Eb_R
     
       CALL GAUSSLEG(0.d0,K2MXB2,NKAPAX,KAPPAX,WGTKPX)
       KAPPAX(:)=DSQRT(KAPPAX(:))
      
       DO IPHI=1,NPHI
        DO IKAPAX=1,NKAPAX
         WGTK(IPHI,IKAPAX)=WGTPHI(IPHI)*WGTKPX(IKAPAX)
        END DO
       END DO

        SILRT=(0.d0,0.d0)
        SILR0=(0.d0,0.d0)
        SILRD=(0.d0,0.d0)
        SIRRT=(0.d0,0.d0)
        SIRR0=(0.d0,0.d0)
        SIRRD=(0.d0,0.d0)

DO ISZ=1,NPTZ
 DO LZ=-N_Z,N_Z
  DO LZP=-N_Z,N_Z
CTMP(LZ,LZP,ISZ)=(CIKZ(LZ)+CIKZ(LZP))*(LVZ(LZ)*LVZ(LZP))*CDEXP((CIKZ(LZ)-CIKZ(LZP))*ZVAL(ISZ,1))*WGTZ(ISZ,1)/LL2(1)   
  END DO
 END DO
END DO


 DO ISZ=1,NPTZ
  DO LZ=-N_Z,N_Z
   DO LZP=-N_Z,N_Z
     DO IPHI1=1,NPHI
      DO IKAPAX1=1,NKAPAX
       DO IPHI2=1,NPHI
        DO IKAPAX2=1,NKAPAX
         DO LY=-N_Y,N_Y
          DO LX=-N_X,N_X

SIRRT(IPHI1,IKAPAX1,IPHI2,IKAPAX2,IE)=SIRRT(IPHI1,IKAPAX1,IPHI2,IKAPAX2,IE)+DCONJG(PSI_R(LX,LY,LZ,IPHI1,IKAPAX1,IE))*PSI_R(LX,LY,LZP,IPHI2,IKAPAX2,IE)*CTMP(LZ,LZP,ISZ)

SIRR0(IPHI1,IKAPAX1,IPHI2,IKAPAX2,IE)=SIRR0(IPHI1,IKAPAX1,IPHI2,IKAPAX2,IE)+DCONJG(PSI0_R(LX,LY,LZ,IPHI1,IKAPAX1,IE))*PSI0_R(LX,LY,LZP,IPHI2,IKAPAX2,IE)*CTMP(LZ,LZP,ISZ)


SILRT(IPHI1,IKAPAX1,IPHI2,IKAPAX2,IE)=SILRT(IPHI1,IKAPAX1,IPHI2,IKAPAX2,IE)+DCONJG(PSI_L(LX,LY,LZ,IPHI1,IKAPAX1,IE))*PSI_R(LX,LY,LZP,IPHI2,IKAPAX2,IE)*CTMP(LZ,LZP,ISZ)

SILR0(IPHI1,IKAPAX1,IPHI2,IKAPAX2,IE)=SILR0(IPHI1,IKAPAX1,IPHI2,IKAPAX2,IE)+DCONJG(PSI0_L(LX,LY,LZ,IPHI1,IKAPAX1,IE))*PSI0_R(LX,LY,LZP,IPHI2,IKAPAX2,IE)*CTMP(LZ,LZP,ISZ)

               END DO
              END DO
            END DO
       END DO
     END DO 
    END DO   
  END DO 
 END DO
END DO
       SIRRD=SIRRT-SIRR0
       SILRD=SILRT-SILR0

           pGIZET=0.d0
           pGIZE0=0.d0
           pGIZED=0.d0
           pGSHTET=0.d0
           pGSHTE0=0.d0
           pGSHTED=0.d0


      DO IPHI1=1,NPHI
       DO IKAPAX1=1,NKAPAX
           TMP1=0.5d0*WGTK(IPHI1,IKAPAX1)*CL
!CURRENT:
pGIZET=pGIZET-DREAL((0.d0,-1.d0)*SIRRT(IPHI1,IKAPAX1,IPHI1,IKAPAX1,IE))*TMP1*WGTE(IE)
pGIZE0=pGIZE0-DREAL((0.d0,-1.d0)*SIRR0(IPHI1,IKAPAX1,IPHI1,IKAPAX1,IE))*TMP1*WGTE(IE)
pGIZED=pGIZED-DREAL((0.d0,-1.d0)*SIRRD(IPHI1,IKAPAX1,IPHI1,IKAPAX1,IE))*TMP1*WGTE(IE)

DO IPHI2=1,NPHI
 DO IKAPAX2=1,NKAPAX

!SHOTNOISE:
TMP2=2.d0*Pi*WGTK(IPHI1,IKAPAX1)*CL*CL*0.5D0*WGTK(IPHI2,IKAPAX2)

           pGSHTZT=pGSHTZT+DREAL(DCONJG(SILRT(IPHI1,IKAPAX1,IPHI2,IKAPAX2,IE))*SILRT(IPHI1,IKAPAX1,IPHI2,IKAPAX2,IE))*TMP2*WGTE(IE)
           pGSHTZ0=pGSHTZ0+DREAL(DCONJG(SILR0(IPHI1,IKAPAX1,IPHI2,IKAPAX2,IE))*SILR0(IPHI1,IKAPAX1,IPHI2,IKAPAX2,IE))*TMP2*WGTE(IE)
           pGSHTZD=pGSHTZD+DREAL(DCONJG(SILRD(IPHI1,IKAPAX1,IPHI2,IKAPAX2,IE))*SILRD(IPHI1,IKAPAX1,IPHI2,IKAPAX2,IE))*TMP2*WGTE(IE)

       END DO
      END DO
   END DO
 END DO
           
!===================================================
         i = 1
         Call MPI_Send(i, 1, MPI_INTEGER, 0, 100, MPI_COMM_WORLD, ierr)
         Call MPI_RECV(IE, 1, MPI_INTEGER, 0, 100, MPI_COMM_WORLD, status, ierr)
END DO !END DO IE DO LOOP
    !*** Send everybody IELIST and NIE
         Call MPI_BCAST(NIE, Size(NIE), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         Call MPI_BCAST(IELIST, Size(IELIST), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         nPerNode = MaXVal(NIE(0:nNodes-1))

         IF(myNode.EQ.0)THEN
         DO i=0,nNodes -1
            WRITE(16,4332)I,NIE(I),(IELIST(I,J),J=1,NIE(I))
4332        FORMAT(' NODE',I3,' NIE=',I2,'  IE=',(20I3))
         End Do
         END IF

!=========================================================
! End the distributed energy loop, communicate distributed data
!=========================================================

      Call MPI_Barrier(MPI_COMM_WORLD,IERR)

        icount=1
       call MPI_Allreduce(pGIZET,GIZT,icount,MPI_REAL8,mpi_sum,&
           MPI_COMM_WORLD,IERR)
       call MPI_Allreduce(pGIZE0,GIZ0,icount,MPI_REAL8,mpi_sum,&
           MPI_COMM_WORLD,IERR)
       call MPI_Allreduce(pGIZED,GIZD,icount,MPI_REAL8,mpi_sum,&
           MPI_COMM_WORLD,IERR)

        icount=1
       call MPI_Allreduce(pGSHTZT,GSHTZT,icount,MPI_REAL8,mpi_sum,&
           MPI_COMM_WORLD,IERR)
       call MPI_Allreduce(pGSHTZ0,GSHTZ0,icount,MPI_REAL8,mpi_sum,&
           MPI_COMM_WORLD,IERR)
       call MPI_Allreduce(pGSHTZD,GSHTZD,icount,MPI_REAL8,mpi_sum,&
           MPI_COMM_WORLD,IERR)

!begin{SHOTNOISE}{WRITE SHOTNOISE INFO TO THE FILE}

       IF(MYNODE.EQ.0) THEN
          WRITE(993,*) "TOTAL SHOT, BARE SHOT, DIFF SHOT"
          WRITE(993,174) GSHTZT,GSHTZ0,GSHTZD
          WRITE(993,*) "TOTAL CURRENT, BARE CURRENT, DIFF CURRENT"
          WRITE(993,174) GIZT,GIZ0,GIZD       
          WRITE(993,*) "TOTAL FANO, METAL FANO, NET FANO"
          WRITE(993,174) GSHTZT/2.D0/GIZT, GSHTZ0/2.D0/GIZ0, GSHTZD/2.D0/GIZD
       END IF
 
174     FORMAT(3E16.7)

        call MPI_Barrier(MPI_COMM_WORLD,IERR)
        Call MPI_FINALIZE(IERR)

END PROGRAM
