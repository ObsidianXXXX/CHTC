        SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1  PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2  KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3  LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
         INCLUDE 'ABA_PARAM.INC'
C
        DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1  SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2  DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3  JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4  PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

C       Initialize Variables
         ! Define Gauss Point of CPS4 and Shape Function Values at Gauss Point
         real*8:: w(4), gpt, s(4), t(4), N(4,4), dNs(4,4), dNt(4,4)
         ! Define Coordinate Vector
         real*8:: xEle(NNODE,1), yEle(NNODE,1)
         ! Define Stiffness Related Parameters
         real*8:: E(NNODE,1), Ele(1,1), nu, thickness, D(3,3), Dele(3,3)
         real*8:: DefE1, alpha, beta
         ! Define Jacobian Matrix
         real*8:: J(2,2), detJ
         ! Define Displacement-Strain Map Matrix B
         real*8:: B(3,NDOFEL)
         ! Output Variables
         ! SVARS(1:3): Strain(E11,E22,E33,E12)
         ! SVARS(1:3): Stress(S11,S22,S33,S12)
         real*8:: strainIgauss(3,1), stressIgauss(3,1)
         real*8:: Strain(16), Stress(16)
         integer:: iDof, jDof, iNode, iGauss, NELE, num
         ! Define Differentiation Parameters
         real*8:: dxds(1,1), dxdt(1,1), dyds(1,1), dydt(1,1)
         real*8:: duds(1,NDOFEL), dudt(1,NDOFEL)
         real*8:: dvds(1,NDOFEL), dvdt(1,NDOFEL)
         real*8:: N_i(1,4), dNs_i(1,4), dNt_i(1,4)
         real*8:: dudx(1,NDOFEL), dudy(1,NDOFEL)
         real*8:: dvdx(1,NDOFEL), dvdy(1,NDOFEL)
         ! Define Filename
C       ! Read Filepath from .ini
         character(len=100):: file
         character(len=120):: filepath
         character(len=200):: datpath,StrainBin,StressBin
         character(len=8):: filenum
         num = JPROPS(2)
         write(filenum,'(I0)') num
         file = 
     1    "D:\CHTC\2D Exponential\Batch\"
         filepath = trim(file)//trim(filenum)//"\"
         datpath = trim(filepath)//trim(filenum)//".ueldat"
         StrainBin = trim(filepath)//trim(filenum)//"Strain.bin"
         StressBin = trim(filepath)//trim(filenum)//"Stress.bin"
         ! Initialize the output file .ueldat
         if((time(1).eq.0.0) .and.(JELEM .eq.1 )) then
            call WriteHead(datpath)
         endif
         ! Gauss Points at Local Coordinates, Counterclockwise from LeftDown, Begining at the Third Quadrant
         w = (/1.0,1.0,1.0,1.0/)
         gpt = 1.0 / sqrt(3.0)
         s = (/-gpt, gpt, gpt, -gpt/)
         t = (/gpt, gpt, -gpt, -gpt/)
         N(1:4,1) = (1.0-s)*(1.0-t) / 4.0
         N(1:4,2) = (1.0+s)*(1.0-t) / 4.0
         N(1:4,3) = (1.0+s)*(1.0+t) / 4.0
         N(1:4,4) = (1.0-s)*(1.0+t) / 4.0
         dNs(1:4,1) = -(1.0-t) / 4.0
         dNs(1:4,2) = (1.0-t) / 4.0
         dNs(1:4,3) = (1.0+t) / 4.0
         dNs(1:4,4) = -(1.0+t) / 4.0
         dNt(1:4,1) = -(1.0-s) / 4.0
         dNt(1:4,2) = -(1.0+s) / 4.0
         dNt(1:4,3) = (1.0+s) / 4.0
         dNt(1:4,4) = (1.0-s) / 4.0
         ! Get the Coordinates Vector
         xEle(1:NNODE,1) = COORDS(1,1:NNODE)
         yEle(1:NNODE,1) = COORDS(2,1:NNODE)
         ! Get Material Matrix
         ! Material Properties
         DefE1 = 1.0
         alpha = PROPS(1)
         beta = PROPS(2)
         nu = PROPS(3)
         thickness = 1.0
         E = DefE1*exp(alpha*xEle+beta*yEle)
         D = 0.0
         D(1,1) = 1.0
         D(1,2) = nu
         D(2,1) = nu
         D(2,2) = 1.0
         D(3,3) = 0.5*(1-nu)
         D = D / (1.0-nu**2)

         ! Initialize Differentiation Parameters
         duds = 0.0
         dudt = 0.0
         dvds = 0.0
         dvdt = 0.0
         ! Loop through GaussPoint
         svarOld = SVARS(1)
         do iGauss = 1, 4
            ! Calculate the Stiffness Matrix at Gauss Point
            N_i(1,:) = N(iGauss,1:4)
            Ele = matmul(N_i,E)
            Dele = Ele(1,1)*D
            ! Calculate Displacement-Strain Map Matrix B
            dNs_i(1,:) = dNs(iGauss,1:4)
            dNt_i(1,:) = dNt(iGauss,1:4)
            dxds = matmul(dNs_i,xEle)
            dxdt = matmul(dNt_i,xEle)
            dyds = matmul(dNs_i,yEle)
            dydt = matmul(dNt_i,yEle)
            detJ = dxds(1,1)*dydt(1,1)-dyds(1,1)*dydt(1,1)
            dxds = dxds / detJ
            dxdt = dxdt / detJ
            dyds = dyds / detJ
            dydt = dydt / detJ
            ! Calculate grad of u and v in s-t coord
            duds(1,1:NDOFEL:2) = dNs_i(1,1:4)
            dudt(1,1:NDOFEL:2) = dNt_i(1,1:4)
            dvds(1,2:NDOFEL:2) = dNs_i(1,1:4)
            dvdt(1,2:NDOFEL:2) = dNt_i(1,1:4)
            ! Calculate grad of u and v in x-y coord
            dudx = dydt(1,1)*duds - dyds(1,1)*dudt
            dudy = -dxdt(1,1)*duds + dxds(1,1)*dudt
            dvdx = dydt(1,1)*dvds - dyds(1,1)*dvdt
            dvdy = -dxdt(1,1)*dvds + dxds(1,1)*dvdt
            ! Get Strain
            B(1,:) = dudx(1,:)
            B(2,:) = dvdy(1,:)
            B(3,:) = dudy(1,:) + dvdx(1,:)
            ! Get Ke
            AMATRX(1:NDOFEL,1:NDOFEL) = AMATRX(1:NDOFEL,1:NDOFEL) +
     1     w(iGauss)*detJ*thickness*matmul(transpose(B),matmul(Dele,B))
            strainIgauss(1:3,1) = matmul(B,U)
            stressIgauss(1:3,1) = matmul(Dele,matmul(B,U))
            ! Assign Strain and Stress to the SVD
            SVARS((iGauss-1)*8+1:(iGauss-1)*8+2) = strainIgauss(1:2,1)
            SVARS((iGauss-1)*8+3) = 0.0
            SVARS((iGauss-1)*8+4) = strainIgauss(3,1)
            SVARS((iGauss-1)*8+5:(iGauss-1)*8+6) = stressIgauss(1:2,1)
            SVARS((iGauss-1)*8+7) = 0.0
            SVARS((iGauss-1)*8+8) = stressIgauss(3,1)
            Strain((iGauss-1)*4+1:(iGauss-1)*4+4) =
     1          SVARS((iGauss-1)*8+1:(iGauss-1)*8+4)
            Stress((iGauss-1)*4+1:(iGauss-1)*4+4) =
     1          SVARS((iGauss-1)*8+5:(iGauss-1)*8+8)
         end do
         RHS(1:NDOFEL,1) = RHS(1:NDOFEL,1) - matmul(AMATRX,U)
         NELE = JPROPS(1)
         call WriteEle(JELEM,TIME,DTIME,SVARS,NELE,
     1        Strain,Stress,datpath,StressBin,StrainBin)
      end subroutine UEL

      subroutine WriteHead(datpath)
         character(len=200):: datpath
         open(unit=10,file=datpath,status='NEW',
     1            position='APPEND')
         write(10,*)
     1    "THE FOLLOWING TABLE IS PRINTED FOR ALL ELEMENTS AT THE
     2          INTEGRATION POINTS"
         write(10,*)
     1    "**OUTPUT FORMAT(Gauss Point Clockwise from [-1,1])"
         write(10,*) "**ELEMENT NUMBER:   NO.ELE=1"
         write(10,*) 
     1    "**iGauss   E11   E22   E33   E12   S11   S22   S33   S12"
         write(10,*) "**"
         write(10,*) "=================================="
         close(10)
         return
      end subroutine WriteHead

      subroutine WriteEle(JELEM,TIME,DTIME,SVARS,NELE,Strain,Stress,
     1                  datpath,StressBin,StrainBin)
         integer:: N_Gauss,JELEM,NELE
         real*8:: TIME(2),DTIME
         real*8:: Strain(16),Stress(16),SVARS(32)
         character(len=200):: datpath,StressBin,StrainBin
         logical:: Isend
         call FindEnd(datpath,IsEnd,TIME,DTIME,SVARS)
         if(Isend) then
            continue
         else
            open(10,file=datpath,status='OLD',position='APPEND')
            write(10,*) "**ELEMENT NUMBER:   NO.ELE=",JELEM,
     1                "Current TIME", TIME(1)
            write(10,"(A8,6A16)") "iGauss","E11","E22","E33","E12",
     1                "S11","S22","S33","S12"
            do N_Gauss = 1, 4
               write(10,"(I8,6E16.6)") N_Gauss,
     1                    Strain((N_Gauss-1)*4+1),
     2                    Strain((N_Gauss-1)*4+2),
     3                    Strain((N_Gauss-1)*4+3),
     4                    Strain((N_Gauss-1)*4+4),
     5                    Stress((N_Gauss-1)*4+1),
     6                    Stress((N_Gauss-1)*4+2),
     7                    Stress((N_Gauss-1)*4+3),
     8                    Strain((N_Gauss-1)*4+4)
            end do
            write(10,*) "**"
            write(10,*) "**"
            if(JELEM.eq.NELE) then
               write(10,*) "END_AT_TIME=",TIME(1)
               write(10,*) "=================================="
            endif
            close(10)
            if((time(1).eq.0.0).and.(JELEM.eq.1)) then
               ! Create & Open Binary File
               open(unit=20,file=StressBin,status='NEW',
     1              access='STREAM',form='UNFORMATTED',
     2              position='APPEND')
               open(unit=30,file=StrainBin,status='NEW',
     1              access='STREAM',form='UNFORMATTED',
     2              position='APPEND')
            else
               ! Open Binary File & Write Data
               open(unit=20,file=StressBin,status='OLD',
     1              access='STREAM',form='UNFORMATTED',
     2              position='APPEND')             
               open(unit=30,file=StrainBin,status='OLD',
     1              access='STREAM',form='UNFORMATTED',
     2              position='APPEND')
            endif
            write(20) Stress
            write(30) Strain
            close(20)
            close(30)
         endif
         return
      end subroutine WriteEle

      subroutine FindEnd(datapath,IsEnd,TIME,DTIME,SVARS)
         character(len=200):: datapath
         character(len=13):: endStr
         real*8:: timeOld, TIME(2),timeNew,DTIME,SVARS(32)
         logical:: Isend
         Isend = .false.
         timeNew = TIME(2)
         if (timeNew.eq.0.0) then
            Isend = .false.
            close(10)
         else if ((timeNew.eq.DTIME)) then
            if (SVARS(1) .eq. 0.0) then
               Isend = .true.
            else
               Isend = .false.
            endif
            close(10)
         else
            open(10,file=datapath,status='OLD',position='APPEND')
            backspace(10)
            read(10,"(A13)") EndStr
            if (EndStr .eq. " ============")  then
               backspace(10)
               backspace(10)
               read(10,"(A13,E16.6)") EndStr, timeOld
               if(timeOld.eq.timeNew) then
                  Isend = .true.
               endif
            endif
            close(10)
         endif
         return
      end subroutine FindEnd