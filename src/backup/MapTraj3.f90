!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Program :                                                           !
!       Mapping the trajectories to numerical PES by RMSD               !
!                                                                       !
!   Input command-line arguments:                                       !
!       $1 = single trajectory with separate forward and reverse parts  !
!               format:                                                 !
!                       # of atom                                       !
!                       comment (TSS has 'runpoint 1' in this line)     !
!                       #atom       x   y   z                           !
!       $2 = *_E.dat                                                    !
!       $3 = *_Struc.xyz                                                !
!       $4 = constrain of reverse direction; 'R', "P1/P2' or 'none'.    !
!       $5 = constrain of forward direction; 'R', "P1/P2' or 'none'.    !
!                                                                       !
!   Output :                                                            !
!           coord.$1                                                    !
!                                                                       !
!   History:                                                            !
! 2018/10/12, Grace                                                     !
! 2019/06/04, Grace, modify the structure of code. Map trajectory from  !
!   the transition-state structurecd .                                  !
! 2019/06/10, Grace, Tune the coordinate of trajectory; add function    !
!   tuneTraj().                                                         !
! 2019/06/11, Grace, Add different constrains in forward and reverse    !
!   direction in order to make sure the well-behaved trajectory.        !
!   Modified the subroutine; findNextPts().                             !
! 2019/07/23, Grace, add the third input argument which record PES in   !
!   mass-weighted coordinate. Shift coordinate of trajectory in one grid! 
!   area and then also output its coordinate in mass-weighted unit;     !
!   Change the coordinate from the amount of grid points to             !
!   mass-weighted coordinate; add subroutine gp2mw().                   !
! 2019/07/26, Grace, rewrite this code since the memory of passing      !
!   variables may have overload problem in get_vec() subroutine. Also   !
!   reduce the memory using in the new version.                         !
!   reference: https://en.wikibooks.org/wiki/Fortran/memory_management  !
! 2019/07/31, Grace, modify the searching method, and add the 4th and   !
!   5th command-line arguments in order to set up constraint the searing!
!   direction. Remove right neighbors for half trajectory which goes to !
!   reactant, and remove left neighbors for half trajectory which goes  !
!   to the product region; for both product 1 and product 2.            !
! 2019/08/22, Grace, modify the last step of searching coordinate; in   !
!   one brick area, tune the cooridnate of trajectory by using its      !
!   gradient in both x and y directions. Replace tuneTraj() as          !
!   tuneTraj2().                                                        !
! 2019/09/03, Grace, modify the boundary condition in the subroutine,   !
!   filter2TwoNeighbors(), which should not pass index equal to zero.   !
! 2019/09/04, Grace, modify tuneTraj2() while tunning the coordinate,   !
!   it is not linear such that if the new coordinate is out of range,   !
!   redefine the boundary again.                                        !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module gloable
    implicit none
    integer(4)  :: NAtoms ! number of atoms 
    integer(4)  :: numPts ! amount of points in 2D-PES 
    real(8),allocatable,dimension(:,:)   :: PESList
    character(len=100)  :: output

    contains 

    subroutine print_purpose()
        implicit none
        character(len=100)   :: filename
    
        call GETARG(1,filename)
        write(*,'(A)') 'output file: coord.'//TRIM(ADJUSTL(filename))
        return
    end subroutine print_purpose

    subroutine get_NAtoms(NAtoms)
        implicit none
        integer(4),intent(out)  :: NAtoms
        character(len=100)  :: struc
        logical :: filestat
    
        call GETARG(1,struc)
        INQUIRE(file=struc,exist=filestat)
        if (filestat) then
            open(10,file=struc,status='old')
        else
            write(*,'(A)') TRIM(struc)//" doesn't exist"
            stop
        end if
        read(10,*) NAtoms
        close(10)
        return
    end subroutine get_NAtoms

    subroutine get_numPts(NAtoms,numPts)
        implicit none
        integer(4),intent(in)   :: NAtoms
        integer(4),intent(out)   :: numPts
        ! local variable
        character(len=100)  :: filename

        call GETARG(2,filename)
        call system('wc -l '//TRIM(filename)//" | awk '{print $1}' > totalLine.dat")
        open(10,file='totalLine.dat',status='old')
        read(10,*)  numPts
        close(10,status='delete')
        return
    end subroutine get_numPts

    subroutine get_PESList(numPts,PESList)
        implicit none
        integer(4),intent(in)   :: numPts
        real(8),dimension(numPts,3),intent(out) :: PESList
        ! local variable
        character(len=100)   :: filename
        logical :: filestat
        integer(4)  :: i
    
        call GETARG(2,filename)
        INQUIRE(file=filename,exist=filestat)
        open(10,file=filename,status='old',action='read')
        do i=1,numPts
            read(10,*) PESList(i,1:3)
        end do
        close(10)
        return
    end subroutine get_PESList

    subroutine get_output(output)
        implicit none
        character(len=100),intent(out)  :: output
        ! local variable
        character(len=100)  :: filename 
        call GETARG(1,filename)
        output = 'coord.'//TRIM(ADJUSTL(filename))
        return 
    end subroutine get_output
end module gloable  

Module extractData
    implicit none
    contains 
    subroutine readStruc(filename,NAtoms,coord)
        implicit none
        character(len=100),intent(in)   :: filename
        integer(4),intent(in)   :: NAtoms
        real(8),dimension(NAtoms,3),intent(out) :: coord
        ! locat variable 
        integer(4)  :: i
        character(len=100)  :: buffer
    
        open(10,file=TRIM(ADJUSTL(filename)),status='old',action='read')
        read(10,*) buffer
        read(10,*) buffer
        do i=1,NAtoms
            read(10,*) buffer,coord(i,1:3)
        end do
        close(10,status='delete')
        return
    end subroutine readStruc

    subroutine pickTssTraj(NAtoms,tssTraj)
        implicit none 
        integer(4),intent(in)   :: NAtoms 
        real(8),dimension(NAtoms,3) :: tssTraj
        ! local variables 
        character(len=100)  :: filename, jobL 
        integer(4)  :: jobNum 

        call GETARG(1,filename)
        jobNum =  NAtoms + 2
        write(jobL,*)jobNum 
        call system('head -n '//TRIM(ADJUSTL(jobL))//' ' &
            //TRIM(ADJUSTL(filename))//' > tssTraj.xyz')
        filename='tssTraj.xyz'
        call readStruc(filename,NAtoms,tssTraj)
        return 
    end subroutine pickTssTraj

    subroutine pickPESstruc(index,NAtoms,coord)
        implicit none
        integer(4),intent(in)   :: index, NAtoms
        real(8),dimension(NAtoms,3),intent(out) :: coord
        ! local variables
        integer(4)  :: i, jobNum, iniNum, finNum
        character(len=100)  :: filename , jobL, iniL, finL 

        jobNum = NAtoms + 2 
        iniNum = 1 + jobNum * ( index - 1 )
        finNum = index * jobNum 
        write(iniL,*) iniNum
        write(finL,*) finNum 

        call GETARG(3, filename)
        call system('sed -n "'//TRIM(ADJUSTL(iniL))//',' &
            //TRIM(ADJUSTL(finL))//' p" '//TRIM(ADJUSTL(filename)) &
            //' > tmp.xyz')
        filename='tmp.xyz'
        call readStruc(filename,NAtoms,coord)
        return
    end subroutine pickPESstruc

    subroutine splitTraj()
        implicit none 
        integer(4)  :: num, totline
        character(len=100)  :: filename, lineNum 

        call GETARG(1,filename)

        call system(" grep -n 'runpoint 1' "//TRIM(ADJUSTL(filename)) & 
            //" | tail -n 1 | cut -d ':' -f 1 > lineNum.dat")
        open(10,file='lineNum.dat',status='old')
        read(10,*) num
        close(10,status='delete')

        call system('wc -l '//TRIM(ADJUSTL(filename)) & 
            //" | awk '{print $1}' > totLine.dat" )
        open(10,file='totLine.dat',status='old')
        read(10,*) totline 
        close(10,status='delete')

        write(lineNum,*) num - 2 
        call system(" head -n "//TRIM(ADJUSTL(lineNum))//' ' & 
            //TRIM(ADJUSTL(filename))//' > half1Traj.xyz')

        write(lineNum,*) totline - ( num - 2 )
        call system(" tail -n "//TRIM(ADJUSTL(lineNum))//' ' & 
            //TRIM(ADJUSTL(filename))//' > half2Traj.xyz')
        return
    end subroutine splitTraj

    subroutine reverseFile(outputFile)
        implicit none
        character(len=100),intent(in)   :: outputFile 
        ! local variables 
        integer(4)  :: totline, i
        character(len=100),allocatable,dimension(:)   :: tmpFile

        call system('wc -l '//TRIM(ADJUSTL(outputFile))//' > totline.dat')
        open(10,file='totline.dat',status='old')
        read(10,*) totline 
        close(10,status='delete')
        allocate(tmpFile(totline))

        open(10,file=outputFile,status='old')
        do i = 1, totline 
            read(10,'(A)') tmpFile(i)
        end do
        close(10)

        open(10,file=outputFile,status='replace')
        do i = totline, 1, -1 
            write(10,'(A)') tmpFile(i)
        end do
        close(10)
        return 
    end subroutine reverseFile 

    subroutine get_RMSD(NAtoms,Coord1,Coord2,rmsd)
        implicit none
        integer(4), intent(in)   :: NAtoms
        real(8), dimension(NAtoms,3), intent(in)  :: Coord1,Coord2
        real(8), intent(out) :: rmsd
        ! local variable
        real(8) :: sumDiff
        integer(4)  :: i,j
        sumDiff = 0.0D0
        do i = 1, NAtoms
            do j = 1,3
                sumDiff = sumDiff + (coord1(i,j)-coord2(i,j))**2
            end do
        end do
        rmsd = SQRT(sumDiff/(NAtoms*1.0D0)) !NAtoms must larger then 0
        return
    end subroutine get_RMSD

    subroutine vaildNeighbors(index,numPts,cstr,numNeighbor,neighbors)
        implicit none
        integer(4),intent(in)    :: index, numPts
        integer(4),intent(in)   :: numNeighbor
        character(len=100),intent(in)   :: cstr
        integer(4),dimension(numNeighbor),intent(out) :: neighbors
        ! local variable
        integer(4)  :: i, num1D
    
        num1D = INT(SQRT(numPts*1.0D0))
    
        neighbors = 0
        neighbors(1) = index + num1D        ! right
        neighbors(2) = index - num1D        ! left
        neighbors(3) = index + 1            ! upper
        neighbors(4) = index - 1            ! lower
        if ( numNeighbor .eq. 8 ) then 
            neighbors(5) = index + num1D + 1    ! upper right corner
            neighbors(6) = index + num1D - 1    ! lower right corner
            neighbors(7) = index - num1D + 1    ! upper left corner
            neighbors(8) = index - num1D - 1    ! lower left corner
        end if 
    
        ! index must stay in this surface; i.e. within the range 1 ~ numPts
        if ( ( index .lt. 1 ) .or. ( index .gt. numPts ) ) then 
            neighbors = 0
        else ! edges and corners conditions 
            ! left edge 
            if ( ( index .ge. 1 ) .and. ( index .le. num1D  ) ) then 
                select case (numNeighbor)
                    case (4)
                        neighbors(2) = 0
                    case (8)
                        neighbors(2) = 0
                        neighbors(7) = 0
                        neighbors(8) = 0
                end select 
           
                ! index has only 2 neighbors - upper and lower left corners
                    if ( index .eq. 1 ) then 
                        ! lower left corner  
                        select case (numNeighbor)
                            case (4)
                                neighbors(4) = 0
                            case (8)
                                neighbors(4) = 0
                                neighbors(6) = 0
                        end select 
    
                    else if ( index .eq. num1D ) then 
                        ! upper left corner 
                        select case (numNeighbor)
                        case (4)
                            neighbors(3) = 0
                        case (8)
                            neighbors(3) = 0
                            neighbors(5) = 0
                        end select 
                    end if
    
            ! right edge 
            else if ( ( index .ge. 1 + (num1D - 1) * num1D) .and. ( index .le. numPts) ) then 
                select case (numNeighbor)
                case (4)
                    neighbors(1) = 0
                case (8)
                    neighbors(1) = 0
                    neighbors(5) = 0
                    neighbors(6) = 0
                end select 
    
                ! index has only 2 neighbors - upper and lower right corners
                    if ( index .eq. 1 + (num1D - 1) * num1D ) then 
                        ! lower right corner
                        select case (numNeighbor)
                        case (4)
                            neighbors(4) = 0
                        case (8)
                            neighbors(4) = 0
                            neighbors(8) = 0
                        end select 
    
                    else if ( index .eq. numPts ) then 
                        ! upper right corner 
                        select case (numNeighbor)
                        case (4)
                            neighbors(3) = 0
                        case (8)
                            neighbors(3) = 0
                            neighbors(7) = 0
                        end select
                    end if 
            end if 
    
            ! upper edge with double count the upper left/right corner 
            if ( ( MODULO(index,num1D) .eq. 0 ) ) then !& 
                !.and. ( index/num1D .ne. 0 ) .and. ( index/num1D .ne. ( num1D - 1 ) ) ) then 
                select case (numNeighbor)
                case (4)
                    neighbors(3) = 0
                case (8)
                    neighbors(3) = 0
                    neighbors(5) = 0
                    neighbors(7)  = 0 
                end select  
                
            ! lower edge with double count the upper left/right corner 
            else if ( ( MODULO(index,num1D) .eq. 1 ) ) then !&
                !.and. ( index/num1D .ne. 1 ) .and. ( index/num1D .ne. num1D ) ) then 
                select case (numNeighbor)
                case (4)
                    neighbors(4) = 0
                case (8)
                    neighbors(4) = 0
                    neighbors(6) = 0
                    neighbors(8)  = 0 
                end select 
            end if 
    
        end if 
    
        ! set up constraint for reverse or forward direction 
        ! if ( numNeighbor .eq. 8 ) then 
        !     select case (cstr)
        !         case ('R')
        !             neighbors(1) = 0
        !             ! neighbors(3) = 0
        !             ! neighbors(4) = 0
        !             neighbors(5) = 0
        !             neighbors(6) = 0
        !         case ('P1')
        !             neighbors(2) = 0
        !             neighbors(7) = 0
        !             neighbors(8) = 0
        !         case ('P2')
        !             neighbors(2) = 0
        !             neighbors(7) = 0
        !             neighbors(8) = 0
        !     end select
        ! end if
        return 
    end subroutine vaildNeighbors

    subroutine filter2TwoNeighbors(NAtoms,numPts,PESList,coordTraj,center,indexX,indexY)
        implicit none 
        integer(4),intent(in)   :: NAtoms, numPts, center
        real(8),dimension(numPts,3),intent(in) :: PESList
        real(8),dimension(NAtoms,3),intent(in) :: coordTraj   
        integer(4),intent(out)  :: indexX, indexY
        ! local variable 
        integer(4)  :: i
        character(len=100)  :: buff 
        integer(4),parameter    :: numNeighbor = 4
        integer(4),dimension(numNeighbor) :: neighbors
        real(8),dimension(numNeighbor,NAtoms,3)   :: coordNeighbors
        real(8),dimension(NAtoms,3) :: matrixTmp1, matrixTmp2
    
        ! 1. Extract the 4 indeice around the input index (i.e. variable: center), 
        !   and also remove the invalide situations, i.e. edge and corner of 2D-PES. 
            call vaildNeighbors(center,numPts,buff,numNeighbor,neighbors)
    
        ! 2. Extract the corresponding structures of the 4 neighbors. After that, 
        !   select only 2 vaild neighbors by comparing the structure difference 
        !    between the coordinate of trajectory and (upper,below) or (left,right).
            coordNeighbors = 0.0D0 
            do i = 1, 4
                if ( neighbors(i) .eq. 0 ) cycle 
                call pickPESstruc(neighbors(i),NAtoms,matrixTmp1)
                !2019/07/27, Grace, cannot pass 3D array, such that create a tmp 2D array 
                coordNeighbors(i,:,:) = matrixTmp1 
            end do
    
            ! Compare right and left entry, and then select one of them 
            indexX = 0
            ! 2019/09/03, Grace, boundary should not be zero
            ! write(*,*) neighbors(1), neighbors(2)
            ! if ( ( neighbors(1) .ne. 0 ) .and. (neighbors(2) .ne. 0 ) ) then  
                matrixTmp1 = coordNeighbors(1,:,:)
                matrixTmp2 = coordNeighbors(2,:,:)
                call TwoNeighbors(NAtoms,coordTraj,neighbors(1),neighbors(2),matrixTmp1,matrixTmp2)
                if ( neighbors(1) .eq. 0 ) then 
                    indexX = neighbors(2)
                else 
                    indexX = neighbors(1)
                end if
                
                ! if ( indexX .eq. 0 ) then 
                !     indexX = center
                ! end if
            ! end if 
    
            ! Compare upper and lower entry, and then select one of them 
            indexY = 0
            ! if ( ( neighbors(3) .ne. 0 ) .and. (neighbors(4) .ne. 0 ) ) then  
                matrixTmp1 = coordNeighbors(3,:,:)
                matrixTmp2 = coordNeighbors(4,:,:)
                call TwoNeighbors(NAtoms,coordTraj,neighbors(3),neighbors(4),matrixTmp1,matrixTmp2)
                if ( neighbors(3) .eq. 0 ) then 
                    indexY = neighbors(4)
                else 
                    indexY = neighbors(3)
                end if
                ! 2019/09/03, Grace, boundary should not be zero
            !     if ( indexY .eq. 0 ) then 
            !         indexY = center
            !     end if
            ! end if  
            ! write(*,*) indexX, indexY 
        return 
    end subroutine filter2TwoNeighbors

    subroutine TwoNeighbors(NAtoms,coordTraj,index1,index2,coord1,coord2)
        implicit none
        integer(4),intent(in)   :: NAtoms
        real(8),dimension(NAtoms,3),intent(in) :: coordTraj, coord1, coord2
        integer(4),intent(inout)    :: index1, index2
        ! local variables
        integer(4)  :: i, j
        real(8) :: length1, length2
        real(8),dimension(NAtoms,3) :: diff1, diff2
    
        diff1 = 0.0D0 
        diff2 = 0.0D0
        do j = 1, 3
            do i = 1, NAtoms 
                diff1(i,j) = coordTraj(i,j) - coord1(i,j)
                diff2(i,j) = coordTraj(i,j) - coord2(i,j)
                length1 = length1 + diff1(i,j) * diff1(i,j)
                length2 = length2 + diff2(i,j) * diff2(i,j)
            end do
        end do
        length1 = SQRT(length1)
        length2 = SQRT(length2)
    
        if ( length1 .lt. length2 ) then 
            index2 = 0 
        else 
            index1 = 0 
        end if
        return 
    end subroutine TwoNeighbors

    subroutine filter1Neighbor(NAtoms,numPts,coordTraj,cstr,inIndex,outIndex)
        implicit none 
        integer(4),intent(in)   :: NAtoms, numPts, inIndex
        real(8),dimension(NAtoms,3),intent(in)  :: coordTraj
        character(len=100),intent(in)   :: cstr
        integer(4),intent(out)  :: outIndex 
        ! local variables 
        integer(4),parameter :: numNeighbor = 8 ! FIXME: test 4 or 8 neighbors 
        integer(4),dimension(numNeighbor) :: neighbors  
        integer(4)  :: i
        real(8) :: rmsd1, rmsd2 
        real(8),dimension(NAtoms,3) :: coordPES
    
        call vaildNeighbors(inIndex,numPts,cstr,numNeighbor,neighbors)
        rmsd1 = 10000.D0 ! to be replaced
        do i = 1, numNeighbor
            if ( neighbors(i) .eq. 0 ) cycle
            call pickPESstruc(neighbors(i),NAtoms,coordPES)
            call get_RMSD(NAtoms,coordTraj,coordPES,rmsd2)
            if ( rmsd1 .gt. rmsd2 ) then 
                rmsd1 = rmsd2
                outIndex = neighbors(i)
            end if 
        end do 
        return 
    end subroutine filter1Neighbor

    
    subroutine get_dRMSD(NAtoms,coordTraj,coordCenter,vec,interval,length,dRMSD)
        ! Implement centered finite difference method, accuracy: O(h^2) 
        ! reference: Computational science and engineering I, MIT, G. Strang, Lecture 2
        ! https://ocw.mit.edu/courses/mathematics/18-085-computational-science-and-engineering-i-fall-2008/video-lectures/lecture-2-differential-eqns-and-difference-eqns/
        implicit none
        integer(4),intent(in)   :: NAtoms
        real(8),dimension(NAtoms,3),intent(in)  :: coordTraj, coordCenter, vec 
        real(8),intent(in)  :: interval, length
        real(8),intent(out) :: dRMSD
        ! local variables 
        real(8), parameter :: pts = 100.0D0 
        real(8) :: dx  ! interval is the length of one grid, and dx is length / pts 
        integer(4)  :: i, j
        real(8) :: RMSD_F, RMSD_B
        real(8),dimension(NAtoms,3) :: coordShift, coordForward, coordBackward

        dx = interval / pts
        coordShift = 0.0D0 
        coordForward = 0.0D0 
        coordBackward = 0.0D0 
        do j = 1, 3 
            do i = 1, NAtoms 
                coordShift(i,j) = coordCenter(i,j) + vec(i,j) * length 
                coordForward(i,j) =  coordShift(i,j) + vec(i,j) / pts
                coordBackward(i,j) =  coordShift(i,j) - vec(i,j) /pts
            end do 
        end do
        call get_RMSD(NAtoms,coordTraj,coordForward,RMSD_F)
        call get_RMSD(NAtoms,coordTraj,coordBackward,RMSD_B)
        ! dRMSD = RMSD_F
        dRMSD = ( RMSD_F - RMSD_B ) / (2 * dx)

        return 
    end subroutine get_dRMSD

end module extractData

Program main
    use gloable
    implicit none
    integer(4)  :: tss_index    ! index of tss which can search the 
                                ! corresponding coordinate and energy 
    integer(4)  :: fileID

    ! Step 1. I/O import trajectory and numerical PES 
        call print_purpose()    
        call get_NAtoms(NAtoms)
        call get_numPts(NAtoms,numPts)
        allocate(PESList(numPts,3))
        call get_PESList(numPts,PESList)
        
    ! Step 2. Find the initial coordinate of the first point; TSS.
        call findFirstPt(NAtoms,numPts,PESList,tss_index)

    ! Step 3. Compare the RMSD to its neighborhood, and find the next point
        call get_output(output) ! output: coord.$1 
        fileID = 100 
        open(fileID,file=output,status='replace')
        
        call findOtherPts(fileID,NAtoms,numPts,PESList,tss_index)

        deallocate(PESList)
        close(fileID)
    stop
End program main

subroutine findFirstPt(NAtoms,numPts,PESList,tss_index)
    use extractData
    implicit none
    integer(4),intent(in)   :: NAtoms, numPts
    real(8),dimension(NAtoms,3),intent(in)   :: PESList
    integer(4),intent(out)  :: tss_index 
    ! local variables 
    integer(4)  :: i, jobNum, iniNum, finNum
    character(len=100)  :: struc, filename , jobL, iniL, finL 
    real(8),dimension(NAtoms,3) :: tssTraj, coordPES
    real(8) :: rmsd1,rmsd2
    real(8), dimension(2)   :: tune_coord

    ! 1. Extract the tss of the trajectory
        call pickTssTraj(NAtoms,tssTraj)

    ! 2. Search the corresponding coordinate in this PES 
        call GETARG(3,struc)
        jobNum =  NAtoms + 2
        write(jobL,*)jobNum 
        call system('head -n '//TRIM(ADJUSTL(jobL))//' ' &
        //TRIM(ADJUSTL(struc))//' > 1ptsPES.xyz')
        filename='1ptsPES.xyz'
        call readStruc(filename,NAtoms,coordPES)
        call get_RMSD(NAtoms,tssTraj,coordPES,rmsd1)

        call GETARG(2,filename)
        ! FIXME: only for test system : H2O 
        ! call system("grep -n '0.97587863 108.25046046' "//TRIM(ADJUSTL(filename)) &
        call system("grep -n '0.0000 0.0000' "//TRIM(ADJUSTL(filename)) &
            //" | cut -d ':' -f 1 > line.dat")

        open(10,file='line.dat',status='old')
        read(10,*) tss_index
        close(10,status='delete')

    return 
end subroutine findFirstPt

subroutine findOtherPts(fileID,NAtoms,numPts,PESList,tss_index)
    use extractData
    implicit none
    integer(4),intent(in)   :: fileID, NAtoms, numPts,tss_index
    real(8),dimension(numPts,3),intent(in)  :: PESList 
    ! local variables
    character(len=100)  :: inputFile, outputFile, finalOutput, buffer 
    character(len=100)  :: cstr_r, cstr_f ! constraint for reverse/forward directions
    integer(4)  :: numLine, i
    

    ! 1. Split trajectory into two parts
        call splitTraj() !output: half1Traj.xyz and half2Traj.xyz 

    ! 2. Assign the final output filename 
        call GETARG(1,inputFile)
        finalOutput='coord.'//TRIM(ADJUSTL(inputFile)) 

    ! 3. Search the coordinate for the first half trajectory, 
    !   and then reverse its order. TSS coordinate is also removed, 
    !   since the later half trajectory already have it.
        inputFile = 'half1Traj.xyz'
        outputFile = 'coord.'//TRIM(ADJUSTL(inputFile))
        call GETARG(4,cstr_r)
        write(*,*) tss_index
        call coordHalfTraj(NAtoms,numPts,PESList,tss_index,cstr_r,inputFile,outputFile)
        write(*,*) tss_index
stop
        call reverseFile(outputFile)

        call system('wc -l '//TRIM(ADJUSTL(outputFile))//' > totline.dat')
        open(10,file='totline.dat',status='old')
        read(10,*) numLine 
        close(10,status='delete')

        open(10,file=outputFile,status='old')
        open(999,file=finalOutput,status='replace')
        do i = 1, numLine - 1 ! remove TSS
            read(10,'(A)') buffer
            write(999,'(A)') TRIM(ADJUSTL(buffer))
        end do
        ! close(10)
        close(10,status='delete')
        write(*,'()')
        write(*,'(A)') 'Half'
        write(*,'()')
        
    ! 4. Search the coordinate for the last half trajectory, print the TSS
    !   coordinate as well.
        inputFile = 'half2Traj.xyz'
        outputFile = 'coord.'//TRIM(ADJUSTL(inputFile))
        call GETARG(5,cstr_f)
        call coordHalfTraj(NAtoms,numPts,PESList,tss_index,cstr_f,inputFile,outputFile)

        call system('wc -l '//TRIM(ADJUSTL(outputFile))//' > totline.dat')
        open(10,file='totline.dat',status='old')
        read(10,*) numLine 
        close(10,status='delete')

        open(10,file=outputFile,status='old')
        do i = 1, numLine 
            read(10,'(A)') buffer
            write(999,'(A)') TRIM(ADJUSTL(buffer))
        end do
        ! close(10)
        close(10,status='delete')

        close(999)
        call system('rm -f half1Traj.xyz half2Traj.xyz' )
        ! write(*,*) cstr_r,cstr_f

    return 
end subroutine findOtherPts

subroutine coordHalfTraj(NAtoms,numPts,PESList,tss_index,cstr,inputFile,outputFile)
    use extractData
    implicit none 
    integer(4),intent(in)   :: NAtoms, numPts, tss_index 
    real(8),dimension(numPts,3),intent(in)  :: PESList 
    character(len=100),intent(in)   :: cstr, inputFile, outputFile 
    ! local variables 
    real(8),dimension(2)    :: tune_coord
    real(8),dimension(NAtoms,3) :: tssTraj, oneCoord 
    character(len=100)  :: charline, charline2, filename 
    integer(4)  :: numline, oneJob, numjobs, inIndex, outIndex
    integer(4)  :: i, iniL, finL

    ! 1.  Count the total amount of points (i.e. numjobs).
        call system("wc -l "//TRIM(ADJUSTL(inputFile)) & 
            //" | awk '{print $1}' > charline.dat")
        open(10,file='charline.dat',status='old')
        read(10,*) numline
        close(10,status='delete')
        oneJob = NAtoms + 2 
        numjobs = numline / oneJob

    ! 2. The coordinate of the first point (TSS) is already founded, 
    !   tune this coordinate, and then save it in outputFile
        inIndex = tss_index 
        call pickTssTraj(NAtoms,tssTraj)
        call tuneTraj2(NAtoms,numPts,PESList,tssTraj,inIndex,tune_coord)
        open(20,file=outputFile,status='replace')
        write(20,"(3(F15.8,1X))") tune_coord(1:2),PESList(inIndex,3)
        ! write(*,'()')
        ! write(*,"(3(F15.8,1X))") tune_coord(1:2),PESList(inIndex,3)
        ! stop
    
    ! 3. Search the coordinate of the following  points.
        inIndex = tss_index 
        ! write(*,*) TRIM(ADJUSTL(inputFile))
        do i = 2, numjobs
            iniL = 1 + oneJob * ( i - 1 ) 
            finL = oneJob * i
            write(charLine,*) iniL
            write(charLine2,*) finL 
            
            call system("sed -n '"//TRIM(ADJUSTL(charLine))//',' &
                //TRIM(ADJUSTL(charline2))//" p' " & 
                //TRIM(ADJUSTL(inputFile))// " > oneCoord.xyz")
            filename = 'oneCoord.xyz'
            call readStruc(filename,NAtoms,oneCoord)

            ! Compare the neighest 8 neighbors, select the one with the smallest RMSD
            call filter1Neighbor(NAtoms,numPts,oneCoord,cstr,inIndex,outIndex)
            
            ! Tune the coordinate in both x and y directions
            call tuneTraj2(NAtoms,numPts,PESList,tssTraj,outIndex,tune_coord)

            inIndex = outIndex 
            write(20,"(3(F15.8,1X))") tune_coord(1:2),PESList(outIndex,3)
            ! write(*,"(3(F15.8,1X))") tune_coord(1:2),PESList(tss_index,3)
        end do 
        close(10)
    return 
end subroutine coordHalfTraj

subroutine tuneTraj(NAtoms,numPts,PESList,coordTraj,center,coord)
    use extractData ! call pickPESstruc
    implicit none
    integer(4),intent(in)   :: NAtoms, numPts, center
    real(8),dimension(numPts,3),intent(in) :: PESList
    real(8),dimension(NAtoms,3),intent(in) :: coordTraj   
    real(8),dimension(2),intent(out)    :: coord
    ! local variables
    integer(4)  :: i, j, indexX, indexY
    real(8),dimension(2)    :: x0y0, horizontal, vertical, interval, dotProd, norm
    real(8),dimension(NAtoms,3) :: coordCenter, coordHor, coordVer
    real(8),dimension(NAtoms,3) :: vecX, vecY, vecTraj_Center
    
    ! x = x0 + dot(v,+dx) * dx 
    ! y = y0 + dot(v,+dy) * dy

    ! 1. Extract x0 and y0
    !   x0 = x0y0(1) and y0 = x0y0(2). 
    !   And then extract the the corresponding structure; i.e. coordCenter(NAtoms, 3)
        do i = 1, 2
            x0y0(i) = PESList(center,i)
        end do  
        call pickPESstruc(center,NAtoms,coordCenter)

    ! 2. Select its valid neighbors; get indexX and indexY. 
    !   Select one horizontal and one vertical neighbors with smaller RMSD from 
    !   four nearest-neighbor; right, left, upper and lower neighbors. 
        call filter2TwoNeighbors(NAtoms,numPts,PESList,coordTraj,center,indexX,indexY)
        do i = 1, 2 
            horizontal(i) = PESList(indexX,i)
            vertical(i) = PESList(indexY,i)
        end do
        call pickPESstruc(indexX,NAtoms,coordHor)
        call pickPESstruc(indexY,NAtoms,coordVer)

        interval = 0.0D0 
        interval(1) = horizontal(1) - x0y0(1)
        interval(2) = vertical(2) - x0y0(2)

    ! 3.  Calculate dot product; get dotProd.
        !   3.1. Calculate the vaild normalized x vector and y vector, 
        !   which are the structure differences between the center and 
        !   its horizontal and vertical neighbors. And then normalize 
        !   these vectors.
            vecX = 0.0D0 
            vecY = 0.0D0 
            norm = 0.0D0 
            do j = 1, 3
                do i = 1, NAtoms 
                    vecX(i,j) = coordHor(i,j) - coordCenter(i,j)
                    vecY(i,j) = coordVer(i,j) - coordCenter(i,j)
                    norm(1) = norm(1) + vecX(i,j) ** 2
                    norm(2) = norm(2) + vecY(i,j) ** 2
                end do 
            end do 
            norm(1) = SQRT(norm(1))
            norm(2) = SQRT(norm(2))
            do j = 1, 3
                do i = 1, NAtoms 
                    vecX(i,j) = vecX(i,j)/norm(1)
                    vecY(i,j) = vecY(i,j)/norm(2)
                end do
            end do

        ! 3.2. Calculate the structure difference between the trajectory and the 
        !   center structure of this surface.
            vecTraj_Center = 0.0D0
            do j = 1, 3
                do i = 1, NAtoms 
                    vecTraj_Center(i,j) = coordTraj(i,j) - coordCenter(i,j)
                end do
            end do

        ! 3.3. Calculate dot products
            dotProd = 0.0D0 
            do j = 1, 3 
                do i = 1, NAtoms 
                    dotProd(1) = dotProd(1) + vecTraj_Center(i,j) * vecX(i,j)
                    dotProd(2) = dotProd(2) + vecTraj_Center(i,j) * vecY(i,j)
                end do
            end do

    ! 4. Calculate the new coordinate; get coord.
        coord = 0.0D0
        do i = 1, 2 
            coord(i) = x0y0(i) + dotProd(i) * interval(i)
        end do 
        
    return 
end subroutine tuneTraj

subroutine tuneTraj2(NAtoms,numPts,PESList,coordTraj,center,coord)
    use extractData
    implicit none
    integer(4),intent(in)   :: NAtoms, numPts, center
    real(8),dimension(numPts,3),intent(in) :: PESList
    real(8),dimension(NAtoms,3),intent(in) :: coordTraj   
    real(8),dimension(2),intent(out)    :: coord
    ! local variables 
    integer(4)  :: i

    ! x = x0 + shift_x
    ! y = y0 + shift_y

        ! 3.2. Calculate new coordinate if its derivative of RMSD is near to zero.
            coord = 0.0D0 
            do i = 1, 2
                call get_zero_dRMSD(NAtoms,numPts,PESList,coordTraj,center,i,coord(i))
            end do 
            ! write(*,*) coord 
    return 
end subroutine tuneTraj2

subroutine get_zero_dRMSD(NAtoms,numPts,PESList,coordTraj,center,XorY,coord)
    use extractData
    implicit none 
    integer(4),intent(in)   :: NAtoms, numPts, XorY
    real(8),dimension(numPts,3),intent(in) :: PESList
    real(8),dimension(NAtoms,3),intent(in)  :: coordTraj
    integer(4),intent(inout)   :: center
    real(8),intent(out)  :: coord
    ! local variables 
    integer(4)  :: iter, i, min_iter, indexX, indexY, tmpX, tmpY
    real(8), parameter  :: zero = 0.001D0
    real(8), parameter  :: MaxDer = 1.0D0 
    real(8) :: interval, coord_ini
    real(8) :: length, minLini, minLfin, df1, df2, Prod1, Prod2
    real(8),dimension(2)    :: intervalXandY
    real(8),dimension(NAtoms,3) :: vec, vecX, vecY, coordCenter

    ! 1. Calculate the derivative of RMSD at the initial point
        call get_shift(NAtoms,numPts,PESList,coordTraj,center,indexX,indexY,intervalXandY,vecX,vecY) 
        !output: indexX, indexY, interval, vecX, vecY
        if ( XorY .eq. 1 ) then 
            coord_ini = PESList(center,1)
            vec = vecX
            interval = intervalXandY(1)
        else 
            coord_ini = PESList(center,2)
            vec = vecY 
            interval = intervalXandY(2)
        end if
        call pickPESstruc(center,NAtoms,coordCenter)

        length = 0.0D0 
        call get_dRMSD(NAtoms,coordTraj,coordCenter,vec,interval,length,df1) !output: df1 = dRMSD

    ! 2. Calculate the derivative of RMSD at the half length 
        length = 0.5D0 
        call get_dRMSD(NAtoms,coordTraj,coordCenter,vec,interval,length,df2) !output: df2 = dRMSD
    
    ! 3. Condition: stop or iterative 
        coord = coord_ini + interval * length
        Prod1 = df1 * df2 
        Prod2 = Prod1
        iter = 1 ! first iter. compares the initial point and the point at 0.5D0 
        min_iter = 0

        ! 2019/09/04, Grace, discard this trajectory if the difference 
        !   between the structure of trajectory and the structure of 
        !   potential is too much. c.f. df1 = dRMSD/dx. 
        if ( ( df1 .lt. - MaxDer ) .or. ( df1 .gt. MaxDer ) ) then 
            write(*,'(A)') 'Fail at mapping trajtory'
            write(*,'(A)') 'Error: get_zero_dRMSD() in MapTraj.f90'
            stop
        end if 

        write(*,'(A)') '    iter, center, indexX, indexY,df1, df2, length, Prod1, Prod2'
        do while ( ( Prod2 .lt. - zero ) .or. ( Prod2 .gt. zero ) )       
            iter = iter + 1 
            ! works after the second loop
            Prod1 = Prod2          

            if ( Prod2 .gt. zero ) then 
                df1 = df2 
                ! move the new point toward forward direction
                length = length + 0.1D0 

                ! 2019/09/05, Grace, change vec since this PES is not linear,
                ! and 0 <= length <= 1.0D0 
                ! if ( length .ge. 1.0D0 ) then 
                !     if ( XorY .eq. 1 ) then 
                !         center = indexX
                !     else 
                !         center = indexY
                ! end if
                ! 2019/09/07, Grace, change the constraint of boundary from one 
                ! direction to two possible directions. 
                ! write(*,*) center,indexX,indexY, length  
                if ( length .ge. 1 ) then 
                    if ( XorY .eq. 1 ) then 
                        ! x-direction
                        call filter2TwoNeighbors(NAtoms,numPts,PESList,coordTraj,indexX,tmpX,tmpY)
                        center = tmpY 
                    else
                        ! y-direction
                        call filter2TwoNeighbors(NAtoms,numPts,PESList,coordTraj,indexY,tmpX,tmpY)
                        center = tmpY
                    end if
                end if
                ! write(*,*) center,indexX,indexY , length  
                call get_shift(NAtoms,numPts,PESList,coordTraj,center,indexX,indexY,intervalXandY,vecX,vecY) 
                ! write(*,*) center,indexX,indexY , length  
                ! stop
                if ( XorY .eq. 1 ) then 
                    ! x-direction
                    coord_ini = PESList(center,1)
                    vec = vecX
                    interval = intervalXandY(1)
                else 
                    ! y-direction
                    coord_ini = PESList(center,2)
                    vec = vecY 
                    interval = intervalXandY(2)
                end if

                call pickPESstruc(center,NAtoms,coordCenter)

                if ( (center .gt. indexX ) .or. ( center .gt. indexY ) ) then 
                    length = length + 0.5D0
                else
                    length = length -0.5D0 
                end if

                call get_dRMSD(NAtoms,coordTraj,coordCenter,vec,interval,length,df2)
                coord = coord_ini + interval * length
            else 
                ! search a local minimum by moving initial and final points
                min_iter = min_iter + 1 
                minLini = length - 0.1D0 + zero * min_iter 
                minLfin = length - zero * min_iter 
                call get_dRMSD(NAtoms,coordTraj,coordCenter,vec,interval,minLini,df1)
                call get_dRMSD(NAtoms,coordTraj,coordCenter,vec,interval,minLfin,df2)
                coord = coord_ini + interval * minLfin
            end if

            Prod2 = df1 * df2 
            write(*,'(I6,3(1X,I4),5(1X,F8.5))') iter, center, indexX, indexY,df1, df2, length, Prod1, Prod2
        end do
        ! write(*,*) iter, min_iter
        ! stop
    return 
end subroutine get_zero_dRMSD

subroutine get_shift(NAtoms,numPts,PESList,coordTraj,center,indexX,indexY,interval,vecX,vecY)
    use extractData
    implicit none 
    integer(4),intent(in)   :: NAtoms, numPts, center
    real(8),dimension(numPts,3),intent(in) :: PESList
    real(8),dimension(NAtoms,3),intent(in) :: coordTraj   
    real(8),dimension(2),intent(out)    :: interval
    integer(4),intent(out)  :: indexX, indexY
    real(8),dimension(NAtoms,3),intent(out) :: vecX, vecY
    ! local variable
    integer(4)  :: i,j
    real(8),dimension(2)    :: x0y0, horizontal, vertical
    real(8),dimension(NAtoms,3) :: coordCenter, coordHor, coordVer


    ! 1. Extract x0 and y0
    !   x0 = x0y0(1) and y0 = x0y0(2). 
    !   And then extract the the corresponding structure; i.e. coordCenter(NAtoms, 3)
        do i = 1, 2
            x0y0(i) = PESList(center,i)
        end do  
        call pickPESstruc(center,NAtoms,coordCenter)

    ! 2. Select its valid neighbors; get indexX and indexY. 
    !   Select one horizontal and one vertical neighbors with smaller RMSD from 
    !   four nearest neighbor; right, left, upper and lower neighbors. 
        call filter2TwoNeighbors(NAtoms,numPts,PESList,coordTraj,center,indexX,indexY)
        do i = 1, 2 
            horizontal(i) = PESList(indexX,i)
            vertical(i) = PESList(indexY,i)
        end do
        call pickPESstruc(indexX,NAtoms,coordHor)
        call pickPESstruc(indexY,NAtoms,coordVer)

        interval = 0.0D0 
        interval(1) = horizontal(1) - x0y0(1)
        interval(2) = vertical(2) - x0y0(2)

    ! 3. Calculate vector without normalization 
        vecX = 0.0D0 
        vecY = 0.0D0 
        do j = 1, 3
            do i = 1, NAtoms 
                vecX(i,j) = coordHor(i,j) - coordCenter(i,j)
                vecY(i,j) = coordVer(i,j) - coordCenter(i,j)
            end do 
        end do 
    return 
end subroutine get_shift