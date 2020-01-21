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
!   the transition-state structurec.                                    !
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
! 2019/09/10, Grace, debug of the subroutine findFirstPt(), add loop    !
!   back which is removed before because of considering the heavy cost. !
!   Modify get_zero_dRMSD(); in a local minimum region, use uneven      !
!   length to search local minimum.                                     !
! 2019/09/26, Grace, expand step2_select1fromNeighbors() form the first !
!   shell to the third shell. It converges at the H3CO system (i.e.     !
!   check the subpages of H3CO at the excel; Summary.xls). For tune     !
!   trajectories, locate a unimodal in 2D surface first, and then       !
!   implement 2D golden section search in tuneTraj().                   !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module gloable
    implicit none
    integer(4)  :: NAtoms ! number of atoms 
    integer(4)  :: numPts ! amount of points in 2D-PES 
    real(8),allocatable,dimension(:,:)   :: PESList

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

end module gloable  

Module extractData
    implicit none
    contains 

    subroutine get_output(output)
        implicit none
        character(len=100),intent(out)  :: output
        ! local variable
        character(len=100)  :: filename 
        call GETARG(1,filename)
        output = 'coord.'//TRIM(ADJUSTL(filename))
        return 
    end subroutine get_output

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
        ! The first point of trajectory is transition-state structure
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
        do j = 1, 3 
            do i = 1, Natoms 
                sumDiff = sumDiff + ( coord1(i,j) - coord2(i,j) )**2
            end do
        end do
        rmsd = SQRT(sumDiff/(NAtoms*1.0D0)) ! NAtoms must larger then 0
        return
    end subroutine get_RMSD

    subroutine vaildNeighbors(index,numPts,cstr,numNeighbor,neighbors)
        implicit none
        integer(4),intent(in)    :: index, numPts
        integer(4),intent(in)   :: numNeighbor
        character(len=100),intent(in)   :: cstr ! constraint for reverse/forward directions
        integer(4),dimension(numNeighbor),intent(out) :: neighbors
        ! local variable
        integer(4)  :: i, num1D
    
        num1D = INT(SQRT(numPts*1.0D0))
    
        ! Step 1. Assign neighbors, which has four classes. 
        ! 1. 4 neighbors; right, left, upper and lower elements
        ! 2. 8 neighbors; the first shell of the center index
        ! 3. 24 neighbors; the second shell
        ! 4. 48 neighbors; the third shell
        if ( ( numNeighbor .eq. 4 ) .or. ( numNeighbor .eq. 8 ) & 
            .or. ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
            neighbors = 0
            neighbors(1) = index + num1D        ! right
            neighbors(2) = index - num1D        ! left
            neighbors(3) = index + 1            ! upper
            neighbors(4) = index - 1            ! lower
            ! the closest neighbors; the first shell
            if ( ( numNeighbor .eq. 8 ) .or. ( numNeighbor .eq. 24 ) &
                .or. ( numNeighbor .eq. 48 ) ) then 
                neighbors(5) = index + num1D + 1    ! upper right corner
                neighbors(6) = index + num1D - 1    ! lower right corner
                neighbors(7) = index - num1D + 1    ! upper left corner
                neighbors(8) = index - num1D - 1    ! lower left corner
            end if 
            ! the second shell 
            if ( ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                do i = -2, 2 
                    neighbors(11 + i) = index + 2 * num1D + i 
                    neighbors(22 + i) = index - 2 * num1D + i 
                end do
                neighbors(14) = neighbors(5) + 1 
                neighbors(15) = neighbors(6) - 1 
                neighbors(16) = neighbors(3) + 1
                neighbors(17) = neighbors(4) - 1
                neighbors(18) = neighbors(7) + 1
                neighbors(19) = neighbors(8) - 1 
            end if
            if ( numNeighbor .eq. 48 ) then 
                do i = -3, 3
                    neighbors(28 + i ) = index + 3 * num1D + i 
                    neighbors(45 + i ) = index + 3 * num1D + i 
                end do 
                neighbors(32) = neighbors(9) + 1 
                neighbors(33) = neighbors(13) - 1 
                neighbors(34) = neighbors(14) + 1
                neighbors(35) = neighbors(15) - 1
                neighbors(36) = neighbors(16) + 1
                neighbors(37) = neighbors(17) - 1 
                neighbors(38) = neighbors(18) + 1 
                neighbors(39) = neighbors(19) - 1 
                neighbors(40) = neighbors(20) + 1
                neighbors(41) = neighbors(21) - 1
            end if
        end if 
    
        ! Step 2. Indeices must stay in this surface; i.e. within the range from 1 to numPts
            do i = 1, numNeighbor
                if ( neighbors(i) .gt. numPts) then 
                    neighbors(i) = 0 
                end if 
            end do 

        ! Step 3. Consider edge and corner conditions of the center index
        if ( ( index .lt. 1 ) .or. ( index .gt. numPts ) ) then 
            neighbors = 0
        else ! edges and corners conditions 
            
            ! Left edge
            if ( ( index .ge. 1 ) .and. ( index .le. num1D  ) ) then 
                ! left edge 
                    if ( ( numNeighbor .eq. 4 ) .or. ( numNeighbor .eq. 8 ) & 
                        .or. ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(2) = 0
                    end if
                    if ( ( numNeighbor .eq. 8 ) .or. ( numNeighbor .eq. 24 ) & 
                        .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(7) = 0
                        neighbors(8) = 0
                    end if
                    if ( ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(18) = 0
                        neighbors(19) = 0
                        do i = 20, 24
                            neighbors(i) = 0
                        end do
                    end if
                    if ( numNeighbor .eq. 48 ) then 
                        do i = 38, 48 
                            neighbors(i) = 0 
                        end do
                    end if
           
                ! upper and lower left corners
                if ( index .eq. 1 ) then 
                ! lower left corner  
                    if ( ( numNeighbor .eq. 4 ) .or. ( numNeighbor .eq. 8 ) & 
                        .or. ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(4) = 0
                    end if
                    if ( ( numNeighbor .eq. 8 ) .or. ( numNeighbor .eq. 24 ) & 
                        .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(6) = 0
                    end if
                    if ( ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(12) = 0
                        neighbors(13) = 0
                        neighbors(15) = 0
                        neighbors(17) = 0
                    end if
                    if ( numNeighbor .eq. 48 ) then 
                        neighbors(29) = 0
                        neighbors(30) = 0
                        neighbors(31) = 0
                        neighbors(33) = 0
                        neighbors(35) = 0
                        neighbors(37) = 0
                    end if
                else if ( index .eq. num1D ) then 
                ! upper left corner 
                    if ( ( numNeighbor .eq. 4 ) .or. ( numNeighbor .eq. 8 ) & 
                        .or. ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(3) = 0
                    end if
                    if ( ( numNeighbor .eq. 8 ) .or. ( numNeighbor .eq. 24 ) & 
                        .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(5) = 0
                    end if
                    if ( ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(9) = 0
                        neighbors(10) = 0
                        neighbors(14) = 0
                        neighbors(16) = 0
                    end if
                    if ( numNeighbor .eq. 48 ) then 
                        neighbors(25) = 0
                        neighbors(26) = 0
                        neighbors(27) = 0
                        neighbors(32) = 0
                        neighbors(34) = 0
                        neighbors(36) = 0
                    end if 
                end if
    
            ! Right edge
            else if ( ( index .ge. 1 + (num1D - 1) * num1D) .and. ( index .le. numPts) ) then 
                ! right edge 
                    if ( ( numNeighbor .eq. 4 ) .or. ( numNeighbor .eq. 8 ) & 
                        .or. ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(1) = 0
                    end if
                    if ( ( numNeighbor .eq. 8 ) .or. ( numNeighbor .eq. 24 ) & 
                        .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(5) = 0
                        neighbors(6) = 0
                    end if
                    if ( ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                        do i = 9, 15
                            neighbors(i) = 0
                        end do
                    end if
                    if ( numNeighbor .eq. 48 ) then 
                        do i = 25, 35
                            neighbors(i) = 0
                        end do
                    end if

                ! upper and lower right corners
                if ( index .eq. 1 + (num1D - 1) * num1D ) then 
                ! lower right corner 
                    if ( ( numNeighbor .eq. 4 ) .or. ( numNeighbor .eq. 8 ) & 
                        .or. ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(4) = 0
                    end if
                    if ( ( numNeighbor .eq. 8 ) .or. ( numNeighbor .eq. 24 ) & 
                        .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(8) = 0
                    end if
                    if ( ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(17) = 0
                        neighbors(19) = 0
                        neighbors(23) = 0
                        neighbors(24) = 0
                    end if
                    if ( numNeighbor .eq. 48 ) then 
                        neighbors(37) = 0
                        neighbors(39) = 0
                        neighbors(41) = 0
                        neighbors(46) = 0
                        neighbors(47) = 0
                        neighbors(48) = 0
                    end if
                else if ( index .eq. numPts ) then 
                ! upper right corner 
                    if ( ( numNeighbor .eq. 4 ) .or. ( numNeighbor .eq. 8 ) & 
                        .or. ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(3) = 0
                    end if
                    if ( ( numNeighbor .eq. 8 ) .or. ( numNeighbor .eq. 24 ) & 
                        .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(7) = 0
                    end if
                    if ( ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                        neighbors(16) = 0
                        neighbors(18) = 0
                        neighbors(20) = 0
                        neighbors(21) = 0
                    end if
                    if ( numNeighbor .eq. 48 ) then 
                        neighbors(36) = 0
                        neighbors(38) = 0
                        neighbors(40) = 0
                        neighbors(42) = 0
                        neighbors(43) = 0
                        neighbors(44) = 0
                    end if
                end if 
            end if 
    
            ! upper edge with double count the upper left/right corner 
            if ( ( MODULO(index,num1D) .eq. 0 ) ) then !& 
                if ( ( numNeighbor .eq. 4 ) .or. ( numNeighbor .eq. 8 ) & 
                    .or. ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                    neighbors(3) = 0
                end if
                if ( ( numNeighbor .eq. 8 ) .or. ( numNeighbor .eq. 24 ) & 
                    .or. ( numNeighbor .eq. 48 ) ) then 
                    neighbors(5) = 0
                    neighbors(7) = 0
                end if
                if ( ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                    neighbors(9) = 0
                    neighbors(10) = 0
                    neighbors(14) = 0 
                    neighbors(16) = 0
                    neighbors(18) = 0
                    neighbors(20) = 0 
                    neighbors(21) = 0 
                end if
                if ( numNeighbor .eq. 48 ) then 
                    do i = 25, 27
                        neighbors(i) = 0
                    end do
                    neighbors(32) = 0
                    neighbors(34) = 0 
                    neighbors(36) = 0
                    neighbors(38) = 0
                    neighbors(40) = 0 
                    do i = 42, 44
                        neighbors(i) = 0
                    end do
                end if
                
            ! lower edge with double count the upper left/right corner 
            else if ( ( MODULO(index,num1D) .eq. 1 ) ) then !&
                if ( ( numNeighbor .eq. 4 ) .or. ( numNeighbor .eq. 8 ) & 
                    .or. ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                    neighbors(4) = 0
                end if
                if ( ( numNeighbor .eq. 8 ) .or. ( numNeighbor .eq. 24 ) & 
                    .or. ( numNeighbor .eq. 48 ) ) then 
                    neighbors(6) = 0
                    neighbors(8) = 0
                end if
                if ( ( numNeighbor .eq. 24 ) .or. ( numNeighbor .eq. 48 ) ) then 
                    neighbors(12) = 0
                    neighbors(13) = 0
                    neighbors(15) = 0 
                    neighbors(17) = 0
                    neighbors(19) = 0
                    neighbors(23) = 0 
                    neighbors(24) = 0 
                end if
                if ( numNeighbor .eq. 48 ) then 
                    do i = 29, 31
                        neighbors(i) = 0
                    end do
                    neighbors(33) = 0
                    neighbors(35) = 0 
                    neighbors(37) = 0
                    neighbors(39) = 0
                    neighbors(41) = 0 
                    do i = 46, 48
                        neighbors(i) = 0
                    end do
                end if
            end if 
        end if 
    
        ! set up constraint for reverse or forward direction 
        ! if ( numNeighbor .eq. 8 ) then 
        !     select case (cstr)
        !         case ('defR.dat')
        !             neighbors(1) = 0
        !             neighbors(3) = 0
        !             neighbors(4) = 0
        !             neighbors(5) = 0
        !             neighbors(6) = 0
        !         case default 
        !             ! neighbors(2) = 0
        !             ! neighbors(7) = 0
        !             ! neighbors(8) = 0
        !     end select
        ! end if
        return 
    end subroutine vaildNeighbors

    subroutine filter2NearestNeighbors(NAtoms,numPts,PESList,strucTraj,center,indexHor,indexVer)
        implicit none 
        integer(4),intent(in)   :: NAtoms, numPts, center
        real(8),dimension(numPts,3),intent(in) :: PESList
        real(8),dimension(NAtoms,3),intent(in) :: strucTraj   
        integer(4),intent(out)  :: indexHor, indexVer
        ! local variable 
        integer(4)  :: i
        character(len=100)  :: buff 
        integer(4),parameter    :: numNeighbor = 4
        integer(4),dimension(numNeighbor) :: neighbors
        real(8),dimension(numNeighbor,NAtoms,3)   :: strucNeighbors
        real(8),dimension(NAtoms,3) :: strucTmp1, strucTmp2
    
        ! 1. Extract the 4 indeice around the input index (i.e. variable: center), 
        !   and also remove the invalide situations, i.e. edge and corner of 2D-PES. 
            call vaildNeighbors(center,numPts,buff,numNeighbor,neighbors)

        ! 2. Extract the corresponding structures of the 4 neighbors. After that, 
        !   select only 2 vaild neighbors by comparing the structure difference 
        !    between the coordinate of trajectory and (upper,below) or (left,right).
            strucNeighbors = 0.0D0 
            do i = 1, 4
                if ( neighbors(i) .eq. 0 ) cycle 
                call pickPESstruc(neighbors(i),NAtoms,strucTmp1)
                !2019/07/27, Grace, cannot pass 3D array, such that create a tmp 2D array 
                strucNeighbors(i,:,:) = strucTmp1 
            end do
    
            ! Compare right and left entry, and then select one of them 
                indexHor = 0
                ! 2019/09/03, Grace, boundary should not be zero
                ! write(*,*) neighbors(1), neighbors(2)
                ! if ( ( neighbors(1) .ne. 0 ) .and. (neighbors(2) .ne. 0 ) ) then  
                strucTmp1 = strucNeighbors(1,:,:)
                strucTmp2 = strucNeighbors(2,:,:)
                call TwoNeighbors(NAtoms,strucTraj,neighbors(1),neighbors(2),strucTmp1,strucTmp2)
                if ( neighbors(1) .eq. 0 ) then 
                    indexHor = neighbors(2)
                else 
                    indexHor = neighbors(1)
                end if
    
            ! Compare upper and lower entry, and then select one of them 
                indexVer = 0
                strucTmp1 = strucNeighbors(3,:,:)
                strucTmp2 = strucNeighbors(4,:,:)
                call TwoNeighbors(NAtoms,strucTraj,neighbors(3),neighbors(4),strucTmp1,strucTmp2)
                if ( neighbors(3) .eq. 0 ) then 
                    indexVer = neighbors(4)
                else 
                    indexVer = neighbors(3)
                end if
        return 
    end subroutine filter2NearestNeighbors

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

    subroutine get_length(indexCen,indexHor,indexVer,XorY,moveL,length)
        implicit none
        integer(4),intent(in)   :: indexCen, indexHor, indexVer , XorY
        real(8),intent(in)  :: moveL 
        real(8),intent(inout)   :: length 
        ! Move the new point toward forward/backward direction.
        ! The step size of one step is $moveL. 

        if ( XorY .eq. 0 ) then 
            ! x coordinate 
            if ( indexCen .lt. indexHor ) then 
                ! initial point locates at the left side of the final point.
                ! such that move to the positive direction.
                length = length + moveL
            else
                ! initial point locates at the right side of the final point.
                ! such that move to the negative direction.
                length = length - moveL
            end if
        else
            ! y coordinate 
            if ( indexCen .lt. indexVer ) then 
                ! initial point locates at the left side of the final point.
                ! such that move to the positive direction.
                length = length + moveL
            else
                ! initial point locates at the right side of the final point.
                ! such that move to the negative direction.
                length = length - moveL
            end if
        end if 
        return 
    end subroutine get_length

    subroutine get_CoordShift(NAtoms,coord0,vec,length,coord1)
        implicit none 
        integer(4),intent(in)   :: NAtoms
        real(8),dimension(NAtoms,3),intent(in)  :: coord0,vec
        real(8),intent(in)  :: length 
        real(8),dimension(NAtoms,3),intent(out) :: coord1
        ! local variables 
        integer(4)  :: i, j

        coord1 = 0.0D0 
        do j = 1, 3 
            do i = 1, NAtoms 
                coord1(i,j) = coord0(i,j) + vec(i,j) * length 
            end do
        end do

        return
    end subroutine 

    subroutine get_iniNeighbor(NAtoms,numPts,PESList,strucTraj,indexCen,indexHor,indexVer,intervalHor,intervalVer,vecHor,vecVer)
        implicit none 
        integer(4),intent(in)   :: NAtoms, numPts, indexCen
        real(8),dimension(numPts,3),intent(in) :: PESList
        real(8),dimension(NAtoms,3),intent(in) :: strucTraj   
        real(8),intent(out)    :: intervalHor, intervalVer
        integer(4), intent(out)  :: indexHor, indexVer
        real(8),dimension(NAtoms,3),intent(out) :: vecHor, vecVer
        ! local variable
        integer(4)  :: i,j
        real(8),dimension(2)    :: coordCen, coordHor, coordVer
        real(8),dimension(NAtoms,3) :: strucCen, strucHor, strucVer
    
    
        ! 1. Extract x0 and y0
        !   extract the the corresponding structure; i.e. strucCenter(NAtoms, 3)
            do i = 1, 2
                coordCen = PESList(indexCen,i)
            end do  
            call pickPESstruc(indexCen,NAtoms,strucCen)
    
        ! 2. Select its valid neighbors; get indexX and indexY. 
        !   Select one horizontal and one vertical neighbors with smaller RMSD from 
        !   four nearest neighbor; right, left, upper and lower neighbors. 
            call filter2NearestNeighbors(NAtoms,numPts,PESList,strucTraj,indexCen,indexHor,indexVer)

            do i = 1, 2 
                coordHor(i) = PESList(indexHor,i)
                coordVer(i) = PESList(indexVer,i)
            end do
            call pickPESstruc(indexHor,NAtoms,strucHor)
            call pickPESstruc(indexVer,NAtoms,strucVer)
    
            intervalHor = coordHor(1) - coordCen(1)
            intervalVer = coordVer(2) - coordCen(2)
    
        ! 3. Calculate vector without normalization 
            vecHor = 0.0D0 
            vecVer = 0.0D0 
            do j = 1, 3
                do i = 1, NAtoms 
                    vecHor(i,j) = strucHor(i,j) - strucCen(i,j)
                    vecVer(i,j) = strucVer(i,j) - strucCen(i,j)
                end do 
            end do 
        return 
    end subroutine get_iniNeighbor

    subroutine get_dRMSD(NAtoms,strucTraj,strucCenter,vec,interval,length,dRMSD)
        ! Implement centered finite difference method, accuracy: O(h^2) 
        ! reference: Computational science and engineering I, MIT, G. Strang, Lecture 2
        ! https://ocw.mit.edu/courses/mathematics/18-085-computational-science-and-engineering-i-fall-2008/video-lectures/lecture-2-differential-eqns-and-difference-eqns/
        implicit none
        integer(4),intent(in)   :: NAtoms
        real(8),dimension(NAtoms,3),intent(in)  :: strucTraj, strucCenter, vec 
        real(8),intent(in)  :: interval, length
        real(8),intent(out) :: dRMSD
        ! local variables 
        real(8), parameter :: pts = 100.0D0 
        real(8) :: dx  ! interval is the length of one grid, and dx is length / pts 
        real(8) :: RMSD_F, RMSD_B
        real(8),dimension(NAtoms,3) :: strucForward, strucBackward ! strucShift

        dx = interval / pts
        ! coordShift = 0.0D0 
        strucForward = 0.0D0 
        strucBackward = 0.0D0 
        ! call get_CoordShift(NAtoms,strucCenter,vec,length,strucShift)

        call get_CoordShift(NAtoms,strucCenter,vec,length + dx,strucForward)
        call get_CoordShift(NAtoms,strucCenter,vec,length - dx,strucBackward)

        call get_RMSD(NAtoms,strucTraj,strucForward,RMSD_F)
        call get_RMSD(NAtoms,strucTraj,strucBackward,RMSD_B)
        ! dRMSD = RMSD_F
        dRMSD = ( RMSD_F - RMSD_B ) / (2 * dx)

        return 
    end subroutine get_dRMSD

    subroutine get_1Dunimodal(XorY,zero, NAtoms,numPts,PESList,strucTraj,indexCen,indexNei,vec,interval,df2)
        implicit none 
        integer(4),intent(in) :: XorY, NAtoms, numPts  ! x: XorY=0, y: XorY=1
        real(8), intent(in)  :: zero  
        real(8),dimension(numPts,3),intent(in) :: PESList
        real(8),dimension(NAtoms,3),intent(in) :: strucTraj   
        integer(4),intent(inout)    :: indexCen
        integer(4),intent(out)  :: indexNei
        real(8),dimension(NAtoms,3),intent(out) :: vec
        real(8),intent(out) :: interval
        real(8),intent(inout)   :: df2
        ! local variables 
        integer(4)  :: i, indexHor, indexVer
        real(8), parameter  :: moveL = 0.5D0
        real(8) :: length, intervalHor, intervalVer
        real(8) :: df1, Prod1, Prod2
        real(8),dimension(NAtoms,3) :: strucCen, strucHor, strucVer
        real(8),dimension(NAtoms,3) :: vecHor, vecVer
        ! write(*,*) strucTraj(1,:)
        ! 1. Calculate dRMSD at the initial point
            call get_iniNeighbor(NAtoms,numPts,PESList,strucTraj,indexCen,indexHor,indexVer,intervalHor,intervalVer,vecHor,vecVer)

            ! output: strucCen
            call pickPESstruc(indexCen,NAtoms,strucCen)

            ! output: vec, interval
            if ( XorY .eq. 0 ) then  ! Horizental direction; XorY = 0 
                vec = vecHor
                interval = intervalHor
            else ! Vertical direction; XorY = 1
                vec = vecVer
                interval = intervalHor
            end if
            
            ! output: df1 = dRMSD
            ! get_dRMSD() calculates the new structure inside the subroutine 
            length = 0.0D0 
            call get_dRMSD(NAtoms,strucTraj,strucCen,vec,interval,length,df1) 

        ! 2. Find unimodal region 
        ! If the product of two derivatives is positive, there is no extreme in between, 
        ! otherwise, if it is negative, an extreme exists. 

            ! 2019/09/10, Grace, If the trajectory is exactly on this surface, 
            !   it should stop at the first while loop. Such that, set df2 = df1. 
            df2 = df1 
           
            Prod1 = df1 * df2 
            Prod2 = Prod1
            i = 0 
            do while ( Prod2 .gt. zero )
                i = i + 1 
                ! write(*,*) 'inside1Dunimodal',i,prod2
                if ( i .ge. numPts ) exit !2019/09/19, Grace
                df1 = df2 
                
                ! Calculate the new $length 
                ! 2019/09/05, Grace, change $vec since this PES is not linear,
                ! and the range of length is 
                ! positive direction: 0.0D0  <= $length <= 1.0D
                ! negative direction: -1.0D0 <= $length <= 0.0D0
                ! If the variable, $length, is not stay in above region, change the relative variables.
                ! If length changes, the relative variables need to replace as well. 
                if (  abs(int(length)) .eq. 0  ) then 
                    ! move in the same direction 
                    call get_length(indexCen,indexHor,indexVer,XorY,moveL,length)
                else ! Move point to extent region, desire ouptut: length, strucCen, vec, interval 

                    ! move the original point to the select one
                    if ( XorY .eq. 0 ) then ! Horizental direction; XorY = 0 
                        indexCen = indexHor
                    else ! Vertical direction; XorY = 1
                        indexCen = indexVer
                    end if

                    ! output: indexHor, indexVer, intervalHor, intervalVer, vecHor, vecVer
                    call get_iniNeighbor(NAtoms,numPts,PESList,strucTraj,indexCen,indexHor,indexVer,intervalHor,intervalVer,vecHor,vecVer)

                    ! output: coord, vec, interval
                    if ( XorY .eq. 0 ) then ! Horizental direction; XorY = 0 
                        vec = vecHor
                        interval = intervalHor
                    else ! Vertical direction; XorY = 1
                        vec = vecVer 
                        interval = intervalVer
                    end if

                    ! output: strucCen
                    call pickPESstruc(indexCen,NAtoms,strucCen)

                    ! output: length 
                    length = 0.0D0
                end if
                
                ! Desire new parameters to calculate a new dRMSD: 
                ! 1. $strucCen, 2. $vec, 3. $interval, 4. $length
                ! output: df2 = dRMSD/dx 
                call get_dRMSD(NAtoms,strucTraj,strucCen,vec,interval,length,df2)

                Prod2 = df1 * df2 
                ! write(*,*) Prod2, length 
            end do

            if ( XorY .eq. 0 ) then 
                ! x-direction
                indexNei = indexHor
            else 
                ! y-direction
                indexNei = indexVer
            end if
        return 
    end subroutine get_1Dunimodal

    subroutine get_2Dunimodal(zero,NAtoms,numPts,PESList,strucTraj,indexCen,vecHor,vecVer,intervalHor,intervalVer,prodHor,prodVer)
        implicit none 
        integer(4),intent(in)   :: NAtoms, numPts
        real(8),intent(in)  :: zero 
        real(8),dimension(numPts,3),intent(in) :: PESList
        real(8),dimension(NAtoms,3),intent(in) :: strucTraj   
        integer(4),intent(inout)    :: indexCen 
        real(8),intent(inout)   :: intervalHor, intervalVer, prodHor, prodVer  
        real(8),dimension(NAtoms,3),intent(inout) :: vecHor, vecVer
        ! local variable
        integer(4)   :: i, XorY, indexHor, indexVer, indexCentmp
        real(8) :: length,  df1, dfHor, dfVer
        real(8),dimension(NAtoms,3) :: strucCen
        
            
            XorY = 0
            call get_1Dunimodal(XorY,zero, NAtoms,numPts,PESList,strucTraj,indexCen,indexHor,vecHor,intervalHor,dfHor)

            indexCentmp = indexCen 
            XorY = 1 
            call get_1Dunimodal(XorY,zero, NAtoms,numPts,PESList,strucTraj,indexCentmp,indexVer,vecVer,intervalVer,dfVer)

            i = 0 
            do while ( indexCen .ne. indexCentmp )
                i = i + 1 
                if ( i .eq. 10 ) exit ! 2019/10/03, Grace
                indexCen = indexCentmp 
                XorY = 0
                call get_1Dunimodal(XorY,zero, NAtoms,numPts,PESList,strucTraj,indexCen,indexHor,vecHor,intervalHor,dfHor)

                indexCentmp = indexCen 
                XorY = 1 
                call get_1Dunimodal(XorY,zero, NAtoms,numPts,PESList,strucTraj,indexCentmp,indexVer,vecVer,intervalVer,dfVer)
                ! write(*,*) indexCen,indexCentmp
            end do 

            call pickPESstruc(indexCen,NAtoms,strucCen)
            length = 0.0D0 
            call get_dRMSD(NAtoms,strucTraj,strucCen,vecHor,intervalHor,length,df1) 

            ! length = 1.0D0 
            ! call get_dRMSD(NAtoms,strucTraj,strucCen,vecHor,intervalHor,length,df2) 
            prodHor = df1 * dfHor 

            ! call get_dRMSD(NAtoms,strucTraj,strucCen,vecVer,intervalVer,length,df2) 
            prodVer = df1 * dfVer
            ! write(*,*) '1Dunimodal',indexCen,df1,dfHor,dfVer
        return 
    end subroutine get_2Dunimodal 

    subroutine get_1DGSS(XorY,NAtoms,numPts,PESList,indexCen,indexNei,vec,interval,strucTraj,coord)
        ! 2019/09/10, Grace, Golden section search
        ! ref: https://www.essie.ufl.edu/~kgurl/Classes/Lect3421/NM6_optim_s02.pdf
        ! Golden ratio: https://en.wikipedia.org/wiki/Golden_ratio
        implicit none
        integer(4),intent(in)   :: XorY, NAtoms, numPts, indexCen, indexNei
        real(8),dimension(NAtoms,3),intent(in)  :: strucTraj,vec
        real(8),dimension(numPts,3),intent(in)  :: PESList
        real(8),intent(in)  :: interval
        real(8),intent(out)  :: coord
        ! local variables,
        real(8) :: x1, x2, rmsdX1,rmsdX2,err, xmin, xmax, length, coordIni, coordNei
        real(8) :: intervalHor, intervalVer
        real(8),parameter   :: R = 0.6183D0, C = 1.0D0 - R ! Golden ratio
        real(8),parameter   :: tol = 0.00001D0 !  tolerance 
        real(8),dimension(NAtoms,3) :: strucCen, strucNei
        real(8),dimension(NAtoms,3) :: strucX1, strucX2, vecHor, vecVer
        integer(4)  :: indexHor, indexVer, i, j 

        ! 1. 
        call pickPESstruc(indexCen,NAtoms,strucCen)

        if ( XorY .eq. 0 ) then  
            ! x coordinate; XorY = 0 
            coordIni = PESList(indexCen,1)
        else 
            ! y coordinate; XorY = 1 
            coordIni = PESList(indexCen,2)
        end if

        xmin = 0.0D0 
        xmax = 1.0D0 
        err = C * ( xmax - xmin )
        do while ( err .gt. tol )
            ! 1. Deivde interval into 3 sections by adding two internal points between ends
            x1 = xmax - R * ( xmax - xmin )
            x2 = xmin + R * ( xmax - xmin )

            ! 2. Evalute the function at the two internal points
            call get_CoordShift(NAtoms,strucCen,vec,x1,strucX1)
            call get_CoordShift(NAtoms,strucCen,vec,x2,strucX2)

            call get_RMSD(NAtoms,strucTraj,strucX1,rmsdX1)
            call get_RMSD(NAtoms,strucTraj,strucX2,rmsdX2)

            if ( rmsdX1 .gt. rmsdX2 ) then
                xmax = x2
            else
                xmin = x1 
            end if 
            err = C * ( xmax - xmin )
            ! write(*,'(4(F5.2,1X))') x1, x2, err, tol 
        end do 

        length = ( x2 - x1 ) / 2.0D0 

        coord = coordIni + interval * length 
        ! write(*,*) coord, coordIni, interval, length 
        ! write(*,*) length 
        ! do i = 0, 10 
        !     x1 = i * 1.0/10.0D0 
        !     call get_CoordShift(NAtoms,strucCen,vec,x1,strucX1)
        !     call get_RMSD(NAtoms,strucTraj,strucX1,rmsdX1)
        !     write(*,*) x1, rmsdX1 
        ! end do 
        ! write(*,'()') 
        return 
    end subroutine get_1DGSS 
end module extractData

Program main
    use gloable
    implicit none
    integer(4)  :: tss_index    ! index of tss which can search the 
                                ! corresponding coordinate and energy 

    ! Step 1. I/O import trajectory and numerical PES 
        call print_purpose()    
        call get_NAtoms(NAtoms)
        call get_numPts(NAtoms,numPts)
        allocate(PESList(numPts,3))
        call get_PESList(numPts,PESList)
        
    ! Step 2. Find the initial coordinate of the first point; TSS.
        call findFirstPt(NAtoms,numPts,tss_index) ! output: tss_index

    ! Step 3. Compare the RMSD to its neighborhood, and find the next point.
    ! Also use tuneTraj2() for the first point and all the following points in 
    ! the subroutine coordHalfTraj().
    ! Output the corresponding coordinate in findOtherPts().

        call findOtherPts(NAtoms,numPts,PESList,tss_index)

        deallocate(PESList)
        
    stop
End program main


subroutine findFirstPt(NAtoms,numPts,tss_index)
    use extractData
    implicit none
    integer(4),intent(in)   :: NAtoms, numPts
    integer(4),intent(out)  :: tss_index 
    ! local variables 
    integer(4)  :: i
    real(8) :: rmsd1, rmsd2
    real(8),dimension(NAtoms,3) :: tssTraj, coordPES

    ! 1. Extract the tss of the trajectory; the 1st point of trajectory
        call pickTssTraj(NAtoms,tssTraj)

    ! 2. Search the corresponding structure on this PES 
        ! TODO: uncomment the do loop 
        ! rmsd1 = 10000.D0 ! to be replaced
        ! do i = 1, numPts
        !     call pickPESstruc(i,NAtoms,coordPES)
        !     call get_RMSD(NAtoms,tssTraj,coordPES,rmsd2)
        !     if ( rmsd1 .gt. rmsd2 ) then 
        !         rmsd1 = rmsd2
        !         tss_index = i
        !     end if 
        ! end do 
        tss_index = 1741 ! 1741 !761 !2381
    return 
end subroutine findFirstPt

subroutine findOtherPts(NAtoms,numPts,PESList,tss_index)
    use extractData
    implicit none
    integer(4),intent(in)   :: NAtoms, numPts,tss_index
    real(8),dimension(numPts,3),intent(in)  :: PESList 
    ! local variables
    character(len=100)  :: inputFile, outputFile, finalOutput, buffer 
    character(len=100)  :: cstr_r, cstr_f ! constraint for reverse/forward directions; defR.dat or defP$n.dat
    integer(4)  :: numLine, i
    

    ! 1. Split trajectory into two parts
        call splitTraj() !output: half1Traj.xyz and half2Traj.xyz 

    ! 2. Assign the final output filename 
        call get_output(finalOutput) ! output: coord.$1 ; $1 = 1st arg. in commend line

    ! 3. Search the coordinate for the first half trajectory, 
    !   and then reverse its order. TSS coordinate is also removed, 
    !   since the later half trajectory already have it.
        inputFile = 'half1Traj.xyz'
        outputFile = 'coord.'//TRIM(ADJUSTL(inputFile))
        call GETARG(4,cstr_r)
        call coordHalfTraj(NAtoms,numPts,PESList,tss_index,cstr_r,inputFile,outputFile)

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
        ! write(*,'()')
        ! write(*,'(A)') 'Half'
        ! write(*,'()')
        
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

    ! 5. Close file ID after finish printing out the coordinate for both backward and 
    !   forward direction. Kill both tmp files as well. 
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
    ! cstr = constraint for reverse/forward directions; defR.dat or defP$n.dat
    ! local variables 
    real(8),dimension(2)    :: tune_coord
    real(8),dimension(NAtoms,3) :: tssTraj, oneStruc
    character(len=100)  :: charline, charline2, filename
    integer(4)  :: numline, oneJob, numTrajPts, inIndex, outIndex
    integer(4)  :: i, iniL, finL

    ! 1.  Count the total amount of points of one trajectory (i.e. numjobs).
        call system("wc -l "//TRIM(ADJUSTL(inputFile)) & 
            //" | awk '{print $1}' > charline.dat")
        open(10,file='charline.dat',status='old')
        read(10,*) numline
        close(10,status='delete')
        oneJob = NAtoms + 2 
        numTrajPts = numline / oneJob

    ! 2. The coordinate of the first point (TSS) is already founded, 
    !   tune this coordinate, and then save it in outputFile
        inIndex = tss_index 
        call pickTssTraj(NAtoms,tssTraj)

        ! Tune the coordinate in both x and y directions
        call tuneTraj(NAtoms,numPts,PESList,tssTraj,inIndex,tune_coord)

        open(20,file=outputFile,status='replace')
        ! write(20,"(3(F15.8,1X))") PESList(inIndex,1:3) ! turn-off tuneTraj2()
        ! write(*,"(3(F15.8,1X))") PESList(inIndex,1:3) ! turn-off tuneTraj2()
        write(20,"(3(F15.8,1X))") tune_coord(1:2),PESList(inIndex,3)

        ! for debug 
            ! write(*,'()')
            ! write(*,"(I4,1X,3(F15.8,1X))") inIndex, PESList(inIndex,1:3)
            ! write(*,"(3(F15.8,1X))") tune_coord(1:2),PESList(inIndex,3)
            ! stop
    
    ! 3. Search the coordinate of the rest  points.
        ! inIndex = tss_index 
        ! write(*,*) TRIM(ADJUSTL(inputFile))
        do i = 2,  numTrajPts 
            iniL = 1 + oneJob * ( i - 1 ) 
            finL = oneJob * i
            write(charLine,*) iniL
            write(charLine2,*) finL 
            
            call system("sed -n '"//TRIM(ADJUSTL(charLine))//',' &
                //TRIM(ADJUSTL(charline2))//" p' " & 
                //TRIM(ADJUSTL(inputFile))// " > oneStruc.xyz")
            filename = 'oneStruc.xyz'
            call readStruc(filename,NAtoms,oneStruc)

            ! Compare the neighest 8 neighbors, select the one with the smallest RMSD
            ! cstr = constraint for reverse/forward directions; defR.dat or defP$n.dat
            call step2_select1fromNeighbors(NAtoms,numPts,oneStruc,cstr,inIndex,outIndex)
            
            ! Tune the coordinate in both x and y directions FIXME:
            ! call tuneTraj(NAtoms,numPts,PESList,oneStruc,outIndex,tune_coord)

            inIndex = outIndex 
            write(20,"(3(F15.8,1X))") PESList(outIndex,1:3) ! turn-off tuneTraj()
            ! write(20,"(3(F15.8,1X))") tune_coord(1:2),PESList(outIndex,3)
            ! write(*,'(A)') '---------------------------'
        end do 
        close(10)
    return 
end subroutine coordHalfTraj

subroutine step2_select1fromNeighbors(NAtoms,numPts,strucTraj,cstr,inIndex,outIndex)
    use extractData
    implicit none 
    integer(4),intent(in)   :: NAtoms, numPts, inIndex
    real(8),dimension(NAtoms,3),intent(in)  :: strucTraj
    character(len=100),intent(in)   :: cstr ! constraint for reverse/forward directions
    integer(4),intent(out)  :: outIndex 
    ! local variables 
    ! numNeighbor = 4, right, left, upper, lower neighbors
    ! numNeighbor = 8, first shell 
    ! numNeighbor = 24, first and second shell 
    ! numNeighbor = 48, first, second and the third shell
    integer(4),parameter :: numNeighbor = 48

    integer(4),dimension(numNeighbor) :: neighbors  
    integer(4)  :: i
    real(8) :: rmsd1, rmsd2 
    real(8),dimension(NAtoms,3) :: strucPES

    call vaildNeighbors(inIndex,numPts,cstr,numNeighbor,neighbors)
    ! write(*,'(9(I4,1X))') inIndex,neighbors

    rmsd1 = 10000.D0 ! to be replaced
    do i = 1, numNeighbor
        if ( neighbors(i) .eq. 0 ) cycle
        call pickPESstruc(neighbors(i),NAtoms,strucPES)
        call get_RMSD(NAtoms,strucTraj,strucPES,rmsd2)
        ! write(*,*) rmsd2
        if ( rmsd1 .gt. rmsd2 ) then 
            rmsd1 = rmsd2
            outIndex = neighbors(i)
        end if 
        ! write(*,*) i,inIndex, outIndex
        ! write(*,*) neighbors(:)
    end do 
    return 
end subroutine step2_select1fromNeighbors

subroutine tuneTraj(NAtoms,numPts,PESList,strucTraj,indexCen,coord)
    use extractData
    implicit none
    integer(4),intent(in)   :: NAtoms, numPts
    real(8),dimension(numPts,3),intent(in) :: PESList
    real(8),dimension(NAtoms,3),intent(in) :: strucTraj   
    integer(4),intent(inout)    :: indexCen
    real(8),dimension(2),intent(out)    :: coord
    ! local variables 
    integer(4)  :: i, XorY, indexHor, indexVer, indexCentmp
    real(8) :: prodHor, prodVer,intervalHor, intervalVer
    real(8),dimension(NAtoms,3) :: vecHor, vecVer
    real(8),parameter  :: zero = 0.00001D0 ! do not change this parameter

        ! Calculate new coordinate if its derivative of RMSD is near to zero.
        ! 1. Locate a unimodal area; constrain a point of trajectory in 
        !   certain region which has a local minimum. 

            ! initialize 
            call get_2Dunimodal(zero,NAtoms,numPts,PESList,strucTraj,indexCen,vecHor,vecVer,intervalHor,intervalVer,prodHor,prodVer)

            i = 0 
            do while ( ( prodHor .gt. zero ) .or. (prodVer .gt. zero ) )
                i = i + 1
                call get_2Dunimodal(zero,NAtoms,numPts,PESList,strucTraj,indexCen,vecHor,vecVer,intervalHor,intervalVer,prodHor,prodVer)
                if ( i .eq. 10 ) exit ! 2019/10/03, Grace
                ! write(*,*) '2Dunimodal',i,prodHor, prodVer
            end do

            ! write(*,*) indexCen 
            ! write(*,*) prodHor, prodVer
            
            ! coord(1) = PESList(indexCen,1)
            ! coord(2) = PESList(indexCen,2)
        ! 2. Golden-section search 
            coord = 0.0D0 
            XorY = 0
            call get_1DGSS(XorY,NAtoms,numPts,PESList,indexCen,indexHor,vecHor,intervalHor,strucTraj,coord(1))
            XorY = 1 
            call get_1DGSS(XorY,NAtoms,numPts,PESList,indexCen,indexVer,vecVer,intervalVer,strucTraj,coord(2))

    return 
end subroutine tuneTraj