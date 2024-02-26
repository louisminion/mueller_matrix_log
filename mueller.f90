! Reads a csv from B23 Mueller polarimeter and performs Mueller matrix decomposition on each matrix in the array.
! There is a suspected bug in the matrix log algorithm due to mismatch with python program results.
! © Louis Minion 2023. Released under GNU GPL V3.
program mueller
    implicit none
    ! Variable declaration
    double precision :: G(4, 4) ! G = diag([1, -1, -1, -1]) Minkowski metric
    double precision, allocatable, dimension( :, :, :) :: expmm, lu, lm, l
    double precision, allocatable, dimension(:, :) :: expmmi, lui, lmi, li
    external :: dgeev, DGETRF, DGETRI
    intrinsic :: min, int, log, transpose, matmul
    ! double precision :: wavelengths, As, LDs, 
    Integer :: open_error, i,j, nlines, nargs
    character(len=256) :: filename, outfilename
    ! Load matrix into memory
    ! Split matrix by wavelength
    ! For each wavelength, take the mueller matrix, find the logarithm
    ! Split the differential matrix into symmetric and asymmetric parts by mat multiplication
    ! Append an array for each of the variables extracted
    type :: Experimental_Mueller_Matrix
     integer :: record_n
     double precision :: wavelength
     double precision :: x_pos
     double precision :: y_pos
     character(19) :: time
     double precision :: m_00
     double precision :: m_01
     double precision :: m_02
     double precision :: m_03
     double precision :: m_10
     double precision :: m_11
     double precision :: m_12
     double precision :: m_13
     double precision :: m_20
     double precision :: m_21
     double precision :: m_22
     double precision :: m_23
     double precision :: m_30
     double precision :: m_31
     double precision :: m_32
     double precision :: m_33
     double precision :: A
     double precision :: LDabs
     double precision :: LDang
     double precision :: CD
     double precision :: LRabs
     double precision :: LRang
     double precision :: CR
     double precision :: LD
     double precision :: LDp
     double precision :: LR
     double precision :: LRp
     double precision :: TR
     double precision :: AA
     double precision :: LDabsA
     double precision :: LDangA
     double precision :: CDA
     double precision :: LRabsA
     double precision :: LRangA
     double precision :: CRA
     double precision :: LDA
     double precision :: LDpA
     double precision :: LRA
     double precision :: LRpA
     double precision :: TRA
     double precision :: N
     double precision :: NAngle
     double precision :: S
     double precision :: SAngle
     double precision :: C
     double precision :: Mraw00
     double precision :: Mraw01
     double precision :: Mraw02
     double precision :: Mraw03
     double precision :: Mraw10
     double precision :: Mraw11
     double precision :: Mraw12
     double precision :: Mraw13
     double precision :: Mraw20
     double precision :: Mraw21
     double precision :: Mraw22
     double precision :: Mraw23
     double precision :: Mraw30
     double precision :: Mraw31
     double precision :: Mraw32
     double precision :: Mraw33
     double precision :: MSTD00
     double precision :: MSTD01
     double precision :: MSTD02
     double precision :: MSTD03
     double precision :: MSTD10
     double precision :: MSTD11
     double precision :: MSTD12
     double precision :: MSTD13
     double precision :: MSTD20
     double precision :: MSTD21
     double precision :: MSTD22
     double precision :: MSTD23
     double precision :: MSTD30
     double precision :: MSTD31
     double precision :: MSTD32
     double precision :: MSTD33
     integer :: WavelengthIndex
     integer :: XIndex
     integer :: YIndex
     integer :: TemperatureIndex
     integer :: RepeatIndex
     double precision :: Temp
    end type Experimental_Mueller_Matrix
    type(Experimental_Mueller_Matrix), allocatable, dimension(:) :: measured_matrices
    nargs = command_argument_count()
    if (nargs == 0) then
        Write(*,*) 'No filename entered. Enter filename below: '
        Read(*,*) filename
    else
        call GET_COMMAND_ARGUMENT(1,filename)
    end if
    filename = TRIM(filename)
    j = LEN_TRIM(filename)
    outfilename = filename(1:j-4) // '_LOGDECOMP.csv'
    G = reshape((/double precision :: 1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,-1/), (/4,4/))
    nlines = 0
    open(1, file = filename)
    do
        read(1,*, end=10)
        nlines = nlines + 1
    end do
    10 close(1)
    write(*,*) "Number of matrices in file:",nlines-17
    allocate( measured_matrices(nlines-17))
    Open(Unit=12, File=filename, Action='Read', Iostat=open_error, status='old')
    If (open_error /= 0) Stop 'Error opening data file'
    Write(*,*) 'Reading from input data...'
    j = 0 
    do i= 1, 100000
        if (j<17) then
         read(12, *)
         j = j +1
         cycle
        end if
        j = j + 1
        read(12, *, end=100) measured_matrices(i-17)
        ! write(*,*) measured_matrices(i-17)%wavelength
    end do
    100 Close(12)
    allocate(expmm(4,4,nlines))
    do i=1,nlines
        expmm(1,1,i) = measured_matrices(i)%m_00
        expmm(1,2,i) = measured_matrices(i)%m_01
        expmm(1,3,i) = measured_matrices(i)%m_02
        expmm(1,4,i) = measured_matrices(i)%m_03
        expmm(2,1,i) = measured_matrices(i)%m_10
        expmm(2,2,i) = measured_matrices(i)%m_11
        expmm(2,3,i) = measured_matrices(i)%m_12
        expmm(2,4,i) = measured_matrices(i)%m_13
        expmm(3,1,i) = measured_matrices(i)%m_20
        expmm(3,2,i) = measured_matrices(i)%m_21
        expmm(3,3,i) = measured_matrices(i)%m_22
        expmm(3,4,i) = measured_matrices(i)%m_23
        expmm(4,1,i) = measured_matrices(i)%m_30
        expmm(4,2,i) = measured_matrices(i)%m_31
        expmm(4,3,i) = measured_matrices(i)%m_32
        expmm(4,4,i) = measured_matrices(i)%m_33
    end do
    allocate(expmmi(4,4))
    allocate(l(4,4,nlines-17))
    allocate(li(4,4))
    allocate(lui(4,4))
    allocate(lmi(4,4))
    allocate(lu(4,4,nlines-17))
    allocate(lm(4,4,nlines-17))
    do i=1, nlines-17
        expmmi = expmm(:,:,i)
        call CALC_LOGM(expmmi,li)
        ! write(*,*) li
        l(:,:,i) = li
        ! Lm = 1/2(L-GL^TG), Lu = 1/2(L+GL^TG)
        lmi = 0.5*(li-matmul(G,matmul(transpose(li),G)))
        lui = 0.5*(li+matmul(G,matmul(transpose(li),G)))
        lm(:,:,i) = lmi
        lu(:,:,i) = lui
    end do
    ! Write data out
    open(unit=15,file=outfilename,Action='Write', Iostat=open_error)
    If (open_error /= 0) Stop 'Error opening output file'
    101 format(1x,*(g0, :, ", "))
    do i=1,nlines-17
     if (i==1) then
      Write(15,101) "LM00","LM01","LM02","LM03","LM10","LM11", &
      "LM12","LM13","LM20","LM21","LM22","LM23","LM30","LM31","LM32","LM33", &
      "LU00","LU01","LU02","LU03","LU10","LU11","LU12","LU13","LU20","LU21", &
      "LU22","LU23","LU30","LU31","LU32","LU33","L00","L01","L02","L03", &
      "L10","L11","L12","L13","L20","L21","L22","L23","L30","L31","L32","L33"
     end if
        Write(15,101) lm(1,1,i) , lm(1,2,i), lm(1,3,i), lm(1,4,i), &
        lm(2,1,i) , lm(2,2,i), lm(2,3,i), lm(2,4,i), &
        lm(3,1,i) , lm(3,2,i), lm(3,3,i), lm(3,4,i), &
        lm(4,1,i) , lm(4,2,i), lm(4,3,i), lm(4,4,i), &
        lu(1,1,i) , lu(1,2,i), lu(1,3,i), lu(1,4,i), &
        lu(2,1,i) , lu(2,2,i), lu(2,3,i), lu(2,4,i), &
        lu(3,1,i) , lu(3,2,i), lu(3,3,i), lu(3,4,i), &
        lu(4,1,i) , lu(4,2,i), lu(4,3,i), lu(4,4,i), &
        l(1,1,i) , l(1,2,i), l(1,3,i), l(1,4,i), &
        l(2,1,i) , l(2,2,i), l(2,3,i), l(2,4,i), &
        l(3,1,i) , l(3,2,i), l(3,3,i), l(3,4,i), &
        l(4,1,i) , l(4,2,i), l(4,3,i), l(4,4,i)
    end do
    Close(15)
        

end program mueller

SUBROUTINE CALC_LOGM(M, LOGM)
    ! M is an input square matrix to find the logarithm of. Must be diagonizable.
    IMPLICIT NONE
    integer :: n,i,j
    Integer :: LWORK, LWMAX, INFO
    Parameter (LWMAX = 1000)
    Parameter (n=4)
    double precision WORK(LWMAX)
    double precision,dimension(n)::WR,WI,ipiv
    double precision, dimension(n,n) :: M, LOGM, VL, VR, logADASH, VRINV
    ! external :: dgeev, dgetri, dgetrf

    LWORK = -1
    ! OPTIMAL WORKSPACE QUERY
    call dgeev('N','V',n,M,n, WR, WI, VL, n, VR, n,WORK,LWORK,INFO)
    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
    call dgeev('N','V',n,M,n, WR, WI, VL, n, VR, n,WORK,LWORK,INFO)
    if( INFO.GT.0 ) then
        Write(*,*)'The algorithm failed to compute eigenvalues.'
        GOTO 57
    end if
    ! Check eigenvalues >0
    do i=1,n
        if (WR(i) .LE. 0.d0) then
            Write(*,*) 'Eigenvalue less than zero, cannot compute '
            GOTO 57
        end if
    end do

    ! initialise array as zero
    do i=1,n
        do j=1,n
            logADASH(i,j) = 0.d0
        end do
    end do

    do i =1,n
        logADASH(i,i) = LOG(WR(i))
    end do

    VRINV = VR
    call DGETRF(n,n,VRINV,n, ipiv, INFO)
    if (info .NE. 0) then
        Write(*,*) 'LU factorisation failed exit code', INFO
        GOTO 57
    end if
    call DGETRI( n, VRINV, n, ipiv, WORK, LWORK, INFO )
    if (info .NE. 0) then
        Write(*,*) 'LU factorisation failed exit code', INFO
        GOTO 57
    end if
    logM = matmul(matmul(VR,logADASH), VRINV)
57 END SUBROUTINE CALC_LOGM
