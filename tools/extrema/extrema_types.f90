MODULE extrema_types
  ! This module sets the types used in the Fortran 90 modules
  ! follows the example of Healpix and Numerical Recipes
  !
  INTEGER, PARAMETER, public :: i4b = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER, public :: i8b = SELECTED_INT_KIND(16)
  INTEGER, PARAMETER, public :: i2b = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER, public :: i1b = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER, public :: sp  = SELECTED_REAL_KIND(5,30)
  INTEGER, PARAMETER, public :: dp  = SELECTED_REAL_KIND(12,200)
  INTEGER, PARAMETER, public :: lgt = KIND(.TRUE.)
  INTEGER, PARAMETER, public :: spc = KIND((1.0_sp, 1.0_sp))
  INTEGER, PARAMETER, public :: dpc = KIND((1.0_dp, 1.0_dp))
  !
  INTEGER(I8B),  PARAMETER, public :: max_i8b = HUGE(1_i8b)
  INTEGER,       PARAMETER, public :: max_i4b = HUGE(1_i4b)
  INTEGER,       PARAMETER, public :: max_i2b = HUGE(1_i2b)
  INTEGER,       PARAMETER, public :: max_i1b = 127
  REAL(kind=sp), PARAMETER, public :: max_sp  = HUGE(1.0_sp)
  REAL(kind=dp), PARAMETER, public :: max_dp  = HUGE(1.0_dp)

  ! Numerical Constant (Double precision)
  REAL(kind=dp), PARAMETER, public :: QUARTPI=0.785398163397448309615660845819875721049_dp
  REAL(kind=dp), PARAMETER, public :: HALFPI= 1.570796326794896619231321691639751442099_dp
  REAL(kind=dp), PARAMETER, public :: PI    = 3.141592653589793238462643383279502884197_dp
  REAL(kind=dp), PARAMETER, public :: TWOPI = 6.283185307179586476925286766559005768394_dp
  REAL(kind=dp), PARAMETER, public :: FOURPI=12.56637061435917295385057353311801153679_dp
  REAL(kind=dp), PARAMETER, public :: SQRT2 = 1.41421356237309504880168872420969807856967_dp
  REAL(kind=dp), PARAMETER, public :: EULER = 0.5772156649015328606065120900824024310422_dp
  REAL(kind=dp), PARAMETER, public :: SQ4PI_INV = 0.2820947917738781434740397257803862929220_dp
  REAL(kind=dp), PARAMETER, public :: TWOTHIRD = 0.6666666666666666666666666666666666666666_dp

  real(kind=DP), parameter, public :: RAD2DEG = 180.0_DP / PI
  real(kind=DP), parameter, public :: DEG2RAD = PI / 180.0_DP
  real(kind=SP), parameter, public :: hpx_sbadval = -1.6375e30_sp
  real(kind=DP), parameter, public :: hpx_dbadval = -1.6375e30_dp

  ! Maximum length of filenames
  integer, parameter :: filenamelen = 1024

  ! Control structure
  type CND_CNTRL_TYPE
     integer      :: STAT=0,  EXACT=1
     integer      :: NPROC=1
     logical      :: justprint=.false.
     logical      :: normalize=.true.
  end type CND_CNTRL_TYPE

  ! Input formats
  type inputformats
     integer :: mpgrafic=1, plain=0
  end type inputformats

  ! Output container for extrema
  type EXT_DATA
     INTEGER(I8B) :: pix
     REAL(DP)     :: pos(3)
     REAL(DP)     :: eig(3)
     REAL(DP)     :: val
     INTEGER(I4B) :: typ
  end type EXT_DATA

    type EXT_DATA_3D
     INTEGER(I8B) :: pix(3)
     REAL(DP)     :: pos(3)
     REAL(DP)     :: eig(3)
     REAL(DP)     :: val
     INTEGER(I4B) :: typ
  end type EXT_DATA_3D
  TYPE(inputformats), public  :: inputforms
  INTEGER(I8B),       public  :: NPIX

END MODULE extrema_types
