module module_read_write

  use mpi
  use module_precision
  contains

    subroutine read_dbl_from_file(filename, dbl)
      implicit none
      ! user defined parameter structure
      real(kind=rk), intent(out)      :: dbl
      ! inifile name
      character(len=*), intent(in)    :: filename
      integer(kind=ik)                :: io_error

      open (unit=99, file=filename, status='old', action='read', iostat=io_error)
      if ( io_error == 0) then
        read(99, *) dbl
      else
        write(*,*) 'error', io_error,' while opening the file', filename
      end if
      close(99)
    end subroutine read_dbl_from_file

    subroutine read_vec_from_file(filename, vector, L)
      implicit none

      integer(kind=ik), intent(in)                           :: L
      ! user defined parameter structure
      real(kind=rk), intent(out), dimension(:), allocatable  :: vector
      ! inifile name
      character(len=*), intent(in)                           :: filename
      integer(kind=ik)                                       :: n, io_error, stat
      real(kind=rk)                                          :: x

      open (unit=99, file=filename, status='old', action='read', iostat=io_error)
      if ( io_error == 0) then
        n    = 0
        stat = 0
        do while(stat == 0)
          n = n + 1
          read(99,*,iostat=stat) x
        enddo
        allocate(vector(n))
        if ( n .ne. L ) then
          write(*,*) 'ERROR! Length of vector in file and parameter struct disagree', n, L
        else
          read(99,*,iostat=stat) vector
        end if
      else
        write(*,*) 'error', io_error,' while opening the file', filename
      end if
      close(99)
    end subroutine read_vec_from_file

    subroutine read_vec(filename, data, Ndata)
      implicit none

      integer(kind=ik), intent(in)                           :: Ndata
      ! user defined parameter structure
      real(kind=rk), intent(out), dimension(:)               :: data
      ! inifile name
      character(len=*), intent(in)                           :: filename
      integer(kind=ik)                                       :: io_error

      open (unit=99, file=filename, status='old', action='read', iostat=io_error)
      if ( io_error == 0) then
        !allocate(data(Ndata))
        read(99,*) data
      else
        write(*,*) 'error', io_error,' while opening the file', filename
      end if
      close(99)
    end subroutine read_vec

    subroutine read_matrix_from_file(filename, matrix, cols, skip)
      implicit none
      ! user defined parameter structure
      real(kind=rk), intent(out), dimension(:,:), allocatable    :: matrix
      ! inifile name
      character(len=*), intent(in)                               :: filename
      integer(kind=ik), intent(in)                               :: cols
      logical, intent(in)                                        :: skip
      real(kind=rk)                                              :: x
      integer(kind=ik)                                           :: n, io_error, stat, k

      open (unit=99, file=filename, status='old', action='read', iostat=io_error)
      if ( io_error == 0) then
        n    = -1
        stat = 0
        do while(stat == 0)
          n = n + 1
          read(99,*,iostat=stat) x
        enddo
        rewind(99)
        if (skip) then
          allocate(matrix(cols,n-1))
          k = -1
          stat = 0
          do while(stat == 0)
            k = k + 1
            if (k==0 .or. k==n) then
              read(99,*,iostat=stat) x
            else
              read(99,*,iostat=stat) matrix(:,k)
            end if
          enddo
        else
          allocate(matrix(cols,n))
          read(99,*) matrix
        end if
      else
        write(*,*) 'error', io_error,' while opening the file', filename
      end if
      close(99)
    end subroutine read_matrix_from_file
    subroutine write_matrix(filename, data)
      implicit none
      character(len=*), intent(in)                           :: filename
      real(kind=rk), intent(in), dimension(:,:)              :: data
      integer(kind=ik)                                       :: io_error,i

      open (unit=99, file=trim(adjustl(filename)), status='replace', action='write', iostat=io_error)
      if ( io_error == 0) then
        do i =1,size(data,2)
          write(99,*) data(:,i)
        end do
      else
        write(*,*) 'error', io_error,' while opening the file', filename
      end if
      close(99)
    end subroutine write_matrix

end module
