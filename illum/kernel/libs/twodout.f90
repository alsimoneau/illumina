subroutine twodout(nbx,nby,filename,bindata) 
! write double precision array in binary 
    integer :: nbx,nby,i,j 
    real :: bindata(512,512) 
    character*72 :: filename 
    open(unit=1,form='unformatted',file=filename,action = 'write') 
    write(1) nbx,nby 
    do j = nby,1,-1 
      do i = 1,nbx 
        write(1) bindata(i,j) 
      end do 
    end do 
    close(unit = 1) 
    return 
end subroutine twodout 
