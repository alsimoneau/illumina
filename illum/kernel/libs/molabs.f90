!  adding molecular absorption 
 
subroutine molabs(lambda,bandw,tabs) 
    real :: lambda,bandw,tabs,wmax,wmin,wla,ta,bwa,bwabs 
    wmin = lambda-bandw/2. 
    wmax = lambda+bandw/2. 
    tabs = 0. 
    wla = 0. 
    bwabs = 0. 
    open(unit=1,file='MolecularAbs.txt',status = 'old') 
    read(1,*) 
    do while (wla < wmax) 
      read(1,*) wla,ta,bwa 
      if (wla > wmin) then 
        tabs = tabs+ta*bwa 
        bwabs = bwabs+bwa 
      end if 
    end do 
    if (bwabs > 0.) tabs = tabs/bwabs 
    close(unit = 1) 
    return 
end subroutine molabs 
