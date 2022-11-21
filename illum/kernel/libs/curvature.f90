!  adding earth curvature to the topography 
 
subroutine curvature(distc,hcur) 
    real :: rearth,distc,hcur 
    Rearth = 6371000. 
    hcur = Rearth-sqrt(Rearth**2.-distc**2.) 
    hcur = -1.*hcur 
    return 
end subroutine curvature 
