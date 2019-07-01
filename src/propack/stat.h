!                                                                       
!     (C) Rasmus Munk Larsen, Stanford University, 2000                 
!                                                                       
                                                                        
      integer nopx, nreorth, nreorthu, nreorthv, ndot, nitref,          &
     &     nrestart, nbsvd, nlandim, nsing                              
      real    tmvopx, tgetu0, tupdmu, tupdnu, tintv, tlanbpro,          &
     &     treorthu, treorthv, telru, telrv, tbsvd, tnorm2,             &
     &     tlansvd, tritzvec, trestart, treorth, tdot                   
                                                                        
                                                                        
      common/timing/ nopx, nreorth, ndot, nreorthu, nreorthv, nitref,   &
     &     nrestart, nbsvd, tmvopx, tgetu0, tupdmu, tupdnu, tintv,      &
     &     tlanbpro, treorth, treorthu, treorthv, telru, telrv, tbsvd,  &
     &     tnorm2, tlansvd, nlandim, tritzvec, trestart, tdot, nsing    
                                                                        
