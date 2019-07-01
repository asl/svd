!                                                                       
!     (C) Rasmus Munk Larsen, Stanford University, 2000                 
!                                                                       
      subroutine clearstat 
      implicit none 
      include 'stat.h' 
      nopx = 0 
      nreorth = 0 
      ndot = 0 
      nitref = 0 
      nbsvd = 0 
      nrestart = 0 
      tmvopx = 0 
      tgetu0 = 0 
      tupdmu = 0 
      tupdnu = 0 
      tintv = 0 
      tlanbpro = 0 
      treorth = 0 
      treorthu = 0 
      treorthv = 0 
      telru = 0 
      telrv = 0 
      tbsvd = 0 
      tnorm2 = 0 
      tdot = 0 
      tlansvd = 0 
      nlandim = 0 
      nsing = 0 
      tritzvec = 0 
      trestart = 0 
      END                                           

      module cfuncs
        use iso_c_binding, only: c_char, c_null_char
        interface
           subroutine printchar0(label) bind(C, name = 'printchar0')
             use iso_c_binding, only: c_char
             character(kind = c_char):: label(*)
           end subroutine printchar0

           subroutine printint0(label, data) bind(C, name = 'printint0')
             use iso_c_binding, only: c_char, c_int
             character(kind = c_char):: label(*)
             integer(kind = c_int), value :: data
           end subroutine printint0

           subroutine printdbl0(label, data) bind(C, name = 'printdbl0')
             use iso_c_binding, only: c_char, c_double
             character(kind = c_char):: label(*)
             real(kind = c_double), value :: data
           end subroutine printdbl0
        end interface
      end module cfuncs

      subroutine printchar(label)
        use cfuncs
        implicit none
        character(*) label 
        call printchar0(label//c_null_char) 
      END                                           
                                                                        
      subroutine printint(label, data) 
        use cfuncs
        implicit none 
        character(*) label 
        integer data 
        call printint0(label//c_null_char, data) 
      END                                           
                                                                        
      subroutine printdbl(label, data) 
        use cfuncs
        implicit none 
        character(*) label 
        real data 
        double precision outd
        outd = data 
        call printdbl0(label//c_null_char, outd) 
      END                                           
                                                                        
      subroutine printstat 
      implicit none 
                                                                        
      include 'stat.h' 
                                                                        
      call printchar(' +------------------------------------------------&
     &-----------+')                                                    
      call printint(' Dimension of Lanczos basis                  = ',  &
     &     nlandim)                                                     
      call printint(' Number of singular values requested         = ',  &
     &     nsing)                                                       
      call printint(' Number of restarts                          = ',  &
     &     nrestart)                                                    
      call printint(' Number of matrix-vector multiplications     = ',  &
     &     nopx)                                                        
      call printint(' Number of reorthogonalizations              = ',  &
     &     nreorth)                                                     
      call printint(' Number of inner products in reorth.         = ',  &
     &     ndot)                                                        
!      print *,'Number of iterative refinement steps        = ',nitref  
      call printint(' Number of bidiagonal SVDs calculated        = ',  &
     &     nbsvd)                                                       
      call printchar('') 
      call printchar('') 
                                                                        
      call printdbl('  Time spent doing matrix-vector multiply    = ',  &
     &     tmvopx)                                                      
      call printdbl('  Time spent generating starting vectors     = ',  &
     &     tgetu0)                                                      
      call printdbl('    Time spent reorthogonalizing U_{j+1}     = ',  &
     &     treorthu)                                                    
      call printdbl('    Time spent reorthogonalizing V_{j}       = ',  &
     &     treorthv)                                                    
      call printdbl('  Time spent reorthogonalizing               = ',  &
     &     treorth)                                                     
      call printdbl(' Total Time spent in LANBPRO                 = ',  &
     &     tlanbpro)                                                    
                                                                        
!      print *                                                          
!      print *,'Time spent updating mu-recurrence           = ',tupdmu  
!      print *,'Time spent updating nu-recurrence           = ',tupdnu  
!      print *,'Time spent on local reorth. on U_{j+1}      = ',telru   
!      print *,'Time spent on local reorth. on V_{j+1}      = ',telrv   
!      print *,'Time spent in PDNORM2                       = ',tnorm2  
!      print *,'Time spent in PDDOT                         = ',tdot    
      call printchar('') 
      call printdbl('  Time spent in LANBPRO                      = ',  &
     &     tlanbpro)                                                    
      call printdbl('  Time spent computing bidiagonal SVDs       = ',  &
     &     tbsvd)                                                       
      call printdbl('  Time spent doing implicit restarts         = ',  &
     &     trestart)                                                    
      call printdbl('  Time spent computing Ritz vectors          = ',  &
     &     tritzvec)                                                    
      call printchar('') 
      call printdbl(' Total Time spent in LANSVD                  = ',  &
     &     tlansvd)                                                     
      call printchar(' +------------------------------------------------&
     &----------+')                                                     
      END                                           
