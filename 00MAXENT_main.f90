!This code do the Monte Carlo simulation to find the most
!uniform high-entropy alloy by the max-entropy theory
!(Entropy 2013, 15, 5536-5548;Computational Materials Science 142 (2018) 332–337))
!Developed by Prof. Shin-Pon Ju 2019/03/25 at UCLA
!
!1. You need to assign the HEA UNIT cell first
!Output file:

!!!test
!!!!
include 'information.f90'
PROGRAM MAXEnt
!    use omp_lib
    USE information
    implicit real*8(a-h,o-z)
    integer nbetter ! better counter
	!!!integer ( kind = 4 ) thread_num          !!!!!!!!!
    real*8 Boltz,r,ini_acceptratio
    integer NN(5)
    !CALL random_seed()
    call init_random_seed()
	!do i=1,10
	!call random_number(r)
	!write(*,*)i,r
	!enddo
    call general_para ! activate all parameters
    call celllist ! get topology of system

    open(101,file='00MaxEntropy_summary.dat',status='unknown')
    open(109,file='00MCkT_summary.dat',status='unknown')
    nadjust = 10 !interval to adjust kT
    ini_acceptratio = 0.7 ! the initial value for acceptration
    kT = 1.e3 ! the initial value for Boltzmann factor

    Entmin = 1.e20 ! initial entropy
    nbetter = 0
   ! call conf_entropy !first evaluation
    do 111 ini=1,inirand
     
        call shuffle
		call conf_entropy
	 	 
        !	do ient =1,natom
        !    write(*,'(a,i5,a,f10.3)')"Atom ",ient," Entropy:",atomentropy(ient)
        !    enddo
        !    pause

        if(mod(ini,1000) .eq. 0) write(*,*)"Initial shuffle iterations: ",ini

        if(confentropy .le. Entmin)then
            nbetter = nbetter + 1
            Entmin = confentropy
            !xmin = x
            !ymin = y
            !zmin = z
			minatkeep = atype
            !atkeep = atype
			katomentropy = atomentropy
            call lmpdata("00Shuffle",nbetter)
             ! thread_num = omp_get_max_threads ( )!!!!!!!!!!!!!!
            write(101,*)"Randomly Shuffle process"
            write(101,*)"The searched configuration No : ",nbetter
			write(101,*)"Shuffle Steps: ",ini			
            write(101,*) "***current normalized configurational entropy: ",confentropy/dble(natom)
            write(101,*)""
            write(*,*) "***current RANDOM searched better configuration No.: ",nbetter
            write(*,*)"Normalized configurational entropy: ",confentropy/dble(natom)
			!write(*,*) '  The number of threads available is:    ', thread_num
            write(*,*)""
 
        endif

    !   call shuffle
	!	call conf_entropy
	        !    call partialshuffle
        !
        111 continue

        write(101,*)"***** The following is MC output"
        write(101,*)""
        ! Use the best configuration after the above random search
        !x = xmin
        !y = ymin
        !z = zmin
		confentropy = Entmin
        atype = minatkeep
        !xkeep = xmin
        !ykeep = ymin
        !zkeep = zmin

        Entkeep = Entmin ! Entkeep for MC, the initial value uses the best Entmin after the above
        atomentropy = katomentropy
        afterconfentropy = confentropy
        nbetter = 0
        naccept = 0

do 112 imc=1,iterMC
            !acceptratio gradually becomes smaller with the increasing MC step, like annealing
            acceptratio = ini_acceptratio - ( dble(imc)/dble(iterMC) )*0.4 !! lowest: 0.5 (from 0.9)
            !if(mod(imc,500) .eq. 0)then
            !  call shuffle
			
			!filterValue = weight(2)*(6-imc+1)
			filterValue = 0.
			!if (filterValue .lt. weight(2))filterValue = weight(2)
!!!!!!! Do while loop             !else



!!!!!!!! end of do while loop            !endif            
!             ishindex = 0
!            if(mod(imc,3) .eq. 0) then
!                      
!			call kshuffle ! pick the highest entropy atoms
!            ishindex = 0
!			elseif(mod(imc,3) .eq. 1) then
!                
!				call normshuffle ! pick one highest and one lowest
!			ishindex = 1	
!          !  elseif(mod(imc,5) .eq. 3) then
!          !      call twoshuffle ! randomly pick two
!            elseif(mod(imc,3) .eq. 3) then
!            
!				call randshuffle ! randomly pick to shuffle
!				ishindex = 2
!          !  else
!          !      call shuffle ! all atoms for shuffle
!            endif
!			
!            call conf_entropy

            !write(*,*)"1"
            !call twoshuffle
            !call kshuffle
            !	call normshuffle
            !write(*,*)"2"
		if(Mod(imc,50).eq.0)then	
		  
	!	   keepNo = 0
	!	  do i11 = 1,natom
	!      if(atomentropy(i11) .ge. filterValue)then
	!	     keepNo = keepNo + 1
	!	
	!	      keepeID(keepNo) = i11
	!	!write(*,*)i,atomentropy(i),keepNo
	!	!pause
	!       endif	
	!!if(keepNo .ge. 200) exit ! only consider the maximal number of 100
    !      enddo
		  
		  
          write(*,*)""
          write(*,*)""		  
          write(*,*)"##********  MC iterations: ",imc
	!	  write(*,*)"##  atom Number with large scoring value:",keepNo
		 
	!	  do ikeep =1,keepNO
	!	    write(*,*)ikeep,keepeID(ikeep),atomentropy(keepeID(ikeep))
	!	  enddo
		  write(*,*)"Ave confentropy:", afterconfentropy/dble(natom)
          write(*,*)"*** CURRENT best:",Entmin/dble(natom)
		  write(*,*) "*** The first three nearest neighbour numbers for current best***"	
          write(*,*) "***N1:",NN(1)/2," ***N2",NN(2)/2," ***N3",NN(3)/2," ***N4",NN(4)/2," ***N5",NN(5)/2			  
		  write(*,*)""
		  !pause
		endif 
	!do i=1,10
	!call random_number(r)
	!write(*,*)"before ten rand",i,r
	!enddo
!**************** by OMP
    	!write(*,*)"before Ave confentropy:", confentropy/dble(natom) 	
		  call SCF ! use each poossile pair to get the best final atype
		afterconfentropy = confentropy
		!write(*,*)"after Ave confentropy:", confentropy/dble(natom)   
	!	  	do i=1,10
	!call random_number(r)
	!write(*,*)iscf,"after ten rand",i,r
	!enddo
    
! output the lowest potential, atype, atomentropy 

            !!! Keep the better one
            if(confentropy .lt. Entmin)then
                nbetter = nbetter + 1
                Entmin = confentropy
                !xmin = x
                !ymin = y
                !zmin = z

                !xkeep = x
                !ykeep = y
                !zkeep = z
				minatkeep = atype
                !atkeep = atype ! keep atom type
				katomentropy = atomentropy ! keep atom potential
!!!!!!!!!!!!!				! get neigbour atom number 
				NN =0
								
				do 1 i=1,natom
				!write(*,*)"atom: ",i," atomentropy: ",atomentropy(i)
					do 2 ine =1,5 !neighbor atom type
					ntemp = CN_No(i,ine)
						do 3 neID = 1,ntemp ! atom ID of a neighbor type
						JID = CN_ID(i,ine,neID)          
						if(atype(i) .eq. atype(JID))then
							NN(ine) = NN(ine) + 1
						endif 
				3 		continue                      
				2   continue  						
				1 continue			
				
                !Entkeep = Entmin ! for MC
                call lmpdata("00MC",nbetter)
                write(101,*)"#MC process at MC step:", imc
                write(101,*)"The searched configuration No : ",nbetter
                write(101,*) "***current normalized configurational entropy: ",confentropy/dble(natom)
                write(101,*) "***N1:",NN(1)/2," ***N2",NN(2)/2," ***N3",NN(3)/2," ***N4",NN(4)/2," ***N5",NN(5)/2	! the corresponding NN atom No.			
                
				!keepNo = 0
		 ! do i11 = 1,natom
	     ! if(atomentropy(i11) .ge. filterValue)then
		 !    keepNo = keepNo + 1
		!
		 !     keepeID(keepNo) = i11
		!!write(*,*)i,atomentropy(i),keepNo
		!!pause
	     !  endif	
	!if(keepNo .ge. 200) exit ! only consider the maximal number of 100
          !enddo
		!	do ikeep =1,keepNO
		  !  write(101,*)ikeep,keepeID(ikeep),atomentropy(keepeID(ikeep))
		  !enddo	
				
				
				write(101,*)""
                write(*,*) "******current MC searched better configuration No.: ",nbetter
				write(*,*) "***current MC setp for this configuration:",imc				
                write(*,*)"Normalized configurational entropy: ",confentropy/dble(natom)
				write(*,*) "***N1:",NN(1)/2," ***N2",NN(2)/2," ***N3",NN(3)/2," ***N4",NN(4)/2," ***N5",NN(5)/2	! the corresponding NN atom No.	
                write(*,*)""

            !else
			!	confentropy = Entmin
			!	atype = minatkeep
			!	atomentropy = katomentropy
            endif
            ! begin Boltzmann evaluation

            Boltz =  Exp(-(confentropy-Entkeep)/(kT*dble(natom)))
            call random_number(r)
            ! write(*,*)"###########Boltz",Boltz,"r",r
			! write(*,*)"Entkeep",Entkeep,"confentropy",confentropy
            if(confentropy .le. Entkeep)then
                !xkeep = x
                !ykeep = y
                !zkeep = z
			!	write(*,*)"******Find a better one without Boltzmann factor at loop: ",imc
			!	write(*,*)"confentropy: ",confentropy/dble(natom),"Entkeep: ",Entkeep/dble(natom)
				!pause
                Entkeep = confentropy
				atkeep = atype
				katomentropy = atomentropy ! keep atom potential
                naccept = naccept + 1
				
         !   elseif(r .lt. Boltz) then
		!		!if (Boltz .lt. 0.01) write(*,*)"*****Using the case within the worst 1% for evolution!!!" 
         !       !xkeep = x
         !       !ykeep = y
         !       !zkeep = z
         !       Entkeep = confentropy
		!		atkeep = atype
		!		katomentropy = atomentropy ! keep atom potential
         !       naccept = naccept + 1
         !    !   write(*,*)"******Boltzmann factor works at loop: ",imc
         !    !   write(*,*)r,Boltz
         !       !call sleep(1)
            else
                !x = xkeep
                !y = ykeep
                !z = zkeep
				confentropy = Entkeep
				atype = atkeep
				atomentropy = katomentropy ! keep atom potential
                
             !   write(*,*)"******Go back to the last structure at loop: ",imc
             !   write(*,*)r,Boltz
                !call sleep(1)
            endif
!pause
           ! adjust kT

            if(mod(imc,nadjust) .eq.0)then

                write(109,*)  "***MC steps: ",imc
                write(109,*)    "Total MC moves: ", nadjust
                write(109,*)     "Accepted MC moves: ",naccept
                temp = dble(naccept)/dble(nadjust)
                write(109,*)"Accepted_ratio: ",temp
                write(109,*)"Desired_ratio: ", acceptratio

                if(temp .lt. acceptratio) then
                    write(109,*)"current_ratio < desired_upper"
                    write(109,*)"Old kT: ", kT
                    kT = kT*2.
                    write(109,*)"Adjusted kT: ", kT
                    write(109,*)""
                elseif(temp .gt. acceptratio) then
                    write(109,*)"current_ratio > desired_upper"
                    write(109,*)"Old kT: ", kT
                    kT = kT/2.
                    write(109,*)"Adjusted kT: ", kT
                    write(109,*)""
                endif
                naccept = 0 ! activate a new counter
            endif
!! generate the next 
            call kshuffle(imc)
            !call lmpdata("00MC",imc)
112 continue

            close(101)
            close(109)
            write(*,*)"MC DONE!!"

END

include 'shuffle.f90' ! all random shuffle
include 'kshuffle.f90' ! pick the highest entropy atoms
include 'normshuffle.f90' ! pick one highest and one lowest
include 'twoshuffle.f90' ! randomly pick two
include 'randshuffle.f90' ! randomly pick two
!include 'partialshuffle.f95'
include 'celllist.f90'
include 'lmpdata.f90'
include 'confentropy.f90'
include 'init_random_seed.f90'
include 'confentropyP.f90'
include 'SCF.f90'
