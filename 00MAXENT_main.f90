!This code do the Monte Carlo simulation to find the most
!uniform high-entropy alloy by the max-entropy theory
!(Entropy 2013, 15, 5536-5548;Computational Materials Science 142 (2018) 332–337))
!Developed by Prof. Shin-Pon Ju 2019/03/25 at UCLA
!major modification by Prof. Shin-Pon Ju 2020/03/12 at NSYSU
!
include 'information.f90'
PROGRAM MAXEnt
!    use omp_lib
    USE information
    implicit real*8(a-h,o-z)
    integer nbetter ! better counter
	!!!integer ( kind = 4 ) thread_num          !!!!!!!!!
    real*8 Boltz,r,ini_acceptratio
    integer NN(5)
    CALL random_seed()
    !call init_random_seed(123)
	call general_para ! activate all parameters
    call celllist ! get topology of system

    open(101,file='00MaxEntropy_summary.dat',status='unknown')
    open(109,file='00MCkT_summary.dat',status='unknown')
    nadjust = 10 !interval to adjust kT
    ini_acceptratio = 0.7 ! the initial value for acceptration
    kT = 1.e3 ! the initial value for Boltzmann factor

    Entmin = 1.e30 ! initial entropy
    nbetter = 0
   ! call conf_entropy !first evaluation
    do 111 ini=1,inirand
     
        call shuffle
		call conf_entropy

        if(mod(ini,1000) .eq. 0) write(*,*)"Initial shuffle iterations: ",ini

        if(confentropy .le. Entmin)then
            nbetter = nbetter + 1
            Entmin = confentropy
			minatkeep = atype
			katomentropy = atomentropy
            call lmpdata("00Shuffle",nbetter)
            write(101,*)"Randomly Shuffle process"
            write(101,*)"The searched configuration No : ",nbetter
			write(101,*)"Shuffle Steps: ",ini			
            write(101,*) "***current normalized configurational entropy: ",confentropy/dble(natom)
            write(101,*)""
            write(*,*) "***current RANDOM searched better configuration No.: ",nbetter
            write(*,*)"Normalized configurational entropy: ",confentropy/dble(natom)
            call conf_entropy
            write(*,*)"call again Normalized configurational entropy: ",confentropy/dble(natom)

            write(*,*)""
            NN =0								
				do 1112 i=1,natom
					do 211 ine =1,5 !neighbor atom type
					ntemp = CN_No(i,ine)
						do 311 neID = 1,ntemp ! atom ID of a neighbor type
						JID = CN_ID(i,ine,neID)          
						if(atype(i) .eq. atype(JID))then
							NN(ine) = NN(ine) + 1
						endif 
				311 		continue                      
				211   continue  						
				1112 continue			
				
                 write(*,*) "***N1:",NN(1)/2," ***N2",NN(2)/2," ***N3",NN(3)/2," ***N4",NN(4)/2," ***N5",NN(5)/2	! the corresponding NN atom No.			
 
        endif
        111 continue
!pause
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
			filterValue = 0.
		if(Mod(imc,50).eq.0)then	
          write(*,*)""
          write(*,*)""		  
          write(*,*)"##********  MC iterations: ",imc
		  write(*,*)"Ave confentropy:", afterconfentropy/dble(natom)
          write(*,*)"*** CURRENT best:",Entmin/dble(natom)
		  write(*,*) "*** The first three nearest neighbour numbers for current best***"	
          write(*,*) "***N1:",NN(1)/2," ***N2",NN(2)/2," ***N3",NN(3)/2," ***N4",NN(4)/2," ***N5",NN(5)/2			  
          write(*,*) "KeepNo in SCF sub:",keepNo,"natom:",natom, "fraction: ",dble(keepNo)/dble(natom)			
		  write(*,*)""
		endif 
!**************** by OMP

		  call SCF ! use each poossile pair to get the best final atype

		afterconfentropy = confentropy
    
! output the lowest potential, atype, atomentropy 

            !!! Keep the better one
            if(confentropy .lt. Entmin)then
                nbetter = nbetter + 1
                Entmin = confentropy
				minatkeep = atype
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
                write(101,*) "KeepNo in SCF sub:",keepNo,"natom:",natom, "fraction: ",dble(keepNo)/dble(natom)			

                if(keepNo .eq. 0 )then !all bad atoms are done
                    write(*,*)"#MC process at MC step:", imc
					write(*,*)"***keepNo = 0, All atoms with larger scoring values have been processed!!***"
					write(*,*)"YOU CAN TERMINATE YOUR CODE or WAIT for BETTER by MC after several more MC runs!!!!!!!!!!!"
					write(*,*) "***N1:",NN(1)/2," ***N2",NN(2)/2," ***N3",NN(3)/2," ***N4",NN(4)/2," ***N5",NN(5)/2	! the corresponding NN atom No.			

					write(101,*)"***keepNo = 0, All atoms with larger scoring values have been processed!!***"
					write(101,*)"YOU CAN TERMINATE YOUR CODE after several more MC runs  !!!!!!!!!!!"
					
				endif 
	
				write(101,*)""
                write(*,*) "******current MC searched better configuration No.: ",nbetter
				write(*,*) "***current MC setp for this configuration:",imc				
                write(*,*)"Normalized configurational entropy: ",confentropy/dble(natom)
				write(*,*) "***N1:",NN(1)/2," ***N2",NN(2)/2," ***N3",NN(3)/2," ***N4",NN(4)/2," ***N5",NN(5)/2	! the corresponding NN atom No.	
				write(*,*) "KeepNo in SCF sub:",keepNo,"natom:",natom, "fraction: ",dble(keepNo)/dble(natom)			

                write(*,*)""

            endif
            ! begin Boltzmann evaluation

            Boltz =  Exp(-(confentropy-Entkeep)/(kT*dble(natom)))
            call random_number(r)
            if(confentropy .le. Entkeep)then
                Entkeep = confentropy
				atkeep = atype
				katomentropy = atomentropy ! keep atom potential
                naccept = naccept + 1
				
           elseif(r .lt. Boltz) then
               Entkeep = confentropy
				atkeep = atype
				katomentropy = atomentropy ! keep atom potential
               naccept = naccept + 1
            else
				confentropy = Entkeep
				atype = atkeep
				atomentropy = katomentropy ! keep atom potential
            endif

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
     !       print*,"kshuffle-1" 
            call kshuffle(imc)
      !      print*,"kshuffle-2" 

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
