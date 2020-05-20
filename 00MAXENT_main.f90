!This code do the Monte Carlo simulation to find the most
!uniform high-entropy alloy by the max-entropy theory
!(Entropy 2013, 15, 5536-5548;Computational Materials Science 142 (2018) 332–337))
!Developed by Prof. Shin-Pon Ju 2019/03/25 at UCLA
!major modification by Prof. Shin-Pon Ju 2020/03/12 at NSYSU
!
include 'information.f90'
PROGRAM MAXEnt
    use omp_lib
    USE information
    implicit real*8(a-h,o-z)
    integer nbetter ! better counter
    real*8 Boltz,r,ini_acceptratio
    integer NN(5)
    !omp
    integer tid,nthreads
    integer,allocatable:: omp_atype(:,:),omp_keepNo(:)
    real*8,allocatable:: omp_confentropy(:),omp_atomentropy(:,:)

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
    nthreads = 1 !how many threads do u want to use
    call omp_set_num_threads(nthreads)

    !open(969,file='omp.dat',status='unknown')
    allocate (omp_confentropy(0:nthreads-1))
    allocate (omp_atomentropy(0:nthreads-1,natom))
    allocate (omp_atype(0:nthreads-1,natom))
    allocate (omp_keepNo(0:nthreads-1))

!call conf_entropy !first evaluation
do 111 ini=1,inirand
!$OMP PARALLEL PRIVATE(tid,ith) shared(natom,weight,CN_No,CN_ID) FIRSTPRIVATE (atype)
tid = omp_get_thread_num()
!$OMP DO
    do 970 ith=1,nthreads       
        !if (.NOT. supercell)then
            !print *,"shuffle",ini
            call shuffle(natom,atype(:))
        !else
            !print *,"kshuffle",ini
            !call kshuffle(inirand)
        !endif
        call conf_entropy(tid,natom,atype(:),omp_confentropy(tid),omp_atomentropy(tid,:),weight(:),CN_No(:,:),CN_ID(:,:,:))
        omp_atype(tid,:)=atype(:) 
        !print*,"tid",tid,"confentropy(",tid,")",omp_confentropy(tid) 
        !print*,"tid",tid,"atomentropy(",tid,")",omp_atomentropy(tid,:)
        !write(*,*)"tid",tid,"omp_atype",atype
    970 continue
!$OMP END DO
!$OMP END PARALLEL

!找目前最小值
tid = minloc(omp_confentropy(:),DIM=1)-1 !tid min number is 0
confentropy = omp_confentropy(tid)
atomentropy = omp_atomentropy(tid,:)
atype = omp_atype(tid,:)

!print*,omp_confentropy(:) 
!print*,"tid",tid," confentropy: ",confentropy
!print*,"tid",tid,"atype",atype
!pause
!if(mod(ini,1000) .eq. 0) write(*,*)"Initial shuffle iterations: ",ini

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
        call conf_entropy(tid,natom,atype(:),confentropy,atomentropy(:),weight(:),CN_No(:,:),CN_ID(:,:,:))
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
    	        311 continue                      
    	    211 continue  						
    	1112 continue    	
         write(*,*) "***N1:",NN(1)/2," ***N2",NN(2)/2," ***N3",NN(3)/2," ***N4",NN(4)/2," ***N5",NN(5)/2	! the corresponding NN atom No. 
    endif
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

!將omp陣列=目前最小值
omp_confentropy(:) = confentropy
do 987 n=0,nthreads-1
    omp_atomentropy(n,:) = omp_atomentropy(tid,:)
    omp_atype(n,:) = omp_atype(tid,:)
987 continue

!print*,"no scf confentropy",confentropy
!print*,"confentropy",omp_confentropy(:)
!print*,"omp_data",omp_atype(0,:)
!print*,' '
!print*,"omp_data",omp_atype(1,:)
!print*,"atomentropy",atomentropy
!print*,"atype",atype
!pause

!FIXME: confentropy
do 112 imc=1,iterMC
    !$OMP PARALLEL PRIVATE(tid,ith,filterValue) shared(natom,weight,CN_No,CN_ID)
    tid = omp_get_thread_num()
    !$OMP DO
    do 971 ith=1,nthreads
        !acceptratio gradually becomes smaller with the increasing MC step, like annealing
        acceptratio = ini_acceptratio - ( dble(imc)/dble(iterMC) )*0.4 !! lowest: 0.5 (from 0.9)
        filterValue = 0.
        
        call SCF(tid,natom,omp_atype(tid,:),omp_confentropy(tid),omp_atomentropy(tid,:),filterValue, &
                 omp_keepNo(tid),weight(:),CN_No(:,:),CN_ID(:,:,:))! use each poossile pair to get the best final atype
        !call SCF(natom,atype(:),confentropy,atomentropy(:),filterValue,keepNo,weight,CN_No(:,:),CN_ID(:,:,:))
        !print*,"tid",tid,"confentropy(",tid,")",omp_confentropy(tid)
    971 continue
    !$OMP END DO
    !$OMP END PARALLEL 
!pause

!找目前最小值
tid = minloc(omp_confentropy(:),DIM=1)-1 !tid min number is 0
confentropy = omp_confentropy(tid)
atomentropy = omp_atomentropy(tid,:)
atype = omp_atype(tid,:)
keepNo = omp_keepNo(tid)
!print*,' '
!print*,"tid",tid,"min_confentropy",omp_confentropy(tid)
!print*,'keepNo',keepNo
!print*,' '

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
    
    afterconfentropy = confentropy
    ! output the lowest potential, atype, atomentropy 
    !!! Keep the better one

    if(confentropy .lt. Entmin)then
        nbetter = nbetter + 1
        Entmin = confentropy
		minatkeep = atype
		katomentropy = atomentropy ! keep atom potential
        !!!!!!!!!!!!!!! get neigbour atom number 
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
		        3 continue                      
		    2 continue  						
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

    !將omp陣列=目前最小值
    omp_confentropy(:) = confentropy
    do 986 n=0,nthreads-1
        omp_atomentropy(n,:) = omp_atomentropy(tid,:)
        omp_atype(n,:) = atype(:)
    986 continue

!print*,omp_confentropy(:)
!pause
    !! generate the next
    !  print*,"kshuffle-1" 
    !$OMP PARALLEL PRIVATE(tid,ith) shared(imc,natom,weight,CN_No,CN_ID,filterValue)
    tid = omp_get_thread_num()
    !$OMP DO
        do 973 ith=1,nthreads
        call kshuffle(tid,imc,natom,omp_atype(tid,:),omp_confentropy(tid),&
                      omp_atomentropy(tid,:),weight(:),CN_No(:,:),CN_ID(:,:,:))
        973 continue
    !$OMP END DO
    !$OMP END PARALLEL
    !print*,"confentropy",omp_confentropy(:)
    !print*,"omp",omp_atype(0,:)
    !print*,"omp",omp_atype(1,:)


    !  print*,"kshuffle-2" 
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
