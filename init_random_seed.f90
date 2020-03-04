SUBROUTINE init_random_seed()
  INTEGER :: i, n, clock
  INTEGER, ALLOCATABLE :: seed(:)

  call RANDOM_SEED(SIZE = n)
  ALLOCATE(seed(n))
  call SYSTEM_CLOCK(COUNT=clock)
  seed = clock - 2047 * (/ (i - 1, i = 1, n) /)
  seed = seed * 1103515245 + 2531011
  !PUT = seed
  !seed = 1000005
  call RANDOM_SEED(PUT = seed)
  DEALLOCATE(seed)
END SUBROUTINE