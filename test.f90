program ocean

    use ncio 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: prec = sp 


    character(len=512) :: file_in 
    character(len=512) :: file_out 
    
    integer :: nx, ny 
    real(prec), allocatable :: xc(:)
    real(prec), allocatable :: yc(:)
    real(prec), allocatable :: f_grnd(:,:) 
    integer,    allocatable :: mask_ocn_ref(:,:) 
    integer,    allocatable :: mask_ocn(:,:) 

    ! Options 
    file_in  = "ANT-32KM_TOPO-RTOPO-2.0.1.nc"
    file_out = "./ocean.nc" 
    
    ! Determine size of input array
    nx = nc_size(file_in,"xc")
    ny = nc_size(file_in,"yc")

    ! Allocate working arrays
    allocate(xc(nx))
    allocate(yc(ny)) 
    allocate(f_grnd(nx,ny))
    allocate(mask_ocn(nx,ny))
    allocate(mask_ocn_ref(nx,ny))

    ! Load domain dimensions from input file 
    call nc_read(file_in,"xc",xc)
    call nc_read(file_in,"yc",yc)
    
    ! Load f_grnd from file 
    call load_f_grnd(f_grnd,file_in)

    ! Add a few holes in f_grnd that are not connected to the ocean 
    ! for testing the algorithm 
    f_grnd(130:135,80:82)   = 0.0 
    f_grnd(110:112,130:137) = 0.0 
    
    ! Define reference ocean mask: this sets the 
    ! seed location(s) of where we know open ocean exists 
    call define_ocn_ref(mask_ocn_ref,xc,yc,domain="Antarctica")
    
    ! Determine open-ocean mask for our domain 
    call find_ocean(mask_ocn,f_grnd,mask_ocn_ref)
    
    ! Write output to file 
    call write_ocn(file_out,f_grnd,mask_ocn_ref,mask_ocn,xc,yc)

contains 

subroutine find_ocean(mask,f_grnd,mask_ref)
    ! Brute-force routine to find all ocean points 
    ! (ie, when f_grnd < 1) connected to
    ! the open ocean as defined in mask_ref. 

    implicit none 

    integer,    intent(OUT) :: mask(:,:) 
    real(prec), intent(IN)  :: f_grnd(:,:)
    integer,    intent(IN)  :: mask_ref(:,:) 

    ! Local variables 
    integer :: i, j, q, nx, ny
    integer :: im1, ip1, jm1, jp1 
    integer :: n_unfilled 
    
    integer, parameter :: qmax = 1000 

    nx = size(mask,1)
    ny = size(mask,2) 

    ! First specify land points (mask==0)
    ! and assume all floating points are 'closed ocean' points (mask==-1)
    mask = 0 
    where (f_grnd .lt. 1.0) mask = -1 

    ! Now populate our ocean mask to be consistent 
    ! with known open ocean points 
    where (mask .eq. -1 .and. mask_ref .eq. 1) mask = 1 

if (.TRUE.) then 
    ! Iteratively fill in open-ocean points that are found
    ! next to known open-ocean points
    do q = 1, qmax 

        n_unfilled = 0 

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices 
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)
            
            if (mask(i,j) .eq. 1) then 
                ! This is an open-ocean point 
                ! Define any neighbor ocean points as open-ocean points
                
                if (mask(im1,j) .eq. -1) then 
                    mask(im1,j) = 1
                    n_unfilled = n_unfilled + 1
                end if 

                if (mask(ip1,j) .eq. -1) then 
                    mask(ip1,j) = 1 
                    n_unfilled = n_unfilled + 1
                end if

                if (mask(i,jm1) .eq. -1) then 
                    mask(i,jm1) = 1 
                    n_unfilled = n_unfilled + 1
                end if 

                if (mask(i,jp1) .eq. -1) then 
                    mask(i,jp1) = 1 
                    n_unfilled = n_unfilled + 1
                end if 

            end if
                
        end do 
        end do  

        write(*,*) q, n_unfilled, count(mask .eq. -1) 

        ! Exit loop if no more open-ocean points are found 
        if (n_unfilled .eq. 0) exit 

    end do 

end if 

    return 

end subroutine find_ocean

subroutine load_f_grnd(f_grnd,filename)

    implicit none 

    real(prec),         intent(OUT) :: f_grnd(:,:) 
    character(len=512), intent(IN)  :: filename 

    ! Local variables 
    integer :: nx, ny 
    real(prec), allocatable :: z_bed(:,:) 
    real(prec), allocatable :: H_ice(:,:)
    
    real(prec), parameter :: rho_ice =  910.0_prec 
    real(prec), parameter :: rho_sw  = 1028.0_prec 
    real(prec), parameter :: z_sl    = 0.0_prec 

    real(prec) :: rho_sw_ice 
    
    nx = size(f_grnd,1)
    ny = size(f_grnd,2) 

    allocate(z_bed(nx,ny))
    allocate(H_ice(nx,ny)) 

    ! Load input data 
    call nc_read(filename,"z_bed",z_bed)
    call nc_read(filename,"H_ice",H_ice)


    ! Calculate grounded fraction from ice thickness overburden 

    rho_sw_ice = rho_sw/rho_ice ! Ratio of density of seawater to ice [--]
        
    ! Calculate new H_grnd (ice thickness overburden)
    f_grnd = H_ice - rho_sw_ice*max(z_sl-z_bed,0.0_prec)

    where(f_grnd .ge. 0.0) f_grnd = 1.0 
    where(f_grnd .lt. 0.0) f_grnd = 0.0 

    return 

end subroutine load_f_grnd

subroutine define_ocn_ref(mask,xc,yc,domain)

    implicit none 

    integer,    intent(OUT) :: mask(:,:) 
    real(prec), intent(IN)  :: xc(:) 
    real(prec), intent(IN)  :: yc(:) 
    character(len=*), intent(IN) :: domain 

    ! Local variables 
    integer :: i, j 
    real(prec) :: x1, y1 

    ! First set mask to zero everywhere == not ocean 
    mask = 0 

    ! Define known ocean points depending on the domain 
    select case(trim(domain))

        case("Antarctica")

            ! Points above this value are open ocean 
            y1 = 2400.0
            j = minloc(abs(yc - y1),1)
            mask(:,j:ny) = 1 

            ! Points below this value are open ocean 
            y1 = -2300.0
            j = minloc(abs(yc - y1),1)
            mask(:,1:j) = 1 

            ! Points to the right of this value are open ocean 
            x1 = 2700.0 
            i = minloc(abs(xc - x1),1)
            mask(i:nx,:) = 1 

            ! Points to the left of this value are open ocean 
            x1 = -2700.0 
            i = minloc(abs(xc - x1),1)
            mask(1:i,:) = 1 

        case("Greenland")


        case DEFAULT 

            write(*,*) "define_ocn_ref:: Error: domain not recognized."
            write(*,*) "domain = ", trim(domain)
            stop 

    end select 


    return 

end subroutine define_ocn_ref

subroutine write_ocn(filename,f_grnd,mask_ref,mask,xc,yc)

    implicit none 

    character(len=512), intent(IN)  :: filename  
    real(prec),         intent(IN)  :: f_grnd(:,:) 
    integer,            intent(IN)  :: mask_ref(:,:)
    integer,            intent(IN)  :: mask(:,:) 
    real(prec),         intent(IN)  :: xc(:) 
    real(prec),         intent(IN)  :: yc(:) 
    
    call nc_create(filename) 

    call nc_write_dim(filename,"xc",x=xc)
    call nc_write_dim(filename,"yc",x=yc)
    
    call nc_write(filename,"f_grnd",  f_grnd,  dim1="xc",dim2="yc") 
    call nc_write(filename,"mask_ref",mask_ref,dim1="xc",dim2="yc") 
    call nc_write(filename,"mask",    mask,    dim1="xc",dim2="yc") 

    return 

end subroutine write_ocn

end program ocean

