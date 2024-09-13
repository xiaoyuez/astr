!+---------------------------------------------------------------------+
!| This module contains subroutines for pre and post-process           |
!+---------------------------------------------------------------------+
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!|  28-05-2021  | Created by J. Fang @ Warrington                      |
!+---------------------------------------------------------------------+
module pp
  !
  use constdef
  use stlaio,  only: get_unit
  !
  implicit none
  !
  interface ProjectP2
    module procedure ProjectP2_2D
    module procedure ProjectP2_3D
  end interface
  !
  interface ProjectP3
    module procedure ProjectP3_2D
    module procedure ProjectP3_3D
  end interface
  !
  interface ProjectPi2
    module procedure ProjectPi2_2D
    module procedure ProjectPi2_3D
  end interface
  !
  interface ProjectPi3
    module procedure ProjectPi3_2D
    module procedure ProjectPi3_3D
  end interface
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is the entrance of post/pre-process.              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine ppentrance
    !
    use cmdefne
    use commvar,         only : ia,ja,ka
    use readwrite,       only : readinput,readic,readinput_serial
    use gridgeneration,  only : gridgen
    use solver,          only : refcal
    use udf_postprocess
    use parallel,        only : mpirank,bcast,mpisize,lio
    use commarray,       only : allocommarray
    !
    ! local data
    character(len=64) :: cmd,casefolder,inputfile,outputfile,viewmode, &
                         flowfieldfile, readmode
    integer :: filenumb
    integer :: ian,jan,kan
    !
    if(mpirank == 0) then
      call readkeyboad(cmd)
      print*,' ** pp command: ',cmd
    endif
    !
    call bcast(cmd)
    !
    if(trim(cmd)=='init') then
      !
      call readkeyboad(casefolder)
      call examplegen(trim(casefolder))
      ! generate an example channel flow case
      !
    elseif(trim(cmd)=='gridgen') then
      call readinput
      !
      call allocommarray
      !
      call gridgen
    elseif(trim(cmd)=='hitgen') then
      call readinput
      !
      if(mpirank == 0) then
        call readkeyboad(readmode) 
      endif
      call bcast(readmode)
      !
      if(trim(readmode)=='serial') then
        !
        if(mpisize .ne. 1) then
          stop ' Serial mode please only use 1 proc'
        endif
        !
        call hitgen
      elseif(trim(readmode)=='parallel') then
        !
        if(mpirank == 0)   print *, ' Parallel used'
        call hitgen_parallel
        !
      else
        print* ,"Readmode is not defined!", readmode
      endif
      !
    elseif(trim(cmd)=='hitgen2d') then
      !
      if(mpirank == 0) then
            call readkeyboad(readmode) 
        endif
      call bcast(readmode)
      !
      if(trim(readmode)=='serial') then
        !
        if(mpisize .ne. 1) then
          stop ' Serial mode please only use 1 proc'
        endif
        !
        call hitgen2d
      elseif(trim(readmode)=='parallel') then
        !
        if(mpirank == 0)   print *, ' Parallel used'
        call hitgen2d_parallel
        !
      else
        print* ,"Readmode is not defined!", readmode
      endif
    elseif(trim(cmd)=='solid') then
      call solidpp
    elseif(trim(cmd)=='datacon') then
      !
      call readkeyboad(inputfile)
      !
      if(trim(inputfile)=='all') then
        call stream2struc('outdat/flowfield')
        call stream2struc('outdat/meanflow')
        call stream2struc('outdat/2ndsta')
        call stream2struc('outdat/3rdsta')
        call stream2struc('outdat/budget')
      else
        call stream2struc(trim(inputfile))
      endif
      !
    elseif(trim(cmd)=='parinfo') then
      !
      call parallelifogen
      !
    elseif(trim(cmd)=='scale') then
      !
      if(mpirank == 0) then
        call readkeyboad(flowfieldfile)
        print*,' ** flowfieldfile command: ',flowfieldfile
      endif
      call bcast(flowfieldfile)
      !
      !
      if(mpirank == 0) then
        call readkeyboad(inputfile) 
        read(inputfile,'(i4)') ian
        call readkeyboad(inputfile) 
        read(inputfile,'(i4)') jan
        call readkeyboad(inputfile) 
        read(inputfile,'(i4)') kan
        print*,' ** Output size: ',ian,'x',jan,'x',kan
      endif
      call bcast(ian)
      call bcast(jan)
      call bcast(kan)
      !
      call fieldscale(trim(flowfieldfile),ian,jan,kan) 
      !
    elseif(trim(cmd)=='view') then
      !
      call readkeyboad(flowfieldfile)
      call readkeyboad(outputfile)
      call readkeyboad(viewmode)
      call readkeyboad(inputfile)
      !
      call fieldview(trim(flowfieldfile),trim(outputfile),trim(viewmode),trim(inputfile))
      !
    elseif(trim(cmd)=='flamegen') then
        !
      call readkeyboad(flowfieldfile)
      call readkeyboad(viewmode)
      call readkeyboad(inputfile)
        !
      call flamegen(trim(flowfieldfile),trim(viewmode),trim(inputfile))
        !
    elseif(trim(cmd)=='udf') then
        !
        call readkeyboad(inputfile)
        !
        call flame2d_pp(trim(inputfile))
    elseif(trim(cmd)=='spectra') then
      !
      !
      if(mpirank == 0) then
          call readkeyboad(readmode) 
      endif
      call bcast(readmode)
      !
      if(trim(readmode)=='instant2D') then
        !
        if(mpirank == 0) then
          call readkeyboad(inputfile) 
          read(inputfile,'(i4)') filenumb
        endif
        call bcast(filenumb)
        call instantspectra2D(filenumb)
        !
      elseif(trim(readmode)=='instant3D') then
        !
        if(mpirank == 0) then
          call readkeyboad(inputfile) 
          read(inputfile,'(i4)') filenumb
        endif
        call bcast(filenumb)
        call instantspectra3D(filenumb)
        !
      elseif(trim(readmode)=='instant2Davg') then
        !
        if(mpirank == 0) then
          call readkeyboad(inputfile) 
          read(inputfile,'(i4)') filenumb
        endif
        call bcast(filenumb)
        call instantspectra2Davg(filenumb)
        !
      elseif(trim(readmode)=='instant3Davg') then
        !
        if(mpirank == 0) then
          call readkeyboad(inputfile) 
          read(inputfile,'(i4)') filenumb
        endif
        call bcast(filenumb)
        call instantspectra3Davg(filenumb)
        !
      elseif(trim(readmode)=='skewness') then
        !
        if(mpirank == 0) then
          call readkeyboad(inputfile) 
          read(inputfile,'(i4)') filenumb
        endif
        call bcast(filenumb)
        call instantspectraskewness(filenumb)
        !
      elseif(trim(readmode)=='triad2D') then
        !
        if(mpirank == 0) then
          call readkeyboad(inputfile) 
          read(inputfile,'(i4)') filenumb
        endif
        call bcast(filenumb)
        call instanttriad2D(filenumb)
        !
      elseif(trim(readmode)=='triad3D') then
        !
        if(mpirank == 0) then
          call readkeyboad(inputfile) 
          read(inputfile,'(i4)') filenumb
        endif
        call bcast(filenumb)
        call instanttriad3D(filenumb)
        !
      elseif(trim(readmode)=='triad2Davg') then
        !
        if(mpirank == 0) then
          call readkeyboad(inputfile) 
          read(inputfile,'(i4)') filenumb
        endif
        call bcast(filenumb)
        call instanttriad2Davg(filenumb)
        !
      elseif(trim(readmode)=='triad3Davg') then
        !
        if(mpirank == 0) then
          call readkeyboad(inputfile) 
          read(inputfile,'(i4)') filenumb
        endif
        call bcast(filenumb)
        call instanttriad3Davg(filenumb)
        !
      else
        print* ,"Readmode is not defined!", readmode
      endif
      !
    elseif(trim(cmd)=='velgradient') then
      !
      if(mpirank == 0) then
        call readkeyboad(readmode)
        print*,' ** readmode command: ',readmode
      endif
      call bcast(readmode)
      !
      if(trim(readmode)=='instant') then
        ! 
        if(mpirank == 0) then
          call readkeyboad(inputfile) 
          read(inputfile,'(i4)') filenumb
          print*,' ** Filenumb: ',filenumb
        endif
        call bcast(filenumb)
        call instantvelgradient(filenumb)
        !
      elseif(trim(readmode)=='PQR') then
        ! TODO: write PQR
        if(mpirank == 0) then
          call readkeyboad(inputfile) 
          read(inputfile,'(i4)') filenumb
          print*,' ** Filenumb: ',filenumb
        endif
        call bcast(filenumb)
        stop ' !! PQR not defined'
        !call PQR
      elseif(trim(readmode)=='ScaleLen') then
        !
        if(mpirank == 0) then
          call readkeyboad(inputfile) 
          read(inputfile,'(i4)') filenumb
          print*,' ** Filenumb: ',filenumb
        endif
        call bcast(filenumb)
        !
        call velgradient_scale_lengths(filenumb) 
        !
      else
        print* ,"Readmode is not defined!", readmode
      endif
      !
    elseif(trim(cmd)=='SGS') then
      !
      !
      if(mpirank == 0) then
        call readkeyboad(readmode)
        print*,' ** readmode command: ',readmode
      endif
      call bcast(readmode)
      !
      if(trim(readmode)=='Pi2Dint') then
        ! 
        if(mpirank == 0) then
          print* ,"Use SGSPi2Dint"
          call readkeyboad(inputfile) 
          read(inputfile,'(i4)') filenumb
          print*,' ** Filenumb: ',filenumb
        endif
        call bcast(filenumb)
        call SGSPi2Dint(filenumb)
        !
      elseif(trim(readmode)=='Pi2Dtot') then
          ! 
          if(mpirank == 0) then
            print* ,"Use SGSPi2Dtot"
            call readkeyboad(inputfile) 
            read(inputfile,'(i4)') filenumb
            print*,' ** Filenumb: ',filenumb
          endif
          call bcast(filenumb)
          call SGSPi2Dtot(filenumb)
          !
      elseif(trim(readmode)=='PiOmgea2D') then
            ! 
            if(mpirank == 0) then
              print* ,"Use SGSPiOmgea2D"
              call readkeyboad(inputfile) 
              read(inputfile,'(i4)') filenumb
              print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSPiOmgea2D(filenumb)
            !
      elseif(trim(readmode)=='Pi2Dlocal') then
          ! 
          if(mpirank == 0) then
            print* ,"Use SGSPi2Dlocal"
            call readkeyboad(inputfile) 
            read(inputfile,'(i4)') filenumb
            print*,' ** Filenumb: ',filenumb
          endif
          call bcast(filenumb)
          call SGSPi2Dlocal(filenumb)
          !
      elseif(trim(readmode)=='Pi3Dint') then
        ! 
        if(mpirank == 0) then
          print* ,"Use SGSPi3Dint"
          call readkeyboad(inputfile) 
          read(inputfile,'(i4)') filenumb
          print*,' ** Filenumb: ',filenumb
        endif
        call bcast(filenumb)
        call SGSPi3Dint(filenumb)
        !
      elseif(trim(readmode)=='Pi3Dtot') then
          ! 
          if(mpirank == 0) then
            print* ,"Use SGSPi3Dtot"
            call readkeyboad(inputfile) 
            read(inputfile,'(i4)') filenumb
            print*,' ** Filenumb: ',filenumb
          endif
          call bcast(filenumb)
          call SGSPi3Dtot(filenumb)
          !
      elseif(trim(readmode)=='Pi3Dlocal') then
          ! 
          if(mpirank == 0) then
            print* ,"Use SGSPi3Dlocal"
            call readkeyboad(inputfile) 
            read(inputfile,'(i4)') filenumb
            print*,' ** Filenumb: ',filenumb
          endif
          call bcast(filenumb)
          call SGSPi3Dlocal(filenumb)
          !
      elseif(trim(readmode)=='LES3D') then
          ! 
          if(mpirank == 0) then
            print* ,"Use SGSLES3D"
            call readkeyboad(inputfile) 
            read(inputfile,'(i4)') filenumb
            print*,' ** Filenumb: ',filenumb
          endif
          call bcast(filenumb)
          call SGSLES3D(filenumb)
          !
      elseif(trim(readmode)=='LES2D') then
          ! 
          if(mpirank == 0) then
            print* ,"Use SGSLES2D"
            call readkeyboad(inputfile) 
            read(inputfile,'(i4)') filenumb
            print*,' ** Filenumb: ',filenumb
          endif
          call bcast(filenumb)
          ! call SGSLES2D(filenumb)
          !
      elseif(trim(readmode)=='T3D') then
        ! 
        if(mpirank == 0) then
          print* ,"Use SGST3D"
          call readkeyboad(inputfile) 
          read(inputfile,'(i4)') filenumb
          print*,' ** Filenumb: ',filenumb
        endif
        call bcast(filenumb)
        call SGST3D(filenumb)
        !
      elseif(trim(readmode)=='stress2D') then
          ! 
          if(mpirank == 0) then
            print* ,"Use SGSstress2D"
            call readkeyboad(inputfile) 
            read(inputfile,'(i4)') filenumb
            print*,' ** Filenumb: ',filenumb
          endif
          call bcast(filenumb)
          call SGSstress2D(filenumb)
          !
      elseif(trim(readmode)=='stress3D') then
        ! 
        if(mpirank == 0) then
          print* ,"Use SGSstress3D"
          call readkeyboad(inputfile) 
          read(inputfile,'(i4)') filenumb
          print*,' ** Filenumb: ',filenumb
        endif
        call bcast(filenumb)
        call SGSstress3D(filenumb)
        !
      else
        print* ,"Readmode is not defined!", readmode
      endif
      !
      elseif(trim(cmd)=='initparam') then
        !
        !
        if(mpirank == 0) then
          call readkeyboad(readmode)
          print*,' ** readmode command: ',readmode
        endif
        call bcast(readmode)
        !
        if(trim(readmode)=='2D') then
          !
          call initparam2D
          !
        elseif(trim(readmode)=='3D') then
          ! 
          call initparam3D
          !
        else
          print* ,"Readmode is not defined!", readmode
        endif
        !
    else
      stop ' !! pp command not defined. @ ppentrance'
    endif
    ! 
    !
  end subroutine ppentrance
  !+-------------------------------------------------------------------+
  !| The end of the subroutine preprocess.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate parallel.info file.                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-04-2022  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine parallelifogen
    !
    use readwrite, only: readinput
    use parallel,  only: mpisizedis,parapp,mpisize,mpirankmax
    !
    call readinput
    !
    mpisize=128*128
    mpirankmax=mpisize-1
    !
    call mpisizedis
    !
    call parapp
    !
  end subroutine parallelifogen
  !+-------------------------------------------------------------------+
  !| The end of the subroutine parallelifogen.                         |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate an example case.              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-05-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine examplegen(folder)
    !
    use parallel, only: mpirank
    !
    character(len=*),intent(in) :: folder
    integer :: fh
    !
    if(mpirank==0) then
      !
      print*,' ** Generating an example case.'
      !
      call system('mkdir '//trim(folder))
      call system('mkdir '//trim(folder)//'/datin')
      !
      fh=get_unit()
      !
      open(fh,file=trim(folder)//'/datin/input.chl',form='formatted')
      write(fh,'(A)')'########################################################################'
      write(fh,'(A)')'#                     input file of ASTR code                          #'
      write(fh,'(A)')'########################################################################'
      write(fh,*)
      write(fh,'(A)')'# flowtype                                              : The type of flow problem'
      write(fh,'(A)')'channel'
      write(fh,*)
      write(fh,'(A)')'# im,jm,km                                              : The size of grid.'
      write(fh,'(A)')'128,128,128'
      write(fh,*)
      write(fh,'(A)')'# lihomo,ljhomo,lkhomo                                  : The homogeneous directions'
      write(fh,'(A)')'t,f,t'
      write(fh,*)
      write(fh,'(A)')'# nondimen,diffterm,lfilter,lreadgrid,lfftz,limmbou     : Parameters'
      write(fh,'(A)')'t,t,t,t,f,f'
      write(fh,*)
      write(fh,'(A)')'# lrestar                                               : start mode'
      write(fh,'(A)')'f'
      write(fh,*)
      write(fh,'(A)')'# alfa_filter, kcutoff                                  : Filter parameters'
      write(fh,'(A)')'0.49d0, 48'
      write(fh,*)
      write(fh,'(A)')'# ref_tem,reynolds,mach                                   : Reference variables'
      write(fh,'(A)')'273.15d0,  3000.d0,  0.5d0 '
      write(fh,*)
      write(fh,'(A)')'# conschm,difschm,rkscheme                              : Numerical scheme'
      write(fh,'(A)')'642c, 642c, rk4 '
      write(fh,*)
      write(fh,'(A)')'# recon_schem, lchardecomp,bfacmpld                     : Parameters for upwind-biased scheme'
      write(fh,'(A)')'0, f, 0.3d0 '
      write(fh,*)
      write(fh,'(A)')'# num_species                                           : number of species'
      write(fh,'(A)')'1'
      write(fh,*)
      write(fh,'(A)')'# turbmode                                              : turbulence model'
      write(fh,'(A)')'none'
      write(fh,*)
      write(fh,'(A)')'# bctype                                                : Boundary condition definition '
      write(fh,'(A)')'1'
      write(fh,'(A)')'1'
      write(fh,'(A)')'41, 1.d0'
      write(fh,'(A)')'41, 1.d0'
      write(fh,'(A)')'1'
      write(fh,'(A)')'1'
      write(fh,*)
      write(fh,'(A)')'# ninit                                                 : Initial method'
      write(fh,'(A)')'1'
      write(fh,*)
      write(fh,'(A)')'# spg_imin,spg_imax,spg_jmin,spg_jmax,spg_kmin,spg_kmax : Sponge layer range'
      write(fh,'(A)')'0, 0, 0, 0, 0, 0'
      write(fh,*)
      write(fh,'(A)')'# gridfile                                              : grid'
      write(fh,'(A)')'./datin/grid.chl'
      write(fh,*)
      write(fh,*)
      write(fh,*)
      write(fh,*)'########################################################################'
      write(fh,*)'# bctype                                                               #'
      write(fh,*)'#   1 : periodic bc,     nothing will be done.                         #'
      write(fh,*)'#  41 : isothermal wall, wall temperature input.                       #'
      write(fh,*)'########################################################################'
      close(fh)
      print*,' << ',trim(folder),'/datin/input.chl'
      !
      call gridchannel(320,128,128,trim(folder))
      !
      open(fh,file=trim(folder)//'/datin/controller',form='formatted')
      write(fh,'(A)')'############################################################'
      write(fh,'(A)')'#               control file for ASTRR code                #'
      write(fh,'(A)')'############################################################'
      write(fh,*)
      write(fh,'(A)')'# lwsequ,lwslic,lavg'
      write(fh,'(A)')'f,f'
      write(fh,*)
      write(fh,'(A)')'# maxstep,nwrite,ninst,nlstep,navg'
      write(fh,'(A)')'1000,100,20,1,20'
      write(fh,*)
      write(fh,'(A)')'# deltat '
      write(fh,'(A)')'1.d-3'
      write(fh,'(A)')'+----------------------------------------------------------+'
      write(fh,'(A)')'| This file will be read each time after checkpoint        |'
      write(fh,'(A)')'+--------------+-------------------------------------------+'
      write(fh,'(A)')'|       lwrite | to write a sequence of flowfield files.   |'
      write(fh,'(A)')'|       lwslic | to write a sequence of slice cut files.   |'
      write(fh,'(A)')'|         lavg | to conduct on-fly stastistics.            |'
      write(fh,'(A)')'+--------------+-------------------------------------------+'
      write(fh,'(A)')'|      maxstep | max step to run.                          |'
      write(fh,'(A)')'|       nwrite | frequency of dumping checkpoint.          |'
      write(fh,'(A)')'|        ninst | frequency of writing slice.               |'
      write(fh,'(A)')'|       nlstep | frequency of listing computing state.     |'
      write(fh,'(A)')'|         navg | frequency of calculating statistics.      |'
      write(fh,'(A)')'+--------------+-------------------------------------------+'
      write(fh,'(A)')'|       deltat | time step.                                |'
      write(fh,'(A)')'+--------------+-------------------------------------------+'
      close(fh)
      print*,' << ',trim(folder),'/datin/controller'
      !
      call system('cp -v ./bin/astr ./'//trim(folder))
      !
      print*,' ** An example case is generated.'
      print*,' ** you can now run a simulation of channel flow by : '
      print*,' mpirun -np 8 ./astr run datin/input.chl'
      !
    endif
    !
  end subroutine examplegen
  !+-------------------------------------------------------------------+
  !| The end of the subroutine examplegen.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| The end of the subroutine examplegen.                             |
  !+-------------------------------------------------------------------+
  !
  subroutine fieldscale(flowfile,ian,jan,kan)
    !
    use hdf5io
    use readwrite, only : readinput
    use commvar,   only : im,jm,km,ia,ja,ka
    use commarray, only : rho,vel,prs,tmp
    use parallel,  only : dataswap, mpisizedis,parapp,parallelini,mpirank, mpistop,mpi_ikgroup,mpi_kgroup
    use solver,    only : refcal
    !
    character(len=*),intent(in) :: flowfile
    integer,intent(in) :: ian,jan,kan
    !
    real(8), allocatable, dimension(:,:,:) :: rhon,prsn,tmpn
    real(8), allocatable, dimension(:,:,:,:) :: veln
    !
    integer :: i,j,k,ratio1,ratio2,ratio3,l,m,n,imn,jmn,kmn
    !
    character(len=1) :: modeio
    character(len=128) :: outfilename
    !
    call readinput
    !
    call mpisizedis
    if(mpirank==0)then
      print*, '** mpisizedis done!'
    endif
    !
    call parapp
    if(mpirank==0)then
      print*, '** parapp done!'
    endif
    !
    call parallelini
    if(mpirank==0)then
      print*, '** parallelini done!'
    endif
    !
    call refcal
    if(mpirank==0)then
      print*, '** refcal done!'
    endif
    !
    modeio='h'
    !
    if(mpirank==0)then
      if(ka==0)then
        print *,"2D, ia:",ia,",ja:",ja
      else
        print *,"3D, ia:",ia,",ja:",ja, ",ka:", ka
      endif
    endif
    !
    !
    !allocate(x(0:im,0:jm,0:km,1:3) )
    allocate(vel(0:im,0:jm,0:km,1:3))
    allocate(rho(0:im,0:jm,0:km),prs(0:im,0:jm,0:km),tmp(0:im,0:jm,0:km))
    !
    !call geomcal
    !
    call h5io_init(filename=flowfile,mode='read')
    !
    call h5read(varname='ro',var=rho(0:im,0:jm,0:km), mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='p',var=prs(0:im,0:jm,0:km), mode = modeio)
    call h5read(varname='t',var=tmp(0:im,0:jm,0:km), mode = modeio)
    !
    call h5io_end
    !
    !!!! Global size verification assumption
    if(mpirank==0)then
      if(mod(ia,ian) .ne. 0)then
        print *,'ia=',ia,'ian=',ian
        stop "Ia should be multiple of new ia"
      endif
      !
      if(mod(ja,jan) .ne. 0)then
        print *,'ja=',ja,'jan=',jan
        stop "Ja should be multiple of new ja"
      endif
      !
      if(ka .ne. 0)then
        if(mod(ka,kan) .ne. 0)then
          print *,'ka=',ka,'kan=',kan
          stop "Ka should be multiple of new ka"
        endif
      endif
      print *, "* Global size verification passed!"
    endif
    !
    !!!! Local size verification assumption
    ratio1 = ia/ian
    ratio2 = ja/jan
    if(ka .ne. 0)then
      ratio3 = ka/kan
    endif
    !
    if(mod(im,ratio1).ne.0)then
      print *,'mpirank=',mpirank,'im=',im,"ratio1=",ratio1
      stop "Not possible for integer division of im by ratio1"
    endif
    !
    if(mod(jm,ratio2).ne.0)then
      print *,'mpirank=',mpirank,'jm=',jm,"ratio2=",ratio2
      stop "Not possible for integer division of jm by ratio2"
    endif
    !
    if(ka .ne. 0)then
      if(mod(km,ratio3).ne.0)then
        print *,'mpirank=',mpirank,'km=',jm,"ratio3=",ratio3
        stop "Not possible for integer division of km by ratio3"
      endif
    endif
    print *, "* Local size verification passed for mpirank=",mpirank
    ! 
    !!!! Do scale
    if(ka == 0)then
      !! 2D case
      imn = im / ratio1
      jmn = jm / ratio2
      !
      allocate(veln(0:imn,0:jmn,0:0,1:3))
      allocate(rhon(0:imn,0:jmn,0:0), &
               prsn(0:imn,0:jmn,0:0), &
               tmpn(0:imn,0:jmn,0:0))
      !
      rhon = 0.d0
      veln = 0.d0
      prsn = 0.d0
      tmpn = 0.d0
      !
      do i=1,imn
      do j=1,jmn
        do m = 0,(ratio1-1)
        do n = 0,(ratio2-1)
          rhon(i,j,0) = rhon(i,j,0) + rho(ratio1*i - m,ratio2*j - n,0)
          veln(i,j,0,1) = veln(i,j,0,1) + vel(ratio1*i - m,ratio2*j - n,0,1)
          veln(i,j,0,2) = veln(i,j,0,2) + vel(ratio1*i - m,ratio2*j - n,0,2)
          prsn(i,j,0) = prsn(i,j,0) + prs(ratio1*i - m,ratio2*j - n,0)
          tmpn(i,j,0) = tmpn(i,j,0) + tmp(ratio1*i - m,ratio2*j - n,0)
        end do
        end do
      end do
      end do
      !
      rhon = rhon / (ratio1*ratio2)
      veln = veln / (ratio1*ratio2)
      prsn = prsn / (ratio1*ratio2)
      tmpn = tmpn / (ratio1*ratio2)
      !
      rhon(0,:,:)  =rhon(imn,:,:)
      veln(0,:,:,:)=veln(imn,:,:,:)
      prsn(0,:,:)  =prsn(imn,:,:)
      tmpn(0,:,:)  =tmpn(imn,:,:)
      !
      rhon(:,0,:)  =rhon(:,jmn,:)
      veln(:,0,:,:)=veln(:,jmn,:,:)
      prsn(:,0,:)  =prsn(:,jmn,:)
      tmpn(:,0,:)  =tmpn(:,jmn,:)
      !
      !
      ia = ia / ratio1
      ja = ja / ratio2
      !
      call parapp
      !
      deallocate(mpi_ikgroup)
      !
      call parallelini
      !
      if((imn .ne. im) .and. (jmn .ne. jm))then
        stop "Problem emerge when redistributing the array"
      end if
      !
      outfilename = flowfile//'_scaled.'//modeio//'5'
      !
      call h5io_init(trim(outfilename),mode='write')
      call h5write(varname='ro',var=rhon(0:im,0:jm,0),  dir='k')
      call h5write(varname='u1',var=veln(0:im,0:jm,0,1),  dir='k')
      call h5write(varname='u2',var=veln(0:im,0:jm,0,2),  dir='k')
      call h5write(varname='p', var=prsn(0:im,0:jm,0),  dir='k')
      call h5write(varname='t', var=tmpn(0:im,0:jm,0),  dir='k')
      call h5io_end
      !
      if(mpirank==0)then
        print*,' <<< ', outfilename, '... done.'
      endif
      !
    else
      !! 3D case
      allocate(veln(0:(im/ratio1),0:(jm/ratio2),0:(km/ratio3),1:3))
      allocate(rhon(0:(im/ratio1),0:(jm/ratio2),0:(km/ratio3)), &
               prsn(0:(im/ratio1),0:(jm/ratio2),0:(km/ratio3)), &
               tmpn(0:(im/ratio1),0:(jm/ratio2),0:(km/ratio3)))
      !
      rhon = 0.d0
      veln = 0.d0
      prsn = 0.d0
      tmpn = 0.d0
      !
      do i=1,(im/ratio1)
      do j=1,(jm/ratio2)
      do k=1,(km/ratio3)
        do m = 0,(ratio1-1)
        do n = 0,(ratio2-1)
        do l = 0,(ratio3-1)
          rhon(i,j,k)   = rhon(i,j,k)   + rho(ratio1*i - m,ratio2*j - n,ratio3*k-l)
          veln(i,j,k,1) = veln(i,j,k,1) + vel(ratio1*i - m,ratio2*j - n,ratio3*k-l,1)
          veln(i,j,k,2) = veln(i,j,k,2) + vel(ratio1*i - m,ratio2*j - n,ratio3*k-l,2)
          veln(i,j,k,3) = veln(i,j,k,3) + vel(ratio1*i - m,ratio2*j - n,ratio3*k-l,3)
          prsn(i,j,k)   = prsn(i,j,k)   + prs(ratio1*i - m,ratio2*j - n,ratio3*k-l)
          tmpn(i,j,k)   = tmpn(i,j,k)   + tmp(ratio1*i - m,ratio2*j - n,ratio3*k-l)
        end do
        end do
        end do
      end do
      end do
      end do
      !
      rhon = rhon / (ratio1*ratio2*ratio3)
      veln = veln / (ratio1*ratio2*ratio3)
      prsn = prsn / (ratio1*ratio2*ratio3)
      tmpn = tmpn / (ratio1*ratio2*ratio3)
      !
      rhon(0,:,:)  =rhon((im/ratio1),:,:)
      veln(0,:,:,:)=veln((im/ratio1),:,:,:)
      prsn(0,:,:)  =prsn((im/ratio1),:,:)
      tmpn(0,:,:)  =tmpn((im/ratio1),:,:)
      !
      rhon(:,0,:)  =rhon(:,(jm/ratio2),:)
      veln(:,0,:,:)=veln(:,(jm/ratio2),:,:)
      prsn(:,0,:)  =prsn(:,(jm/ratio2),:)
      tmpn(:,0,:)  =tmpn(:,(jm/ratio2),:)
      !
      rhon(:,:,0)  =rhon(:,:,(km/ratio3))
      veln(:,:,0,:)=veln(:,:,(km/ratio3),:)
      prsn(:,:,0)  =prsn(:,:,(km/ratio3))
      tmpn(:,:,0)  =tmpn(:,:,(km/ratio3))
      !
      ia = ia / ratio1
      ja = ja / ratio2
      ka = ka / ratio3
      !
      call parapp
      !
      deallocate(mpi_ikgroup)
      !
      call parallelini
      !
      if((imn .ne. im) .and. (jmn .ne. jm) .and. (kmn .ne. km))then
        stop "Problem emerge when redistributing the array"
      end if
      !
      outfilename = flowfile//'_scaled.'//modeio//'5'
      !
      call h5io_init(trim(outfilename),mode='write')
      call h5write(varname='ro',var=rhon(0:im,0:jm,0:km), mode = modeio)
      call h5write(varname='u1',var=veln(0:im,0:jm,0:km,1), mode = modeio)
      call h5write(varname='u2',var=veln(0:im,0:jm,0:km,2), mode = modeio)
      call h5write(varname='u3',var=veln(0:im,0:jm,0:km,3), mode = modeio)
      call h5write(varname='p', var=prsn(0:im,0:jm,0:km), mode = modeio)
      call h5write(varname='t', var=tmpn(0:im,0:jm,0:km), mode = modeio)
      call h5io_end
      !
      if(mpirank==0)then
        print*,' <<< ', outfilename, '... done.'
      endif
      !
    endif
    !
    call mpistop
    !
    deallocate(rho,vel,prs,tmp)
    deallocate(veln,rhon,prsn,tmpn)
    !
  end subroutine fieldscale
  !+-------------------------------------------------------------------+
  !| This subroutine is used to convert a streammed data to a          |
  !| structured data.                                                  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 31-03-2022  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine stream2struc(fname)
    !
    use hdf5io
    !
    character(len=*),intent(in) :: fname
    !
    integer :: fh,isize,jsize,ksize,rankmax,jrank,vai1,irp,jrp,krp,    &
               asizemax,ima,jma,kma,i,j,k,n,nstep
    real(8) :: time
    integer,allocatable :: im(:),jm(:),km(:),i0(:),j0(:),k0(:),        &
                           offset(:),asize(:)
    !
    real(8),allocatable,dimension(:) :: data_1d
    real(8),allocatable,dimension(:,:,:) :: data_3d
    !
    fh=get_unit()
    open(fh,file='datin/parallel.info',form='formatted',action='read')
    read(fh,*)
    read(fh,*)isize,jsize,ksize
    rankmax=isize*jsize*ksize-1
    print*,' ** rankmax=',rankmax
    allocate(im(0:rankmax),jm(0:rankmax),km(0:rankmax),                &
             i0(0:rankmax),j0(0:rankmax),k0(0:rankmax),                &
             offset(0:rankmax),asize(0:rankmax))
    read(fh,*)
    do jrank=0,rankmax
      read(fh,*)vai1,irp,jrp,krp,im(jrank),jm(jrank),km(jrank),        &
                i0(jrank),j0(jrank),k0(jrank)
    enddo
    close(fh)
    print*,' >> datin/parallel.info'
    !
    do jrank=0,rankmax
      asize(jrank)=(im(jrank)+1)*(jm(jrank)+1)*(km(jrank)+1)
    enddo
    !
    offset(0)=0
    do jrank=1,rankmax
      offset(jrank)=offset(jrank-1)+asize(jrank-1)
    enddo
    !
    asizemax=offset(rankmax)+asize(rankmax)
    !
    ima=i0(rankmax)+im(rankmax)
    jma=j0(rankmax)+jm(rankmax)
    kma=k0(rankmax)+km(rankmax)
    !
    print*,' **     stream data size:',asizemax
    print*,' ** structured data size:',ima,jma,kma
    ! !
    allocate( data_1d(1:asizemax),data_3d(0:ima,0:jma,0:kma) )
    !
    call h5sread(varname='nstep',var=nstep,filename=fname//'.s5')
    call h5srite(var=nstep,varname='nstep',filename=fname//'.h5',newfile=.true.)
    !
    if(fname=='outdat/flowfield') then
      !
      call h5sread(varname='time',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='time',filename=fname//'.h5')
      !
      call data3d_rcw(varname='ro')
      call data3d_rcw(varname='u1')
      call data3d_rcw(varname='u2')
      call data3d_rcw(varname='u3')
      call data3d_rcw(varname='p')
      call data3d_rcw(varname='t')
      !
    elseif(fname=='outdat/meanflow') then
      !
      call h5sread(varname='nsamples',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='nsamples',filename=fname//'.h5')
      !
      call h5sread(varname='nstep_sbeg',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='nstep_sbeg',filename=fname//'.h5')
      !
      call h5sread(varname='time_sbeg',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='time_sbeg',filename=fname//'.h5')
      !
      call data3d_rcw(varname='rom')
      call data3d_rcw(varname='u1m')
      call data3d_rcw(varname='u2m')
      call data3d_rcw(varname='u3m')
      call data3d_rcw(varname='pm')
      call data3d_rcw(varname='tm')
      !
    elseif(fname=='outdat/2ndsta') then
      !
      call h5sread(varname='nsamples',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='nsamples',filename=fname//'.h5')
      !
      call h5sread(varname='nstep_sbeg',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='nstep_sbeg',filename=fname//'.h5')
      !
      call h5sread(varname='time_sbeg',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='time_sbeg',filename=fname//'.h5')
      !
      call data3d_rcw(varname= 'pp')
      call data3d_rcw(varname= 'tt')
      call data3d_rcw(varname='tu1')
      call data3d_rcw(varname='tu2')
      call data3d_rcw(varname='tu3')
      call data3d_rcw(varname='u11')
      call data3d_rcw(varname='u12')
      call data3d_rcw(varname='u13')
      call data3d_rcw(varname='u22')
      call data3d_rcw(varname='u23')
      call data3d_rcw(varname='u33')
      !
    elseif(fname=='outdat/3rdsta') then
      !
      call h5sread(varname='nsamples',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='nsamples',filename=fname//'.h5')
      !
      call h5sread(varname='nstep_sbeg',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='nstep_sbeg',filename=fname//'.h5')
      !
      call h5sread(varname='time_sbeg',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='time_sbeg',filename=fname//'.h5')
      !
      call data3d_rcw(varname='u111')
      call data3d_rcw(varname='u112')
      call data3d_rcw(varname='u113')
      call data3d_rcw(varname='u122')
      call data3d_rcw(varname='u123')
      call data3d_rcw(varname='u133')
      call data3d_rcw(varname='u222')
      call data3d_rcw(varname='u223')
      call data3d_rcw(varname='u233')
      call data3d_rcw(varname='u333')
      !
    elseif(fname=='outdat/budget') then
      !
      call h5sread(varname='nsamples',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='nsamples',filename=fname//'.h5')
      !
      call h5sread(varname='nstep_sbeg',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='nstep_sbeg',filename=fname//'.h5')
      !
      call h5sread(varname='time_sbeg',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='time_sbeg',filename=fname//'.h5')
      !
      call data3d_rcw(varname= 'predil')
      call data3d_rcw(varname=    'pu1')
      call data3d_rcw(varname=    'pu2')
      call data3d_rcw(varname=    'pu3')
      call data3d_rcw(varname='sgmam11')
      call data3d_rcw(varname='sgmam12')
      call data3d_rcw(varname='sgmam13')
      call data3d_rcw(varname='sgmam22')
      call data3d_rcw(varname='sgmam23')
      call data3d_rcw(varname='sgmam33')
      call data3d_rcw(varname=  'u1rem')
      call data3d_rcw(varname=  'u2rem')
      call data3d_rcw(varname=  'u3rem')
      call data3d_rcw(varname='visdif1')
      call data3d_rcw(varname='visdif2')
      call data3d_rcw(varname='visdif3')
      !
    else
      stop ' !! file name error @ stream2struc'
    endif
    !
    contains
    !
    subroutine data3d_rcw(varname)
      !
      character(len=*),intent(in) :: varname
      !
      call h5sread(varname=varname,var= data_1d,dim=asizemax,filename=fname//'.s5')
      !
      write(*,'(A)',advance='no')'  ** data converting ... '
      do jrank=0,rankmax
        !
        n=offset(jrank)
        !
        do k=k0(jrank),k0(jrank)+km(jrank)
        do j=j0(jrank),j0(jrank)+jm(jrank)
        do i=i0(jrank),i0(jrank)+im(jrank)
          !
          n=n+1
          !
          data_3d(i,j,k)=data_1d(n)
          !
        enddo
        enddo
        enddo
        !
        write(*,'(1A1,A,I4,A,$)')char(13),'  ** data converting ... ', &
                                                100*jrank/rankmax,'  % '
        !
      enddo
      write(*,*)''
      !
      call h5srite(var=data_3d,varname=varname,filename=fname//'.h5')
      !
    end subroutine data3d_rcw
    !
  end subroutine stream2struc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine stream2struc.                           |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate an example grid.              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-05-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine gridchannel(im,jm,km,folder)
    !
    use tecio
    use hdf5io
    use constdef
    use commfunc, only : argtanh
    !
    integer,intent(in) :: im,jm,km
    character(len=*),intent(in) :: folder
    !
    real(8),allocatable,dimension(:,:,:) :: x,y,z
    integer :: npa,n,i,j,k,jmm,fh
    real(8) :: varc,var1,var2,Retau,dx
    !
    allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))
    !
    Retau=185.d0
    jmm=jm/2
    !
    varc=1.075d0
    do j=0,jm
      !
      var1=argtanh(1.d0/varc)
      var2=2.d0*(j)/(jm)*1.d0-1.d0
      !
      y(0,j,0)=1.d0*(1.d0+varc*dtanh(var1*var2))
      !
    end do
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      x(i,j,k)=4.d0*pi/im*i
      y(i,j,k)=y(0,j,0)
      !
      if(km==0) then
        z(i,j,k)=0.d0
      else
        z(i,j,k)=num4d3*pi/km*k
      endif
    end do
    end do
    end do
    print*,' ** y1+=',y(0,1,0)*Retau
    print*,' ** ym+=',(y(0,jmm,0)-y(0,jmm-1,0))*Retau
    print*,' ** dx+=',(x(1,0,0)-x(0,0,0))*Retau,(x(1,0,0)-x(0,0,0))
    print*,' ** dz+=',(z(0,0,1)-z(0,0,0))*Retau,(z(0,0,1)-z(0,0,0))
    print*,' ** lx+=',(x(im,0,0)-x(0,0,0))*Retau
    print*,' ** lz+=',(z(0,0,km)-z(0,0,0))*Retau
    !
    fh=get_unit()
    open(fh,file='dy.dat')
    write(fh,*)(j*1.d0,y(0,j,0),y(0,j,0)-y(0,j-1,0),j=1,jm)
    close(fh)
    print*,' << dy.dat'
    !
    ! call writetecbin('Testout/tecgrid.plt',x,'x',y,'y',z,'z')
    !
    call h5srite(x,'x',folder//'/datin/grid.chl',newfile=.true.)
    call h5srite(y,'y',folder//'/datin/grid.chl')
    call h5srite(z,'z',folder//'/datin/grid.chl')
    !
    !
    deallocate(x,y,z)
    !
  end subroutine gridchannel
  !+-------------------------------------------------------------------+
  !| The end of the subroutine gridchannel.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to process solid file to improve its accuracy. |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solidpp
    !
    use commtype, only : solid,triangle
    use readwrite,only : readsolid
    use commvar,  only : nsolid,immbody
    use geom,     only : solidrange,solidresc,solidshif,solidimpro,solidrota
    use tecio,    only : tecsolid
    use stlaio,   only : stla_write
    use cmdefne,  only : readkeyboad
    !
    ! arguments
    !
    ! local data
    integer :: js
    character(len=64) :: inputfile
    character(len=4) :: cmd
    real(8) :: resc_fact
    real(8) :: rot_vec(3),rot_theta,shift_cor(3)
    !
    call readkeyboad(cmd)
    !
    print*,cmd
    !
    if(cmd=='sgen') then
      call solidgen_circle
      !
      shift_cor=(/5.d0,2.5d0,0.d0/)
      !
    elseif(cmd=='proc') then
      call readkeyboad(inputfile)
      !
      call readsolid(trim(inputfile))
      !
      ! resc_fact=0.015d0
      resc_fact=0.005d0
      rot_vec=(/1.d0,0.d0,0.d0/)
      rot_theta=-90.d0
      shift_cor=(/5.d0,5.d0,5.d0/)
    else
      !
      print*,' !! ERROR 1, cmd not defined !!'
      print*,' ** cmd=',cmd
      print*,' ** inputfile',inputfile
      !
      stop
      !
    endif
    !
    ! 
    ! call solidgen_circle
    ! call solidgen_cub
    ! call solidgen_triagnle
    ! call solidgen_airfoil
    !
    do js=1,nsolid
      call solidrange(immbody(js))
      !
      ! call solidresc(immbody(js),resc_fact)
      ! call solidrota(immbody(js),rot_theta,rot_vec)
      call solidshif(immbody(js),x=shift_cor(1)-immbody(js)%xcen(1),  &
                                 y=shift_cor(2)-immbody(js)%xcen(2),  &
                                 z=shift_cor(3)-immbody(js)%xcen(3))
      !
    enddo
    !
    !
    !
    !
    ! do js=1,nsolid
    !   call solidimpro(immbody(js))
    ! enddo
    !
    call tecsolid('tecsolid.plt',immbody)
    !
    call stla_write('solid_in.stl',immbody)
    !
  end subroutine solidpp
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidpp.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate a solid usign stl file format.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solidgen_airfoil
    !
    use commtype, only : solid,triangle
    use readwrite,only : readsolid
    use commvar,  only : nsolid,immbody
    use geom,     only : solidrange,solidresc,solidshif,solidimpro
    use tecio,     only : tecsolid
    use stlaio,    only : stla_write
    use commfunc,  only : cross_product,argtanh
    !
    ! local data
    integer :: i,j,k,jf,km,jm,im,nface
    real(8) :: dthe,dphi,theter1,theter2,phi1,phi2,radi, &
               x1(3),x2(3),x3(3),x4(3),x5(3),x6(3),norm1(3),var1,var2, &
               xc,yc
    real(8) :: epsilon,varc,theter
    integer :: map
    character(len=4) :: nacaname
    real(8),allocatable :: xap(:),yap(:),xin(:)
    type(triangle),allocatable :: tempface(:)
    !
    epsilon=1.d-12
    !
    print*,' ** generating solid'
    !
    map=1024
    nacaname='4412'
    allocate(xin(0:map),xap(0:map),yap(0:map))
    !
    varc=1.02d0
    do i=0,map/2
      !
      var1=argtanh(1.d0/varc)
      var2=2.d0*(i)/(map)*1.d0-1.d0
      !
      xin(i)=1.d0*(1.d0+varc*dtanh(var1*var2))
      call naca4digit(xin(i),xap(i),yap(i),nacaname,'upper')
      !
    end do
    !
    do i=map/2+1,map
      xin(i)=xin(map-i)
      call naca4digit(xin(i),xap(i),yap(i),nacaname,'lower')
    enddo
    !
    print*,' ** input angle of attack in degree'
    read(*,*)theter
    !
    theter=-theter/180.d0*pi
    !
    do i=0,map
      var1=xap(i)*cos(theter)-yap(i)*sin(theter)
      var2=xap(i)*sin(theter)+yap(i)*cos(theter)
      xap(i)=var1
      yap(i)=var2
    enddo
    !
    open(18,file='naca'//nacaname//'.dat')
    do i=0,map
      write(18,*)xap(i),yap(i)
    enddo
    close(18)
    print*,' << naca',nacaname,'.dat'
    !
    ! open(12,file='naca4412.dat')
    ! read(12,*)map
    ! allocate(xap(map),yap(map))
    ! do i=1,map
    !   read(12,*)xap(i),yap(i)
    ! enddo
    ! close(12)
    ! print*,' << naca4412.dat'
    !
    nface=0
    !
    xc=0.d0
    yc=0.d0
    do i=1,map-1
      xc=xc+xap(i)
      yc=yc+yap(i)
    enddo
    !
    xc=xc/dble(map-1)
    yc=yc/dble(map-1)
    !
    nsolid=1
    allocate(immbody(nsolid))
    immbody(1)%name='cube'
    !
    allocate(tempface(map*16))
    !
    do i=0,map-1
      x1(1)=xap(i)
      x1(2)=yap(i)
      x1(3)=0.d0
      !
      x2(1)=xap(i+1)
      x2(2)=yap(i+1)
      x2(3)=0.d0
      !
      x3(1)=xap(i+1)
      x3(2)=yap(i+1)
      x3(3)=1.d0
      !
      x4(1)=xap(i)
      x4(2)=yap(i)
      x4(3)=1.d0
      !
      x5(1)=xc
      x5(2)=yc
      x5(3)=0.d0
      !
      x6(1)=xc
      x6(2)=yc
      x6(3)=1.d0
      !
      nface=nface+1
      tempface(nface)%a=x1
      tempface(nface)%b=x2
      tempface(nface)%c=x3
      tempface(nface)%normdir=cross_product(x3-x1,x2-x1)
      !
      nface=nface+1
      tempface(nface)%a=x3
      tempface(nface)%b=x4
      tempface(nface)%c=x1
      tempface(nface)%normdir=cross_product(x1-x3,x4-x3)
      !
      nface=nface+1
      tempface(nface)%a=x1
      tempface(nface)%b=x2
      tempface(nface)%c=x5
      tempface(nface)%normdir=(/0.d0,0.d0,-1.d0/)
      !
      nface=nface+1
      tempface(nface)%a=x4
      tempface(nface)%b=x3
      tempface(nface)%c=x6
      tempface(nface)%normdir=(/0.d0,0.d0,1.d0/)
      !
    enddo
    !
    immbody(1)%num_face=nface
    call immbody(1)%alloface()
    immbody(1)%face(1:nface)=tempface(1:nface)
    !
    call tecsolid('tecsolid.plt',immbody)
    !
    print*,' ** solid generated'
    !
  end subroutine solidgen_airfoil
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidgen_airfoil.                       |
  !+-------------------------------------------------------------------+
  !
  
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate a mvg solid usign stl file format. |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 22-10-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solidgen_mvg
    !
    use commtype, only : solid,triangle
    use readwrite,only : readsolid
    use commvar,  only : nsolid,immbody
    use geom,     only : solidrange,solidresc,solidshif,solidimpro
    use tecio,     only : tecsolid
    use stlaio,    only : stla_write
    use commfunc,  only : cross_product
    !
    ! local data
    integer :: i,j,k,jf,km,jm,im,nface
    real(8) :: dthe,dphi,theter1,theter2,phi1,phi2,radi, &
               x1(3),x2(3),x3(3),x4(3),norm1(3),var1,var2
    real(8) :: epsilon,delta,x0,lx,lz,h,sz,chord
    type(triangle),allocatable :: tempface(:)
    !
    epsilon=1.d-12
    delta=1.d0
    x0=10.d0
    lx=2.523597432d0
    lz=2.25d0
    chord=sqrt(lx**2+0.25d0*lz**2)
    sz=0.25d0
    h =0.3838d0
    nface=0
    !
    print*,' ** generating solid'
    !
    nsolid=2
    allocate(immbody(nsolid))
    immbody(1)%name='mvg1'
    !
    allocate(tempface(12))
    !
    x1(1)=x0
    x1(2)=0.d0
    x1(3)=0.5d0*sz
    !
    x2(1)=x0
    x2(2)=0.d0
    x2(3)=0.5d0*sz+lz
    !
    x3(1)=x0+lx
    x3(2)=0.d0
    x3(3)=0.5d0*sz+0.5d0*lz
    !
    x4(1)=x0+lx
    x4(2)=0.d0+h
    x4(3)=0.5d0*sz+0.5d0*lz
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/0.d0,-1.d0,0.d0/)
    !
    var1=0.5d0*lz/chord
    var2=      lx/chord
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x3
    tempface(nface)%c=x4
    tempface(nface)%normdir=(/var1,0.d0,-var2/)
    !
    var1=0.5d0*lz/chord
    var2=      lx/chord
    nface=nface+1
    tempface(nface)%a=x2
    tempface(nface)%b=x3
    tempface(nface)%c=x4
    tempface(nface)%normdir=(/var1,0.d0,var2/)
    !
    var1=h /sqrt(lx**2+h**2)
    var2=lx/sqrt(lx**2+h**2)
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x4
    tempface(nface)%normdir=(/-var1,var2,0.d0/)
    !
    immbody(1)%num_face=nface
    call immbody(1)%alloface()
    immbody(1)%face(1:nface)=tempface(1:nface)
    !
    nface=0
    !
    x1(3)=x1(3)+sz+lz
    x2(3)=x2(3)+sz+lz
    x3(3)=x3(3)+sz+lz
    x4(3)=x4(3)+sz+lz
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/0.d0,-1.d0,0.d0/)
    !
    var1=0.5d0*lz/chord
    var2=      lx/chord
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x3
    tempface(nface)%c=x4
    tempface(nface)%normdir=(/var1,0.d0,-var2/)
    !
    var1=0.5d0*lz/chord
    var2=      lx/chord
    nface=nface+1
    tempface(nface)%a=x2
    tempface(nface)%b=x3
    tempface(nface)%c=x4
    tempface(nface)%normdir=(/var1,0.d0, var2/)
    !
    var1=h /sqrt(lx**2+h**2)
    var2=lx/sqrt(lx**2+h**2)
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x4
    tempface(nface)%normdir=(/-var1,var2,0.d0/)
    !
    immbody(2)%num_face=nface
    call immbody(2)%alloface()
    immbody(2)%face(1:nface)=tempface(1:nface)
    !
    call tecsolid('tecsolid.plt',immbody)
    !
    print*,' ** solid generated'
    !
  end subroutine solidgen_mvg
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate a solid usign stl file format.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solidgen_triagnle
    !
    use commtype, only : solid,triangle
    use readwrite,only : readsolid
    use commvar,  only : nsolid,immbody
    use geom,     only : solidrange,solidresc,solidshif,solidimpro
    use tecio,     only : tecsolid
    use stlaio,    only : stla_write
    use commfunc,  only : cross_product
    !
    ! local data
    integer :: i,j,k,jf,km,jm,im,nface
    real(8) :: dthe,dphi,theter1,theter2,phi1,phi2,radi, &
               x1(3),x2(3),x3(3),x4(3),norm1(3),var1
    real(8) :: epsilon
    type(triangle),allocatable :: tempface(:)
    !
    epsilon=1.d-12
    !
    print*,' ** generating solid'
    !
    nsolid=1
    allocate(immbody(nsolid))
    immbody(1)%name='cube'
    !
    allocate(tempface(12))
    !
    nface=0
    !
    x1(1)=0.d0
    x1(2)=0.d0
    x1(3)=0.d0
    !
    x2(1)=1.d0
    x2(2)=-0.5d0
    x2(3)=0.d0
    !
    x3(1)=1.d0
    x3(2)=0.5d0
    x3(3)=0.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/0.d0,0.d0,-1.d0/)
    !
    x1(1)=0.d0
    x1(2)=0.d0
    x1(3)=1.d0
    !
    x2(1)=1.d0
    x2(2)=-0.5d0
    x2(3)=1.d0
    !
    x3(1)=1.d0
    x3(2)=0.5d0
    x3(3)=1.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/0.d0,0.d0, 1.d0/)
    !
    x1(1)=1.d0
    x1(2)=-0.5d0
    x1(3)=0.d0
    !
    x2(1)=1.d0
    x2(2)=-0.5d0
    x2(3)=1.d0
    !
    x3(1)=1.d0
    x3(2)=0.5d0
    x3(3)=1.d0
    !
    x4(1)=1.d0
    x4(2)=0.5d0
    x4(3)=0.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/1.d0,0.d0,0.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    tempface(nface)%normdir=(/1.d0,0.d0,0.d0/)
    !
    x1(1)=0.d0
    x1(2)=0.d0
    x1(3)=0.d0
    !
    x2(1)=0.d0
    x2(2)=0.d0
    x2(3)=1.d0
    !
    x3(1)=1.d0
    x3(2)=0.5d0
    x3(3)=1.d0
    !
    x4(1)=1.d0
    x4(2)=0.5d0
    x4(3)=0.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/-1.d0/sqrt(5.d0),2.d0/sqrt(5.d0),0.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    tempface(nface)%normdir=(/-1.d0/sqrt(5.d0),2.d0/sqrt(5.d0),0.d0/)
    !
    x1(1)=0.d0
    x1(2)=0.d0
    x1(3)=0.d0
    !
    x2(1)=0.d0
    x2(2)=0.d0
    x2(3)=1.d0
    !
    x3(1)=1.d0
    x3(2)=-0.5d0
    x3(3)=1.d0
    !
    x4(1)=1.d0
    x4(2)=-0.5d0
    x4(3)=0.d0
    !
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/-1.d0/sqrt(5.d0),-2.d0/sqrt(5.d0),0.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    tempface(nface)%normdir=(/-1.d0/sqrt(5.d0),-2.d0/sqrt(5.d0),0.d0/)
    !
    !
    immbody(1)%num_face=nface
    call immbody(1)%alloface()
    immbody(1)%face(1:nface)=tempface(1:nface)
    !
    call tecsolid('tecsolid.plt',immbody)
    !
    print*,' ** solid generated'
    !
  end subroutine solidgen_triagnle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidgen_triagnle.                      |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate a solid usign stl file format.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solidgen_cub
    !
    use commtype, only : solid,triangle
    use readwrite,only : readsolid
    use commvar,  only : nsolid,immbody
    use geom,     only : solidrange,solidresc,solidshif,solidimpro
    use tecio,     only : tecsolid
    use stlaio,    only : stla_write
    use commfunc,  only : cross_product
    !
    ! local data
    integer :: i,j,k,jf,km,jm,im,nface
    real(8) :: dthe,dphi,theter1,theter2,phi1,phi2,radi, &
               x1(3),x2(3),x3(3),x4(3),norm1(3),var1
    real(8) :: epsilon
    type(triangle),allocatable :: tempface(:)
    !
    epsilon=1.d-12
    !
    print*,' ** generating solid'
    !
    nsolid=1
    allocate(immbody(nsolid))
    immbody(1)%name='cube'
    !
    allocate(tempface(12))
    nface=0
    !
    x1(1)=0.d0
    x1(2)=0.d0
    x1(3)=0.d0
    !
    x2(1)=1.d0
    x2(2)=0.d0
    x2(3)=0.d0
    !
    x3(1)=1.d0
    x3(2)=1.d0
    x3(3)=0.d0
    !
    x4(1)=0.d0
    x4(2)=1.d0
    x4(3)=0.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/0.d0,0.d0,-1.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    tempface(nface)%normdir=(/0.d0,0.d0,-1.d0/)
    !
    x1(1)=1.d0
    x1(2)=0.d0
    x1(3)=0.d0
    !
    x2(1)=1.d0
    x2(2)=0.d0
    x2(3)=1.d0
    !
    x3(1)=1.d0
    x3(2)=1.d0
    x3(3)=1.d0
    !
    x4(1)=1.d0
    x4(2)=1.d0
    x4(3)=0.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/1.d0,0.d0,0.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    tempface(nface)%normdir=(/1.d0,0.d0,0.d0/)
    !
    x1(1)=0.d0
    x1(2)=0.d0
    x1(3)=0.d0
    !
    x2(1)=0.d0
    x2(2)=0.d0
    x2(3)=1.d0
    !
    x3(1)=0.d0
    x3(2)=1.d0
    x3(3)=1.d0
    !
    x4(1)=0.d0
    x4(2)=1.d0
    x4(3)=0.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/-1.d0,0.d0,0.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    tempface(nface)%normdir=(/-1.d0,0.d0,0.d0/)
    !
    x1(1)=0.d0
    x1(2)=0.d0
    x1(3)=1.d0
    !
    x2(1)=1.d0
    x2(2)=0.d0
    x2(3)=1.d0
    !
    x3(1)=1.d0
    x3(2)=1.d0
    x3(3)=1.d0
    !
    x4(1)=0.d0
    x4(2)=1.d0
    x4(3)=1.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/0.d0,0.d0,1.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    tempface(nface)%normdir=(/0.d0,0.d0,1.d0/)
    !
    x1(1)=0.d0
    x1(2)=0.d0
    x1(3)=0.d0
    !
    x2(1)=1.d0
    x2(2)=0.d0
    x2(3)=0.d0
    !
    x3(1)=1.d0
    x3(2)=0.d0
    x3(3)=1.d0
    !
    x4(1)=0.d0
    x4(2)=0.d0
    x4(3)=1.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/0.d0,-1.d0,0.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    tempface(nface)%normdir=(/0.d0,-1.d0,0.d0/)
    !
    x1(1)=0.d0
    x1(2)=1.d0
    x1(3)=0.d0
    !
    x2(1)=1.d0
    x2(2)=1.d0
    x2(3)=0.d0
    !
    x3(1)=1.d0
    x3(2)=1.d0
    x3(3)=1.d0
    !
    x4(1)=0.d0
    x4(2)=1.d0
    x4(3)=1.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/0.d0,1.d0,0.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    
    tempface(nface)%normdir=(/0.d0,1.d0,0.d0/)
    !
    immbody(1)%num_face=nface
    call immbody(1)%alloface()
    immbody(1)%face(1:nface)=tempface(1:nface)
    !
    call tecsolid('tecsolid.plt',immbody)
    !
    print*,' ** solid generated'
    !
  end subroutine solidgen_cub
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidgen_cub.                           |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate a sphere STL file with many        |
  !| equilateral triangle.                                             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solidgen_sphere_tri
    !
    use commtype, only : solid,triangle
    use readwrite,only : readsolid
    use commvar,  only : nsolid,immbody
    use geom,     only : solidrange,solidresc,solidshif,solidimpro
    use tecio,     only : tecsolid
    use stlaio,    only : stla_write
    use commfunc,  only : cross_product,dis2point
    !
    ! local data
    integer :: i,j,k,jf,km,jm,im,nface,nface2,jface,ntrimax,m,n,l
    real(8) :: dthe,dphi,theter1,theter2,phi1,phi2,radi,   &
               x1(3),x2(3),x3(3),x4(3),norm1(3),norm2(3),var1,var2, &
               ab(3),ac(3),bc(3),abc(3),xico(12,3),disab,disac,disbc
    integer :: lmnsav(120,3)
    real(8) :: epsilon
    type(triangle),allocatable :: tempface(:),tempface2(:)
    !
    epsilon=1.d-12
    ntrimax=10000
    !
    print*,' ** generating solid'
    !
    nsolid=1
    allocate(immbody(nsolid))
    immbody(1)%name='sphere'
    !
    allocate(tempface(ntrimax),tempface2(2*ntrimax))
    !
    ! first build a icosahedron
    !
    var1=0.5d0*(sqrt(5.d0)+1.d0)
    !
    xico(1,1)=0.d0; xico(2,1)= 0.d0; xico(3,1)= 0.d0; xico(4,1)= 0.d0
    xico(1,2)=1.d0; xico(2,2)= 1.d0; xico(3,2)=-1.d0; xico(4,2)=-1.d0
    xico(1,3)=var1; xico(2,3)=-var1; xico(3,3)=-var1; xico(4,3)= var1
    !
    xico(5,1)=1.d0; xico(6,1)= 1.d0; xico(7,1)=-1.d0; xico(8,1)=-1.d0
    xico(5,2)=var1; xico(6,2)=-var1; xico(7,2)=-var1; xico(8,2)= var1
    xico(5,3)=0.d0; xico(6,3)= 0.d0; xico(7,3)= 0.d0; xico(8,3)= 0.d0
    !
    xico(9,1)=var1; xico(10,1)=-var1; xico(11,1)=-var1; xico(12,1)= var1
    xico(9,2)=0.d0; xico(10,2)= 0.d0; xico(11,2)= 0.d0; xico(12,2)= 0.d0
    xico(9,3)=1.d0; xico(10,3)= 1.d0; xico(11,3)=-1.d0; xico(12,3)=-1.d0
    !
    radi=sqrt(xico(1,1)**2+xico(1,2)**2+xico(1,3)**2)
    print*,' ** radi= ',radi
    xico=xico/radi
    var2=2.d0/radi
    !
    nface=0
    lmnsav=0.d0
    do l=1,12
    do m=1,12
    do n=1,12
      !
      disab=abs(dis2point(xico(m,:),xico(n,:))-var2)
      disac=abs(dis2point(xico(l,:),xico(n,:))-var2)
      disbc=abs(dis2point(xico(m,:),xico(l,:))-var2)
      !
      if(disab<epsilon .and. disac<epsilon .and. disbc<epsilon) then
        !
        do i=1,nface
          !
          if( (lmnsav(i,1)==l .and. lmnsav(i,2)==m .and. lmnsav(i,3)==n) .or. &
              (lmnsav(i,1)==l .and. lmnsav(i,2)==n .and. lmnsav(i,3)==m) .or. &
              (lmnsav(i,1)==m .and. lmnsav(i,2)==l .and. lmnsav(i,3)==n) .or. &
              (lmnsav(i,1)==m .and. lmnsav(i,2)==n .and. lmnsav(i,3)==l) .or. &
              (lmnsav(i,1)==n .and. lmnsav(i,2)==l .and. lmnsav(i,3)==m) .or. &
              (lmnsav(i,1)==n .and. lmnsav(i,2)==m .and. lmnsav(i,3)==l) ) then
            exit
          endif
        enddo
        !
        if(i==nface+1) then
          !
          nface=nface+1
          !
          print*,' ** nface= ',nface
          !
          lmnsav(nface,1)=l
          lmnsav(nface,2)=m
          lmnsav(nface,3)=n
          !
          tempface(nface)%a=xico(l,:)
          tempface(nface)%b=xico(m,:)
          tempface(nface)%c=xico(n,:)
          !
          norm1=cross_product(tempface(nface)%b-tempface(nface)%a,         &
                              tempface(nface)%c-tempface(nface)%a )
          norm2=tempface(nface)%a
          !
          if(dot_product(norm1,norm2)<0.d0) then
            !
            tempface(nface)%a=xico(l,:)
            tempface(nface)%b=xico(n,:)
            tempface(nface)%c=xico(m,:)
            !
            norm1=cross_product(tempface(nface)%b-tempface(nface)%a,         &
                                tempface(nface)%c-tempface(nface)%a )
            norm2=tempface(nface)%a
            !
          endif
          !
          !
        endif
        !
      endif
      !
    enddo
    enddo
    enddo
    !
    do while(nface<ntrimax)
      !
      nface2=0
      !
      do jface=1,nface
        !
        ab=0.5d0*(tempface(jface)%a+tempface(jface)%b)
        ac=0.5d0*(tempface(jface)%a+tempface(jface)%c)
        bc=0.5d0*(tempface(jface)%b+tempface(jface)%c)
        abc=num1d3*(tempface(jface)%a+tempface(jface)%b+tempface(jface)%c)
        !
        var1=sqrt(ab(1)**2+ab(2)**2+ab(3)**2)
        ab=ab/var1
        var1=sqrt(ac(1)**2+ac(2)**2+ac(3)**2)
        ac=ac/var1
        var1=sqrt(bc(1)**2+bc(2)**2+bc(3)**2)
        bc=bc/var1
        var1=sqrt(abc(1)**2+abc(2)**2+abc(3)**2)
        abc=abc/var1
        !
        nface2=nface2+1
        tempface2(nface2)%a=tempface(jface)%a
        tempface2(nface2)%b=ab
        tempface2(nface2)%c=ac
        !
        nface2=nface2+1
        tempface2(nface2)%a=tempface(jface)%b
        tempface2(nface2)%b=ab
        tempface2(nface2)%c=bc
        !
        nface2=nface2+1
        tempface2(nface2)%a=tempface(jface)%c
        tempface2(nface2)%b=ac
        tempface2(nface2)%c=bc
        !
        nface2=nface2+1
        tempface2(nface2)%a=ab
        tempface2(nface2)%b=ac
        tempface2(nface2)%c=bc
        !
        if(nface2>ntrimax) exit
        !
      enddo
      !
      if(nface2>ntrimax) exit
      !
      nface=nface2
      tempface=tempface2
      !
      print*,' ** nface=',nface
      !
    enddo
    !
    !
    immbody(1)%num_face=nface
    call immbody(1)%alloface()
    immbody(1)%face(1:nface)=tempface(1:nface)
    !
    do jf=1,immbody(1)%num_face
      !
      norm1=cross_product(immbody(1)%face(jf)%a-immbody(1)%face(jf)%b,         &
                          immbody(1)%face(jf)%c-immbody(1)%face(jf)%b )
      norm2=immbody(1)%face(jf)%a
      !
      if(dot_product(norm1,norm2)<0.d0) then
        x1=immbody(1)%face(jf)%c
        immbody(1)%face(jf)%c=immbody(1)%face(jf)%a
        immbody(1)%face(jf)%a=x1
        !
        norm1=cross_product(immbody(1)%face(jf)%a-immbody(1)%face(jf)%b,       &
                            immbody(1)%face(jf)%c-immbody(1)%face(jf)%b )
      endif
      !
      var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
      !
      if(abs(var1)<epsilon) then
        !
        norm1=cross_product(immbody(1)%face(jf)%C-immbody(1)%face(jf)%a,         &
                            immbody(1)%face(jf)%b-immbody(1)%face(jf)%a )
        !
        var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
        !
      endif
      !
      if(abs(var1)<epsilon) then
        !
        norm1=cross_product(immbody(1)%face(jf)%a-immbody(1)%face(jf)%c,         &
                            immbody(1)%face(jf)%b-immbody(1)%face(jf)%c )
        !
        var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
        !
      endif
      !
      immbody(1)%face(jf)%normdir=norm1/var1
      !
      if(isnan(immbody(1)%face(jf)%normdir(1))) then
        print*,jf
        print*,immbody(1)%face(jf)%a
        print*,immbody(1)%face(jf)%b
        print*,immbody(1)%face(jf)%c
        print*,norm1,'|',var1
        stop
      endif
      !
    enddo
    !
    call tecsolid('tecsolid.plt',immbody)
    !
    print*,' ** solid generated'
    !
  end subroutine solidgen_sphere_tri
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidgen_sphere_tri.                    |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate a solid usign stl file format.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solidgen_sphere
    !
    use commtype, only : solid,triangle
    use readwrite,only : readsolid
    use commvar,  only : nsolid,immbody
    use geom,     only : solidrange,solidresc,solidshif,solidimpro
    use tecio,     only : tecsolid
    use stlaio,    only : stla_write
    use commfunc,  only : cross_product
    !
    ! local data
    integer :: i,j,k,jf,km,jm,im,nface
    real(8) :: dthe,dphi,theter1,theter2,phi1,phi2,radi, &
               x1(3),x2(3),x3(3),x4(3),norm1(3),var1
    real(8) :: epsilon
    type(triangle),allocatable :: tempface(:)
    !
    epsilon=1.d-12
    !
    print*,' ** generating solid'
    !
    nsolid=1
    allocate(immbody(nsolid))
    immbody(1)%name='sphere'
    !
    km=16
    jm=32
    !
    allocate(tempface(km*jm*2))
    !
    dthe=2.d0*pi/km
    dphi=pi/jm
    radi=0.5d0
    !
    nface=0
    !
    do k=0,km-1
    do j=0,jm-1
      !
      theter1=dthe*k
      theter2=dthe*(k+1)
      phi1=dphi*j
      phi2=dphi*(j+1)
      !
      x1(1)=radi*sin(phi1)*cos(theter1)
      x1(2)=radi*sin(phi1)*sin(theter1)
      x1(3)=radi*cos(phi1)
      !
      x2(1)=radi*sin(phi2)*cos(theter1)
      x2(2)=radi*sin(phi2)*sin(theter1)
      x2(3)=radi*cos(phi2)
      !
      x3(1)=radi*sin(phi2)*cos(theter2)
      x3(2)=radi*sin(phi2)*sin(theter2)
      x3(3)=radi*cos(phi2)
      !
      x4(1)=radi*sin(phi1)*cos(theter2)
      x4(2)=radi*sin(phi1)*sin(theter2)
      x4(3)=radi*cos(phi1)
      !
      if(abs(sin(phi1))<epsilon) then
        !
        nface=nface+1
        tempface(nface)%a=x1
        tempface(nface)%b=x2
        tempface(nface)%c=x3
        !
      elseif(abs(sin(phi2))<epsilon) then
        !
        nface=nface+1
        tempface(nface)%a=x3
        tempface(nface)%b=x4
        tempface(nface)%c=x1
        !
      else
        !
        nface=nface+1
        tempface(nface)%a=x1
        tempface(nface)%b=x2
        tempface(nface)%c=x3
        !
        nface=nface+1
        tempface(nface)%a=x3
        tempface(nface)%b=x4
        tempface(nface)%c=x1
        !
      endif
      !
    enddo
    enddo
    !
    immbody(1)%num_face=nface
    call immbody(1)%alloface()
    immbody(1)%face(1:nface)=tempface(1:nface)
    !
    do jf=1,immbody(1)%num_face
      !
      norm1=cross_product(immbody(1)%face(jf)%b-immbody(1)%face(jf)%a,         &
                          immbody(1)%face(jf)%c-immbody(1)%face(jf)%a )
      !
      var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
      !
      if(abs(var1)<epsilon) then
        !
        norm1=cross_product(immbody(1)%face(jf)%a-immbody(1)%face(jf)%b,         &
                            immbody(1)%face(jf)%c-immbody(1)%face(jf)%b )
        !
        var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
        !
      endif
      !
      if(abs(var1)<epsilon) then
        !
        norm1=cross_product(immbody(1)%face(jf)%a-immbody(1)%face(jf)%c,         &
                            immbody(1)%face(jf)%b-immbody(1)%face(jf)%c )
        !
        var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
        !
      endif
      !
      immbody(1)%face(jf)%normdir=norm1/var1
      !
      if(isnan(immbody(1)%face(jf)%normdir(1))) then
        print*,jf
        print*,immbody(1)%face(jf)%a
        print*,immbody(1)%face(jf)%b
        print*,immbody(1)%face(jf)%c
        print*,norm1,'|',var1
        stop
      endif
      !
    enddo
    !
    call tecsolid('tecsolid.plt',immbody)
    !
    print*,' ** solid generated'
    !
  end subroutine solidgen_sphere
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidgen_sphere.                        |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate a solid usign stl file format.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 06-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solidgen_circle
    !
    use commtype, only : solid,triangle
    use readwrite,only : readsolid
    use commvar,  only : nsolid,immbody
    use geom,     only : solidrange,solidresc,solidshif,solidimpro
    use tecio,     only : tecsolid
    use stlaio,    only : stla_write
    use commfunc,  only : cross_product
    !
    ! local data
    integer :: i,j,k,jf,km,jm,im,nface
    real(8) :: dthe,dphi,theter1,theter2,phi1,phi2,radi, &
               x1(3),x2(3),x3(3),x4(3),x5(3),x6(3),norm1(3),var1,xcent(3)
    real(8) :: epsilon
    type(triangle),allocatable :: tempface(:)
    !
    epsilon=1.d-12
    !
    print*,' ** generating solid'
    !
    nsolid=1
    allocate(immbody(nsolid))
    immbody(1)%name='circle'
    !
    km=360
    !
    allocate(tempface(km*4))
    !
    dthe=2.d0*pi/km
    radi=0.5d0
    !
    nface=0
    !
    do k=0,km-1
      !
      theter1=dthe*k
      theter2=dthe*(k+1)
      !
      x1(1)=radi*cos(theter1)
      x1(2)=radi*sin(theter1)
      x1(3)=-0.01d0
      !
      x2(1)=radi*cos(theter2)
      x2(2)=radi*sin(theter2)
      x2(3)=-0.01d0
      !
      x3(1)=radi*cos(theter2)
      x3(2)=radi*sin(theter2)
      x3(3)=0.01d0
      !
      x4(1)=radi*cos(theter1)
      x4(2)=radi*sin(theter1)
      x4(3)=0.01d0
      !
      x5(1)=0.d0
      x5(2)=0.d0
      x5(3)=-0.01d0
      !
      x6(1)=0.d0
      x6(2)=0.d0
      x6(3)=0.01d0
      !
      nface=nface+1
      tempface(nface)%a=x1
      tempface(nface)%b=x2
      tempface(nface)%c=x3
      !
      nface=nface+1
      tempface(nface)%a=x3
      tempface(nface)%b=x4
      tempface(nface)%c=x1
      !
      nface=nface+1
      xcent=(/0.d0,0.d0,-0.01d0/)
      tempface(nface)%a=x1
      tempface(nface)%b=x2
      tempface(nface)%c=xcent
      ! !
      nface=nface+1
      xcent=(/0.d0,0.d0,0.01d0/)
      tempface(nface)%a=x3
      tempface(nface)%b=x4
      tempface(nface)%c=xcent
      !
    enddo
    !
    immbody(1)%num_face=nface
    call immbody(1)%alloface()
    immbody(1)%face(1:nface)=tempface(1:nface)
    !
    do jf=1,immbody(1)%num_face
      !
      norm1=cross_product(immbody(1)%face(jf)%b-immbody(1)%face(jf)%a,         &
                          immbody(1)%face(jf)%c-immbody(1)%face(jf)%a )
      !
      var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
      !
      if(abs(var1)<epsilon) then
        !
        norm1=cross_product(immbody(1)%face(jf)%a-immbody(1)%face(jf)%b,         &
                            immbody(1)%face(jf)%c-immbody(1)%face(jf)%b )
        !
        var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
        !
      endif
      !
      if(abs(var1)<epsilon) then
        !
        norm1=cross_product(immbody(1)%face(jf)%a-immbody(1)%face(jf)%c,         &
                            immbody(1)%face(jf)%b-immbody(1)%face(jf)%c )
        !
        var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
        !
      endif
      !
      immbody(1)%face(jf)%normdir=norm1/var1
      !
      if(isnan(immbody(1)%face(jf)%normdir(1))) then
        print*,jf
        print*,immbody(1)%face(jf)%a
        print*,immbody(1)%face(jf)%b
        print*,immbody(1)%face(jf)%c
        print*,norm1,'|',var1
        stop
      endif
      !
    enddo
    !
    call tecsolid('tecsolid.plt',immbody)
    !
    print*,' ** solid generated'
    !
  end subroutine solidgen_circle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidgen_circle.                        |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is to generate a  4-digit NACA airfoil.             |
  !+-------------------------------------------------------------------+
  !| ref: https://en.wikipedia.org/wiki/NACA_airfoil                   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine naca4digit(xin,x,y,name,surface)
    !
    ! arguments
    real(8),intent(in) :: xin
    real(8),intent(out) :: x,y
    character(len=4),intent(in) :: name
    character(len=5),intent(in) :: surface
    !
    ! local data 
    real(8) :: yt,t,yc,theter,mr,pr
    integer :: m,p
    !
    read(name(1:1),*)m
    read(name(2:2),*)p
    read(name(3:4),*)t
    !
    t=t/100.d0
    !
    yt=5.d0*t*(0.2969d0*sqrt(xin)-0.1260d0*xin-0.3516d0*xin*xin+       &
               0.2843d0*xin**3-0.1036d0*xin**4)
    !
    if(p==0 .and. m==0) then
      ! NACA-00XX, a symmetrical 4-digit NACA airfoil
      !
      if(surface=='upper') then
        x=xin
        y=yt
      elseif(surface=='lower') then
        x=xin
        y=-yt
      else
        print*,' surface=',surface
        stop '!! ERROR 1 @ naca4digit'
      endif
      !
    else
      !
      mr=0.01d0*dble(m)
      pr=0.1d0*dble(p)
      !
      if(xin>=0.d0 .and. xin<=pr) then
        yc=mr/pr/pr*(2.d0*pr*xin-xin*xin)
        theter=atan(2.d0*mr/pr/pr*(pr-xin))
      elseif(xin>=pr .and. xin<=1.d0) then
        yc=mr/(1-pr)**2*(1.d0-2.d0*pr+2.d0*pr*xin-xin**2)
        theter=atan(2.d0*mr/(1-pr)**2*(pr-xin))
      else
        stop '!! ERROR 2 @ naca4digit'
      endif
      !
      if(surface=='upper') then
        x=xin-yt*sin(theter)
        y= yc+yt*cos(theter)
      elseif(surface=='lower') then
        x=xin+yt*sin(theter)
        y= yc-yt*cos(theter)
      else
        print*,' surface=',surface
        stop '!! ERROR 3 @ naca4digit'
      endif
      !
    endif
    !
  end subroutine naca4digit
  !+-------------------------------------------------------------------+
  !| The end of the subroutine naca4digit.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is to view flow by postprocess data.                |
  !+-------------------------------------------------------------------+
  !| ref: https://en.wikipedia.org/wiki/NACA_airfoil                   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine fieldview(flowfile,outputfile,viewmode,inputfile)
    !
    use hdf5io
    use tecio
    use WriteVTK
    !
    ! arguments
    character(len=*),intent(in) :: flowfile,outputfile,viewmode,inputfile
    !
    ! local data
    integer :: im,jm,km,num_species,jsp,i
    logical :: lihomo,ljhomo,lkhomo
    real(8) :: ref_tem,reynolds,mach
    character(len=32) :: gridfile
    !
    real(8),allocatable,dimension(:) :: x_1d,y_1d,z_1d,ro_1d,u_1d,   &
                                        v_1d,w_1d,p_1d,t_1d
    real(8),allocatable,dimension(:,:) :: spc_1d
    real(8),allocatable,dimension(:,:) :: x_2d,y_2d,z_2d,ro_2d,u_2d,   &
                                          v_2d,w_2d,p_2d,t_2d
    real(8),allocatable,dimension(:,:,:) :: x,y,z,ro,u,v,w,p,t
    real(8),allocatable :: var1d(:)
    !
    character(len=3) :: spname
    !
    print*,' ==========================readinput=========================='
    !
    open(11,file=inputfile,form='formatted',status='old')
    read(11,'(///////)')
    read(11,*)im,jm,km
    read(11,"(/)")
    read(11,*)lihomo,ljhomo,lkhomo
    read(11,'(//////////)')
    read(11,*)ref_tem,reynolds,mach
    read(11,'(///////)')
    read(11,*)num_species
    print*,' ** num_species: ',num_species
    read(11,'(//////////////////)')
    read(11,'(A)')gridfile
    close(11)
    print*,' >> ',inputfile
    !
    print*,' ** grid file: ',trim(gridfile)
    !
    if(viewmode=='xy') then
      allocate(x_2d(0:im,0:jm),y_2d(0:im,0:jm))
      call H5ReadSubset(x_2d,im,jm,km,'x',gridfile,kslice=0)
      call H5ReadSubset(y_2d,im,jm,km,'y',gridfile,kslice=0)
      !
      allocate(ro_2d(0:im,0:jm),u_2d(0:im,0:jm), v_2d(0:im,0:jm),      &
                w_2d(0:im,0:jm),p_2d(0:im,0:jm), t_2d(0:im,0:jm)       )
      !
      call h5_read2dfrom3d(ro_2d,im,jm,km,'ro',flowfile,kslice=0)
      call h5_read2dfrom3d( u_2d,im,jm,km,'u1',flowfile,kslice=0)
      call h5_read2dfrom3d( v_2d,im,jm,km,'u2',flowfile,kslice=0)
      call h5_read2dfrom3d( w_2d,im,jm,km,'u3',flowfile,kslice=0)
      call h5_read2dfrom3d( p_2d,im,jm,km, 'p',flowfile,kslice=0)
      call h5_read2dfrom3d( t_2d,im,jm,km, 't',flowfile,kslice=0)
      !
      call writeprvbin(outputfile,x_2d,'x',y_2d,'y',ro_2d,'ro',u_2d,'u',   &
                                     v_2d,'v',p_2d,'p',t_2d,'t',im,jm)
       ! call tecbin(outputfile,x_2d,'x',y_2d,'y',ro_2d,'ro',u_2d,'u',   &
       !                                   v_2d,'v',p_2d,'p',t_2d,'t')
    elseif(viewmode=='3d') then
      allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))
    elseif(viewmode=='1d') then
      !
      allocate(x_1d(0:im),ro_1d(0:im),u_1d(0:im),p_1d(0:im),t_1d(0:im))
      !
      call H5ReadSubset( x_1d,im,jm,km, 'x',gridfile,jslice=jm/2-1,kslice=-1)
      !
      call H5ReadSubset(ro_1d,im,jm,km,'ro',flowfile,jslice=jm/2-1,kslice=-1)
      call H5ReadSubset( u_1d,im,jm,km,'u1',flowfile,jslice=jm/2-1,kslice=-1)
      call H5ReadSubset( p_1d,im,jm,km, 'p',flowfile,jslice=jm/2-1,kslice=-1)
      call H5ReadSubset( t_1d,im,jm,km, 't',flowfile,jslice=jm/2-1,kslice=-1)
      !
      allocate(spc_1d(0:im,1:num_species))
      allocate(var1d(0:im))
      do jsp=1,num_species
        write(spname,'(i3.3)')jsp
        call H5ReadSubset(var1d,im,jm,km,'sp'//spname,flowfile,jslice=jm/2-1,kslice=-1)
        spc_1d(:,jsp)=var1d
      enddo
      !
      open(18,file=outputfile)
      write(18,"(8(1X,A15))")'x','ro','u','p','t','Y1','Y2','Y3'
      write(18,"(8(1X,E15.7E3))")(x_1d(i),ro_1d(i),u_1d(i),p_1d(i), &
                  t_1d(i),spc_1d(i,1),spc_1d(i,2),spc_1d(i,3),i=0,im)
      close(18)
      print*,' << ',outputfile
      !
      call h5srite(var=ro_1d,varname='ro',filename='flowini1d.h5',explicit=.true.,newfile=.true.)
      call h5srite(var= u_1d,varname='u1',filename='flowini1d.h5',explicit=.true.,newfile=.false.)
      call h5srite(var=t_1d, varname= 't',filename='flowini1d.h5',explicit=.true.,newfile=.false.)
      call h5srite(var=p_1d, varname= 'p',filename='flowini1d.h5',explicit=.true.,newfile=.false.)
      !
      do jsp=1,num_species
        write(spname,'(i3.3)')jsp
        call h5srite(var=spc_1d(:,jsp),varname='sp'//spname,filename='flowini1d.h5',explicit=.true.,newfile=.false.)
      enddo
    else
      print*,viewmode
      stop ' !! mode is not defined @ fieldview'
    endif
    !
  end subroutine fieldview
  !+-------------------------------------------------------------------+
  !| The end of the subroutine flowfieldview.                          |
  !+-------------------------------------------------------------------+
  !+-------------------------------------------------------------------+
  !| This function is to generate specific flame with postprocess data.|
  !+-------------------------------------------------------------------+
  !| ref: https://en.wikipedia.org/wiki/NACA_airfoil                   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 11-Aug-2021: Created by Yifan Xu @ Peking University              |
  !+-------------------------------------------------------------------+
  subroutine flamegen(flowfile,viewmode,inputfile)
    !
    use hdf5io
    use tecio
    use WriteVTK
    !
    ! arguments
    character(len=*),intent(in) :: flowfile,viewmode,inputfile
    !
    ! local data
    integer :: im,jm,km,num_species,jsp,i
    integer :: iflame_in,iflame_out
    logical :: lihomo,ljhomo,lkhomo
    real(8) :: ref_tem,reynolds,mach
    character(len=32) :: gridfile
    !
    real(8),allocatable,dimension(:) :: x_1d,y_1d,z_1d,ro_1d,u_1d,   &
                                        v_1d,w_1d,p_1d,t_1d
    real(8),allocatable,dimension(:,:) :: spc_1d
    real(8),allocatable,dimension(:,:) :: x_2d,y_2d,z_2d,ro_2d,u_2d,   &
                                          v_2d,w_2d,p_2d,t_2d
    real(8),allocatable,dimension(:,:,:) :: spc_2d
    real(8),allocatable,dimension(:,:,:) :: x,y,z,ro,u,v,w,p,t
    real(8),allocatable,dimension(:,:,:,:) :: spc
    real(8),allocatable :: var1d(:)
    !
    character(len=3) :: spname
    !
    print*,' ==========================readinput=========================='
    !
    open(11,file=inputfile,form='formatted',status='old')
    read(11,'(///////)')
    read(11,*)im,jm,km
    read(11,"(/)")
    read(11,*)lihomo,ljhomo,lkhomo
    read(11,'(//////////)')
    read(11,*)ref_tem,reynolds,mach
    read(11,'(///////)')
    read(11,*)num_species
    print*,' ** num_species: ',num_species
    read(11,'(//////////////////)')
    read(11,'(A)')gridfile
    close(11)
    print*,' >> ',inputfile
    !
    print*,' ** grid file: ',trim(gridfile)
    !
    iflame_out=im/4
    !
    ! if(viewmode=='3d') then
    !   allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))
    ! elseif(viewmode=='1d') then
    !   !
    allocate(x_1d(0:im),ro_1d(0:im),u_1d(0:im),p_1d(0:im),t_1d(0:im))
    !
    call H5ReadSubset( x_1d,im,jm,km, 'x',gridfile,jslice=jm/2-1,kslice=-1)
    !
    call H5ReadSubset(ro_1d,im,jm,km,'ro',flowfile,jslice=jm/2-1,kslice=-1)
    call H5ReadSubset( u_1d,im,jm,km,'u1',flowfile,jslice=jm/2-1,kslice=-1)
    call H5ReadSubset( p_1d,im,jm,km, 'p',flowfile,jslice=jm/2-1,kslice=-1)
    call H5ReadSubset( t_1d,im,jm,km, 't',flowfile,jslice=jm/2-1,kslice=-1)
    !
    allocate(spc_1d(0:im,1:num_species))
    allocate(var1d(0:im))
    do jsp=1,num_species
      write(spname,'(i3.3)')jsp
      call H5ReadSubset(var1d,im,jm,km,'sp'//spname,flowfile,jslice=jm/2-1,kslice=-1)
      spc_1d(:,jsp)=var1d
    enddo

    if(viewmode=='3d') then
      !
      allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km),ro(0:im,0:jm,0:km),u(0:im,0:jm,0:km), &
               v(0:im,0:jm,0:km),w(0:im,0:jm,0:km),p(0:im,0:jm,0:km),t(0:im,0:jm,0:km))
      allocate(spc(0:im,0:jm,0:km,1:num_species))
      !
      call H5io_init(filename=gridfile,mode='read')
      call H5Read(varname='x',var=x(0:im,0:jm,0:km),mode='h')
      call H5Read(varname='y',var=y(0:im,0:jm,0:km),mode='h')
      call H5Read(varname='z',var=z(0:im,0:jm,0:km),mode='h')
      !
      do i=1,im
        !
        if((t_1d(i-1)<=400.d0 .and. t_1d(i)>=400.d0) .or. &
            (t_1d(i-1)>=400.d0 .and. t_1d(i)<=400.d0)) then
          iflame_in=i
          exit
        endif
      enddo

      do i=-im/8,im/8
        ro(i+iflame_out,:,:)=ro_1d(i+iflame_in)
        u(i+iflame_out,:,:) =u_1d(i+iflame_in)
        v(i+iflame_out,:,:) =0.d0
        w(i+iflame_out,:,:) =0.d0
        p(i+iflame_out,:,:) =p_1d(i+iflame_in)
        t(i+iflame_out,:,:) =t_1d(i+iflame_in)
        do jsp=1,num_species
          spc(i+iflame_out,:,:,jsp)=spc_1d(i+iflame_in,jsp)
        enddo
      enddo

      do i=0,iflame_out-im/8-1
        ro(i,:,:)=ro(iflame_out-im/8,:,:)
        u(i,:,:) =u(iflame_out-im/8,:,:)
        v(i,:,:) =v(iflame_out-im/8,:,:)
        w(i,:,:) =w(iflame_out-im/8,:,:)
        p(i,:,:) =p(iflame_out-im/8,:,:)
        t(i,:,:) =t(iflame_out-im/8,:,:)
        do jsp=1,num_species
          spc(i,:,:,jsp)=spc(iflame_out-im/8,:,:,jsp)
        enddo
      enddo

      do i=iflame_out+im/8+1,im
        ro(i,:,:)=ro(iflame_out+im/8,:,:)
        u(i,:,:) =u(iflame_out+im/8,:,:)
        v(i,:,:) =v(iflame_out+im/8,:,:)
        w(i,:,:) =w(iflame_out+im/8,:,:)
        p(i,:,:) =p(iflame_out+im/8,:,:)
        t(i,:,:) =t(iflame_out+im/8,:,:)
        do jsp=1,num_species
          spc(i,:,:,jsp)=spc(iflame_out+im/8,:,:,jsp)
        enddo
      enddo
    !
      call h5srite(var= x,varname= 'x',filename='datin/flowini3d.h5',explicit=.true.,newfile=.true.)
      call h5srite(var= y,varname= 'y',filename='datin/flowini3d.h5',explicit=.true.,newfile=.false.)
      call h5srite(var= z,varname= 'z',filename='datin/flowini3d.h5',explicit=.true.,newfile=.false.)
      call h5srite(var=ro,varname='ro',filename='datin/flowini3d.h5',explicit=.true.,newfile=.false.)
      call h5srite(var= u,varname='u1',filename='datin/flowini3d.h5',explicit=.true.,newfile=.false.)
      call h5srite(var= v,varname='u2',filename='datin/flowini3d.h5',explicit=.true.,newfile=.false.)
      call h5srite(var= w,varname='u3',filename='datin/flowini3d.h5',explicit=.true.,newfile=.false.)
      call h5srite(var= t,varname= 't',filename='datin/flowini3d.h5',explicit=.true.,newfile=.false.)
      call h5srite(var= p,varname= 'p',filename='datin/flowini3d.h5',explicit=.true.,newfile=.false.)
    !
      do jsp=1,num_species
        write(spname,'(i3.3)')jsp
        call h5srite(var=spc(:,:,:,jsp),varname='sp'//spname,filename='datin/flowini3d.h5',explicit=.true.,newfile=.false.)
      enddo
    !  
    else
      print*,viewmode
      stop ' !! mode is not defined @ fieldview'
    endif
    !
  end subroutine flamegen
  !+-------------------------------------------------------------------+
  !| The end of the subroutine flamegen.                          |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate a random field for homogeneous|
  !| isotropic turbulence.                                             |
  !+-------------------------------------------------------------------+
  !| Ref: Blaisdell, G. A., Numerical simulation of compressible       |
  !|      homogeneous turbulence, Phd, 1991, Stanford University       |
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 25-04-2023: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine hitgen
    !
    use cmdefne,   only : readkeyboad
    use commvar,   only : gridfile,im,jm,km,ia,ja,ka,hm,Mach,Reynolds, &
                          tinf,roinf,spcinf,num_species,nondimen
    use bc,        only : twall
    use readwrite, only : readgrid, readic, readinput
    use commarray, only : x,vel,rho,tmp,prs,dvel,spc
    use solver,    only : refcal
    use parallel,  only : mpisizedis,parapp,parallelini
    use geom,      only : geomcal
    use fludyna,   only :  thermal
    use hdf5io,    only : h5srite,h5sread
    use gridgeneration
    use tecio
    !
    real(8) :: urms,kenergy,ufmx,roav,uav,vav,wav,tav,pav
    integer :: i,j,k
    character(len=4) :: genmethod
    !
    print*,' ** cmd to generate a box turbulence: astr pp hitgen <input> . gen/read'
    !
    !
    call readinput
    !
    call readic
    !
    im=ka ! ensure a cubic box
    !
    jm=ka
    !
    km=ka
    !
    call mpisizedis
    !
    call parapp
    !
    call parallelini
    !
    call refcal
    !
    allocate(   x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    allocate( vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    allocate(rho(0:im,0:jm,0:km),tmp(0:im,0:jm,0:km),prs(0:im,0:jm,0:km))
#ifdef COMB
    allocate(spc(0:im,0:jm,0:km,1:num_species))
#endif
    !
    ! call gridhitflame(mode='cubic')
    call gridcube(2.d0*pi,2.d0*pi,2.d0*pi)
    ! call readgrid(trim(gridfile))
    !
    call geomcal
    !
    call readkeyboad(genmethod)
    !
    if(trim(genmethod)=='gen') then
      !
      call div_free_gen(im,jm,km,vel(0:im,0:jm,0:km,1),   &
                                 vel(0:im,0:jm,0:km,2),   &
                                 vel(0:im,0:jm,0:km,3) )
      !
      urms=1.d0 
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        rho(i,j,k)  = roinf
        tmp(i,j,k)  = tinf 
        if(nondimen) then
            prs(i,j,k)  = thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
        else
            spc(i,j,k,:)= spcinf(:)
            prs(i,j,k)  = thermal(density=rho(i,j,k),temperature=tmp(i,j,k),species=spc(i,j,k,:))
        endif
        vel(i,j,k,1)= urms*vel(i,j,k,1)
        vel(i,j,k,2)= urms*vel(i,j,k,2)
        vel(i,j,k,3)= urms*vel(i,j,k,3)
        !
      enddo
      enddo
      enddo
      !
      call div_test(vel,dvel)
      !
      call hitsta
      !
    elseif(trim(genmethod)=='read') then
      !
      call h5sread(vel(0:im,0:jm,0:km,1),'u1',im,jm,km,'outdat/flowfield.h5')
      call h5sread(vel(0:im,0:jm,0:km,2),'u2',im,jm,km,'outdat/flowfield.h5')
      call h5sread(vel(0:im,0:jm,0:km,3),'u3',im,jm,km,'outdat/flowfield.h5')
      call h5sread(rho(0:im,0:jm,0:km),  'ro',im,jm,km,'outdat/flowfield.h5')
      call h5sread(tmp(0:im,0:jm,0:km),  't', im,jm,km,'outdat/flowfield.h5')
      !
      call div_test(vel,dvel)
      !
      call hitsta
      !
    else
      print*,genmethod
      stop ' !! genmethod error '
    endif
    ! 
    call h5srite(var=rho,                  varname='ro',filename='flowini3d.h5',explicit=.true.,newfile=.true.) 
    call h5srite(var=vel(0:im,0:jm,0:km,1),varname='u1',filename='flowini3d.h5',explicit=.true.)
    call h5srite(var=vel(0:im,0:jm,0:km,2),varname='u2',filename='flowini3d.h5',explicit=.true.)
    call h5srite(var=vel(0:im,0:jm,0:km,3),varname='u3',filename='flowini3d.h5',explicit=.true.)
    call h5srite(var=prs,                  varname='p', filename='flowini3d.h5',explicit=.true.)
    call h5srite(var=tmp,                  varname='t', filename='flowini3d.h5',explicit=.true.)
    !
    roav=0.d0
    uav=0.d0
    vav=0.d0
    wav=0.d0
    tav=0.d0
    pav=0.d0
    !
    do k=0,km
    do j=0,jm
    do i=0,im 
      if(nondimen) then
        prs(i,j,k)  = thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
      else
        spc(i,j,k,:)= spcinf(:)
        prs(i,j,k)  = thermal(density=rho(i,j,k),temperature=tmp(i,j,k),species=spc(i,j,k,:))
      endif
      !
      if(i==0 .or. j==0 .or. k==0) cycle
      !
      roav=roav+rho(i,j,k)
      uav=uav+vel(i,j,k,1)
      vav=vav+vel(i,j,k,2)
      wav=wav+vel(i,j,k,3)
      tav=tav+tmp(i,j,k)
      pav=pav+prs(i,j,k)
      !
    enddo
    enddo
    enddo
    !
    roav=roav/dble(im*jm*km)
    uav = uav/dble(im*jm*km)
    vav = vav/dble(im*jm*km)
    wav = wav/dble(im*jm*km)
    tav = tav/dble(im*jm*km)
    pav = pav/dble(im*jm*km)
    !
    print*, '** mean density    :',roav
    print*, '** mean velocity   :',uav,vav,wav
    print*, '** mean temperature:',tav
    print*, '** mean pressure   :',pav
    !
    rho(:,:,:)   = rho(:,:,:)  -roav
    vel(:,:,:,1) = vel(:,:,:,1)-uav
    vel(:,:,:,2) = vel(:,:,:,2)-vav
    vel(:,:,:,3) = vel(:,:,:,3)-wav
    tmp(:,:,:)   = tmp(:,:,:)  -tav
    prs(:,:,:)   = prs(:,:,:)  -pav
    !
    call h5srite(var=x(0:im,0,0,1),        varname='x', filename='datin/flowin.h5',explicit=.true.,newfile=.true.) 
    call h5srite(var=rho,                  varname='ro',filename='datin/flowin.h5',explicit=.true.)
    call h5srite(var=vel(0:im,0:jm,0:km,1),varname='u1',filename='datin/flowin.h5',explicit=.true.)
    call h5srite(var=vel(0:im,0:jm,0:km,2),varname='u2',filename='datin/flowin.h5',explicit=.true.)
    call h5srite(var=vel(0:im,0:jm,0:km,3),varname='u3',filename='datin/flowin.h5',explicit=.true.)
    call h5srite(var=prs,                  varname='p', filename='datin/flowin.h5',explicit=.true.)
    call h5srite(var=tmp,                  varname='t', filename='datin/flowin.h5',explicit=.true.)
    
    ! call tecbin('techit.plt',x(0:im,0:jm,0:km,1),'x', &
    !                          x(0:im,0:jm,0:km,2),'y', &
    !                          x(0:im,0:jm,0:km,3),'z', &
    !                        vel(0:im,0:jm,0:km,1),'u', &
    !                        vel(0:im,0:jm,0:km,2),'v', &
    !                        vel(0:im,0:jm,0:km,3),'w' )
    !
  end subroutine hitgen
  !+-------------------------------------------------------------------+
  !| The end of the subroutine hitgen.                                 |
  !+-------------------------------------------------------------------+
  !
  subroutine hitgen_parallel
    !
    use, intrinsic :: iso_c_binding
    use fftwlink
    use commvar,   only : gridfile,im,jm,km,ia,ja,ka,hm,Mach,Reynolds, &
                          tinf,roinf,spcinf,num_species,nondimen,&
                          ickmax,iomode,icamplitude,icsolenoidal,icdilatational
    use bc,        only : twall
    use readwrite, only : readgrid, readic, readinput
    use commarray, only : x,vel,rho,tmp,prs,dvel,spc
    use solver,    only : refcal
    use parallel,  only : parallelini,mpistop,mpi_sizeof,mpirank, psum, pmax, dataswap
    use geom,      only : geomcal
    use fludyna,   only : thermal
    use hdf5io
    use gridgeneration
    use tecio
    include 'fftw3-mpi.f03'
    !
    real(8) :: urms,kenergy,ufmx,roav,uav,vav,wav,tav,pav,vmax
    integer :: i,j,k,n,clock,irandom,total_m,proc_m,m
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
    integer,allocatable :: seed(:)
    real(8) :: wn1, wn2,wn3, wn12, wna, var1, var2
    real(8) :: ran1, ran2,ran3,ran4,rn3
    complex(8) :: vac1, vac2, vac3, crn1, crn2, crn4
    real(8) :: ISEA
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: u1c,u2c,u3c
    real(C_DOUBLE), pointer, dimension(:,:,:) :: u1r,u2r,u3r
    type(C_PTR) ::  backward_plan, c_u1c, c_u2c, c_u3c, c_u1r, c_u2r, c_u3r
    character(len=1) :: modeio
    modeio = 'h'
    !
    !
    call readinput
    !
    call readic
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
    !
    call mpisizedis_half_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    call refcal
    if(mpirank==0)  print*, '** refcal done!'
    !
    allocate( vel(-hm:2*im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    allocate(rho(0:(2*im),0:jm,0:km),tmp(0:(2*im),0:jm,0:km),prs(0:(2*im),0:jm,0:km))
    !
    !
    ! Generate field
    ISEA=1.d0/224.7699d0
    !
    !! random seed
    call random_seed(size=n)
    allocate(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock  +  37  *  (/ (irandom  -  1, irandom = 1, n) /)
    call random_seed(put=seed)
    deallocate(seed)
    !
    !
    !! wavenumber generation
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      ! if(jm .ne. ja)then
      !   stop "error! jm /= ja"
      ! endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j,k) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j,k) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if(j <= (ja/2+1)) then
        k2(i,j,k) = real(j-1,8)
      else if(j<=ja) then
        k2(i,j,k) = real(j-ja-1,8)
      else
        print *,"Error, no wave number possible, j must smaller than ja-1 !"
      end if
      !
      !
      if((k+k0) <= (ka/2+1)) then
        k3(i,j,k) = real(k+k0-1,8)
      else if((k+k0)<=(ka)) then
        k3(i,j,k) = real(k+k0-ka-1,8)
      else
        print *,"Error, no wave number possible, (k+k0) must smaller than ka-1 !"
      end if
      !
      !
    end do
    end do
    end do
    !
    !! complex speed allocation
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw,kmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw,kmfftw])
    c_u3c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3c, u3c, [imfftw,jmfftw,kmfftw])
    !
    c_u1r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u1r, u1r, [2*imfftw,jmfftw,kmfftw])
    c_u2r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u2r, u2r, [2*imfftw,jmfftw,kmfftw])
    c_u3r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u3r, u3r, [2*imfftw,jmfftw,kmfftw])
    !
    backward_plan = fftw_mpi_plan_dft_c2r_3d(kafftw,jafftw,iafftw, u1c,u1r, MPI_COMM_WORLD,FFTW_MEASURE)
    !
    !! half spectral generation
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      if(k1(i,j,k)==0 .and. k2(i,j,k)==0) then
        u1c(i,j,k)=0.d0
        u2c(i,j,k)=0.d0
        u3c(i,j,k)=0.d0
      else
        !
        call random_number(ran1)
        call random_number(ran2)
        call random_number(ran3)
        call random_number(ran4)
        !
        crn1=ran1*2.d0*pi*(0.d0,1.d0)
        crn2=ran2*2.d0*pi*(0.d0,1.d0)
        rn3 =ran3*2.d0*pi
        crn4=ran4*2.d0*pi*(0.d0,1.d0)
        !
        ! Calculate the modul of the wavenumber in each direction
        wn1=real(k1(i,j,k))
        wn2=real(k2(i,j,k))
        wn3=real(k3(i,j,k))
        wn12=sqrt(wn1**2+wn2**2)
        wna=sqrt(wn1**2+wn2**2+wn3**2)
        !
        var1=IniEnergDis(ISEA*2.d0,ickmax*1.d0,wna)
        var2=sqrt(var1/4.d0/pi/wna**2)
        !
        vac1=var2*cdexp(crn1)*dcos(rn3)
        vac2=var2*cdexp(crn2)*dsin(rn3)
        vac3=var2*cdexp(crn4)
        !
        if(k1(i,j,k)==0 .and. k2(i,j,k)==0) then
          u1c(i,j,k)=vac1*icsolenoidal + vac3*wn1/wna*icdilatational
          u2c(i,j,k)=vac2*icsolenoidal + vac3*wn2/wna*icdilatational
          u3c(i,j,k)=vac3*wn3/wna*icdilatational
        else
          u1c(i,j,k)=(vac1*wna*wn2+vac2*wn1*wn3)/(wna*wn12)*icsolenoidal + vac3*wn1/wna*icdilatational
          u2c(i,j,k)=(vac2*wn2*wn3-vac1*wna*wn1)/(wna*wn12)*icsolenoidal + vac3*wn2/wna*icdilatational
          u3c(i,j,k)=-vac2*wn12/wna*icsolenoidal + vac3*wn3/wna*icdilatational
        end if
        !
      endif
    enddo
    enddo
    enddo
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !
    if(mpirank==0)  print*, '** field generated!'
    call fftw_mpi_execute_dft_c2r(backward_plan,u1c,u1r)
    call fftw_mpi_execute_dft_c2r(backward_plan,u2c,u2r)
    call fftw_mpi_execute_dft_c2r(backward_plan,u3c,u3r)
    !
    if(mpirank==0)  print*,' ** project to physical space. '
    !
    im = im*2-2
    !
    do k=1,km
    do j=1,jm
    do i=1,im
        ! 
        vel(i,j,k,1)=u1r(i,j,k)
        vel(i,j,k,2)=u2r(i,j,k)
        vel(i,j,k,3)=u3r(i,j,k)
        !
    end do
    end do
    end do
    !
    vel(0,1:jm,1:km,1)=vel(im,1:jm,1:km,1)
    vel(0,1:jm,1:km,2)=vel(im,1:jm,1:km,2)
    vel(0,1:jm,1:km,3)=vel(im,1:jm,1:km,3)
    !
    vel(0:im,0,1:km,1)=vel(0:im,jm,1:km,1)
    vel(0:im,0,1:km,2)=vel(0:im,jm,1:km,2)
    vel(0:im,0,1:km,3)=vel(0:im,jm,1:km,3)
    ! !
    vel(0:im,0:jm,0:(km-1),1)=vel(0:im,0:jm,1:km,1)
    vel(0:im,0:jm,0:(km-1),2)=vel(0:im,0:jm,1:km,2)
    vel(0:im,0:jm,0:(km-1),3)=vel(0:im,0:jm,1:km,3)
    !
    call dataswap(vel)
    !
    roav=0.d0
    uav=0.d0
    vav=0.d0
    wav=0.d0
    tav=0.d0
    pav=0.d0
    vmax=0.d0
    !
    do k=0,km
    do j=0,jm
    do i=0,im 
      rho(i,j,k)  = roinf
      tmp(i,j,k)  = tinf 
      if(nondimen) then
        prs(i,j,k)  = thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
      else
        spc(i,j,k,:)= spcinf(:)
        prs(i,j,k)  = thermal(density=rho(i,j,k),temperature=tmp(i,j,k),species=spc(i,j,k,:))
      endif
      !
      vel(i,j,k,1)= icamplitude*vel(i,j,k,1)
      vel(i,j,k,2)= icamplitude*vel(i,j,k,2)
      vel(i,j,k,3)= icamplitude*vel(i,j,k,3)
      !
      if(i==0 .or. j==0 .or. k==0) cycle
      !
      roav=roav+rho(i,j,k)
      uav=uav+vel(i,j,k,1)
      vav=vav+vel(i,j,k,2)
      wav=wav+vel(i,j,k,3)
      tav=tav+tmp(i,j,k)
      pav=pav+prs(i,j,k)
      vmax = max(vmax,sqrt(vel(i,j,k,1)**2+vel(i,j,k,2)**2+vel(i,j,k,3)**2))
      !
    enddo
    enddo
    enddo
    !
    roav= psum(roav)/dble(ia*ja*ka)
    uav = psum(uav) /dble(ia*ja*ka)
    vav = psum(vav) /dble(ia*ja*ka)
    wav = psum(wav) /dble(ia*ja*ka)
    tav = psum(tav) /dble(ia*ja*ka)
    pav = psum(pav) /dble(ia*ja*ka)
    vmax= pmax(vmax)
    !
    if(mpirank==0)  then
      print*, '** mean density    :',roav
      print*, '** mean velocity   :',uav,vav,wav
      print*, '** mean temperature:',tav
      print*, '** mean pressure   :',pav
      print*, '** max velocity    :',vmax
    endif
    !
    vel(:,:,:,1) = vel(:,:,:,1)-uav
    vel(:,:,:,2) = vel(:,:,:,2)-vav
    vel(:,:,:,3) = vel(:,:,:,3)-wav
    ! 
    call mpi_barrier(mpi_comm_world,ierr)
    !
    call h5io_init(trim('datin/flowini3d.h5'),mode='write')
    !
    if((k0+km)==ka)then
      call h5write(var=rho(0:im,0:jm,0:km),  varname='ro', mode = modeio) 
      call h5write(var=vel(0:im,0:jm,0:km,1),varname='u1', mode = modeio)
      call h5write(var=vel(0:im,0:jm,0:km,2),varname='u2', mode = modeio)
      call h5write(var=vel(0:im,0:jm,0:km,3),varname='u3', mode = modeio)
      call h5write(var=prs(0:im,0:jm,0:km),  varname='p', mode = modeio)
      call h5write(var=tmp(0:im,0:jm,0:km),  varname='t', mode = modeio)
    else
      call h5write(var=rho(0:im,0:jm,0:(km-1)),  varname='ro', mode = modeio) 
      call h5write(var=vel(0:im,0:jm,0:(km-1),1),varname='u1', mode = modeio)
      call h5write(var=vel(0:im,0:jm,0:(km-1),2),varname='u2', mode = modeio)
      call h5write(var=vel(0:im,0:jm,0:(km-1),3),varname='u3', mode = modeio)
      call h5write(var=prs(0:im,0:jm,0:(km-1)),  varname='p', mode = modeio)
      call h5write(var=tmp(0:im,0:jm,0:(km-1)),  varname='t', mode = modeio)
    end if
    call h5io_end
    !
    allocate(   x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    call gridcube(2.d0*pi,2.d0*pi,2.d0*pi)
    !
    call geomcal
    !
    call div_test(vel,dvel)
    ! !
    ! call hitsta!
    !
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1c)
    call fftw_free(c_u2c)
    call fftw_free(c_u3c)
    call fftw_free(c_u1r)
    call fftw_free(c_u2r)
    call fftw_free(c_u3r)
    call mpistop
    !
  end subroutine hitgen_parallel
  !+-------------------------------------------------------------------+
  !| The end of the subroutine hitgen.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used calcualted the statistics of hit.         |
  !+-------------------------------------------------------------------+
  !| Ref: Samtaney, R., Pullin, D. I., & Kosović, B. (2001). Direct    |
  !|      numerical simulation of decaying compressible turbulence and |
  !|      shocklet statistics. Physics of Fluids, 13(5), 1415-1430.    |
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17 Nov. 2023: Created by C.S. Luo @ Beihang University            |
  !+-------------------------------------------------------------------+
  subroutine hitgen2d
    !
    use commvar,   only : gridfile,im,jm,km,ia,ja,ka,hm,Mach,Reynolds, &
                          tinf,roinf,spcinf,num_species,nondimen,&
                          ickmax,iomode,icamplitude,icsolenoidal,icdilatational
    use bc,        only : twall
    use readwrite, only : readgrid,readinput, readic
    use commarray, only : x,vel,rho,tmp,prs,dvel,spc
    use solver,    only : refcal
    use parallel,  only : mpisizedis,parapp,parallelini
    use geom,      only : geomcal
    use fludyna,   only : thermal
    use hdf5io
    use gridgeneration
    use tecio
    !
    integer :: i,j
    !
    call readinput
    call readic
    !
    call mpisizedis
    print*, '** mpisizedis done!'
    !
    call parapp
    print*, '** parapp done!'
    !
    call parallelini
    print*, '** parallelini done!'
    !
    call refcal
    print*, '** refcal done!'
    !
    allocate(   x(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3) )
    allocate( vel(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3) )
    allocate(rho(0:im,0:jm,0:0),tmp(0:im,0:jm,0:0),prs(0:im,0:jm,0:0))
    print*, '** allocation finished!'
    !
    call gridcube(2.d0*pi,2.d0*pi,0.d0)
    ! call readgrid(trim(gridfile))
    !
    ! generate field
    !call div_free_2d_gen(im,jm,ickmax,vel(0:im,0:jm,0,1), vel(0:im,0:jm,0,2))
    call solenoidal_dilatational_2d_gen(im,jm,ickmax,icsolenoidal,icdilatational,vel(0:im,0:jm,0,1),vel(0:im,0:jm,0,2))
    print*, '** field generated!'
    !
    vel(0:im,0:jm,0,3)=0.d0
    !
    do j=0,jm
    do i=0,im
      !
      rho(i,j,0)  = roinf
      tmp(i,j,0)  = tinf 
      if(nondimen) then
          prs(i,j,0)  = thermal(density=rho(i,j,0),temperature=tmp(i,j,0))
      else
          spc(i,j,0,:)= spcinf(:)
          prs(i,j,0)  = thermal(density=rho(i,j,0),temperature=tmp(i,j,0),species=spc(i,j,0,:))
      endif
      vel(i,j,0,1)= icamplitude*vel(i,j,0,1)
      vel(i,j,0,2)= icamplitude*vel(i,j,0,2)
      !
    enddo
    enddo
    !
    !
    call h5io_init(trim('datin/flowini2d.h5'),mode='write')
    !
    call h5wa2d_r8(varname='ro',var=rho(0:im,0:jm,0),  dir='k')
    !
    call h5wa2d_r8(varname='u1',var=vel(0:im,0:jm,0,1),dir='k')
    call h5wa2d_r8(varname='u2',var=vel(0:im,0:jm,0,2),dir='k')
    call h5wa2d_r8(varname='p', var=prs(0:im,0:jm,0),  dir='k')
    call h5wa2d_r8(varname='t', var=tmp(0:im,0:jm,0),  dir='k')
    call h5io_end
    !
    !
    ! test
    print*, '-- Test!'
    !
    call geomcal
    ! 
    call div_test_2d(vel,dvel)
    !
    call hitsta2d
    !
    call mpistop
    !
  end subroutine hitgen2d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine hitgen2d.                               |
  !+-------------------------------------------------------------------+
  !
  !
  subroutine hitgen2d_parallel
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput, readgrid, readic
    use fftwlink
    use commvar,   only : gridfile,im,jm,km,ia,ja,ka,hm,Mach,Reynolds, &
                          tinf,roinf,spcinf,num_species,nondimen,&
                          ickmax,iomode,icamplitude,icsolenoidal,icdilatational
    use bc,        only : twall
    use commarray, only : x,vel,rho,tmp,prs,dvel,spc
    use solver,    only : refcal
    use parallel,  only : parallelini,mpistop,mpi_sizeof,mpirank, psum, pmax, dataswap
    use geom,      only : geomcal
    use fludyna,   only : thermal
    use hdf5io
    use gridgeneration
    use tecio
    include 'fftw3-mpi.f03'
    !
    integer :: i,j,n,clock,irandom,total_m,proc_m,m
    real(8), allocatable, dimension(:,:) :: k1,k2
    integer,allocatable :: seed(:)
    real(8) :: wn1, wn2, wna, var1, var2, ran1, ran2
    real(8) :: dudi,lambda,ke0,en0,lint,tau,eta0,vmax
    complex(8) :: vac1, vac2, crn1, crn2
    real(8) :: Kenergy,Enstropy,ITGscale,LETT,KolmLength,urms,ufmx,ISEA
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1c,u2c
    real(C_DOUBLE), pointer, dimension(:,:) :: u1r,u2r
    type(C_PTR) ::  backward_plan, c_u1c, c_u2c, c_u1r, c_u2r
    !
    call readinput
    call readic
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja
    !
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    call mpisizedis_half_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    call refcal
    if(mpirank==0)  print*, '** refcal done!'
    !
    !
    allocate( vel(-hm:2*im+hm,-hm:jm+hm,-hm:hm,1:3) )
    allocate(rho(0:(2*im),0:jm,0:0),tmp(0:(2*im),0:jm,0:0),prs(0:(2*im),0:jm,0:0))
    if(mpirank==0)  print*, '** allocation finished!'
    !
    ! call readgrid(trim(gridfile))
    !
    ! Generate field
    ISEA=1.d0/224.7699d0
    !
    !! random seed
    call random_seed(size=n)
    allocate(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock  +  37  *  (/ (irandom  -  1, irandom = 1, n) /)
    call random_seed(put=seed)
    deallocate(seed)
    !
    !! wavenumber generation
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    do j=1,jm
    do i=1,im
      !
      ! if(jm .ne. ja)then
      !   stop "error! jm /= ja"
      ! endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if((j+j0) <= (ja/2+1)) then
        k2(i,j) = real(j+j0-1,8)
      else if((j+j0)<=(ja)) then
        k2(i,j) = real(j+j0-ja-1,8)
      else
        print *,"Error, no wave number possible, (j+jm) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    !
    !! complex speed allocation
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw])
    !
    c_u1r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u1r, u1r, [2*imfftw,jmfftw])
    c_u2r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u2r, u2r, [2*imfftw,jmfftw])
    !
    backward_plan = fftw_mpi_plan_dft_c2r_2d(jafftw,iafftw, u1c,u1r, MPI_COMM_WORLD,FFTW_MEASURE)
    !
    !! half spectral generation
    !
    do j=1,jm
    do i=1,im
      if(k1(i,j)==0 .and. k2(i,j)==0) then
        u1c(i,j)=0.d0
        u2c(i,j)=0.d0
      else
        call random_number(ran1)
        call random_number(ran2)
        ! ran1: random number distributied in (0,1)
        !
        crn1=ran1*2.d0*pi*(0.d0,1.d0)
        crn2=ran2*2.d0*pi*(0.d0,1.d0)
        !
        ! Calculate the modul of the wavenumber in each direction
        wn1=real(k1(i,j))
        wn2=real(k2(i,j))
        wna=sqrt(wn1**2+wn2**2)
        !
        var1=IniEnergDis(ISEA*2.d0,ickmax*1.d0,wna)
        var2=sqrt(var1/2.d0/pi/wna)
        !
        vac1=var2*cdexp(crn1)
        vac2=var2*cdexp(crn2)
        !
        u1c(i,j)=vac1*wn2/wna*icsolenoidal + vac2*wn1/wna*icdilatational
        u2c(i,j)=-vac1*wn1/wna*icsolenoidal + vac2*wn2/wna*icdilatational
      end if
    enddo
    enddo
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !
    if(mpirank==0)  print*, '** field generated!'
    !
    call fftw_mpi_execute_dft_c2r(backward_plan,u1c,u1r)
    call fftw_mpi_execute_dft_c2r(backward_plan,u2c,u2r)
    !
    if(mpirank==0)  print*,' ** project to physical space. '
    !
    im = im*2-2
    !
    do j=1,jm
    do i=1,im
      ! 
      vel(i,j,0,1)=u1r(i,j)
      vel(i,j,0,2)=u2r(i,j)
      !
    end do
    end do
    !
    vel(0,1:jm,0,1)=vel(im,1:jm,0,1)
    vel(0,1:jm,0,2)=vel(im,1:jm,0,2)
    !
    vel(0:im,0:(jm-1),0,1)=vel(0:im,1:jm,0,1)
    vel(0:im,0:(jm-1),0,2)=vel(0:im,1:jm,0,2)
    !
    ! Prepare other field
    vel(0:im,0:jm,0,3)=0.d0
    !
    call dataswap(vel)
    !
    do j=0,jm
    do i=0,im
      !
      rho(i,j,0)  = roinf
      tmp(i,j,0)  = tinf 
      if(nondimen) then
          prs(i,j,0)  = thermal(density=rho(i,j,0),temperature=tmp(i,j,0))
      else
          spc(i,j,0,:)= spcinf(:)
          prs(i,j,0)  = thermal(density=rho(i,j,0),temperature=tmp(i,j,0),species=spc(i,j,0,:))
      endif
      vel(i,j,0,1)= icamplitude*vel(i,j,0,1)
      vel(i,j,0,2)= icamplitude*vel(i,j,0,2)
      !
      vmax = max(vmax,sqrt(vel(i,j,0,1)**2+vel(i,j,0,2)**2))
      !
    enddo
    enddo
    !
    vmax= pmax(vmax)
    !
    if(mpirank==0)  then
      print*, '** max velocity    :',vmax
    endif
    !
    ! Output
    call h5io_init(trim('datin/flowini2d.h5'),mode='write')
    !
    if((j0+jm)==ja)then
      call h5wa2d_r8(varname='ro',var=rho(0:im,0:jm,0),  dir='k')
      call h5wa2d_r8(varname='u1',var=vel(0:im,0:jm,0,1),dir='k')
      call h5wa2d_r8(varname='u2',var=vel(0:im,0:jm,0,2),dir='k')
      call h5wa2d_r8(varname='p', var=prs(0:im,0:jm,0),  dir='k')
      call h5wa2d_r8(varname='t', var=tmp(0:im,0:jm,0),  dir='k')
    else
      call h5wa2d_r8(varname='ro',var=rho(0:im,0:(jm-1),0),  dir='k')
      call h5wa2d_r8(varname='u1',var=vel(0:im,0:(jm-1),0,1),dir='k')
      call h5wa2d_r8(varname='u2',var=vel(0:im,0:(jm-1),0,2),dir='k')
      call h5wa2d_r8(varname='p', var=prs(0:im,0:(jm-1),0),  dir='k')
      call h5wa2d_r8(varname='t', var=tmp(0:im,0:(jm-1),0),  dir='k')
    endif
    !
    call h5io_end
    !
    ! Param
    !
    ! if(mpirank == 0) then
    !   ke0=3.d0*ISEA/64.d0*sqrt(2.d0*pi)*dble(ickmax**5)
    !   en0=15.d0*ISEA/256.d0*sqrt(2.d0*pi)*dble(ickmax**7)
    !   lint=sqrt(2.d0*pi)/ke0
    !   tau =sqrt(32.d0/ISEA*sqrt(2.d0*pi))/sqrt(dble(ickmax**7))
    !   eta0=1.d0/sqrt(sqrt(2.d0*en0*Reynolds**2))
    !   !
    !   print*,' ---------------------------------------------------------------'
    !   print*,'        statistics according to the initial energy spectrum     '
    !   print*,' --------------------------+------------------------------------'
    !   print*,'                   kenergy |',ke0
    !   print*,'                 enstrophy |',en0
    !   print*,'           integral length |',lint
    !   print*,'  large-eddy-turnover time |',tau
    !   print*,'         kolmogorov length |',eta0
    !   print*,' --------------------------+------------------------------------'
    ! endif
    !
    ! Test
    if(mpirank==0)  print*, '-- Test!'
    ! 
    allocate(x(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3) )
    call gridcube(2.d0*pi,2.d0*pi,0.d0)
    call geomcal
    !
    call div_test_2d(vel,dvel)
    !
    !call hitsta2d
    !
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1c)
    call fftw_free(c_u2c)
    call fftw_free(c_u1r)
    call fftw_free(c_u2r)
    call mpistop
    ! 
    !
  end subroutine hitgen2d_parallel
  !+-------------------------------------------------------------------+
  !| The end of the subroutine hitgen2d_parallel.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used calcualted the statistics of hit.         |
  !+-------------------------------------------------------------------+
  subroutine hitsta
    !
    use constdef
    use commvar,   only : reynolds,mach,im,jm,km
    use commarray, only : vel,dvel,rho,tmp,x
    use fludyna,   only : miucal
    !
    real(8) :: roav,kenergy,enstrophy,dissp,miuav,urms,kolmlength,ukolm, &
               timekolm,taylorlength,retaylor,tav,machrms
    !
    integer :: i,j,k
    real(8) :: var1,vort1,vort2,vort3,vorts,s11,s22,s33,s12,s13,s23, &
               div,dudx2,miu,ufmx
    !
      roav=0.d0
      tav=0.d0
      !
      kenergy=0.d0
      !
      enstrophy=0.d0
      !
      dissp=0.d0
      !
      miuav=0.d0
      urms=0.d0
      dudx2=0.d0
      !
      ufmx=0.d0
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        roav=roav+rho(i,j,k)
        tav=tav+tmp(i,j,k)
        !
        var1=vel(i,j,k,1)**2+vel(i,j,k,2)**2+vel(i,j,k,3)**2
        kenergy=kenergy+rho(i,j,k)*var1
        !
        urms=urms+var1
        !
        vort1=dvel(i,j,k,2,3)-dvel(i,j,k,3,2)
        vort2=dvel(i,j,k,3,1)-dvel(i,j,k,1,3)
        vort3=dvel(i,j,k,1,2)-dvel(i,j,k,2,1)
        vorts=vort1*vort1+vort2*vort2+vort3*vort3
        !
        enstrophy=enstrophy+rho(i,j,k)*vorts
        !
        miu=miucal(tmp(i,j,k))/reynolds
        miuav=miuav+miu
        !
        s11=dvel(i,j,k,1,1)
        s22=dvel(i,j,k,2,2)
        s33=dvel(i,j,k,3,3)
        div=s11+s22+s33
        s12=0.5d0*(dvel(i,j,k,2,1)+dvel(i,j,k,1,2))
        s13=0.5d0*(dvel(i,j,k,3,1)+dvel(i,j,k,1,3))
        s23=0.5d0*(dvel(i,j,k,3,2)+dvel(i,j,k,2,3))
        dissp=dissp+2.d0*miu*(s11**2+s22**2+s33**2+2.d0*(s12**2+s13**2+s23**2)-num1d3*div**2)
        
        !
        dudx2=dudx2+dvel(i,j,k,1,1)**2+dvel(i,j,k,2,2)**2+dvel(i,j,k,3,3)**2
        !
        ufmx=max(ufmx,abs(vel(i,j,k,1)),abs(vel(i,j,k,2)),abs(vel(i,j,k,3)))
      end do
      end do
      end do
      !
      var1=dble(im*jm*km)

      roav      =roav     /var1   
      tav       =tav     /var1     
      kenergy   =kenergy  /var1*0.5d0
      enstrophy =enstrophy/var1*0.5d0
      dissp     =dissp    /var1
      !
      miuav     =miuav    /var1
      !
      urms      =urms /var1
      dudx2     =dudx2/var1
      !
      kolmlength=sqrt(sqrt((miuav/roav)**3/dissp))
      !
      ukolm=sqrt(sqrt(dissp*miuav/roav))
      !
      timekolm=sqrt(miuav/roav/dissp)
      !
      taylorlength=sqrt(urms/dudx2)
      retaylor=taylorlength*roav*urms/miuav/1.7320508075688773d0
      !
      machrms=urms/sqrt(tav)*mach
      !
      print*,' ---------------------------------------------------------------'
      print*,'              statistics according to actual field              '
      print*,' --------------------------+------------------------------------'
      print*,'                      urms |',urms
      print*,'                   machrms |',machrms
      print*,'                   kenergy |',kenergy
      print*,'           max fluctuation |',ufmx
      print*,'                 enstrophy |',Enstrophy
      print*,'             Kolmlength, η |',kolmlength
      print*,'                       η/Δ |',kolmlength/(x(1,0,0,1)-x(0,0,0,1))
      print*,'                     ukolm |',ukolm
      print*,'                     tkolm |',kolmlength/ukolm
      print*,'              Taylorlength |',taylorlength
      print*,'                  Retaylor |',retaylor
      print*,' --------------------------+------------------------------------'
      !
      !
  end subroutine hitsta
  !+-------------------------------------------------------------------+
  !| The end of the subroutine hitsta.                                 |
  !+-------------------------------------------------------------------+
  !
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used calcualted the statistics of hit.         |
  !+-------------------------------------------------------------------+
  subroutine hitsta2d
    !
    use constdef
    use commvar,   only : reynolds,mach,im,jm,ia,ja
    use commarray, only : vel,dvel,rho,tmp,x
    use fludyna,   only : miucal
    use parallel,  only : psum
    !
    real(8) :: roav,kenergy,enstrophy,dissp,miuav,urms,kolmlength,ukolm, &
               timekolm,taylorlength,retaylor,tav,machrms
    !
    integer :: i,j,k
    real(8) :: var1,vort3,vorts,s11,s22,s12,div,dudx2,miu,ufmx
    !
      roav=0.d0
      tav=0.d0
      !
      kenergy=0.d0
      !
      enstrophy=0.d0
      !
      dissp=0.d0
      !
      miuav=0.d0
      urms=0.d0
      dudx2=0.d0
      !
      ufmx=0.d0
      !
      do j=1,jm
      do i=1,im
        !
        roav=roav+rho(i,j,0)
        tav=tav+tmp(i,j,0)
        !
        var1=vel(i,j,0,1)**2+vel(i,j,0,2)**2
        kenergy=kenergy+rho(i,j,0)*var1
        !
        urms=urms+var1
        !
        vort3=dvel(i,j,0,1,2)-dvel(i,j,0,2,1)
        vorts=vort3*vort3
        !
        enstrophy=enstrophy+rho(i,j,0)*vorts
        !
        miu=miucal(tmp(i,j,0))/reynolds
        miuav=miuav+miu
        !
        s11=dvel(i,j,0,1,1)
        s22=dvel(i,j,0,2,2)
        div=s11+s22
        s12=0.5d0*(dvel(i,j,0,2,1)+dvel(i,j,0,1,2))
        dissp=dissp+2.d0*miu*(s11**2+s22**2+2.d0*s12**2-num1d3*div**2)
        
        !
        dudx2=dudx2+dvel(i,j,0,1,1)**2+dvel(i,j,0,2,2)**2+dvel(i,j,0,3,3)**2
        !
        ufmx=max(ufmx,abs(vel(i,j,0,1)),abs(vel(i,j,0,2)))
      end do
      end do
      !
      var1=dble(ia*ja)
      !
      roav      =psum(roav)     /var1   
      tav       =psum(tav)     /var1     
      kenergy   =psum(kenergy)  /var1*0.5d0
      enstrophy =psum(enstrophy)/var1*0.5d0
      dissp     =psum(dissp)   /var1
      !
      miuav     =psum(miuav)    /var1
      !
      urms      =sqrt(psum(urms) /var1)
      dudx2     =psum(dudx2)/var1
      !
      kolmlength=sqrt(sqrt((miuav/roav)**3/dissp))
      !
      ukolm=sqrt(sqrt(dissp*miuav/roav))
      !
      timekolm=sqrt(miuav/roav/dissp)
      !
      taylorlength=urms/sqrt(dudx2)
      retaylor=taylorlength*roav*urms/miuav/1.7320508075688773d0
      !
      machrms=urms/sqrt(tav)*mach
      !
      print*,' ---------------------------------------------------------------'
      print*,'              statistics according to actual field              '
      print*,' --------------------------+------------------------------------'
      print*,'                      urms |',urms
      print*,'                   machrms |',machrms
      print*,'                   kenergy |',kenergy
      print*,'           max fluctuation |',ufmx
      print*,'                 enstrophy |',Enstrophy
      print*,'             Kolmlength, η |',kolmlength
      print*,'                       η/Δ |',kolmlength/(x(1,0,0,1)-x(0,0,0,1))
      print*,'                     ukolm |',ukolm
      print*,'                     tkolm |',kolmlength/ukolm
      print*,'              Taylorlength |',taylorlength
      print*,'                  Retaylor |',retaylor
      print*,' --------------------------+------------------------------------'
      !
      !
  end subroutine hitsta2d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine hitsta2d.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used calcualted the divergence.                |
  !+-------------------------------------------------------------------+
  subroutine div_test(u,du)
    !
    use commvar,   only : im,jm,km,hm,ia,ja,ka
    use parallel,  only : dataswap, psum, mpirank
    use comsolver, only : solvrinit,grad
    !
    real(8),intent(inout) :: u(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3)
    !
    real(8),allocatable,intent(out),dimension(:,:,:,:,:) :: du
    !
    integer :: i,j,k
    real(8) :: div,div2
    !
    allocate( du(0:im,0:jm,0:km,1:3,1:3))
    !
    call dataswap(u)
    !
    call solvrinit

    du(:,:,:,:,1)=grad(u(:,:,:,1))
    du(:,:,:,:,2)=grad(u(:,:,:,2))
    du(:,:,:,:,3)=grad(u(:,:,:,3))
    !
    div=0.d0
    div2=0.d0
    do k=1,km
    do j=1,jm
    do i=1,im
      div =div +du(i,j,k,1,1)+du(i,j,k,2,2)+du(i,j,k,3,3)
      div2=div2+(du(i,j,k,1,1)+du(i,j,k,2,2)+du(i,j,k,3,3))**2
    enddo
    enddo
    enddo
    !
    div = psum(div) /dble(ia*ja*ka)
    div2= psum(div2)/dble(ia*ja*ka)
    !
    if(mpirank == 0) then
      print*,' ** averaged div is:',div
      print*,' ** variance div is:',div2
    endif
    !
  end subroutine div_test
  !+-------------------------------------------------------------------+
  !| The end of the subroutine div_test.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used calcualted the divergence.                |
  !+-------------------------------------------------------------------+
  subroutine div_test_2d(u,du)
    !
    use commvar,   only : im,jm,km,hm,ia,ja,ka
    use parallel,  only : dataswap, psum, mpirank
    use comsolver, only : solvrinit,grad
    !
    real(8),intent(inout) :: u(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3)
    !
    real(8),allocatable,intent(out),dimension(:,:,:,:,:) :: du
    !
    integer :: i,j,k
    real(8) :: div,div2
    !
    allocate( du(0:im,0:jm,0:0,1:3,1:3))
    !
    call dataswap(u)
    !
    call solvrinit
    !
    du(:,:,:,:,1)=grad(u(:,:,:,1))
    du(:,:,:,:,2)=grad(u(:,:,:,2))
    du(:,:,:,:,3)=0.0d0
    !
    div=0.d0
    div2=0.d0
    do j=1,jm
    do i=1,im
      div =div +du(i,j,0,1,1)+du(i,j,0,2,2)
      div2=div2+(du(i,j,0,1,1)+du(i,j,0,2,2))**2
    enddo
    enddo
    !
    div = psum(div) /dble(ia*ja)
    div2= psum(div2)/dble(ia*ja)
    !
    if(mpirank == 0) then
      print*,' ** averaged div is:',div
      print*,' ** variance div is:',div2
    endif
    !
  end subroutine div_test_2d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine div_test_2d.                            |
  !+-------------------------------------------------------------------+
  !
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate a divergence free fluctuation.|
  !+-------------------------------------------------------------------+
  !| Ref: Blaisdell, G. A., Numerical simulation of compressible       |
  !|      homogeneous turbulence, Phd, 1991, Stanford University       |
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 26-09-2022: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine div_free_gen(idim,jdim,kdim,u1,u2,u3)
    !
    use singleton
    use commvar,only : Reynolds
    !
    integer,intent(in) :: idim,jdim,kdim
    real(8),intent(out),dimension(0:idim,0:jdim,0:kdim) :: u1,u2,u3
    !
    ! local data
    integer :: kmi,kmj,kmk,kmax
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! kmi: maximal wavenumber in i direction
    ! kmj: maximal wavenumber in j direction
    ! kmk: maximal wavenumber in k direction
    ! kmax: maximal wavenumber in all direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) :: wn1,wn2,wn3,wn12,wna
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! wn1: modul of wavenumber in i direction
    ! wn2: modul of wavenumber in j direction
    ! wn3: modul of wavenumber in k direction
    ! wn12: wn12=sqrt(wn1**2+wn2**2)
    ! wna: modul of wavenumber in all direction
    ! (k0*1.d0): the wavenumber at maximum given 
    !     spectrum
    ! Ac: the intensity of given spectrum
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, dimension(:,:,:) :: u1tp,u2tp,u3tp
    !
    complex(8), allocatable, dimension(:,:,:) :: u1c,u2c,u3c,u1ct,u2ct,u3ct,u4ct
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Egv: the given initial energy spectrum
    ! u1c: the spectral velocity in k1 direction
    ! u2c: the spectral velocity in k2 direction
    ! u3c: the spectral velocity in k3 direction
    ! uct: the spectrl variable in (1~*2km)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) :: Kenergy,Enstropy,ITGscale,LETT,KolmLength,urms,ufmx
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Kenergy: initial volume averaged turbulent 
    !          kinetic energy
    ! Enstropy: initial volume averaged e
    !           nstrophy
    ! ITGscale: initial integral length scale
    ! LETT: initial large-eddy-turnover time
    ! KolmLength: initial Kolmogorov scale
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: k1,k2,k3,k0,i,j,k
    real(8) :: ran1,ran2,ran3,rn1,rn2,rn3,var1,var2,var3,ISEA
    complex(8) :: vac1,vac2,vac3,vac4,crn1,crn2
    real(8) :: dudi,lambda,ke0,en0,lint,tau,eta0
    !
    kmi=idim/2
    kmj=jdim/2
    kmk=kdim/2
    kmax=idnint(sqrt((kmi**2+kmj**2+kmk**2)*1.d0))+1
    !
    allocate(u1c(-kmi:kmi,-kmj:kmj,-kmk:kmk),                        &
             u2c(-kmi:kmi,-kmj:kmj,-kmk:kmk),                        &
             u3c(-kmi:kmi,-kmj:kmj,-kmk:kmk)                         )
    allocate(u1ct(1:idim,1:jdim,1:kdim),u2ct(1:idim,1:jdim,1:kdim),              &  
             u3ct(1:idim,1:jdim,1:kdim) )
    allocate(u1tp(1:idim,1:jdim,1:kdim),u2tp(1:idim,1:jdim,1:kdim),              &
             u3tp(1:idim,1:jdim,1:kdim),u4ct(1:idim,1:jdim,1:kdim)               )
    !
    ! Give the inital energy spectrum.
    ISEA=1.d0/224.7699d0
    k0=4
    !
    ! Generate the random velocity field according the given energy 
    ! spectrum
    ! Blaisdell, G. A. 1991 took E(k)=Integer(ui*uicoj*dA(k)). 
    ! This program takes E(k)=Integer(0.5*ui*uicoj*dA(k)). 
    ! Therefor, we take the Ek as twice of that from Blaisdell.
    print*,' ** Generate the random velocity field according the given energy spectrum'
    do k1=0,kmi
    do k2=-kmj,kmj
    do k3=-kmk,kmk
      !
      call random_number(ran1)
      call random_number(ran2)
      call random_number(ran3)
      ! ran1,ran2,ran3: random number distributied in (0,1)
      !
      rn1=ran1*2.d0*pi
      rn2=ran2*2.d0*pi
      rn3=ran3*2.d0*pi
      !
      ! Calculate the modul of the wavenumber in each direction
      wn1=real(k1,8)
      wn2=real(k2,8)
      wn3=real(k3,8)
      wn12=sqrt(wn1**2+wn2**2)
      wna=sqrt(wn1**2+wn2**2+wn3**2)
      !
      ! Calculate the initidiml energy spectral
      if(k1==0 .and. k2==0 .and. k3==0) then
        var1=0.d0
        var2=0.d0
      else
        var1=IniEnergDis(ISEA*2.d0,K0*1.d0,wna)
        var2=sqrt(var1/4.d0/pi/wna**2)
        ! var2=1.d0
      end if
      !
      ! Gererate the velocity spectrum in half-wavenumber space.
      crn1=rn1*(0.d0,1.d0)
      crn2=rn2*(0.d0,1.d0)
      !
      vac1=var2*cdexp(crn1)*dcos(rn3)
      vac2=var2*cdexp(crn2)*dsin(rn3)
      !
      if(k1==0 .and. k2==0 .and. k3==0) then
        u1c(k1,k2,k3)=0.d0
        u2c(k1,k2,k3)=0.d0
        u3c(k1,k2,k3)=0.d0
      elseif(k1==0 .and. k2==0) then
        u1c(k1,k2,k3)=vac1
        u2c(k1,k2,k3)=vac2
        u3c(k1,k2,k3)=0.d0
      else
        u1c(k1,k2,k3)=(vac1*wna*wn2+vac2*wn1*wn3)/(wna*wn12)
        u2c(k1,k2,k3)=(vac2*wn2*wn3-vac1*wna*wn1)/(wna*wn12)
        u3c(k1,k2,k3)=-vac2*wn12/wna
      end if
      !
    end do
    end do
    end do
    !
    print*,' ** Generate the velocity spectrum in another half-wavenumber space '
    ! Generate the velocity spectrum in another half-wavenumber space
    ! by using conjunction relation
    do k1=-kmi,-1
    do k2=-kmj,kmj
    do k3=-kmk,kmk
      u1c(k1,k2,k3)=conjg(u1c(-k1,-k2,-k3))
      u2c(k1,k2,k3)=conjg(u2c(-k1,-k2,-k3))
      u3c(k1,k2,k3)=conjg(u3c(-k1,-k2,-k3))
    end do
    end do
    end do
    ! !
    ! Transform the spectrum from (-N/2+1,N/2) to (1,N) fo rthe 
    ! convenience of using external FFT subroutine
    print*,' ** Transform the spectrum from (-N/2+1,N/2) to (1,N)  '
    !
    do k=1,kdim
    do j=1,jdim
    do i=1,idim
      if(i<=idim/2+1) then
        k1=i-1
      else
        k1=i-idim-1
      end if
      if(j<=jdim/2+1) then
        k2=j-1
      else
        k2=j-jdim-1
      end if
      if(k<=kdim/2+1) then
        k3=k-1
      else
        k3=k-kdim-1
      end if
      !
      u1ct(i,j,k)=u1c(k1,k2,k3)
      u2ct(i,j,k)=u2c(k1,k2,k3)
      u3ct(i,j,k)=u3c(k1,k2,k3)
    end do
    end do
    end do
    ! !
    u1ct=FFT(u1ct,inv=.true.)
    u2ct=FFT(u2ct,inv=.true.)
    u3ct=FFT(u3ct,inv=.true.)
    !
    print*,' ** project to physical space. '
    !
    do k=1,kdim
    do j=1,jdim
    do i=1,idim
      ! multiply sqrt(NxNyNz) for return standard FFT
      u1(i,j,k)=real(u1ct(i,j,k),8)*sqrt(real(idim*jdim*kdim,8))
      u2(i,j,k)=real(u2ct(i,j,k),8)*sqrt(real(idim*jdim*kdim,8))
      u3(i,j,k)=real(u3ct(i,j,k),8)*sqrt(real(idim*jdim*kdim,8))
      !
    end do
    end do
    end do
    !
    u1(0,1:jdim,1:kdim)=u1(idim,1:jdim,1:kdim)
    u2(0,1:jdim,1:kdim)=u2(idim,1:jdim,1:kdim)
    u3(0,1:jdim,1:kdim)=u3(idim,1:jdim,1:kdim)
    !
    u1(0:idim,0,1:kdim)=u1(0:idim,jdim,1:kdim)
    u2(0:idim,0,1:kdim)=u2(0:idim,jdim,1:kdim)
    u3(0:idim,0,1:kdim)=u3(0:idim,jdim,1:kdim)
    !
    u1(0:idim,0:jdim,0)=u1(0:idim,0:jdim,kdim)
    u2(0:idim,0:jdim,0)=u2(0:idim,0:jdim,kdim)
    u3(0:idim,0:jdim,0)=u3(0:idim,0:jdim,kdim)
    !
    ! urms=0.d0
    ! ! ufmx=0.d0
    ! ! Kenergy=0.d0
    ! do k=1,kdim
    ! do j=1,jdim
    ! do i=1,idim
    !   Kenergy=Kenergy+0.5d0*(u1(i,j,k)**2+u2(i,j,k)**2+u3(i,j,k)**2)
    !   urms=urms+u1(i,j,k)**2+u2(i,j,k)**2+u3(i,j,k)**2
    !   ufmx=max(ufmx,dabs(u1(i,j,k)),dabs(u2(i,j,k)),dabs(u3(i,j,k)))
    ! end do
    ! end do
    ! end do
    ! urms=sqrt(urms/real(idim*jdim*kdim,8))
    ! Kenergy=Kenergy/real(idim*jdim*kdim,8)
    ! !
    ! u1=u1/urms
    ! u2=u2/urms
    ! u3=u3/urms
    ! Kenergy=Kenergy/urms/urms
    ! urms=urms/urms
    ! !
    ke0=3.d0*ISEA/64.d0*sqrt(2.d0*pi)*dble(k0**5)
    en0=15.d0*ISEA/256.d0*sqrt(2.d0*pi)*dble(k0**7)
    lint=sqrt(2.d0*pi)/ke0
    tau =sqrt(32.d0/ISEA*sqrt(2.d0*pi))/sqrt(dble(k0**7))
    eta0=1.d0/sqrt(sqrt(2.d0*en0*Reynolds**2))
    !
    print*,' ---------------------------------------------------------------'
    print*,'        statistics according to the initial energy spectrum     '
    print*,' --------------------------+------------------------------------'
    print*,'                   kenergy |',ke0
    print*,'                 enstrophy |',en0
    print*,'           integral length |',lint
    print*,'  large-eddy-turnover time |',tau
    print*,'         kolmogorov length |',eta0
    print*,' --------------------------+------------------------------------'
    ! !
    ! call h5srite(var=u1,varname='u1',filename='velocity.h5',explicit=.true.,newfile=.true.)
    ! call h5srite(var=u2,varname='u2',filename='velocity.h5',explicit=.true.)
    ! call h5srite(var=u3,varname='u3',filename='velocity.h5',explicit=.true.)
    !
    deallocate(u1c,u2c,u3c)
    deallocate(u1ct,u2ct,u3ct)
    deallocate(u1tp,u2tp,u3tp)
    !
  end subroutine div_free_gen
  !+-------------------------------------------------------------------+
  !| The end of the function div_free_gen.                             |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate a divergence free fluctuation |
  !| in 2D.                                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16 Oct. 2023: Created by C.S. Luo @ Beihang University            |
  !+-------------------------------------------------------------------+
  subroutine div_free_2d_gen(idim,jdim,k0,u1,u2)
    !
    use singleton
    use commvar,only : Reynolds
    !
    integer,intent(in) :: idim,jdim,k0
    real(8),intent(out),dimension(0:idim,0:jdim,0:0) :: u1,u2
    !
    ! local data
    integer :: kmi,kmj,kmk,kmax
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! kmi: maximal wavenumber in i direction
    ! kmj: maximal wavenumber in j direction
    ! kmk: maximal wavenumber in k direction
    ! kmax: maximal wavenumber in all direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) :: wn1,wn2,wna, Ac
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! wn1: modul of wavenumber in i direction
    ! wn2: modul of wavenumber in j direction
    ! wna: modul of wavenumber in all direction
    ! (k0*1.d0): the wavenumber at maximum given 
    !     spectrum
    ! Ac: the intensity of given spectrum
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    complex(8), allocatable, dimension(:,:) :: u1c,u2c,u1ct,u2ct
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Egv: the given initial energy spectrum
    ! u1c: the spectral velocity in k1 direction
    ! u2c: the spectral velocity in k2 direction
    ! u3c: the spectral velocity in k3 direction
    ! uct: the spectrl variable in (1~*2km)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) :: Kenergy,Enstropy,ITGscale,LETT,KolmLength,urms,ufmx,ISEA
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Kenergy: initial volume averaged turbulent 
    !          kinetic energy
    ! Enstropy: initial volume averaged e
    !           nstrophy
    ! ITGscale: initial integral length scale
    ! LETT: initial large-eddy-turnover time
    ! KolmLength: initial Kolmogorov scale
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: k1,k2,i,j,n,clock, irandom
    !
    real(8) :: ran1,rn1,var1,var2
    complex(8) :: crn1,vac1
    real(8) :: ke0,en0,lint,tau,eta0
    integer,allocatable :: seed(:)
    !
    kmi=idim/2
    kmj=jdim/2
    kmax=idnint(sqrt((kmi**2+kmj**2)*1.d0))+1
    !
    allocate(u1c(-kmi:kmi,-kmj:kmj),u2c(-kmi:kmi,-kmj:kmj))
    allocate(u1ct(1:idim,1:jdim),u2ct(1:idim,1:jdim))
    !
    ! Give the inital energy spectrum.
    ISEA=1.d0/224.7699d0
    !
    ! Generate the random velocity field according the given energy spectrum
    ! Blaisdell, G. A. 1991 took E(k)=Integer(ui*uicoj*dA(k)). 
    ! This program takes E(k)=Integer(0.5*ui*uicoj*dA(k)). 
    ! Therefore, we take the Ek as twice of that from Blaisdell.
    print*,' ** Generate the random velocity field according the given energy spectrum'
    !
    ! Insert random seed
    call random_seed(size=n)
    allocate(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock  +  37  *  (/ (irandom  -  1, irandom = 1, n) /)
    call random_seed(put=seed)
    deallocate(seed)
    !
    do k1=0,kmi
    do k2=-kmj,kmj
      !
      !
      call random_number(ran1)
      ! ran1: random number distributied in (0,1)
      !
      rn1=ran1*2.d0*pi
      !
      ! Calculate the modul of the wavenumber in each direction
      wn1=real(k1,8)
      wn2=real(k2,8)
      wna=sqrt(wn1**2+wn2**2)
      !
      ! Calculate the initidiml energy spectral
      if(k1==0 .and. k2==0) then
        var1=0.d0
        var2=0.d0
      else
        var1=IniEnergDis(ISEA*2.d0,K0*1.d0,wna)
        var2=sqrt(var1/4.d0/pi/wna**2)
        ! var2=1.d0
      end if
      !
      ! Gererate the velocity spectrum in half-wavenumber space.
      !
      crn1=rn1*(0.d0,1.d0)
      !
      vac1=var2*cdexp(crn1)
      !
      if(k1==0 .and. k2==0) then
        u1c(k1,k2)=0.d0
        u2c(k1,k2)=0.d0
      else
        u1c(k1,k2)=vac1*wn2/wna
        u2c(k1,k2)=-vac1*wn1/wna
      end if
      !
    end do
    end do
    !
    print*,' ** Generate the velocity spectrum in another half-wavenumber space '
    ! Generate the velocity spectrum in another half-wavenumber space
    ! by using conjunction relation
    do k1=-kmi,-1
    do k2=-kmj,kmj
      u1c(k1,k2)=conjg(u1c(-k1,-k2))
      u2c(k1,k2)=conjg(u2c(-k1,-k2))
    end do
    end do
    ! !
    ! Transform the spectrum from (-N/2+1,N/2) to (1,N) for the 
    ! convenience of using external FFT subroutine
    print*,' ** Transform the spectrum from (-N/2+1,N/2) to (1,N)  '
    !
    do j=1,jdim
    do i=1,idim
      if(i<=idim/2+1) then
        k1=i-1
      else
        k1=i-idim-1
      end if
      if(j<=jdim/2+1) then
        k2=j-1
      else
        k2=j-jdim-1
      end if
      !
      u1ct(i,j)=u1c(k1,k2)
      u2ct(i,j)=u2c(k1,k2)
    end do
    end do
    ! !
    u1ct=FFT(u1ct,inv=.true.)
    u2ct=FFT(u2ct,inv=.true.)
    !
    print*,' ** project to physical space. '
    !
    do j=1,jdim
    do i=1,idim
      ! multiply sqrt(NxNy) for return standard FFT
      u1(i,j,0)=real(u1ct(i,j),8)*sqrt(real(idim*jdim,8))
      u2(i,j,0)=real(u2ct(i,j),8)*sqrt(real(idim*jdim,8))
      !
    end do
    end do
    !
    u1(0,1:jdim,0)=u1(idim,1:jdim,0)
    u2(0,1:jdim,0)=u2(idim,1:jdim,0)
    !
    u1(0:idim,0,0)=u1(0:idim,jdim,0)
    u2(0:idim,0,0)=u2(0:idim,jdim,0)
    !
    ! !
    ke0=3.d0*ISEA/64.d0*sqrt(2.d0*pi)*dble(k0**5)
    en0=15.d0*ISEA/256.d0*sqrt(2.d0*pi)*dble(k0**7)
    lint=sqrt(2.d0*pi)/ke0
    tau =sqrt(32.d0/ISEA*sqrt(2.d0*pi))/sqrt(dble(k0**7))
    eta0=1.d0/sqrt(sqrt(2.d0*en0*Reynolds**2))
    !
    print*,' ---------------------------------------------------------------'
    print*,'        statistics according to the initial energy spectrum     '
    print*,' --------------------------+------------------------------------'
    print*,'                   kenergy |',ke0
    print*,'                 enstrophy |',en0
    print*,'           integral length |',lint
    print*,'  large-eddy-turnover time |',tau
    print*,'         kolmogorov length |',eta0
    print*,' --------------------------+------------------------------------'
    ! !
    !
    deallocate(u1c,u2c)
    deallocate(u1ct,u2ct)
    !
  end subroutine div_free_2d_gen
  !+-------------------------------------------------------------------+
  !| The end of the function div_free_2d_gen.                             |
  !+-------------------------------------------------------------------+
  !!
    !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate a solenoidal-dilatational     !
  !| fluctuation in 2D.                                                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17 Nov. 2023: Created by C.S. Luo @ Beihang University            |
  !+-------------------------------------------------------------------+
  subroutine solenoidal_dilatational_2d_gen(idim,jdim,k0,icsolenoidal,icdilatational,u1,u2)
    !
    use singleton
    use commvar,only : Reynolds
    !
    integer,intent(in) :: idim,jdim,k0
    real(8),intent(in) :: icsolenoidal,icdilatational
    real(8),intent(out),dimension(0:idim,0:jdim,0:0) :: u1,u2
    !
    ! local data
    integer :: kmi,kmj,kmk,kmax
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! kmi: maximal wavenumber in i direction
    ! kmj: maximal wavenumber in j direction
    ! kmk: maximal wavenumber in k direction
    ! kmax: maximal wavenumber in all direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) :: wn1,wn2,wna, Ac
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! wn1: modul of wavenumber in i direction
    ! wn2: modul of wavenumber in j direction
    ! wna: modul of wavenumber in all direction
    ! (k0*1.d0): the wavenumber at maximum given 
    !     spectrum
    ! Ac: the intensity of given spectrum
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    complex(8), allocatable, dimension(:,:) :: u1c,u2c,u1ct,u2ct
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Egv: the given initial energy spectrum
    ! u1c: the spectral velocity in k1 direction
    ! u2c: the spectral velocity in k2 direction
    ! u3c: the spectral velocity in k3 direction
    ! uct: the spectrl variable in (1~*2km)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) :: Kenergy,Enstropy,ITGscale,LETT,KolmLength,urms,ufmx,ISEA
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Kenergy: initial volume averaged turbulent 
    !          kinetic energy
    ! Enstropy: initial volume averaged e
    !           nstrophy
    ! ITGscale: initial integral length scale
    ! LETT: initial large-eddy-turnover time
    ! KolmLength: initial Kolmogorov scale
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: k1,k2,i,j,n,clock, irandom
    !
    real(8) :: ran1,ran2,rn1,rn2,var1,var2
    complex(8) :: crn1,crn2,vac1,vac2
    real(8) :: ke0,en0,lint,tau,eta0
    integer,allocatable :: seed(:)
    !
    kmi=idim/2
    kmj=jdim/2
    kmax=idnint(sqrt((kmi**2+kmj**2)*1.d0))+1
    !
    allocate(u1c(-kmi:kmi,-kmj:kmj),u2c(-kmi:kmi,-kmj:kmj))
    allocate(u1ct(1:idim,1:jdim),u2ct(1:idim,1:jdim))
    !
    ! Give the inital energy spectrum.
    ISEA=1.d0/224.7699d0
    !
    ! Generate the random velocity field according the given energy spectrum
    ! Blaisdell, G. A. 1991 took E(k)=Integer(ui*uicoj*dA(k)). 
    ! This program takes E(k)=Integer(0.5*ui*uicoj*dA(k)). 
    ! Therefore, we take the Ek as twice of that from Blaisdell.
    print*,' ** Generate the random velocity field according the given energy spectrum'
    !
    ! Insert random seed
    call random_seed(size=n)
    allocate(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock  +  37  *  (/ (irandom  -  1, irandom = 1, n) /)
    call random_seed(put=seed)
    deallocate(seed)
    !
    do k1=0,kmi
    do k2=-kmj,kmj
      !
      !
      call random_number(ran1)
      call random_number(ran2)
      ! ran1: random number distributied in (0,1)
      !
      rn1=ran1*2.d0*pi
      rn2=ran2*2.d0*pi
      !
      ! Calculate the modul of the wavenumber in each direction
      wn1=real(k1,8)
      wn2=real(k2,8)
      wna=sqrt(wn1**2+wn2**2)
      !
      ! Calculate the initidiml energy spectral
      if(k1==0 .and. k2==0) then
        var1=0.d0
        var2=0.d0
      else
        var1=IniEnergDis(ISEA*2.d0,K0*1.d0,wna)
        var2=sqrt(var1/4.d0/pi/wna**2)
        ! var2=1.d0
      end if
      !
      ! Gererate the velocity spectrum in half-wavenumber space.
      !
      crn1=rn1*(0.d0,1.d0)
      crn2=rn2*(0.d0,1.d0)
      !
      vac1=var2*cdexp(crn1)
      vac2=var2*cdexp(crn2)
      !
      if(k1==0 .and. k2==0) then
        u1c(k1,k2)=0.d0
        u2c(k1,k2)=0.d0
      else
        u1c(k1,k2)=vac1*wn2/wna*icsolenoidal + vac2*wn1/wna*icdilatational
        u2c(k1,k2)=-vac1*wn1/wna*icsolenoidal + vac2*wn2/wna*icdilatational
      end if
      !
    end do
    end do
    !
    print*,' ** Generate the velocity spectrum in another half-wavenumber space '
    ! Generate the velocity spectrum in another half-wavenumber space
    ! by using conjunction relation
    do k1=-kmi,-1
    do k2=-kmj,kmj
      u1c(k1,k2)=conjg(u1c(-k1,-k2))
      u2c(k1,k2)=conjg(u2c(-k1,-k2))
    end do
    end do
    ! !
    ! Transform the spectrum from (-N/2+1,N/2) to (1,N) for the 
    ! convenience of using external FFT subroutine
    print*,' ** Transform the spectrum from (-N/2+1,N/2) to (1,N)  '
    !
    do j=1,jdim
    do i=1,idim
      if(i<=idim/2+1) then
        k1=i-1
      else
        k1=i-idim-1
      end if
      if(j<=jdim/2+1) then
        k2=j-1
      else
        k2=j-jdim-1
      end if
      !
      u1ct(i,j)=u1c(k1,k2)
      u2ct(i,j)=u2c(k1,k2)
    end do
    end do
    ! !
    u1ct=FFT(u1ct,inv=.true.)
    u2ct=FFT(u2ct,inv=.true.)
    !
    print*,' ** project to physical space. '
    !
    do j=1,jdim
    do i=1,idim
      ! multiply sqrt(NxNy) for return standard FFT
      u1(i,j,0)=real(u1ct(i,j),8)*sqrt(real(idim*jdim,8))
      u2(i,j,0)=real(u2ct(i,j),8)*sqrt(real(idim*jdim,8))
      !
    end do
    end do
    !
    u1(0,1:jdim,0)=u1(idim,1:jdim,0)
    u2(0,1:jdim,0)=u2(idim,1:jdim,0)
    !
    u1(0:idim,0,0)=u1(0:idim,jdim,0)
    u2(0:idim,0,0)=u2(0:idim,jdim,0)
    !
    ! !
    ke0=3.d0*ISEA/64.d0*sqrt(2.d0*pi)*dble(k0**5)
    en0=15.d0*ISEA/256.d0*sqrt(2.d0*pi)*dble(k0**7)
    lint=sqrt(2.d0*pi)/ke0
    tau =sqrt(32.d0/ISEA*sqrt(2.d0*pi))/sqrt(dble(k0**7))
    eta0=1.d0/sqrt(sqrt(2.d0*en0*Reynolds**2))
    !
    print*,' ---------------------------------------------------------------'
    print*,'        statistics according to the initial energy spectrum     '
    print*,' --------------------------+------------------------------------'
    print*,'                   kenergy |',ke0
    print*,'                 enstrophy |',en0
    print*,'           integral length |',lint
    print*,'  large-eddy-turnover time |',tau
    print*,'         kolmogorov length |',eta0
    print*,' --------------------------+------------------------------------'
    ! !
    !
    deallocate(u1c,u2c)
    deallocate(u1ct,u2ct)
    !
  end subroutine solenoidal_dilatational_2d_gen
  !+-------------------------------------------------------------------+
  !| The end of the function solenoidal_dilatational_2d_gen.           |
  !+-------------------------------------------------------------------+
  !!
  !
  !
  subroutine instantspectra2Davg(thefilenumb)
    !
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km,ia,ja,ka
    use commarray, only: vel, rho, prs
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    integer,intent(in) :: thefilenumb
    character(len=128) :: infilename
    character(len=4) :: stepname
    real(8) :: u1mean,u2mean,rhomean,prsmean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1spe,u2spe,pspe
    real(8), allocatable, dimension(:) :: ES,EC,Ecount,Eall,Puc
    complex(8), allocatable, dimension(:,:) :: usspe,ucspe
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1c,u2c,u1s,u2s
    ! real(8), allocatable, dimension(:,:) :: usspeR,usspeI,ucspeR,ucspeI
    real(8), allocatable, dimension(:,:) :: u1cR,u2cR,u1sR,u2sR,ucuc,ucus,usus,ReUsConjUd
    real(8), allocatable, dimension(:,:) :: k1,k2
    integer :: allkmax
    real(8) :: k,dk !wave number
    real(8) :: Ecspe,Esspe,Pucspe,Ecphy,Esphy,ucusphy,roav,Ecmax
    real(8) :: ReUsConjUdmax,ReUsConjUdmin,k2Es,k2Ec,IntLengthAbove
    character(len=128) :: outfilename
    integer :: hand_a,hand_b
    character(len=1) :: modeio
    integer :: i,j,n
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_pspe, c_u1c, c_u2c, c_u1s, c_u2s
    !
    call readinput
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    allkmax=ceiling(sqrt(2.d0)/3*min(ia,ja))
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,"allkmax:",allkmax
    if(ka .ne. 0) stop 'Please use instantspectra3D'
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !
    allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km), prs(0:im,0:jm,0:km))
    !
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='p',  var=prs(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    ! Calculate average
    u1mean = 0.0d0
    u2mean = 0.0d0
    rhomean = 0.0d0
    prsmean = 0.0d0
    !
    do i=1,im
      do j=1,jm
        u1mean = u1mean + vel(i,j,0,1)
        u2mean = u2mean + vel(i,j,0,2)
        rhomean = rhomean + rho(i,j,0)
        prsmean = prsmean + prs(i,j,0)
      enddo
    enddo
    rhomean = psum(rhomean) / (1.0d0*ia*ja)
    u1mean = psum(u1mean) / (1.d0*ia*ja)
    u2mean = psum(u2mean) / (1.d0*ia*ja)
    prsmean = psum(prsmean) / (1.d0*ia*ja)
    if(mpirank==0) print *, 'u1mean=',u1mean, 'u2mean=',u2mean, 'prsmean=',prsmean
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw])
    c_pspe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_pspe, pspe, [imfftw,jmfftw])
    !
    ! Planning
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=CMPLX(vel(i,j,0,1)-u1mean,0.d0,C_INTPTR_T);
      u2spe(i,j)=CMPLX(vel(i,j,0,2)-u2mean,0.d0,C_INTPTR_T);
      pspe(i,j)=CMPLX(prs(i,j,0)-prsmean,0.d0,C_INTPTR_T);
      !
    end do
    end do
    !
    !!!! Do 2d FFT
    !
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    call fftw_mpi_execute_dft(forward_plan,pspe,pspe)

    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=u1spe(i,j)/(1.d0*ia*ja)
      u2spe(i,j)=u2spe(i,j)/(1.d0*ia*ja)
      pspe(i,j)=pspe(i,j)/(1.d0*ia*ja)
      !
    end do
    end do
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    do j=1,jm
    do i=1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if((j+j0) <= (ja/2+1)) then
        k2(i,j) = real(j+j0-1,8)
      else if((j+j0)<=(ja)) then
        k2(i,j) = real(j+j0-ja-1,8)
      else
        print *,"Error, no wave number possible, (j+jm) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    !
    !!!! Do S-C decomposition
    allocate(usspe(1:im,1:jm),ucspe(1:im,1:jm))
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw])
    c_u1s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1s, u1s, [imfftw,jmfftw])
    c_u2s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2s, u2s, [imfftw,jmfftw])
    allocate(u1cR(1:im,1:jm),u2cR(1:im,1:jm),u1sR(1:im,1:jm),u2sR(1:im,1:jm),&
            ucuc(1:im,1:jm),ucus(1:im,1:jm),usus(1:im,1:jm))
    allocate(ReUsConjUd(1:im,1:jm))
    !
    ! 
    !
    do j=1,jm
    do i=1,im
        k=dsqrt(k1(i,j)**2+k2(i,j)**2+1.d-15)
        usspe(i,j) = u1spe(i,j)*k2(i,j)/k - u2spe(i,j)*k1(i,j)/k
        ucspe(i,j) = u1spe(i,j)*k1(i,j)/k + u2spe(i,j)*k2(i,j)/k
        !
        u1c(i,j)=  ucspe(i,j)*k1(i,j)/k 
        u2c(i,j)=  ucspe(i,j)*k2(i,j)/k
        u1s(i,j)=  usspe(i,j)*k2(i,j)/k 
        u2s(i,j)= -usspe(i,j)*k1(i,j)/k
        !
        ReUsConjUd(i,j) = real(conjg(usspe(i,j))*ucspe(i,j))
        !
      end do
    end do
    !
    if(mpirank==0)  print*, '** spectral decomposition calculation finish'
    !
    !
    !!!! Give S-C spectra and spectral energy
    !
    dk = 0.5d0
    allocate(ES(0:allkmax),EC(0:allkmax),Ecount(0:allkmax),Eall(0:allkmax),Puc(0:allkmax))
    do i=0,allkmax
      ES(i) = 0.0d0
      EC(i) = 0.0d0
      Ecount(i) = 0.0d0
      Eall(i) = 0.0d0
      Puc(i) = 0.0d0
    end do
    !
    Ecspe = 0.0d0
    Esspe = 0.0d0
    Pucspe = 0.0d0
    k2Es = 0.d0
    k2Ec = 0.d0
    !
    do j=1,jm
    do i=1,im
        k=dsqrt(k1(i,j)**2+k2(i,j)**2+1.d-15)
        if (abs(k-nint(k))<=dk .and. kint(k,dk)<=allkmax) then
          ES(kint(k,dk)) = ES(kint(k,dk)) + usspe(i,j)*conjg(usspe(i,j))*k
          EC(kint(k,dk)) = EC(kint(k,dk)) + ucspe(i,j)*conjg(ucspe(i,j))*k
          Eall(kint(k,dk)) = Eall(kint(k,dk)) + u1spe(i,j)*conjg(u1spe(i,j))*k + &
                              u2spe(i,j)*conjg(u2spe(i,j))*k
          Puc(kint(k,dk)) = Puc(kint(k,dk)) + dimag(pspe(i,j)*dconjg(ucspe(i,j))*k)*k
          Ecount(kint(k,dk)) = Ecount(kint(k,dk)) + 1
        end if
        Ecspe = Ecspe + (ucspe(i,j)*dconjg(ucspe(i,j)))/2
        Esspe = Esspe + (usspe(i,j)*dconjg(usspe(i,j)))/2
        IntLengthAbove = IntLengthAbove + usspe(i,j)*conjg(usspe(i,j))/2/k + &
                        (ucspe(i,j)*dconjg(ucspe(i,j)))/2/k
        Pucspe = Pucspe + dimag(pspe(i,j)*dconjg(ucspe(i,j))*k)/2
        k2Es = k2Es + k**2 * (usspe(i,j)*dconjg(usspe(i,j)))/2
        k2Ec = k2Ec + k**2 * (ucspe(i,j)*dconjg(ucspe(i,j)))/2
      end do
    end do
    !
    do i=0,allkmax
      ES(i) = psum(ES(i))
      EC(i) = psum(EC(i))
      Eall(i) = psum(Eall(i))
      Ecount(i) = psum(Ecount(i))
      Puc(i) = psum(Puc(i))
    end do
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    do i=0,allkmax
      if(Ecount(i) .ne. 0)then
        ES(i) = ES(i)/Ecount(i)*pi
        EC(i) = EC(i)/Ecount(i)*pi
        Eall(i) = Eall(i)/Ecount(i)*pi
        Ecount(i) = Ecount(i)
        Puc(i) = - 2*Puc(i)/Ecount(i)*pi
      endif
    end do
    !
    !
    Ecspe = psum(Ecspe)
    Esspe = psum(Esspe)
    Pucspe = psum(Pucspe)
    IntLengthAbove = psum(IntLengthAbove)
    k2Es = psum(k2Es)
    k2Ec = psum(k2Ec)
    !
    if(mpirank==0)  print*, '** spectra calculation finish'
    !
    !!!! Do inverse FFT
    !
    call fftw_mpi_execute_dft(backward_plan,u1c,u1c)
    call fftw_mpi_execute_dft(backward_plan,u2c,u2c)
    call fftw_mpi_execute_dft(backward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(backward_plan,u2s,u2s)
    !
    !!!! Give S-C physical energy
    Ecphy = 0.0d0
    Esphy = 0.0d0
    ucusphy = 0.0d0
    roav = 0.0d0
    Ecmax = 0.0d0
    ReUsConjUdmax = -100.d0
    ReUsConjUdmin = 100.d0
    do j=1,jm
    do i=1,im
        u1cR(i,j) = real(u1c(i,j))
        u2cR(i,j) = real(u2c(i,j))
        u1sR(i,j) = real(u1s(i,j))
        u2sR(i,j) = real(u2s(i,j))
        ucuc(i,j) = u1cR(i,j)*u1cR(i,j)+u2cR(i,j)*u2cR(i,j)
        ucus(i,j) = u1cR(i,j)*u1sR(i,j)+u2cR(i,j)*u2sR(i,j)
        usus(i,j) = u1sR(i,j)*u1sR(i,j)+u2sR(i,j)*u2sR(i,j)
        Ecmax = max(Ecmax,ucuc(i,j))
        Ecphy = Ecphy + ucuc(i,j)
        Esphy = Esphy + usus(i,j)
        ucusphy = ucusphy + ucus(i,j)
        roav = roav + rho(i,j,0)
        ReUsConjUdmax = max(ReUsConjUdmax,ReUsConjUd(i,j))
        ReUsConjUdmin = min(ReUsConjUdmin,ReUsConjUd(i,j))
      end do
    end do
    !
    roav = psum(roav)/(1.d0*ia*ja)
    Ecphy = psum(Ecphy)/(2.d0*ia*ja)
    Esphy = psum(Esphy)/(2.d0*ia*ja)
    ucusphy = psum(ucusphy)/(1.d0*ia*ja)
    Ecmax = pmax(Ecmax)
    ReUsConjUdmax = pmax(ReUsConjUdmax)
    ReUsConjUdmin = pmin(ReUsConjUdmin)
    !
    if(mpirank==0)  print*, '** physical calculation finish'
    !
    if(mpirank == 0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec'//stepname//'.dat'
      else
        outfilename = 'pp/Espec.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_a, &
                        firstline='nstep time k ES EC Eall Puc')
      do i=0,allkmax
        call listwrite(hand_a,dble(i),ES(i),EC(i),Eall(i),Puc(i))
      end do
      !
      print*,' <<< '//outfilename//'... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec_aux'//stepname//'.dat'
      else
        outfilename = 'pp/Espec_aux.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
            firstline='nstep time Ecphy Esphy ucusphy Ephyall Ecspe Esspe Pucspe Espeall Ecmax k2Es k2Ec IntLen')
      call listwrite(hand_b,Ecphy,Esphy,ucusphy,Ecphy+Esphy+ucusphy&
                    ,Ecspe,Esspe,Pucspe,Ecspe+Esspe,Ecmax,k2Es,k2Ec,IntLengthAbove/(Ecspe+Esspe))
      !
      print*,' <<< '//outfilename//'... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec_ReUsConjUd'//stepname//'.dat'
      else
        outfilename = 'pp/Espec_ReUsConjUd.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
                    firstline='nstep time ReUsConjUdmax ReUsConjUdmin')
      call listwrite(hand_b,ReUsConjUdmax,ReUsConjUdmin)
      !
      print*,' <<< '//outfilename//'... done.'
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_pspe)
    call fftw_free(c_u1c)
    call fftw_free(c_u1s)
    call fftw_free(c_u2c)
    call fftw_free(c_u2s)
    call mpistop
    
    deallocate(ES,EC,Ecount,Eall,Puc)
    deallocate(usspe,ucspe)
    deallocate(u1cR,u2cR,u1sR,u2sR,ucuc,ucus,usus,ReUsConjUd)
    deallocate(k1,k2)
    !
  end subroutine instantspectra2Davg
  !
  !
  subroutine instantspectra3Davg(thefilenumb)
    !
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km,ia,ja,ka
    use commarray, only: vel, rho, prs
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    integer,intent(in) :: thefilenumb
    character(len=128) :: infilename
    character(len=4) :: stepname
    real(8) :: u1mean,u2mean,u3mean,rhomean,prsmean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: u1spe,u2spe,u3spe,pspe
    real(8), allocatable, dimension(:) :: ES,EC,Ecount,Eall,Puc
    complex(8), allocatable, dimension(:,:,:) :: ucspe
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: u1c,u2c,u3c,u1s,u2s,u3s
    real(8), allocatable, dimension(:,:,:) :: ucuc,ucus,usus,ReUsConjUd
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
    integer :: allkmax
    real(8) :: kk,dk !wave number
    real(8) :: Ecspe,Esspe,Pucspe,Ecphy,Esphy,ucusphy,roav,Ecmax,ReUsConjUdmax,ReUsConjUdmin
    real(8) :: IntLengthAbove
    character(len=128) :: outfilename
    integer :: hand_a,hand_b
    character(len=1) :: modeio
    integer :: i,j,k,n
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_u3spe, c_pspe, c_u1c, c_u2c, c_u3c, c_u1s, c_u2s, c_u3s
    !
    call readinput
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    allkmax=ceiling(sqrt(2.d0)/3*min(min(ia,ja),ka))
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:", ka,"allkmax:",allkmax
    if(ka == 0) stop 'Please use instantspectra2D'
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !
    allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km), prs(0:im,0:jm,0:km))
    !
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='p',  var=prs(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    ! Calculate average
    u1mean = 0.0d0
    u2mean = 0.0d0
    u3mean = 0.0d0
    rhomean = 0.0d0
    prsmean = 0.0d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      u1mean = u1mean + vel(i,j,k,1)
      u2mean = u2mean + vel(i,j,k,2)
      u3mean = u3mean + vel(i,j,k,3)
      rhomean = rhomean + rho(i,j,k)
      prsmean = prsmean + prs(i,j,k)
    enddo
    enddo
    enddo
    !
    rhomean = psum(rhomean) / (1.0d0*ia*ja*ka)
    u1mean = psum(u1mean) / (1.d0*ia*ja*ka)
    u2mean = psum(u2mean) / (1.d0*ia*ja*ka)
    u3mean = psum(u3mean) / (1.d0*ia*ja*ka)
    prsmean = psum(prsmean) / (1.d0*ia*ja*ka)
    if(mpirank==0) print *, 'u1mean=',u1mean, 'u2mean=',u2mean, 'u3mean=', u3mean, 'prsmean=',prsmean
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw,kmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw,kmfftw])
    c_u3spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3spe, u3spe, [imfftw,jmfftw,kmfftw])
    c_pspe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_pspe, pspe, [imfftw,jmfftw,kmfftw])
    !
    ! Planning
    forward_plan = fftw_mpi_plan_dft_3d(kafftw, jafftw, iafftw, u1spe,u1spe, &
                    MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_3d(kafftw, jafftw, iafftw, u1spe,u1spe, &
                    MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j,k)=CMPLX(vel(i,j,k,1)-u1mean,0.d0,C_INTPTR_T);
      u2spe(i,j,k)=CMPLX(vel(i,j,k,2)-u2mean,0.d0,C_INTPTR_T);
      u3spe(i,j,k)=CMPLX(vel(i,j,k,3)-u3mean,0.d0,C_INTPTR_T);
      pspe(i,j,k)=CMPLX(prs(i,j,k)-prsmean,0.d0,C_INTPTR_T);
      !
    enddo
    enddo
    enddo
    !
    !!!! Do FFT
    !
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    call fftw_mpi_execute_dft(forward_plan,u3spe,u3spe)
    call fftw_mpi_execute_dft(forward_plan,pspe,pspe)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j,k)=u1spe(i,j,k)/(1.d0*ia*ja*ka)
      u2spe(i,j,k)=u2spe(i,j,k)/(1.d0*ia*ja*ka)
      u3spe(i,j,k)=u3spe(i,j,k)/(1.d0*ia*ja*ka)
      pspe(i,j,k)=pspe(i,j,k)/(1.d0*ia*ja*ka)
      !
    enddo
    enddo
    enddo
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j,k) = real(i-1,8)
      else if(i<=ia) then
        k1(i,j,k) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if(j <= (ja/2+1)) then
        k2(i,j,k) = real(j-1,8)
      else if(j<=ja) then
        k2(i,j,k) = real(j-ja-1,8)
      else
        print *,"Error, no wave number possible, j must smaller than ja-1 !"
      end if
      !
      if((k+k0) <= (ka/2+1)) then
        k3(i,j,k) = real(k+k0-1,8)
      else if((k+k0)<=ka) then
        k3(i,j,k) = real(k+k0-ka-1,8)
      else
        print *,"Error, no wave number possible, (k+k0) must smaller than ka-1 !"
      end if
      !
    enddo
    enddo
    enddo
    !
    !!!! Do S-C decomposition
    allocate(ucspe(1:im,1:jm,1:km))
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw,kmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw,kmfftw])
    c_u3c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3c, u3c, [imfftw,jmfftw,kmfftw])
    c_u1s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1s, u1s, [imfftw,jmfftw,kmfftw])
    c_u2s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2s, u2s, [imfftw,jmfftw,kmfftw])
    c_u3s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3s, u3s, [imfftw,jmfftw,kmfftw])
    allocate(ucuc(1:im,1:jm,1:km),ucus(1:im,1:jm,1:km),usus(1:im,1:jm,1:km))
    allocate(ReUsConjUd(1:im,1:jm,1:km))
    !
    ! 
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      kk=k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2+1.d-15
      !
      ucspe(i,j,k) = k1(i,j,k)/kk * u1spe(i,j,k) + k2(i,j,k)/kk * u2spe(i,j,k) + k3(i,j,k)/kk * u3spe(i,j,k)
      u1c(i,j,k)=  k1(i,j,k)*k1(i,j,k)/kk * u1spe(i,j,k) + k1(i,j,k)*k2(i,j,k)/kk * u2spe(i,j,k) &
                + k1(i,j,k)*k3(i,j,k)/kk * u3spe(i,j,k)
      u2c(i,j,k)=  k2(i,j,k)*k1(i,j,k)/kk * u1spe(i,j,k) + k2(i,j,k)*k2(i,j,k)/kk * u2spe(i,j,k) &
                + k2(i,j,k)*k3(i,j,k)/kk * u3spe(i,j,k)
      u3c(i,j,k)=  k3(i,j,k)*k1(i,j,k)/kk * u1spe(i,j,k) + k3(i,j,k)*k2(i,j,k)/kk * u2spe(i,j,k) &
                + k3(i,j,k)*k3(i,j,k)/kk * u3spe(i,j,k)
      u1s(i,j,k)=  u1spe(i,j,k) - u1c(i,j,k)
      u2s(i,j,k)=  u2spe(i,j,k) - u2c(i,j,k)
      u3s(i,j,k)=  u3spe(i,j,k) - u3c(i,j,k)
      !
      ReUsConjUd(i,j,k) = real(conjg(u1s(i,j,k))*ucspe(i,j,k)) + real(conjg(u2s(i,j,k))*ucspe(i,j,k))
      !
    enddo
    enddo
    enddo
    !
    !
    if(mpirank==0)  print*, '** spectral decomposition calculation finish'
    !
    !!!! Give S-C spectra and spectral energy
    !
    dk = 0.5d0
    allocate(ES(0:allkmax),EC(0:allkmax),Ecount(0:allkmax),Eall(0:allkmax),Puc(0:allkmax))
    !
    do i=0,allkmax
      ES(i) = 0.0d0
      EC(i) = 0.0d0
      Ecount(i) = 0.0d0
      Eall(i) = 0.0d0
      Puc(i) = 0.0d0
    end do
    !
    Ecspe = 0.0d0
    Esspe = 0.0d0
    Pucspe = 0.0d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      kk=dsqrt(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2+1.d-15)
      if (abs(kk-nint(kk))<=dk .and. kint(kk,dk)<=allkmax) then
        ES(kint(kk,dk)) = ES(kint(kk,dk)) + u1s(i,j,k)*conjg(u1s(i,j,k))*(kk**2) + &
                          u2s(i,j,k)*conjg(u2s(i,j,k))*(kk**2) + &
                          u3s(i,j,k)*conjg(u3s(i,j,k))*(kk**2)
        EC(kint(kk,dk)) = EC(kint(kk,dk)) + u1c(i,j,k)*conjg(u1c(i,j,k))*(kk**2) + &
                          u2c(i,j,k)*conjg(u2c(i,j,k))*(kk**2) + &
                          u3c(i,j,k)*conjg(u3c(i,j,k))*(kk**2)
        Eall(kint(kk,dk)) = Eall(kint(kk,dk)) + u1spe(i,j,k)*conjg(u1spe(i,j,k))*(kk**2) + &
                            u2spe(i,j,k)*conjg(u2spe(i,j,k))*(kk**2) + &
                            u3spe(i,j,k)*conjg(u3spe(i,j,k))*(kk**2)
        Puc(kint(kk,dk)) = Puc(kint(kk,dk)) + dimag(pspe(i,j,k)*dconjg(ucspe(i,j,k))*kk)*(kk**2)
        Ecount(kint(kk,dk)) = Ecount(kint(kk,dk)) + 1
      end if
      Ecspe = Ecspe + u1c(i,j,k)*conjg(u1c(i,j,k))/2 + &
              u2c(i,j,k)*conjg(u2c(i,j,k))/2 + &
              u3c(i,j,k)*conjg(u3c(i,j,k))/2
      Esspe = Esspe + u1s(i,j,k)*conjg(u1s(i,j,k))/2 + &
              u2s(i,j,k)*conjg(u2s(i,j,k))/2 + &
              u3s(i,j,k)*conjg(u3s(i,j,k))/2
      IntLengthAbove = IntLengthAbove + u1c(i,j,k)*conjg(u1c(i,j,k))/2/kk + &
              u2c(i,j,k)*conjg(u2c(i,j,k))/2/kk + &
              u3c(i,j,k)*conjg(u3c(i,j,k))/2/kk + &
              u1s(i,j,k)*conjg(u1s(i,j,k))/2/kk + &
              u2s(i,j,k)*conjg(u2s(i,j,k))/2/kk + &
              u3s(i,j,k)*conjg(u3s(i,j,k))/2/kk
      Pucspe = Pucspe + dimag(pspe(i,j,k)*dconjg(ucspe(i,j,k))*(kk**2))/2
    enddo
    enddo
    enddo
    !
    do i=0,allkmax
      ES(i) = psum(ES(i))
      EC(i) = psum(EC(i))
      Eall(i) = psum(Eall(i))
      Ecount(i) = psum(Ecount(i))
      Puc(i) = psum(Puc(i))
    end do
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    do i=0,allkmax
      if(Ecount(i) .ne. 0)then
        ES(i) = ES(i)/Ecount(i)*2*pi
        EC(i) = EC(i)/Ecount(i)*2*pi
        Eall(i) = Eall(i)/Ecount(i)*2*pi
        Puc(i) = - Puc(i)/Ecount(i)*4*pi
      endif
    enddo
    !
    Ecspe = psum(Ecspe)
    Esspe = psum(Esspe)
    Pucspe = psum(Pucspe)
    IntLengthAbove = psum(IntLengthAbove)
    !
    if(mpirank==0)  print*, '** spectra calculation finish'
    !
    !!!! Do inverse FFT
    !
    call fftw_mpi_execute_dft(backward_plan,u1c,u1c)
    call fftw_mpi_execute_dft(backward_plan,u2c,u2c)
    call fftw_mpi_execute_dft(backward_plan,u3c,u3c)
    call fftw_mpi_execute_dft(backward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(backward_plan,u2s,u2s)
    call fftw_mpi_execute_dft(backward_plan,u3s,u3s)
    !
    !!!! Give S-C physical energy
    Ecphy = 0.0d0
    Esphy = 0.0d0
    ucusphy = 0.0d0
    roav = 0.0d0
    Ecmax = 0.0d0
    ReUsConjUdmax = -100.d0
    ReUsConjUdmin = 100.d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      ucuc(i,j,k) = real(u1c(i,j,k))*real(u1c(i,j,k)) + real(u2c(i,j,k))*real(u2c(i,j,k)) + real(u3c(i,j,k))*real(u3c(i,j,k))
      ucus(i,j,k) = real(u1c(i,j,k))*real(u1s(i,j,k)) + real(u2c(i,j,k))*real(u2s(i,j,k)) + real(u3c(i,j,k))*real(u3s(i,j,k))
      usus(i,j,k) = real(u1s(i,j,k))*real(u1s(i,j,k)) + real(u2s(i,j,k))*real(u2s(i,j,k)) + real(u3s(i,j,k))*real(u3s(i,j,k))
      Ecmax = max(Ecmax,ucuc(i,j,k))
      Ecphy = Ecphy + ucuc(i,j,k)
      Esphy = Esphy + usus(i,j,k)
      ucusphy = ucusphy + ucus(i,j,k)
      roav = roav + rho(i,j,k)
      ReUsConjUdmax = max(ReUsConjUdmax,ReUsConjUd(i,j,k))
      ReUsConjUdmin = min(ReUsConjUdmin,ReUsConjUd(i,j,k))
    enddo
    enddo
    enddo
    !
    roav = psum(roav)/(1.d0*ia*ja*ka)
    Ecphy = psum(Ecphy)/(2.d0*ia*ja*ka)
    Esphy = psum(Esphy)/(2.d0*ia*ja*ka)
    ucusphy = psum(ucusphy)/(1.d0*ia*ja*ka)
    Ecmax = pmax(Ecmax)
    ReUsConjUdmax = pmax(ReUsConjUdmax)
    ReUsConjUdmin = pmin(ReUsConjUdmin)
    !
    if(mpirank==0)  print*, '** physical calculation finish'
    !
    if(mpirank == 0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec'//stepname//'.dat'
      else
        outfilename = 'pp/Espec.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_a, &
                        firstline='nstep time k ES EC Eall Puc')
      do i=0,allkmax
        call listwrite(hand_a,dble(i),ES(i),EC(i),Eall(i),Puc(i))
      end do
      !
      print*,' <<< '//outfilename//'... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec_aux'//stepname//'.dat'
      else
        outfilename = 'pp/Espec_aux.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
                    firstline='nstep time Ecphy Esphy ucusphy Ephyall Ecspe Esspe Pucspe Espeall Ecmax IntLen')
      call listwrite(hand_b,Ecphy,Esphy,ucusphy,Ecphy+Esphy+ucusphy&
                    ,Ecspe,Esspe,Pucspe,Ecspe+Esspe,Ecmax,IntLengthAbove/(Ecspe+Esspe))
      !
      print*,' <<< '//outfilename//'... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec_ReUsConjUd'//stepname//'.dat'
      else
        outfilename = 'pp/Espec_ReUsConjUd.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
                    firstline='nstep time ReUsConjUdmax ReUsConjUdmin')
      call listwrite(hand_b,ReUsConjUdmax,ReUsConjUdmin)
      !
      print*,' <<< '//outfilename//'... done.'
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_u3spe)
    call fftw_free(c_pspe)
    call fftw_free(c_u1c)
    call fftw_free(c_u1s)
    call fftw_free(c_u2c)
    call fftw_free(c_u2s)
    call fftw_free(c_u3c)
    call fftw_free(c_u3s)
    call mpistop
    
    deallocate(ES,EC,Ecount,Eall,Puc)
    deallocate(ucspe)
    deallocate(ucuc,ucus,usus,ReUsConjUd)
    deallocate(k1,k2,k3)
    !
  end subroutine instantspectra3Davg
  !
  subroutine instantspectra2D(thefilenumb)
    !
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km,ia,ja,ka
    use commarray, only: vel, rho, prs
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    integer,intent(in) :: thefilenumb
    character(len=128) :: infilename
    character(len=4) :: stepname
    real(8) :: u1mean,u2mean,rhomean,prsmean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1spe,u2spe,pspe
    real(8), allocatable, dimension(:) :: ES,EC,Ecount,Eall,Puc
    complex(8), allocatable, dimension(:,:) :: usspe,ucspe
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1c,u2c,u1s,u2s
    ! real(8), allocatable, dimension(:,:) :: usspeR,usspeI,ucspeR,ucspeI
    real(8), allocatable, dimension(:,:) :: u1cR,u2cR,u1sR,u2sR,ucuc,ucus,usus,ReUsConjUd
    real(8), allocatable, dimension(:,:) :: k1,k2
    integer :: allkmax
    real(8) :: k,dk !wave number
    real(8) :: Ecspe,Esspe,Pucspe,Ecphy,Esphy,ucusphy,roav,Ecmax
    real(8) :: ReUsConjUdmax,ReUsConjUdmin,k2Es,k2Ec,IntLengthAbove
    character(len=128) :: outfilename
    integer :: hand_a,hand_b
    character(len=1) :: modeio
    integer :: i,j,n
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_pspe, c_u1c, c_u2c, c_u1s, c_u2s
    !
    call readinput
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    allkmax=ceiling(sqrt(2.d0)/3*min(ia,ja))
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,"allkmax:",allkmax
    if(ka .ne. 0) stop 'Please use instantspectra3D'
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !
    allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km), prs(0:im,0:jm,0:km))
    !
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='p',  var=prs(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    ! Calculate average
    u1mean = 0.0d0
    u2mean = 0.0d0
    rhomean = 0.0d0
    prsmean = 0.0d0
    !
    do i=1,im
      do j=1,jm
        u1mean = u1mean + vel(i,j,0,1)
        u2mean = u2mean + vel(i,j,0,2)
        rhomean = rhomean + rho(i,j,0)
        prsmean = prsmean + prs(i,j,0)
      enddo
    enddo
    rhomean = psum(rhomean) / (1.0d0*ia*ja)
    u1mean = psum(u1mean) / (1.d0*ia*ja)
    u2mean = psum(u2mean) / (1.d0*ia*ja)
    prsmean = psum(prsmean) / (1.d0*ia*ja)
    if(mpirank==0) print *, 'u1mean=',u1mean, 'u2mean=',u2mean, 'prsmean=',prsmean
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw])
    c_pspe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_pspe, pspe, [imfftw,jmfftw])
    !
    ! Planning
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=CMPLX(vel(i,j,0,1)-u1mean,0.d0,C_INTPTR_T);
      u2spe(i,j)=CMPLX(vel(i,j,0,2)-u2mean,0.d0,C_INTPTR_T);
      pspe(i,j)=CMPLX(prs(i,j,0)-prsmean,0.d0,C_INTPTR_T);
      !
    end do
    end do
    !
    !!!! Do 2d FFT
    !
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    call fftw_mpi_execute_dft(forward_plan,pspe,pspe)

    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=u1spe(i,j)/(1.d0*ia*ja)
      u2spe(i,j)=u2spe(i,j)/(1.d0*ia*ja)
      pspe(i,j)=pspe(i,j)/(1.d0*ia*ja)
      !
    end do
    end do
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    do j=1,jm
    do i=1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if((j+j0) <= (ja/2+1)) then
        k2(i,j) = real(j+j0-1,8)
      else if((j+j0)<=(ja)) then
        k2(i,j) = real(j+j0-ja-1,8)
      else
        print *,"Error, no wave number possible, (j+jm) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    !
    !!!! Do S-C decomposition
    allocate(usspe(1:im,1:jm),ucspe(1:im,1:jm))
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw])
    c_u1s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1s, u1s, [imfftw,jmfftw])
    c_u2s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2s, u2s, [imfftw,jmfftw])
    allocate(u1cR(1:im,1:jm),u2cR(1:im,1:jm),u1sR(1:im,1:jm),u2sR(1:im,1:jm),&
            ucuc(1:im,1:jm),ucus(1:im,1:jm),usus(1:im,1:jm))
    allocate(ReUsConjUd(1:im,1:jm))
    !
    ! 
    !
    do j=1,jm
    do i=1,im
        k=dsqrt(k1(i,j)**2+k2(i,j)**2+1.d-15)
        usspe(i,j) = u1spe(i,j)*k2(i,j)/k - u2spe(i,j)*k1(i,j)/k
        ucspe(i,j) = u1spe(i,j)*k1(i,j)/k + u2spe(i,j)*k2(i,j)/k
        !
        u1c(i,j)=  ucspe(i,j)*k1(i,j)/k 
        u2c(i,j)=  ucspe(i,j)*k2(i,j)/k
        u1s(i,j)=  usspe(i,j)*k2(i,j)/k 
        u2s(i,j)= -usspe(i,j)*k1(i,j)/k
        !
        ReUsConjUd(i,j) = real(conjg(usspe(i,j))*ucspe(i,j))
        !
      end do
    end do
    !
    if(mpirank==0)  print*, '** spectral decomposition calculation finish'
    !
    !
    !!!! Give S-C spectra and spectral energy
    !
    dk = 0.5d0
    allocate(ES(0:allkmax),EC(0:allkmax),Ecount(0:allkmax),Eall(0:allkmax),Puc(0:allkmax))
    do i=0,allkmax
      ES(i) = 0.0d0
      EC(i) = 0.0d0
      Ecount(i) = 0.0d0
      Eall(i) = 0.0d0
      Puc(i) = 0.0d0
    end do
    !
    Ecspe = 0.0d0
    Esspe = 0.0d0
    Pucspe = 0.0d0
    k2Es = 0.d0
    k2Ec = 0.d0
    !
    do j=1,jm
    do i=1,im
        k=dsqrt(k1(i,j)**2+k2(i,j)**2+1.d-15)
        if (abs(k-nint(k))<=dk .and. kint(k,dk)<=allkmax) then
          ES(kint(k,dk)) = ES(kint(k,dk)) + usspe(i,j)*conjg(usspe(i,j))/2
          EC(kint(k,dk)) = EC(kint(k,dk)) + ucspe(i,j)*conjg(ucspe(i,j))/2
          Eall(kint(k,dk)) = Eall(kint(k,dk)) + u1spe(i,j)*conjg(u1spe(i,j))/2 + &
                              u2spe(i,j)*conjg(u2spe(i,j))/2
          Puc(kint(k,dk)) = Puc(kint(k,dk)) + dimag(pspe(i,j)*dconjg(ucspe(i,j))*k)
          Ecount(kint(k,dk)) = Ecount(kint(k,dk)) + 1
        end if
        Ecspe = Ecspe + (ucspe(i,j)*dconjg(ucspe(i,j)))/2
        Esspe = Esspe + (usspe(i,j)*dconjg(usspe(i,j)))/2
        IntLengthAbove = IntLengthAbove + usspe(i,j)*conjg(usspe(i,j))/2/k + &
                        (ucspe(i,j)*dconjg(ucspe(i,j)))/2/k
        Pucspe = Pucspe + dimag(pspe(i,j)*dconjg(ucspe(i,j))*k)/2
        k2Es = k2Es + k**2 * (usspe(i,j)*dconjg(usspe(i,j)))/2
        k2Ec = k2Ec + k**2 * (ucspe(i,j)*dconjg(ucspe(i,j)))/2
      end do
    end do
    !
    do i=0,allkmax
      ES(i) = psum(ES(i))
      EC(i) = psum(EC(i))
      Eall(i) = psum(Eall(i))
      Ecount(i) = psum(Ecount(i))
      Puc(i) = - psum(Puc(i))
    end do
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !
    !
    Ecspe = psum(Ecspe)
    Esspe = psum(Esspe)
    Pucspe = psum(Pucspe)
    IntLengthAbove = psum(IntLengthAbove)
    k2Es = psum(k2Es)
    k2Ec = psum(k2Ec)
    !
    if(mpirank==0)  print*, '** spectra calculation finish'
    !
    !!!! Do inverse FFT
    !
    call fftw_mpi_execute_dft(backward_plan,u1c,u1c)
    call fftw_mpi_execute_dft(backward_plan,u2c,u2c)
    call fftw_mpi_execute_dft(backward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(backward_plan,u2s,u2s)
    !
    !!!! Give S-C physical energy
    Ecphy = 0.0d0
    Esphy = 0.0d0
    ucusphy = 0.0d0
    roav = 0.0d0
    Ecmax = 0.0d0
    ReUsConjUdmax = -100.d0
    ReUsConjUdmin = 100.d0
    do j=1,jm
    do i=1,im
        u1cR(i,j) = real(u1c(i,j))
        u2cR(i,j) = real(u2c(i,j))
        u1sR(i,j) = real(u1s(i,j))
        u2sR(i,j) = real(u2s(i,j))
        ucuc(i,j) = u1cR(i,j)*u1cR(i,j)+u2cR(i,j)*u2cR(i,j)
        ucus(i,j) = u1cR(i,j)*u1sR(i,j)+u2cR(i,j)*u2sR(i,j)
        usus(i,j) = u1sR(i,j)*u1sR(i,j)+u2sR(i,j)*u2sR(i,j)
        Ecmax = max(Ecmax,ucuc(i,j))
        Ecphy = Ecphy + ucuc(i,j)
        Esphy = Esphy + usus(i,j)
        ucusphy = ucusphy + ucus(i,j)
        roav = roav + rho(i,j,0)
        ReUsConjUdmax = max(ReUsConjUdmax,ReUsConjUd(i,j))
        ReUsConjUdmin = min(ReUsConjUdmin,ReUsConjUd(i,j))
      end do
    end do
    !
    roav = psum(roav)/(1.d0*ia*ja)
    Ecphy = psum(Ecphy)/(2.d0*ia*ja)
    Esphy = psum(Esphy)/(2.d0*ia*ja)
    ucusphy = psum(ucusphy)/(1.d0*ia*ja)
    Ecmax = pmax(Ecmax)
    ReUsConjUdmax = pmax(ReUsConjUdmax)
    ReUsConjUdmin = pmin(ReUsConjUdmin)
    !
    if(mpirank==0)  print*, '** physical calculation finish'
    !
    if(mpirank == 0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec'//stepname//'.dat'
      else
        outfilename = 'pp/Espec.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_a, &
                        firstline='nstep time k ES EC Eall Puc')
      do i=0,allkmax
        call listwrite(hand_a,dble(i),ES(i),EC(i),Eall(i),Puc(i))
      end do
      !
      print*,' <<< '//outfilename//'... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec_aux'//stepname//'.dat'
      else
        outfilename = 'pp/Espec_aux.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
            firstline='nstep time Ecphy Esphy ucusphy Ephyall Ecspe Esspe Pucspe Espeall Ecmax k2Es k2Ec IntLen')
      call listwrite(hand_b,Ecphy,Esphy,ucusphy,Ecphy+Esphy+ucusphy&
                    ,Ecspe,Esspe,Pucspe,Ecspe+Esspe,Ecmax,k2Es,k2Ec,IntLengthAbove/(Ecspe+Esspe))
      !
      print*,' <<< '//outfilename//'... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec_ReUsConjUd'//stepname//'.dat'
      else
        outfilename = 'pp/Espec_ReUsConjUd.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
                    firstline='nstep time ReUsConjUdmax ReUsConjUdmin')
      call listwrite(hand_b,ReUsConjUdmax,ReUsConjUdmin)
      !
      print*,' <<< '//outfilename//'... done.'
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_pspe)
    call fftw_free(c_u1c)
    call fftw_free(c_u1s)
    call fftw_free(c_u2c)
    call fftw_free(c_u2s)
    call mpistop
    
    deallocate(ES,EC,Ecount,Eall,Puc)
    deallocate(usspe,ucspe)
    deallocate(u1cR,u2cR,u1sR,u2sR,ucuc,ucus,usus,ReUsConjUd)
    deallocate(k1,k2)
    !
  end subroutine instantspectra2D
  !
  !
  subroutine instantspectra3D(thefilenumb)
    !
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km,ia,ja,ka
    use commarray, only: vel, rho, prs
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    integer,intent(in) :: thefilenumb
    character(len=128) :: infilename
    character(len=4) :: stepname
    real(8) :: u1mean,u2mean,u3mean,rhomean,prsmean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: u1spe,u2spe,u3spe,pspe
    real(8), allocatable, dimension(:) :: ES,EC,Ecount,Eall,Puc
    complex(8), allocatable, dimension(:,:,:) :: ucspe
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: u1c,u2c,u3c,u1s,u2s,u3s
    real(8), allocatable, dimension(:,:,:) :: ucuc,ucus,usus,ReUsConjUd
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
    integer :: allkmax
    real(8) :: kk,dk !wave number
    real(8) :: Ecspe,Esspe,Pucspe,Ecphy,Esphy,ucusphy,roav,Ecmax,ReUsConjUdmax,ReUsConjUdmin
    real(8) :: IntLengthAbove
    character(len=128) :: outfilename
    integer :: hand_a,hand_b
    character(len=1) :: modeio
    integer :: i,j,k,n
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_u3spe, c_pspe, c_u1c, c_u2c, c_u3c, c_u1s, c_u2s, c_u3s
    !
    call readinput
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    allkmax=ceiling(sqrt(2.d0)/3*min(min(ia,ja),ka))
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:", ka,"allkmax:",allkmax
    if(ka == 0) stop 'Please use instantspectra2D'
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !
    allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km), prs(0:im,0:jm,0:km))
    !
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='p',  var=prs(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    ! Calculate average
    u1mean = 0.0d0
    u2mean = 0.0d0
    u3mean = 0.0d0
    rhomean = 0.0d0
    prsmean = 0.0d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      u1mean = u1mean + vel(i,j,k,1)
      u2mean = u2mean + vel(i,j,k,2)
      u3mean = u3mean + vel(i,j,k,3)
      rhomean = rhomean + rho(i,j,k)
      prsmean = prsmean + prs(i,j,k)
    enddo
    enddo
    enddo
    !
    rhomean = psum(rhomean) / (1.0d0*ia*ja*ka)
    u1mean = psum(u1mean) / (1.d0*ia*ja*ka)
    u2mean = psum(u2mean) / (1.d0*ia*ja*ka)
    u3mean = psum(u3mean) / (1.d0*ia*ja*ka)
    prsmean = psum(prsmean) / (1.d0*ia*ja*ka)
    if(mpirank==0) print *, 'u1mean=',u1mean, 'u2mean=',u2mean, 'u3mean=', u3mean, 'prsmean=',prsmean
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw,kmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw,kmfftw])
    c_u3spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3spe, u3spe, [imfftw,jmfftw,kmfftw])
    c_pspe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_pspe, pspe, [imfftw,jmfftw,kmfftw])
    !
    ! Planning
    forward_plan = fftw_mpi_plan_dft_3d(kafftw, jafftw, iafftw, u1spe,u1spe, &
                    MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_3d(kafftw, jafftw, iafftw, u1spe,u1spe, &
                    MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j,k)=CMPLX(vel(i,j,k,1)-u1mean,0.d0,C_INTPTR_T);
      u2spe(i,j,k)=CMPLX(vel(i,j,k,2)-u2mean,0.d0,C_INTPTR_T);
      u3spe(i,j,k)=CMPLX(vel(i,j,k,3)-u3mean,0.d0,C_INTPTR_T);
      pspe(i,j,k)=CMPLX(prs(i,j,k)-prsmean,0.d0,C_INTPTR_T);
      !
    enddo
    enddo
    enddo
    !
    !!!! Do FFT
    !
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    call fftw_mpi_execute_dft(forward_plan,u3spe,u3spe)
    call fftw_mpi_execute_dft(forward_plan,pspe,pspe)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j,k)=u1spe(i,j,k)/(1.d0*ia*ja*ka)
      u2spe(i,j,k)=u2spe(i,j,k)/(1.d0*ia*ja*ka)
      u3spe(i,j,k)=u3spe(i,j,k)/(1.d0*ia*ja*ka)
      pspe(i,j,k)=pspe(i,j,k)/(1.d0*ia*ja*ka)
      !
    enddo
    enddo
    enddo
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j,k) = real(i-1,8)
      else if(i<=ia) then
        k1(i,j,k) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if(j <= (ja/2+1)) then
        k2(i,j,k) = real(j-1,8)
      else if(j<=ja) then
        k2(i,j,k) = real(j-ja-1,8)
      else
        print *,"Error, no wave number possible, j must smaller than ja-1 !"
      end if
      !
      if((k+k0) <= (ka/2+1)) then
        k3(i,j,k) = real(k+k0-1,8)
      else if((k+k0)<=ka) then
        k3(i,j,k) = real(k+k0-ka-1,8)
      else
        print *,"Error, no wave number possible, (k+k0) must smaller than ka-1 !"
      end if
      !
    enddo
    enddo
    enddo
    !
    !!!! Do S-C decomposition
    allocate(ucspe(1:im,1:jm,1:km))
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw,kmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw,kmfftw])
    c_u3c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3c, u3c, [imfftw,jmfftw,kmfftw])
    c_u1s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1s, u1s, [imfftw,jmfftw,kmfftw])
    c_u2s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2s, u2s, [imfftw,jmfftw,kmfftw])
    c_u3s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3s, u3s, [imfftw,jmfftw,kmfftw])
    allocate(ucuc(1:im,1:jm,1:km),ucus(1:im,1:jm,1:km),usus(1:im,1:jm,1:km))
    allocate(ReUsConjUd(1:im,1:jm,1:km))
    !
    ! 
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      kk=k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2+1.d-15
      !
      ucspe(i,j,k) = k1(i,j,k)/kk * u1spe(i,j,k) + k2(i,j,k)/kk * u2spe(i,j,k) + k3(i,j,k)/kk * u3spe(i,j,k)
      u1c(i,j,k)=  k1(i,j,k)*k1(i,j,k)/kk * u1spe(i,j,k) + k1(i,j,k)*k2(i,j,k)/kk * u2spe(i,j,k) &
                + k1(i,j,k)*k3(i,j,k)/kk * u3spe(i,j,k)
      u2c(i,j,k)=  k2(i,j,k)*k1(i,j,k)/kk * u1spe(i,j,k) + k2(i,j,k)*k2(i,j,k)/kk * u2spe(i,j,k) &
                + k2(i,j,k)*k3(i,j,k)/kk * u3spe(i,j,k)
      u3c(i,j,k)=  k3(i,j,k)*k1(i,j,k)/kk * u1spe(i,j,k) + k3(i,j,k)*k2(i,j,k)/kk * u2spe(i,j,k) &
                + k3(i,j,k)*k3(i,j,k)/kk * u3spe(i,j,k)
      u1s(i,j,k)=  u1spe(i,j,k) - u1c(i,j,k)
      u2s(i,j,k)=  u2spe(i,j,k) - u2c(i,j,k)
      u3s(i,j,k)=  u3spe(i,j,k) - u3c(i,j,k)
      !
      ReUsConjUd(i,j,k) = real(conjg(u1s(i,j,k))*ucspe(i,j,k)) + real(conjg(u2s(i,j,k))*ucspe(i,j,k))
      !
    enddo
    enddo
    enddo
    !
    !
    if(mpirank==0)  print*, '** spectral decomposition calculation finish'
    !
    !!!! Give S-C spectra and spectral energy
    !
    dk = 0.5d0
    allocate(ES(0:allkmax),EC(0:allkmax),Ecount(0:allkmax),Eall(0:allkmax),Puc(0:allkmax))
    !
    do i=0,allkmax
      ES(i) = 0.0d0
      EC(i) = 0.0d0
      Ecount(i) = 0.0d0
      Eall(i) = 0.0d0
      Puc(i) = 0.0d0
    end do
    !
    Ecspe = 0.0d0
    Esspe = 0.0d0
    Pucspe = 0.0d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      kk=dsqrt(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2+1.d-15)
      if (abs(kk-nint(kk))<=dk .and. kint(kk,dk)<=allkmax) then
        ES(kint(kk,dk)) = ES(kint(kk,dk)) + u1s(i,j,k)*conjg(u1s(i,j,k))/2 + &
                          u2s(i,j,k)*conjg(u2s(i,j,k))/2 + &
                          u3s(i,j,k)*conjg(u3s(i,j,k))/2
        EC(kint(kk,dk)) = EC(kint(kk,dk)) + u1c(i,j,k)*conjg(u1c(i,j,k))/2 + &
                          u2c(i,j,k)*conjg(u2c(i,j,k))/2 + &
                          u3c(i,j,k)*conjg(u3c(i,j,k))/2
        Eall(kint(kk,dk)) = Eall(kint(kk,dk)) + u1spe(i,j,k)*conjg(u1spe(i,j,k))/2 + &
                            u2spe(i,j,k)*conjg(u2spe(i,j,k))/2 + &
                            u3spe(i,j,k)*conjg(u3spe(i,j,k))/2
        Puc(kint(kk,dk)) = Puc(kint(kk,dk)) + dimag(pspe(i,j,k)*dconjg(ucspe(i,j,k))*kk)
        Ecount(kint(kk,dk)) = Ecount(kint(kk,dk)) + 1
      end if
      Ecspe = Ecspe + u1c(i,j,k)*conjg(u1c(i,j,k))/2 + &
              u2c(i,j,k)*conjg(u2c(i,j,k))/2 + &
              u3c(i,j,k)*conjg(u3c(i,j,k))/2
      Esspe = Esspe + u1s(i,j,k)*conjg(u1s(i,j,k))/2 + &
              u2s(i,j,k)*conjg(u2s(i,j,k))/2 + &
              u3s(i,j,k)*conjg(u3s(i,j,k))/2
      IntLengthAbove = IntLengthAbove + u1c(i,j,k)*conjg(u1c(i,j,k))/2/kk + &
              u2c(i,j,k)*conjg(u2c(i,j,k))/2/kk + &
              u3c(i,j,k)*conjg(u3c(i,j,k))/2/kk + &
              u1s(i,j,k)*conjg(u1s(i,j,k))/2/kk + &
              u2s(i,j,k)*conjg(u2s(i,j,k))/2/kk + &
              u3s(i,j,k)*conjg(u3s(i,j,k))/2/kk
      Pucspe = Pucspe + dimag(pspe(i,j,k)*dconjg(ucspe(i,j,k))*(kk**2))/2
    enddo
    enddo
    enddo
    !
    do i=0,allkmax
      ES(i) = psum(ES(i))
      EC(i) = psum(EC(i))
      Eall(i) = psum(Eall(i))
      Ecount(i) = psum(Ecount(i))
      Puc(i) = - psum(Puc(i))
    end do
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !
    !
    Ecspe = psum(Ecspe)
    Esspe = psum(Esspe)
    Pucspe = psum(Pucspe)
    IntLengthAbove = psum(IntLengthAbove)
    !
    if(mpirank==0)  print*, '** spectra calculation finish'
    !
    !!!! Do inverse FFT
    !
    call fftw_mpi_execute_dft(backward_plan,u1c,u1c)
    call fftw_mpi_execute_dft(backward_plan,u2c,u2c)
    call fftw_mpi_execute_dft(backward_plan,u3c,u3c)
    call fftw_mpi_execute_dft(backward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(backward_plan,u2s,u2s)
    call fftw_mpi_execute_dft(backward_plan,u3s,u3s)
    !
    !!!! Give S-C physical energy
    Ecphy = 0.0d0
    Esphy = 0.0d0
    ucusphy = 0.0d0
    roav = 0.0d0
    Ecmax = 0.0d0
    ReUsConjUdmax = -100.d0
    ReUsConjUdmin = 100.d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      ucuc(i,j,k) = real(u1c(i,j,k))*real(u1c(i,j,k)) + real(u2c(i,j,k))*real(u2c(i,j,k)) + real(u3c(i,j,k))*real(u3c(i,j,k))
      ucus(i,j,k) = real(u1c(i,j,k))*real(u1s(i,j,k)) + real(u2c(i,j,k))*real(u2s(i,j,k)) + real(u3c(i,j,k))*real(u3s(i,j,k))
      usus(i,j,k) = real(u1s(i,j,k))*real(u1s(i,j,k)) + real(u2s(i,j,k))*real(u2s(i,j,k)) + real(u3s(i,j,k))*real(u3s(i,j,k))
      Ecmax = max(Ecmax,ucuc(i,j,k))
      Ecphy = Ecphy + ucuc(i,j,k)
      Esphy = Esphy + usus(i,j,k)
      ucusphy = ucusphy + ucus(i,j,k)
      roav = roav + rho(i,j,k)
      ReUsConjUdmax = max(ReUsConjUdmax,ReUsConjUd(i,j,k))
      ReUsConjUdmin = min(ReUsConjUdmin,ReUsConjUd(i,j,k))
    enddo
    enddo
    enddo
    !
    roav = psum(roav)/(1.d0*ia*ja*ka)
    Ecphy = psum(Ecphy)/(2.d0*ia*ja*ka)
    Esphy = psum(Esphy)/(2.d0*ia*ja*ka)
    ucusphy = psum(ucusphy)/(1.d0*ia*ja*ka)
    Ecmax = pmax(Ecmax)
    ReUsConjUdmax = pmax(ReUsConjUdmax)
    ReUsConjUdmin = pmin(ReUsConjUdmin)
    !
    if(mpirank==0)  print*, '** physical calculation finish'
    !
    if(mpirank == 0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec'//stepname//'.dat'
      else
        outfilename = 'pp/Espec.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_a, &
                        firstline='nstep time k ES EC Eall Puc')
      do i=0,allkmax
        call listwrite(hand_a,dble(i),ES(i),EC(i),Eall(i),Puc(i))
      end do
      !
      print*,' <<< '//outfilename//'... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec_aux'//stepname//'.dat'
      else
        outfilename = 'pp/Espec_aux.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
                    firstline='nstep time Ecphy Esphy ucusphy Ephyall Ecspe Esspe Pucspe Espeall Ecmax IntLen')
      call listwrite(hand_b,Ecphy,Esphy,ucusphy,Ecphy+Esphy+ucusphy&
                    ,Ecspe,Esspe,Pucspe,Ecspe+Esspe,Ecmax,IntLengthAbove/(Ecspe+Esspe))
      !
      print*,' <<< '//outfilename//'... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec_ReUsConjUd'//stepname//'.dat'
      else
        outfilename = 'pp/Espec_ReUsConjUd.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
                    firstline='nstep time ReUsConjUdmax ReUsConjUdmin')
      call listwrite(hand_b,ReUsConjUdmax,ReUsConjUdmin)
      !
      print*,' <<< '//outfilename//'... done.'
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_u3spe)
    call fftw_free(c_pspe)
    call fftw_free(c_u1c)
    call fftw_free(c_u1s)
    call fftw_free(c_u2c)
    call fftw_free(c_u2s)
    call fftw_free(c_u3c)
    call fftw_free(c_u3s)
    call mpistop
    
    deallocate(ES,EC,Ecount,Eall,Puc)
    deallocate(ucspe)
    deallocate(ucuc,ucus,usus,ReUsConjUd)
    deallocate(k1,k2,k3)
    !
  end subroutine instantspectra3D
  !
  subroutine initparam2D
    !
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput,readic
    use fftwlink
    use commvar,only : im,jm,km,ia,ja,ka,ickmax
    use commarray, only: vel
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    character(len=128) :: infilename
    character(len=4) :: stepname
    integer :: allkmax
    character(len=1) :: modeio
    !
    real(8) :: u1mean,u2mean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1spe,u2spe
    type(C_PTR) :: forward_plan, c_u1spe, c_u2spe
    real(8), allocatable, dimension(:,:) :: k1,k2
    !
    real(8) :: k,dk !wave number
    real(8) :: Espeall,TauAbove,urms,tau,L
    integer :: i,j,n
    
    !
    call readinput
    !
    call readic
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    allkmax=ceiling(sqrt(2.d0)/3*min(ia,ja))
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,"allkmax:",allkmax
    if(ka .ne. 0) stop 'Please use instantspectra3D'
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !
    allocate(vel(0:im,0:jm,0:km,1:2))
    !
    !
    infilename='datin/flowini2d.h5'
    !
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    ! Calculate average
    u1mean = 0.0d0
    u2mean = 0.0d0
    !
    do i=1,im
      do j=1,jm
        u1mean = u1mean + vel(i,j,0,1)
        u2mean = u2mean + vel(i,j,0,2)
      enddo
    enddo
    u1mean = psum(u1mean) / (1.d0*ia*ja)
    u2mean = psum(u2mean) / (1.d0*ia*ja)
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw])
    !
    ! Planning
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    !
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=CMPLX(vel(i,j,0,1)-u1mean,0.d0,C_INTPTR_T);
      u2spe(i,j)=CMPLX(vel(i,j,0,2)-u2mean,0.d0,C_INTPTR_T);
      !
    end do
    end do
    !
    !!!! Do 2d FFT
    !
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    !
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=u1spe(i,j)/(1.d0*ia*ja)
      u2spe(i,j)=u2spe(i,j)/(1.d0*ia*ja)
      !
    end do
    end do
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    !
    do j=1,jm
    do i=1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if((j+j0) <= (ja/2+1)) then
        k2(i,j) = real(j+j0-1,8)
      else if((j+j0)<=(ja)) then
        k2(i,j) = real(j+j0-ja-1,8)
      else
        print *,"Error, no wave number possible, (j+jm) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    !
    !
    dk = 0.5d0
    !
    Espeall = 0.d0
    TauAbove = 0.d0
    !
    do j=1,jm
    do i=1,im
        k=dsqrt(k1(i,j)**2+k2(i,j)**2+1.d-15)
        Espeall = Espeall + (u1spe(i,j)*dconjg(u1spe(i,j)))/2 + &
                            (u2spe(i,j)*dconjg(u2spe(i,j)))/2
        TauAbove = TauAbove + (u1spe(i,j)*dconjg(u1spe(i,j)))/k/2 + &
                              (u2spe(i,j)*dconjg(u2spe(i,j)))/k/2
      end do
    end do
    !
    !
    TauAbove = psum(TauAbove)
    Espeall = psum(Espeall)
    !
    urms = sqrt(Espeall)
    L = TauAbove/urms/urms*pi/2
    tau = L/urms
    !
    if(mpirank==0)  print*, '** spectra calculation finish'
    !
    !
    if(mpirank == 0) then
      !
      print *, 'urms:', urms
      print *, 'L:', L
      print *, 'tau:', tau
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call mpistop
    !
    deallocate(k1,k2)
    !
  end subroutine initparam2D
  !
  !
  subroutine initparam3D
    !
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput,readic
    use fftwlink
    use commvar,only : im,jm,km,ia,ja,ka,ickmax
    use commarray, only: vel
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    character(len=128) :: infilename
    character(len=4) :: stepname
    integer :: allkmax
    character(len=1) :: modeio
    !
    real(8) :: u1mean,u2mean,u3mean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: u1spe,u2spe,u3spe
    type(C_PTR) :: forward_plan, c_u1spe, c_u2spe, c_u3spe
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
    !
    real(8) :: kk,dk !wave number
    real(8) :: Espeall,TauAbove,urms,tau,L
    integer :: i,j,k,n
    !
    call readinput
    !
    call readic
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    allkmax=ceiling(sqrt(2.d0)/3*min(min(ia,ja),ka))
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:", ka,"allkmax:",allkmax
    if(ka == 0) stop 'Please use instantspectra2D'
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !
    allocate(vel(0:im,0:jm,0:km,1:3))
    !
    !
    !
    infilename='datin/flowini3d.h5'
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    ! Calculate average
    u1mean = 0.0d0
    u2mean = 0.0d0
    u3mean = 0.0d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      u1mean = u1mean + vel(i,j,k,1)
      u2mean = u2mean + vel(i,j,k,2)
      u3mean = u3mean + vel(i,j,k,3)
    enddo
    enddo
    enddo
    !
    u1mean = psum(u1mean) / (1.d0*ia*ja*ka)
    u2mean = psum(u2mean) / (1.d0*ia*ja*ka)
    u3mean = psum(u3mean) / (1.d0*ia*ja*ka)
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw,kmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw,kmfftw])
    c_u3spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3spe, u3spe, [imfftw,jmfftw,kmfftw])
    !
    ! Planning
    forward_plan = fftw_mpi_plan_dft_3d(kafftw, jafftw, iafftw, u1spe,u1spe, &
                    MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j,k)=CMPLX(vel(i,j,k,1)-u1mean,0.d0,C_INTPTR_T);
      u2spe(i,j,k)=CMPLX(vel(i,j,k,2)-u2mean,0.d0,C_INTPTR_T);
      u3spe(i,j,k)=CMPLX(vel(i,j,k,3)-u3mean,0.d0,C_INTPTR_T);
      !
    enddo
    enddo
    enddo
    !
    !!!! Do FFT
    !
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    call fftw_mpi_execute_dft(forward_plan,u3spe,u3spe)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j,k)=u1spe(i,j,k)/(1.d0*ia*ja*ka)
      u2spe(i,j,k)=u2spe(i,j,k)/(1.d0*ia*ja*ka)
      u3spe(i,j,k)=u3spe(i,j,k)/(1.d0*ia*ja*ka)
      !
    enddo
    enddo
    enddo
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j,k) = real(i-1,8)
      else if(i<=ia) then
        k1(i,j,k) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if(j <= (ja/2+1)) then
        k2(i,j,k) = real(j-1,8)
      else if(j<=ja) then
        k2(i,j,k) = real(j-ja-1,8)
      else
        print *,"Error, no wave number possible, j must smaller than ja-1 !"
      end if
      !
      if((k+k0) <= (ka/2+1)) then
        k3(i,j,k) = real(k+k0-1,8)
      else if((k+k0)<=ka) then
        k3(i,j,k) = real(k+k0-ka-1,8)
      else
        print *,"Error, no wave number possible, (k+k0) must smaller than ka-1 !"
      end if
      !
    enddo
    enddo
    enddo
    !
    dk = 0.5d0
    !
    Espeall = 0.d0
    TauAbove = 0.d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      kk=dsqrt(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2+1.d-15)
      !
      Espeall = Espeall + (u1spe(i,j,k)*dconjg(u1spe(i,j,k)))/2 + &
                            (u2spe(i,j,k)*dconjg(u2spe(i,j,k)))/2+&
                            (u3spe(i,j,k)*dconjg(u3spe(i,j,k)))/2
      TauAbove = TauAbove + (u1spe(i,j,k)*dconjg(u1spe(i,j,k)))/2/kk + &
                            (u2spe(i,j,k)*dconjg(u2spe(i,j,k)))/2/kk + &
                            (u3spe(i,j,k)*dconjg(u3spe(i,j,k)))/2/kk
      !
    enddo
    enddo
    enddo
    !
    !
    if(mpirank==0)  print*, '** spectral decomposition calculation finish'
    !
    !!!! Give S-C spectra and spectral energy
    !
    TauAbove = psum(TauAbove)
    Espeall = psum(Espeall)
    !
    urms = sqrt(2.d0*Espeall/3.d0)
    L = TauAbove/urms/urms*pi/2
    tau = L/urms
    !
    if(mpirank==0)  print*, '** spectra calculation finish'
    !
    !
    if(mpirank == 0) then
      !
      print *, 'urms:', urms
      print *, 'L1:', sqrt(2*pi)/ickmax
      print *, 'L:', L
      print *, 'tau:', tau
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_u3spe)
    call mpistop
    !
    deallocate(k1,k2,k3)
    !
  end subroutine initparam3D
  !
  subroutine instanttriad2D(thefilenumb)
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km
    use commarray, only: vel, rho, prs
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    integer,intent(in) :: thefilenumb
    character(len=128) :: infilename
    character(len=4) :: stepname
    real(8) :: u1mean,u2mean,rhomean,prsmean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1spe,u2spe,pspe
    complex(8), allocatable, dimension(:,:) :: usspe,ucspe
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: theta,u1c,u2c,u1s,u2s
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: T11,T12,T21,T22,Tstheta1,Tstheta2,Tdtheta1,Tdtheta2
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: T11s,T12s,T21s,T22s,T11d,T12d,T21d,T22d
    real(8), allocatable, dimension(:) :: Ts,Tstheta,Td,Tdd,Tcount,Tss,Tdstheta,Tddtheta
    real(8), allocatable, dimension(:,:) :: k1,k2,u1,u2
    integer :: allkmax
    real(8) :: k,dk,p,q,kx,ky !wave number
    character(len=128) :: outfilename
    integer :: hand_a,hand_b
    character(len=1) :: modeio
    integer :: i,j,n,s,t
    real(8) :: k2Ts, k2Tc,Tsall,Tcall, Tssall, Tccall, k2Tss, k2Tcc
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_pspe, c_theta,c_u1c,c_u2c,c_u1s,c_u2s
    type(C_PTR) :: c_T11,c_T12,c_T21,c_T22,c_Tstheta1,c_Tstheta2,c_Tdtheta1,c_Tdtheta2
    type(C_PTR) :: c_T11s,c_T12s,c_T21s,c_T22s,c_T11d,c_T12d,c_T21d,c_T22d
    call readinput
    !
    modeio='h'
    dk = 0.5d0
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    allkmax=ceiling(1.d0*min(ia,ja))
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,"allkmax:",allkmax
    if(ka .ne. 0) stop 'Please use instantspectra3D'
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km), prs(0:im,0:jm,0:km))
    !
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='p',  var=prs(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    allocate(u1(1:im,1:jm),u2(1:im,1:jm))
    ! Calculate favre average
    u1mean = 0.0d0
    u2mean = 0.0d0
    rhomean = 0.0d0
    prsmean = 0.0d0
    !
    do i=1,im
      do j=1,jm
        u1mean = u1mean + vel(i,j,0,1)
        u2mean = u2mean + vel(i,j,0,2)
        rhomean = rhomean + rho(i,j,0)
        prsmean = prsmean + prs(i,j,0)
      enddo
    enddo
    !
    rhomean = psum(rhomean) / (1.0d0*ia*ja)
    u1mean = psum(u1mean) / (1.d0*ia*ja)
    u2mean = psum(u2mean) / (1.d0*ia*ja)
    prsmean = psum(prsmean) / (1.d0*ia*ja)
    if(mpirank==0) print *, 'u1mean=',u1mean, 'u2mean=',u2mean, 'prsmean=',prsmean
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw])
    c_pspe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_pspe, pspe, [imfftw,jmfftw])
    !
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    !!!! Convert to complex
    do j=1,jm
    do i=1,im
      !
      u1(i,j) = vel(i,j,0,1)-u1mean
      u2(i,j) = vel(i,j,0,2)-u2mean
      u1spe(i,j)=CMPLX(u1(i,j),0.d0,C_INTPTR_T)
      u2spe(i,j)=CMPLX(u2(i,j),0.d0,C_INTPTR_T)
      pspe(i,j)=CMPLX(prs(i,j,0)-prsmean,0.d0,C_INTPTR_T)
      !
    end do
    end do
    !
    !!!! Do 2d FFT
    !
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    call fftw_mpi_execute_dft(forward_plan,pspe,pspe)
    !
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=u1spe(i,j)/(1.d0*ia*ja)
      u2spe(i,j)=u2spe(i,j)/(1.d0*ia*ja)
      pspe(i,j)=pspe(i,j)/(1.d0*ia*ja)
      !
    end do
    end do
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    do j=1,jm
    do i=1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if((j+j0) <= (ja/2+1)) then
        k2(i,j) = real(j+j0-1,8)
      else if((j+j0)<=(ja)) then
        k2(i,j) = real(j+j0-ja-1,8)
      else
        print *,"Error, no wave number possible, (j+j0) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    !
    allocate(usspe(1:im,1:jm),ucspe(1:im,1:jm))
    !
    c_theta = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_theta, theta, [imfftw,jmfftw])
    c_u1s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1s, u1s, [imfftw,jmfftw])
    c_u2s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2s, u2s, [imfftw,jmfftw])
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw])
    !
    !!!! Do S-C decomposition
    do j=1,jm
      do i=1,im
        !
        k=dsqrt(k1(i,j)**2+k2(i,j)**2+1.d-15)
        usspe(i,j) = u1spe(i,j)*k2(i,j)/k - u2spe(i,j)*k1(i,j)/k
        ucspe(i,j) = u1spe(i,j)*k1(i,j)/k + u2spe(i,j)*k2(i,j)/k
        !
        u1c(i,j)=  ucspe(i,j)*k1(i,j)/k 
        u2c(i,j)=  ucspe(i,j)*k2(i,j)/k
        u1s(i,j)=  usspe(i,j)*k2(i,j)/k 
        u2s(i,j)= -usspe(i,j)*k1(i,j)/k
        !
        theta(i,j)= CMPLX(0.d0,1.d0,C_INTPTR_T) * (ucspe(i,j) * k)
        !
      end do
    end do
    !
    !!!! Do inverse FFT
    !
    call fftw_mpi_execute_dft(backward_plan,theta,theta)
    call fftw_mpi_execute_dft(backward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(backward_plan,u2s,u2s)
    call fftw_mpi_execute_dft(backward_plan,u1c,u1c)
    call fftw_mpi_execute_dft(backward_plan,u2c,u2c)
    !
    c_T11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11,T11, [imfftw,jmfftw])
    c_T12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12,T12, [imfftw,jmfftw])
    c_T21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21,T21, [imfftw,jmfftw])
    c_T22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22,T22, [imfftw,jmfftw])
    !
    c_T11s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11s,T11s, [imfftw,jmfftw])
    c_T12s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12s,T12s, [imfftw,jmfftw])
    c_T21s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21s,T21s, [imfftw,jmfftw])
    c_T22s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22s,T22s, [imfftw,jmfftw])
    !
    c_T11d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11d,T11d, [imfftw,jmfftw])
    c_T12d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12d,T12d, [imfftw,jmfftw])
    c_T21d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21d,T21d, [imfftw,jmfftw])
    c_T22d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22d,T22d, [imfftw,jmfftw])
    !
    c_Tstheta1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tstheta1,Tstheta1, [imfftw,jmfftw])
    c_Tstheta2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tstheta2,Tstheta2, [imfftw,jmfftw])
    c_Tdtheta1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tdtheta1,Tdtheta1, [imfftw,jmfftw])
    c_Tdtheta2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tdtheta2,Tdtheta2, [imfftw,jmfftw])
    !
    ! This part is in physical space
    do j=1,jm
    do i=1,im
      T11(i,j) = u1(i,j)*u1(i,j)
      T12(i,j) = u1(i,j)*u2(i,j)
      T21(i,j) = u2(i,j)*u1(i,j)
      T22(i,j) = u2(i,j)*u2(i,j)
      !
      T11s(i,j) = u1s(i,j)*u1s(i,j)
      T12s(i,j) = u1s(i,j)*u2s(i,j)
      T21s(i,j) = u2s(i,j)*u1s(i,j)
      T22s(i,j) = u2s(i,j)*u2s(i,j)
      !
      T11d(i,j) = u1c(i,j)*u1c(i,j)
      T12d(i,j) = u1c(i,j)*u2c(i,j)
      T21d(i,j) = u2c(i,j)*u1c(i,j)
      T22d(i,j) = u2c(i,j)*u2c(i,j)
      !
      Tstheta1(i,j) = theta(i,j)*u1s(i,j)
      Tstheta2(i,j) = theta(i,j)*u2s(i,j)
      Tdtheta1(i,j) = theta(i,j)*u1c(i,j)
      Tdtheta2(i,j) = theta(i,j)*u2c(i,j)
    enddo
    enddo
    !
    call fftw_mpi_execute_dft(forward_plan,T11,T11)
    call fftw_mpi_execute_dft(forward_plan,T12,T12)
    call fftw_mpi_execute_dft(forward_plan,T21,T21)
    call fftw_mpi_execute_dft(forward_plan,T22,T22)
    !
    call fftw_mpi_execute_dft(forward_plan,T11s,T11s)
    call fftw_mpi_execute_dft(forward_plan,T12s,T12s)
    call fftw_mpi_execute_dft(forward_plan,T21s,T21s)
    call fftw_mpi_execute_dft(forward_plan,T22s,T22s)
    !
    call fftw_mpi_execute_dft(forward_plan,T11d,T11d)
    call fftw_mpi_execute_dft(forward_plan,T12d,T12d)
    call fftw_mpi_execute_dft(forward_plan,T21d,T21d)
    call fftw_mpi_execute_dft(forward_plan,T22d,T22d)
    !
    call fftw_mpi_execute_dft(forward_plan,Tstheta1,Tstheta1)
    call fftw_mpi_execute_dft(forward_plan,Tstheta2,Tstheta2)
    call fftw_mpi_execute_dft(forward_plan,Tdtheta1,Tdtheta1)
    call fftw_mpi_execute_dft(forward_plan,Tdtheta2,Tdtheta2)
    !
    call fftw_mpi_execute_dft(forward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(forward_plan,u2s,u2s)
    call fftw_mpi_execute_dft(forward_plan,u1c,u1c)
    call fftw_mpi_execute_dft(forward_plan,u2c,u2c)
    !
    do j=1,jm
    do i=1,im
      T11(i,j)=T11(i,j)/(1.d0*ia*ja)
      T12(i,j)=T12(i,j)/(1.d0*ia*ja)
      T21(i,j)=T21(i,j)/(1.d0*ia*ja)
      T22(i,j)=T22(i,j)/(1.d0*ia*ja)
      !
      T11s(i,j)=T11s(i,j)/(1.d0*ia*ja)
      T12s(i,j)=T12s(i,j)/(1.d0*ia*ja)
      T21s(i,j)=T21s(i,j)/(1.d0*ia*ja)
      T22s(i,j)=T22s(i,j)/(1.d0*ia*ja)
      !
      T11d(i,j)=T11d(i,j)/(1.d0*ia*ja)
      T12d(i,j)=T12d(i,j)/(1.d0*ia*ja)
      T21d(i,j)=T21d(i,j)/(1.d0*ia*ja)
      T22d(i,j)=T22d(i,j)/(1.d0*ia*ja)
      !
      Tstheta1(i,j)=Tstheta1(i,j)/(1.d0*ia*ja)
      Tstheta2(i,j)=Tstheta2(i,j)/(1.d0*ia*ja)
      Tdtheta1(i,j)=Tdtheta1(i,j)/(1.d0*ia*ja)
      Tdtheta2(i,j)=Tdtheta2(i,j)/(1.d0*ia*ja)
      !
      u1s(i,j)=u1s(i,j)/(1.d0*ia*ja)
      u2s(i,j)=u2s(i,j)/(1.d0*ia*ja)
      u1c(i,j)=u1c(i,j)/(1.d0*ia*ja)
      u2c(i,j)=u2c(i,j)/(1.d0*ia*ja)
    enddo
    enddo
    !
    !
    allocate(Ts(0:allkmax),Tstheta(0:allkmax),Td(0:allkmax),Tdstheta(0:allkmax))
    allocate(Tddtheta(0:allkmax),Tcount(0:allkmax),Tss(0:allkmax),Tdd(0:allkmax))
    !
    do i=0,allkmax
      Ts(i)= 0.0d0
      Tstheta(i)= 0.0d0
      Td(i)= 0.0d0
      Tddtheta(i)= 0.0d0
      Tdstheta(i)= 0.0d0
      Tcount(i)=0
      Tss(i)=0.d0
      Tdd(i)=0.d0
    end do
    !
    k2Ts = 0.d0
    k2Tc = 0.d0
    Tsall = 0.d0
    Tcall = 0.d0
    Tssall = 0.d0
    Tccall = 0.d0
    k2Tss = 0.d0
    k2Tcc = 0.d0
    !
    do j = 1,jm
    do i = 1,im
      kx = k1(i,j)
      ky = k2(i,j)
      k=dsqrt(kx**2+ky**2+1.d-15)
      !
      if(kint(k,dk)==0)then
        print *, u1s(i,j), u2s(i,j), T11s(i,j), T12s(i,j), T21s(i,j), T22s(i,j)
      endif
      !
      if(kint(k,dk)<= allkmax)then
        Ts(kint(k,dk)) = Ts(kint(k,dk)) + &
              ProjectP3(1,1,1,kx,ky) * dimag(u1s(i,j)*conjg(T11(i,j))) + &
              ProjectP3(2,1,1,kx,ky) * dimag(u2s(i,j)*conjg(T11(i,j))) + &
              ProjectP3(1,1,2,kx,ky) * dimag(u1s(i,j)*conjg(T12(i,j))) + &
              ProjectP3(2,1,2,kx,ky) * dimag(u2s(i,j)*conjg(T12(i,j))) + &
              ProjectP3(1,2,1,kx,ky) * dimag(u1s(i,j)*conjg(T21(i,j))) + &
              ProjectP3(2,2,1,kx,ky) * dimag(u2s(i,j)*conjg(T21(i,j))) + &
              ProjectP3(1,2,2,kx,ky) * dimag(u1s(i,j)*conjg(T22(i,j))) + &
              ProjectP3(2,2,2,kx,ky) * dimag(u2s(i,j)*conjg(T22(i,j)))

        Tss(kint(k,dk)) = Tss(kint(k,dk)) + &
              ProjectP3(1,1,1,kx,ky) * dimag(u1s(i,j)*conjg(T11s(i,j))) + &
              ProjectP3(2,1,1,kx,ky) * dimag(u2s(i,j)*conjg(T11s(i,j))) + &
              ProjectP3(1,1,2,kx,ky) * dimag(u1s(i,j)*conjg(T12s(i,j))) + &
              ProjectP3(2,1,2,kx,ky) * dimag(u2s(i,j)*conjg(T12s(i,j))) + &
              ProjectP3(1,2,1,kx,ky) * dimag(u1s(i,j)*conjg(T21s(i,j))) + &
              ProjectP3(2,2,1,kx,ky) * dimag(u2s(i,j)*conjg(T21s(i,j))) + &
              ProjectP3(1,2,2,kx,ky) * dimag(u1s(i,j)*conjg(T22s(i,j))) + &
              ProjectP3(2,2,2,kx,ky) * dimag(u2s(i,j)*conjg(T22s(i,j)))
            !
        Tstheta(kint(k,dk)) = Tstheta(kint(k,dk)) + &
              ProjectP2(1,1,kx,ky) * dreal(u1s(i,j)*conjg(Tstheta1(i,j)+Tdtheta1(i,j))) + &
              ProjectP2(1,2,kx,ky) * dreal(u1s(i,j)*conjg(Tstheta2(i,j)+Tdtheta2(i,j))) + &
              ProjectP2(2,1,kx,ky) * dreal(u2s(i,j)*conjg(Tstheta1(i,j)+Tstheta1(i,j))) + &
              ProjectP2(2,2,kx,ky) * dreal(u2s(i,j)*conjg(Tstheta2(i,j)+Tdtheta2(i,j)))
            !
        Td(kint(k,dk)) = Td(kint(k,dk)) + &
              ProjectPi3(1,1,1,kx,ky) * dimag(u1c(i,j)*conjg(T11(i,j))) + &
              ProjectPi3(2,1,1,kx,ky) * dimag(u2c(i,j)*conjg(T11(i,j))) + &
              ProjectPi3(1,1,2,kx,ky) * dimag(u1c(i,j)*conjg(T12(i,j))) + &
              ProjectPi3(2,1,2,kx,ky) * dimag(u2c(i,j)*conjg(T12(i,j))) + &
              ProjectPi3(1,2,1,kx,ky) * dimag(u1c(i,j)*conjg(T21(i,j))) + &
              ProjectPi3(2,2,1,kx,ky) * dimag(u2c(i,j)*conjg(T21(i,j))) + &
              ProjectPi3(1,2,2,kx,ky) * dimag(u1c(i,j)*conjg(T22(i,j))) + &
              ProjectPi3(2,2,2,kx,ky) * dimag(u2c(i,j)*conjg(T22(i,j)))
            !
        Tdd(kint(k,dk)) = Tdd(kint(k,dk)) + &
              ProjectPi3(1,1,1,kx,ky) * dimag(u1c(i,j)*conjg(T11d(i,j))) + &
              ProjectPi3(2,1,1,kx,ky) * dimag(u2c(i,j)*conjg(T11d(i,j))) + &
              ProjectPi3(1,1,2,kx,ky) * dimag(u1c(i,j)*conjg(T12d(i,j))) + &
              ProjectPi3(2,1,2,kx,ky) * dimag(u2c(i,j)*conjg(T12d(i,j))) + &
              ProjectPi3(1,2,1,kx,ky) * dimag(u1c(i,j)*conjg(T21d(i,j))) + &
              ProjectPi3(2,2,1,kx,ky) * dimag(u2c(i,j)*conjg(T21d(i,j))) + &
              ProjectPi3(1,2,2,kx,ky) * dimag(u1c(i,j)*conjg(T22d(i,j))) + &
              ProjectPi3(2,2,2,kx,ky) * dimag(u2c(i,j)*conjg(T22d(i,j)))
            !
        Tdstheta(kint(k,dk)) = Tdstheta(kint(k,dk)) + &
              ProjectPi2(1,1,kx,ky) * dreal(u1c(i,j)*conjg(Tstheta1(i,j))) + &
              ProjectPi2(1,2,kx,ky) * dreal(u1c(i,j)*conjg(Tstheta2(i,j))) + &
              ProjectPi2(2,1,kx,ky) * dreal(u2c(i,j)*conjg(Tstheta1(i,j))) + &
              ProjectPi2(2,2,kx,ky) * dreal(u2c(i,j)*conjg(Tstheta2(i,j)))
              !
        Tddtheta(kint(k,dk)) = Tddtheta(kint(k,dk)) + &
              ProjectPi2(1,1,kx,ky) * dreal(u1c(i,j)*conjg(Tdtheta1(i,j))) + &
              ProjectPi2(1,2,kx,ky) * dreal(u1c(i,j)*conjg(Tdtheta2(i,j))) + &
              ProjectPi2(2,1,kx,ky) * dreal(u2c(i,j)*conjg(Tdtheta1(i,j))) + &
              ProjectPi2(2,2,kx,ky) * dreal(u2c(i,j)*conjg(Tdtheta2(i,j)))
        Tcount(kint(k,dk)) = Tcount(kint(k,dk))+1
      endif
      Tsall = Tsall + ProjectP3(1,1,1,kx,ky) * dimag(u1s(i,j)*conjg(T11(i,j))) + &
                      ProjectP3(2,1,1,kx,ky) * dimag(u2s(i,j)*conjg(T11(i,j))) + &
                      ProjectP3(1,1,2,kx,ky) * dimag(u1s(i,j)*conjg(T12(i,j))) + &
                      ProjectP3(2,1,2,kx,ky) * dimag(u2s(i,j)*conjg(T12(i,j))) + &
                      ProjectP3(1,2,1,kx,ky) * dimag(u1s(i,j)*conjg(T21(i,j))) + &
                      ProjectP3(2,2,1,kx,ky) * dimag(u2s(i,j)*conjg(T21(i,j))) + &
                      ProjectP3(1,2,2,kx,ky) * dimag(u1s(i,j)*conjg(T22(i,j))) + &
                      ProjectP3(2,2,2,kx,ky) * dimag(u2s(i,j)*conjg(T22(i,j)))
      Tssall=Tssall + ProjectP3(1,1,1,kx,ky) * dimag(u1s(i,j)*conjg(T11s(i,j))) + &
                      ProjectP3(2,1,1,kx,ky) * dimag(u2s(i,j)*conjg(T11s(i,j))) + &
                      ProjectP3(1,1,2,kx,ky) * dimag(u1s(i,j)*conjg(T12s(i,j))) + &
                      ProjectP3(2,1,2,kx,ky) * dimag(u2s(i,j)*conjg(T12s(i,j))) + &
                      ProjectP3(1,2,1,kx,ky) * dimag(u1s(i,j)*conjg(T21s(i,j))) + &
                      ProjectP3(2,2,1,kx,ky) * dimag(u2s(i,j)*conjg(T21s(i,j))) + &
                      ProjectP3(1,2,2,kx,ky) * dimag(u1s(i,j)*conjg(T22s(i,j))) + &
                      ProjectP3(2,2,2,kx,ky) * dimag(u2s(i,j)*conjg(T22s(i,j)))
      Tcall = Tcall + ProjectPi3(1,1,1,kx,ky) * dimag(u1c(i,j)*conjg(T11(i,j))) + &
                      ProjectPi3(2,1,1,kx,ky) * dimag(u2c(i,j)*conjg(T11(i,j))) + &
                      ProjectPi3(1,1,2,kx,ky) * dimag(u1c(i,j)*conjg(T12(i,j))) + &
                      ProjectPi3(2,1,2,kx,ky) * dimag(u2c(i,j)*conjg(T12(i,j))) + &
                      ProjectPi3(1,2,1,kx,ky) * dimag(u1c(i,j)*conjg(T21(i,j))) + &
                      ProjectPi3(2,2,1,kx,ky) * dimag(u2c(i,j)*conjg(T21(i,j))) + &
                      ProjectPi3(1,2,2,kx,ky) * dimag(u1c(i,j)*conjg(T22(i,j))) + &
                      ProjectPi3(2,2,2,kx,ky) * dimag(u2c(i,j)*conjg(T22(i,j)))
      Tccall=Tccall + ProjectPi3(1,1,1,kx,ky) * dimag(u1c(i,j)*conjg(T11d(i,j))) + &
                      ProjectPi3(2,1,1,kx,ky) * dimag(u2c(i,j)*conjg(T11d(i,j))) + &
                      ProjectPi3(1,1,2,kx,ky) * dimag(u1c(i,j)*conjg(T12d(i,j))) + &
                      ProjectPi3(2,1,2,kx,ky) * dimag(u2c(i,j)*conjg(T12d(i,j))) + &
                      ProjectPi3(1,2,1,kx,ky) * dimag(u1c(i,j)*conjg(T21d(i,j))) + &
                      ProjectPi3(2,2,1,kx,ky) * dimag(u2c(i,j)*conjg(T21d(i,j))) + &
                      ProjectPi3(1,2,2,kx,ky) * dimag(u1c(i,j)*conjg(T22d(i,j))) + &
                      ProjectPi3(2,2,2,kx,ky) * dimag(u2c(i,j)*conjg(T22d(i,j)))
      k2Ts = k2Ts + k**2 * ProjectP3(1,1,1,kx,ky) * dimag(u1s(i,j)*conjg(T11(i,j))) + &
              k**2 * ProjectP3(2,1,1,kx,ky) * dimag(u2s(i,j)*conjg(T11(i,j))) + &
              k**2 * ProjectP3(1,1,2,kx,ky) * dimag(u1s(i,j)*conjg(T12(i,j))) + &
              k**2 * ProjectP3(2,1,2,kx,ky) * dimag(u2s(i,j)*conjg(T12(i,j))) + &
              k**2 * ProjectP3(1,2,1,kx,ky) * dimag(u1s(i,j)*conjg(T21(i,j))) + &
              k**2 * ProjectP3(2,2,1,kx,ky) * dimag(u2s(i,j)*conjg(T21(i,j))) + &
              k**2 * ProjectP3(1,2,2,kx,ky) * dimag(u1s(i,j)*conjg(T22(i,j))) + &
              k**2 * ProjectP3(2,2,2,kx,ky) * dimag(u2s(i,j)*conjg(T22(i,j)))
      k2Tss= k2Tss+ k**2 * ProjectP3(1,1,1,kx,ky) * dimag(u1s(i,j)*conjg(T11s(i,j))) + &
              k**2 * ProjectP3(2,1,1,kx,ky) * dimag(u2s(i,j)*conjg(T11s(i,j))) + &
              k**2 * ProjectP3(1,1,2,kx,ky) * dimag(u1s(i,j)*conjg(T12s(i,j))) + &
              k**2 * ProjectP3(2,1,2,kx,ky) * dimag(u2s(i,j)*conjg(T12s(i,j))) + &
              k**2 * ProjectP3(1,2,1,kx,ky) * dimag(u1s(i,j)*conjg(T21s(i,j))) + &
              k**2 * ProjectP3(2,2,1,kx,ky) * dimag(u2s(i,j)*conjg(T21s(i,j))) + &
              k**2 * ProjectP3(1,2,2,kx,ky) * dimag(u1s(i,j)*conjg(T22s(i,j))) + &
              k**2 * ProjectP3(2,2,2,kx,ky) * dimag(u2s(i,j)*conjg(T22s(i,j)))
      k2Tc =  k2Tc + k**2 * ProjectPi3(1,1,1,kx,ky) * dimag(u1c(i,j)*conjg(T11(i,j))) + &
              k**2 * ProjectPi3(2,1,1,kx,ky) * dimag(u2c(i,j)*conjg(T11(i,j))) + &
              k**2 * ProjectPi3(1,1,2,kx,ky) * dimag(u1c(i,j)*conjg(T12(i,j))) + &
              k**2 * ProjectPi3(2,1,2,kx,ky) * dimag(u2c(i,j)*conjg(T12(i,j))) + &
              k**2 * ProjectPi3(1,2,1,kx,ky) * dimag(u1c(i,j)*conjg(T21(i,j))) + &
              k**2 * ProjectPi3(2,2,1,kx,ky) * dimag(u2c(i,j)*conjg(T21(i,j))) + &
              k**2 * ProjectPi3(1,2,2,kx,ky) * dimag(u1c(i,j)*conjg(T22(i,j))) + &
              k**2 * ProjectPi3(2,2,2,kx,ky) * dimag(u2c(i,j)*conjg(T22(i,j)))
      k2Tcc=  k2Tcc+ k**2 * ProjectPi3(1,1,1,kx,ky) * dimag(u1c(i,j)*conjg(T11d(i,j))) + &
              k**2 * ProjectPi3(2,1,1,kx,ky) * dimag(u2c(i,j)*conjg(T11d(i,j))) + &
              k**2 * ProjectPi3(1,1,2,kx,ky) * dimag(u1c(i,j)*conjg(T12d(i,j))) + &
              k**2 * ProjectPi3(2,1,2,kx,ky) * dimag(u2c(i,j)*conjg(T12d(i,j))) + &
              k**2 * ProjectPi3(1,2,1,kx,ky) * dimag(u1c(i,j)*conjg(T21d(i,j))) + &
              k**2 * ProjectPi3(2,2,1,kx,ky) * dimag(u2c(i,j)*conjg(T21d(i,j))) + &
              k**2 * ProjectPi3(1,2,2,kx,ky) * dimag(u1c(i,j)*conjg(T22d(i,j))) + &
              k**2 * ProjectPi3(2,2,2,kx,ky) * dimag(u2c(i,j)*conjg(T22d(i,j)))
    enddo
    enddo
    !
    do i=0,allkmax
      Ts(i) = -psum(Ts(i))
      Tstheta(i) = psum(Tstheta(i))
      Td(i) = -psum(Td(i))
      Tdstheta(i) = psum(Tdstheta(i))
      Tddtheta(i) = psum(Tddtheta(i))
      Tss(i) = -psum(Tss(i))
      Tdd(i) = -psum(Tdd(i))
      Tcount(i) = psum(Tcount(i))
    end do
    !
    !
    Tsall = psum(Tsall)
    Tcall = psum(Tcall)
    k2Ts = psum(k2Ts)
    k2Tc = psum(k2Tc)
    Tssall = psum(Tssall)
    Tccall = psum(Tccall)
    k2Tss = psum(k2Tss)
    k2Tcc = psum(k2Tcc)
    !
    if(mpirank == 0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Tspec'//stepname//'.dat'
      else
        outfilename = 'pp/Tspec.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_a, &
                        firstline='nstep time k Ts Tss Tstheta Td Tdd Tdstheta Tddtheta')
      do i=0,allkmax
        call listwrite(hand_a,dble(i),Ts(i),Tss(i),Tstheta(i),Td(i),Tdd(i),Tdstheta(i),Tddtheta(i))
      end do
      !
      print*,' <<< pp/Tspec'//stepname//'.dat ... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Tspec_aux'//stepname//'.dat'
      else
        outfilename = 'pp/Tspec_aux.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
            firstline='nstep time k2Ts k2Tc Tsall Tcall Tssall Tccall k2Tss k2Tcc')
      call listwrite(hand_b,k2Ts,k2Tc,Tsall,Tcall,Tssall,Tccall,k2Tss,k2Tcc)
      !
      print*,' <<< '//outfilename//'... done.'
    endif
    !
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_pspe)
    call fftw_free(c_theta)
    call fftw_free(c_u1s)
    call fftw_free(c_u2s)
    call fftw_free(c_u1c)
    call fftw_free(c_u2c)
    call fftw_free(c_T11)
    call fftw_free(c_T12)
    call fftw_free(c_T21)
    call fftw_free(c_T22)
    call fftw_free(c_T11s)
    call fftw_free(c_T12s)
    call fftw_free(c_T21s)
    call fftw_free(c_T22s)
    call fftw_free(c_T11d)
    call fftw_free(c_T12d)
    call fftw_free(c_T21d)
    call fftw_free(c_T22d)
    call fftw_free(c_Tstheta1)
    call fftw_free(c_Tstheta2)
    call fftw_free(c_Tdtheta1)
    call fftw_free(c_Tdtheta2)
    call mpistop
    !
  end subroutine instanttriad2D
  !
  subroutine instanttriad3D(thefilenumb)
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,   only : time,nstep,im,jm,km
    use commarray, only : vel, rho, prs
    use hdf5io
    use utility,   only : listinit,listwrite
    use parallel,  only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    integer,intent(in) :: thefilenumb
    character(len=128) :: infilename
    character(len=4) :: stepname
    real(8) :: u1mean,u2mean,u3mean,rhomean,prsmean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: u1spe,u2spe,u3spe,pspe
    complex(8), allocatable, dimension(:,:,:) :: ucspe
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: theta,u1c,u2c,u3c,u1s,u2s,u3s
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: T11,T12,T13,T21,T22,T23,T31,T32,T33
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: Tstheta1,Tstheta2,Tstheta3,Tdtheta1,Tdtheta2,Tdtheta3
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: T11s,T12s,T13s,T21s,T22s,T23s,T31s,T32s,T33s
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: T11d,T12d,T13d,T21d,T22d,T23d,T31d,T32d,T33d
    real(8), allocatable, dimension(:) :: Ts,Tstheta,Td,Tddtheta,Tdstheta,Tcount,Tss,Tdd
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3,u1,u2,u3
    integer :: allkmax
    real(8) :: kk,dk,p,q,kx,ky,kz !wave number
    character(len=128) :: outfilename
    integer :: hand_a,hand_b
    character(len=1) :: modeio
    integer :: i,j,n,s,t,k
    real(8) :: k2Ts, k2Tc, Tsall, Tcall, Tssall, Tccall, k2Tss, k2Tcc
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_u3spe, c_pspe
    type(C_PTR) :: c_theta,c_u1c,c_u2c,c_u3c,c_u1s,c_u2s,c_u3s
    type(C_PTR) :: c_T11,c_T12,c_T13,c_T21,c_T22,c_T23,c_T31,c_T32,c_T33
    type(C_PTR) :: c_Tstheta1,c_Tstheta2,c_Tstheta3,c_Tdtheta1,c_Tdtheta2,c_Tdtheta3
    type(C_PTR) :: c_T11s,c_T12s,c_T13s,c_T21s,c_T22s,c_T23s,c_T31s,c_T32s,c_T33s
    type(C_PTR) :: c_T11d,c_T12d,c_T13d,c_T21d,c_T22d,c_T23d,c_T31d,c_T32d,c_T33d
    call readinput
    !
    modeio='h'
    dk = 0.5d0
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    allkmax=ceiling(1.d0*min(ia,ja,ka))
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka,"allkmax:",allkmax
    if(ka==0) stop 'Please use instantspectra2D'
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km), prs(0:im,0:jm,0:km))
    !
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='p',  var=prs(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    allocate(u1(1:im,1:jm,1:km),u2(1:im,1:jm,1:km),u3(1:im,1:jm,1:km))
    ! Calculate favre average
    u1mean = 0.0d0
    u2mean = 0.0d0
    u3mean = 0.0d0
    rhomean = 0.0d0
    prsmean = 0.0d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      u1mean = u1mean + vel(i,j,k,1)
      u2mean = u2mean + vel(i,j,k,2)
      u3mean = u3mean + vel(i,j,k,3)
      rhomean = rhomean + rho(i,j,k)
      prsmean = prsmean + prs(i,j,k)
    enddo
    enddo
    enddo
    !
    rhomean = psum(rhomean) / (1.0d0*ia*ja*ka)
    u1mean = psum(u1mean) / (1.d0*ia*ja*ka)
    u2mean = psum(u2mean) / (1.d0*ia*ja*ka)
    u3mean = psum(u3mean) / (1.d0*ia*ja*ka)
    prsmean = psum(prsmean) / (1.d0*ia*ja*ka)
    if(mpirank==0) print *, 'u1mean=',u1mean, 'u2mean=',u2mean, 'u3mean=',u3mean, 'prsmean=',prsmean
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw,kmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw,kmfftw])
    c_u3spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3spe, u3spe, [imfftw,jmfftw,kmfftw])
    c_pspe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_pspe, pspe, [imfftw,jmfftw,kmfftw])
    !
    forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    !!!! Convert to complex
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1(i,j,k) = vel(i,j,k,1)-u1mean
      u2(i,j,k) = vel(i,j,k,2)-u2mean
      u3(i,j,k) = vel(i,j,k,3)-u3mean
      u1spe(i,j,k)=CMPLX(u1(i,j,k),0.d0,C_INTPTR_T)
      u2spe(i,j,k)=CMPLX(u2(i,j,k),0.d0,C_INTPTR_T)
      u3spe(i,j,k)=CMPLX(u3(i,j,k),0.d0,C_INTPTR_T)
      pspe(i,j,k)=CMPLX(prs(i,j,k)-prsmean,0.d0,C_INTPTR_T)
      !
    end do
    end do
    end do
    !
    !!!! Do FFT
    !
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    call fftw_mpi_execute_dft(forward_plan,u3spe,u3spe)
    call fftw_mpi_execute_dft(forward_plan,pspe,pspe)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j,k)=u1spe(i,j,k)/(1.d0*ia*ja*ka)
      u2spe(i,j,k)=u2spe(i,j,k)/(1.d0*ia*ja*ka)
      u3spe(i,j,k)=u3spe(i,j,k)/(1.d0*ia*ja*ka)
      pspe(i,j,k)=pspe(i,j,k)/(1.d0*ia*ja*ka)
      !
    end do
    end do
    end do
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j,k) = real(i-1,8)
      else if(i<=ia) then
        k1(i,j,k) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if(j <= (ja/2+1)) then
        k2(i,j,k) = real(j-1,8)
      else if(j<=ja) then
        k2(i,j,k) = real(j-ja-1,8)
      else
        print *,"Error, no wave number possible, j must smaller than ja-1 !"
      end if
      !
      if((k+k0) <= (ka/2+1)) then
        k3(i,j,k) = real(k+k0-1,8)
      else if((k+k0)<=ka) then
        k3(i,j,k) = real(k+k0-ka-1,8)
      else
        print *,"Error, no wave number possible, (k+k0) must smaller than ka-1 !"
      end if
      !
    enddo
    enddo
    enddo
    !
    allocate(ucspe(1:im,1:jm,1:km))
    !
    c_theta = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_theta, theta, [imfftw,jmfftw,kmfftw])
    c_u1s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1s, u1s, [imfftw,jmfftw,kmfftw])
    c_u2s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2s, u2s, [imfftw,jmfftw,kmfftw])
    c_u3s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3s, u3s, [imfftw,jmfftw,kmfftw])
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw,kmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw,kmfftw])
    c_u3c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3c, u3c, [imfftw,jmfftw,kmfftw])
    !
    !!!! Do S-C decomposition
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      kk=dsqrt(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2+1.d-15)
      !
      ucspe(i,j,k) = k1(i,j,k)/kk * u1spe(i,j,k) + k2(i,j,k)/kk * u2spe(i,j,k) + k3(i,j,k)/kk * u3spe(i,j,k)
      u1c(i,j,k)=  k1(i,j,k)*k1(i,j,k)/(kk**2) * u1spe(i,j,k) + k1(i,j,k)*k2(i,j,k)/(kk**2) * u2spe(i,j,k) &
                + k1(i,j,k)*k3(i,j,k)/(kk**2) * u3spe(i,j,k)
      u2c(i,j,k)=  k2(i,j,k)*k1(i,j,k)/(kk**2) * u1spe(i,j,k) + k2(i,j,k)*k2(i,j,k)/(kk**2) * u2spe(i,j,k) &
                + k2(i,j,k)*k3(i,j,k)/(kk**2) * u3spe(i,j,k)
      u3c(i,j,k)=  k3(i,j,k)*k1(i,j,k)/(kk**2) * u1spe(i,j,k) + k3(i,j,k)*k2(i,j,k)/(kk**2) * u2spe(i,j,k) &
                + k3(i,j,k)*k3(i,j,k)/(kk**2) * u3spe(i,j,k)
      u1s(i,j,k)=  u1spe(i,j,k) - u1c(i,j,k)
      u2s(i,j,k)=  u2spe(i,j,k) - u2c(i,j,k)
      u3s(i,j,k)=  u3spe(i,j,k) - u3c(i,j,k)
      !
      theta(i,j,k)= CMPLX(0.d0,1.d0,C_INTPTR_T) * (ucspe(i,j,k) * kk)
      !
    end do
    end do
    enddo
    !
    !!!! Do inverse FFT
    !
    call fftw_mpi_execute_dft(backward_plan,theta,theta)
    call fftw_mpi_execute_dft(backward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(backward_plan,u2s,u2s)
    call fftw_mpi_execute_dft(backward_plan,u3s,u3s)
    call fftw_mpi_execute_dft(backward_plan,u1c,u1c)
    call fftw_mpi_execute_dft(backward_plan,u2c,u2c)
    call fftw_mpi_execute_dft(backward_plan,u3c,u3c)
    !
    c_T11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11,T11, [imfftw,jmfftw,kmfftw])
    c_T12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12,T12, [imfftw,jmfftw,kmfftw])
    c_T13 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T13,T13, [imfftw,jmfftw,kmfftw])
    c_T21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21,T21, [imfftw,jmfftw,kmfftw])
    c_T22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22,T22, [imfftw,jmfftw,kmfftw])
    c_T23 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T23,T23, [imfftw,jmfftw,kmfftw])
    c_T31 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T31,T31, [imfftw,jmfftw,kmfftw])
    c_T32 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T32,T32, [imfftw,jmfftw,kmfftw])
    c_T33 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T33,T33, [imfftw,jmfftw,kmfftw])
    !
    c_T11s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11s,T11s, [imfftw,jmfftw,kmfftw])
    c_T12s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12s,T12s, [imfftw,jmfftw,kmfftw])
    c_T13s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T13s,T13s, [imfftw,jmfftw,kmfftw])
    c_T21s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21s,T21s, [imfftw,jmfftw,kmfftw])
    c_T22s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22s,T22s, [imfftw,jmfftw,kmfftw])
    c_T23s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T23s,T23s, [imfftw,jmfftw,kmfftw])
    c_T31s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T31s,T31s, [imfftw,jmfftw,kmfftw])
    c_T32s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T32s,T32s, [imfftw,jmfftw,kmfftw])
    c_T33s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T33s,T33s, [imfftw,jmfftw,kmfftw])
    !
    c_T11d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11d,T11d, [imfftw,jmfftw,kmfftw])
    c_T12d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12d,T12d, [imfftw,jmfftw,kmfftw])
    c_T13d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T13d,T13d, [imfftw,jmfftw,kmfftw])
    c_T21d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21d,T21d, [imfftw,jmfftw,kmfftw])
    c_T22d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22d,T22d, [imfftw,jmfftw,kmfftw])
    c_T23d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T23d,T23d, [imfftw,jmfftw,kmfftw])
    c_T31d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T31d,T31d, [imfftw,jmfftw,kmfftw])
    c_T32d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T32d,T32d, [imfftw,jmfftw,kmfftw])
    c_T33d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T33d,T33d, [imfftw,jmfftw,kmfftw])
    !
    c_Tstheta1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tstheta1,Tstheta1, [imfftw,jmfftw,kmfftw])
    c_Tstheta2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tstheta2,Tstheta2, [imfftw,jmfftw,kmfftw])
    c_Tstheta3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tstheta3,Tstheta3, [imfftw,jmfftw,kmfftw])
    c_Tdtheta1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tdtheta1,Tdtheta1, [imfftw,jmfftw,kmfftw])
    c_Tdtheta2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tdtheta2,Tdtheta2, [imfftw,jmfftw,kmfftw])
    c_Tdtheta3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tdtheta3,Tdtheta3, [imfftw,jmfftw,kmfftw])
    !
    ! This part is in physical space
    do k=1,km
    do j=1,jm
    do i=1,im
      T11(i,j,k) = u1(i,j,k)*u1(i,j,k)
      T12(i,j,k) = u1(i,j,k)*u2(i,j,k)
      T13(i,j,k) = u1(i,j,k)*u3(i,j,k)
      T21(i,j,k) = u2(i,j,k)*u1(i,j,k)
      T22(i,j,k) = u2(i,j,k)*u2(i,j,k)
      T23(i,j,k) = u2(i,j,k)*u3(i,j,k)
      T31(i,j,k) = u3(i,j,k)*u1(i,j,k)
      T32(i,j,k) = u3(i,j,k)*u2(i,j,k)
      T33(i,j,k) = u3(i,j,k)*u3(i,j,k)
      !
      T11s(i,j,k) = u1s(i,j,k)*u1s(i,j,k)
      T12s(i,j,k) = u1s(i,j,k)*u2s(i,j,k)
      T13s(i,j,k) = u1s(i,j,k)*u3s(i,j,k)
      T21s(i,j,k) = u2s(i,j,k)*u1s(i,j,k)
      T22s(i,j,k) = u2s(i,j,k)*u2s(i,j,k)
      T23s(i,j,k) = u2s(i,j,k)*u3s(i,j,k)
      T31s(i,j,k) = u3s(i,j,k)*u1s(i,j,k)
      T32s(i,j,k) = u3s(i,j,k)*u2s(i,j,k)
      T33s(i,j,k) = u3s(i,j,k)*u3s(i,j,k)
      !
      T11d(i,j,k) = u1c(i,j,k)*u1c(i,j,k)
      T12d(i,j,k) = u1c(i,j,k)*u2c(i,j,k)
      T13d(i,j,k) = u1c(i,j,k)*u3c(i,j,k)
      T21d(i,j,k) = u2c(i,j,k)*u1c(i,j,k)
      T22d(i,j,k) = u2c(i,j,k)*u2c(i,j,k)
      T23d(i,j,k) = u2c(i,j,k)*u3c(i,j,k)
      T31d(i,j,k) = u3c(i,j,k)*u1c(i,j,k)
      T32d(i,j,k) = u3c(i,j,k)*u2c(i,j,k)
      T33d(i,j,k) = u3c(i,j,k)*u3c(i,j,k)
      !
      Tstheta1(i,j,k) = theta(i,j,k)*u1s(i,j,k)
      Tstheta2(i,j,k) = theta(i,j,k)*u2s(i,j,k)
      Tstheta3(i,j,k) = theta(i,j,k)*u3s(i,j,k)
      Tdtheta1(i,j,k) = theta(i,j,k)*u1c(i,j,k)
      Tdtheta2(i,j,k) = theta(i,j,k)*u2c(i,j,k)
      Tdtheta3(i,j,k) = theta(i,j,k)*u3c(i,j,k)
    enddo
    enddo
    enddo
    !
    call fftw_mpi_execute_dft(forward_plan,T11,T11)
    call fftw_mpi_execute_dft(forward_plan,T12,T12)
    call fftw_mpi_execute_dft(forward_plan,T13,T13)
    call fftw_mpi_execute_dft(forward_plan,T21,T21)
    call fftw_mpi_execute_dft(forward_plan,T22,T22)
    call fftw_mpi_execute_dft(forward_plan,T23,T23)
    call fftw_mpi_execute_dft(forward_plan,T31,T31)
    call fftw_mpi_execute_dft(forward_plan,T32,T32)
    call fftw_mpi_execute_dft(forward_plan,T33,T33)
    !
    call fftw_mpi_execute_dft(forward_plan,T11s,T11s)
    call fftw_mpi_execute_dft(forward_plan,T12s,T12s)
    call fftw_mpi_execute_dft(forward_plan,T13s,T13s)
    call fftw_mpi_execute_dft(forward_plan,T21s,T21s)
    call fftw_mpi_execute_dft(forward_plan,T22s,T22s)
    call fftw_mpi_execute_dft(forward_plan,T23s,T23s)
    call fftw_mpi_execute_dft(forward_plan,T31s,T31s)
    call fftw_mpi_execute_dft(forward_plan,T32s,T32s)
    call fftw_mpi_execute_dft(forward_plan,T33s,T33s)
    !
    call fftw_mpi_execute_dft(forward_plan,T11d,T11d)
    call fftw_mpi_execute_dft(forward_plan,T12d,T12d)
    call fftw_mpi_execute_dft(forward_plan,T13d,T13d)
    call fftw_mpi_execute_dft(forward_plan,T21d,T21d)
    call fftw_mpi_execute_dft(forward_plan,T22d,T22d)
    call fftw_mpi_execute_dft(forward_plan,T23d,T23d)
    call fftw_mpi_execute_dft(forward_plan,T31d,T31d)
    call fftw_mpi_execute_dft(forward_plan,T32d,T32d)
    call fftw_mpi_execute_dft(forward_plan,T33d,T33d)
    !
    call fftw_mpi_execute_dft(forward_plan,Tstheta1,Tstheta1)
    call fftw_mpi_execute_dft(forward_plan,Tstheta2,Tstheta2)
    call fftw_mpi_execute_dft(forward_plan,Tstheta3,Tstheta3)
    call fftw_mpi_execute_dft(forward_plan,Tdtheta1,Tdtheta1)
    call fftw_mpi_execute_dft(forward_plan,Tdtheta2,Tdtheta2)
    call fftw_mpi_execute_dft(forward_plan,Tdtheta3,Tdtheta3)
    !
    call fftw_mpi_execute_dft(forward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(forward_plan,u2s,u2s)
    call fftw_mpi_execute_dft(forward_plan,u3s,u3s)
    call fftw_mpi_execute_dft(forward_plan,u1c,u1c)
    call fftw_mpi_execute_dft(forward_plan,u2c,u2c)
    call fftw_mpi_execute_dft(forward_plan,u3c,u3c)
    !
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      T11(i,j,k)=T11(i,j,k)/(1.d0*ia*ja*ka)
      T12(i,j,k)=T12(i,j,k)/(1.d0*ia*ja*ka)
      T13(i,j,k)=T13(i,j,k)/(1.d0*ia*ja*ka)
      T21(i,j,k)=T21(i,j,k)/(1.d0*ia*ja*ka)
      T22(i,j,k)=T22(i,j,k)/(1.d0*ia*ja*ka)
      T23(i,j,k)=T23(i,j,k)/(1.d0*ia*ja*ka)
      T31(i,j,k)=T31(i,j,k)/(1.d0*ia*ja*ka)
      T32(i,j,k)=T32(i,j,k)/(1.d0*ia*ja*ka)
      T33(i,j,k)=T33(i,j,k)/(1.d0*ia*ja*ka)
      !
      T11s(i,j,k)=T11s(i,j,k)/(1.d0*ia*ja*ka)
      T12s(i,j,k)=T12s(i,j,k)/(1.d0*ia*ja*ka)
      T13s(i,j,k)=T13s(i,j,k)/(1.d0*ia*ja*ka)
      T21s(i,j,k)=T21s(i,j,k)/(1.d0*ia*ja*ka)
      T22s(i,j,k)=T22s(i,j,k)/(1.d0*ia*ja*ka)
      T23s(i,j,k)=T23s(i,j,k)/(1.d0*ia*ja*ka)
      T31s(i,j,k)=T31s(i,j,k)/(1.d0*ia*ja*ka)
      T32s(i,j,k)=T32s(i,j,k)/(1.d0*ia*ja*ka)
      T33s(i,j,k)=T33s(i,j,k)/(1.d0*ia*ja*ka)
      !
      T11d(i,j,k)=T11d(i,j,k)/(1.d0*ia*ja*ka)
      T12d(i,j,k)=T12d(i,j,k)/(1.d0*ia*ja*ka)
      T13d(i,j,k)=T13d(i,j,k)/(1.d0*ia*ja*ka)
      T21d(i,j,k)=T21d(i,j,k)/(1.d0*ia*ja*ka)
      T22d(i,j,k)=T22d(i,j,k)/(1.d0*ia*ja*ka)
      T23d(i,j,k)=T23d(i,j,k)/(1.d0*ia*ja*ka)
      T31d(i,j,k)=T31d(i,j,k)/(1.d0*ia*ja*ka)
      T32d(i,j,k)=T32d(i,j,k)/(1.d0*ia*ja*ka)
      T33d(i,j,k)=T33d(i,j,k)/(1.d0*ia*ja*ka)
      !
      Tstheta1(i,j,k)=Tstheta1(i,j,k)/(1.d0*ia*ja*ka)
      Tstheta2(i,j,k)=Tstheta2(i,j,k)/(1.d0*ia*ja*ka)
      Tstheta3(i,j,k)=Tstheta3(i,j,k)/(1.d0*ia*ja*ka)
      Tdtheta1(i,j,k)=Tdtheta1(i,j,k)/(1.d0*ia*ja*ka)
      Tdtheta2(i,j,k)=Tdtheta2(i,j,k)/(1.d0*ia*ja*ka)
      Tdtheta3(i,j,k)=Tdtheta3(i,j,k)/(1.d0*ia*ja*ka)
      !
      u1s(i,j,k)=u1s(i,j,k)/(1.d0*ia*ja*ka)
      u2s(i,j,k)=u2s(i,j,k)/(1.d0*ia*ja*ka)
      u3s(i,j,k)=u3s(i,j,k)/(1.d0*ia*ja*ka)
      u1c(i,j,k)=u1c(i,j,k)/(1.d0*ia*ja*ka)
      u2c(i,j,k)=u2c(i,j,k)/(1.d0*ia*ja*ka)
      u3c(i,j,k)=u3c(i,j,k)/(1.d0*ia*ja*ka)
      !
    enddo
    enddo
    enddo
    !
    allocate(Ts(0:allkmax),Tstheta(0:allkmax),Td(0:allkmax),Tdstheta(0:allkmax))
    allocate(Tddtheta(0:allkmax),Tcount(0:allkmax),Tss(0:allkmax),Tdd(0:allkmax))
    !
    do i=0,allkmax
      Ts(i)= 0.0d0
      Tstheta(i)= 0.0d0
      Td(i)= 0.0d0
      Tddtheta(i)= 0.0d0
      Tdstheta(i)= 0.0d0
      Tcount(i)=0
      Tss(i)=0.d0
      Tdd(i)=0.d0
    end do
    !
    k2Ts = 0.d0
    k2Tc = 0.d0
    Tsall = 0.d0
    Tcall = 0.d0
    Tssall = 0.d0
    Tccall = 0.d0
    k2Tss = 0.d0
    k2Tcc = 0.d0
    !
    do k = 1,km
    do j = 1,jm
    do i = 1,im
      kx = k1(i,j,k)
      ky = k2(i,j,k)
      kz = k3(i,j,k)
      kk=dsqrt(kx**2+ky**2+kz**2+1.d-15)
      !
      !
      if(kint(kk,dk)<= allkmax)then
        Ts(kint(kk,dk)) = Ts(kint(kk,dk)) + &
                          ProjectP3(1,1,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T11(i,j,k))) + &
                          ProjectP3(2,1,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T11(i,j,k))) + &
                          ProjectP3(3,1,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T11(i,j,k))) + &
                          ProjectP3(1,1,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T12(i,j,k))) + &
                          ProjectP3(2,1,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T12(i,j,k))) + &
                          ProjectP3(3,1,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T12(i,j,k))) + &
                          ProjectP3(1,1,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T13(i,j,k))) + &
                          ProjectP3(2,1,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T13(i,j,k))) + &
                          ProjectP3(3,1,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T13(i,j,k))) + &
                          ProjectP3(1,2,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T21(i,j,k))) + &
                          ProjectP3(2,2,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T21(i,j,k))) + &
                          ProjectP3(3,2,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T21(i,j,k))) + &
                          ProjectP3(1,2,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T22(i,j,k))) + &
                          ProjectP3(2,2,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T22(i,j,k))) + &
                          ProjectP3(3,2,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T22(i,j,k))) + &
                          ProjectP3(1,2,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T23(i,j,k))) + &
                          ProjectP3(2,2,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T23(i,j,k))) + &
                          ProjectP3(3,2,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T23(i,j,k))) + &
                          ProjectP3(1,3,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T31(i,j,k))) + &
                          ProjectP3(2,3,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T31(i,j,k))) + &
                          ProjectP3(3,3,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T31(i,j,k))) + &
                          ProjectP3(1,3,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T32(i,j,k))) + &
                          ProjectP3(2,3,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T32(i,j,k))) + &
                          ProjectP3(3,3,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T32(i,j,k))) + &
                          ProjectP3(1,3,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T33(i,j,k))) + &
                          ProjectP3(2,3,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T33(i,j,k))) + &
                          ProjectP3(3,3,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T33(i,j,k)))
            !
        Tss(kint(kk,dk)) = Tss(kint(kk,dk)) + &
                            ProjectP3(1,1,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T11s(i,j,k))) + &
                            ProjectP3(2,1,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T11s(i,j,k))) + &
                            ProjectP3(3,1,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T11s(i,j,k))) + &
                            ProjectP3(1,1,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T12s(i,j,k))) + &
                            ProjectP3(2,1,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T12s(i,j,k))) + &
                            ProjectP3(3,1,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T12s(i,j,k))) + &
                            ProjectP3(1,1,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T13s(i,j,k))) + &
                            ProjectP3(2,1,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T13s(i,j,k))) + &
                            ProjectP3(3,1,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T13s(i,j,k))) + &
                            ProjectP3(1,2,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T21s(i,j,k))) + &
                            ProjectP3(2,2,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T21s(i,j,k))) + &
                            ProjectP3(3,2,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T21s(i,j,k))) + &
                            ProjectP3(1,2,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T22s(i,j,k))) + &
                            ProjectP3(2,2,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T22s(i,j,k))) + &
                            ProjectP3(3,2,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T22s(i,j,k))) + &
                            ProjectP3(1,2,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T23s(i,j,k))) + &
                            ProjectP3(2,2,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T23s(i,j,k))) + &
                            ProjectP3(3,2,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T23s(i,j,k))) + &
                            ProjectP3(1,3,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T31s(i,j,k))) + &
                            ProjectP3(2,3,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T31s(i,j,k))) + &
                            ProjectP3(3,3,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T31s(i,j,k))) + &
                            ProjectP3(1,3,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T32s(i,j,k))) + &
                            ProjectP3(2,3,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T32s(i,j,k))) + &
                            ProjectP3(3,3,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T32s(i,j,k))) + &
                            ProjectP3(1,3,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T33s(i,j,k))) + &
                            ProjectP3(2,3,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T33s(i,j,k))) + &
                            ProjectP3(3,3,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T33s(i,j,k)))
              !
        Tstheta(kint(kk,dk)) = Tstheta(kint(kk,dk)) + &
                                ProjectP2(1,1,kx,ky,kz) * dreal(u1s(i,j,k)*conjg(Tstheta1(i,j,k)+Tdtheta1(i,j,k))) + &
                                ProjectP2(1,2,kx,ky,kz) * dreal(u1s(i,j,k)*conjg(Tstheta2(i,j,k)+Tdtheta2(i,j,k))) + &
                                ProjectP2(1,3,kx,ky,kz) * dreal(u1s(i,j,k)*conjg(Tstheta3(i,j,k)+Tdtheta3(i,j,k))) + &
                                ProjectP2(2,1,kx,ky,kz) * dreal(u2s(i,j,k)*conjg(Tstheta1(i,j,k)+Tdtheta1(i,j,k))) + &
                                ProjectP2(2,2,kx,ky,kz) * dreal(u2s(i,j,k)*conjg(Tstheta2(i,j,k)+Tdtheta2(i,j,k))) + &
                                ProjectP2(2,3,kx,ky,kz) * dreal(u2s(i,j,k)*conjg(Tstheta3(i,j,k)+Tdtheta3(i,j,k))) + &
                                ProjectP2(3,1,kx,ky,kz) * dreal(u3s(i,j,k)*conjg(Tstheta1(i,j,k)+Tdtheta1(i,j,k))) + &
                                ProjectP2(3,2,kx,ky,kz) * dreal(u3s(i,j,k)*conjg(Tstheta2(i,j,k)+Tdtheta2(i,j,k))) + &
                                ProjectP2(3,3,kx,ky,kz) * dreal(u3s(i,j,k)*conjg(Tstheta3(i,j,k)+Tdtheta3(i,j,k)))
            !
        Td(kint(kk,dk)) = Td(kint(kk,dk)) + &
                          ProjectPi3(1,1,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T11(i,j,k))) + &
                          ProjectPi3(2,1,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T11(i,j,k))) + &
                          ProjectPi3(3,1,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T11(i,j,k))) + &
                          ProjectPi3(1,1,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T12(i,j,k))) + &
                          ProjectPi3(2,1,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T12(i,j,k))) + &
                          ProjectPi3(3,1,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T12(i,j,k))) + &
                          ProjectPi3(1,1,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T13(i,j,k))) + &
                          ProjectPi3(2,1,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T13(i,j,k))) + &
                          ProjectPi3(3,1,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T13(i,j,k))) + &
                          ProjectPi3(1,2,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T21(i,j,k))) + &
                          ProjectPi3(2,2,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T21(i,j,k))) + &
                          ProjectPi3(3,2,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T21(i,j,k))) + &
                          ProjectPi3(1,2,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T22(i,j,k))) + &
                          ProjectPi3(2,2,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T22(i,j,k))) + &
                          ProjectPi3(3,2,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T22(i,j,k))) + &
                          ProjectPi3(1,2,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T23(i,j,k))) + &
                          ProjectPi3(2,2,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T23(i,j,k))) + &
                          ProjectPi3(3,2,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T23(i,j,k))) + &
                          ProjectPi3(1,3,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T31(i,j,k))) + &
                          ProjectPi3(2,3,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T31(i,j,k))) + &
                          ProjectPi3(3,3,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T31(i,j,k))) + &
                          ProjectPi3(1,3,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T32(i,j,k))) + &
                          ProjectPi3(2,3,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T32(i,j,k))) + &
                          ProjectPi3(3,3,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T32(i,j,k))) + &
                          ProjectPi3(1,3,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T33(i,j,k))) + &
                          ProjectPi3(2,3,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T33(i,j,k))) + &
                          ProjectPi3(3,3,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T33(i,j,k)))
          !
        Tdd(kint(kk,dk)) = Tdd(kint(kk,dk)) + &
                            ProjectPi3(1,1,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T11d(i,j,k))) + &
                            ProjectPi3(2,1,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T11d(i,j,k))) + &
                            ProjectPi3(3,1,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T11d(i,j,k))) + &
                            ProjectPi3(1,1,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T12d(i,j,k))) + &
                            ProjectPi3(2,1,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T12d(i,j,k))) + &
                            ProjectPi3(3,1,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T12d(i,j,k))) + &
                            ProjectPi3(1,1,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T13d(i,j,k))) + &
                            ProjectPi3(2,1,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T13d(i,j,k))) + &
                            ProjectPi3(3,1,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T13d(i,j,k))) + &
                            ProjectPi3(1,2,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T21d(i,j,k))) + &
                            ProjectPi3(2,2,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T21d(i,j,k))) + &
                            ProjectPi3(3,2,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T21d(i,j,k))) + &
                            ProjectPi3(1,2,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T22d(i,j,k))) + &
                            ProjectPi3(2,2,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T22d(i,j,k))) + &
                            ProjectPi3(3,2,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T22d(i,j,k))) + &
                            ProjectPi3(1,2,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T23d(i,j,k))) + &
                            ProjectPi3(2,2,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T23d(i,j,k))) + &
                            ProjectPi3(3,2,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T23d(i,j,k))) + &
                            ProjectPi3(1,3,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T31d(i,j,k))) + &
                            ProjectPi3(2,3,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T31d(i,j,k))) + &
                            ProjectPi3(3,3,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T31d(i,j,k))) + &
                            ProjectPi3(1,3,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T32d(i,j,k))) + &
                            ProjectPi3(2,3,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T32d(i,j,k))) + &
                            ProjectPi3(3,3,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T32d(i,j,k))) + &
                            ProjectPi3(1,3,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T33d(i,j,k))) + &
                            ProjectPi3(2,3,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T33d(i,j,k))) + &
                            ProjectPi3(3,3,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T33d(i,j,k)))
            !
        Tdstheta(kint(kk,dk)) = Tdstheta(kint(kk,dk)) + &
                                ProjectPi2(1,1,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tstheta1(i,j,k))) + &
                                ProjectPi2(1,2,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tstheta2(i,j,k))) + &
                                ProjectPi2(1,3,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tstheta3(i,j,k))) + &
                                ProjectPi2(2,1,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tstheta1(i,j,k))) + &
                                ProjectPi2(2,2,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tstheta2(i,j,k))) + &
                                ProjectPi2(2,3,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tstheta3(i,j,k))) + &
                                ProjectPi2(3,1,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tstheta1(i,j,k))) + &
                                ProjectPi2(3,2,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tstheta2(i,j,k))) + &
                                ProjectPi2(3,3,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tstheta3(i,j,k)))
                                !
        Tddtheta(kint(kk,dk)) = Tddtheta(kint(kk,dk)) + &
                                ProjectPi2(1,1,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tdtheta1(i,j,k))) + &
                                ProjectPi2(1,2,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tdtheta2(i,j,k))) + &
                                ProjectPi2(1,3,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tdtheta3(i,j,k))) + &
                                ProjectPi2(2,1,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tdtheta1(i,j,k))) + &
                                ProjectPi2(2,2,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tdtheta2(i,j,k))) + &
                                ProjectPi2(2,3,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tdtheta3(i,j,k))) + &
                                ProjectPi2(3,1,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tdtheta1(i,j,k))) + &
                                ProjectPi2(3,2,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tdtheta2(i,j,k))) + &
                                ProjectPi2(3,3,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tdtheta3(i,j,k)))
        Tcount(kint(kk,dk)) = Tcount(kint(kk,dk))+1
      endif
      Tsall = Tsall + ProjectP3(1,1,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T11(i,j,k))) + &
                      ProjectP3(2,1,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T11(i,j,k))) + &
                      ProjectP3(3,1,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T11(i,j,k))) + &
                      ProjectP3(1,1,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T12(i,j,k))) + &
                      ProjectP3(2,1,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T12(i,j,k))) + &
                      ProjectP3(3,1,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T12(i,j,k))) + &
                      ProjectP3(1,1,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T13(i,j,k))) + &
                      ProjectP3(2,1,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T13(i,j,k))) + &
                      ProjectP3(3,1,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T13(i,j,k))) + &
                      ProjectP3(1,2,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T21(i,j,k))) + &
                      ProjectP3(2,2,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T21(i,j,k))) + &
                      ProjectP3(3,2,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T21(i,j,k))) + &
                      ProjectP3(1,2,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T22(i,j,k))) + &
                      ProjectP3(2,2,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T22(i,j,k))) + &
                      ProjectP3(3,2,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T22(i,j,k))) + &
                      ProjectP3(1,2,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T23(i,j,k))) + &
                      ProjectP3(2,2,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T23(i,j,k))) + &
                      ProjectP3(3,2,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T23(i,j,k))) + &
                      ProjectP3(1,3,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T31(i,j,k))) + &
                      ProjectP3(2,3,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T31(i,j,k))) + &
                      ProjectP3(3,3,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T31(i,j,k))) + &
                      ProjectP3(1,3,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T32(i,j,k))) + &
                      ProjectP3(2,3,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T32(i,j,k))) + &
                      ProjectP3(3,3,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T32(i,j,k))) + &
                      ProjectP3(1,3,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T33(i,j,k))) + &
                      ProjectP3(2,3,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T33(i,j,k))) + &
                      ProjectP3(3,3,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T33(i,j,k)))
      !
      Tssall=Tssall + ProjectP3(1,1,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T11s(i,j,k))) + &
                      ProjectP3(2,1,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T11s(i,j,k))) + &
                      ProjectP3(3,1,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T11s(i,j,k))) + &
                      ProjectP3(1,1,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T12s(i,j,k))) + &
                      ProjectP3(2,1,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T12s(i,j,k))) + &
                      ProjectP3(3,1,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T12s(i,j,k))) + &
                      ProjectP3(1,1,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T13s(i,j,k))) + &
                      ProjectP3(2,1,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T13s(i,j,k))) + &
                      ProjectP3(3,1,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T13s(i,j,k))) + &
                      ProjectP3(1,2,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T21s(i,j,k))) + &
                      ProjectP3(2,2,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T21s(i,j,k))) + &
                      ProjectP3(3,2,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T21s(i,j,k))) + &
                      ProjectP3(1,2,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T22s(i,j,k))) + &
                      ProjectP3(2,2,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T22s(i,j,k))) + &
                      ProjectP3(3,2,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T22s(i,j,k))) + &
                      ProjectP3(1,2,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T23s(i,j,k))) + &
                      ProjectP3(2,2,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T23s(i,j,k))) + &
                      ProjectP3(3,2,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T23s(i,j,k))) + &
                      ProjectP3(1,3,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T31s(i,j,k))) + &
                      ProjectP3(2,3,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T31s(i,j,k))) + &
                      ProjectP3(3,3,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T31s(i,j,k))) + &
                      ProjectP3(1,3,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T32s(i,j,k))) + &
                      ProjectP3(2,3,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T32s(i,j,k))) + &
                      ProjectP3(3,3,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T32s(i,j,k))) + &
                      ProjectP3(1,3,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T33s(i,j,k))) + &
                      ProjectP3(2,3,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T33s(i,j,k))) + &
                      ProjectP3(3,3,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T33s(i,j,k)))
      !
      Tcall = Tcall + ProjectPi3(1,1,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T11(i,j,k))) + &
                      ProjectPi3(2,1,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T11(i,j,k))) + &
                      ProjectPi3(3,1,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T11(i,j,k))) + &
                      ProjectPi3(1,1,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T12(i,j,k))) + &
                      ProjectPi3(2,1,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T12(i,j,k))) + &
                      ProjectPi3(3,1,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T12(i,j,k))) + &
                      ProjectPi3(1,1,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T13(i,j,k))) + &
                      ProjectPi3(2,1,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T13(i,j,k))) + &
                      ProjectPi3(3,1,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T13(i,j,k))) + &
                      ProjectPi3(1,2,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T21(i,j,k))) + &
                      ProjectPi3(2,2,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T21(i,j,k))) + &
                      ProjectPi3(3,2,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T21(i,j,k))) + &
                      ProjectPi3(1,2,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T22(i,j,k))) + &
                      ProjectPi3(2,2,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T22(i,j,k))) + &
                      ProjectPi3(3,2,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T22(i,j,k))) + &
                      ProjectPi3(1,2,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T23(i,j,k))) + &
                      ProjectPi3(2,2,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T23(i,j,k))) + &
                      ProjectPi3(3,2,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T23(i,j,k))) + &
                      ProjectPi3(1,3,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T31(i,j,k))) + &
                      ProjectPi3(2,3,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T31(i,j,k))) + &
                      ProjectPi3(3,3,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T31(i,j,k))) + &
                      ProjectPi3(1,3,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T32(i,j,k))) + &
                      ProjectPi3(2,3,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T32(i,j,k))) + &
                      ProjectPi3(3,3,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T32(i,j,k))) + &
                      ProjectPi3(1,3,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T33(i,j,k))) + &
                      ProjectPi3(2,3,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T33(i,j,k))) + &
                      ProjectPi3(3,3,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T33(i,j,k)))
      !
      Tccall=Tccall + ProjectPi3(1,1,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T11d(i,j,k))) + &
                      ProjectPi3(2,1,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T11d(i,j,k))) + &
                      ProjectPi3(3,1,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T11d(i,j,k))) + &
                      ProjectPi3(1,1,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T12d(i,j,k))) + &
                      ProjectPi3(2,1,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T12d(i,j,k))) + &
                      ProjectPi3(3,1,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T12d(i,j,k))) + &
                      ProjectPi3(1,1,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T13d(i,j,k))) + &
                      ProjectPi3(2,1,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T13d(i,j,k))) + &
                      ProjectPi3(3,1,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T13d(i,j,k))) + &
                      ProjectPi3(1,2,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T21d(i,j,k))) + &
                      ProjectPi3(2,2,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T21d(i,j,k))) + &
                      ProjectPi3(3,2,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T21d(i,j,k))) + &
                      ProjectPi3(1,2,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T22d(i,j,k))) + &
                      ProjectPi3(2,2,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T22d(i,j,k))) + &
                      ProjectPi3(3,2,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T22d(i,j,k))) + &
                      ProjectPi3(1,2,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T23d(i,j,k))) + &
                      ProjectPi3(2,2,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T23d(i,j,k))) + &
                      ProjectPi3(3,2,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T23d(i,j,k))) + &
                      ProjectPi3(1,3,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T31d(i,j,k))) + &
                      ProjectPi3(2,3,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T31d(i,j,k))) + &
                      ProjectPi3(3,3,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T31d(i,j,k))) + &
                      ProjectPi3(1,3,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T32d(i,j,k))) + &
                      ProjectPi3(2,3,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T32d(i,j,k))) + &
                      ProjectPi3(3,3,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T32d(i,j,k))) + &
                      ProjectPi3(1,3,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T33d(i,j,k))) + &
                      ProjectPi3(2,3,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T33d(i,j,k))) + &
                      ProjectPi3(3,3,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T33d(i,j,k)))
      !
      k2Ts = k2Ts + kk**2 * ProjectP3(1,1,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T11(i,j,k))) + &
                    kk**2 * ProjectP3(2,1,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T11(i,j,k))) + &
                    kk**2 * ProjectP3(3,1,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T11(i,j,k))) + &
                    kk**2 * ProjectP3(1,1,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T12(i,j,k))) + &
                    kk**2 * ProjectP3(2,1,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T12(i,j,k))) + &
                    kk**2 * ProjectP3(3,1,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T12(i,j,k))) + &
                    kk**2 * ProjectP3(1,1,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T13(i,j,k))) + &
                    kk**2 * ProjectP3(2,1,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T13(i,j,k))) + &
                    kk**2 * ProjectP3(3,1,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T13(i,j,k))) + &
                    kk**2 * ProjectP3(1,2,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T21(i,j,k))) + &
                    kk**2 * ProjectP3(2,2,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T21(i,j,k))) + &
                    kk**2 * ProjectP3(3,2,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T21(i,j,k))) + &
                    kk**2 * ProjectP3(1,2,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T22(i,j,k))) + &
                    kk**2 * ProjectP3(2,2,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T22(i,j,k))) + &
                    kk**2 * ProjectP3(3,2,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T22(i,j,k))) + &
                    kk**2 * ProjectP3(1,2,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T23(i,j,k))) + &
                    kk**2 * ProjectP3(2,2,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T23(i,j,k))) + &
                    kk**2 * ProjectP3(3,2,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T23(i,j,k))) + &
                    kk**2 * ProjectP3(1,3,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T31(i,j,k))) + &
                    kk**2 * ProjectP3(2,3,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T31(i,j,k))) + &
                    kk**2 * ProjectP3(3,3,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T31(i,j,k))) + &
                    kk**2 * ProjectP3(1,3,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T32(i,j,k))) + &
                    kk**2 * ProjectP3(2,3,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T32(i,j,k))) + &
                    kk**2 * ProjectP3(3,3,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T32(i,j,k))) + &
                    kk**2 * ProjectP3(1,3,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T33(i,j,k))) + &
                    kk**2 * ProjectP3(2,3,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T33(i,j,k))) + &
                    kk**2 * ProjectP3(3,3,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T33(i,j,k)))
      !
      k2Tss=k2Tss + kk**2 * ProjectP3(1,1,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T11s(i,j,k))) + &
                    kk**2 * ProjectP3(2,1,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T11s(i,j,k))) + &
                    kk**2 * ProjectP3(3,1,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T11s(i,j,k))) + &
                    kk**2 * ProjectP3(1,1,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T12s(i,j,k))) + &
                    kk**2 * ProjectP3(2,1,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T12s(i,j,k))) + &
                    kk**2 * ProjectP3(3,1,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T12s(i,j,k))) + &
                    kk**2 * ProjectP3(1,1,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T13s(i,j,k))) + &
                    kk**2 * ProjectP3(2,1,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T13s(i,j,k))) + &
                    kk**2 * ProjectP3(3,1,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T13s(i,j,k))) + &
                    kk**2 * ProjectP3(1,2,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T21s(i,j,k))) + &
                    kk**2 * ProjectP3(2,2,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T21s(i,j,k))) + &
                    kk**2 * ProjectP3(3,2,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T21s(i,j,k))) + &
                    kk**2 * ProjectP3(1,2,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T22s(i,j,k))) + &
                    kk**2 * ProjectP3(2,2,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T22s(i,j,k))) + &
                    kk**2 * ProjectP3(3,2,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T22s(i,j,k))) + &
                    kk**2 * ProjectP3(1,2,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T23s(i,j,k))) + &
                    kk**2 * ProjectP3(2,2,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T23s(i,j,k))) + &
                    kk**2 * ProjectP3(3,2,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T23s(i,j,k))) + &
                    kk**2 * ProjectP3(1,3,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T31s(i,j,k))) + &
                    kk**2 * ProjectP3(2,3,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T31s(i,j,k))) + &
                    kk**2 * ProjectP3(3,3,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T31s(i,j,k))) + &
                    kk**2 * ProjectP3(1,3,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T32s(i,j,k))) + &
                    kk**2 * ProjectP3(2,3,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T32s(i,j,k))) + &
                    kk**2 * ProjectP3(3,3,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T32s(i,j,k))) + &
                    kk**2 * ProjectP3(1,3,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T33s(i,j,k))) + &
                    kk**2 * ProjectP3(2,3,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T33s(i,j,k))) + &
                    kk**2 * ProjectP3(3,3,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T33s(i,j,k)))

      k2Tc =  k2Tc+ kk**2 * ProjectPi3(1,1,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T11(i,j,k))) + &
                    kk**2 * ProjectPi3(2,1,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T11(i,j,k))) + &
                    kk**2 * ProjectPi3(3,1,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T11(i,j,k))) + &
                    kk**2 * ProjectPi3(1,1,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T12(i,j,k))) + &
                    kk**2 * ProjectPi3(2,1,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T12(i,j,k))) + &
                    kk**2 * ProjectPi3(3,1,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T12(i,j,k))) + &
                    kk**2 * ProjectPi3(1,1,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T13(i,j,k))) + &
                    kk**2 * ProjectPi3(2,1,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T13(i,j,k))) + &
                    kk**2 * ProjectPi3(3,1,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T13(i,j,k))) + &
                    kk**2 * ProjectPi3(1,2,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T21(i,j,k))) + &
                    kk**2 * ProjectPi3(2,2,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T21(i,j,k))) + &
                    kk**2 * ProjectPi3(3,2,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T21(i,j,k))) + &
                    kk**2 * ProjectPi3(1,2,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T22(i,j,k))) + &
                    kk**2 * ProjectPi3(2,2,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T22(i,j,k))) + &
                    kk**2 * ProjectPi3(3,2,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T22(i,j,k))) + &
                    kk**2 * ProjectPi3(1,2,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T23(i,j,k))) + &
                    kk**2 * ProjectPi3(2,2,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T23(i,j,k))) + &
                    kk**2 * ProjectPi3(3,2,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T23(i,j,k))) + &
                    kk**2 * ProjectPi3(1,3,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T31(i,j,k))) + &
                    kk**2 * ProjectPi3(2,3,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T31(i,j,k))) + &
                    kk**2 * ProjectPi3(3,3,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T31(i,j,k))) + &
                    kk**2 * ProjectPi3(1,3,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T32(i,j,k))) + &
                    kk**2 * ProjectPi3(2,3,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T32(i,j,k))) + &
                    kk**2 * ProjectPi3(3,3,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T32(i,j,k))) + &
                    kk**2 * ProjectPi3(1,3,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T33(i,j,k))) + &
                    kk**2 * ProjectPi3(2,3,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T33(i,j,k))) + &
                    kk**2 * ProjectPi3(3,3,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T33(i,j,k)))
      !
      k2Tcc=k2Tcc + kk**2 * ProjectPi3(1,1,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T11d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,1,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T11d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,1,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T11d(i,j,k))) + &
                    kk**2 * ProjectPi3(1,1,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T12d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,1,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T12d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,1,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T12d(i,j,k))) + &
                    kk**2 * ProjectPi3(1,1,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T13d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,1,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T13d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,1,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T13d(i,j,k))) + &
                    kk**2 * ProjectPi3(1,2,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T21d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,2,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T21d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,2,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T21d(i,j,k))) + &
                    kk**2 * ProjectPi3(1,2,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T22d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,2,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T22d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,2,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T22d(i,j,k))) + &
                    kk**2 * ProjectPi3(1,2,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T23d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,2,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T23d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,2,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T23d(i,j,k))) + &
                    kk**2 * ProjectPi3(1,3,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T31d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,3,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T31d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,3,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T31d(i,j,k))) + &
                    kk**2 * ProjectPi3(1,3,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T32d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,3,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T32d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,3,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T32d(i,j,k))) + &
                    kk**2 * ProjectPi3(1,3,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T33d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,3,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T33d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,3,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T33d(i,j,k)))
    enddo
    enddo
    enddo
    !
    do i=0,allkmax
      Ts(i) = - psum(Ts(i))
      Tstheta(i) = psum(Tstheta(i))
      Td(i) = -psum(Td(i))
      Tdstheta(i) = psum(Tdstheta(i))
      Tddtheta(i) = psum(Tddtheta(i))
      Tss(i) = -psum(Tss(i))
      Tdd(i) = -psum(Tdd(i))
      Tcount(i) = psum(Tcount(i))
    end do
    !
    !
    Tsall = psum(Tsall)
    Tcall = psum(Tcall)
    k2Ts = psum(k2Ts)
    k2Tc = psum(k2Tc)
    Tssall = psum(Tssall)
    Tccall = psum(Tccall)
    k2Tss = psum(k2Tss)
    k2Tcc = psum(k2Tcc)
    !
    if(mpirank == 0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Tspec'//stepname//'.dat'
      else
        outfilename = 'pp/Tspec.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_a, &
                        firstline='nstep time k Ts Tss Tstheta Td Tdd Tdstheta Tddtheta')
      do i=0,allkmax
        call listwrite(hand_a,dble(i),Ts(i), Tss(i),Tstheta(i),Td(i),Tdd(i),Tdstheta(i),Tddtheta(i))
      end do
      !
      print*,' <<< pp/Tspec'//stepname//'.dat ... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Tspec_aux'//stepname//'.dat'
      else
        outfilename = 'pp/Tspec_aux.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
            firstline='nstep time k2Ts k2Tc Tsall Tcall Tssall Tccall k2Tss k2Tcc')
      call listwrite(hand_b,k2Ts,k2Tc,Tsall,Tcall,Tssall,Tccall,k2Tss,k2Tcc)
      !
      print*,' <<< '//outfilename//'... done.'
    endif
    !
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_u3spe)
    call fftw_free(c_pspe)
    call fftw_free(c_theta)
    call fftw_free(c_u1s)
    call fftw_free(c_u2s)
    call fftw_free(c_u3s)
    call fftw_free(c_u1c)
    call fftw_free(c_u2c)
    call fftw_free(c_u3c)
    call fftw_free(c_T11s)
    call fftw_free(c_T12s)
    call fftw_free(c_T13s)
    call fftw_free(c_T21s)
    call fftw_free(c_T22s)
    call fftw_free(c_T23s)
    call fftw_free(c_T31s)
    call fftw_free(c_T32s)
    call fftw_free(c_T33s)
    call fftw_free(c_T11d)
    call fftw_free(c_T12d)
    call fftw_free(c_T13d)
    call fftw_free(c_T21d)
    call fftw_free(c_T22d)
    call fftw_free(c_T23d)
    call fftw_free(c_T31d)
    call fftw_free(c_T32d)
    call fftw_free(c_T33d)
    call fftw_free(c_T11)
    call fftw_free(c_T12)
    call fftw_free(c_T13)
    call fftw_free(c_T21)
    call fftw_free(c_T22)
    call fftw_free(c_T23)
    call fftw_free(c_T31)
    call fftw_free(c_T32)
    call fftw_free(c_T33)
    call fftw_free(c_Tstheta1)
    call fftw_free(c_Tstheta2)
    call fftw_free(c_Tstheta3)
    call fftw_free(c_Tdtheta1)
    call fftw_free(c_Tdtheta2)
    call fftw_free(c_Tdtheta3)
    call mpistop
    !
  end subroutine instanttriad3D
  !
  subroutine instanttriad2Davg(thefilenumb)
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km
    use commarray, only: vel, rho, prs
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    integer,intent(in) :: thefilenumb
    character(len=128) :: infilename
    character(len=4) :: stepname
    real(8) :: u1mean,u2mean,rhomean,prsmean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1spe,u2spe,pspe
    complex(8), allocatable, dimension(:,:) :: usspe,ucspe
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: theta,u1c,u2c,u1s,u2s
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: T11,T12,T21,T22,Tstheta1,Tstheta2,Tdtheta1,Tdtheta2
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: T11s,T12s,T21s,T22s,T11d,T12d,T21d,T22d
    real(8), allocatable, dimension(:) :: Ts,Tstheta,Td,Tdd,Tcount,Tss,Tdstheta,Tddtheta
    real(8), allocatable, dimension(:,:) :: k1,k2,u1,u2
    integer :: allkmax
    real(8) :: k,dk,p,q,kx,ky !wave number
    character(len=128) :: outfilename
    integer :: hand_a,hand_b
    character(len=1) :: modeio
    integer :: i,j,n,s,t
    real(8) :: k2Ts, k2Tc,Tsall,Tcall, Tssall, Tccall, k2Tss, k2Tcc
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_pspe, c_theta,c_u1c,c_u2c,c_u1s,c_u2s
    type(C_PTR) :: c_T11,c_T12,c_T21,c_T22,c_Tstheta1,c_Tstheta2,c_Tdtheta1,c_Tdtheta2
    type(C_PTR) :: c_T11s,c_T12s,c_T21s,c_T22s,c_T11d,c_T12d,c_T21d,c_T22d
    call readinput
    !
    modeio='h'
    dk = 0.5d0
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    allkmax=ceiling(1.d0*min(ia,ja))
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,"allkmax:",allkmax
    if(ka .ne. 0) stop 'Please use instantspectra3D'
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km), prs(0:im,0:jm,0:km))
    !
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='p',  var=prs(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    allocate(u1(1:im,1:jm),u2(1:im,1:jm))
    ! Calculate favre average
    u1mean = 0.0d0
    u2mean = 0.0d0
    rhomean = 0.0d0
    prsmean = 0.0d0
    !
    do i=1,im
      do j=1,jm
        u1mean = u1mean + vel(i,j,0,1)
        u2mean = u2mean + vel(i,j,0,2)
        rhomean = rhomean + rho(i,j,0)
        prsmean = prsmean + prs(i,j,0)
      enddo
    enddo
    !
    rhomean = psum(rhomean) / (1.0d0*ia*ja)
    u1mean = psum(u1mean) / (1.d0*ia*ja)
    u2mean = psum(u2mean) / (1.d0*ia*ja)
    prsmean = psum(prsmean) / (1.d0*ia*ja)
    if(mpirank==0) print *, 'u1mean=',u1mean, 'u2mean=',u2mean, 'prsmean=',prsmean
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw])
    c_pspe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_pspe, pspe, [imfftw,jmfftw])
    !
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    !!!! Convert to complex
    do j=1,jm
    do i=1,im
      !
      u1(i,j) = vel(i,j,0,1)-u1mean
      u2(i,j) = vel(i,j,0,2)-u2mean
      u1spe(i,j)=CMPLX(u1(i,j),0.d0,C_INTPTR_T)
      u2spe(i,j)=CMPLX(u2(i,j),0.d0,C_INTPTR_T)
      pspe(i,j)=CMPLX(prs(i,j,0)-prsmean,0.d0,C_INTPTR_T)
      !
    end do
    end do
    !
    !!!! Do 2d FFT
    !
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    call fftw_mpi_execute_dft(forward_plan,pspe,pspe)
    !
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=u1spe(i,j)/(1.d0*ia*ja)
      u2spe(i,j)=u2spe(i,j)/(1.d0*ia*ja)
      pspe(i,j)=pspe(i,j)/(1.d0*ia*ja)
      !
    end do
    end do
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    do j=1,jm
    do i=1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if((j+j0) <= (ja/2+1)) then
        k2(i,j) = real(j+j0-1,8)
      else if((j+j0)<=(ja)) then
        k2(i,j) = real(j+j0-ja-1,8)
      else
        print *,"Error, no wave number possible, (j+j0) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    !
    allocate(usspe(1:im,1:jm),ucspe(1:im,1:jm))
    !
    c_theta = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_theta, theta, [imfftw,jmfftw])
    c_u1s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1s, u1s, [imfftw,jmfftw])
    c_u2s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2s, u2s, [imfftw,jmfftw])
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw])
    !
    !!!! Do S-C decomposition
    do j=1,jm
      do i=1,im
        !
        k=dsqrt(k1(i,j)**2+k2(i,j)**2+1.d-15)
        usspe(i,j) = u1spe(i,j)*k2(i,j)/k - u2spe(i,j)*k1(i,j)/k
        ucspe(i,j) = u1spe(i,j)*k1(i,j)/k + u2spe(i,j)*k2(i,j)/k
        !
        u1c(i,j)=  ucspe(i,j)*k1(i,j)/k 
        u2c(i,j)=  ucspe(i,j)*k2(i,j)/k
        u1s(i,j)=  usspe(i,j)*k2(i,j)/k 
        u2s(i,j)= -usspe(i,j)*k1(i,j)/k
        !
        theta(i,j)= CMPLX(0.d0,1.d0,C_INTPTR_T) * (ucspe(i,j) * k)
        !
      end do
    end do
    !
    !!!! Do inverse FFT
    !
    call fftw_mpi_execute_dft(backward_plan,theta,theta)
    call fftw_mpi_execute_dft(backward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(backward_plan,u2s,u2s)
    call fftw_mpi_execute_dft(backward_plan,u1c,u1c)
    call fftw_mpi_execute_dft(backward_plan,u2c,u2c)
    !
    c_T11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11,T11, [imfftw,jmfftw])
    c_T12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12,T12, [imfftw,jmfftw])
    c_T21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21,T21, [imfftw,jmfftw])
    c_T22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22,T22, [imfftw,jmfftw])
    !
    c_T11s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11s,T11s, [imfftw,jmfftw])
    c_T12s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12s,T12s, [imfftw,jmfftw])
    c_T21s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21s,T21s, [imfftw,jmfftw])
    c_T22s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22s,T22s, [imfftw,jmfftw])
    !
    c_T11d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11d,T11d, [imfftw,jmfftw])
    c_T12d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12d,T12d, [imfftw,jmfftw])
    c_T21d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21d,T21d, [imfftw,jmfftw])
    c_T22d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22d,T22d, [imfftw,jmfftw])
    !
    c_Tstheta1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tstheta1,Tstheta1, [imfftw,jmfftw])
    c_Tstheta2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tstheta2,Tstheta2, [imfftw,jmfftw])
    c_Tdtheta1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tdtheta1,Tdtheta1, [imfftw,jmfftw])
    c_Tdtheta2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tdtheta2,Tdtheta2, [imfftw,jmfftw])
    !
    ! This part is in physical space
    do j=1,jm
    do i=1,im
      T11(i,j) = u1(i,j)*u1(i,j)
      T12(i,j) = u1(i,j)*u2(i,j)
      T21(i,j) = u2(i,j)*u1(i,j)
      T22(i,j) = u2(i,j)*u2(i,j)
      !
      T11s(i,j) = u1s(i,j)*u1s(i,j)
      T12s(i,j) = u1s(i,j)*u2s(i,j)
      T21s(i,j) = u2s(i,j)*u1s(i,j)
      T22s(i,j) = u2s(i,j)*u2s(i,j)
      !
      T11d(i,j) = u1c(i,j)*u1c(i,j)
      T12d(i,j) = u1c(i,j)*u2c(i,j)
      T21d(i,j) = u2c(i,j)*u1c(i,j)
      T22d(i,j) = u2c(i,j)*u2c(i,j)
      !
      Tstheta1(i,j) = theta(i,j)*u1s(i,j)
      Tstheta2(i,j) = theta(i,j)*u2s(i,j)
      Tdtheta1(i,j) = theta(i,j)*u1c(i,j)
      Tdtheta2(i,j) = theta(i,j)*u2c(i,j)
    enddo
    enddo
    !
    call fftw_mpi_execute_dft(forward_plan,T11,T11)
    call fftw_mpi_execute_dft(forward_plan,T12,T12)
    call fftw_mpi_execute_dft(forward_plan,T21,T21)
    call fftw_mpi_execute_dft(forward_plan,T22,T22)
    !
    call fftw_mpi_execute_dft(forward_plan,T11s,T11s)
    call fftw_mpi_execute_dft(forward_plan,T12s,T12s)
    call fftw_mpi_execute_dft(forward_plan,T21s,T21s)
    call fftw_mpi_execute_dft(forward_plan,T22s,T22s)
    !
    call fftw_mpi_execute_dft(forward_plan,T11d,T11d)
    call fftw_mpi_execute_dft(forward_plan,T12d,T12d)
    call fftw_mpi_execute_dft(forward_plan,T21d,T21d)
    call fftw_mpi_execute_dft(forward_plan,T22d,T22d)
    !
    call fftw_mpi_execute_dft(forward_plan,Tstheta1,Tstheta1)
    call fftw_mpi_execute_dft(forward_plan,Tstheta2,Tstheta2)
    call fftw_mpi_execute_dft(forward_plan,Tdtheta1,Tdtheta1)
    call fftw_mpi_execute_dft(forward_plan,Tdtheta2,Tdtheta2)
    !
    call fftw_mpi_execute_dft(forward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(forward_plan,u2s,u2s)
    call fftw_mpi_execute_dft(forward_plan,u1c,u1c)
    call fftw_mpi_execute_dft(forward_plan,u2c,u2c)
    !
    do j=1,jm
    do i=1,im
      T11(i,j)=T11(i,j)/(1.d0*ia*ja)
      T12(i,j)=T12(i,j)/(1.d0*ia*ja)
      T21(i,j)=T21(i,j)/(1.d0*ia*ja)
      T22(i,j)=T22(i,j)/(1.d0*ia*ja)
      !
      T11s(i,j)=T11s(i,j)/(1.d0*ia*ja)
      T12s(i,j)=T12s(i,j)/(1.d0*ia*ja)
      T21s(i,j)=T21s(i,j)/(1.d0*ia*ja)
      T22s(i,j)=T22s(i,j)/(1.d0*ia*ja)
      !
      T11d(i,j)=T11d(i,j)/(1.d0*ia*ja)
      T12d(i,j)=T12d(i,j)/(1.d0*ia*ja)
      T21d(i,j)=T21d(i,j)/(1.d0*ia*ja)
      T22d(i,j)=T22d(i,j)/(1.d0*ia*ja)
      !
      Tstheta1(i,j)=Tstheta1(i,j)/(1.d0*ia*ja)
      Tstheta2(i,j)=Tstheta2(i,j)/(1.d0*ia*ja)
      Tdtheta1(i,j)=Tdtheta1(i,j)/(1.d0*ia*ja)
      Tdtheta2(i,j)=Tdtheta2(i,j)/(1.d0*ia*ja)
      !
      u1s(i,j)=u1s(i,j)/(1.d0*ia*ja)
      u2s(i,j)=u2s(i,j)/(1.d0*ia*ja)
      u1c(i,j)=u1c(i,j)/(1.d0*ia*ja)
      u2c(i,j)=u2c(i,j)/(1.d0*ia*ja)
    enddo
    enddo
    !
    !
    allocate(Ts(0:allkmax),Tstheta(0:allkmax),Td(0:allkmax),Tdstheta(0:allkmax))
    allocate(Tddtheta(0:allkmax),Tcount(0:allkmax),Tss(0:allkmax),Tdd(0:allkmax))
    !
    do i=0,allkmax
      Ts(i)= 0.0d0
      Tstheta(i)= 0.0d0
      Td(i)= 0.0d0
      Tddtheta(i)= 0.0d0
      Tdstheta(i)= 0.0d0
      Tcount(i)=0
      Tss(i)=0.d0
      Tdd(i)=0.d0
    end do
    !
    k2Ts = 0.d0
    k2Tc = 0.d0
    Tsall = 0.d0
    Tcall = 0.d0
    Tssall = 0.d0
    Tccall = 0.d0
    k2Tss = 0.d0
    k2Tcc = 0.d0
    !
    do j = 1,jm
    do i = 1,im
      kx = k1(i,j)
      ky = k2(i,j)
      k=dsqrt(kx**2+ky**2+1.d-15)
      !
      if(kint(k,dk)==0)then
        print *, u1s(i,j), u2s(i,j), T11s(i,j), T12s(i,j), T21s(i,j), T22s(i,j)
      endif
      !
      if(kint(k,dk)<= allkmax)then
        Ts(kint(k,dk)) = Ts(kint(k,dk)) + &
        k * ProjectP3(1,1,1,kx,ky) * dimag(u1s(i,j)*conjg(T11(i,j))) + &
        k * ProjectP3(2,1,1,kx,ky) * dimag(u2s(i,j)*conjg(T11(i,j))) + &
        k * ProjectP3(1,1,2,kx,ky) * dimag(u1s(i,j)*conjg(T12(i,j))) + &
        k * ProjectP3(2,1,2,kx,ky) * dimag(u2s(i,j)*conjg(T12(i,j))) + &
        k * ProjectP3(1,2,1,kx,ky) * dimag(u1s(i,j)*conjg(T21(i,j))) + &
        k * ProjectP3(2,2,1,kx,ky) * dimag(u2s(i,j)*conjg(T21(i,j))) + &
        k * ProjectP3(1,2,2,kx,ky) * dimag(u1s(i,j)*conjg(T22(i,j))) + &
        k * ProjectP3(2,2,2,kx,ky) * dimag(u2s(i,j)*conjg(T22(i,j)))

        Tss(kint(k,dk)) = Tss(kint(k,dk)) + &
        k * ProjectP3(1,1,1,kx,ky) * dimag(u1s(i,j)*conjg(T11s(i,j))) + &
        k * ProjectP3(2,1,1,kx,ky) * dimag(u2s(i,j)*conjg(T11s(i,j))) + &
        k * ProjectP3(1,1,2,kx,ky) * dimag(u1s(i,j)*conjg(T12s(i,j))) + &
        k * ProjectP3(2,1,2,kx,ky) * dimag(u2s(i,j)*conjg(T12s(i,j))) + &
        k * ProjectP3(1,2,1,kx,ky) * dimag(u1s(i,j)*conjg(T21s(i,j))) + &
        k * ProjectP3(2,2,1,kx,ky) * dimag(u2s(i,j)*conjg(T21s(i,j))) + &
        k * ProjectP3(1,2,2,kx,ky) * dimag(u1s(i,j)*conjg(T22s(i,j))) + &
        k * ProjectP3(2,2,2,kx,ky) * dimag(u2s(i,j)*conjg(T22s(i,j)))
            !
        Tstheta(kint(k,dk)) = Tstheta(kint(k,dk)) + &
        k * ProjectP2(1,1,kx,ky) * dreal(u1s(i,j)*conjg(Tstheta1(i,j)+Tdtheta1(i,j))) + &
        k * ProjectP2(1,2,kx,ky) * dreal(u1s(i,j)*conjg(Tstheta2(i,j)+Tdtheta2(i,j))) + &
        k * ProjectP2(2,1,kx,ky) * dreal(u2s(i,j)*conjg(Tstheta1(i,j)+Tstheta1(i,j))) + &
        k * ProjectP2(2,2,kx,ky) * dreal(u2s(i,j)*conjg(Tstheta2(i,j)+Tdtheta2(i,j)))
            !
        Td(kint(k,dk)) = Td(kint(k,dk)) + &
        k * ProjectPi3(1,1,1,kx,ky) * dimag(u1c(i,j)*conjg(T11(i,j))) + &
        k * ProjectPi3(2,1,1,kx,ky) * dimag(u2c(i,j)*conjg(T11(i,j))) + &
        k * ProjectPi3(1,1,2,kx,ky) * dimag(u1c(i,j)*conjg(T12(i,j))) + &
        k * ProjectPi3(2,1,2,kx,ky) * dimag(u2c(i,j)*conjg(T12(i,j))) + &
        k * ProjectPi3(1,2,1,kx,ky) * dimag(u1c(i,j)*conjg(T21(i,j))) + &
        k * ProjectPi3(2,2,1,kx,ky) * dimag(u2c(i,j)*conjg(T21(i,j))) + &
        k * ProjectPi3(1,2,2,kx,ky) * dimag(u1c(i,j)*conjg(T22(i,j))) + &
        k * ProjectPi3(2,2,2,kx,ky) * dimag(u2c(i,j)*conjg(T22(i,j)))
            !
        Tdd(kint(k,dk)) = Tdd(kint(k,dk)) + &
        k * ProjectPi3(1,1,1,kx,ky) * dimag(u1c(i,j)*conjg(T11d(i,j))) + &
        k * ProjectPi3(2,1,1,kx,ky) * dimag(u2c(i,j)*conjg(T11d(i,j))) + &
        k * ProjectPi3(1,1,2,kx,ky) * dimag(u1c(i,j)*conjg(T12d(i,j))) + &
        k * ProjectPi3(2,1,2,kx,ky) * dimag(u2c(i,j)*conjg(T12d(i,j))) + &
        k * ProjectPi3(1,2,1,kx,ky) * dimag(u1c(i,j)*conjg(T21d(i,j))) + &
        k * ProjectPi3(2,2,1,kx,ky) * dimag(u2c(i,j)*conjg(T21d(i,j))) + &
        k * ProjectPi3(1,2,2,kx,ky) * dimag(u1c(i,j)*conjg(T22d(i,j))) + &
        k * ProjectPi3(2,2,2,kx,ky) * dimag(u2c(i,j)*conjg(T22d(i,j)))
            !
        Tdstheta(kint(k,dk)) = Tdstheta(kint(k,dk)) + &
        k * ProjectPi2(1,1,kx,ky) * dreal(u1c(i,j)*conjg(Tstheta1(i,j))) + &
        k * ProjectPi2(1,2,kx,ky) * dreal(u1c(i,j)*conjg(Tstheta2(i,j))) + &
        k * ProjectPi2(2,1,kx,ky) * dreal(u2c(i,j)*conjg(Tstheta1(i,j))) + &
        k * ProjectPi2(2,2,kx,ky) * dreal(u2c(i,j)*conjg(Tstheta2(i,j)))
              !
        Tddtheta(kint(k,dk)) = Tddtheta(kint(k,dk)) + &
        k * ProjectPi2(1,1,kx,ky) * dreal(u1c(i,j)*conjg(Tdtheta1(i,j))) + &
        k * ProjectPi2(1,2,kx,ky) * dreal(u1c(i,j)*conjg(Tdtheta2(i,j))) + &
        k * ProjectPi2(2,1,kx,ky) * dreal(u2c(i,j)*conjg(Tdtheta1(i,j))) + &
        k * ProjectPi2(2,2,kx,ky) * dreal(u2c(i,j)*conjg(Tdtheta2(i,j)))
        Tcount(kint(k,dk)) = Tcount(kint(k,dk))+1
      endif
      Tsall = Tsall + ProjectP3(1,1,1,kx,ky) * dimag(u1s(i,j)*conjg(T11(i,j))) + &
                      ProjectP3(2,1,1,kx,ky) * dimag(u2s(i,j)*conjg(T11(i,j))) + &
                      ProjectP3(1,1,2,kx,ky) * dimag(u1s(i,j)*conjg(T12(i,j))) + &
                      ProjectP3(2,1,2,kx,ky) * dimag(u2s(i,j)*conjg(T12(i,j))) + &
                      ProjectP3(1,2,1,kx,ky) * dimag(u1s(i,j)*conjg(T21(i,j))) + &
                      ProjectP3(2,2,1,kx,ky) * dimag(u2s(i,j)*conjg(T21(i,j))) + &
                      ProjectP3(1,2,2,kx,ky) * dimag(u1s(i,j)*conjg(T22(i,j))) + &
                      ProjectP3(2,2,2,kx,ky) * dimag(u2s(i,j)*conjg(T22(i,j)))
      Tssall=Tssall + ProjectP3(1,1,1,kx,ky) * dimag(u1s(i,j)*conjg(T11s(i,j))) + &
                      ProjectP3(2,1,1,kx,ky) * dimag(u2s(i,j)*conjg(T11s(i,j))) + &
                      ProjectP3(1,1,2,kx,ky) * dimag(u1s(i,j)*conjg(T12s(i,j))) + &
                      ProjectP3(2,1,2,kx,ky) * dimag(u2s(i,j)*conjg(T12s(i,j))) + &
                      ProjectP3(1,2,1,kx,ky) * dimag(u1s(i,j)*conjg(T21s(i,j))) + &
                      ProjectP3(2,2,1,kx,ky) * dimag(u2s(i,j)*conjg(T21s(i,j))) + &
                      ProjectP3(1,2,2,kx,ky) * dimag(u1s(i,j)*conjg(T22s(i,j))) + &
                      ProjectP3(2,2,2,kx,ky) * dimag(u2s(i,j)*conjg(T22s(i,j)))
      Tcall = Tcall + ProjectPi3(1,1,1,kx,ky) * dimag(u1c(i,j)*conjg(T11(i,j))) + &
                      ProjectPi3(2,1,1,kx,ky) * dimag(u2c(i,j)*conjg(T11(i,j))) + &
                      ProjectPi3(1,1,2,kx,ky) * dimag(u1c(i,j)*conjg(T12(i,j))) + &
                      ProjectPi3(2,1,2,kx,ky) * dimag(u2c(i,j)*conjg(T12(i,j))) + &
                      ProjectPi3(1,2,1,kx,ky) * dimag(u1c(i,j)*conjg(T21(i,j))) + &
                      ProjectPi3(2,2,1,kx,ky) * dimag(u2c(i,j)*conjg(T21(i,j))) + &
                      ProjectPi3(1,2,2,kx,ky) * dimag(u1c(i,j)*conjg(T22(i,j))) + &
                      ProjectPi3(2,2,2,kx,ky) * dimag(u2c(i,j)*conjg(T22(i,j)))
      Tccall=Tccall + ProjectPi3(1,1,1,kx,ky) * dimag(u1c(i,j)*conjg(T11d(i,j))) + &
                      ProjectPi3(2,1,1,kx,ky) * dimag(u2c(i,j)*conjg(T11d(i,j))) + &
                      ProjectPi3(1,1,2,kx,ky) * dimag(u1c(i,j)*conjg(T12d(i,j))) + &
                      ProjectPi3(2,1,2,kx,ky) * dimag(u2c(i,j)*conjg(T12d(i,j))) + &
                      ProjectPi3(1,2,1,kx,ky) * dimag(u1c(i,j)*conjg(T21d(i,j))) + &
                      ProjectPi3(2,2,1,kx,ky) * dimag(u2c(i,j)*conjg(T21d(i,j))) + &
                      ProjectPi3(1,2,2,kx,ky) * dimag(u1c(i,j)*conjg(T22d(i,j))) + &
                      ProjectPi3(2,2,2,kx,ky) * dimag(u2c(i,j)*conjg(T22d(i,j)))
      k2Ts = k2Ts + k**2 * ProjectP3(1,1,1,kx,ky) * dimag(u1s(i,j)*conjg(T11(i,j))) + &
              k**2 * ProjectP3(2,1,1,kx,ky) * dimag(u2s(i,j)*conjg(T11(i,j))) + &
              k**2 * ProjectP3(1,1,2,kx,ky) * dimag(u1s(i,j)*conjg(T12(i,j))) + &
              k**2 * ProjectP3(2,1,2,kx,ky) * dimag(u2s(i,j)*conjg(T12(i,j))) + &
              k**2 * ProjectP3(1,2,1,kx,ky) * dimag(u1s(i,j)*conjg(T21(i,j))) + &
              k**2 * ProjectP3(2,2,1,kx,ky) * dimag(u2s(i,j)*conjg(T21(i,j))) + &
              k**2 * ProjectP3(1,2,2,kx,ky) * dimag(u1s(i,j)*conjg(T22(i,j))) + &
              k**2 * ProjectP3(2,2,2,kx,ky) * dimag(u2s(i,j)*conjg(T22(i,j)))
      k2Tss= k2Tss+ k**2 * ProjectP3(1,1,1,kx,ky) * dimag(u1s(i,j)*conjg(T11s(i,j))) + &
              k**2 * ProjectP3(2,1,1,kx,ky) * dimag(u2s(i,j)*conjg(T11s(i,j))) + &
              k**2 * ProjectP3(1,1,2,kx,ky) * dimag(u1s(i,j)*conjg(T12s(i,j))) + &
              k**2 * ProjectP3(2,1,2,kx,ky) * dimag(u2s(i,j)*conjg(T12s(i,j))) + &
              k**2 * ProjectP3(1,2,1,kx,ky) * dimag(u1s(i,j)*conjg(T21s(i,j))) + &
              k**2 * ProjectP3(2,2,1,kx,ky) * dimag(u2s(i,j)*conjg(T21s(i,j))) + &
              k**2 * ProjectP3(1,2,2,kx,ky) * dimag(u1s(i,j)*conjg(T22s(i,j))) + &
              k**2 * ProjectP3(2,2,2,kx,ky) * dimag(u2s(i,j)*conjg(T22s(i,j)))
      k2Tc =  k2Tc + k**2 * ProjectPi3(1,1,1,kx,ky) * dimag(u1c(i,j)*conjg(T11(i,j))) + &
              k**2 * ProjectPi3(2,1,1,kx,ky) * dimag(u2c(i,j)*conjg(T11(i,j))) + &
              k**2 * ProjectPi3(1,1,2,kx,ky) * dimag(u1c(i,j)*conjg(T12(i,j))) + &
              k**2 * ProjectPi3(2,1,2,kx,ky) * dimag(u2c(i,j)*conjg(T12(i,j))) + &
              k**2 * ProjectPi3(1,2,1,kx,ky) * dimag(u1c(i,j)*conjg(T21(i,j))) + &
              k**2 * ProjectPi3(2,2,1,kx,ky) * dimag(u2c(i,j)*conjg(T21(i,j))) + &
              k**2 * ProjectPi3(1,2,2,kx,ky) * dimag(u1c(i,j)*conjg(T22(i,j))) + &
              k**2 * ProjectPi3(2,2,2,kx,ky) * dimag(u2c(i,j)*conjg(T22(i,j)))
      k2Tcc=  k2Tcc+ k**2 * ProjectPi3(1,1,1,kx,ky) * dimag(u1c(i,j)*conjg(T11d(i,j))) + &
              k**2 * ProjectPi3(2,1,1,kx,ky) * dimag(u2c(i,j)*conjg(T11d(i,j))) + &
              k**2 * ProjectPi3(1,1,2,kx,ky) * dimag(u1c(i,j)*conjg(T12d(i,j))) + &
              k**2 * ProjectPi3(2,1,2,kx,ky) * dimag(u2c(i,j)*conjg(T12d(i,j))) + &
              k**2 * ProjectPi3(1,2,1,kx,ky) * dimag(u1c(i,j)*conjg(T21d(i,j))) + &
              k**2 * ProjectPi3(2,2,1,kx,ky) * dimag(u2c(i,j)*conjg(T21d(i,j))) + &
              k**2 * ProjectPi3(1,2,2,kx,ky) * dimag(u1c(i,j)*conjg(T22d(i,j))) + &
              k**2 * ProjectPi3(2,2,2,kx,ky) * dimag(u2c(i,j)*conjg(T22d(i,j)))
    enddo
    enddo
    !
    do i=0,allkmax
      Ts(i) = psum(Ts(i))
      Tstheta(i) = psum(Tstheta(i))
      Td(i) = psum(Td(i))
      Tdstheta(i) = psum(Tdstheta(i))
      Tddtheta(i) = psum(Tddtheta(i))
      Tss(i) = psum(Tss(i))
      Tdd(i) = psum(Tdd(i))
      Tcount(i) = psum(Tcount(i))
    end do
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    do i=0,allkmax
      if(Tcount(i) .ne. 0) then
        Ts(i) = -Ts(i)/Tcount(i)*2*pi
        Tstheta(i) = Tstheta(i)/Tcount(i)*2*pi
        Td(i) = -Td(i)/Tcount(i)*2*pi
        Tdstheta(i) = Tdstheta(i)/Tcount(i)*2*pi
        Tddtheta(i) = Tddtheta(i)/Tcount(i)*2*pi
        Tss(i) = -Tss(i)/Tcount(i)*2*pi
        Tdd(i) = -Tdd(i)/Tcount(i)*2*pi
      endif
    end do
    !
    Tsall = psum(Tsall)
    Tcall = psum(Tcall)
    k2Ts = psum(k2Ts)
    k2Tc = psum(k2Tc)
    Tssall = psum(Tssall)
    Tccall = psum(Tccall)
    k2Tss = psum(k2Tss)
    k2Tcc = psum(k2Tcc)
    !
    if(mpirank == 0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Tspec'//stepname//'.dat'
      else
        outfilename = 'pp/Tspec.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_a, &
                        firstline='nstep time k Ts Tss Tstheta Td Tdd Tdstheta Tddtheta')
      do i=0,allkmax
        call listwrite(hand_a,dble(i),Ts(i),Tss(i),Tstheta(i),Td(i),Tdd(i),Tdstheta(i),Tddtheta(i))
      end do
      !
      print*,' <<< pp/Tspec'//stepname//'.dat ... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Tspec_aux'//stepname//'.dat'
      else
        outfilename = 'pp/Tspec_aux.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
            firstline='nstep time k2Ts k2Tc Tsall Tcall Tssall Tccall k2Tss k2Tcc')
      call listwrite(hand_b,k2Ts,k2Tc,Tsall,Tcall,Tssall,Tccall,k2Tss,k2Tcc)
      !
      print*,' <<< '//outfilename//'... done.'
    endif
    !
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_pspe)
    call fftw_free(c_theta)
    call fftw_free(c_u1s)
    call fftw_free(c_u2s)
    call fftw_free(c_u1c)
    call fftw_free(c_u2c)
    call fftw_free(c_T11)
    call fftw_free(c_T12)
    call fftw_free(c_T21)
    call fftw_free(c_T22)
    call fftw_free(c_T11s)
    call fftw_free(c_T12s)
    call fftw_free(c_T21s)
    call fftw_free(c_T22s)
    call fftw_free(c_T11d)
    call fftw_free(c_T12d)
    call fftw_free(c_T21d)
    call fftw_free(c_T22d)
    call fftw_free(c_Tstheta1)
    call fftw_free(c_Tstheta2)
    call fftw_free(c_Tdtheta1)
    call fftw_free(c_Tdtheta2)
    call mpistop
    !
  end subroutine instanttriad2Davg
  !
  subroutine instanttriad3Davg(thefilenumb)
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,   only : time,nstep,im,jm,km
    use commarray, only : vel, rho, prs
    use hdf5io
    use utility,   only : listinit,listwrite
    use parallel,  only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    integer,intent(in) :: thefilenumb
    character(len=128) :: infilename
    character(len=4) :: stepname
    real(8) :: u1mean,u2mean,u3mean,rhomean,prsmean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: u1spe,u2spe,u3spe,pspe
    complex(8), allocatable, dimension(:,:,:) :: ucspe
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: theta,u1c,u2c,u3c,u1s,u2s,u3s
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: T11,T12,T13,T21,T22,T23,T31,T32,T33
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: Tstheta1,Tstheta2,Tstheta3,Tdtheta1,Tdtheta2,Tdtheta3
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: T11s,T12s,T13s,T21s,T22s,T23s,T31s,T32s,T33s
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: T11d,T12d,T13d,T21d,T22d,T23d,T31d,T32d,T33d
    real(8), allocatable, dimension(:) :: Ts,Tstheta,Td,Tddtheta,Tdstheta,Tcount,Tss,Tdd
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3,u1,u2,u3
    integer :: allkmax
    real(8) :: kk,dk,p,q,kx,ky,kz !wave number
    character(len=128) :: outfilename
    integer :: hand_a,hand_b
    character(len=1) :: modeio
    integer :: i,j,n,s,t,k
    real(8) :: k2Ts, k2Tc, Tsall, Tcall, Tssall, Tccall, k2Tss, k2Tcc
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_u3spe, c_pspe
    type(C_PTR) :: c_theta,c_u1c,c_u2c,c_u3c,c_u1s,c_u2s,c_u3s
    type(C_PTR) :: c_T11,c_T12,c_T13,c_T21,c_T22,c_T23,c_T31,c_T32,c_T33
    type(C_PTR) :: c_Tstheta1,c_Tstheta2,c_Tstheta3,c_Tdtheta1,c_Tdtheta2,c_Tdtheta3
    type(C_PTR) :: c_T11s,c_T12s,c_T13s,c_T21s,c_T22s,c_T23s,c_T31s,c_T32s,c_T33s
    type(C_PTR) :: c_T11d,c_T12d,c_T13d,c_T21d,c_T22d,c_T23d,c_T31d,c_T32d,c_T33d
    call readinput
    !
    modeio='h'
    dk = 0.5d0
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    allkmax=ceiling(1.d0*min(ia,ja,ka))
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka,"allkmax:",allkmax
    if(ka==0) stop 'Please use instantspectra2D'
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km), prs(0:im,0:jm,0:km))
    !
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='p',  var=prs(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    allocate(u1(1:im,1:jm,1:km),u2(1:im,1:jm,1:km),u3(1:im,1:jm,1:km))
    ! Calculate favre average
    u1mean = 0.0d0
    u2mean = 0.0d0
    u3mean = 0.0d0
    rhomean = 0.0d0
    prsmean = 0.0d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      u1mean = u1mean + vel(i,j,k,1)
      u2mean = u2mean + vel(i,j,k,2)
      u3mean = u3mean + vel(i,j,k,3)
      rhomean = rhomean + rho(i,j,k)
      prsmean = prsmean + prs(i,j,k)
    enddo
    enddo
    enddo
    !
    rhomean = psum(rhomean) / (1.0d0*ia*ja*ka)
    u1mean = psum(u1mean) / (1.d0*ia*ja*ka)
    u2mean = psum(u2mean) / (1.d0*ia*ja*ka)
    u3mean = psum(u3mean) / (1.d0*ia*ja*ka)
    prsmean = psum(prsmean) / (1.d0*ia*ja*ka)
    if(mpirank==0) print *, 'u1mean=',u1mean, 'u2mean=',u2mean, 'u3mean=',u3mean, 'prsmean=',prsmean
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw,kmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw,kmfftw])
    c_u3spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3spe, u3spe, [imfftw,jmfftw,kmfftw])
    c_pspe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_pspe, pspe, [imfftw,jmfftw,kmfftw])
    !
    forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    !!!! Convert to complex
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1(i,j,k) = vel(i,j,k,1)-u1mean
      u2(i,j,k) = vel(i,j,k,2)-u2mean
      u3(i,j,k) = vel(i,j,k,3)-u3mean
      u1spe(i,j,k)=CMPLX(u1(i,j,k),0.d0,C_INTPTR_T)
      u2spe(i,j,k)=CMPLX(u2(i,j,k),0.d0,C_INTPTR_T)
      u3spe(i,j,k)=CMPLX(u3(i,j,k),0.d0,C_INTPTR_T)
      pspe(i,j,k)=CMPLX(prs(i,j,k)-prsmean,0.d0,C_INTPTR_T)
      !
    end do
    end do
    end do
    !
    !!!! Do FFT
    !
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    call fftw_mpi_execute_dft(forward_plan,u3spe,u3spe)
    call fftw_mpi_execute_dft(forward_plan,pspe,pspe)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j,k)=u1spe(i,j,k)/(1.d0*ia*ja*ka)
      u2spe(i,j,k)=u2spe(i,j,k)/(1.d0*ia*ja*ka)
      u3spe(i,j,k)=u3spe(i,j,k)/(1.d0*ia*ja*ka)
      pspe(i,j,k)=pspe(i,j,k)/(1.d0*ia*ja*ka)
      !
    end do
    end do
    end do
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j,k) = real(i-1,8)
      else if(i<=ia) then
        k1(i,j,k) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if(j <= (ja/2+1)) then
        k2(i,j,k) = real(j-1,8)
      else if(j<=ja) then
        k2(i,j,k) = real(j-ja-1,8)
      else
        print *,"Error, no wave number possible, j must smaller than ja-1 !"
      end if
      !
      if((k+k0) <= (ka/2+1)) then
        k3(i,j,k) = real(k+k0-1,8)
      else if((k+k0)<=ka) then
        k3(i,j,k) = real(k+k0-ka-1,8)
      else
        print *,"Error, no wave number possible, (k+k0) must smaller than ka-1 !"
      end if
      !
    enddo
    enddo
    enddo
    !
    allocate(ucspe(1:im,1:jm,1:km))
    !
    c_theta = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_theta, theta, [imfftw,jmfftw,kmfftw])
    c_u1s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1s, u1s, [imfftw,jmfftw,kmfftw])
    c_u2s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2s, u2s, [imfftw,jmfftw,kmfftw])
    c_u3s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3s, u3s, [imfftw,jmfftw,kmfftw])
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw,kmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw,kmfftw])
    c_u3c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3c, u3c, [imfftw,jmfftw,kmfftw])
    !
    !!!! Do S-C decomposition
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      kk=dsqrt(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2+1.d-15)
      !
      ucspe(i,j,k) = k1(i,j,k)/kk * u1spe(i,j,k) + k2(i,j,k)/kk * u2spe(i,j,k) + k3(i,j,k)/kk * u3spe(i,j,k)
      u1c(i,j,k)=  k1(i,j,k)*k1(i,j,k)/(kk**2) * u1spe(i,j,k) + k1(i,j,k)*k2(i,j,k)/(kk**2) * u2spe(i,j,k) &
                + k1(i,j,k)*k3(i,j,k)/(kk**2) * u3spe(i,j,k)
      u2c(i,j,k)=  k2(i,j,k)*k1(i,j,k)/(kk**2) * u1spe(i,j,k) + k2(i,j,k)*k2(i,j,k)/(kk**2) * u2spe(i,j,k) &
                + k2(i,j,k)*k3(i,j,k)/(kk**2) * u3spe(i,j,k)
      u3c(i,j,k)=  k3(i,j,k)*k1(i,j,k)/(kk**2) * u1spe(i,j,k) + k3(i,j,k)*k2(i,j,k)/(kk**2) * u2spe(i,j,k) &
                + k3(i,j,k)*k3(i,j,k)/(kk**2) * u3spe(i,j,k)
      u1s(i,j,k)=  u1spe(i,j,k) - u1c(i,j,k)
      u2s(i,j,k)=  u2spe(i,j,k) - u2c(i,j,k)
      u3s(i,j,k)=  u3spe(i,j,k) - u3c(i,j,k)
      !
      theta(i,j,k)= CMPLX(0.d0,1.d0,C_INTPTR_T) * (ucspe(i,j,k) * kk)
      !
    end do
    end do
    enddo
    !
    !!!! Do inverse FFT
    !
    call fftw_mpi_execute_dft(backward_plan,theta,theta)
    call fftw_mpi_execute_dft(backward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(backward_plan,u2s,u2s)
    call fftw_mpi_execute_dft(backward_plan,u3s,u3s)
    call fftw_mpi_execute_dft(backward_plan,u1c,u1c)
    call fftw_mpi_execute_dft(backward_plan,u2c,u2c)
    call fftw_mpi_execute_dft(backward_plan,u3c,u3c)
    !
    c_T11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11,T11, [imfftw,jmfftw,kmfftw])
    c_T12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12,T12, [imfftw,jmfftw,kmfftw])
    c_T13 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T13,T13, [imfftw,jmfftw,kmfftw])
    c_T21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21,T21, [imfftw,jmfftw,kmfftw])
    c_T22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22,T22, [imfftw,jmfftw,kmfftw])
    c_T23 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T23,T23, [imfftw,jmfftw,kmfftw])
    c_T31 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T31,T31, [imfftw,jmfftw,kmfftw])
    c_T32 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T32,T32, [imfftw,jmfftw,kmfftw])
    c_T33 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T33,T33, [imfftw,jmfftw,kmfftw])
    !
    c_T11s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11s,T11s, [imfftw,jmfftw,kmfftw])
    c_T12s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12s,T12s, [imfftw,jmfftw,kmfftw])
    c_T13s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T13s,T13s, [imfftw,jmfftw,kmfftw])
    c_T21s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21s,T21s, [imfftw,jmfftw,kmfftw])
    c_T22s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22s,T22s, [imfftw,jmfftw,kmfftw])
    c_T23s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T23s,T23s, [imfftw,jmfftw,kmfftw])
    c_T31s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T31s,T31s, [imfftw,jmfftw,kmfftw])
    c_T32s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T32s,T32s, [imfftw,jmfftw,kmfftw])
    c_T33s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T33s,T33s, [imfftw,jmfftw,kmfftw])
    !
    c_T11d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11d,T11d, [imfftw,jmfftw,kmfftw])
    c_T12d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12d,T12d, [imfftw,jmfftw,kmfftw])
    c_T13d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T13d,T13d, [imfftw,jmfftw,kmfftw])
    c_T21d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21d,T21d, [imfftw,jmfftw,kmfftw])
    c_T22d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22d,T22d, [imfftw,jmfftw,kmfftw])
    c_T23d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T23d,T23d, [imfftw,jmfftw,kmfftw])
    c_T31d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T31d,T31d, [imfftw,jmfftw,kmfftw])
    c_T32d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T32d,T32d, [imfftw,jmfftw,kmfftw])
    c_T33d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T33d,T33d, [imfftw,jmfftw,kmfftw])
    !
    c_Tstheta1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tstheta1,Tstheta1, [imfftw,jmfftw,kmfftw])
    c_Tstheta2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tstheta2,Tstheta2, [imfftw,jmfftw,kmfftw])
    c_Tstheta3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tstheta3,Tstheta3, [imfftw,jmfftw,kmfftw])
    c_Tdtheta1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tdtheta1,Tdtheta1, [imfftw,jmfftw,kmfftw])
    c_Tdtheta2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tdtheta2,Tdtheta2, [imfftw,jmfftw,kmfftw])
    c_Tdtheta3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tdtheta3,Tdtheta3, [imfftw,jmfftw,kmfftw])
    !
    ! This part is in physical space
    do k=1,km
    do j=1,jm
    do i=1,im
      T11(i,j,k) = u1(i,j,k)*u1(i,j,k)
      T12(i,j,k) = u1(i,j,k)*u2(i,j,k)
      T13(i,j,k) = u1(i,j,k)*u3(i,j,k)
      T21(i,j,k) = u2(i,j,k)*u1(i,j,k)
      T22(i,j,k) = u2(i,j,k)*u2(i,j,k)
      T23(i,j,k) = u2(i,j,k)*u3(i,j,k)
      T31(i,j,k) = u3(i,j,k)*u1(i,j,k)
      T32(i,j,k) = u3(i,j,k)*u2(i,j,k)
      T33(i,j,k) = u3(i,j,k)*u3(i,j,k)
      !
      T11s(i,j,k) = u1s(i,j,k)*u1s(i,j,k)
      T12s(i,j,k) = u1s(i,j,k)*u2s(i,j,k)
      T13s(i,j,k) = u1s(i,j,k)*u3s(i,j,k)
      T21s(i,j,k) = u2s(i,j,k)*u1s(i,j,k)
      T22s(i,j,k) = u2s(i,j,k)*u2s(i,j,k)
      T23s(i,j,k) = u2s(i,j,k)*u3s(i,j,k)
      T31s(i,j,k) = u3s(i,j,k)*u1s(i,j,k)
      T32s(i,j,k) = u3s(i,j,k)*u2s(i,j,k)
      T33s(i,j,k) = u3s(i,j,k)*u3s(i,j,k)
      !
      T11d(i,j,k) = u1c(i,j,k)*u1c(i,j,k)
      T12d(i,j,k) = u1c(i,j,k)*u2c(i,j,k)
      T13d(i,j,k) = u1c(i,j,k)*u3c(i,j,k)
      T21d(i,j,k) = u2c(i,j,k)*u1c(i,j,k)
      T22d(i,j,k) = u2c(i,j,k)*u2c(i,j,k)
      T23d(i,j,k) = u2c(i,j,k)*u3c(i,j,k)
      T31d(i,j,k) = u3c(i,j,k)*u1c(i,j,k)
      T32d(i,j,k) = u3c(i,j,k)*u2c(i,j,k)
      T33d(i,j,k) = u3c(i,j,k)*u3c(i,j,k)
      !
      Tstheta1(i,j,k) = theta(i,j,k)*u1s(i,j,k)
      Tstheta2(i,j,k) = theta(i,j,k)*u2s(i,j,k)
      Tstheta3(i,j,k) = theta(i,j,k)*u3s(i,j,k)
      Tdtheta1(i,j,k) = theta(i,j,k)*u1c(i,j,k)
      Tdtheta2(i,j,k) = theta(i,j,k)*u2c(i,j,k)
      Tdtheta3(i,j,k) = theta(i,j,k)*u3c(i,j,k)
    enddo
    enddo
    enddo
    !
    call fftw_mpi_execute_dft(forward_plan,T11,T11)
    call fftw_mpi_execute_dft(forward_plan,T12,T12)
    call fftw_mpi_execute_dft(forward_plan,T13,T13)
    call fftw_mpi_execute_dft(forward_plan,T21,T21)
    call fftw_mpi_execute_dft(forward_plan,T22,T22)
    call fftw_mpi_execute_dft(forward_plan,T23,T23)
    call fftw_mpi_execute_dft(forward_plan,T31,T31)
    call fftw_mpi_execute_dft(forward_plan,T32,T32)
    call fftw_mpi_execute_dft(forward_plan,T33,T33)
    !
    call fftw_mpi_execute_dft(forward_plan,T11s,T11s)
    call fftw_mpi_execute_dft(forward_plan,T12s,T12s)
    call fftw_mpi_execute_dft(forward_plan,T13s,T13s)
    call fftw_mpi_execute_dft(forward_plan,T21s,T21s)
    call fftw_mpi_execute_dft(forward_plan,T22s,T22s)
    call fftw_mpi_execute_dft(forward_plan,T23s,T23s)
    call fftw_mpi_execute_dft(forward_plan,T31s,T31s)
    call fftw_mpi_execute_dft(forward_plan,T32s,T32s)
    call fftw_mpi_execute_dft(forward_plan,T33s,T33s)
    !
    call fftw_mpi_execute_dft(forward_plan,T11d,T11d)
    call fftw_mpi_execute_dft(forward_plan,T12d,T12d)
    call fftw_mpi_execute_dft(forward_plan,T13d,T13d)
    call fftw_mpi_execute_dft(forward_plan,T21d,T21d)
    call fftw_mpi_execute_dft(forward_plan,T22d,T22d)
    call fftw_mpi_execute_dft(forward_plan,T23d,T23d)
    call fftw_mpi_execute_dft(forward_plan,T31d,T31d)
    call fftw_mpi_execute_dft(forward_plan,T32d,T32d)
    call fftw_mpi_execute_dft(forward_plan,T33d,T33d)
    !
    call fftw_mpi_execute_dft(forward_plan,Tstheta1,Tstheta1)
    call fftw_mpi_execute_dft(forward_plan,Tstheta2,Tstheta2)
    call fftw_mpi_execute_dft(forward_plan,Tstheta3,Tstheta3)
    call fftw_mpi_execute_dft(forward_plan,Tdtheta1,Tdtheta1)
    call fftw_mpi_execute_dft(forward_plan,Tdtheta2,Tdtheta2)
    call fftw_mpi_execute_dft(forward_plan,Tdtheta3,Tdtheta3)
    !
    call fftw_mpi_execute_dft(forward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(forward_plan,u2s,u2s)
    call fftw_mpi_execute_dft(forward_plan,u3s,u3s)
    call fftw_mpi_execute_dft(forward_plan,u1c,u1c)
    call fftw_mpi_execute_dft(forward_plan,u2c,u2c)
    call fftw_mpi_execute_dft(forward_plan,u3c,u3c)
    !
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      T11(i,j,k)=T11(i,j,k)/(1.d0*ia*ja*ka)
      T12(i,j,k)=T12(i,j,k)/(1.d0*ia*ja*ka)
      T13(i,j,k)=T13(i,j,k)/(1.d0*ia*ja*ka)
      T21(i,j,k)=T21(i,j,k)/(1.d0*ia*ja*ka)
      T22(i,j,k)=T22(i,j,k)/(1.d0*ia*ja*ka)
      T23(i,j,k)=T23(i,j,k)/(1.d0*ia*ja*ka)
      T31(i,j,k)=T31(i,j,k)/(1.d0*ia*ja*ka)
      T32(i,j,k)=T32(i,j,k)/(1.d0*ia*ja*ka)
      T33(i,j,k)=T33(i,j,k)/(1.d0*ia*ja*ka)
      !
      T11s(i,j,k)=T11s(i,j,k)/(1.d0*ia*ja*ka)
      T12s(i,j,k)=T12s(i,j,k)/(1.d0*ia*ja*ka)
      T13s(i,j,k)=T13s(i,j,k)/(1.d0*ia*ja*ka)
      T21s(i,j,k)=T21s(i,j,k)/(1.d0*ia*ja*ka)
      T22s(i,j,k)=T22s(i,j,k)/(1.d0*ia*ja*ka)
      T23s(i,j,k)=T23s(i,j,k)/(1.d0*ia*ja*ka)
      T31s(i,j,k)=T31s(i,j,k)/(1.d0*ia*ja*ka)
      T32s(i,j,k)=T32s(i,j,k)/(1.d0*ia*ja*ka)
      T33s(i,j,k)=T33s(i,j,k)/(1.d0*ia*ja*ka)
      !
      T11d(i,j,k)=T11d(i,j,k)/(1.d0*ia*ja*ka)
      T12d(i,j,k)=T12d(i,j,k)/(1.d0*ia*ja*ka)
      T13d(i,j,k)=T13d(i,j,k)/(1.d0*ia*ja*ka)
      T21d(i,j,k)=T21d(i,j,k)/(1.d0*ia*ja*ka)
      T22d(i,j,k)=T22d(i,j,k)/(1.d0*ia*ja*ka)
      T23d(i,j,k)=T23d(i,j,k)/(1.d0*ia*ja*ka)
      T31d(i,j,k)=T31d(i,j,k)/(1.d0*ia*ja*ka)
      T32d(i,j,k)=T32d(i,j,k)/(1.d0*ia*ja*ka)
      T33d(i,j,k)=T33d(i,j,k)/(1.d0*ia*ja*ka)
      !
      Tstheta1(i,j,k)=Tstheta1(i,j,k)/(1.d0*ia*ja*ka)
      Tstheta2(i,j,k)=Tstheta2(i,j,k)/(1.d0*ia*ja*ka)
      Tstheta3(i,j,k)=Tstheta3(i,j,k)/(1.d0*ia*ja*ka)
      Tdtheta1(i,j,k)=Tdtheta1(i,j,k)/(1.d0*ia*ja*ka)
      Tdtheta2(i,j,k)=Tdtheta2(i,j,k)/(1.d0*ia*ja*ka)
      Tdtheta3(i,j,k)=Tdtheta3(i,j,k)/(1.d0*ia*ja*ka)
      !
      u1s(i,j,k)=u1s(i,j,k)/(1.d0*ia*ja*ka)
      u2s(i,j,k)=u2s(i,j,k)/(1.d0*ia*ja*ka)
      u3s(i,j,k)=u3s(i,j,k)/(1.d0*ia*ja*ka)
      u1c(i,j,k)=u1c(i,j,k)/(1.d0*ia*ja*ka)
      u2c(i,j,k)=u2c(i,j,k)/(1.d0*ia*ja*ka)
      u3c(i,j,k)=u3c(i,j,k)/(1.d0*ia*ja*ka)
      !
    enddo
    enddo
    enddo
    !
    allocate(Ts(0:allkmax),Tstheta(0:allkmax),Td(0:allkmax),Tdstheta(0:allkmax))
    allocate(Tddtheta(0:allkmax),Tcount(0:allkmax),Tss(0:allkmax),Tdd(0:allkmax))
    !
    do i=0,allkmax
      Ts(i)= 0.0d0
      Tstheta(i)= 0.0d0
      Td(i)= 0.0d0
      Tddtheta(i)= 0.0d0
      Tdstheta(i)= 0.0d0
      Tcount(i)=0
      Tss(i)=0.d0
      Tdd(i)=0.d0
    end do
    !
    k2Ts = 0.d0
    k2Tc = 0.d0
    Tsall = 0.d0
    Tcall = 0.d0
    Tssall = 0.d0
    Tccall = 0.d0
    k2Tss = 0.d0
    k2Tcc = 0.d0
    !
    do k = 1,km
    do j = 1,jm
    do i = 1,im
      kx = k1(i,j,k)
      ky = k2(i,j,k)
      kz = k3(i,j,k)
      kk=dsqrt(kx**2+ky**2+kz**2+1.d-15)
      !
      !
      if(kint(kk,dk)<= allkmax)then
        Ts(kint(kk,dk)) = Ts(kint(kk,dk)) + &
                          kk**2 * ProjectP3(1,1,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T11(i,j,k))) + &
                          kk**2 * ProjectP3(2,1,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T11(i,j,k))) + &
                          kk**2 * ProjectP3(3,1,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T11(i,j,k))) + &
                          kk**2 * ProjectP3(1,1,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T12(i,j,k))) + &
                          kk**2 * ProjectP3(2,1,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T12(i,j,k))) + &
                          kk**2 * ProjectP3(3,1,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T12(i,j,k))) + &
                          kk**2 * ProjectP3(1,1,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T13(i,j,k))) + &
                          kk**2 * ProjectP3(2,1,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T13(i,j,k))) + &
                          kk**2 * ProjectP3(3,1,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T13(i,j,k))) + &
                          kk**2 * ProjectP3(1,2,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T21(i,j,k))) + &
                          kk**2 * ProjectP3(2,2,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T21(i,j,k))) + &
                          kk**2 * ProjectP3(3,2,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T21(i,j,k))) + &
                          kk**2 * ProjectP3(1,2,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T22(i,j,k))) + &
                          kk**2 * ProjectP3(2,2,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T22(i,j,k))) + &
                          kk**2 * ProjectP3(3,2,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T22(i,j,k))) + &
                          kk**2 * ProjectP3(1,2,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T23(i,j,k))) + &
                          kk**2 * ProjectP3(2,2,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T23(i,j,k))) + &
                          kk**2 * ProjectP3(3,2,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T23(i,j,k))) + &
                          kk**2 * ProjectP3(1,3,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T31(i,j,k))) + &
                          kk**2 * ProjectP3(2,3,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T31(i,j,k))) + &
                          kk**2 * ProjectP3(3,3,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T31(i,j,k))) + &
                          kk**2 * ProjectP3(1,3,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T32(i,j,k))) + &
                          kk**2 * ProjectP3(2,3,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T32(i,j,k))) + &
                          kk**2 * ProjectP3(3,3,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T32(i,j,k))) + &
                          kk**2 * ProjectP3(1,3,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T33(i,j,k))) + &
                          kk**2 * ProjectP3(2,3,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T33(i,j,k))) + &
                          kk**2 * ProjectP3(3,3,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T33(i,j,k)))
            !
        Tss(kint(kk,dk)) = Tss(kint(kk,dk)) + &
        kk**2 * ProjectP3(1,1,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T11s(i,j,k))) + &
        kk**2 * ProjectP3(2,1,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T11s(i,j,k))) + &
        kk**2 * ProjectP3(3,1,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T11s(i,j,k))) + &
        kk**2 * ProjectP3(1,1,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T12s(i,j,k))) + &
        kk**2 * ProjectP3(2,1,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T12s(i,j,k))) + &
        kk**2 * ProjectP3(3,1,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T12s(i,j,k))) + &
        kk**2 * ProjectP3(1,1,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T13s(i,j,k))) + &
        kk**2 * ProjectP3(2,1,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T13s(i,j,k))) + &
        kk**2 * ProjectP3(3,1,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T13s(i,j,k))) + &
        kk**2 * ProjectP3(1,2,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T21s(i,j,k))) + &
        kk**2 * ProjectP3(2,2,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T21s(i,j,k))) + &
        kk**2 * ProjectP3(3,2,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T21s(i,j,k))) + &
        kk**2 * ProjectP3(1,2,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T22s(i,j,k))) + &
        kk**2 * ProjectP3(2,2,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T22s(i,j,k))) + &
        kk**2 * ProjectP3(3,2,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T22s(i,j,k))) + &
        kk**2 * ProjectP3(1,2,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T23s(i,j,k))) + &
        kk**2 * ProjectP3(2,2,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T23s(i,j,k))) + &
        kk**2 * ProjectP3(3,2,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T23s(i,j,k))) + &
        kk**2 * ProjectP3(1,3,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T31s(i,j,k))) + &
        kk**2 * ProjectP3(2,3,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T31s(i,j,k))) + &
        kk**2 * ProjectP3(3,3,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T31s(i,j,k))) + &
        kk**2 * ProjectP3(1,3,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T32s(i,j,k))) + &
        kk**2 * ProjectP3(2,3,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T32s(i,j,k))) + &
        kk**2 * ProjectP3(3,3,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T32s(i,j,k))) + &
        kk**2 * ProjectP3(1,3,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T33s(i,j,k))) + &
        kk**2 * ProjectP3(2,3,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T33s(i,j,k))) + &
        kk**2 * ProjectP3(3,3,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T33s(i,j,k)))
              !
        Tstheta(kint(kk,dk)) = Tstheta(kint(kk,dk)) + &
        kk**2 * ProjectP2(1,1,kx,ky,kz) * dreal(u1s(i,j,k)*conjg(Tstheta1(i,j,k)+Tdtheta1(i,j,k))) + &
        kk**2 * ProjectP2(1,2,kx,ky,kz) * dreal(u1s(i,j,k)*conjg(Tstheta2(i,j,k)+Tdtheta2(i,j,k))) + &
        kk**2 * ProjectP2(1,3,kx,ky,kz) * dreal(u1s(i,j,k)*conjg(Tstheta3(i,j,k)+Tdtheta3(i,j,k))) + &
        kk**2 * ProjectP2(2,1,kx,ky,kz) * dreal(u2s(i,j,k)*conjg(Tstheta1(i,j,k)+Tdtheta1(i,j,k))) + &
        kk**2 * ProjectP2(2,2,kx,ky,kz) * dreal(u2s(i,j,k)*conjg(Tstheta2(i,j,k)+Tdtheta2(i,j,k))) + &
        kk**2 * ProjectP2(2,3,kx,ky,kz) * dreal(u2s(i,j,k)*conjg(Tstheta3(i,j,k)+Tdtheta3(i,j,k))) + &
        kk**2 * ProjectP2(3,1,kx,ky,kz) * dreal(u3s(i,j,k)*conjg(Tstheta1(i,j,k)+Tdtheta1(i,j,k))) + &
        kk**2 * ProjectP2(3,2,kx,ky,kz) * dreal(u3s(i,j,k)*conjg(Tstheta2(i,j,k)+Tdtheta2(i,j,k))) + &
        kk**2 * ProjectP2(3,3,kx,ky,kz) * dreal(u3s(i,j,k)*conjg(Tstheta3(i,j,k)+Tdtheta3(i,j,k)))
            !
        Td(kint(kk,dk)) = Td(kint(kk,dk)) + &
        kk**2 * ProjectPi3(1,1,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T11(i,j,k))) + &
        kk**2 * ProjectPi3(2,1,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T11(i,j,k))) + &
        kk**2 * ProjectPi3(3,1,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T11(i,j,k))) + &
        kk**2 * ProjectPi3(1,1,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T12(i,j,k))) + &
        kk**2 * ProjectPi3(2,1,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T12(i,j,k))) + &
        kk**2 * ProjectPi3(3,1,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T12(i,j,k))) + &
        kk**2 * ProjectPi3(1,1,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T13(i,j,k))) + &
        kk**2 * ProjectPi3(2,1,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T13(i,j,k))) + &
        kk**2 * ProjectPi3(3,1,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T13(i,j,k))) + &
        kk**2 * ProjectPi3(1,2,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T21(i,j,k))) + &
        kk**2 * ProjectPi3(2,2,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T21(i,j,k))) + &
        kk**2 * ProjectPi3(3,2,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T21(i,j,k))) + &
        kk**2 * ProjectPi3(1,2,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T22(i,j,k))) + &
        kk**2 * ProjectPi3(2,2,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T22(i,j,k))) + &
        kk**2 * ProjectPi3(3,2,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T22(i,j,k))) + &
        kk**2 * ProjectPi3(1,2,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T23(i,j,k))) + &
        kk**2 * ProjectPi3(2,2,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T23(i,j,k))) + &
        kk**2 * ProjectPi3(3,2,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T23(i,j,k))) + &
        kk**2 * ProjectPi3(1,3,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T31(i,j,k))) + &
        kk**2 * ProjectPi3(2,3,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T31(i,j,k))) + &
        kk**2 * ProjectPi3(3,3,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T31(i,j,k))) + &
        kk**2 * ProjectPi3(1,3,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T32(i,j,k))) + &
        kk**2 * ProjectPi3(2,3,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T32(i,j,k))) + &
        kk**2 * ProjectPi3(3,3,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T32(i,j,k))) + &
        kk**2 * ProjectPi3(1,3,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T33(i,j,k))) + &
        kk**2 * ProjectPi3(2,3,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T33(i,j,k))) + &
        kk**2 * ProjectPi3(3,3,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T33(i,j,k)))
          !
        Tdd(kint(kk,dk)) = Tdd(kint(kk,dk)) + &
        kk**2 * ProjectPi3(1,1,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T11d(i,j,k))) + &
        kk**2 * ProjectPi3(2,1,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T11d(i,j,k))) + &
        kk**2 * ProjectPi3(3,1,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T11d(i,j,k))) + &
        kk**2 * ProjectPi3(1,1,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T12d(i,j,k))) + &
        kk**2 * ProjectPi3(2,1,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T12d(i,j,k))) + &
        kk**2 * ProjectPi3(3,1,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T12d(i,j,k))) + &
        kk**2 * ProjectPi3(1,1,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T13d(i,j,k))) + &
        kk**2 * ProjectPi3(2,1,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T13d(i,j,k))) + &
        kk**2 * ProjectPi3(3,1,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T13d(i,j,k))) + &
        kk**2 * ProjectPi3(1,2,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T21d(i,j,k))) + &
        kk**2 * ProjectPi3(2,2,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T21d(i,j,k))) + &
        kk**2 * ProjectPi3(3,2,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T21d(i,j,k))) + &
        kk**2 * ProjectPi3(1,2,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T22d(i,j,k))) + &
        kk**2 * ProjectPi3(2,2,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T22d(i,j,k))) + &
        kk**2 * ProjectPi3(3,2,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T22d(i,j,k))) + &
        kk**2 * ProjectPi3(1,2,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T23d(i,j,k))) + &
        kk**2 * ProjectPi3(2,2,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T23d(i,j,k))) + &
        kk**2 * ProjectPi3(3,2,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T23d(i,j,k))) + &
        kk**2 * ProjectPi3(1,3,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T31d(i,j,k))) + &
        kk**2 * ProjectPi3(2,3,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T31d(i,j,k))) + &
        kk**2 * ProjectPi3(3,3,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T31d(i,j,k))) + &
        kk**2 * ProjectPi3(1,3,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T32d(i,j,k))) + &
        kk**2 * ProjectPi3(2,3,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T32d(i,j,k))) + &
        kk**2 * ProjectPi3(3,3,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T32d(i,j,k))) + &
        kk**2 * ProjectPi3(1,3,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T33d(i,j,k))) + &
        kk**2 * ProjectPi3(2,3,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T33d(i,j,k))) + &
        kk**2 * ProjectPi3(3,3,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T33d(i,j,k)))
            !
        Tdstheta(kint(kk,dk)) = Tdstheta(kint(kk,dk)) + &
        kk**2 * ProjectPi2(1,1,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tstheta1(i,j,k))) + &
        kk**2 * ProjectPi2(1,2,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tstheta2(i,j,k))) + &
        kk**2 * ProjectPi2(1,3,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tstheta3(i,j,k))) + &
        kk**2 * ProjectPi2(2,1,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tstheta1(i,j,k))) + &
        kk**2 * ProjectPi2(2,2,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tstheta2(i,j,k))) + &
        kk**2 * ProjectPi2(2,3,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tstheta3(i,j,k))) + &
        kk**2 * ProjectPi2(3,1,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tstheta1(i,j,k))) + &
        kk**2 * ProjectPi2(3,2,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tstheta2(i,j,k))) + &
        kk**2 * ProjectPi2(3,3,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tstheta3(i,j,k)))
                                !
        Tddtheta(kint(kk,dk)) = Tddtheta(kint(kk,dk)) + &
        kk**2 * ProjectPi2(1,1,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tdtheta1(i,j,k))) + &
        kk**2 * ProjectPi2(1,2,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tdtheta2(i,j,k))) + &
        kk**2 * ProjectPi2(1,3,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tdtheta3(i,j,k))) + &
        kk**2 * ProjectPi2(2,1,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tdtheta1(i,j,k))) + &
        kk**2 * ProjectPi2(2,2,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tdtheta2(i,j,k))) + &
        kk**2 * ProjectPi2(2,3,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tdtheta3(i,j,k))) + &
        kk**2 * ProjectPi2(3,1,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tdtheta1(i,j,k))) + &
        kk**2 * ProjectPi2(3,2,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tdtheta2(i,j,k))) + &
        kk**2 * ProjectPi2(3,3,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tdtheta3(i,j,k)))
        !
        Tcount(kint(kk,dk)) = Tcount(kint(kk,dk))+1
      endif
      Tsall = Tsall + ProjectP3(1,1,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T11(i,j,k))) + &
                      ProjectP3(2,1,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T11(i,j,k))) + &
                      ProjectP3(3,1,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T11(i,j,k))) + &
                      ProjectP3(1,1,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T12(i,j,k))) + &
                      ProjectP3(2,1,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T12(i,j,k))) + &
                      ProjectP3(3,1,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T12(i,j,k))) + &
                      ProjectP3(1,1,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T13(i,j,k))) + &
                      ProjectP3(2,1,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T13(i,j,k))) + &
                      ProjectP3(3,1,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T13(i,j,k))) + &
                      ProjectP3(1,2,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T21(i,j,k))) + &
                      ProjectP3(2,2,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T21(i,j,k))) + &
                      ProjectP3(3,2,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T21(i,j,k))) + &
                      ProjectP3(1,2,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T22(i,j,k))) + &
                      ProjectP3(2,2,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T22(i,j,k))) + &
                      ProjectP3(3,2,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T22(i,j,k))) + &
                      ProjectP3(1,2,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T23(i,j,k))) + &
                      ProjectP3(2,2,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T23(i,j,k))) + &
                      ProjectP3(3,2,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T23(i,j,k))) + &
                      ProjectP3(1,3,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T31(i,j,k))) + &
                      ProjectP3(2,3,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T31(i,j,k))) + &
                      ProjectP3(3,3,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T31(i,j,k))) + &
                      ProjectP3(1,3,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T32(i,j,k))) + &
                      ProjectP3(2,3,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T32(i,j,k))) + &
                      ProjectP3(3,3,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T32(i,j,k))) + &
                      ProjectP3(1,3,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T33(i,j,k))) + &
                      ProjectP3(2,3,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T33(i,j,k))) + &
                      ProjectP3(3,3,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T33(i,j,k)))
      !
      Tssall=Tssall + ProjectP3(1,1,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T11s(i,j,k))) + &
                      ProjectP3(2,1,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T11s(i,j,k))) + &
                      ProjectP3(3,1,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T11s(i,j,k))) + &
                      ProjectP3(1,1,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T12s(i,j,k))) + &
                      ProjectP3(2,1,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T12s(i,j,k))) + &
                      ProjectP3(3,1,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T12s(i,j,k))) + &
                      ProjectP3(1,1,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T13s(i,j,k))) + &
                      ProjectP3(2,1,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T13s(i,j,k))) + &
                      ProjectP3(3,1,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T13s(i,j,k))) + &
                      ProjectP3(1,2,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T21s(i,j,k))) + &
                      ProjectP3(2,2,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T21s(i,j,k))) + &
                      ProjectP3(3,2,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T21s(i,j,k))) + &
                      ProjectP3(1,2,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T22s(i,j,k))) + &
                      ProjectP3(2,2,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T22s(i,j,k))) + &
                      ProjectP3(3,2,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T22s(i,j,k))) + &
                      ProjectP3(1,2,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T23s(i,j,k))) + &
                      ProjectP3(2,2,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T23s(i,j,k))) + &
                      ProjectP3(3,2,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T23s(i,j,k))) + &
                      ProjectP3(1,3,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T31s(i,j,k))) + &
                      ProjectP3(2,3,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T31s(i,j,k))) + &
                      ProjectP3(3,3,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T31s(i,j,k))) + &
                      ProjectP3(1,3,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T32s(i,j,k))) + &
                      ProjectP3(2,3,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T32s(i,j,k))) + &
                      ProjectP3(3,3,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T32s(i,j,k))) + &
                      ProjectP3(1,3,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T33s(i,j,k))) + &
                      ProjectP3(2,3,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T33s(i,j,k))) + &
                      ProjectP3(3,3,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T33s(i,j,k)))
      !
      Tcall = Tcall + ProjectPi3(1,1,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T11(i,j,k))) + &
                      ProjectPi3(2,1,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T11(i,j,k))) + &
                      ProjectPi3(3,1,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T11(i,j,k))) + &
                      ProjectPi3(1,1,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T12(i,j,k))) + &
                      ProjectPi3(2,1,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T12(i,j,k))) + &
                      ProjectPi3(3,1,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T12(i,j,k))) + &
                      ProjectPi3(1,1,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T13(i,j,k))) + &
                      ProjectPi3(2,1,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T13(i,j,k))) + &
                      ProjectPi3(3,1,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T13(i,j,k))) + &
                      ProjectPi3(1,2,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T21(i,j,k))) + &
                      ProjectPi3(2,2,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T21(i,j,k))) + &
                      ProjectPi3(3,2,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T21(i,j,k))) + &
                      ProjectPi3(1,2,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T22(i,j,k))) + &
                      ProjectPi3(2,2,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T22(i,j,k))) + &
                      ProjectPi3(3,2,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T22(i,j,k))) + &
                      ProjectPi3(1,2,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T23(i,j,k))) + &
                      ProjectPi3(2,2,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T23(i,j,k))) + &
                      ProjectPi3(3,2,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T23(i,j,k))) + &
                      ProjectPi3(1,3,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T31(i,j,k))) + &
                      ProjectPi3(2,3,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T31(i,j,k))) + &
                      ProjectPi3(3,3,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T31(i,j,k))) + &
                      ProjectPi3(1,3,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T32(i,j,k))) + &
                      ProjectPi3(2,3,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T32(i,j,k))) + &
                      ProjectPi3(3,3,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T32(i,j,k))) + &
                      ProjectPi3(1,3,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T33(i,j,k))) + &
                      ProjectPi3(2,3,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T33(i,j,k))) + &
                      ProjectPi3(3,3,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T33(i,j,k)))
      !
      Tccall=Tccall + ProjectPi3(1,1,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T11d(i,j,k))) + &
                      ProjectPi3(2,1,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T11d(i,j,k))) + &
                      ProjectPi3(3,1,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T11d(i,j,k))) + &
                      ProjectPi3(1,1,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T12d(i,j,k))) + &
                      ProjectPi3(2,1,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T12d(i,j,k))) + &
                      ProjectPi3(3,1,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T12d(i,j,k))) + &
                      ProjectPi3(1,1,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T13d(i,j,k))) + &
                      ProjectPi3(2,1,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T13d(i,j,k))) + &
                      ProjectPi3(3,1,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T13d(i,j,k))) + &
                      ProjectPi3(1,2,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T21d(i,j,k))) + &
                      ProjectPi3(2,2,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T21d(i,j,k))) + &
                      ProjectPi3(3,2,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T21d(i,j,k))) + &
                      ProjectPi3(1,2,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T22d(i,j,k))) + &
                      ProjectPi3(2,2,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T22d(i,j,k))) + &
                      ProjectPi3(3,2,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T22d(i,j,k))) + &
                      ProjectPi3(1,2,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T23d(i,j,k))) + &
                      ProjectPi3(2,2,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T23d(i,j,k))) + &
                      ProjectPi3(3,2,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T23d(i,j,k))) + &
                      ProjectPi3(1,3,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T31d(i,j,k))) + &
                      ProjectPi3(2,3,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T31d(i,j,k))) + &
                      ProjectPi3(3,3,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T31d(i,j,k))) + &
                      ProjectPi3(1,3,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T32d(i,j,k))) + &
                      ProjectPi3(2,3,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T32d(i,j,k))) + &
                      ProjectPi3(3,3,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T32d(i,j,k))) + &
                      ProjectPi3(1,3,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T33d(i,j,k))) + &
                      ProjectPi3(2,3,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T33d(i,j,k))) + &
                      ProjectPi3(3,3,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T33d(i,j,k)))
      !
      k2Ts = k2Ts + kk**2 * ProjectP3(1,1,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T11(i,j,k))) + &
                    kk**2 * ProjectP3(2,1,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T11(i,j,k))) + &
                    kk**2 * ProjectP3(3,1,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T11(i,j,k))) + &
                    kk**2 * ProjectP3(1,1,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T12(i,j,k))) + &
                    kk**2 * ProjectP3(2,1,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T12(i,j,k))) + &
                    kk**2 * ProjectP3(3,1,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T12(i,j,k))) + &
                    kk**2 * ProjectP3(1,1,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T13(i,j,k))) + &
                    kk**2 * ProjectP3(2,1,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T13(i,j,k))) + &
                    kk**2 * ProjectP3(3,1,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T13(i,j,k))) + &
                    kk**2 * ProjectP3(1,2,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T21(i,j,k))) + &
                    kk**2 * ProjectP3(2,2,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T21(i,j,k))) + &
                    kk**2 * ProjectP3(3,2,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T21(i,j,k))) + &
                    kk**2 * ProjectP3(1,2,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T22(i,j,k))) + &
                    kk**2 * ProjectP3(2,2,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T22(i,j,k))) + &
                    kk**2 * ProjectP3(3,2,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T22(i,j,k))) + &
                    kk**2 * ProjectP3(1,2,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T23(i,j,k))) + &
                    kk**2 * ProjectP3(2,2,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T23(i,j,k))) + &
                    kk**2 * ProjectP3(3,2,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T23(i,j,k))) + &
                    kk**2 * ProjectP3(1,3,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T31(i,j,k))) + &
                    kk**2 * ProjectP3(2,3,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T31(i,j,k))) + &
                    kk**2 * ProjectP3(3,3,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T31(i,j,k))) + &
                    kk**2 * ProjectP3(1,3,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T32(i,j,k))) + &
                    kk**2 * ProjectP3(2,3,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T32(i,j,k))) + &
                    kk**2 * ProjectP3(3,3,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T32(i,j,k))) + &
                    kk**2 * ProjectP3(1,3,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T33(i,j,k))) + &
                    kk**2 * ProjectP3(2,3,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T33(i,j,k))) + &
                    kk**2 * ProjectP3(3,3,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T33(i,j,k)))
      !
      k2Tss=k2Tss + kk**2 * ProjectP3(1,1,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T11s(i,j,k))) + &
                    kk**2 * ProjectP3(2,1,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T11s(i,j,k))) + &
                    kk**2 * ProjectP3(3,1,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T11s(i,j,k))) + &
                    kk**2 * ProjectP3(1,1,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T12s(i,j,k))) + &
                    kk**2 * ProjectP3(2,1,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T12s(i,j,k))) + &
                    kk**2 * ProjectP3(3,1,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T12s(i,j,k))) + &
                    kk**2 * ProjectP3(1,1,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T13s(i,j,k))) + &
                    kk**2 * ProjectP3(2,1,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T13s(i,j,k))) + &
                    kk**2 * ProjectP3(3,1,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T13s(i,j,k))) + &
                    kk**2 * ProjectP3(1,2,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T21s(i,j,k))) + &
                    kk**2 * ProjectP3(2,2,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T21s(i,j,k))) + &
                    kk**2 * ProjectP3(3,2,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T21s(i,j,k))) + &
                    kk**2 * ProjectP3(1,2,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T22s(i,j,k))) + &
                    kk**2 * ProjectP3(2,2,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T22s(i,j,k))) + &
                    kk**2 * ProjectP3(3,2,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T22s(i,j,k))) + &
                    kk**2 * ProjectP3(1,2,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T23s(i,j,k))) + &
                    kk**2 * ProjectP3(2,2,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T23s(i,j,k))) + &
                    kk**2 * ProjectP3(3,2,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T23s(i,j,k))) + &
                    kk**2 * ProjectP3(1,3,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T31s(i,j,k))) + &
                    kk**2 * ProjectP3(2,3,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T31s(i,j,k))) + &
                    kk**2 * ProjectP3(3,3,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T31s(i,j,k))) + &
                    kk**2 * ProjectP3(1,3,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T32s(i,j,k))) + &
                    kk**2 * ProjectP3(2,3,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T32s(i,j,k))) + &
                    kk**2 * ProjectP3(3,3,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T32s(i,j,k))) + &
                    kk**2 * ProjectP3(1,3,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T33s(i,j,k))) + &
                    kk**2 * ProjectP3(2,3,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T33s(i,j,k))) + &
                    kk**2 * ProjectP3(3,3,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T33s(i,j,k)))

      k2Tc =  k2Tc+ kk**2 * ProjectPi3(1,1,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T11(i,j,k))) + &
                    kk**2 * ProjectPi3(2,1,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T11(i,j,k))) + &
                    kk**2 * ProjectPi3(3,1,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T11(i,j,k))) + &
                    kk**2 * ProjectPi3(1,1,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T12(i,j,k))) + &
                    kk**2 * ProjectPi3(2,1,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T12(i,j,k))) + &
                    kk**2 * ProjectPi3(3,1,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T12(i,j,k))) + &
                    kk**2 * ProjectPi3(1,1,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T13(i,j,k))) + &
                    kk**2 * ProjectPi3(2,1,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T13(i,j,k))) + &
                    kk**2 * ProjectPi3(3,1,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T13(i,j,k))) + &
                    kk**2 * ProjectPi3(1,2,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T21(i,j,k))) + &
                    kk**2 * ProjectPi3(2,2,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T21(i,j,k))) + &
                    kk**2 * ProjectPi3(3,2,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T21(i,j,k))) + &
                    kk**2 * ProjectPi3(1,2,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T22(i,j,k))) + &
                    kk**2 * ProjectPi3(2,2,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T22(i,j,k))) + &
                    kk**2 * ProjectPi3(3,2,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T22(i,j,k))) + &
                    kk**2 * ProjectPi3(1,2,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T23(i,j,k))) + &
                    kk**2 * ProjectPi3(2,2,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T23(i,j,k))) + &
                    kk**2 * ProjectPi3(3,2,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T23(i,j,k))) + &
                    kk**2 * ProjectPi3(1,3,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T31(i,j,k))) + &
                    kk**2 * ProjectPi3(2,3,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T31(i,j,k))) + &
                    kk**2 * ProjectPi3(3,3,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T31(i,j,k))) + &
                    kk**2 * ProjectPi3(1,3,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T32(i,j,k))) + &
                    kk**2 * ProjectPi3(2,3,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T32(i,j,k))) + &
                    kk**2 * ProjectPi3(3,3,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T32(i,j,k))) + &
                    kk**2 * ProjectPi3(1,3,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T33(i,j,k))) + &
                    kk**2 * ProjectPi3(2,3,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T33(i,j,k))) + &
                    kk**2 * ProjectPi3(3,3,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T33(i,j,k)))
      !
      k2Tcc=k2Tcc + kk**2 * ProjectPi3(1,1,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T11d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,1,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T11d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,1,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T11d(i,j,k))) + &
                    kk**2 * ProjectPi3(1,1,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T12d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,1,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T12d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,1,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T12d(i,j,k))) + &
                    kk**2 * ProjectPi3(1,1,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T13d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,1,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T13d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,1,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T13d(i,j,k))) + &
                    kk**2 * ProjectPi3(1,2,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T21d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,2,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T21d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,2,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T21d(i,j,k))) + &
                    kk**2 * ProjectPi3(1,2,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T22d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,2,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T22d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,2,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T22d(i,j,k))) + &
                    kk**2 * ProjectPi3(1,2,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T23d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,2,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T23d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,2,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T23d(i,j,k))) + &
                    kk**2 * ProjectPi3(1,3,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T31d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,3,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T31d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,3,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T31d(i,j,k))) + &
                    kk**2 * ProjectPi3(1,3,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T32d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,3,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T32d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,3,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T32d(i,j,k))) + &
                    kk**2 * ProjectPi3(1,3,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T33d(i,j,k))) + &
                    kk**2 * ProjectPi3(2,3,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T33d(i,j,k))) + &
                    kk**2 * ProjectPi3(3,3,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T33d(i,j,k)))
    enddo
    enddo
    enddo
    !
    do i=0,allkmax
      Ts(i) = psum(Ts(i))
      Tstheta(i) = psum(Tstheta(i))
      Td(i) = psum(Td(i))
      Tdstheta(i) = psum(Tdstheta(i))
      Tddtheta(i) = psum(Tddtheta(i))
      Tss(i) = psum(Tss(i))
      Tdd(i) = psum(Tdd(i))
      Tcount(i) = psum(Tcount(i))
    end do
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    do i=0,allkmax
      if(Tcount(i) .ne. 0) then
        Ts(i) = -Ts(i)/Tcount(i)*4*pi
        Tss(i) = -Tss(i)/Tcount(i)*4*pi
        Tstheta(i) = Tstheta(i)/Tcount(i)*4*pi
        Td(i) = - Td(i)/Tcount(i)*4*pi
        Tdd(i) = - Tdd(i)/Tcount(i)*4*pi
        Tddtheta(i) = Tddtheta(i)/Tcount(i)*4*pi
        Tdstheta(i) = Tdstheta(i)/Tcount(i)*4*pi
      endif
    enddo
    !
    Tsall = psum(Tsall)
    Tcall = psum(Tcall)
    k2Ts = psum(k2Ts)
    k2Tc = psum(k2Tc)
    Tssall = psum(Tssall)
    Tccall = psum(Tccall)
    k2Tss = psum(k2Tss)
    k2Tcc = psum(k2Tcc)
    !
    if(mpirank == 0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Tspec'//stepname//'.dat'
      else
        outfilename = 'pp/Tspec.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_a, &
                        firstline='nstep time k Ts Tss Tstheta Td Tdd Tdstheta Tddtheta')
      do i=0,allkmax
        call listwrite(hand_a,dble(i),Ts(i), Tss(i),Tstheta(i),Td(i),Tdd(i),Tdstheta(i),Tddtheta(i))
      end do
      !
      print*,' <<< pp/Tspec'//stepname//'.dat ... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Tspec_aux'//stepname//'.dat'
      else
        outfilename = 'pp/Tspec_aux.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
            firstline='nstep time k2Ts k2Tc Tsall Tcall Tssall Tccall k2Tss k2Tcc')
      call listwrite(hand_b,k2Ts,k2Tc,Tsall,Tcall,Tssall,Tccall,k2Tss,k2Tcc)
      !
      print*,' <<< '//outfilename//'... done.'
    endif
    !
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_u3spe)
    call fftw_free(c_pspe)
    call fftw_free(c_theta)
    call fftw_free(c_u1s)
    call fftw_free(c_u2s)
    call fftw_free(c_u3s)
    call fftw_free(c_u1c)
    call fftw_free(c_u2c)
    call fftw_free(c_u3c)
    call fftw_free(c_T11s)
    call fftw_free(c_T12s)
    call fftw_free(c_T13s)
    call fftw_free(c_T21s)
    call fftw_free(c_T22s)
    call fftw_free(c_T23s)
    call fftw_free(c_T31s)
    call fftw_free(c_T32s)
    call fftw_free(c_T33s)
    call fftw_free(c_T11d)
    call fftw_free(c_T12d)
    call fftw_free(c_T13d)
    call fftw_free(c_T21d)
    call fftw_free(c_T22d)
    call fftw_free(c_T23d)
    call fftw_free(c_T31d)
    call fftw_free(c_T32d)
    call fftw_free(c_T33d)
    call fftw_free(c_T11)
    call fftw_free(c_T12)
    call fftw_free(c_T13)
    call fftw_free(c_T21)
    call fftw_free(c_T22)
    call fftw_free(c_T23)
    call fftw_free(c_T31)
    call fftw_free(c_T32)
    call fftw_free(c_T33)
    call fftw_free(c_Tstheta1)
    call fftw_free(c_Tstheta2)
    call fftw_free(c_Tstheta3)
    call fftw_free(c_Tdtheta1)
    call fftw_free(c_Tdtheta2)
    call fftw_free(c_Tdtheta3)
    call mpistop
    !
  end subroutine instanttriad3Davg
  !
  subroutine instantspectraskewness(thefilenumb)
    !
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,   only : time,nstep,im,jm,km,ia,ja,ka
    use commarray, only : vel, rho
    use hdf5io
    use utility,   only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    integer,intent(in) :: thefilenumb
    character(len=128) :: infilename
    character(len=4) :: stepname
    real(8) :: u1mean,u2mean,rhomean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1spe,u2spe
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u11,u12,u21,u22
    real(8), allocatable, dimension(:,:) :: k1,k2
    integer :: allkmax
    real(8) :: u11c,u11s,u22c,u22s,thetac,thetas,omegas,psi2theta,phi2theta,thetamax
    character(len=128) :: outfilename
    integer :: hand_a,hand_b
    character(len=1) :: modeio
    integer :: i,j,n
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_u11, c_u12, c_u21, c_u22
    !
    call readinput
    !
    modeio='h'
    !
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    allkmax=ceiling(sqrt(2.d0)/3*min(ia,ja))
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,"allkmax:",allkmax
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km))
    !
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    ! Calculate favre average
    u1mean = 0.0d0
    u2mean = 0.0d0
    rhomean = 0.0d0
    !
    do i=1,im
      do j=1,jm
        u1mean = u1mean + rho(i,j,0) * vel(i,j,0,1)
        u2mean = u2mean + rho(i,j,0) * vel(i,j,0,2)
        rhomean = rhomean + rho(i,j,0)
      enddo
    enddo
    rhomean = psum(rhomean) / (1.0d0*ia*ja)
    u1mean = psum(u1mean) / (1.d0*ia*ja*rhomean)
    u2mean = psum(u2mean) / (1.d0*ia*ja*rhomean)
    !
    if(mpirank==0) print *, 'u1mean=',u1mean, 'u2mean=',u2mean
    !
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw])
    !
    ! planning
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    !!!! Convert to complex
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=CMPLX(vel(i,j,0,1)-u1mean,0.d0,C_INTPTR_T);
      u2spe(i,j)=CMPLX(vel(i,j,0,2)-u2mean,0.d0,C_INTPTR_T);
      !
    end do
    end do
    !
    !!!! Do 2d FFT
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    !
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=u1spe(i,j)/(1.d0*ia*ja)
      u2spe(i,j)=u2spe(i,j)/(1.d0*ia*ja)
      !
    end do
    end do
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    do j=1,jm
    do i=1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if((j+j0) <= (ja/2+1)) then
        k2(i,j) = real(j+j0-1,8)
      else if((j+j0)<=(ja)) then
        k2(i,j) = real(j+j0-ja-1,8)
      else
        print *,"Error, no wave number possible, (j+jm) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    !
    c_u11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u11, u11, [imfftw,jmfftw])
    c_u21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u21, u21, [imfftw,jmfftw])
    c_u12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u12, u12, [imfftw,jmfftw])
    c_u22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u22, u22, [imfftw,jmfftw])
    !
    !
    do j=1,jm
      do i=1,im
        u11(i,j) = u1spe(i,j) * CMPLX(0.d0,k1(i,j),C_INTPTR_T)
        u12(i,j) = u1spe(i,j) * CMPLX(0.d0,k2(i,j),C_INTPTR_T)
        u21(i,j) = u2spe(i,j) * CMPLX(0.d0,k1(i,j),C_INTPTR_T)
        u22(i,j) = u2spe(i,j) * CMPLX(0.d0,k2(i,j),C_INTPTR_T)
      end do
    end do
    !
    !!!! Do inverse FFT
    call fftw_mpi_execute_dft(backward_plan,u11,u11)
    call fftw_mpi_execute_dft(backward_plan,u21,u21)
    call fftw_mpi_execute_dft(backward_plan,u12,u12)
    call fftw_mpi_execute_dft(backward_plan,u22,u22)
    !
    u11c = 0.d0
    u11s = 0.d0
    u22c = 0.d0
    u22s = 0.d0
    !
    thetac = 0.d0
    thetas = 0.d0
    omegas = 0.d0
    psi2theta = 0.d0
    phi2theta = 0.d0
    thetamax = 0.d0
    !
    do j=1,jm
      do i=1,im
        u11c = u11c + real(u11(i,j))**3
        u11s = u11s + real(u11(i,j))**2
        u22c = u22c + real(u22(i,j))**3
        u22s = u22s + real(u22(i,j))**2
        !
        thetac = thetac + real(u11(i,j)+u22(i,j))**3
        thetas = thetas + real(u11(i,j)+u22(i,j))**2
        omegas = omegas + real(u21(i,j)-u12(i,j))**2
        thetamax = max(thetamax,abs(u11(i,j)+u22(i,j)))
        !
        psi2theta = psi2theta + real(u12(i,j)+u21(i,j))**2*real(u11(i,j)+u22(i,j))
        phi2theta = phi2theta + real(u11(i,j)-u22(i,j))**2*real(u11(i,j)+u22(i,j))
        !
      end do
    end do
    !
    u11c      = psum(u11c)     /(1.d0*ia*ja)
    u11s      = psum(u11s)     /(1.d0*ia*ja)
    u22c      = psum(u22c)     /(1.d0*ia*ja)
    u22s      = psum(u22s)     /(1.d0*ia*ja)
    thetac    = psum(thetac)   /(1.d0*ia*ja)
    thetas    = psum(thetas)   /(1.d0*ia*ja)
    omegas    = psum(omegas)   /(1.d0*ia*ja)
    psi2theta = psum(psi2theta)/(1.d0*ia*ja)
    phi2theta = psum(phi2theta)/(1.d0*ia*ja)
    thetamax  = pmax(thetamax)
    !
    if(mpirank == 0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/SkewnessSpec'//stepname//'.dat'
      else
        outfilename = 'pp/SkewnessSpec.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
                    firstline='nstep time u11c u11s s1 u22c u22s s2 ullc ulls skew thetamax')
      call listwrite(hand_b,u11c,u11s,u11c/sqrt(u11s**3),&
              u22c,u22s,u22c/sqrt(u22s**3),&
              u11c+u22c,u11s+u22s,(u11c+u22c)/2.d0/sqrt((u11s/2.d0+u22s/2.d0)**3),thetamax)
      print*,' <<< pp/SkewnessSpec'//stepname//'.dat ... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/SkewmomSpec'//stepname//'.dat'
      else
        outfilename = 'pp/SkewmomSpec.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
                    firstline='nstep time thetac thetas omegas psi2theta phi2theta Sa Sb')
      call listwrite(hand_b,thetac,thetas,omegas,psi2theta,phi2theta,&
                    5.d0/16.d0*thetac/sqrt((omegas/8.d0+3.d0*thetas/8.d0)**3),&
                    (thetac/8.d0+3.d0*psi2theta/16.d0+3.d0*phi2theta/16.d0)/sqrt((omegas/8.d0+3.d0*thetas/8.d0)**3))
      print*,' <<< pp/SkewmomSpec'//stepname//'.dat ... done.'
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_u11)
    call fftw_free(c_u12)
    call fftw_free(c_u21)
    call fftw_free(c_u22)
    call mpistop
    !
    deallocate(k1,k2)
    !
  end subroutine instantspectraskewness
  !
  subroutine instantPtheta2D(thefilenumb)
 
    use singleton
    use readwrite, only : readinput
    use commvar,only : time,nstep,im,jm,km,hm,ia,ja,ka
    use commarray, only : x,vel,dvel
    use hdf5io
    use parallel,  only : dataswap, mpisizedis,parapp,parallelini,mpistop,mpirank
    use comsolver, only : solvrinit,grad
    use solver,    only : refcal
    use geom,      only : geomcal
    use gridgeneration
    use fludyna, only : miucal
    use utility, only : listinit, listwrite
    ! arguments
    integer,intent(in) :: thefilenumb
    character(len=128) :: infilename
    character(len=4) :: stepname
    character(len=128) :: outfilename
    character(len=1) :: modeio
    integer :: i,j
    integer :: hand_a, hand_b
    !
    real(8), allocatable, dimension(:,:,:) :: omega,theta,psi,phi,prstheta,tmp,rho
    real(8):: sumprstheta,disspc, dissps,miu,omegarms,rhorms,Ek,rho_0
    call readinput
    !
    call mpisizedis
    if(mpirank == 0) then
      print*, '** mpisizedis done!'
    endif
    !
    call parapp
    if(mpirank == 0) then
      print*, '** parapp done!'
    endif
    !
    call parallelini
    if(mpirank == 0) then
      print*, '** parallelini done!'
    endif
    !
    call refcal
    if(mpirank == 0) then
      print*, '** refcal done!'
    endif
    !
    !
    modeio='h'
    !
    if(mpirank == 0) then
      print *, "ia:",ia,",ja:",ja,'ka:',ka
    endif 
    !
    allocate(x(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3) )
    allocate(vel(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3))
    allocate(dvel(0:im,0:jm,0:0,1:3,1:2))
    allocate(omega(1:im,1:jm,0:0),theta(1:im,1:jm,0:0),psi(1:im,1:jm,0:0),phi(1:im,1:jm,0:0))
    allocate(prs(-hm:im+hm,-hm:jm+hm,-hm:hm))
    allocate(tmp(-hm:im+hm,-hm:jm+hm,-hm:hm))
    allocate(rho(-hm:im+hm,-hm:jm+hm,-hm:hm))
    allocate(prstheta(1:im,i:jm,0:0))
    
    !
    call gridcube(2.d0*pi,2.d0*pi,0.d0)
    !
    call geomcal
    !
    write(stepname,'(i4.4)')thefilenumb
    infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='p', var=prs(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='t',var=tmp(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='ro',var=rho(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    vel(0:im,0:jm,0:km,3) = 0.d0
    sumprstheta = 0.d0
    disspc = 0.d0
    dissps = 0.d0
    Etheta = 0.d0
    Ek = 0.d0
    omegarms = 0.d0
    rhorms = 0.d0
    rho_0 = 1.d0
    !
    !
    if(mpirank == 0) then
      print *, "Swap velocity"
    endif
    !
    call dataswap(vel)
    !
    call solvrinit
    !
    if(mpirank == 0) then
      print *, "Calculate gradient"
    endif
    dvel(:,:,:,:,1)=grad(vel(:,:,:,1))
    dvel(:,:,:,:,2)=grad(vel(:,:,:,2))
    !
    omega(:,:,0) = dvel(:,:,:,2,1) - dvel(:,:,:,1,2)
    theta(:,:,0) = dvel(:,:,:,1,1) + dvel(:,:,:,2,2)
    prstheta(:,:,0) = prs(:,:,0) * theta(:,:,0)
    do i=1,im
    do j=1,jm

    miu = miucal(tmp(i,j,0))/reynolds
    disspc = disspc + 4/3 * miu * theta(i,j,0)**2
    dissps = dissps + miu * abs(omega(i,j,0))**2
    omegarms = omegarms + abs(omega(i,j,0))**2
    rhorms = rhorms + (rho(i,j,0)-rho_0)**2
    !
    enddo
    enddo
    
    sumprstheta = psum(sumprstheta)/(ia*ja*1.d0)
    disspc = psum(disspc)/(ia*ja*1.d0)
    dissps = psum(dissps)/(ia*ja*1.d0)
    rhorms = psum(rhorms)/(ia*ja*1.d0)
    omegarms = psum(omegarms)/(ia*ja*1.d0)
    if (mpirank == 0) then
      outfilename = 'pp/velgradnew'//stepname//'.dat'
      call listinit(filename=outfilename,handle=hand_a, &
                        firstline='nstep time sumprstheta disspc dissps omegarms rhorms')
      call listwrite(hand_a,sumprstheta,disspc,dissps,omegarms, rhorms)

      print*,' <<< pp/velgradnew'//stepname//'.dat ... done.'
    endif

    call mpistop
    deallocate(x,vel,dvel,prs,tmp,rho,omega,theta,psi,phi,prstheta)
  end subroutine instantPtheta2D

  !
  subroutine instantvelgradient(thefilenumb)
    !
    !
    use singleton
    use readwrite, only : readinput
    use commvar,only : time,nstep,im,jm,km,hm,ia,ja,ka
    use commarray, only : x,vel,dvel
    use hdf5io
    use parallel,  only : dataswap, mpisizedis,parapp,parallelini,mpistop,mpirank
    use comsolver, only : solvrinit,grad
    use solver,    only : refcal
    use geom,      only : geomcal
    use gridgeneration
    !
    ! arguments
    integer,intent(in) :: thefilenumb
    character(len=128) :: infilename
    character(len=4) :: stepname
    character(len=128) :: outfilename
    character(len=1) :: modeio
    integer :: i,j
    !
    real(8), allocatable, dimension(:,:,:) :: omega,theta,psi,phi
    real(8), allocatable, dimension(:,:,:) :: alpha,beta,s12,s13,s23,omega1,omega2,omega3
    !
    call readinput
    !
    call mpisizedis
    if(mpirank==0)then
      print*, '** mpisizedis done!'
    endif
    !
    call parapp
    if(mpirank==0)then
      print*, '** parapp done!'
    endif
    !
    call parallelini
    if(mpirank==0)then
      print*, '** parallelini done!'
    endif
    !
    call refcal
    if(mpirank==0)then
      print*, '** refcal done!'
    endif
    !
    modeio='h'
    !
    if(mpirank==0)then
      if(ka==0)then
        print *,"2D, ia:",ia,",ja:",ja
      else
        print *,"3D, ia:",ia,",ja:",ja, ",ka:", ka
      endif
    endif
    !
    allocate(x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    allocate(vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3))
    allocate(dvel(0:im,0:jm,0:km,1:3,1:3))
    !
    if(ka==0)then
      call gridsquare(2.d0*pi,2.d0*pi)
      allocate(omega(0:im,0:jm,0:0),theta(0:im,0:jm,0:0),psi(0:im,0:jm,0:0),phi(0:im,0:jm,0:0))
    else
      call gridcube(2.d0*pi,2.d0*pi,2.d0*pi)
      allocate(theta(0:im,0:jm,0:km),alpha(0:im,0:jm,0:km),beta(0:im,0:jm,0:km), &
                s12(0:im,0:jm,0:km),s13(0:im,0:jm,0:km),s23(0:im,0:jm,0:km),     &
                omega1(0:im,0:jm,0:km),omega2(0:im,0:jm,0:km),omega3(0:im,0:jm,0:km))
    endif
    !
    call geomcal
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    if(mpirank==0)then
      print *, "Swap velocity"
    endif
    !
    call dataswap(vel)
    !
    call solvrinit
    !
    if(mpirank==0)then
      print *, "Calculate gradient"
    endif
    !
    dvel(:,:,:,1,:)=grad(vel(:,:,:,1))
    dvel(:,:,:,2,:)=grad(vel(:,:,:,2))
    if(ka==0)then
      dvel(:,:,:,3,:)=grad(vel(:,:,:,3))
    endif
    !
    if (thefilenumb .ne. 0) then
      outfilename = 'pp/velgrad'//stepname//'.'//modeio//'5'
    else
      outfilename = 'pp/velgrad.'//modeio//'5'
    endif
    
    ! !
    if(ka==0)then
      !
      theta(:,:,:) = dvel(:,:,:,1,1) + dvel(:,:,:,2,2)
      omega(:,:,:) = dvel(:,:,:,2,1) - dvel(:,:,:,1,2)
      psi(:,:,:)   = dvel(:,:,:,2,1) + dvel(:,:,:,1,2)
      phi(:,:,:)   = dvel(:,:,:,1,1) - dvel(:,:,:,2,2)
      !
      call h5io_init(trim(outfilename),mode='write')
      !
      call h5write(varname='m11',var=dvel(1:im,1:jm,0:0,1,1),mode=modeio)
      call h5write(varname='m12',var=dvel(1:im,1:jm,0:0,1,2),mode=modeio)
      call h5write(varname='m21',var=dvel(1:im,1:jm,0:0,2,1),mode=modeio)
      call h5write(varname='m22',var=dvel(1:im,1:jm,0:0,2,2),mode=modeio)
      !
      call h5write(varname='theta',var=theta(1:im,1:jm,0:0),mode=modeio)
      call h5write(varname='omega',var=omega(1:im,1:jm,0:0),mode=modeio)
      call h5write(varname='psi',var=psi(1:im,1:jm,0:0),mode=modeio)
      call h5write(varname='phi',var=phi(1:im,1:jm,0:0),mode=modeio)
      !
      call h5io_end
      !
      if(mpirank==0) then
        call h5srite(varname='time',var=time,filename=trim(outfilename))
        call h5srite(varname='nstep',var=nstep,filename=trim(outfilename))
      endif
      !
      !
    else
      !
      theta(:,:,:) = dvel(:,:,:,1,1) + dvel(:,:,:,2,2) + dvel(:,:,:,3,3)
      alpha(:,:,:) = dvel(:,:,:,1,1) - dvel(:,:,:,2,2)
      beta(:,:,:)  = dvel(:,:,:,1,1) - dvel(:,:,:,3,3)
      s12(:,:,:)   = dvel(:,:,:,1,2) + dvel(:,:,:,2,1)
      s13(:,:,:)   = dvel(:,:,:,1,3) + dvel(:,:,:,3,1)
      s23(:,:,:)   = dvel(:,:,:,2,3) + dvel(:,:,:,3,2)
      omega1(:,:,:)= dvel(:,:,:,3,2) - dvel(:,:,:,2,3)
      omega2(:,:,:)= dvel(:,:,:,1,3) - dvel(:,:,:,3,1)
      omega3(:,:,:)= dvel(:,:,:,2,1) - dvel(:,:,:,1,2)
      !
      call h5io_init(trim(outfilename),mode='write')
      !
      call h5write(varname='m11',var=dvel(1:im,1:jm,1:km,1,1),mode=modeio)
      call h5write(varname='m12',var=dvel(1:im,1:jm,1:km,1,2),mode=modeio)
      call h5write(varname='m13',var=dvel(1:im,1:jm,1:km,1,3),mode=modeio)
      call h5write(varname='m21',var=dvel(1:im,1:jm,1:km,2,1),mode=modeio)
      call h5write(varname='m22',var=dvel(1:im,1:jm,1:km,2,2),mode=modeio)
      call h5write(varname='m23',var=dvel(1:im,1:jm,1:km,2,3),mode=modeio)
      call h5write(varname='m31',var=dvel(1:im,1:jm,1:km,3,1),mode=modeio)
      call h5write(varname='m32',var=dvel(1:im,1:jm,1:km,3,2),mode=modeio)
      call h5write(varname='m33',var=dvel(1:im,1:jm,1:km,3,3),mode=modeio)
      !
      call h5write(varname='theta',var=theta(1:im,1:jm,1:km),mode=modeio)
      call h5write(varname='alpha',var=alpha(1:im,1:jm,1:km),mode=modeio)
      call h5write(varname='beta',var=beta(1:im,1:jm,1:km),mode=modeio)
      call h5write(varname='s12',var=s12(1:im,1:jm,1:km),mode=modeio)
      call h5write(varname='s23',var=s23(1:im,1:jm,1:km),mode=modeio)
      call h5write(varname='s13',var=s13(1:im,1:jm,1:km),mode=modeio)
      call h5write(varname='omega1',var=omega1(1:im,1:jm,1:km),mode=modeio)
      call h5write(varname='omega2',var=omega2(1:im,1:jm,1:km),mode=modeio)
      call h5write(varname='omega3',var=omega3(1:im,1:jm,1:km),mode=modeio)
      !
      !
      call h5io_end
      !
      if(mpirank==0) then
        call h5srite(varname='time',var=time,filename=trim(outfilename))
        call h5srite(varname='nstep',var=nstep,filename=trim(outfilename))
      endif
      !
    endif
    !
    call mpistop
    !
    deallocate(x,vel,dvel)
    if(ka==0)then
      deallocate(omega,theta,psi,phi)
    else
      deallocate(theta,alpha,beta,s12,s13,s23,omega1,omega2,omega3)
    endif
    !
  end subroutine instantvelgradient
  !
  subroutine velgradient_scale_lengths(thefilenumb)
    !
    !
    use singleton
    use readwrite, only : readinput
    use commvar,   only : time,nstep,im,jm,km,hm,ia,ja,ka,Reynolds
    use commarray, only : x,vel,dvel,tmp,rho
    use hdf5io
    use parallel,  only : dataswap, mpisizedis,parapp,parallelini,mpistop,mpirank,psum
    use fludyna,   only : miucal
    use comsolver, only : solvrinit,grad
    use solver,    only : refcal
    use geom,      only : geomcal
    use gridgeneration
    !
    ! arguments
    integer,intent(in) :: thefilenumb
    character(len=128) :: infilename
    character(len=4) :: stepname
    character(len=128) :: outfilename
    character(len=1) :: modeio
    !
    real(8), allocatable, dimension(:,:,:) :: omega,theta,psi,phi,miu
    real(8), allocatable, dimension(:,:,:) :: alpha,beta,s12,s13,s23,omega1,omega2,omega3
    real(8), allocatable, dimension(:,:,:) :: ds, eta, domega
    real(8) :: etaavg,nuavg,disspavg, rsamples
    integer ::i,j,k
    !
    call readinput
    !
    call mpisizedis
    if(mpirank==0)then
      print*, '** mpisizedis done!'
    endif
    !
    call parapp
    if(mpirank==0)then
      print*, '** parapp done!'
    endif
    !
    call parallelini
    if(mpirank==0)then
      print*, '** parallelini done!'
    endif
    !
    call refcal
    if(mpirank==0)then
      print*, '** refcal done!'
    endif
    !
    modeio='h'
    !
    if(mpirank==0)then
      if(ka==0)then
        print *,"2D, ia:",ia,",ja:",ja
      else
        print *,"3D, ia:",ia,",ja:",ja, ",ka:", ka
      endif
    endif
    !
    allocate(x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    allocate(vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),tmp(0:im,0:jm,0:km),rho(0:im,0:jm,0:km),miu(0:im,0:jm,0:km))
    allocate(dvel(0:im,0:jm,0:km,1:3,1:3))
    !
    if(ka==0)then
      call gridsquare(2.d0*pi,2.d0*pi)
      allocate(omega(0:im,0:jm,0:0),theta(0:im,0:jm,0:0),psi(0:im,0:jm,0:0),phi(0:im,0:jm,0:0))
      allocate(ds(0:im,0:jm,0:0),eta(0:im,0:jm,0:0),domega(0:im,0:jm,0:0))
    else
      call gridcube(2.d0*pi,2.d0*pi,2.d0*pi)
      allocate(theta(0:im,0:jm,0:km),alpha(0:im,0:jm,0:km),beta(0:im,0:jm,0:km), &
                s12(0:im,0:jm,0:km),s13(0:im,0:jm,0:km),s23(0:im,0:jm,0:km),     &
                omega1(0:im,0:jm,0:km),omega2(0:im,0:jm,0:km),omega3(0:im,0:jm,0:km))
      allocate(ds(0:im,0:jm,0:km),eta(0:im,0:jm,0:km),domega(0:im,0:jm,0:km))
    endif
    !
    call geomcal
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='t',  var=tmp(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='ro',var=rho(0:im,0:jm,0:km), mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    if(mpirank==0)then
      print *, "Swap velocity"
    endif
    !
    call dataswap(vel)
    !
    call solvrinit
    !
    if(mpirank==0)then
      print *, "Calculate gradient"
    endif
    !
    dvel(:,:,:,1,:)=grad(vel(:,:,:,1))
    dvel(:,:,:,2,:)=grad(vel(:,:,:,2))
    if(ka==0)then
      dvel(:,:,:,3,:)=grad(vel(:,:,:,3))
    endif
    !
    if (thefilenumb .ne. 0) then
      outfilename = 'pp/velgrad_scale'//stepname//'.'//modeio//'5'
    else
      outfilename = 'pp/velgrad_scale.'//modeio//'5'
    endif
    
    ! !
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      miu(i,j,k) = miucal(tmp(i,j,k))/Reynolds
    enddo
    enddo
    enddo
    !
    if(ka==0)then
      !
      theta(:,:,0) = dvel(:,:,0,1,1) + dvel(:,:,0,2,2)
      omega(:,:,0) = dvel(:,:,0,2,1) - dvel(:,:,0,1,2)
      ds = sqrt(miu/rho/abs(theta))
      eta = sqrt(miu/rho/sqrt(theta**2+omega**2))
      domega = sqrt(miu/rho/abs(omega))
      !
      nuavg = 0.d0
      disspavg = 0.d0
      !
      k = 0
      do i=1,im
      do j=1,jm
        nuavg = nuavg + miu(i,j,k)/rho(i,j,k)
        disspavg = disspavg + miu(i,j,k)/rho(i,j,k) * (theta(i,j,k)**2 + omega(i,j,k)**2)
      enddo
      enddo
      !
      rsamples = dble(ia*ja)
      nuavg = psum(nuavg)/rsamples
      disspavg = psum(disspavg)/rsamples
      !
      etaavg = sqrt(sqrt(nuavg**3/disspavg))
      !
      call h5io_init(trim(outfilename),mode='write')
      !
      call h5write(varname='ds',var=ds(1:im,1:jm,0:0),mode=modeio)
      call h5write(varname='eta',var=eta(1:im,1:jm,0:0),mode=modeio)
      call h5write(varname='domega',var=domega(1:im,1:jm,0:0),mode=modeio)
      !
      call h5io_end
      !
      if(mpirank==0) then
        call h5srite(varname='time',var=time,filename=trim(outfilename))
        call h5srite(varname='nstep',var=nstep,filename=trim(outfilename))
        call h5srite(varname='etaavg',var=etaavg,filename=trim(outfilename))
      endif
      !
      !
    else
      !
      theta(:,:,:) = dvel(:,:,:,1,1) + dvel(:,:,:,2,2) + dvel(:,:,:,3,3)
      omega1(:,:,:)= dvel(:,:,:,3,2) - dvel(:,:,:,2,3)
      omega2(:,:,:)= dvel(:,:,:,1,3) - dvel(:,:,:,3,1)
      omega3(:,:,:)= dvel(:,:,:,2,1) - dvel(:,:,:,1,2)
      !
      ds = sqrt(miu/rho/abs(theta))
      eta = sqrt(miu/rho/sqrt(theta**2+omega1**2+omega2**2+omega3**2))
      domega = sqrt(miu/rho/sqrt(omega1**2+omega2**2+omega3**2))
      !
      nuavg = 0.d0
      disspavg = 0.d0
      !
      do i=1,im
      do j=1,jm
      do k=1,km
        nuavg = nuavg + miu(i,j,k)/rho(i,j,k)
        disspavg = disspavg + miu(i,j,k)/rho(i,j,k) * (theta(i,j,k)**2 + omega1(i,j,k)**2 + omega2(i,j,k)**2 + omega3(i,j,k)**2)
      enddo
      enddo
      enddo
      !
      rsamples = dble(ia*ja*ka)
      nuavg = psum(nuavg)/rsamples
      disspavg = psum(disspavg)/rsamples
      !
      etaavg = sqrt(sqrt(nuavg**3/disspavg))
      !
      call h5io_init(trim(outfilename),mode='write')
      !
      call h5write(varname='ds',var=ds(1:im,1:jm,1:km),mode=modeio)
      call h5write(varname='eta',var=eta(1:im,1:jm,1:km),mode=modeio)
      call h5write(varname='domega',var=domega(1:im,1:jm,1:km),mode=modeio)
      !
      call h5io_end
      !
      if(mpirank==0) then
        call h5srite(varname='time',var=time,filename=trim(outfilename))
        call h5srite(varname='nstep',var=nstep,filename=trim(outfilename))
        call h5srite(varname='etaavg',var=etaavg,filename=trim(outfilename))
      endif
      !
    endif
    !
    call mpistop
    !
    deallocate(x,vel,dvel,tmp,rho,miu)
    if(ka==0)then
      deallocate(omega,theta,psi,phi)
    else
      deallocate(theta,alpha,beta,s12,s13,s23,omega1,omega2,omega3)
    endif
    deallocate(ds,eta,domega)
    !
  end subroutine velgradient_scale_lengths
  !
  subroutine SGSPiOmgea2D(thefilenumb)
    ! 
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km,reynolds
    use commarray, only: vel,tmp,rho
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
    use fludyna,   only : miucal
    use solver,    only : refcal
    include 'fftw3-mpi.f03'
    !
    integer,intent(in) :: thefilenumb
    integer :: fh
    integer :: i,j,m,n
    character(len=128) :: infilename,outfilename,outfilename2
    character(len=4) :: stepname,mname
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1,u2,w,wu1,wu2
    real(8), allocatable, dimension(:,:) :: k1,k2
    complex(8) :: imag
    real(8),allocatable,dimension(:) :: l_lim
    integer :: num_l,num_alpha,num_alphamin
    integer :: hand_a,hand_b
    real(8) :: l_min, ratio_max, ratio_min
    real(8) :: Gl, beta,roav,miu,miuav
    real(8), allocatable, dimension(:) :: Pi_omega
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1_filted,u2_filted,w_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: wu1_filted,wu2_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: wx1_filted,wx2_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: wx1,wx2
    !
    !
    type(C_PTR) :: c_u1,c_u2,c_w,forward_plan,backward_plan
    type(C_PTR) :: c_wx1,c_wx2,c_wu1,c_wu2
    type(C_PTR) :: c_u1_filted,c_u2_filted,c_w_filted
    type(C_PTR) :: c_wu1_filted,c_wu2_filted
    type(C_PTR) :: c_wx1_filted,c_wx2_filted
    !
    integer,dimension(8) :: value
    character(len=1) :: modeio
    logical :: loutput
    !
    call readinput
    call refcal
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja, 'Reynolds:',reynolds
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !!!! Read velocity and density field
    allocate(vel(0:im,0:jm,0:km,1:2),rho(0:im,0:jm,0:km),tmp(0:im,0:jm,0:km))
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='t', var=tmp(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    !! wavenumber
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    do j = 1,jm
    do i = 1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if((j+j0) <= (ja/2+1)) then
        k2(i,j) = real(j+j0-1,8)
      else if((j+j0)<=(ja)) then
        k2(i,j) = real(j+j0-ja-1,8)
      else
        print *,"Error, no wave number possible, (j+j0) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    !
    !! Imaginary number prepare
    imag = CMPLX(0.d0,1.d0,8)
    !
    !!!! Prepare initial field in Fourier space
    !! velocity
    c_u1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1, u1, [imfftw,jmfftw])
    c_u2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2, u2, [imfftw,jmfftw])
    c_w = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w, w, [imfftw,jmfftw])
    c_wu1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_wu1, wu1, [imfftw,jmfftw])
    c_wu2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_wu2, wu2, [imfftw,jmfftw])
    c_wx1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_wx1, wx1, [imfftw,jmfftw])
    c_wx2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_wx2, wx2, [imfftw,jmfftw])
    !
    !
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1,u1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1,u1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    do j=1,jm
    do i=1,im
      !
      u1(i,j)=CMPLX(vel(i,j,0,1),0.d0,C_INTPTR_T);
      u2(i,j)=CMPLX(vel(i,j,0,2),0.d0,C_INTPTR_T);
      !
    end do
    end do
    !
    !After this bloc, u1,u2,w are in spectral space
    call fftw_mpi_execute_dft(forward_plan,u1,u1)
    call fftw_mpi_execute_dft(forward_plan,u2,u2)
    !
    do j=1,jm
    do i=1,im
      !
      u1(i,j)=u1(i,j)/(1.d0*ia*ja)
      u2(i,j)=u2(i,j)/(1.d0*ia*ja)
      !
      w(i,j)=imag*k1(i,j)*u2(i,j)-imag*k2(i,j)*u1(i,j)
      wx1(i,j)=imag*k1(i,j)*w(i,j)
      wx2(i,j)=imag*k2(i,j)*w(i,j)
      !
    end do
    end do
    !
    !After this bloc,u1,u2,w,wx1,wx2 are in physical space
    call fftw_mpi_execute_dft(backward_plan,u1,u1)
    call fftw_mpi_execute_dft(backward_plan,u2,u2)
    call fftw_mpi_execute_dft(backward_plan,w,w)
    call fftw_mpi_execute_dft(backward_plan,wx1,wx1)
    call fftw_mpi_execute_dft(backward_plan,wx2,wx2)
    !
    beta = 0.d0
    roav = 0.d0
    miu = 0.d0
    miuav = 0.d0
    !!! 
    do j=1,jm
    do i=1,im
      !
      !
      wu1(i,j)=w(i,j)*u1(i,j)
      wu2(i,j)=w(i,j)*u2(i,j)
      !
      !
      miu=miucal(tmp(i,j,0))/reynolds
      beta = beta + miu/rho(i,j,0) * (wx1(i,j)**2 + wx2(i,j)**2)
      !
      roav=roav+rho(i,j,0)
      miuav=miuav+miu
      !
    end do
    end do
    !
    beta  = psum(beta) / (ia*ja)
    roav  = psum(roav) / (ia*ja)
    miuav = psum(miuav)/ (ia*ja)
    !
    !After this bloc, u1,u2,w,wu1,wu2 are in spectral space
    call fftw_mpi_execute_dft(forward_plan,u1,u1)
    call fftw_mpi_execute_dft(forward_plan,u2,u2)
    call fftw_mpi_execute_dft(forward_plan,w,w)
    call fftw_mpi_execute_dft(forward_plan,wu1,wu1)
    call fftw_mpi_execute_dft(forward_plan,wu2,wu2)
    !
    do j=1,jm
    do i=1,im
      !
      u1(i,j)=u1(i,j)/(1.d0*ia*ja)
      u2(i,j)=u2(i,j)/(1.d0*ia*ja)
      w(i,j)=w(i,j)/(1.d0*ia*ja)
      wu1(i,j)=wu1(i,j)/(1.d0*ia*ja)
      wu2(i,j)=wu2(i,j)/(1.d0*ia*ja)
      !
    end do
    end do
    !
    !
    if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
    !
    !!!! Prepare l,alpha and others
    call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
    l_min = 2*pi/im
    allocate(l_lim(1:num_l))
    !
    do i=1,num_l
      l_lim(i) = exp(log(ratio_min)+(i-1) * (log(ratio_max)-log(ratio_min)) / (num_l-1)) * l_min
    enddo
    !
    if(mpirank==0)  print *, "Integrate point allocated"
    !
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !!!!
    allocate(Pi_omega(1:num_l))
    !
    c_u1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1_filted, u1_filted, [imfftw,jmfftw])
    c_u2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2_filted, u2_filted, [imfftw,jmfftw])
    c_w_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w_filted, w_filted, [imfftw,jmfftw])
    c_wu1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_wu1_filted, wu1_filted, [imfftw,jmfftw])
    c_wu2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_wu2_filted, wu2_filted, [imfftw,jmfftw])
    c_wx1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_wx1_filted, wx1_filted, [imfftw,jmfftw])
    c_wx2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_wx2_filted, wx2_filted, [imfftw,jmfftw])
    !
    Pi_omega = 0.d0
    !
    if(mpirank==0)  print *, "Array allocated and initialized"
    !
    do m=1,num_l
      !
      !!!!!! Filter to get Sij filted by l
      if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
      !
      !
      !!!! Velocity Favre average and density average
      ! After this bloc, u1_filted is in spectral space
      do j=1,jm
      do i=1,im
        Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) !  Filtre scale :l
        !
        u1_filted(i,j) = u1(i,j) *Gl
        u2_filted(i,j) = u2(i,j) *Gl
        w_filted(i,j)  = w(i,j)  *Gl
        wu1_filted(i,j)= wu1(i,j)*Gl
        wu2_filted(i,j)= wu2(i,j)*Gl
        wx1_filted(i,j)= w(i,j)  *Gl*imag*k1(i,j)
        wx2_filted(i,j)= w(i,j)  *Gl*imag*k2(i,j)
        !
      enddo
      enddo
      !
      ! After this bloc, u1_filted is in physical space
      call fftw_mpi_execute_dft(backward_plan,u1_filted,u1_filted)
      call fftw_mpi_execute_dft(backward_plan,u2_filted,u2_filted)
      call fftw_mpi_execute_dft(backward_plan,w_filted,w_filted)
      call fftw_mpi_execute_dft(backward_plan,wu1_filted,wu1_filted)
      call fftw_mpi_execute_dft(backward_plan,wu2_filted,wu2_filted)
      call fftw_mpi_execute_dft(backward_plan,wx1_filted,wx1_filted)
      call fftw_mpi_execute_dft(backward_plan,wx2_filted,wx2_filted)
      !
      ! 
      do j=1,jm
      do i=1,im
        !
        Pi_omega(m) = Pi_omega(m) + &
            dreal(wu1_filted(i,j) - w_filted(i,j) * u1_filted(i,j))*dreal(wx1_filted(i,j)) + &
            dreal(wu2_filted(i,j) - w_filted(i,j) * u2_filted(i,j))*dreal(wx2_filted(i,j)) 
        !
      enddo
      enddo
      !
      !
      if(mpirank==0)  print *, '** l filted!'
      !
      Pi_omega(m) =	 psum(Pi_omega(m)) / (ia*ja)
      !
      !
    enddo
    !
    !
    if(mpirank==0)  print *, 'Job finish'
    !
    if(mpirank==0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/SGS_Piomega_'//stepname//'.dat'
      else
        outfilename = 'pp/SGS_Piomega.dat'
      endif
      
      call listinit(filename=outfilename,handle=hand_a, &
                    firstline='nstep time lOlmin piomega beta nuav lens')
      do m=1,num_l
        call listwrite(hand_a,l_lim(m)/l_min, Pi_omega(m), beta,miuav/roav,((miuav/roav)**3/beta)**(1.d0/6.d0))
      enddo
      !
      print *, '>>>>', outfilename
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1)
    call fftw_free(c_u2)
    call fftw_free(c_w)
    call fftw_free(c_wx1)
    call fftw_free(c_wx2)
    call fftw_free(c_wu1)
    call fftw_free(c_wu2)
    call fftw_free(c_u1_filted)
    call fftw_free(c_u2_filted)
    call fftw_free(c_w_filted)
    call fftw_free(c_wu1_filted)
    call fftw_free(c_wu2_filted)
    call fftw_free(c_wx1_filted)
    call fftw_free(c_wx2_filted)
    call mpistop
    deallocate(k1,k2)
    deallocate(l_lim)
    deallocate(Pi_omega)
    !
  end subroutine SGSPiOmgea2D
  !
  subroutine SGSPi2Dtot(thefilenumb)
    ! 
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km
    use commarray, only: vel, rho
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
    include 'fftw3-mpi.f03'
    !
    integer,intent(in) :: thefilenumb
    integer :: fh
    integer :: i,j,m,n
    character(len=128) :: infilename,outfilename,outfilename2
    character(len=4) :: stepname,mname
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1,w2,rhocom
    real(8), allocatable, dimension(:,:) :: k1,k2
    complex(8) :: imag
    real(8),allocatable,dimension(:) :: l_lim
    integer :: num_l,num_alpha,num_alphamin
    integer :: hand_a,hand_b
    real(8) :: l_min, ratio_max, ratio_min
    real(8) :: Gl
    real(8), allocatable, dimension(:) :: Pi_tot
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1_filted,w2_filted,rho_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1w1,w1w2,w2w1,w2w2
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1w1_filted,w1w2_filted,w2w1_filted,w2w2_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: A11_filted,A12_filted,A21_filted,A22_filted
    real(8), allocatable, dimension(:,:,:,:) :: tau ! 1:2,1:2,1:im,1:jm
    !
    !
    type(C_PTR) :: c_w1,c_w2,c_rhocom,forward_plan,backward_plan
    type(C_PTR) :: c_w1_filted,c_w2_filted,c_rho_filted
    type(C_PTR) :: c_w1w1,c_w1w2,c_w2w1,c_w2w2
    type(C_PTR) :: c_w1w1_filted,c_w1w2_filted,c_w2w1_filted,c_w2w2_filted
    type(C_PTR) :: c_A11_filted,c_A12_filted,c_A21_filted,c_A22_filted
    !
    integer,dimension(8) :: value
    character(len=1) :: modeio
    logical :: loutput
    !
    call readinput
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !!!! Read velocity and density field
    allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km))
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    !!!! Prepare initial field in Fourier space
    !! velocity
    c_w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1, w1, [imfftw,jmfftw])
    c_w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2, w2, [imfftw,jmfftw])
    c_rhocom = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw])
    !
    c_w1w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w1, w1w1, [imfftw,jmfftw])
    c_w1w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w2, w1w2, [imfftw,jmfftw])
    c_w2w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w1, w2w1, [imfftw,jmfftw])
    c_w2w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w2, w2w2, [imfftw,jmfftw])
    !
    c_w1w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w1_filted, w1w1_filted, [imfftw,jmfftw])
    c_w1w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w2_filted, w1w2_filted, [imfftw,jmfftw])
    c_w2w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w1_filted, w2w1_filted, [imfftw,jmfftw])
    c_w2w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w2_filted, w2w2_filted, [imfftw,jmfftw])
    !
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    do j=1,jm
    do i=1,im
      !
      w1(i,j)=CMPLX(vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
      w2(i,j)=CMPLX(vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
      rhocom(i,j)=CMPLX(rho(i,j,0),0.d0,C_INTPTR_T);
      w1w1(i,j)=CMPLX(vel(i,j,0,1)*vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
      w1w2(i,j)=CMPLX(vel(i,j,0,1)*vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
      w2w1(i,j)=CMPLX(vel(i,j,0,2)*vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
      w2w2(i,j)=CMPLX(vel(i,j,0,2)*vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
      !
    end do
    end do
    !
    !After this bloc, w1 is (rho*u1) in spectral space
    call fftw_mpi_execute_dft(forward_plan,w1,w1)
    call fftw_mpi_execute_dft(forward_plan,w2,w2)
    call fftw_mpi_execute_dft(forward_plan,w1w1,w1w1)
    call fftw_mpi_execute_dft(forward_plan,w1w2,w1w2)
    call fftw_mpi_execute_dft(forward_plan,w2w1,w2w1)
    call fftw_mpi_execute_dft(forward_plan,w2w2,w2w2)
    call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
    do j=1,jm
    do i=1,im
      !
      w1(i,j)=w1(i,j)/(1.d0*ia*ja)
      w2(i,j)=w2(i,j)/(1.d0*ia*ja)
      !
      w1w1(i,j)=w1w1(i,j)/(1.d0*ia*ja)
      w1w2(i,j)=w1w2(i,j)/(1.d0*ia*ja)
      w2w1(i,j)=w2w1(i,j)/(1.d0*ia*ja)
      w2w2(i,j)=w2w2(i,j)/(1.d0*ia*ja)
      !
      rhocom(i,j)=rhocom(i,j)/(1.d0*ia*ja)
      !
    end do
    end do

    !
    !
    !! wavenumber
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    do j = 1,jm
    do i = 1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if((j+j0) <= (ja/2+1)) then
        k2(i,j) = real(j+j0-1,8)
      else if((j+j0)<=(ja)) then
        k2(i,j) = real(j+j0-ja-1,8)
      else
        print *,"Error, no wave number possible, (j+j0) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    !
    !! Imaginary number prepare
    imag = CMPLX(0.d0,1.d0,8)
    !
    allocate(tau(1:2,1:2,1:im,1:jm))
    !
    if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
    !!!! Prepare l,alpha and others
    call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
    l_min = 2*pi/im
    allocate(l_lim(1:num_l))
    !
    do i=1,num_l
      l_lim(i) = exp(log(ratio_min)+(i-1) * (log(ratio_max)-log(ratio_min)) / (num_l-1)) * l_min
    enddo
    !
    if(mpirank==0)  print *, "Integrate point allocated"
    !
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !!!!
    allocate(Pi_tot(1:num_l))
    !
    c_w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw])
    c_w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw])
    c_rho_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw])
    !
    c_A11_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A11_filted, A11_filted, [imfftw,jmfftw])
    c_A12_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A12_filted, A12_filted, [imfftw,jmfftw])
    c_A21_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A21_filted, A21_filted, [imfftw,jmfftw])
    c_A22_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A22_filted, A22_filted, [imfftw,jmfftw])
    !
    Pi_tot = 0.d0
    !
    if(mpirank==0)  print *, "Array allocated and initialized"
    !
    do m=1,num_l
      !
      !!!!!! Filter to get Sij filted by l
      if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
      !
      !
      !!!! Velocity Favre average and density average
      ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
      do j=1,jm
      do i=1,im
        Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) !  Filtre scale :l
        !
        w1_filted(i,j)    = w1(i,j)   *Gl
        w2_filted(i,j)    = w2(i,j)   *Gl
        !
        w1w1_filted(i,j)  = w1w1(i,j) *Gl
        w1w2_filted(i,j)  = w1w2(i,j) *Gl
        w2w1_filted(i,j)  = w2w1(i,j) *Gl
        w2w2_filted(i,j)  = w2w2(i,j) *Gl
        !
        rho_filted(i,j)   = rhocom(i,j)*Gl
      enddo
      enddo
      !
      ! After this bloc, w1_filted is (rho*u1)_filted in physical space
      call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
      call fftw_mpi_execute_dft(backward_plan,w1w1_filted,w1w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w1w2_filted,w1w2_filted)
      call fftw_mpi_execute_dft(backward_plan,w2w1_filted,w2w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w2w2_filted,w2w2_filted)
      call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
      !
      ! After this bloc, w1_filted is u1_filted in physical space
      do j=1,jm
      do i=1,im
        w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
        w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
        !
        tau(1,1,i,j) = dreal(w1w1_filted(i,j) - rho_filted(i,j) * w1_filted(i,j) * w1_filted(i,j))
        tau(1,2,i,j) = dreal(w1w2_filted(i,j) - rho_filted(i,j) * w1_filted(i,j) * w2_filted(i,j))
        tau(2,1,i,j) = dreal(w2w1_filted(i,j) - rho_filted(i,j) * w2_filted(i,j) * w1_filted(i,j))
        tau(2,2,i,j) = dreal(w2w2_filted(i,j) - rho_filted(i,j) * w2_filted(i,j) * w2_filted(i,j))
        !
      enddo
      enddo
      !
      ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
      call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
      call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
      !
      do j=1,jm
      do i=1,im
        !
        w1_filted(i,j)  = w1_filted(i,j)/(1.d0*ia*ja)
        w2_filted(i,j)  = w2_filted(i,j)/(1.d0*ia*ja)
        !
        A11_filted(i,j) = imag*w1_filted(i,j)*k1(i,j)
        A21_filted(i,j) = imag*w2_filted(i,j)*k1(i,j)
        A12_filted(i,j) = imag*w1_filted(i,j)*k2(i,j)
        A22_filted(i,j) = imag*w2_filted(i,j)*k2(i,j)
        !
      end do
      end do
      !
      !
      !
      ! After this bloc, A11_filted is A11_filted in physical space
      call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
      call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
      call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
      call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
      !
      do j=1,jm
      do i=1,im
        !
        Pi_tot(m) = Pi_tot(m) + tau(1,1,i,j) * A11_filted(i,j) + &
                                tau(1,2,i,j) * A12_filted(i,j) + &
                                tau(2,1,i,j) * A21_filted(i,j) + &
                                tau(2,2,i,j) * A22_filted(i,j)
        !
      end do
      end do
      !
      if(mpirank==0)  print *, '** l filted!'
      !
      Pi_tot(m) =	 psum(Pi_tot(m)) / (ia*ja)
      !
      !
    enddo
    if(mpirank==0)  print *, 'Job finish'
    !
    if(mpirank==0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/SGS_Pitot_'//stepname//'.dat'
      else
        outfilename = 'pp/SGS_Pitot.dat'
      endif
      
      call listinit(filename=outfilename,handle=hand_a, &
                    firstline='nstep time lOlmin pitot')
      do m=1,num_l
        call listwrite(hand_a,l_lim(m)/l_min, Pi_tot(m))
      enddo
      !
      print *, '>>>>', outfilename
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_w1)
    call fftw_free(c_w2)
    call fftw_free(c_rhocom)
    call fftw_free(c_w1_filted)
    call fftw_free(c_w2_filted)
    call fftw_free(c_w1w1)
    call fftw_free(c_w1w2)
    call fftw_free(c_w2w1)
    call fftw_free(c_w2w2)
    call fftw_free(c_w1w1_filted)
    call fftw_free(c_w1w2_filted)
    call fftw_free(c_w2w1_filted)
    call fftw_free(c_w2w2_filted)
    call fftw_free(c_rho_filted)
    call fftw_free(c_A11_filted)
    call fftw_free(c_A12_filted)
    call fftw_free(c_A21_filted)
    call fftw_free(c_A22_filted)
    call mpistop
    deallocate(k1,k2,tau)
    deallocate(l_lim)
    deallocate(Pi_tot)
    !
  end subroutine SGSPi2Dtot
  !
  subroutine SGSPi2Dlocal(thefilenumb)
    !
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km
    use commarray, only: vel, rho
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
    include 'fftw3-mpi.f03'
    !
    integer,intent(in) :: thefilenumb
    integer :: fh
    integer :: i,j,m,n
    character(len=128) :: infilename,outfilename,outfilename2
    character(len=4) :: stepname,mname
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1,w2,rhocom
    real(8), allocatable, dimension(:,:) :: k1,k2
    complex(8) :: imag
    real(8),allocatable,dimension(:) :: l_lim
    integer :: num_l,num_alpha,num_alphamin
    integer :: hand_a,hand_b
    real(8) :: l_min, ratio_max, ratio_min
    real(8) :: Gl
    real(8), allocatable, dimension(:) :: Pis1,Pis2,Pim2,Pim3,Pid
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1f,w2f,rhof
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: A11,A12,A21,A22
    !
    complex(8), allocatable, dimension(:,:) :: All
    complex(8), allocatable, dimension(:,:) :: S11,S12,S21,S22
    complex(8), allocatable, dimension(:,:) :: W12,W21
    !
    type(C_PTR) :: c_w1,c_w2,c_rhocom,forward_plan,backward_plan
    type(C_PTR) :: c_w1f,c_w2f,c_rho
    type(C_PTR) :: c_A11,c_A12,c_A21,c_A22
    !
    integer,dimension(8) :: value
    character(len=1) :: modeio
    logical :: loutput
    !
    call readinput
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !!!! Read velocity and density field
    allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km))
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    !!!! Prepare initial field in Fourier space
    !! velocity
    c_w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1, w1, [imfftw,jmfftw])
    c_w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2, w2, [imfftw,jmfftw])
    c_rhocom = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw])
    !
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    allocate(All(1:im,1:jm),&
            S11(1:im,1:jm),S12(1:im,1:jm),&
            S21(1:im,1:jm),S22(1:im,1:jm),&
            W12(1:im,1:jm),W21(1:im,1:jm))
    !
    do j=1,jm
    do i=1,im
      !
      w1(i,j)=CMPLX(vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
      w2(i,j)=CMPLX(vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
      rhocom(i,j)=CMPLX(rho(i,j,0),0.d0,C_INTPTR_T);
      !
    end do
    end do
    !
    !After this bloc, w1 is (rho*u1) in spectral space
    call fftw_mpi_execute_dft(forward_plan,w1,w1)
    call fftw_mpi_execute_dft(forward_plan,w2,w2)
    call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
    !
    do j=1,jm
    do i=1,im
      !
      w1(i,j)=w1(i,j)/(1.d0*ia*ja)
      w2(i,j)=w2(i,j)/(1.d0*ia*ja)
      !
      rhocom(i,j)=rhocom(i,j)/(1.d0*ia*ja)
      !
    end do
    end do
    !
    !
    !! wavenumber
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    do j = 1,jm
    do i = 1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if((j+j0) <= (ja/2+1)) then
        k2(i,j) = real(j+j0-1,8)
      else if((j+j0)<=(ja)) then
        k2(i,j) = real(j+j0-ja-1,8)
      else
        print *,"Error, no wave number possible, (j+j0) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    !
    !! Imaginary number prepare
    imag = CMPLX(0.d0,1.d0,8)
    !
    !
    if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
    !!!! Prepare l,alpha and others
    call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
    l_min = 2*pi/im
    allocate(l_lim(1:num_l))
    !
    do i=1,num_l
      l_lim(i) = exp(log(ratio_min)+(i-1) * (log(ratio_max)-log(ratio_min)) / (num_l-1)) * l_min
    enddo
    !
    if(mpirank==0)  print *, "Integrate point allocated"
    !
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !!!!
    allocate(Pis1(1:num_l),Pis2(1:num_l),Pim2(1:num_l),Pim3(1:num_l),Pid(1:num_l))
    !
    c_w1f = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1f, w1f,  [imfftw,jmfftw])
    c_w2f = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2f, w2f,  [imfftw,jmfftw])
    c_rho = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rho, rhof,[imfftw,jmfftw])
    !
    c_A11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A11, A11,[imfftw,jmfftw])
    c_A12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A12, A12,[imfftw,jmfftw])
    c_A21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A21, A21,[imfftw,jmfftw])
    c_A22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A22, A22,[imfftw,jmfftw])
    !
    Pis1=0.d0
    Pis2=0.d0
    Pim2=0.d0
    Pim3=0.d0
    Pid=0.d0
    !
    if(mpirank==0)  print *, "Array allocated and initialized"
    !
    do m=1,num_l
      !
      !!!!!! Filter to get Sij filted by l
      if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
      !
      !
      !!!! Velocity Favre average and density average
      ! After this bloc, w1 is (rho*u1) in spectral space
      do j=1,jm
      do i=1,im
        Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
        !
        w1f(i,j)    = w1(i,j)    *Gl
        w2f(i,j)    = w2(i,j)    *Gl
        !
        rhof(i,j)   = rhocom(i,j)*Gl
      enddo
      enddo
      !
      ! After this bloc, w1 is (rho*u1) in physical space
      call fftw_mpi_execute_dft(backward_plan,w1f,w1f)
      call fftw_mpi_execute_dft(backward_plan,w2f,w2f)
      call fftw_mpi_execute_dft(backward_plan,rhof,rhof)
      !
      ! After this bloc, w1 is u1 in physical space
      do j=1,jm
      do i=1,im
        w1f(i,j) = w1f(i,j)/rhof(i,j)
        w2f(i,j) = w2f(i,j)/rhof(i,j)
        !
        !
      enddo
      enddo
      !
      ! After this bloc, w1 is u1 in fourier space, A11 is A11 in fourier space
      call fftw_mpi_execute_dft(forward_plan,w1f,w1f)
      call fftw_mpi_execute_dft(forward_plan,w2f,w2f)
      !
      do j=1,jm
      do i=1,im
        !
        w1f(i,j)  = w1f(i,j)/(1.d0*ia*ja)
        w2f(i,j)  = w2f(i,j)/(1.d0*ia*ja)
        !
        A11(i,j) = imag*w1f(i,j)*k1(i,j)
        A21(i,j) = imag*w2f(i,j)*k1(i,j)
        A12(i,j) = imag*w1f(i,j)*k2(i,j)
        A22(i,j) = imag*w2f(i,j)*k2(i,j)
        !
      end do
      end do
      !
      !
      !
      ! After this bloc, A11 is A11 in physical space
      call fftw_mpi_execute_dft(backward_plan,A11,A11)
      call fftw_mpi_execute_dft(backward_plan,A21,A21)
      call fftw_mpi_execute_dft(backward_plan,A12,A12)
      call fftw_mpi_execute_dft(backward_plan,A22,A22)
      !
      do j=1,jm
      do i=1,im
        !
        All(i,j) = A11(i,j)+A22(i,j)
        !
        S11(i,j) = A11(i,j) - 1.d0/2.d0 * All(i,j)
        S22(i,j) = A22(i,j) - 1.d0/2.d0 * All(i,j)
        S12(i,j) = (A12(i,j) + A21(i,j))*0.5d0
        S21(i,j) = S12(i,j)
        !
        W12(i,j) = (A12(i,j)-A21(i,j))*0.5d0
        W21(i,j) = -1.d0*W12(i,j)
        !
      end do
      end do
      !
      do j=1,jm
      do i=1,im
        !
        Pis1(m) = Pis1(m) + rhof(i,j) *(S11(i,j)*S11(i,j)*S11(i,j)+ &
                                        S12(i,j)*S12(i,j)*S11(i,j)+ &
                                        S11(i,j)*S21(i,j)*S12(i,j)+ &
                                        S12(i,j)*S22(i,j)*S12(i,j)+ &
                                        S21(i,j)*S11(i,j)*S21(i,j)+ &
                                        S22(i,j)*S12(i,j)*S21(i,j)+ &
                                        S21(i,j)*S21(i,j)*S22(i,j)+ &
                                        S22(i,j)*S22(i,j)*S22(i,j)) * &
                                       l_lim(m) * l_lim(m)
        Pis2(m) = Pis2(m) + rhof(i,j) *(W12(i,j)*W21(i,j)*S11(i,j)+ &
                                          W21(i,j)*W12(i,j)*S22(i,j)) * &
                                        l_lim(m) * l_lim(m)
        Pim2(m) = Pim2(m) + rhof(i,j) *(S11(i,j)*S11(i,j)*All(i,j)+ &
                                          S12(i,j)*S12(i,j)*All(i,j)+ &
                                          S21(i,j)*S21(i,j)*All(i,j)+ &
                                          S22(i,j)*S22(i,j)*All(i,j)) * &
                                        l_lim(m) * l_lim(m)
        Pim3(m) = Pim3(m) + rhof(i,j) *(W12(i,j)*W21(i,j)*All(i,j)+ &
                                          W21(i,j)*W12(i,j)*All(i,j)) * &
                                      l_lim(m) * l_lim(m)
        Pid(m) = Pid(m) + rhof(i,j) * All(i,j) * All(i,j) * All(i,j) * &
                                    l_lim(m) * l_lim(m)
        !
      end do
      end do
      !
      if(mpirank==0)  print *, '** l filted!'
      !
      Pis1(m) =	 psum(Pis1(m)) / (ia*ja)
      Pis2(m) =	 - psum(Pis2(m)) / (ia*ja)
      Pim2(m) =	 psum(Pim2(m)) / (ia*ja) /2
      Pim3(m) =	 - psum(Pim3(m)) / (ia*ja) /2
      Pid(m) =	 psum(Pid(m)) / (ia*ja) /4
      !
      !
    enddo
    !
    if(mpirank==0)  print *, 'Job finish'
    !
    if(mpirank==0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/SGS_Pilocal_'//stepname//'.dat'
      else
        outfilename = 'pp/SGS_Pilocal.dat'
      endif
      
      call listinit(filename=outfilename,handle=hand_a, &
                    firstline='nstep time lOlmin pis1 pis2 pim2 pim3 pid')
      do m=1,num_l
        call listwrite(hand_a,l_lim(m)/l_min, Pis1(m), Pis2(m), Pim2(m), Pim3(m), Pid(m))
      enddo
      !
      print *, '>>>>', outfilename
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_w1)
    call fftw_free(c_w2)
    call fftw_free(c_rhocom)
    call fftw_free(c_w1f)
    call fftw_free(c_w2f)
    call fftw_free(c_rho)
    call fftw_free(c_A11)
    call fftw_free(c_A12)
    call fftw_free(c_A21)
    call fftw_free(c_A22)
    call mpistop
    deallocate(k1,k2)
    deallocate(l_lim)
    deallocate(All,S11,S12,S21,S22,W12,W21)
    deallocate(Pis1,Pis2,Pim2,Pim3,Pid)
    !
  end subroutine SGSPi2Dlocal
  !
  subroutine SGSPi2Dint(thefilenumb)
    !
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km,ia,ja,ka
    use commarray, only: vel, rho
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    include 'fftw3-mpi.f03'
    !
    integer,intent(in) :: thefilenumb
    integer :: fh
    integer :: i,j,m,n
    character(len=128) :: infilename,outfilename,outfilename2
    character(len=4) :: stepname,mname
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1,w2,rhocom
    real(8), allocatable, dimension(:,:) :: k1,k2
    complex(8) :: imag
    real(8),allocatable,dimension(:) :: l_lim
    real(8),allocatable,dimension(:,:) :: l_sqrtalpha,l_phi,dl_alpha
    integer,allocatable,dimension(:) :: num_alphas
    integer :: num_l,num_alpha,num_alphamin
    integer :: hand_a,hand_b
    real(8) :: l_min, ratio_max, ratio_min
    real(8) :: Gl,Galpha,Gphi
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1_filted,w2_filted,rho_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: A11_filted,A12_filted,A21_filted,A22_filted
    complex(8), allocatable, dimension(:,:) :: All_filted,S11_filted,S12_filted,S21_filted,S22_filted,W12_filted,W21_filted
    complex(8), allocatable, dimension(:,:) :: All_filted_l,S11_filted_l,S12_filted_l,S21_filted_l,S22_filted_l
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: term1_11,term1_12,term1_21,term1_22,term2, &
                                                          term3_11,term3_12,term3_21,term3_22, &
                                                          term4_11,term4_22,term5,term6_11,term6_12, &
                                                          term6_21,term6_22,term7
    real(8) :: vxr_D1,vxr_D2,vxr_D3,vxr_D4,vxr_D5,vxr_D6,vxr_D7
    type(C_PTR) :: c_w1,c_w2,c_rhocom,forward_plan,backward_plan
    type(C_PTR) :: c_w1_filted,c_w2_filted,c_rho_filted
    type(C_PTR) :: c_A11_filted,c_A12_filted,c_A21_filted,c_A22_filted
    type(C_PTR) :: c_term1_11,c_term1_12,c_term1_21,c_term1_22,c_term2,c_term3_11,c_term3_12,c_term3_21,c_term3_22, &
                                               c_term4_11,c_term4_22,c_term5,c_term6_11,c_term6_12,c_term6_21,c_term6_22,c_term7
    real(8), allocatable, dimension(:) :: Pi1,Pi2,Pi3,Pi4,Pi5,Pi6,Pi7
    integer,dimension(8) :: value
    character(len=1) :: modeio
    logical :: loutput
    !
    call readinput
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !!!! Read velocity and density field
    allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km))
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    !!!! Prepare initial field in Fourier space
    !! velocity
    c_w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1, w1, [imfftw,jmfftw])
    c_w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2, w2, [imfftw,jmfftw])
    c_rhocom = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw])
    !
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    do j=1,jm
    do i=1,im
      !
      w1(i,j)=CMPLX(vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
      w2(i,j)=CMPLX(vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
      rhocom(i,j)=CMPLX(rho(i,j,0),0.d0,C_INTPTR_T);
      !
    end do
    end do
    !
    !After this bloc, w1 is (rho*u1) in spectral space
    call fftw_mpi_execute_dft(forward_plan,w1,w1)
    call fftw_mpi_execute_dft(forward_plan,w2,w2)
    call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
    do j=1,jm
    do i=1,im
      !
      w1(i,j)=w1(i,j)/(1.d0*ia*ja)
      w2(i,j)=w2(i,j)/(1.d0*ia*ja)
      rhocom(i,j)=rhocom(i,j)/(1.d0*ia*ja)
      !
    end do
    end do
    !
    !! wavenumber
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    do j = 1,jm
    do i = 1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if((j+j0) <= (ja/2+1)) then
        k2(i,j) = real(j+j0-1,8)
      else if((j+j0)<=(ja)) then
        k2(i,j) = real(j+j0-ja-1,8)
      else
        print *,"Error, no wave number possible, (j+j0) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    !
    !! Imaginary number prepare
    imag = CMPLX(0.d0,1.d0,8)
    !
    if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
    !!!! Prepare l,alpha and others
    call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
    l_min = 2*pi/im
    allocate(l_lim(1:num_l),num_alphas(1:num_l),l_sqrtalpha(1:num_l,1:num_alpha))
    allocate(l_phi(1:num_l,1:num_alpha),dl_alpha(1:num_l,1:num_alpha))
    !
    do i=1,num_l
      l_lim(i) = exp(log(ratio_min)+(i-1) * (log(ratio_max)-log(ratio_min)) / (num_l-1)) * l_min
      num_alphas(i) = num_alphamin + int(sqrt(l_lim(i)/ratio_max/l_min)*real(num_alpha-num_alphamin))
      !
      do j=1,num_alphas(i)
        l_sqrtalpha(i,j) = sqrt( exp(log(0.2**2) + &
            (log((l_lim(i)/l_min)**2) - log(0.2**2))* (j-1) / (num_alphas(i)-1)) )*l_min
        l_phi(i,j) = sqrt( abs(l_lim(i)**2 - l_sqrtalpha(i,j)**2) )
      enddo
    enddo
    !
    !
    do i=1,num_l
      dl_alpha(i,1) = l_sqrtalpha(i,1)**2 

      do j=2,num_alphas(i)
        dl_alpha(i,j) = l_sqrtalpha(i,j)**2 -l_sqrtalpha(i,j-1)**2 
      enddo
    enddo
    !
    if(mpirank==0)  print *, "Integrate point allocated"
    !
    if(mpirank==0) then
      open(fh,file='pp/SGSintegral.info',form='formatted')
      write(fh,"(2(A9,1x))")'NumL','NumAlpha'
      write(fh,"(2(I9,1x))")num_l,num_alpha
      write(fh,"(2(A9,1x),2(A15,1x))")'i','j','l_lim','l_sqrtalpha'
      do i=1,num_l
        do j=1,num_alphas(i)
        ! Output file of rank information.
          write(fh,"(2(I9,1x),2(E15.7E3,1x))")i,j,l_lim(i),l_sqrtalpha(i,j)
        enddo
      enddo
      !
      close(fh)
      print*,' << SGSintegral.info ... done !'
    endif
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !!!! allocation
    allocate(Pi1(1:num_l), Pi2(1:num_l), Pi3(1:num_l), Pi4(1:num_l), Pi5(1:num_l), Pi6(1:num_l), Pi7(1:num_l))
    !
    c_w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1_filted, w1_filted, [imfftw,jmfftw])
    c_w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2_filted, w2_filted, [imfftw,jmfftw])
    c_rho_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rho_filted, rho_filted, [imfftw,jmfftw])
    !
    c_A11_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A11_filted, A11_filted, [imfftw,jmfftw])
    c_A12_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A12_filted, A12_filted, [imfftw,jmfftw])
    c_A21_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A21_filted, A21_filted, [imfftw,jmfftw])
    c_A22_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A22_filted, A22_filted, [imfftw,jmfftw])
    !
    allocate(All_filted(1:im,1:jm),S11_filted(1:im,1:jm),S12_filted(1:im,1:jm),S21_filted(1:im,1:jm),S22_filted(1:im,1:jm), &
             W12_filted(1:im,1:jm),W21_filted(1:im,1:jm))
    !
    allocate(All_filted_l(1:im,1:jm),S11_filted_l(1:im,1:jm),S12_filted_l(1:im,1:jm),S21_filted_l(1:im,1:jm),  &
            S22_filted_l(1:im,1:jm))
    !
    c_term1_11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term1_11, term1_11, [imfftw,jmfftw])
    c_term1_12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term1_12, term1_12, [imfftw,jmfftw])
    c_term1_21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term1_21, term1_21, [imfftw,jmfftw])
    c_term1_22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term1_22, term1_22, [imfftw,jmfftw])
    c_term2    = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term2, term2, [imfftw,jmfftw])
    c_term3_11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term3_11, term3_11, [imfftw,jmfftw])
    c_term3_12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term3_12, term3_12, [imfftw,jmfftw])
    c_term3_21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term3_21, term3_21, [imfftw,jmfftw])
    c_term3_22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term3_22, term3_22, [imfftw,jmfftw])
    c_term4_11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term4_11, term4_11, [imfftw,jmfftw])
    c_term4_22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term4_22, term4_22, [imfftw,jmfftw])
    c_term5    = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term5, term5, [imfftw,jmfftw])
    c_term6_11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term6_11, term6_11, [imfftw,jmfftw])
    c_term6_12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term6_12, term6_12, [imfftw,jmfftw])
    c_term6_21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term6_21, term6_21, [imfftw,jmfftw])
    c_term6_22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term6_22, term6_22, [imfftw,jmfftw])
    c_term7    = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term7, term7, [imfftw,jmfftw])
    !
    !
    Pi1 = 0.d0
    Pi2 =	0.d0
    Pi3 =	0.d0
    Pi4 =	0.d0
    Pi5 =	0.d0
    Pi6 =	0.d0
    Pi7 =	0.d0
    !
    if(mpirank==0)  print *, "Array allocated and initialized"
    !
    do m=1,num_l
      !
      !!!!!! Filter to get Sij filted by l
      if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
      !
      if(mpirank == 0) then
        write(mname,'(i4.4)')m
        if (thefilenumb .ne. 0) then
          outfilename2 = 'pp/SGS_Pi_precise_'//stepname//'_'//mname//'.dat'
        else
          outfilename2 = 'pp/SGS_Pi_precise_'//mname//'.dat'
        endif
        call listinit(filename=outfilename2,handle=hand_b, &
                    firstline='nstep time sqrtalpha pi1 pi2 pi3 pi4 pi5 pi6 pi7')
      endif
      !!!! Velocity Favre average and density average
      ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
      do j=1,jm
      do i=1,im
        Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) ! Filtre scale :1
        w1_filted(i,j) = w1(i,j)     *Gl
        w2_filted(i,j) = w2(i,j)     *Gl
        rho_filted(i,j) = rhocom(i,j)*Gl
      enddo
      enddo
      !
      ! After this bloc, w1_filted is (rho*u1)_filted in physical space
      call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
      call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
      !
      ! After this bloc, w1_filted is u1_filted in physical space
      do j=1,jm
      do i=1,im
        w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
        w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
      enddo
      enddo
      !
      ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
      call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
      call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
      do j=1,jm
      do i=1,im
        !
        w1_filted(i,j)  = w1_filted(i,j)/(1.d0*ia*ja)
        w2_filted(i,j)  = w2_filted(i,j)/(1.d0*ia*ja)
        A11_filted(i,j) = imag*w1_filted(i,j)*k1(i,j)
        A21_filted(i,j) = imag*w2_filted(i,j)*k1(i,j)
        A12_filted(i,j) = imag*w1_filted(i,j)*k2(i,j)
        A22_filted(i,j) = imag*w2_filted(i,j)*k2(i,j)
        !
      end do
      end do
      !
      !
      ! After this bloc, A11_filted is A11_filted in physical space
      call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
      call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
      call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
      call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
      !
      do j=1,jm
      do i=1,im
        All_filted_l(i,j) = A11_filted(i,j)+A22_filted(i,j)
      !
        S11_filted_l(i,j) = A11_filted(i,j) - 1.d0/2.d0 * All_filted_l(i,j)
        S12_filted_l(i,j) = (A12_filted(i,j) + A21_filted(i,j))*0.5d0
        S21_filted_l(i,j) = S12_filted_l(i,j)
        S22_filted_l(i,j) = A22_filted(i,j) - 1.d0/2.d0 * All_filted_l(i,j)
      !
      end do
      end do
      !
      if(mpirank==0)  print *, '** l filted!'
      !
      !!!!!! Begin integral
      !
      do n=1,num_alphas(m)
        !
        call date_and_time(values=value) 
        !
        if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alphas(m),',now is ',&
                                value(5), ':', value(6),':',value(7)
        !!!! Velocity Favre average and density average
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        do i=1,im
        do j=1,jm
          Galpha = exp(-(k1(i,j)**2+k2(i,j)**2)*l_sqrtalpha(m,n)**2/2.d0) ! Filtre scale :sqrtalpha
          w1_filted(i,j)  = w1(i,j)    *Galpha
          w2_filted(i,j)  = w2(i,j)    *Galpha
          rho_filted(i,j) = rhocom(i,j)*Galpha
        enddo
        enddo
        !
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        do i=1,im
        do j=1,jm
          w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
          w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
        enddo
        enddo
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j)  = w1_filted(i,j)/(1.d0*ia*ja)
          w2_filted(i,j)  = w2_filted(i,j)/(1.d0*ia*ja)
          A11_filted(i,j) = imag*w1_filted(i,j)*k1(i,j)
          A21_filted(i,j) = imag*w2_filted(i,j)*k1(i,j)
          A12_filted(i,j) = imag*w1_filted(i,j)*k2(i,j)
          A22_filted(i,j) = imag*w2_filted(i,j)*k2(i,j)
          !
        end do
        end do
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
        call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
        call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
        call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
        !
        !
        do j=1,jm
        do i=1,im
          All_filted(i,j) = A11_filted(i,j)+A22_filted(i,j)
          ! 
          S11_filted(i,j) = A11_filted(i,j) - 1.d0/2.d0 * All_filted(i,j)
          S12_filted(i,j) = (A12_filted(i,j) + A21_filted(i,j))*0.5d0
          S21_filted(i,j) = S12_filted(i,j)
          S22_filted(i,j) = A22_filted(i,j) - 1.d0/2.d0 * All_filted(i,j)
          !
          W12_filted(i,j) = (A12_filted(i,j)-A21_filted(i,j))*0.5d0
          W21_filted(i,j) = -1.d0*W12_filted(i,j)
        enddo
        enddo
        !
        !!!! Pi terms
        !
        do j=1,jm
        do i=1,im
          !term1_IJ = rho_filted*SI1_filted*SJ1_filted + rho_filted*SI2_filted*SJ2_filted
          term1_11(i,j) = rho_filted(i,j)*S11_filted(i,j)*S11_filted(i,j) + rho_filted(i,j)*S12_filted(i,j)*S12_filted(i,j)
          term1_21(i,j) = rho_filted(i,j)*S21_filted(i,j)*S11_filted(i,j) + rho_filted(i,j)*S22_filted(i,j)*S12_filted(i,j)
          term1_12(i,j) = rho_filted(i,j)*S11_filted(i,j)*S21_filted(i,j) + rho_filted(i,j)*S12_filted(i,j)*S22_filted(i,j)
          term1_22(i,j) = rho_filted(i,j)*S21_filted(i,j)*S21_filted(i,j) + rho_filted(i,j)*S22_filted(i,j)*S22_filted(i,j)
          ! 
          !
          term2(i,j) = rho_filted(i,j)*S11_filted(i,j)*S11_filted(i,j)+rho_filted(i,j)*S12_filted(i,j)*S12_filted(i,j)+ &
                  rho_filted(i,j)*S21_filted(i,j)*S21_filted(i,j)+rho_filted(i,j)*S22_filted(i,j)*S22_filted(i,j)
          !
          ! term3_IJ = rho_filted*All_filted*SIJ_filted
          term3_11(i,j) = rho_filted(i,j)*All_filted(i,j)*S11_filted(i,j)
          term3_12(i,j) = rho_filted(i,j)*All_filted(i,j)*S12_filted(i,j)
          term3_21(i,j) = rho_filted(i,j)*All_filted(i,j)*S21_filted(i,j)
          term3_22(i,j) = rho_filted(i,j)*All_filted(i,j)*S22_filted(i,j)
          !
          !term4_IJ = rho_filted*WI1_filted*W1J_filted+rho_filted*WI2_filted*W2J_filted
          term4_11(i,j) = rho_filted(i,j)*W12_filted(i,j)*W21_filted(i,j)
          term4_22(i,j) = rho_filted(i,j)*W21_filted(i,j)*W12_filted(i,j)
          !
          term5(i,j) = rho_filted(i,j)*W12_filted(i,j)*W21_filted(i,j) + rho_filted(i,j)*W21_filted(i,j)*W12_filted(i,j)
          !
          !term6_IJ= rho_filted*(S1J_filted*WI1_filted-SI1_filted*W1J_filted) + rho_filted*(S2J_filted*WI2_filted-SI2_filted*W2J_filted)
          term6_11(i,j)=   rho_filted(i,j)*(S21_filted(i,j)*W12_filted(i,j)-S12_filted(i,j)*W21_filted(i,j))
          term6_12(i,j)= - rho_filted(i,j)*(S11_filted(i,j)*W12_filted(i,j)) + rho_filted(i,j)*(S22_filted(i,j)*W12_filted(i,j))
          term6_21(i,j)=   rho_filted(i,j)*(S11_filted(i,j)*W21_filted(i,j)) - rho_filted(i,j)*(S22_filted(i,j)*W21_filted(i,j))
          term6_22(i,j)=   rho_filted(i,j)*(S12_filted(i,j)*W21_filted(i,j)-S21_filted(i,j)*W12_filted(i,j))
          !
          term7(i,j) = rho_filted(i,j)*All_filted(i,j)*All_filted(i,j)
        enddo
        enddo
        !
        ! Do filter phi:
        ! F -> product -> F inverse
        call fftw_mpi_execute_dft(forward_plan,term1_11,term1_11)
        call fftw_mpi_execute_dft(forward_plan,term1_12,term1_12)
        call fftw_mpi_execute_dft(forward_plan,term1_21,term1_21)
        call fftw_mpi_execute_dft(forward_plan,term1_22,term1_22)
        call fftw_mpi_execute_dft(forward_plan,term2   ,term2   )
        call fftw_mpi_execute_dft(forward_plan,term3_11,term3_11)
        call fftw_mpi_execute_dft(forward_plan,term3_12,term3_12)
        call fftw_mpi_execute_dft(forward_plan,term3_21,term3_21)
        call fftw_mpi_execute_dft(forward_plan,term3_22,term3_22)
        call fftw_mpi_execute_dft(forward_plan,term4_11,term4_11)
        call fftw_mpi_execute_dft(forward_plan,term4_22,term4_22)
        call fftw_mpi_execute_dft(forward_plan,term5   ,term5   )
        call fftw_mpi_execute_dft(forward_plan,term6_11,term6_11)
        call fftw_mpi_execute_dft(forward_plan,term6_12,term6_12)
        call fftw_mpi_execute_dft(forward_plan,term6_21,term6_21)
        call fftw_mpi_execute_dft(forward_plan,term6_22,term6_22)
        call fftw_mpi_execute_dft(forward_plan,term7   ,term7   )
        !
        do j=1,jm
        do i=1,im
          Gphi = exp(-(k1(i,j)**2+k2(i,j)**2)*l_phi(m,n)**2/2.d0) ! Filtre scale :phi
          term1_11(i,j) = term1_11(i,j)*Gphi/(1.d0*ia*ja)
          term1_12(i,j) = term1_12(i,j)*Gphi/(1.d0*ia*ja)
          term1_21(i,j) = term1_21(i,j)*Gphi/(1.d0*ia*ja)
          term1_22(i,j) = term1_22(i,j)*Gphi/(1.d0*ia*ja)
          term2(i,j)    = term2(i,j)   *Gphi/(1.d0*ia*ja)
          term3_11(i,j) = term3_11(i,j)*Gphi/(1.d0*ia*ja)
          term3_12(i,j) = term3_12(i,j)*Gphi/(1.d0*ia*ja)
          term3_21(i,j) = term3_21(i,j)*Gphi/(1.d0*ia*ja)
          term3_22(i,j) = term3_22(i,j)*Gphi/(1.d0*ia*ja)
          term4_11(i,j) = term4_11(i,j)*Gphi/(1.d0*ia*ja)
          term4_22(i,j) = term4_22(i,j)*Gphi/(1.d0*ia*ja)
          term5(i,j)    = term5(i,j)   *Gphi/(1.d0*ia*ja)
          term6_11(i,j) = term6_11(i,j)*Gphi/(1.d0*ia*ja)
          term6_12(i,j) = term6_12(i,j)*Gphi/(1.d0*ia*ja)
          term6_21(i,j) = term6_21(i,j)*Gphi/(1.d0*ia*ja)
          term6_22(i,j) = term6_22(i,j)*Gphi/(1.d0*ia*ja)
          term7(i,j)    = term7(i,j)   *Gphi/(1.d0*ia*ja)
          !
        enddo
        enddo
        !
        !
        call fftw_mpi_execute_dft(backward_plan,term1_11,term1_11)
        call fftw_mpi_execute_dft(backward_plan,term1_12,term1_12)
        call fftw_mpi_execute_dft(backward_plan,term1_21,term1_21)
        call fftw_mpi_execute_dft(backward_plan,term1_22,term1_22)
        call fftw_mpi_execute_dft(backward_plan,term2   ,term2   )
        call fftw_mpi_execute_dft(backward_plan,term3_11,term3_11)
        call fftw_mpi_execute_dft(backward_plan,term3_12,term3_12)
        call fftw_mpi_execute_dft(backward_plan,term3_21,term3_21)
        call fftw_mpi_execute_dft(backward_plan,term3_22,term3_22)
        call fftw_mpi_execute_dft(backward_plan,term4_11,term4_11)
        call fftw_mpi_execute_dft(backward_plan,term4_22,term4_22)
        call fftw_mpi_execute_dft(backward_plan,term5   ,term5   )
        call fftw_mpi_execute_dft(backward_plan,term6_11,term6_11)
        call fftw_mpi_execute_dft(backward_plan,term6_12,term6_12)
        call fftw_mpi_execute_dft(backward_plan,term6_21,term6_21)
        call fftw_mpi_execute_dft(backward_plan,term6_22,term6_22)
        call fftw_mpi_execute_dft(backward_plan,term7   ,term7   )
        !
        !
        do j=1,jm
        do i=1,im
          vxr_D1 = term1_11(i,j) * S11_filted_l(i,j) + term1_12(i,j) * S12_filted_l(i,j) + &
                  term1_21(i,j) * S21_filted_l(i,j) + term1_22(i,j) * S22_filted_l(i,j)
          Pi1(m) = Pi1(m) + vxr_D1 * dl_alpha(m,n)
          !
          vxr_D2 = term2(i,j) * All_filted_l(i,j)
          Pi2(m) = Pi2(m) + vxr_D2 * dl_alpha(m,n) * 0.5d0 
          !
          vxr_D3 = term3_11(i,j) * S11_filted_l(i,j) + term3_12(i,j) * S12_filted_l(i,j) + &
                  term3_21(i,j) * S21_filted_l(i,j) + term3_22(i,j) * S22_filted_l(i,j)
          Pi3(m) = Pi3(m) + vxr_D3 * dl_alpha(m,n)
          !
          vxr_D4 = term4_11(i,j) * S11_filted_l(i,j) + term4_22(i,j) * S22_filted_l(i,j)
          Pi4(m) = Pi4(m) - vxr_D4 * dl_alpha(m,n)
          !
          vxr_D5 = term5(i,j) * All_filted_l(i,j)
          Pi5(m) = Pi5(m) - vxr_D5 * dl_alpha(m,n) * 0.5d0
          !
          vxr_D6 = term6_11(i,j) * S11_filted_l(i,j) + term6_12(i,j) * S12_filted_l(i,j) + &
                  term6_21(i,j) * S21_filted_l(i,j) + term6_22(i,j) * S22_filted_l(i,j)
          Pi6(m) = Pi6(m) + vxr_D6 * dl_alpha(m,n)
          !
          vxr_D7 = term7(i,j) * All_filted_l(i,j)
          Pi7(m) = Pi7(m) + vxr_D7 * dl_alpha(m,n) * 0.25d0
        enddo
        enddo
        !
        vxr_D1 = psum(vxr_D1)
        vxr_D2 = psum(vxr_D2)
        vxr_D3 = psum(vxr_D3)
        vxr_D4 = psum(vxr_D4)
        vxr_D5 = psum(vxr_D5)
        vxr_D6 = psum(vxr_D6)
        vxr_D7 = psum(vxr_D7)
        !
        if(mpirank==0) then
          call listwrite(hand_b,l_sqrtalpha(m,n),vxr_D1 * dl_alpha(m,n), &
                        vxr_D2 * dl_alpha(m,n) / 2.d0, &
                        vxr_D3 * dl_alpha(m,n) , &
                        - vxr_D4 * dl_alpha(m,n), &
                        - vxr_D5 * dl_alpha(m,n) * 0.5d0, &
                        vxr_D6 * dl_alpha(m,n), &
                        vxr_D7 * dl_alpha(m,n) * 0.25d0)
        endif
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
      enddo
      !
      Pi1(m) =	 psum(Pi1(m)) / (ia*ja)
      Pi2(m) =	 psum(Pi2(m)) / (ia*ja)
      Pi3(m) =	 psum(Pi3(m)) / (ia*ja)
      Pi4(m) =	 psum(Pi4(m)) / (ia*ja)
      Pi5(m) =	 psum(Pi5(m)) / (ia*ja)
      Pi6(m) =	 psum(Pi6(m)) / (ia*ja)
      Pi7(m) =	 psum(Pi7(m)) / (ia*ja)
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
    enddo
    if(mpirank==0)  print *, 'Job finish'
    !
    if(mpirank==0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/SGS_Pi_'//stepname//'.dat'
      else
        outfilename = 'pp/SGS_Pi.dat'
      endif
      
      call listinit(filename=outfilename,handle=hand_a, &
                    firstline='nstep time lOlmin pi1 pi2 pi3 pi4 pi5 pi6 pi7')
      do m=1,num_l
        call listwrite(hand_a,l_lim(m)/l_min, Pi1(m), Pi2(m), &
                        Pi3(m), Pi4(m), Pi5(m),     &
                        Pi6(m), Pi7(m))
      end do
      !
      print *, '>>>>', outfilename
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_w1)
    call fftw_free(c_w2)
    call fftw_free(c_rhocom)
    call fftw_free(c_w1_filted)
    call fftw_free(c_w2_filted)
    call fftw_free(c_rho_filted)
    call fftw_free(c_A11_filted)
    call fftw_free(c_A12_filted)
    call fftw_free(c_A21_filted)
    call fftw_free(c_A22_filted)
    call fftw_free(c_term1_11)
    call fftw_free(c_term1_12)
    call fftw_free(c_term1_21)
    call fftw_free(c_term1_22)
    call fftw_free(c_term2)
    call fftw_free(c_term3_11)
    call fftw_free(c_term3_12)
    call fftw_free(c_term3_21)
    call fftw_free(c_term3_22)
    call fftw_free(c_term4_11)
    call fftw_free(c_term4_22)
    call fftw_free(c_term5)
    call fftw_free(c_term6_11)
    call fftw_free(c_term6_12)
    call fftw_free(c_term6_21)
    call fftw_free(c_term6_22)
    call fftw_free(c_term7)
    call mpistop
    deallocate(All_filted,S11_filted,S12_filted,S21_filted,S22_filted,W12_filted,W21_filted)
    deallocate(All_filted_l,S11_filted_l,S12_filted_l,S21_filted_l,S22_filted_l)
    deallocate(k1,k2)
    deallocate(l_lim,l_sqrtalpha,l_phi,dl_alpha)
    deallocate(Pi1,Pi2,Pi3,Pi4,Pi5,Pi6,Pi7)
    !
  end subroutine SGSPi2Dint
  !
  subroutine SGSPi3Dtot(thefilenumb)
    ! 
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km
    use commarray, only: vel, rho
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
    include 'fftw3-mpi.f03'
    !
    integer,intent(in) :: thefilenumb
    integer :: fh
    integer :: i,j,k,m,n
    character(len=128) :: infilename,outfilename,outfilename2
    character(len=4) :: stepname,mname
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
    complex(8) :: imag
    real(8),allocatable,dimension(:) :: l_lim
    integer :: num_l,num_alpha,num_alphamin
    integer :: hand_a,hand_b
    real(8) :: l_min, ratio_max, ratio_min
    real(8) :: Gl
    real(8), allocatable, dimension(:) :: Pi_tot
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1w1,w1w2,w1w3
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w2w1,w2w2,w2w3
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w3w1,w3w2,w3w3
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1w1_filted,w1w2_filted,w1w3_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w2w1_filted,w2w2_filted,w2w3_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w3w1_filted,w3w2_filted,w3w3_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11_filted,A12_filted,A13_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21_filted,A22_filted,A23_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31_filted,A32_filted,A33_filted
    real(8), allocatable, dimension(:,:,:,:,:) :: tau ! 1:3,1:3,1:im,1:jm,1:km
    !
    !
    type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,forward_plan,backward_plan
    type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted
    type(C_PTR) :: c_w1w1,c_w1w2,c_w1w3
    type(C_PTR) :: c_w2w1,c_w2w2,c_w2w3
    type(C_PTR) :: c_w3w1,c_w3w2,c_w3w3
    type(C_PTR) :: c_w1w1_filted,c_w1w2_filted,c_w1w3_filted
    type(C_PTR) :: c_w2w1_filted,c_w2w2_filted,c_w2w3_filted
    type(C_PTR) :: c_w3w1_filted,c_w3w2_filted,c_w3w3_filted
    type(C_PTR) :: c_A11_filted,c_A12_filted,c_A13_filted
    type(C_PTR) :: c_A21_filted,c_A22_filted,c_A23_filted
    type(C_PTR) :: c_A31_filted,c_A32_filted,c_A33_filted
    !
    integer,dimension(8) :: value
    character(len=1) :: modeio
    logical :: loutput
    !
    call readinput
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !!!! Read velocity and density field
    allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    !!!! Prepare initial field in Fourier space
    !! velocity
    c_w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
    c_w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
    c_w3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
    c_rhocom = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
    !
    c_w1w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w1, w1w1, [imfftw,jmfftw,kmfftw])
    c_w1w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w2, w1w2, [imfftw,jmfftw,kmfftw])
    c_w1w3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w3, w1w3, [imfftw,jmfftw,kmfftw])
    c_w2w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w1, w2w1, [imfftw,jmfftw,kmfftw])
    c_w2w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w2, w2w2, [imfftw,jmfftw,kmfftw])
    c_w2w3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w3, w2w3, [imfftw,jmfftw,kmfftw])
    c_w3w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3w1, w3w1, [imfftw,jmfftw,kmfftw])
    c_w3w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3w2, w3w2, [imfftw,jmfftw,kmfftw])
    c_w3w3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3w3, w3w3, [imfftw,jmfftw,kmfftw])
    !
    c_w1w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w1_filted, w1w1_filted, [imfftw,jmfftw,kmfftw])
    c_w1w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w2_filted, w1w2_filted, [imfftw,jmfftw,kmfftw])
    c_w1w3_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w3_filted, w1w3_filted, [imfftw,jmfftw,kmfftw])
    c_w2w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w1_filted, w2w1_filted, [imfftw,jmfftw,kmfftw])
    c_w2w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w2_filted, w2w2_filted, [imfftw,jmfftw,kmfftw])
    c_w2w3_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w3_filted, w2w3_filted, [imfftw,jmfftw,kmfftw])
    c_w3w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3w1_filted, w3w1_filted, [imfftw,jmfftw,kmfftw])
    c_w3w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3w2_filted, w3w2_filted, [imfftw,jmfftw,kmfftw])
    c_w3w3_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3w3_filted, w3w3_filted, [imfftw,jmfftw,kmfftw])
    !
    forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
      w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
      w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
      rhocom(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T);
      w1w1(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
      w1w2(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
      w1w3(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
      w2w1(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
      w2w2(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
      w2w3(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
      w3w1(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
      w3w2(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
      w3w3(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
      !
    end do
    end do
    end do
    !
    !After this bloc, w1 is (rho*u1) in spectral space
    call fftw_mpi_execute_dft(forward_plan,w1,w1)
    call fftw_mpi_execute_dft(forward_plan,w2,w2)
    call fftw_mpi_execute_dft(forward_plan,w3,w3)
    call fftw_mpi_execute_dft(forward_plan,w1w1,w1w1)
    call fftw_mpi_execute_dft(forward_plan,w1w2,w1w2)
    call fftw_mpi_execute_dft(forward_plan,w1w3,w1w3)
    call fftw_mpi_execute_dft(forward_plan,w2w1,w2w1)
    call fftw_mpi_execute_dft(forward_plan,w2w2,w2w2)
    call fftw_mpi_execute_dft(forward_plan,w2w3,w2w3)
    call fftw_mpi_execute_dft(forward_plan,w3w1,w3w1)
    call fftw_mpi_execute_dft(forward_plan,w3w2,w3w2)
    call fftw_mpi_execute_dft(forward_plan,w3w3,w3w3)
    call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
      w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
      w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
      !
      w1w1(i,j,k)=w1w1(i,j,k)/(1.d0*ia*ja*ka)
      w1w2(i,j,k)=w1w2(i,j,k)/(1.d0*ia*ja*ka)
      w1w3(i,j,k)=w1w3(i,j,k)/(1.d0*ia*ja*ka)
      w2w1(i,j,k)=w2w1(i,j,k)/(1.d0*ia*ja*ka)
      w2w2(i,j,k)=w2w2(i,j,k)/(1.d0*ia*ja*ka)
      w2w3(i,j,k)=w2w3(i,j,k)/(1.d0*ia*ja*ka)
      w3w1(i,j,k)=w3w1(i,j,k)/(1.d0*ia*ja*ka)
      w3w2(i,j,k)=w3w2(i,j,k)/(1.d0*ia*ja*ka)
      w3w3(i,j,k)=w3w3(i,j,k)/(1.d0*ia*ja*ka)
      !
      rhocom(i,j,k)=rhocom(i,j,k)/(1.d0*ia*ja*ka)
      !
    end do
    end do
    end do

    !
    !
    !! wavenumber
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    do k = 1,km
    do j = 1,jm
    do i = 1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j,k) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j,k) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if(j <= (ja/2+1)) then
        k2(i,j,k) = real(j-1,8)
      else if(i<=(ia)) then
        k2(i,j,k) = real(j-ja-1,8)
      else
        print *,"Error, no wave number possible, j must smaller than ja-1 !"
      end if
      !
      if((k+k0) <= (ka/2+1)) then
        k3(i,j,k) = real(k+k0-1,8)
      else if((k+k0)<=(ka)) then
        k3(i,j,k) = real(k+k0-ka-1,8)
      else
        print *,"Error, no wave number possible, (k+k0) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    end do
    !
    !! Imaginary number prepare
    imag = CMPLX(0.d0,1.d0,8)
    !
    allocate(tau(1:3,1:3,1:im,1:jm,1:km))
    !
    if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
    !!!! Prepare l,alpha and others
    call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
    l_min = 2*pi/im
    allocate(l_lim(1:num_l))
    !
    do i=1,num_l
      l_lim(i) = exp(log(ratio_min)+(i-1) * (log(ratio_max)-log(ratio_min)) / (num_l-1)) * l_min
    enddo
    !
    if(mpirank==0)  print *, "Integrate point allocated"
    !
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !!!!
    allocate(Pi_tot(1:num_l))
    !
    c_w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw,kmfftw])
    c_w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw,kmfftw])
    c_w3_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3_filted, w3_filted,  [imfftw,jmfftw,kmfftw])
    c_rho_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw,kmfftw])
    !
    c_A11_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A11_filted, A11_filted,[imfftw,jmfftw,kmfftw])
    c_A12_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A12_filted, A12_filted,[imfftw,jmfftw,kmfftw])
    c_A13_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A13_filted, A13_filted,[imfftw,jmfftw,kmfftw])
    c_A21_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A21_filted, A21_filted,[imfftw,jmfftw,kmfftw])
    c_A22_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A22_filted, A22_filted,[imfftw,jmfftw,kmfftw])
    c_A23_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A23_filted, A23_filted,[imfftw,jmfftw,kmfftw])
    c_A31_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A31_filted, A31_filted,[imfftw,jmfftw,kmfftw])
    c_A32_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A32_filted, A32_filted,[imfftw,jmfftw,kmfftw])
    c_A33_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A33_filted, A33_filted,[imfftw,jmfftw,kmfftw])
    !
    Pi_tot = 0.d0
    !
    if(mpirank==0)  print *, "Array allocated and initialized"
    !
    do m=1,num_l
      !
      !!!!!! Filter to get Sij filted by l
      if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
      !
      !
      !!!! Velocity Favre average and density average
      ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
      do k=1,km
      do j=1,jm
      do i=1,im
        Gl = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
        !
        w1_filted(i,j,k)    = w1(i,j,k)    *Gl
        w2_filted(i,j,k)    = w2(i,j,k)    *Gl
        w3_filted(i,j,k)    = w3(i,j,k)    *Gl
        !
        w1w1_filted(i,j,k)  = w1w1(i,j,k) *Gl
        w1w2_filted(i,j,k)  = w1w2(i,j,k) *Gl
        w1w3_filted(i,j,k)  = w1w3(i,j,k) *Gl
        w2w1_filted(i,j,k)  = w2w1(i,j,k) *Gl
        w2w2_filted(i,j,k)  = w2w2(i,j,k) *Gl
        w2w3_filted(i,j,k)  = w2w3(i,j,k) *Gl
        w3w1_filted(i,j,k)  = w3w1(i,j,k) *Gl
        w3w2_filted(i,j,k)  = w3w2(i,j,k) *Gl
        w3w3_filted(i,j,k)  = w3w3(i,j,k) *Gl
        !
        rho_filted(i,j,k)   = rhocom(i,j,k)*Gl
      enddo
      enddo
      enddo
      !
      ! After this bloc, w1_filted is (rho*u1)_filted in physical space
      call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
      call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
      call fftw_mpi_execute_dft(backward_plan,w1w1_filted,w1w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w1w2_filted,w1w2_filted)
      call fftw_mpi_execute_dft(backward_plan,w1w3_filted,w1w3_filted)
      call fftw_mpi_execute_dft(backward_plan,w2w1_filted,w2w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w2w2_filted,w2w2_filted)
      call fftw_mpi_execute_dft(backward_plan,w2w3_filted,w2w3_filted)
      call fftw_mpi_execute_dft(backward_plan,w3w1_filted,w3w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w3w2_filted,w3w2_filted)
      call fftw_mpi_execute_dft(backward_plan,w3w3_filted,w3w3_filted)
      call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
      !
      ! After this bloc, w1_filted is u1_filted in physical space
      do k=1,km
      do j=1,jm
      do i=1,im
        w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
        w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
        w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
        !
        tau(1,1,i,j,k) = dreal(w1w1_filted(i,j,k) - rho_filted(i,j,k) * w1_filted(i,j,k) * w1_filted(i,j,k))
        tau(1,2,i,j,k) = dreal(w1w2_filted(i,j,k) - rho_filted(i,j,k) * w1_filted(i,j,k) * w2_filted(i,j,k))
        tau(1,3,i,j,k) = dreal(w1w3_filted(i,j,k) - rho_filted(i,j,k) * w1_filted(i,j,k) * w3_filted(i,j,k))
        tau(2,1,i,j,k) = dreal(w2w1_filted(i,j,k) - rho_filted(i,j,k) * w2_filted(i,j,k) * w1_filted(i,j,k))
        tau(2,2,i,j,k) = dreal(w2w2_filted(i,j,k) - rho_filted(i,j,k) * w2_filted(i,j,k) * w2_filted(i,j,k))
        tau(2,3,i,j,k) = dreal(w2w3_filted(i,j,k) - rho_filted(i,j,k) * w2_filted(i,j,k) * w3_filted(i,j,k))
        tau(3,1,i,j,k) = dreal(w3w1_filted(i,j,k) - rho_filted(i,j,k) * w3_filted(i,j,k) * w1_filted(i,j,k))
        tau(3,2,i,j,k) = dreal(w3w2_filted(i,j,k) - rho_filted(i,j,k) * w3_filted(i,j,k) * w2_filted(i,j,k))
        tau(3,3,i,j,k) = dreal(w3w3_filted(i,j,k) - rho_filted(i,j,k) * w3_filted(i,j,k) * w3_filted(i,j,k))
        !
      enddo
      enddo
      enddo
      !
      ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
      call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
      call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
      call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1_filted(i,j,k)  = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
        w2_filted(i,j,k)  = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
        w3_filted(i,j,k)  = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
        !
        A11_filted(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
        A21_filted(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
        A31_filted(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
        A12_filted(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
        A22_filted(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
        A32_filted(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
        A13_filted(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
        A23_filted(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
        A33_filted(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
        !
      end do
      end do
      end do
      !
      !
      !
      ! After this bloc, A11_filted is A11_filted in physical space
      call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
      call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
      call fftw_mpi_execute_dft(backward_plan,A31_filted,A31_filted)
      call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
      call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
      call fftw_mpi_execute_dft(backward_plan,A32_filted,A32_filted)
      call fftw_mpi_execute_dft(backward_plan,A13_filted,A13_filted)
      call fftw_mpi_execute_dft(backward_plan,A23_filted,A23_filted)
      call fftw_mpi_execute_dft(backward_plan,A33_filted,A33_filted)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        Pi_tot(m) = Pi_tot(m) + tau(1,1,i,j,k) * A11_filted(i,j,k) + &
                                tau(1,2,i,j,k) * A12_filted(i,j,k) + &
                                tau(1,3,i,j,k) * A13_filted(i,j,k) + &
                                tau(2,1,i,j,k) * A21_filted(i,j,k) + &
                                tau(2,2,i,j,k) * A22_filted(i,j,k) + &
                                tau(2,3,i,j,k) * A23_filted(i,j,k) + &
                                tau(3,1,i,j,k) * A31_filted(i,j,k) + &
                                tau(3,2,i,j,k) * A32_filted(i,j,k) + &
                                tau(3,3,i,j,k) * A33_filted(i,j,k)
        !
      end do
      end do
      end do
      !
      if(mpirank==0)  print *, '** l filted!'
      !
      Pi_tot(m) =	 psum(Pi_tot(m)) / (ia*ja*ka)
      !
      !
    enddo
    if(mpirank==0)  print *, 'Job finish'
    !
    if(mpirank==0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/SGS_Pitot_'//stepname//'.dat'
      else
        outfilename = 'pp/SGS_Pitot.dat'
      endif
      
      call listinit(filename=outfilename,handle=hand_a, &
                    firstline='nstep time lOlmin pitot')
      do m=1,num_l
        call listwrite(hand_a,l_lim(m)/l_min, Pi_tot(m))
      enddo
      !
      print *, '>>>>', outfilename
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_w1)
    call fftw_free(c_w2)
    call fftw_free(c_w3)
    call fftw_free(c_rhocom)
    call fftw_free(c_w1_filted)
    call fftw_free(c_w2_filted)
    call fftw_free(c_w3_filted)
    call fftw_free(c_w1w1)
    call fftw_free(c_w1w2)
    call fftw_free(c_w1w3)
    call fftw_free(c_w2w1)
    call fftw_free(c_w2w2)
    call fftw_free(c_w2w3)
    call fftw_free(c_w3w1)
    call fftw_free(c_w3w2)
    call fftw_free(c_w3w3)
    call fftw_free(c_w1w1_filted)
    call fftw_free(c_w1w2_filted)
    call fftw_free(c_w1w3_filted)
    call fftw_free(c_w2w1_filted)
    call fftw_free(c_w2w2_filted)
    call fftw_free(c_w2w3_filted)
    call fftw_free(c_w3w1_filted)
    call fftw_free(c_w3w2_filted)
    call fftw_free(c_w3w3_filted)
    call fftw_free(c_rho_filted)
    call fftw_free(c_A11_filted)
    call fftw_free(c_A12_filted)
    call fftw_free(c_A13_filted)
    call fftw_free(c_A21_filted)
    call fftw_free(c_A22_filted)
    call fftw_free(c_A23_filted)
    call fftw_free(c_A31_filted)
    call fftw_free(c_A32_filted)
    call fftw_free(c_A33_filted)
    call mpistop
    deallocate(k1,k2,k3,tau)
    deallocate(l_lim)
    deallocate(Pi_tot)
    !
  end subroutine SGSPi3Dtot
  !
  subroutine SGSPi3Dlocal(thefilenumb)
    ! 
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km
    use commarray, only: vel, rho
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
    include 'fftw3-mpi.f03'
    !
    integer,intent(in) :: thefilenumb
    integer :: fh
    integer :: i,j,k,m,n
    character(len=128) :: infilename,outfilename,outfilename2
    character(len=4) :: stepname,mname
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
    complex(8) :: imag
    real(8),allocatable,dimension(:) :: l_lim
    integer :: num_l,num_alpha,num_alphamin
    integer :: hand_a,hand_b
    real(8) :: l_min, ratio_max, ratio_min
    real(8) :: Gl
    real(8), allocatable, dimension(:) :: Pis1,Pis2,Pim2,Pim3,Pid
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1f,w2f,w3f,rhof
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11,A12,A13
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21,A22,A23
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31,A32,A33
    !
    complex(8), allocatable, dimension(:,:,:) :: All
    complex(8), allocatable, dimension(:,:,:) :: S11,S12,S13
    complex(8), allocatable, dimension(:,:,:) :: S21,S22,S23
    complex(8), allocatable, dimension(:,:,:) :: S31,S32,S33
    complex(8), allocatable, dimension(:,:,:) :: W12,W21
    complex(8), allocatable, dimension(:,:,:) :: W13,W31
    complex(8), allocatable, dimension(:,:,:) :: W23,W32
    !
    type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,forward_plan,backward_plan
    type(C_PTR) :: c_w1f,c_w2f,c_w3f,c_rho
    type(C_PTR) :: c_A11,c_A12,c_A13
    type(C_PTR) :: c_A21,c_A22,c_A23
    type(C_PTR) :: c_A31,c_A32,c_A33
    !
    integer,dimension(8) :: value
    character(len=1) :: modeio
    logical :: loutput
    !
    call readinput
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !!!! Read velocity and density field
    allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    !!!! Prepare initial field in Fourier space
    !! velocity
    c_w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
    c_w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
    c_w3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
    c_rhocom = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
    !
    forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    allocate(All(1:im,1:jm,1:km),&
            S11(1:im,1:jm,1:km),S12(1:im,1:jm,1:km),S13(1:im,1:jm,1:km),&
            S21(1:im,1:jm,1:km),S22(1:im,1:jm,1:km),S23(1:im,1:jm,1:km),&
            S31(1:im,1:jm,1:km),S32(1:im,1:jm,1:km),S33(1:im,1:jm,1:km),&
            W12(1:im,1:jm,1:km),W21(1:im,1:jm,1:km),&
            W13(1:im,1:jm,1:km),W31(1:im,1:jm,1:km),&
            W23(1:im,1:jm,1:km),W32(1:im,1:jm,1:km))
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
      w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
      w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
      rhocom(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T);
      !
    end do
    end do
    end do
    !
    !After this bloc, w1 is (rho*u1) in spectral space
    call fftw_mpi_execute_dft(forward_plan,w1,w1)
    call fftw_mpi_execute_dft(forward_plan,w2,w2)
    call fftw_mpi_execute_dft(forward_plan,w3,w3)
    call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
      w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
      w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
      !
      rhocom(i,j,k)=rhocom(i,j,k)/(1.d0*ia*ja*ka)
      !
    end do
    end do
    end do

    !
    !
    !! wavenumber
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    do k = 1,km
    do j = 1,jm
    do i = 1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j,k) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j,k) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if(j <= (ja/2+1)) then
        k2(i,j,k) = real(j-1,8)
      else if(i<=(ia)) then
        k2(i,j,k) = real(j-ja-1,8)
      else
        print *,"Error, no wave number possible, j must smaller than ja-1 !"
      end if
      !
      if((k+k0) <= (ka/2+1)) then
        k3(i,j,k) = real(k+k0-1,8)
      else if((k+k0)<=(ka)) then
        k3(i,j,k) = real(k+k0-ka-1,8)
      else
        print *,"Error, no wave number possible, (k+k0) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    end do
    !
    !! Imaginary number prepare
    imag = CMPLX(0.d0,1.d0,8)
    !
    !
    if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
    !!!! Prepare l,alpha and others
    call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
    l_min = 2*pi/im
    allocate(l_lim(1:num_l))
    !
    do i=1,num_l
      l_lim(i) = exp(log(ratio_min)+(i-1) * (log(ratio_max)-log(ratio_min)) / (num_l-1)) * l_min
    enddo
    !
    if(mpirank==0)  print *, "Integrate point allocated"
    !
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !!!!
    allocate(Pis1(1:num_l),Pis2(1:num_l),Pim2(1:num_l),Pim3(1:num_l),Pid(1:num_l))
    !
    c_w1f = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1f, w1f,  [imfftw,jmfftw,kmfftw])
    c_w2f = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2f, w2f,  [imfftw,jmfftw,kmfftw])
    c_w3f = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3f, w3f,  [imfftw,jmfftw,kmfftw])
    c_rho = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rho, rhof,[imfftw,jmfftw,kmfftw])
    !
    c_A11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A11, A11,[imfftw,jmfftw,kmfftw])
    c_A12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A12, A12,[imfftw,jmfftw,kmfftw])
    c_A13 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A13, A13,[imfftw,jmfftw,kmfftw])
    c_A21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A21, A21,[imfftw,jmfftw,kmfftw])
    c_A22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A22, A22,[imfftw,jmfftw,kmfftw])
    c_A23 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A23, A23,[imfftw,jmfftw,kmfftw])
    c_A31 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A31, A31,[imfftw,jmfftw,kmfftw])
    c_A32 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A32, A32,[imfftw,jmfftw,kmfftw])
    c_A33 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A33, A33,[imfftw,jmfftw,kmfftw])
    !
    Pis1=0.d0
    Pis2=0.d0
    Pim2=0.d0
    Pim3=0.d0
    Pid=0.d0
    !
    if(mpirank==0)  print *, "Array allocated and initialized"
    !
    do m=1,num_l
      !
      !!!!!! Filter to get Sij filted by l
      if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
      !
      !
      !!!! Velocity Favre average and density average
      ! After this bloc, w1 is (rho*u1) in spectral space
      do k=1,km
      do j=1,jm
      do i=1,im
        Gl = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
        !
        w1f(i,j,k)    = w1(i,j,k)    *Gl
        w2f(i,j,k)    = w2(i,j,k)    *Gl
        w3f(i,j,k)    = w3(i,j,k)    *Gl
        !
        rhof(i,j,k)   = rhocom(i,j,k)*Gl
      enddo
      enddo
      enddo
      !
      ! After this bloc, w1 is (rho*u1) in physical space
      call fftw_mpi_execute_dft(backward_plan,w1f,w1f)
      call fftw_mpi_execute_dft(backward_plan,w2f,w2f)
      call fftw_mpi_execute_dft(backward_plan,w3f,w3f)
      call fftw_mpi_execute_dft(backward_plan,rhof,rhof)
      !
      ! After this bloc, w1 is u1 in physical space
      do k=1,km
      do j=1,jm
      do i=1,im
        w1f(i,j,k) = w1f(i,j,k)/rhof(i,j,k)
        w2f(i,j,k) = w2f(i,j,k)/rhof(i,j,k)
        w3f(i,j,k) = w3f(i,j,k)/rhof(i,j,k)
        !
        !
      enddo
      enddo
      enddo
      !
      ! After this bloc, w1 is u1 in fourier space, A11 is A11 in fourier space
      call fftw_mpi_execute_dft(forward_plan,w1f,w1f)
      call fftw_mpi_execute_dft(forward_plan,w2f,w2f)
      call fftw_mpi_execute_dft(forward_plan,w3f,w3f)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1f(i,j,k)  = w1f(i,j,k)/(1.d0*ia*ja*ka)
        w2f(i,j,k)  = w2f(i,j,k)/(1.d0*ia*ja*ka)
        w3f(i,j,k)  = w3f(i,j,k)/(1.d0*ia*ja*ka)
        !
        A11(i,j,k) = imag*w1f(i,j,k)*k1(i,j,k)
        A21(i,j,k) = imag*w2f(i,j,k)*k1(i,j,k)
        A31(i,j,k) = imag*w3f(i,j,k)*k1(i,j,k)
        A12(i,j,k) = imag*w1f(i,j,k)*k2(i,j,k)
        A22(i,j,k) = imag*w2f(i,j,k)*k2(i,j,k)
        A32(i,j,k) = imag*w3f(i,j,k)*k2(i,j,k)
        A13(i,j,k) = imag*w1f(i,j,k)*k3(i,j,k)
        A23(i,j,k) = imag*w2f(i,j,k)*k3(i,j,k)
        A33(i,j,k) = imag*w3f(i,j,k)*k3(i,j,k)
        !
      end do
      end do
      end do
      !
      !
      !
      ! After this bloc, A11 is A11 in physical space
      call fftw_mpi_execute_dft(backward_plan,A11,A11)
      call fftw_mpi_execute_dft(backward_plan,A21,A21)
      call fftw_mpi_execute_dft(backward_plan,A31,A31)
      call fftw_mpi_execute_dft(backward_plan,A12,A12)
      call fftw_mpi_execute_dft(backward_plan,A22,A22)
      call fftw_mpi_execute_dft(backward_plan,A32,A32)
      call fftw_mpi_execute_dft(backward_plan,A13,A13)
      call fftw_mpi_execute_dft(backward_plan,A23,A23)
      call fftw_mpi_execute_dft(backward_plan,A33,A33)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        All(i,j,k) = A11(i,j,k)+A22(i,j,k)+A33(i,j,k)
        !
        S11(i,j,k) = A11(i,j,k) - 1.d0/3.d0 * All(i,j,k)
        S22(i,j,k) = A22(i,j,k) - 1.d0/3.d0 * All(i,j,k)
        S33(i,j,k) = A33(i,j,k) - 1.d0/3.d0 * All(i,j,k)
        S12(i,j,k) = (A12(i,j,k) + A21(i,j,k))*0.5d0
        S21(i,j,k) = S12(i,j,k)
        S13(i,j,k) = (A13(i,j,k) + A31(i,j,k))*0.5d0
        S31(i,j,k) = S13(i,j,k)
        S23(i,j,k) = (A23(i,j,k) + A32(i,j,k))*0.5d0
        S32(i,j,k) = S23(i,j,k)
        !
        W12(i,j,k) = (A12(i,j,k)-A21(i,j,k))*0.5d0
        W21(i,j,k) = -1.d0*W12(i,j,k)
        W13(i,j,k) = (A13(i,j,k)-A31(i,j,k))*0.5d0
        W31(i,j,k) = -1.d0*W13(i,j,k)
        W23(i,j,k) = (A23(i,j,k)-A32(i,j,k))*0.5d0
        W32(i,j,k) = -1.d0*W23(i,j,k)
        !
      end do
      end do
      end do
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        Pis1(m) = Pis1(m) + rhof(i,j,k) *(S11(i,j,k)*S11(i,j,k)*S11(i,j,k)+ &
                                        S12(i,j,k)*S12(i,j,k)*S11(i,j,k)+ &
                                        S13(i,j,k)*S13(i,j,k)*S11(i,j,k)+ &
                                        S11(i,j,k)*S21(i,j,k)*S12(i,j,k)+ &
                                        S12(i,j,k)*S22(i,j,k)*S12(i,j,k)+ &
                                        S13(i,j,k)*S23(i,j,k)*S12(i,j,k)+ &
                                        S11(i,j,k)*S31(i,j,k)*S13(i,j,k)+ &
                                        S12(i,j,k)*S32(i,j,k)*S13(i,j,k)+ &
                                        S13(i,j,k)*S33(i,j,k)*S13(i,j,k)+ &
                                        S21(i,j,k)*S11(i,j,k)*S21(i,j,k)+ &
                                        S22(i,j,k)*S12(i,j,k)*S21(i,j,k)+ &
                                        S23(i,j,k)*S13(i,j,k)*S21(i,j,k)+ &
                                        S21(i,j,k)*S21(i,j,k)*S22(i,j,k)+ &
                                        S22(i,j,k)*S22(i,j,k)*S22(i,j,k)+ &
                                        S23(i,j,k)*S23(i,j,k)*S22(i,j,k)+ &
                                        S21(i,j,k)*S31(i,j,k)*S23(i,j,k)+ &
                                        S22(i,j,k)*S32(i,j,k)*S23(i,j,k)+ &
                                        S23(i,j,k)*S33(i,j,k)*S23(i,j,k)+ &
                                        S31(i,j,k)*S11(i,j,k)*S31(i,j,k)+ &
                                        S32(i,j,k)*S12(i,j,k)*S31(i,j,k)+ &
                                        S33(i,j,k)*S13(i,j,k)*S31(i,j,k)+ &
                                        S31(i,j,k)*S21(i,j,k)*S32(i,j,k)+ &
                                        S32(i,j,k)*S22(i,j,k)*S32(i,j,k)+ &
                                        S33(i,j,k)*S23(i,j,k)*S32(i,j,k)+ &
                                        S31(i,j,k)*S31(i,j,k)*S33(i,j,k)+ &
                                        S32(i,j,k)*S32(i,j,k)*S33(i,j,k)+ &
                                        S33(i,j,k)*S33(i,j,k)*S33(i,j,k)) * &
                                       l_lim(m) * l_lim(m)
        Pis2(m) = Pis2(m) + rhof(i,j,k) *(W12(i,j,k)*W21(i,j,k)*S11(i,j,k)+ &
                                          W13(i,j,k)*W31(i,j,k)*S11(i,j,k)+ &
                                          W13(i,j,k)*W32(i,j,k)*S12(i,j,k)+ &
                                          W12(i,j,k)*W23(i,j,k)*S13(i,j,k)+ &
                                          W23(i,j,k)*W31(i,j,k)*S21(i,j,k)+ &
                                          W21(i,j,k)*W12(i,j,k)*S22(i,j,k)+ &
                                          W23(i,j,k)*W32(i,j,k)*S22(i,j,k)+ &
                                          W21(i,j,k)*W13(i,j,k)*S23(i,j,k)+ &
                                          W32(i,j,k)*W21(i,j,k)*S31(i,j,k)+ &
                                          W31(i,j,k)*W12(i,j,k)*S32(i,j,k)+ &
                                          W31(i,j,k)*W13(i,j,k)*S33(i,j,k)+ &
                                          W32(i,j,k)*W23(i,j,k)*S33(i,j,k)) * &
                                        l_lim(m) * l_lim(m)
        Pim2(m) = Pim2(m) + rhof(i,j,k) *(S11(i,j,k)*S11(i,j,k)*All(i,j,k)+ &
                                          S12(i,j,k)*S12(i,j,k)*All(i,j,k)+ &
                                          S13(i,j,k)*S13(i,j,k)*All(i,j,k)+ &
                                          S21(i,j,k)*S21(i,j,k)*All(i,j,k)+ &
                                          S22(i,j,k)*S22(i,j,k)*All(i,j,k)+ &
                                          S23(i,j,k)*S23(i,j,k)*All(i,j,k)+ &
                                          S31(i,j,k)*S31(i,j,k)*All(i,j,k)+ &
                                          S32(i,j,k)*S32(i,j,k)*All(i,j,k)+ &
                                          S33(i,j,k)*S33(i,j,k)*All(i,j,k)) * &
                                        l_lim(m) * l_lim(m)
        Pim3(m) = Pim3(m) + rhof(i,j,k) *(W12(i,j,k)*W21(i,j,k)*All(i,j,k)+ &
                                          W13(i,j,k)*W31(i,j,k)*All(i,j,k)+ &
                                          W21(i,j,k)*W12(i,j,k)*All(i,j,k)+ &
                                          W23(i,j,k)*W32(i,j,k)*All(i,j,k)+ &
                                          W31(i,j,k)*W13(i,j,k)*All(i,j,k)+ &
                                          W32(i,j,k)*W23(i,j,k)*All(i,j,k)) * &
                                      l_lim(m) * l_lim(m)
        Pid(m) = Pid(m) + rhof(i,j,k) * All(i,j,k) * All(i,j,k) * All(i,j,k) * &
                                    l_lim(m) * l_lim(m)
        !
      end do
      end do
      end do
      !
      if(mpirank==0)  print *, '** l filted!'
      !
      Pis1(m) =	 psum(Pis1(m)) / (ia*ja*ka)
      Pis2(m) =	 - psum(Pis2(m)) / (ia*ja*ka)
      Pim2(m) =	 psum(Pim2(m)) / (ia*ja*ka) /3
      Pim3(m) =	 - psum(Pim3(m)) / (ia*ja*ka) /3
      Pid(m) =	 psum(Pid(m)) / (ia*ja*ka) /9
      !
      !
    enddo
    if(mpirank==0)  print *, 'Job finish'
    !
    if(mpirank==0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/SGS_Pilocal_'//stepname//'.dat'
      else
        outfilename = 'pp/SGS_Pilocal.dat'
      endif
      
      call listinit(filename=outfilename,handle=hand_a, &
                    firstline='nstep time lOlmin pis1 pis2 pim2 pim3 pid')
      do m=1,num_l
        call listwrite(hand_a,l_lim(m)/l_min, Pis1(m), Pis2(m), Pim2(m), Pim3(m), Pid(m))
      enddo
      !
      print *, '>>>>', outfilename
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_w1)
    call fftw_free(c_w2)
    call fftw_free(c_w3)
    call fftw_free(c_rhocom)
    call fftw_free(c_w1f)
    call fftw_free(c_w2f)
    call fftw_free(c_w3f)
    call fftw_free(c_rho)
    call fftw_free(c_A11)
    call fftw_free(c_A12)
    call fftw_free(c_A13)
    call fftw_free(c_A21)
    call fftw_free(c_A22)
    call fftw_free(c_A23)
    call fftw_free(c_A31)
    call fftw_free(c_A32)
    call fftw_free(c_A33)
    call mpistop
    deallocate(k1,k2,k3)
    deallocate(l_lim)
    deallocate(All,S11,S12,S13,S21,S22,S23,S31,S32,S33)
    deallocate(W12,W13,W21,W23,W31,W32)
    deallocate(Pis1,Pis2,Pim2,Pim3,Pid)
    !
  end subroutine SGSPi3Dlocal
  !
  subroutine SGSLES3D(thefilenumb)
    ! 
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km
    use commarray, only: vel, rho
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
    include 'fftw3-mpi.f03'
    !
    integer,intent(in) :: thefilenumb
    integer :: fh
    integer :: i,j,k,m,n
    character(len=128) :: infilename,outfilename,outfilename2
    character(len=4) :: stepname,mname
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
    complex(8) :: imag
    real(8),allocatable,dimension(:) :: l_lim
    integer :: num_l,num_alpha,num_alphamin
    integer :: hand_a,hand_b
    real(8) :: l_min, ratio_max, ratio_min
    real(8) :: Gl
    real(8), allocatable, dimension(:) :: AllAll,SijSij
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11_filted,A12_filted,A13_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21_filted,A22_filted,A23_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31_filted,A32_filted,A33_filted
    complex(8), allocatable, dimension(:,:,:) :: All_filted
    complex(8), allocatable, dimension(:,:,:) :: S11_filted,S12_filted,S13_filted
    complex(8), allocatable, dimension(:,:,:) :: S21_filted,S22_filted,S23_filted
    complex(8), allocatable, dimension(:,:,:) :: S31_filted,S32_filted,S33_filted
    !
    !
    type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,forward_plan,backward_plan
    type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted
    type(C_PTR) :: c_A11_filted,c_A12_filted,c_A13_filted
    type(C_PTR) :: c_A21_filted,c_A22_filted,c_A23_filted
    type(C_PTR) :: c_A31_filted,c_A32_filted,c_A33_filted
    !
    integer,dimension(8) :: value
    character(len=1) :: modeio
    logical :: loutput
    !
    call readinput
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !!!! Read velocity and density field
    allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    !!!! Prepare initial field in Fourier space
    !! velocity
    c_w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
    c_w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
    c_w3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
    c_rhocom = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
    !
    forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
      w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
      w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
      rhocom(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T);
      !
    end do
    end do
    end do
    !
    !After this bloc, w1 is (rho*u1) in spectral space
    call fftw_mpi_execute_dft(forward_plan,w1,w1)
    call fftw_mpi_execute_dft(forward_plan,w2,w2)
    call fftw_mpi_execute_dft(forward_plan,w3,w3)
    call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
      w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
      w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
      !
      rhocom(i,j,k)=rhocom(i,j,k)/(1.d0*ia*ja*ka)
      !
    end do
    end do
    end do

    !
    !
    !! wavenumber
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    do k = 1,km
    do j = 1,jm
    do i = 1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j,k) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j,k) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if(j <= (ja/2+1)) then
        k2(i,j,k) = real(j-1,8)
      else if(i<=(ia)) then
        k2(i,j,k) = real(j-ja-1,8)
      else
        print *,"Error, no wave number possible, j must smaller than ja-1 !"
      end if
      !
      if((k+k0) <= (ka/2+1)) then
        k3(i,j,k) = real(k+k0-1,8)
      else if((k+k0)<=(ka)) then
        k3(i,j,k) = real(k+k0-ka-1,8)
      else
        print *,"Error, no wave number possible, (k+k0) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    end do
    !
    !! Imaginary number prepare
    imag = CMPLX(0.d0,1.d0,8)
    !
    if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
    !!!! Prepare l,alpha and others
    call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
    l_min = 2*pi/im
    allocate(l_lim(1:num_l))
    !
    do i=1,num_l
      l_lim(i) = exp(log(ratio_min)+(i-1) * (log(ratio_max)-log(ratio_min)) / (num_l-1)) * l_min
    enddo
    !
    if(mpirank==0)  print *, "Integrate point allocated"
    !
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !!!!
    !
    c_w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw,kmfftw])
    c_w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw,kmfftw])
    c_w3_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3_filted, w3_filted,  [imfftw,jmfftw,kmfftw])
    c_rho_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw,kmfftw])
    !
    c_A11_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A11_filted, A11_filted,[imfftw,jmfftw,kmfftw])
    c_A12_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A12_filted, A12_filted,[imfftw,jmfftw,kmfftw])
    c_A13_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A13_filted, A13_filted,[imfftw,jmfftw,kmfftw])
    c_A21_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A21_filted, A21_filted,[imfftw,jmfftw,kmfftw])
    c_A22_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A22_filted, A22_filted,[imfftw,jmfftw,kmfftw])
    c_A23_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A23_filted, A23_filted,[imfftw,jmfftw,kmfftw])
    c_A31_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A31_filted, A31_filted,[imfftw,jmfftw,kmfftw])
    c_A32_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A32_filted, A32_filted,[imfftw,jmfftw,kmfftw])
    c_A33_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A33_filted, A33_filted,[imfftw,jmfftw,kmfftw])
    !
    !
    allocate(All_filted(1:im,1:jm,1:km),&
            S11_filted(1:im,1:jm,1:km),S12_filted(1:im,1:jm,1:km),S13_filted(1:im,1:jm,1:km),&
            S21_filted(1:im,1:jm,1:km),S22_filted(1:im,1:jm,1:km),S23_filted(1:im,1:jm,1:km),&
            S31_filted(1:im,1:jm,1:km),S32_filted(1:im,1:jm,1:km),S33_filted(1:im,1:jm,1:km))
    !
    allocate(AllAll(1:num_l),SijSij(1:num_l))
    !
    if(mpirank==0)  print *, "Array allocated and initialized"
    !
    AllAll = 0.d0
    SijSij = 0.d0
    !
    do m=1,num_l
      !
      !!!!!! Filter to get Sij filted by l
      if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
      !
      !
      !!!! Velocity Favre average and density average
      ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
      do k=1,km
      do j=1,jm
      do i=1,im
        Gl = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
        !
        w1_filted(i,j,k)    = w1(i,j,k)    *Gl
        w2_filted(i,j,k)    = w2(i,j,k)    *Gl
        w3_filted(i,j,k)    = w3(i,j,k)    *Gl
        rho_filted(i,j,k)   = rhocom(i,j,k)*Gl
      enddo
      enddo
      enddo
      !
      ! After this bloc, w1_filted is (rho*u1)_filted in physical space
      call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
      call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
      call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
      !
      ! After this bloc, w1_filted is u1_filted in physical space
      do k=1,km
      do j=1,jm
      do i=1,im
        w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
        w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
        w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
        !
        !
      enddo
      enddo
      enddo
      !
      ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
      call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
      call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
      call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1_filted(i,j,k)  = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
        w2_filted(i,j,k)  = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
        w3_filted(i,j,k)  = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
        !
        A11_filted(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
        A21_filted(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
        A31_filted(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
        A12_filted(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
        A22_filted(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
        A32_filted(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
        A13_filted(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
        A23_filted(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
        A33_filted(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
        !
      end do
      end do
      end do
      !
      !
      !
      ! After this bloc, A11_filted is A11_filted in physical space
      call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
      call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
      call fftw_mpi_execute_dft(backward_plan,A31_filted,A31_filted)
      call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
      call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
      call fftw_mpi_execute_dft(backward_plan,A32_filted,A32_filted)
      call fftw_mpi_execute_dft(backward_plan,A13_filted,A13_filted)
      call fftw_mpi_execute_dft(backward_plan,A23_filted,A23_filted)
      call fftw_mpi_execute_dft(backward_plan,A33_filted,A33_filted)
      !
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        All_filted(i,j,k) = A11_filted(i,j,k)+A22_filted(i,j,k)+A33_filted(i,j,k)
        !
        S11_filted(i,j,k) = A11_filted(i,j,k) - 1.d0/3.d0 * All_filted(i,j,k)
        S22_filted(i,j,k) = A22_filted(i,j,k) - 1.d0/3.d0 * All_filted(i,j,k)
        S33_filted(i,j,k) = A33_filted(i,j,k) - 1.d0/3.d0 * All_filted(i,j,k)
        S12_filted(i,j,k) = (A12_filted(i,j,k) + A21_filted(i,j,k))*0.5d0
        S21_filted(i,j,k) = S12_filted(i,j,k)
        S13_filted(i,j,k) = (A13_filted(i,j,k) + A31_filted(i,j,k))*0.5d0
        S31_filted(i,j,k) = S13_filted(i,j,k)
        S23_filted(i,j,k) = (A23_filted(i,j,k) + A32_filted(i,j,k))*0.5d0
        S32_filted(i,j,k) = S23_filted(i,j,k)
        !
        AllAll(m) = AllAll(m) + dreal(All_filted(i,j,k))*dreal(All_filted(i,j,k))
        SijSij(m) = SijSij(m) + dreal(S11_filted(i,j,k))* dreal(S11_filted(i,j,k)) + &
                                dreal(S12_filted(i,j,k))* dreal(S12_filted(i,j,k)) + &
                                dreal(S13_filted(i,j,k))* dreal(S13_filted(i,j,k)) + &
                                dreal(S21_filted(i,j,k))* dreal(S21_filted(i,j,k)) + &
                                dreal(S22_filted(i,j,k))* dreal(S22_filted(i,j,k)) + &
                                dreal(S23_filted(i,j,k))* dreal(S23_filted(i,j,k)) + &
                                dreal(S31_filted(i,j,k))* dreal(S31_filted(i,j,k)) + &
                                dreal(S32_filted(i,j,k))* dreal(S32_filted(i,j,k)) + &
                                dreal(S33_filted(i,j,k))* dreal(S33_filted(i,j,k))
        !
      end do
      end do
      end do
      !
      if(mpirank==0)  print *, '** l filted!'
      !
      AllAll(m) =	 psum(AllAll(m)) / (ia*ja*ka)
      SijSij(m) =	 psum(SijSij(m)) / (ia*ja*ka)
      !
      !
    enddo
    if(mpirank==0)  print *, 'Job finish'
    !
    if(mpirank==0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/SGS_LES_'//stepname//'.dat'
      else
        outfilename = 'pp/SGS_LES.dat'
      endif
      
      call listinit(filename=outfilename,handle=hand_a, &
                    firstline='nstep time lOlmin AllAll SijSij')
      do m=1,num_l
        call listwrite(hand_a,l_lim(m)/l_min, AllAll(m),SijSij(m))
      enddo
      !
      print *, '>>>>', outfilename
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_w1)
    call fftw_free(c_w2)
    call fftw_free(c_w3)
    call fftw_free(c_rhocom)
    call fftw_free(c_w1_filted)
    call fftw_free(c_w2_filted)
    call fftw_free(c_w3_filted)
    call fftw_free(c_rho_filted)
    call fftw_free(c_A11_filted)
    call fftw_free(c_A12_filted)
    call fftw_free(c_A13_filted)
    call fftw_free(c_A21_filted)
    call fftw_free(c_A22_filted)
    call fftw_free(c_A23_filted)
    call fftw_free(c_A31_filted)
    call fftw_free(c_A32_filted)
    call fftw_free(c_A33_filted)
    call mpistop
    deallocate(All_filted,S11_filted,S12_filted,S13_filted,&
                S21_filted,S22_filted,S23_filted,&
                S31_filted,S32_filted,S33_filted)
    deallocate(k1,k2,k3)
    deallocate(l_lim)
    deallocate(AllAll,SijSij)
    !
  end subroutine SGSLES3D
  !
  subroutine SGSPi3Dint(thefilenumb)
    ! 
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km
    use commarray, only: vel, rho
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
    include 'fftw3-mpi.f03'
    !
    integer,intent(in) :: thefilenumb
    integer :: fh
    integer :: i,j,k,m,n
    character(len=128) :: infilename,outfilename,outfilename2
    character(len=4) :: stepname,mname
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
    complex(8) :: imag
    real(8),allocatable,dimension(:) :: l_lim
    real(8),allocatable,dimension(:,:) :: l_sqrtalpha,l_phi,dl_alpha
    integer,allocatable,dimension(:) :: num_alphas
    integer :: num_l,num_alpha,num_alphamin
    integer :: hand_a,hand_b
    real(8) :: l_min, ratio_max, ratio_min
    real(8) :: Gl,Galpha,Gphi
    real(8), allocatable, dimension(:) :: Pi1,Pi2,Pi3,Pi4,Pi5,Pi6,Pi7
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11_filted,A12_filted,A13_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21_filted,A22_filted,A23_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31_filted,A32_filted,A33_filted
    complex(8), allocatable, dimension(:,:,:) :: All_filted_l
    complex(8), allocatable, dimension(:,:,:) :: S11_filted_l,S12_filted_l,S13_filted_l
    complex(8), allocatable, dimension(:,:,:) :: S21_filted_l,S22_filted_l,S23_filted_l
    complex(8), allocatable, dimension(:,:,:) :: S31_filted_l,S32_filted_l,S33_filted_l
    complex(8), allocatable, dimension(:,:,:) :: All_filted
    complex(8), allocatable, dimension(:,:,:) :: S11_filted,S12_filted,S13_filted
    complex(8), allocatable, dimension(:,:,:) :: S21_filted,S22_filted,S23_filted
    complex(8), allocatable, dimension(:,:,:) :: S31_filted,S32_filted,S33_filted
    complex(8), allocatable, dimension(:,:,:) :: W12_filted,W21_filted
    complex(8), allocatable, dimension(:,:,:) :: W13_filted,W31_filted
    complex(8), allocatable, dimension(:,:,:) :: W23_filted,W32_filted
    !
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: term1_11,term1_12,term1_13,&
                                                            term1_21,term1_22,term1_23,&
                                                            term1_31,term1_32,term1_33,&
                                                            term2,term5,term7,&
                                                            term3_11,term3_12,term3_13,&
                                                            term3_21,term3_22,term3_23,&
                                                            term3_31,term3_32,term3_33,&
                                                            term4_11,term4_12,term4_13,&
                                                            term4_21,term4_22,term4_23,&
                                                            term4_31,term4_32,term4_33,&
                                                            term6_11,term6_12,term6_13,&
                                                            term6_21,term6_22,term6_23,&
                                                            term6_31,term6_32,term6_33
    real(8) :: vxr_D1,vxr_D2,vxr_D3,vxr_D4,vxr_D5,vxr_D6,vxr_D7
    !
    type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,forward_plan,backward_plan
    type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted
    type(C_PTR) :: c_A11_filted,c_A12_filted,c_A13_filted
    type(C_PTR) :: c_A21_filted,c_A22_filted,c_A23_filted
    type(C_PTR) :: c_A31_filted,c_A32_filted,c_A33_filted
    type(C_PTR) :: c_term1_11,c_term1_12,c_term1_13,c_term1_21,c_term1_22,c_term1_23,c_term1_31,c_term1_32,c_term1_33
    type(C_PTR) :: c_term2,c_term5,c_term7
    type(C_PTR) :: c_term3_11,c_term3_12,c_term3_13,c_term3_21,c_term3_22,c_term3_23,c_term3_31,c_term3_32,c_term3_33
    type(C_PTR) :: c_term4_11,c_term4_12,c_term4_13,c_term4_21,c_term4_22,c_term4_23,c_term4_31,c_term4_32,c_term4_33
    type(C_PTR) :: c_term6_11,c_term6_12,c_term6_13,c_term6_21,c_term6_22,c_term6_23,c_term6_31,c_term6_32,c_term6_33
    !
    integer,dimension(8) :: value
    character(len=1) :: modeio
    logical :: loutput
    !
    call readinput
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !!!! Read velocity and density field
    allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    !!!! Prepare initial field in Fourier space
    !! velocity
    c_w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
    c_w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
    c_w3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
    c_rhocom = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
    !
    forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
      w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
      w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
      rhocom(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T);
      !
    end do
    end do
    end do
    !
    !After this bloc, w1 is (rho*u1) in spectral space
    call fftw_mpi_execute_dft(forward_plan,w1,w1)
    call fftw_mpi_execute_dft(forward_plan,w2,w2)
    call fftw_mpi_execute_dft(forward_plan,w3,w3)
    call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
      w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
      w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
      !
      rhocom(i,j,k)=rhocom(i,j,k)/(1.d0*ia*ja*ka)
      !
    end do
    end do
    end do

    !
    !
    !! wavenumber
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    do k = 1,km
    do j = 1,jm
    do i = 1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j,k) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j,k) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if(j <= (ja/2+1)) then
        k2(i,j,k) = real(j-1,8)
      else if(i<=(ia)) then
        k2(i,j,k) = real(j-ja-1,8)
      else
        print *,"Error, no wave number possible, j must smaller than ja-1 !"
      end if
      !
      if((k+k0) <= (ka/2+1)) then
        k3(i,j,k) = real(k+k0-1,8)
      else if((k+k0)<=(ka)) then
        k3(i,j,k) = real(k+k0-ka-1,8)
      else
        print *,"Error, no wave number possible, (k+k0) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    end do
    !
    !! Imaginary number prepare
    imag = CMPLX(0.d0,1.d0,8)
    !
    if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
    !!!! Prepare l,alpha and others
    call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
    l_min = 2*pi/im
    allocate(l_lim(1:num_l),num_alphas(1:num_l),l_sqrtalpha(1:num_l,1:num_alpha))
    allocate(l_phi(1:num_l,1:num_alpha),dl_alpha(1:num_l,1:num_alpha))
    !
    do i=1,num_l
      l_lim(i) = exp(log(ratio_min)+(i-1) * (log(ratio_max)-log(ratio_min)) / (num_l-1)) * l_min
      num_alphas(i) = num_alphamin + int(sqrt(l_lim(i)/ratio_max/l_min)*real(num_alpha-num_alphamin))
      !
      do j=1,num_alphas(i)
        l_sqrtalpha(i,j) = sqrt( exp(log(0.2**2) + &
            (log((l_lim(i)/l_min)**2) - log(0.2**2))* (j-1) / (num_alphas(i)-1)) )*l_min
        l_phi(i,j) = sqrt( abs(l_lim(i)**2 - l_sqrtalpha(i,j)**2) )
      enddo
    enddo
    !
    do i=1,num_l
      dl_alpha(i,1) = l_sqrtalpha(i,1)**2 

      do j=2,num_alphas(i)
        dl_alpha(i,j) = l_sqrtalpha(i,j)**2 -l_sqrtalpha(i,j-1)**2 
      enddo
    enddo
    !
    if(mpirank==0)  print *, "Integrate point allocated"
    !
    if(mpirank==0) then
      open(fh,file='pp/SGSintegral.info',form='formatted')
      write(fh,"(2(A9,1x))")'NumL','NumAlpha'
      write(fh,"(2(I9,1x))")num_l,num_alpha
      write(fh,"(2(A9,1x),2(A15,1x))")'i','j','l_lim','l_sqrtalpha'
      do i=1,num_l
        do j=1,num_alphas(i)
        ! Output file of rank information.
          write(fh,"(2(I9,1x),2(E15.7E3,1x))")i,j,l_lim(i),l_sqrtalpha(i,j)
        enddo
      enddo
      !
      close(fh)
      print*,' << SGSintegral.info ... done !'
    endif
    !
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !!!!
    allocate(Pi1(1:num_l), Pi2(1:num_l), Pi3(1:num_l), Pi4(1:num_l), Pi5(1:num_l), Pi6(1:num_l), Pi7(1:num_l))
    !
    c_w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw,kmfftw])
    c_w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw,kmfftw])
    c_w3_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3_filted, w3_filted,  [imfftw,jmfftw,kmfftw])
    c_rho_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw,kmfftw])
    !
    c_A11_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A11_filted, A11_filted,[imfftw,jmfftw,kmfftw])
    c_A12_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A12_filted, A12_filted,[imfftw,jmfftw,kmfftw])
    c_A13_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A13_filted, A13_filted,[imfftw,jmfftw,kmfftw])
    c_A21_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A21_filted, A21_filted,[imfftw,jmfftw,kmfftw])
    c_A22_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A22_filted, A22_filted,[imfftw,jmfftw,kmfftw])
    c_A23_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A23_filted, A23_filted,[imfftw,jmfftw,kmfftw])
    c_A31_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A31_filted, A31_filted,[imfftw,jmfftw,kmfftw])
    c_A32_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A32_filted, A32_filted,[imfftw,jmfftw,kmfftw])
    c_A33_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A33_filted, A33_filted,[imfftw,jmfftw,kmfftw])
    !
    allocate(All_filted_l(1:im,1:jm,1:km),&
            S11_filted_l(1:im,1:jm,1:km),S12_filted_l(1:im,1:jm,1:km),S13_filted_l(1:im,1:jm,1:km),&
            S21_filted_l(1:im,1:jm,1:km),S22_filted_l(1:im,1:jm,1:km),S23_filted_l(1:im,1:jm,1:km),&
            S31_filted_l(1:im,1:jm,1:km),S32_filted_l(1:im,1:jm,1:km),S33_filted_l(1:im,1:jm,1:km))
    !
    allocate(All_filted(1:im,1:jm,1:km),&
            S11_filted(1:im,1:jm,1:km),S12_filted(1:im,1:jm,1:km),S13_filted(1:im,1:jm,1:km),&
            S21_filted(1:im,1:jm,1:km),S22_filted(1:im,1:jm,1:km),S23_filted(1:im,1:jm,1:km),&
            S31_filted(1:im,1:jm,1:km),S32_filted(1:im,1:jm,1:km),S33_filted(1:im,1:jm,1:km),&
            W12_filted(1:im,1:jm,1:km),W21_filted(1:im,1:jm,1:km),&
            W13_filted(1:im,1:jm,1:km),W31_filted(1:im,1:jm,1:km),&
            W23_filted(1:im,1:jm,1:km),W32_filted(1:im,1:jm,1:km))
    !
    c_term1_11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term1_11, term1_11, [imfftw,jmfftw,kmfftw])
    c_term1_12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term1_12, term1_12, [imfftw,jmfftw,kmfftw])
    c_term1_13 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term1_13, term1_13, [imfftw,jmfftw,kmfftw])
    c_term1_21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term1_21, term1_21, [imfftw,jmfftw,kmfftw])
    c_term1_22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term1_22, term1_22, [imfftw,jmfftw,kmfftw])
    c_term1_23 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term1_23, term1_23, [imfftw,jmfftw,kmfftw])
    c_term1_31 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term1_31, term1_31, [imfftw,jmfftw,kmfftw])
    c_term1_32 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term1_32, term1_32, [imfftw,jmfftw,kmfftw])
    c_term1_33 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term1_33, term1_33, [imfftw,jmfftw,kmfftw])
    c_term2    = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term2,    term2,    [imfftw,jmfftw,kmfftw])
    c_term3_11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term3_11, term3_11, [imfftw,jmfftw,kmfftw])
    c_term3_12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term3_12, term3_12, [imfftw,jmfftw,kmfftw])
    c_term3_13 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term3_13, term3_13, [imfftw,jmfftw,kmfftw])
    c_term3_21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term3_21, term3_21, [imfftw,jmfftw,kmfftw])
    c_term3_22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term3_22, term3_22, [imfftw,jmfftw,kmfftw])
    c_term3_23 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term3_23, term3_23, [imfftw,jmfftw,kmfftw])
    c_term3_31 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term3_31, term3_31, [imfftw,jmfftw,kmfftw])
    c_term3_32 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term3_32, term3_32, [imfftw,jmfftw,kmfftw])
    c_term3_33 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term3_33, term3_33, [imfftw,jmfftw,kmfftw])
    c_term4_11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term4_11, term4_11, [imfftw,jmfftw,kmfftw])
    c_term4_12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term4_12, term4_12, [imfftw,jmfftw,kmfftw])
    c_term4_13 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term4_13, term4_13, [imfftw,jmfftw,kmfftw])
    c_term4_21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term4_21, term4_21, [imfftw,jmfftw,kmfftw])
    c_term4_22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term4_22, term4_22, [imfftw,jmfftw,kmfftw])
    c_term4_23 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term4_23, term4_23, [imfftw,jmfftw,kmfftw])
    c_term4_31 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term4_31, term4_31, [imfftw,jmfftw,kmfftw])
    c_term4_32 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term4_32, term4_32, [imfftw,jmfftw,kmfftw])
    c_term4_33 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term4_33, term4_33, [imfftw,jmfftw,kmfftw])
    c_term5    = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term5,    term5,    [imfftw,jmfftw,kmfftw])
    c_term6_11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term6_11, term6_11, [imfftw,jmfftw,kmfftw])
    c_term6_12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term6_12, term6_12, [imfftw,jmfftw,kmfftw])
    c_term6_13 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term6_13, term6_13, [imfftw,jmfftw,kmfftw])
    c_term6_21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term6_21, term6_21, [imfftw,jmfftw,kmfftw])
    c_term6_22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term6_22, term6_22, [imfftw,jmfftw,kmfftw])
    c_term6_23 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term6_23, term6_23, [imfftw,jmfftw,kmfftw])
    c_term6_31 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term6_31, term6_31, [imfftw,jmfftw,kmfftw])
    c_term6_32 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term6_32, term6_32, [imfftw,jmfftw,kmfftw])
    c_term6_33 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term6_33, term6_33, [imfftw,jmfftw,kmfftw])
    c_term7    = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_term7,    term7,    [imfftw,jmfftw,kmfftw])
    !
    !
    Pi1 = 0.d0
    Pi2 =	0.d0
    Pi3 =	0.d0
    Pi4 =	0.d0
    Pi5 =	0.d0
    Pi6 =	0.d0
    Pi7 =	0.d0
    !
    if(mpirank==0)  print *, "Array allocated and initialized"
    !
    do m=1,num_l
      !
      !!!!!! Filter to get Sij filted by l
      if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
      !
      if(mpirank == 0) then
        write(mname,'(i4.4)')m
        if (thefilenumb .ne. 0) then
          outfilename2 = 'pp/SGS_Pi_precise_'//stepname//'_'//mname//'.dat'
        else
          outfilename2 = 'pp/SGS_Pi_precise_'//mname//'.dat'
        endif
        call listinit(filename=outfilename2,handle=hand_b, &
                    firstline='nstep time sqrtalpha pi1 pi2 pi3 pi4 pi5 pi6 pi7')
      endif
      !
      !!!! Velocity Favre average and density average
      ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
      do k=1,km
      do j=1,jm
      do i=1,im
        Gl = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
        !
        w1_filted(i,j,k)    = w1(i,j,k)    *Gl
        w2_filted(i,j,k)    = w2(i,j,k)    *Gl
        w3_filted(i,j,k)    = w3(i,j,k)    *Gl
        !
        rho_filted(i,j,k)   = rhocom(i,j,k)*Gl
      enddo
      enddo
      enddo
      !
      ! After this bloc, w1_filted is (rho*u1)_filted in physical space
      call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
      call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
      call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
      !
      ! After this bloc, w1_filted is u1_filted in physical space
      do k=1,km
      do j=1,jm
      do i=1,im
        w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
        w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
        w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
      enddo
      enddo
      enddo
      !
      ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
      call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
      call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
      call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1_filted(i,j,k)  = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
        w2_filted(i,j,k)  = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
        w3_filted(i,j,k)  = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
        !
        A11_filted(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
        A21_filted(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
        A31_filted(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
        A12_filted(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
        A22_filted(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
        A32_filted(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
        A13_filted(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
        A23_filted(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
        A33_filted(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
        !
      end do
      end do
      end do
      !
      !
      !
      ! After this bloc, A11_filted is A11_filted in physical space
      call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
      call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
      call fftw_mpi_execute_dft(backward_plan,A31_filted,A31_filted)
      call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
      call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
      call fftw_mpi_execute_dft(backward_plan,A32_filted,A32_filted)
      call fftw_mpi_execute_dft(backward_plan,A13_filted,A13_filted)
      call fftw_mpi_execute_dft(backward_plan,A23_filted,A23_filted)
      call fftw_mpi_execute_dft(backward_plan,A33_filted,A33_filted)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        All_filted_l(i,j,k) = A11_filted(i,j,k)+A22_filted(i,j,k)+A33_filted(i,j,k)
        !
        S11_filted_l(i,j,k) = A11_filted(i,j,k) - 1.d0/3.d0 * All_filted_l(i,j,k)
        S22_filted_l(i,j,k) = A22_filted(i,j,k) - 1.d0/3.d0 * All_filted_l(i,j,k)
        S33_filted_l(i,j,k) = A33_filted(i,j,k) - 1.d0/3.d0 * All_filted_l(i,j,k)
        S12_filted_l(i,j,k) = (A12_filted(i,j,k) + A21_filted(i,j,k))*0.5d0
        S21_filted_l(i,j,k) = S12_filted_l(i,j,k)
        S13_filted_l(i,j,k) = (A13_filted(i,j,k) + A31_filted(i,j,k))*0.5d0
        S31_filted_l(i,j,k) = S13_filted_l(i,j,k)
        S23_filted_l(i,j,k) = (A23_filted(i,j,k) + A32_filted(i,j,k))*0.5d0
        S32_filted_l(i,j,k) = S23_filted_l(i,j,k)
        !
      end do
      end do
      end do
      !
      if(mpirank==0)  print *, '** l filted!'
      !
      !!!!!! Begin integral
      !
      do n=1,num_alphas(m)
        !
        call date_and_time(values=value) 
        !
        if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alphas(m),',now is ',&
                                value(5), ':', value(6),':',value(7)
        !!!! Velocity Favre average and density average
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        do k=1,km
        do j=1,jm
        do i=1,im
          Galpha = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_sqrtalpha(m,n)**2/2.d0) ! Filtre scale :sqrtalpha
          w1_filted(i,j,k)  = w1(i,j,k)    *Galpha
          w2_filted(i,j,k)  = w2(i,j,k)    *Galpha
          w3_filted(i,j,k)  = w3(i,j,k)    *Galpha
          rho_filted(i,j,k) = rhocom(i,j,k)*Galpha
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        do k=1,km
        do j=1,jm
        do i=1,im
          w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
          w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
          w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j,k)  = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
          w2_filted(i,j,k)  = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
          w3_filted(i,j,k)  = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
          !
          A11_filted(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
          A21_filted(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
          A31_filted(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
          A12_filted(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
          A22_filted(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
          A32_filted(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
          A13_filted(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
          A23_filted(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
          A33_filted(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
          !
        end do
        end do
        end do
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
        call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
        call fftw_mpi_execute_dft(backward_plan,A31_filted,A31_filted)
        call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
        call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
        call fftw_mpi_execute_dft(backward_plan,A32_filted,A32_filted)
        call fftw_mpi_execute_dft(backward_plan,A13_filted,A13_filted)
        call fftw_mpi_execute_dft(backward_plan,A23_filted,A23_filted)
        call fftw_mpi_execute_dft(backward_plan,A33_filted,A33_filted)
        !
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          All_filted(i,j,k) = A11_filted(i,j,k)+A22_filted(i,j,k)+A33_filted(i,j,k)
          !
          S11_filted(i,j,k) = A11_filted(i,j,k) - 1.d0/3.d0 * All_filted(i,j,k)
          S22_filted(i,j,k) = A22_filted(i,j,k) - 1.d0/3.d0 * All_filted(i,j,k)
          S33_filted(i,j,k) = A33_filted(i,j,k) - 1.d0/3.d0 * All_filted(i,j,k)
          S12_filted(i,j,k) = (A12_filted(i,j,k) + A21_filted(i,j,k))*0.5d0
          S21_filted(i,j,k) = S12_filted(i,j,k)
          S13_filted(i,j,k) = (A13_filted(i,j,k) + A31_filted(i,j,k))*0.5d0
          S31_filted(i,j,k) = S13_filted(i,j,k)
          S23_filted(i,j,k) = (A23_filted(i,j,k) + A32_filted(i,j,k))*0.5d0
          S32_filted(i,j,k) = S23_filted(i,j,k)
          !
          W12_filted(i,j,k) = (A12_filted(i,j,k)-A21_filted(i,j,k))*0.5d0
          W21_filted(i,j,k) = -1.d0*W12_filted(i,j,k)
          W13_filted(i,j,k) = (A13_filted(i,j,k)-A31_filted(i,j,k))*0.5d0
          W31_filted(i,j,k) = -1.d0*W13_filted(i,j,k)
          W23_filted(i,j,k) = (A23_filted(i,j,k)-A32_filted(i,j,k))*0.5d0
          W32_filted(i,j,k) = -1.d0*W23_filted(i,j,k)
          !
        end do
        end do
        end do
        !
        !!!! Pi terms
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !term1_IJ = rho_filted*SI1_filted*SJ1_filted + rho_filted*SI2_filted*SJ2_filted + rho_filted*SI3_filted*SJ3_filted
          term1_11(i,j,k) = rho_filted(i,j,k)*S11_filted(i,j,k)*S11_filted(i,j,k) + &
                            rho_filted(i,j,k)*S12_filted(i,j,k)*S12_filted(i,j,k) + &
                            rho_filted(i,j,k)*S13_filted(i,j,k)*S13_filted(i,j,k)
          term1_12(i,j,k) = rho_filted(i,j,k)*S11_filted(i,j,k)*S21_filted(i,j,k) + &
                            rho_filted(i,j,k)*S12_filted(i,j,k)*S22_filted(i,j,k) + &
                            rho_filted(i,j,k)*S13_filted(i,j,k)*S23_filted(i,j,k)
          term1_13(i,j,k) = rho_filted(i,j,k)*S11_filted(i,j,k)*S31_filted(i,j,k) + &
                            rho_filted(i,j,k)*S12_filted(i,j,k)*S32_filted(i,j,k) + &
                            rho_filted(i,j,k)*S13_filted(i,j,k)*S33_filted(i,j,k)
          term1_21(i,j,k) = rho_filted(i,j,k)*S21_filted(i,j,k)*S11_filted(i,j,k) + &
                            rho_filted(i,j,k)*S22_filted(i,j,k)*S12_filted(i,j,k) + &
                            rho_filted(i,j,k)*S23_filted(i,j,k)*S13_filted(i,j,k)
          term1_22(i,j,k) = rho_filted(i,j,k)*S21_filted(i,j,k)*S21_filted(i,j,k) + &
                            rho_filted(i,j,k)*S22_filted(i,j,k)*S22_filted(i,j,k) + &
                            rho_filted(i,j,k)*S23_filted(i,j,k)*S23_filted(i,j,k)
          term1_23(i,j,k) = rho_filted(i,j,k)*S21_filted(i,j,k)*S31_filted(i,j,k) + &
                            rho_filted(i,j,k)*S22_filted(i,j,k)*S32_filted(i,j,k) + &
                            rho_filted(i,j,k)*S23_filted(i,j,k)*S33_filted(i,j,k)
          term1_31(i,j,k) = rho_filted(i,j,k)*S31_filted(i,j,k)*S11_filted(i,j,k) + &
                            rho_filted(i,j,k)*S32_filted(i,j,k)*S12_filted(i,j,k) + &
                            rho_filted(i,j,k)*S33_filted(i,j,k)*S13_filted(i,j,k)
          term1_32(i,j,k) = rho_filted(i,j,k)*S31_filted(i,j,k)*S21_filted(i,j,k) + &
                            rho_filted(i,j,k)*S32_filted(i,j,k)*S22_filted(i,j,k) + &
                            rho_filted(i,j,k)*S33_filted(i,j,k)*S23_filted(i,j,k)
          term1_33(i,j,k) = rho_filted(i,j,k)*S31_filted(i,j,k)*S31_filted(i,j,k) + &
                            rho_filted(i,j,k)*S32_filted(i,j,k)*S32_filted(i,j,k) + &
                            rho_filted(i,j,k)*S33_filted(i,j,k)*S33_filted(i,j,k)
          ! 
          ! term2 
          term2(i,j,k) =  rho_filted(i,j,k)*S11_filted(i,j,k)*S11_filted(i,j,k)+&
                          rho_filted(i,j,k)*S12_filted(i,j,k)*S12_filted(i,j,k)+&
                          rho_filted(i,j,k)*S13_filted(i,j,k)*S13_filted(i,j,k)+& ! i=1,k=1,2,3
                          rho_filted(i,j,k)*S21_filted(i,j,k)*S21_filted(i,j,k)+&
                          rho_filted(i,j,k)*S22_filted(i,j,k)*S22_filted(i,j,k)+&
                          rho_filted(i,j,k)*S23_filted(i,j,k)*S23_filted(i,j,k)+& ! i=2,k=1,2,3
                          rho_filted(i,j,k)*S31_filted(i,j,k)*S31_filted(i,j,k)+&
                          rho_filted(i,j,k)*S32_filted(i,j,k)*S32_filted(i,j,k)+&
                          rho_filted(i,j,k)*S33_filted(i,j,k)*S33_filted(i,j,k) ! i=3,k=1,2,3
          !
          ! term3_IJ = rho_filted*All_filted*SIJ_filted
          term3_11(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S11_filted(i,j,k)
          term3_12(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S12_filted(i,j,k)
          term3_13(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S13_filted(i,j,k)
          !
          term3_21(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S21_filted(i,j,k)
          term3_22(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S22_filted(i,j,k)
          term3_23(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S23_filted(i,j,k)
          !
          term3_31(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S31_filted(i,j,k)
          term3_32(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S32_filted(i,j,k)
          term3_33(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S33_filted(i,j,k)
          !
          !term4_IJ = rho_filted*WI1_filted*W1J_filted+rho_filted*WI2_filted*W2J_filted + &
          !rho_filted*WI3_filted*W3J_filted 
          term4_11(i,j,k) = rho_filted(i,j,k)*W12_filted(i,j,k)*W21_filted(i,j,k) + &
                            rho_filted(i,j,k)*W13_filted(i,j,k)*W31_filted(i,j,k) 
          term4_21(i,j,k) = rho_filted(i,j,k)*W23_filted(i,j,k)*W31_filted(i,j,k) 
          term4_31(i,j,k) = rho_filted(i,j,k)*W32_filted(i,j,k)*W21_filted(i,j,k)
          !
          term4_12(i,j,k) = rho_filted(i,j,k)*W13_filted(i,j,k)*W32_filted(i,j,k) 
          term4_22(i,j,k) = rho_filted(i,j,k)*W21_filted(i,j,k)*W12_filted(i,j,k) + &
                            rho_filted(i,j,k)*W23_filted(i,j,k)*W32_filted(i,j,k) 
          term4_32(i,j,k) = rho_filted(i,j,k)*W31_filted(i,j,k)*W12_filted(i,j,k)
          !
          term4_13(i,j,k) = rho_filted(i,j,k)*W12_filted(i,j,k)*W23_filted(i,j,k)
          term4_23(i,j,k) = rho_filted(i,j,k)*W21_filted(i,j,k)*W13_filted(i,j,k)
          term4_33(i,j,k) = rho_filted(i,j,k)*W31_filted(i,j,k)*W13_filted(i,j,k) + &
                            rho_filted(i,j,k)*W32_filted(i,j,k)*W23_filted(i,j,k)
          !
          ! term5
          term5(i,j,k) = 2.d0*rho_filted(i,j,k)*(W12_filted(i,j,k)*W21_filted(i,j,k) + &
                          W13_filted(i,j,k)*W31_filted(i,j,k) + W23_filted(i,j,k)*W32_filted(i,j,k))
          !
          !term6_IJ= rho_filted*(S1J_filted*WI1_filted-SI1_filted*W1J_filted) + &
          !          rho_filted*(S2J_filted*WI2_filted-SI2_filted*W2J_filted) + &
          !          rho_filted*(S3J_filted*WI3_filted-SI3_filted*W3J_filted)
          term6_11(i,j,k)= rho_filted(i,j,k)*S21_filted(i,j,k)*W12_filted(i,j,k) &
                          -rho_filted(i,j,k)*S12_filted(i,j,k)*W21_filted(i,j,k) &
                          +rho_filted(i,j,k)*S31_filted(i,j,k)*W13_filted(i,j,k) &
                          -rho_filted(i,j,k)*S13_filted(i,j,k)*W31_filted(i,j,k)
          term6_21(i,j,k)= rho_filted(i,j,k)*S11_filted(i,j,k)*W21_filted(i,j,k) &
                          -rho_filted(i,j,k)*S22_filted(i,j,k)*W21_filted(i,j,k) &
                          +rho_filted(i,j,k)*S31_filted(i,j,k)*W23_filted(i,j,k) &
                          -rho_filted(i,j,k)*S23_filted(i,j,k)*W31_filted(i,j,k)
          term6_31(i,j,k)= rho_filted(i,j,k)*S11_filted(i,j,k)*W31_filted(i,j,k) &
                          +rho_filted(i,j,k)*S21_filted(i,j,k)*W32_filted(i,j,k) &
                          -rho_filted(i,j,k)*S32_filted(i,j,k)*W21_filted(i,j,k) &
                          -rho_filted(i,j,k)*S33_filted(i,j,k)*W31_filted(i,j,k) 
          !
          term6_12(i,j,k)=-rho_filted(i,j,k)*S11_filted(i,j,k)*W12_filted(i,j,k) &
                          +rho_filted(i,j,k)*S22_filted(i,j,k)*W12_filted(i,j,k) &
                          +rho_filted(i,j,k)*S32_filted(i,j,k)*W13_filted(i,j,k) &
                          -rho_filted(i,j,k)*S13_filted(i,j,k)*W32_filted(i,j,k)
          term6_22(i,j,k)= rho_filted(i,j,k)*S12_filted(i,j,k)*W21_filted(i,j,k) &
                          -rho_filted(i,j,k)*S21_filted(i,j,k)*W12_filted(i,j,k) &
                          +rho_filted(i,j,k)*S32_filted(i,j,k)*W23_filted(i,j,k) &
                          -rho_filted(i,j,k)*S23_filted(i,j,k)*W32_filted(i,j,k)
          term6_32(i,j,k)= rho_filted(i,j,k)*S12_filted(i,j,k)*W31_filted(i,j,k) &
                          -rho_filted(i,j,k)*S31_filted(i,j,k)*W12_filted(i,j,k) &
                          +rho_filted(i,j,k)*S22_filted(i,j,k)*W32_filted(i,j,k) &
                          -rho_filted(i,j,k)*S33_filted(i,j,k)*W32_filted(i,j,k)
          !
          term6_13(i,j,k)=-rho_filted(i,j,k)*S11_filted(i,j,k)*W13_filted(i,j,k) &
                          +rho_filted(i,j,k)*S23_filted(i,j,k)*W12_filted(i,j,k) &
                          -rho_filted(i,j,k)*S12_filted(i,j,k)*W23_filted(i,j,k) &
                          +rho_filted(i,j,k)*S33_filted(i,j,k)*W13_filted(i,j,k)
          term6_23(i,j,k)= rho_filted(i,j,k)*S13_filted(i,j,k)*W21_filted(i,j,k) &
                          -rho_filted(i,j,k)*S21_filted(i,j,k)*W13_filted(i,j,k) &
                          -rho_filted(i,j,k)*S22_filted(i,j,k)*W23_filted(i,j,k) &
                          +rho_filted(i,j,k)*S33_filted(i,j,k)*W23_filted(i,j,k)
          term6_33(i,j,k)= rho_filted(i,j,k)*S13_filted(i,j,k)*W31_filted(i,j,k) &
                          -rho_filted(i,j,k)*S31_filted(i,j,k)*W13_filted(i,j,k) &
                          +rho_filted(i,j,k)*S23_filted(i,j,k)*W32_filted(i,j,k) &
                          -rho_filted(i,j,k)*S32_filted(i,j,k)*W23_filted(i,j,k)
          !
          ! term7
          term7(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*All_filted(i,j,k)
        enddo
        enddo
        enddo
        !
        ! Do filter phi:
        ! F -> product -> F inverse
        call fftw_mpi_execute_dft(forward_plan,term1_11,term1_11)
        call fftw_mpi_execute_dft(forward_plan,term1_12,term1_12)
        call fftw_mpi_execute_dft(forward_plan,term1_13,term1_13)
        call fftw_mpi_execute_dft(forward_plan,term1_21,term1_21)
        call fftw_mpi_execute_dft(forward_plan,term1_22,term1_22)
        call fftw_mpi_execute_dft(forward_plan,term1_23,term1_23)
        call fftw_mpi_execute_dft(forward_plan,term1_31,term1_31)
        call fftw_mpi_execute_dft(forward_plan,term1_32,term1_32)
        call fftw_mpi_execute_dft(forward_plan,term1_33,term1_33)
        !
        call fftw_mpi_execute_dft(forward_plan,term2   ,term2   )
        !
        call fftw_mpi_execute_dft(forward_plan,term3_11,term3_11)
        call fftw_mpi_execute_dft(forward_plan,term3_12,term3_12)
        call fftw_mpi_execute_dft(forward_plan,term3_13,term3_13)
        call fftw_mpi_execute_dft(forward_plan,term3_21,term3_21)
        call fftw_mpi_execute_dft(forward_plan,term3_22,term3_22)
        call fftw_mpi_execute_dft(forward_plan,term3_23,term3_23)
        call fftw_mpi_execute_dft(forward_plan,term3_31,term3_31)
        call fftw_mpi_execute_dft(forward_plan,term3_32,term3_32)
        call fftw_mpi_execute_dft(forward_plan,term3_33,term3_33)
        !
        call fftw_mpi_execute_dft(forward_plan,term4_11,term4_11)
        call fftw_mpi_execute_dft(forward_plan,term4_12,term4_12)
        call fftw_mpi_execute_dft(forward_plan,term4_13,term4_13)
        call fftw_mpi_execute_dft(forward_plan,term4_21,term4_21)
        call fftw_mpi_execute_dft(forward_plan,term4_22,term4_22)
        call fftw_mpi_execute_dft(forward_plan,term4_23,term4_23)
        call fftw_mpi_execute_dft(forward_plan,term4_31,term4_31)
        call fftw_mpi_execute_dft(forward_plan,term4_32,term4_32)
        call fftw_mpi_execute_dft(forward_plan,term4_33,term4_33)
        !
        call fftw_mpi_execute_dft(forward_plan,term5   ,term5   )
        !
        call fftw_mpi_execute_dft(forward_plan,term6_11,term6_11)
        call fftw_mpi_execute_dft(forward_plan,term6_12,term6_12)
        call fftw_mpi_execute_dft(forward_plan,term6_13,term6_13)
        call fftw_mpi_execute_dft(forward_plan,term6_21,term6_21)
        call fftw_mpi_execute_dft(forward_plan,term6_22,term6_22)
        call fftw_mpi_execute_dft(forward_plan,term6_23,term6_23)
        call fftw_mpi_execute_dft(forward_plan,term6_31,term6_31)
        call fftw_mpi_execute_dft(forward_plan,term6_32,term6_32)
        call fftw_mpi_execute_dft(forward_plan,term6_33,term6_33)
        !
        call fftw_mpi_execute_dft(forward_plan,term7   ,term7   )
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          Gphi = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_phi(m,n)**2/2.d0) ! Filtre scale :phi
          term1_11(i,j,k) = term1_11(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term1_12(i,j,k) = term1_12(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term1_13(i,j,k) = term1_13(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term1_21(i,j,k) = term1_21(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term1_22(i,j,k) = term1_22(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term1_23(i,j,k) = term1_23(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term1_31(i,j,k) = term1_31(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term1_32(i,j,k) = term1_32(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term1_33(i,j,k) = term1_33(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          !
          term2(i,j,k)    = term2(i,j,k)   *Gphi/(1.d0*ia*ja*ka)
          !
          term3_11(i,j,k) = term3_11(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term3_12(i,j,k) = term3_12(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term3_13(i,j,k) = term3_13(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term3_21(i,j,k) = term3_21(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term3_22(i,j,k) = term3_22(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term3_23(i,j,k) = term3_23(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term3_31(i,j,k) = term3_31(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term3_32(i,j,k) = term3_32(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term3_33(i,j,k) = term3_33(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          !
          term4_11(i,j,k) = term4_11(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term4_12(i,j,k) = term4_12(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term4_13(i,j,k) = term4_13(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term4_21(i,j,k) = term4_21(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term4_22(i,j,k) = term4_22(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term4_23(i,j,k) = term4_23(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term4_31(i,j,k) = term4_31(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term4_32(i,j,k) = term4_32(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term4_33(i,j,k) = term4_33(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          !
          term5(i,j,k)    = term5(i,j,k)   *Gphi/(1.d0*ia*ja*ka)
          !
          term6_11(i,j,k) = term6_11(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term6_12(i,j,k) = term6_12(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term6_13(i,j,k) = term6_13(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term6_21(i,j,k) = term6_21(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term6_22(i,j,k) = term6_22(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term6_23(i,j,k) = term6_23(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term6_31(i,j,k) = term6_31(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term6_32(i,j,k) = term6_32(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          term6_33(i,j,k) = term6_33(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          !
          term7(i,j,k)    = term7(i,j,k)   *Gphi/(1.d0*ia*ja*ka)
          !
        enddo
        enddo
        enddo
        !
        !
        call fftw_mpi_execute_dft(backward_plan,term1_11,term1_11)
        call fftw_mpi_execute_dft(backward_plan,term1_12,term1_12)
        call fftw_mpi_execute_dft(backward_plan,term1_13,term1_13)
        call fftw_mpi_execute_dft(backward_plan,term1_21,term1_21)
        call fftw_mpi_execute_dft(backward_plan,term1_22,term1_22)
        call fftw_mpi_execute_dft(backward_plan,term1_23,term1_23)
        call fftw_mpi_execute_dft(backward_plan,term1_31,term1_31)
        call fftw_mpi_execute_dft(backward_plan,term1_32,term1_32)
        call fftw_mpi_execute_dft(backward_plan,term1_33,term1_33)
        !
        call fftw_mpi_execute_dft(backward_plan,term2   ,term2   )
        !
        call fftw_mpi_execute_dft(backward_plan,term3_11,term3_11)
        call fftw_mpi_execute_dft(backward_plan,term3_12,term3_12)
        call fftw_mpi_execute_dft(backward_plan,term3_13,term3_13)
        call fftw_mpi_execute_dft(backward_plan,term3_21,term3_21)
        call fftw_mpi_execute_dft(backward_plan,term3_22,term3_22)
        call fftw_mpi_execute_dft(backward_plan,term3_23,term3_23)
        call fftw_mpi_execute_dft(backward_plan,term3_31,term3_31)
        call fftw_mpi_execute_dft(backward_plan,term3_32,term3_32)
        call fftw_mpi_execute_dft(backward_plan,term3_33,term3_33)
        !
        call fftw_mpi_execute_dft(backward_plan,term4_11,term4_11)
        call fftw_mpi_execute_dft(backward_plan,term4_12,term4_12)
        call fftw_mpi_execute_dft(backward_plan,term4_13,term4_13)
        call fftw_mpi_execute_dft(backward_plan,term4_21,term4_21)
        call fftw_mpi_execute_dft(backward_plan,term4_22,term4_22)
        call fftw_mpi_execute_dft(backward_plan,term4_23,term4_23)
        call fftw_mpi_execute_dft(backward_plan,term4_31,term4_31)
        call fftw_mpi_execute_dft(backward_plan,term4_32,term4_32)
        call fftw_mpi_execute_dft(backward_plan,term4_33,term4_33)
        !
        call fftw_mpi_execute_dft(backward_plan,term5   ,term5   )
        !
        call fftw_mpi_execute_dft(backward_plan,term6_11,term6_11)
        call fftw_mpi_execute_dft(backward_plan,term6_12,term6_12)
        call fftw_mpi_execute_dft(backward_plan,term6_13,term6_13)
        call fftw_mpi_execute_dft(backward_plan,term6_21,term6_21)
        call fftw_mpi_execute_dft(backward_plan,term6_22,term6_22)
        call fftw_mpi_execute_dft(backward_plan,term6_23,term6_23)
        call fftw_mpi_execute_dft(backward_plan,term6_31,term6_31)
        call fftw_mpi_execute_dft(backward_plan,term6_32,term6_32)
        call fftw_mpi_execute_dft(backward_plan,term6_33,term6_33)
        !
        call fftw_mpi_execute_dft(backward_plan,term7   ,term7   )
        !
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          vxr_D1 = real(term1_11(i,j,k) * S11_filted_l(i,j,k) + &
                  term1_12(i,j,k) * S12_filted_l(i,j,k) + &
                  term1_13(i,j,k) * S13_filted_l(i,j,k) + &
                  term1_21(i,j,k) * S21_filted_l(i,j,k) + &
                  term1_22(i,j,k) * S22_filted_l(i,j,k) + &
                  term1_23(i,j,k) * S23_filted_l(i,j,k) + &
                  term1_31(i,j,k) * S31_filted_l(i,j,k) + &
                  term1_32(i,j,k) * S32_filted_l(i,j,k) + &
                  term1_33(i,j,k) * S33_filted_l(i,j,k))
          Pi1(m) = Pi1(m) + vxr_D1 * dl_alpha(m,n)
          !
          vxr_D2 = real(term2(i,j,k) * All_filted_l(i,j,k))
          Pi2(m) = Pi2(m) + vxr_D2 * dl_alpha(m,n) / 3.d0
          !
          vxr_D3 = real(term3_11(i,j,k) * S11_filted_l(i,j,k) + &
                  term3_12(i,j,k) * S12_filted_l(i,j,k) + &
                  term3_13(i,j,k) * S13_filted_l(i,j,k) + &
                  term3_21(i,j,k) * S21_filted_l(i,j,k) + &
                  term3_22(i,j,k) * S22_filted_l(i,j,k) + &
                  term3_23(i,j,k) * S23_filted_l(i,j,k) + &
                  term3_31(i,j,k) * S31_filted_l(i,j,k) + &
                  term3_32(i,j,k) * S32_filted_l(i,j,k) + &
                  term3_33(i,j,k) * S33_filted_l(i,j,k))
          Pi3(m) = Pi3(m) + vxr_D3 * dl_alpha(m,n) * 2.d0/3.d0
          !
          vxr_D4 = real(term4_11(i,j,k) * S11_filted_l(i,j,k) + &
                  term4_12(i,j,k) * S12_filted_l(i,j,k) + &
                  term4_13(i,j,k) * S13_filted_l(i,j,k) + &
                  term4_21(i,j,k) * S21_filted_l(i,j,k) + &
                  term4_22(i,j,k) * S22_filted_l(i,j,k) + &
                  term4_23(i,j,k) * S23_filted_l(i,j,k) + &
                  term4_31(i,j,k) * S31_filted_l(i,j,k) + &
                  term4_32(i,j,k) * S32_filted_l(i,j,k) + &
                  term4_33(i,j,k) * S33_filted_l(i,j,k))
          Pi4(m) = Pi4(m) - vxr_D4 * dl_alpha(m,n)
          !
          vxr_D5 = real(term5(i,j,k) * All_filted_l(i,j,k))
          Pi5(m) = Pi5(m) - vxr_D5 * dl_alpha(m,n) / 3.d0
          !
          vxr_D6 = real(term6_11(i,j,k) * S11_filted_l(i,j,k) + &
                  term6_12(i,j,k) * S12_filted_l(i,j,k) + &
                  term6_13(i,j,k) * S13_filted_l(i,j,k) + &
                  term6_21(i,j,k) * S21_filted_l(i,j,k) + &
                  term6_22(i,j,k) * S22_filted_l(i,j,k) + &
                  term6_23(i,j,k) * S23_filted_l(i,j,k) + &
                  term6_31(i,j,k) * S31_filted_l(i,j,k) + &
                  term6_32(i,j,k) * S32_filted_l(i,j,k) + &
                  term6_33(i,j,k) * S33_filted_l(i,j,k))
          Pi6(m) = Pi6(m) + vxr_D6 * dl_alpha(m,n)
          !
          vxr_D7 = term7(i,j,k) * All_filted_l(i,j,k)
          Pi7(m) = Pi7(m) + vxr_D7 * dl_alpha(m,n) / 9.d0
        enddo
        enddo
        enddo
        !
        vxr_D1 = psum(vxr_D1)
        vxr_D2 = psum(vxr_D2)
        vxr_D3 = psum(vxr_D3)
        vxr_D4 = psum(vxr_D4)
        vxr_D5 = psum(vxr_D5)
        vxr_D6 = psum(vxr_D6)
        vxr_D7 = psum(vxr_D7)
        !
        if(mpirank==0) then
          call listwrite(hand_b,l_sqrtalpha(m,n),vxr_D1 * dl_alpha(m,n), &
                        vxr_D2 * dl_alpha(m,n) / 3.d0, &
                        vxr_D3 * dl_alpha(m,n) * 2.d0/3.d0, &
                        - vxr_D4 * dl_alpha(m,n), &
                        - vxr_D5 * dl_alpha(m,n) / 3.d0, &
                        vxr_D6 * dl_alpha(m,n), &
                        vxr_D7 * dl_alpha(m,n) / 9.d0)
        endif
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
      enddo
      !
      Pi1(m) =	 psum(Pi1(m)) / (ia*ja*ka)
      Pi2(m) =	 psum(Pi2(m)) / (ia*ja*ka)
      Pi3(m) =	 psum(Pi3(m)) / (ia*ja*ka)
      Pi4(m) =	 psum(Pi4(m)) / (ia*ja*ka)
      Pi5(m) =	 psum(Pi5(m)) / (ia*ja*ka)
      Pi6(m) =	 psum(Pi6(m)) / (ia*ja*ka)
      Pi7(m) =	 psum(Pi7(m)) / (ia*ja*ka)
      !
      if(mpirank==0)then
        print *, '>>>>', outfilename2
      endif
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
    enddo
    if(mpirank==0)  print *, 'Job finish'
    !
    if(mpirank==0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/SGS_Pi_'//stepname//'.dat'
      else
        outfilename = 'pp/SGS_Pi.dat'
      endif
      
      call listinit(filename=outfilename,handle=hand_a, &
                    firstline='nstep time lOlmin pi1 pi2 pi3 pi4 pi5 pi6 pi7')
      do m=1,num_l
        call listwrite(hand_a,l_lim(m)/l_min, Pi1(m), Pi2(m), &
                        Pi3(m), Pi4(m), Pi5(m),     &
                        Pi6(m), Pi7(m))
      enddo
      !
      print *, '>>>>', outfilename
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_w1)
    call fftw_free(c_w2)
    call fftw_free(c_w3)
    call fftw_free(c_rhocom)
    call fftw_free(c_w1_filted)
    call fftw_free(c_w2_filted)
    call fftw_free(c_w3_filted)
    call fftw_free(c_rho_filted)
    call fftw_free(c_A11_filted)
    call fftw_free(c_A12_filted)
    call fftw_free(c_A13_filted)
    call fftw_free(c_A21_filted)
    call fftw_free(c_A22_filted)
    call fftw_free(c_A23_filted)
    call fftw_free(c_A31_filted)
    call fftw_free(c_A32_filted)
    call fftw_free(c_A33_filted)
    call fftw_free(c_term1_11)
    call fftw_free(c_term1_12)
    call fftw_free(c_term1_13)
    call fftw_free(c_term1_21)
    call fftw_free(c_term1_22)
    call fftw_free(c_term1_23)
    call fftw_free(c_term1_31)
    call fftw_free(c_term1_32)
    call fftw_free(c_term1_33)
    call fftw_free(c_term2)
    call fftw_free(c_term5)
    call fftw_free(c_term7)
    call fftw_free(c_term3_11)
    call fftw_free(c_term3_12)
    call fftw_free(c_term3_13)
    call fftw_free(c_term3_21)
    call fftw_free(c_term3_22)
    call fftw_free(c_term3_23)
    call fftw_free(c_term3_31)
    call fftw_free(c_term3_32)
    call fftw_free(c_term3_33)
    call fftw_free(c_term4_11)
    call fftw_free(c_term4_12)
    call fftw_free(c_term4_13)
    call fftw_free(c_term4_21)
    call fftw_free(c_term4_22)
    call fftw_free(c_term4_23)
    call fftw_free(c_term4_31)
    call fftw_free(c_term4_32)
    call fftw_free(c_term4_33)
    call fftw_free(c_term6_11)
    call fftw_free(c_term6_12)
    call fftw_free(c_term6_13)
    call fftw_free(c_term6_21)
    call fftw_free(c_term6_22)
    call fftw_free(c_term6_23)
    call fftw_free(c_term6_31)
    call fftw_free(c_term6_32)
    call fftw_free(c_term6_33)
    call mpistop
    deallocate(All_filted_l,S11_filted_l,S12_filted_l,S13_filted_l)
    deallocate(S21_filted_l,S22_filted_l,S23_filted_l)
    deallocate(S31_filted_l,S32_filted_l,S33_filted_l)
    deallocate(All_filted,S11_filted,S12_filted,S13_filted)
    deallocate(S21_filted,S22_filted,S23_filted)
    deallocate(S31_filted,S32_filted,S33_filted)
    deallocate(W12_filted,W21_filted,W13_filted,W31_filted,W23_filted,W32_filted)
    deallocate(k1,k2,k3)
    deallocate(l_lim,l_sqrtalpha,l_phi,dl_alpha)
    deallocate(Pi1,Pi2,Pi3,Pi4,Pi5,Pi6,Pi7)
    !
  end subroutine SGSPi3Dint
  !
  subroutine SGSstress2D(thefilenumb)
    ! 
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km
    use commarray, only: vel, rho
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
    include 'fftw3-mpi.f03'
    !
    integer,intent(in) :: thefilenumb
    integer :: fh
    integer :: i,j,m,n,p,q
    character(len=128) :: infilename,outfilename
    character(len=4) :: stepname,mname
    character(len=10):: termname
    real(8), allocatable, dimension(:,:) :: k1,k2
    complex(8) :: imag
    real(8),allocatable,dimension(:) :: l_lim
    real(8),allocatable,dimension(:,:) :: l_sqrtalpha,l_phi,dl_alpha
    integer,allocatable,dimension(:) :: num_alphas
    integer :: num_l,num_alpha,num_alphamin
    integer :: hand_a
    real(8) :: l_min, ratio_max, ratio_min
    real(8) :: Gl,Galpha,Gphi
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1,w2,rhol
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1_filted,w2_filted,rho_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1w1,w1w2,w2w1,w2w2
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1w1_filted,w1w2_filted,w2w1_filted,w2w2_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: A11,A12,A21,A22
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: tau11_term,tau12_term,tau21_term,tau22_term
    real(8), allocatable, dimension(:,:,:,:) :: tau ! 1:2,1:2,1:im,1:jm,1:km
    real(8), allocatable, dimension(:,:,:,:) :: tau_bis ! 1:2,1:2,1:im,1:jm,1:km
    real(8), allocatable, dimension(:,:) :: errormax,erroravg, errorgtr10,errorgtr100
    real(8) :: result,norm2,norm2bis
    real(8) :: errornorm2max,errornorm2avg,errornorm2gtr10,errornorm2gtr100
    !
    !
    type(C_PTR) :: forward_plan,backward_plan
    type(C_PTR) :: c_w1,c_w2,c_rhol
    type(C_PTR) :: c_w1_filted,c_w2_filted,c_rho_filted
    type(C_PTR) :: c_w1w1,c_w1w2,c_w2w1,c_w2w2
    type(C_PTR) :: c_w1w1_filted,c_w1w2_filted,c_w2w1_filted,c_w2w2_filted
    type(C_PTR) :: c_A11,c_A12,c_A21,c_A22
    type(C_PTR) :: c_tau11_term,c_tau12_term,c_tau21_term,c_tau22_term
    !
    integer,dimension(8) :: value
    character(len=1) :: modeio
    logical :: loutput
    !
    call readinput
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
    !
    allocate(erroravg(1:3,1:3),errormax(1:3,1:3),errorgtr10(1:3,1:3),errorgtr100(1:3,1:3))
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !!!! Read velocity and density field
    allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km))
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    !! wavenumber
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    do j = 1,jm
    do i = 1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if((j+j0) <= (ja/2+1)) then
        k2(i,j) = real(j+j0-1,8)
      else if((j+j0)<=(ja)) then
        k2(i,j) = real(j+j0-ja-1,8)
      else
        print *,"Error, no wave number possible, (j+j0) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    !
    allocate(tau(1:2,1:2,1:im,1:jm),tau_bis(1:2,1:2,1:im,1:jm))
    !
    !! Imaginary number prepare
    imag = CMPLX(0.d0,1.d0,8)
    !
    if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
    !!!! Prepare l,alpha and others
    call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
    l_min = 2*pi/im
    allocate(l_lim(1:num_l),num_alphas(1:num_l),l_sqrtalpha(1:num_l,1:num_alpha))
    allocate(l_phi(1:num_l,1:num_alpha),dl_alpha(1:num_l,1:num_alpha))
    !
    do i=1,num_l
      l_lim(i) = exp(log(ratio_min)+(i-1) * (log(ratio_max)-log(ratio_min)) / (num_l-1)) * l_min
      num_alphas(i) = num_alphamin + int(sqrt(l_lim(i)/ratio_max/l_min)*real(num_alpha-num_alphamin))
      !
      do j=1,num_alphas(i)
        l_sqrtalpha(i,j) = sqrt( exp(log(0.2**2) + &
            (log((l_lim(i)/l_min)**2) - log(0.2**2))* (j-1) / (num_alphas(i)-1)) )*l_min
        l_phi(i,j) = sqrt( abs(l_lim(i)**2 - l_sqrtalpha(i,j)**2) )
      enddo
    enddo
    !
    do i=1,num_l
      dl_alpha(i,1) = l_sqrtalpha(i,1)**2 
      !
      do j=2,num_alphas(i)
        dl_alpha(i,j) = l_sqrtalpha(i,j)**2 -l_sqrtalpha(i,j-1)**2 
      enddo
    enddo
    !
    if(mpirank==0)  print *, "Integrate point allocated"
    !
    if(mpirank==0) then
      open(fh,file='pp/SGSintegral.info',form='formatted')
      write(fh,"(2(A9,1x))")'NumL','NumAlpha'
      write(fh,"(2(I9,1x))")num_l,num_alpha
      write(fh,"(2(A9,1x),2(A15,1x))")'i','j','l_lim','l_sqrtalpha'
      do i=1,num_l
        do j=1,num_alphas(i)
        ! Output file of rank information.
          write(fh,"(2(I9,1x),2(E15.7E3,1x))")i,j,l_lim(i),l_sqrtalpha(i,j)
        enddo
      enddo
      !
      close(fh)
      print*,' << SGSintegral.info ... done !'
    endif
    !
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    c_w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1, w1, [imfftw,jmfftw])
    c_w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2, w2, [imfftw,jmfftw])
    c_rhol = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rhol, rhol, [imfftw,jmfftw])
    !
    c_w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1_filted, w1_filted, [imfftw,jmfftw])
    c_w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2_filted, w2_filted, [imfftw,jmfftw])
    c_rho_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rho_filted, rho_filted, [imfftw,jmfftw])
    !
    c_w1w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w1, w1w1, [imfftw,jmfftw])
    c_w1w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w2, w1w2, [imfftw,jmfftw])
    c_w2w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w1, w2w1, [imfftw,jmfftw])
    c_w2w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w2, w2w2, [imfftw,jmfftw])
    !
    c_w1w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w1_filted, w1w1_filted, [imfftw,jmfftw])
    c_w1w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w2_filted, w1w2_filted, [imfftw,jmfftw])
    c_w2w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w1_filted, w2w1_filted, [imfftw,jmfftw])
    c_w2w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w2_filted, w2w2_filted, [imfftw,jmfftw])
    !
    c_A11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A11, A11, [imfftw,jmfftw])
    c_A12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A12, A12, [imfftw,jmfftw])
    c_A21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A21, A21, [imfftw,jmfftw])
    c_A22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A22, A22, [imfftw,jmfftw])
    !
    c_tau11_term = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_tau11_term, tau11_term, [imfftw,jmfftw])
    c_tau12_term = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_tau12_term, tau12_term, [imfftw,jmfftw])
    c_tau21_term = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_tau21_term, tau21_term, [imfftw,jmfftw])
    c_tau22_term = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_tau22_term, tau22_term, [imfftw,jmfftw])
    !
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    ! wi,wiwj,rhol in physical space, not filted
    do j=1,jm
    do i=1,im
      !
      w1(i,j)=CMPLX(vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
      w2(i,j)=CMPLX(vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
      w1w1(i,j)=CMPLX(vel(i,j,0,1)*vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
      w1w2(i,j)=CMPLX(vel(i,j,0,1)*vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
      w2w1(i,j)=CMPLX(vel(i,j,0,2)*vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
      w2w2(i,j)=CMPLX(vel(i,j,0,2)*vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
      rhol(i,j)=CMPLX(rho(i,j,0),0.d0,C_INTPTR_T);
      !
    end do
    end do
    !
    ! wi=rho ui,wiwj =rho ui uj,rhol in spectral space, not filted
    call fftw_mpi_execute_dft(forward_plan,w1,w1)
    call fftw_mpi_execute_dft(forward_plan,w2,w2)
    call fftw_mpi_execute_dft(forward_plan,w1w1,w1w1)
    call fftw_mpi_execute_dft(forward_plan,w1w2,w1w2)
    call fftw_mpi_execute_dft(forward_plan,w2w1,w2w1)
    call fftw_mpi_execute_dft(forward_plan,w2w2,w2w2)
    call fftw_mpi_execute_dft(forward_plan,rhol,rhol)
    !
    do j=1,jm
    do i=1,im
      !
      w1(i,j)=w1(i,j)/(1.d0*ia*ja)
      w2(i,j)=w2(i,j)/(1.d0*ia*ja)
      !
      w1w1(i,j)=w1w1(i,j)/(1.d0*ia*ja)
      w1w2(i,j)=w1w2(i,j)/(1.d0*ia*ja)
      w2w1(i,j)=w2w1(i,j)/(1.d0*ia*ja)
      w2w2(i,j)=w2w2(i,j)/(1.d0*ia*ja)
      !
      rhol(i,j)=rhol(i,j)/(1.d0*ia*ja)
      !
    end do
    end do
    !
    do m=1,num_l
      !
      if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
      !
      ! Method 1: tauij = rho uiuj - rho ui uj
      ! Filter scale: l
      !
      ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
      ! wi_filted=rho ui,wiwj=rho ui uj,rho_filted in spectral space, filted by l
      do j=1,jm
      do i=1,im
        Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) ! 
        !
        w1_filted(i,j)    = w1(i,j)   *Gl
        w2_filted(i,j)    = w2(i,j)   *Gl
        !
        w1w1_filted(i,j)  = w1w1(i,j) *Gl
        w1w2_filted(i,j)  = w1w2(i,j) *Gl
        w2w1_filted(i,j)  = w2w1(i,j) *Gl
        w2w2_filted(i,j)  = w2w2(i,j) *Gl
        !
        rho_filted(i,j)   = rhol(i,j) *Gl
      enddo
      enddo
      !
      ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
      ! wi_filted=(rho ui)_filted,wiwj=(rho ui uj)_filted,rho_filted in physical space, filted by l 
      call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
      call fftw_mpi_execute_dft(backward_plan,w1w1_filted,w1w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w1w2_filted,w1w2_filted)
      call fftw_mpi_execute_dft(backward_plan,w2w1_filted,w2w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w2w2_filted,w2w2_filted)
      call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
      !
      ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
      ! wi_filted=(ui)~filted,wiwj=(rho ui uj)_filted,rho_filted in physical space, filted by l 
      do j=1,jm
      do i=1,im
        !
        w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
        w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
        !
        tau(1,1,i,j) = dreal(w1w1_filted(i,j) - rho_filted(i,j) * w1_filted(i,j) * w1_filted(i,j))
        tau(1,2,i,j) = dreal(w1w2_filted(i,j) - rho_filted(i,j) * w1_filted(i,j) * w2_filted(i,j))
        tau(2,1,i,j) = dreal(w2w1_filted(i,j) - rho_filted(i,j) * w2_filted(i,j) * w1_filted(i,j))
        tau(2,2,i,j) = dreal(w2w2_filted(i,j) - rho_filted(i,j) * w2_filted(i,j) * w2_filted(i,j))
        !
        do p=1,2
        do q=1,2
          tau_bis(p,q,i,j) = 0.d0
        enddo
        enddo
        !
      enddo
      enddo
      !
      ! Method 2: tauij = int_0^l2 （rho_√α  Aik_√α  Ajk_√α_√l2-α 
      ! Filter scale: l
      !
      do n=1,num_alphas(m)
        !
        call date_and_time(values=value) 
        !
        if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alphas(m),',now is ',&
                                value(5), ':', value(6),':',value(7)
        !
        ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
        ! wi_filted=rho ui,rho_filted in spectral space, filted by sqrtalpha
        do j=1,jm
        do i=1,im
          Galpha = exp(-(k1(i,j)**2+k2(i,j)**2)*l_sqrtalpha(m,n)**2/2.d0) ! Filtre scale :sqrtalpha
          w1_filted(i,j)  = w1(i,j)  *Galpha
          w2_filted(i,j)  = w2(i,j)  *Galpha
          rho_filted(i,j) = rhol(i,j)*Galpha
        enddo
        enddo
        !
        ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
        ! wi_filted=rho ui,rho_filted in physical space, filted by sqrtalpha
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! wi_filted=(ui)~filted, rho_filted in physical space, filted by sqrtalpha 
        do j=1,jm
        do i=1,im
          w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
          w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
        enddo
        enddo
        !
        ! wi_filted=(ui)~filted, Aij = Aij~filted in spectral space, filted by sqrtalpha 
        ! rho_filted in physical space, filted by sqrtalpha 
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j)  = w1_filted(i,j)/(1.d0*ia*ja)
          w2_filted(i,j)  = w2_filted(i,j)/(1.d0*ia*ja)
          !
          A11(i,j) = imag*w1_filted(i,j)*k1(i,j)
          A21(i,j) = imag*w2_filted(i,j)*k1(i,j)
          A12(i,j) = imag*w1_filted(i,j)*k2(i,j)
          A22(i,j) = imag*w2_filted(i,j)*k2(i,j)
          !
        end do
        end do
        !
        ! wi_filted=(ui)~filted in spectral space, filted by sqrtalpha 
        ! rho_filted, Aij = Aij~filted in physical space, filted by sqrtalpha 
        call fftw_mpi_execute_dft(backward_plan,A11,A11)
        call fftw_mpi_execute_dft(backward_plan,A21,A21)
        call fftw_mpi_execute_dft(backward_plan,A12,A12)
        call fftw_mpi_execute_dft(backward_plan,A22,A22)
        !
        do j=1,jm
        do i=1,im
          !
          tau11_term(i,j) = rho_filted(i,j) * (A11(i,j)*A11(i,j)+A12(i,j)*A12(i,j))
          tau12_term(i,j) = rho_filted(i,j) * (A11(i,j)*A21(i,j)+A12(i,j)*A22(i,j))
          tau21_term(i,j) = rho_filted(i,j) * (A21(i,j)*A11(i,j)+A22(i,j)*A12(i,j))
          tau22_term(i,j) = rho_filted(i,j) * (A21(i,j)*A21(i,j)+A22(i,j)*A22(i,j))
          !
        end do
        end do
        !
        ! Do filter phi:
        ! F -> product -> F inverse
        call fftw_mpi_execute_dft(forward_plan,tau11_term,tau11_term)
        call fftw_mpi_execute_dft(forward_plan,tau12_term,tau12_term)
        call fftw_mpi_execute_dft(forward_plan,tau21_term,tau21_term)
        call fftw_mpi_execute_dft(forward_plan,tau22_term,tau22_term)
        !
        do j=1,jm
        do i=1,im
          Gphi = exp(-(k1(i,j)**2+k2(i,j)**2)*l_phi(m,n)**2/2.d0) ! Filtre scale :phi
          tau11_term(i,j) = tau11_term(i,j)*Gphi/(1.d0*ia*ja)
          tau12_term(i,j) = tau12_term(i,j)*Gphi/(1.d0*ia*ja)
          tau21_term(i,j) = tau21_term(i,j)*Gphi/(1.d0*ia*ja)
          tau22_term(i,j) = tau22_term(i,j)*Gphi/(1.d0*ia*ja)
          !
        enddo
        enddo
        !
        !
        call fftw_mpi_execute_dft(backward_plan,tau11_term,tau11_term)
        call fftw_mpi_execute_dft(backward_plan,tau12_term,tau12_term)
        call fftw_mpi_execute_dft(backward_plan,tau21_term,tau21_term)
        call fftw_mpi_execute_dft(backward_plan,tau22_term,tau22_term)
        !
        !
        do j=1,jm
        do i=1,im
          !
          tau_bis(1,1,i,j) = tau_bis(1,1,i,j) + dreal(tau11_term(i,j)) * dl_alpha(m,n)
          tau_bis(1,2,i,j) = tau_bis(1,2,i,j) + dreal(tau12_term(i,j)) * dl_alpha(m,n)
          tau_bis(2,1,i,j) = tau_bis(2,1,i,j) + dreal(tau21_term(i,j)) * dl_alpha(m,n)
          tau_bis(2,2,i,j) = tau_bis(2,2,i,j) + dreal(tau22_term(i,j)) * dl_alpha(m,n)
          !
        enddo
        enddo
        !
      enddo ! loop of integral (alpha)
      !
      do p=1,2
      do q=1,2
        erroravg(p,q)=0.d0
        errormax(p,q)=0.d0
        errorgtr10(p,q)=0.d0
        errorgtr100(p,q)=0.d0
        errornorm2max=0.d0
        errornorm2avg=0.d0
        errornorm2gtr10=0.d0
        errornorm2gtr100=0.d0
      enddo
      enddo
      !
      ! Output comparaison results
      do j=1,jm
      do i=1,im
        !
        norm2 = 0.d0
        norm2bis = 0.d0
        do p=1,2
        do q=1,2
          if(abs(tau(p,q,i,j))>1.d-6)then
            result = abs(tau_bis(p,q,i,j)-tau(p,q,i,j))/(abs(tau_bis(p,q,i,j))+abs(tau(p,q,i,j)))
            norm2 = norm2 + tau(p,q,i,j)**2
            norm2bis = norm2bis + tau_bis(p,q,i,j)**2
            errormax(p,q) = max(errormax(p,q),result)
            erroravg(p,q) = erroravg(p,q) + result
            if(result > 0.1)then
              errorgtr10(p,q) = errorgtr10(p,q) + 1.d0
            endif
            if(result > 1)then
              errorgtr100(p,q) = errorgtr100(p,q) + 1.d0
            endif
          endif
        enddo
        enddo
        result = abs(norm2bis-norm2)/(abs(norm2)+abs(norm2bis))
        errornorm2max = max(errornorm2max,result)
        errornorm2avg = errornorm2avg + result
        if(result > 0.1)then
          errornorm2gtr10 = errornorm2gtr10 + 1.d0
        endif
        if(result > 1)then
          errornorm2gtr100 = errornorm2gtr100 + 1.d0
        endif
        !
      enddo
      enddo
      !
      !
      do p=1,2
      do q=1,2
        errormax(p,q) = pmax(errormax(p,q))
        erroravg(p,q) = psum(erroravg(p,q))/(1.d0*ia*ja)
        errorgtr10(p,q) = psum(errorgtr10(p,q))/(1.d0*ia*ja)
        errorgtr100(p,q) = psum(errorgtr100(p,q))/(1.d0*ia*ja)
        errornorm2max = psum(errornorm2max)
        errornorm2avg = psum(errornorm2avg)/(1.d0*ia*ja)
        errornorm2gtr10 = psum(errornorm2gtr10)/(1.d0*ia*ja)
        errornorm2gtr100 = psum(errornorm2gtr100)/(1.d0*ia*ja)
      enddo
      enddo
      !
      if(mpirank==0) then
        write(mname,'(i4.4)')m
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_stress_relative_error_'//stepname//'_'//mname//'.dat'
        else
          outfilename = 'pp/SGS_stress_relative_error_'//mname//'.dat'
        endif
        !
        if(mpirank == 0)then
          open(fh,file=outfilename,form='formatted')
          write(fh,"(A7,1x,2(A20,1x))")'nstep','time','l'
          write(fh,"(I7,1x,2(E20.13E2,1x))")nstep,time,l_lim(m)
          write(fh,"(A8,1x,5(A20,1x))")'type','tau11','tau12','tau21','tau22','norm2'
          write(fh,"(A8,1x,5(E20.13E2,1x))")'max',((errormax(p,q),q=1,2),p=1,2),errornorm2max
          write(fh,"(A8,1x,5(E20.13E2,1x))")'avg',((erroravg(p,q),q=1,2),p=1,2),errornorm2avg
          write(fh,"(A8,1x,5(E20.13E2,1x))")'gtr0.1',((errorgtr10(p,q),q=1,2),p=1,2),errornorm2gtr10
          write(fh,"(A8,1x,5(E20.13E2,1x))")'gtr1',((errorgtr100(p,q),q=1,2),p=1,2),errornorm2gtr100
          close(fh)
          print *, '>>>>', outfilename
        endif
        !
        !
        !
      endif
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(loutput)then
        !
        write(mname,'(i4.4)')m
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_stress_'//stepname//'_'//mname//'.h5'
        else
          outfilename = 'pp/SGS_stress_'//mname//'.h5'
        endif
        !
        call h5io_init(trim(outfilename),mode='write')
        !
        do p=1,2
        do q= 1,2
          write (termname, "(A3,I1,I1)") "tau",p,q
          call h5wa2d_r8(varname=termname,var=tau(p,q,1:im,1:jm),    dir='k')
          write (termname, "(A3,I1,I1,A3)") "tau",p,q,"bis"
          call h5wa2d_r8(varname=termname,var=tau_bis(p,q,1:im,1:jm),dir='k')
        enddo
        enddo
        !
        call h5io_end
        !
      endif
      !
    enddo ! loop of filter point l
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_w1)
    call fftw_free(c_w2)
    call fftw_free(c_rhol)
    call fftw_free(c_w1_filted)
    call fftw_free(c_w2_filted)
    call fftw_free(c_rho_filted)
    call fftw_free(c_w1w1)
    call fftw_free(c_w1w2)
    call fftw_free(c_w2w1)
    call fftw_free(c_w2w2)
    call fftw_free(c_w1w1_filted)
    call fftw_free(c_w1w2_filted)
    call fftw_free(c_w2w1_filted)
    call fftw_free(c_w2w2_filted)
    call fftw_free(c_A11)
    call fftw_free(c_A12)
    call fftw_free(c_A21)
    call fftw_free(c_A22)
    call fftw_free(c_tau11_term)
    call fftw_free(c_tau12_term)
    call fftw_free(c_tau21_term)
    call fftw_free(c_tau22_term)
    call mpistop
    deallocate(tau,tau_bis,erroravg,errormax,errorgtr10,errorgtr100)
    !
  end subroutine SGSstress2D
  !
  subroutine SGSstress3D(thefilenumb)
    ! 
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km
    use commarray, only: vel, rho
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
    include 'fftw3-mpi.f03'
    !
    integer,intent(in) :: thefilenumb
    integer :: fh
    integer :: i,j,k,m,n,p,q
    character(len=128) :: infilename,outfilename
    character(len=4) :: stepname,mname
    character(len=10):: termname
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
    complex(8) :: imag
    real(8),allocatable,dimension(:) :: l_lim
    real(8),allocatable,dimension(:,:) :: l_sqrtalpha,l_phi,dl_alpha
    integer,allocatable,dimension(:) :: num_alphas
    integer :: num_l,num_alpha,num_alphamin
    integer :: hand_a
    real(8) :: l_min, ratio_max, ratio_min
    real(8) :: Gl,Galpha,Gphi
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhol
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1w1,w1w2,w1w3
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w2w1,w2w2,w2w3
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w3w1,w3w2,w3w3
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1w1_filted,w1w2_filted,w1w3_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w2w1_filted,w2w2_filted,w2w3_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w3w1_filted,w3w2_filted,w3w3_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11,A12,A13
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21,A22,A23
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31,A32,A33
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: tau11_term,tau12_term,tau13_term
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: tau21_term,tau22_term,tau23_term
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: tau31_term,tau32_term,tau33_term
    real(8), allocatable, dimension(:,:,:,:,:) :: tau ! 1:3,1:3,1:im,1:jm,1:km
    real(8), allocatable, dimension(:,:,:,:,:) :: tau_bis ! 1:3,1:3,1:im,1:jm,1:km
    real(8), allocatable, dimension(:,:) :: errormax,erroravg, errorgtr10,errorgtr100
    real(8) :: result,norm2,norm2bis
    real(8) :: errornorm2max,errornorm2avg,errornorm2gtr10,errornorm2gtr100
    !
    !
    type(C_PTR) :: forward_plan,backward_plan
    type(C_PTR) :: c_w1,c_w2,c_w3,c_rhol
    type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted
    type(C_PTR) :: c_w1w1,c_w1w2,c_w1w3
    type(C_PTR) :: c_w2w1,c_w2w2,c_w2w3
    type(C_PTR) :: c_w3w1,c_w3w2,c_w3w3
    type(C_PTR) :: c_w1w1_filted,c_w1w2_filted,c_w1w3_filted
    type(C_PTR) :: c_w2w1_filted,c_w2w2_filted,c_w2w3_filted
    type(C_PTR) :: c_w3w1_filted,c_w3w2_filted,c_w3w3_filted
    type(C_PTR) :: c_A11,c_A12,c_A13
    type(C_PTR) :: c_A21,c_A22,c_A23
    type(C_PTR) :: c_A31,c_A32,c_A33
    type(C_PTR) :: c_tau11_term,c_tau12_term,c_tau13_term,c_tau21_term,c_tau22_term
    type(C_PTR) :: c_tau23_term,c_tau31_term,c_tau32_term,c_tau33_term
    !
    integer,dimension(8) :: value
    character(len=1) :: modeio
    logical :: loutput
    !
    call readinput
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
    !
    allocate(erroravg(1:3,1:3),errormax(1:3,1:3),errorgtr10(1:3,1:3),errorgtr100(1:3,1:3))
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !!!! Read velocity and density field
    allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    !! wavenumber
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    do k = 1,km
    do j = 1,jm
    do i = 1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j,k) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j,k) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if(j <= (ja/2+1)) then
        k2(i,j,k) = real(j-1,8)
      else if(i<=(ia)) then
        k2(i,j,k) = real(j-ja-1,8)
      else
        print *,"Error, no wave number possible, j must smaller than ja-1 !"
      end if
      !
      if((k+k0) <= (ka/2+1)) then
        k3(i,j,k) = real(k+k0-1,8)
      else if((k+k0)<=(ka)) then
        k3(i,j,k) = real(k+k0-ka-1,8)
      else
        print *,"Error, no wave number possible, (k+k0) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    end do
    !
    allocate(tau(1:3,1:3,1:im,1:jm,1:km),tau_bis(1:3,1:3,1:im,1:jm,1:km))
    !
    !! Imaginary number prepare
    imag = CMPLX(0.d0,1.d0,8)
    !
    if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
    !!!! Prepare l,alpha and others
    call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
    l_min = 2*pi/im
    allocate(l_lim(1:num_l),num_alphas(1:num_l),l_sqrtalpha(1:num_l,1:num_alpha))
    allocate(l_phi(1:num_l,1:num_alpha),dl_alpha(1:num_l,1:num_alpha))
    !
    do i=1,num_l
      l_lim(i) = exp(log(ratio_min)+(i-1) * (log(ratio_max)-log(ratio_min)) / (num_l-1)) * l_min
      num_alphas(i) = num_alphamin + int(sqrt(l_lim(i)/ratio_max/l_min)*real(num_alpha-num_alphamin))
      !
      do j=1,num_alphas(i)
        l_sqrtalpha(i,j) = sqrt( exp(log(0.2**2) + &
            (log((l_lim(i)/l_min)**2) - log(0.2**2))* (j-1) / (num_alphas(i)-1)) )*l_min
        l_phi(i,j) = sqrt( abs(l_lim(i)**2 - l_sqrtalpha(i,j)**2) )
      enddo
    enddo
    !
    do i=1,num_l
      dl_alpha(i,1) = l_sqrtalpha(i,1)**2 
      !
      do j=2,num_alphas(i)
        dl_alpha(i,j) = l_sqrtalpha(i,j)**2 -l_sqrtalpha(i,j-1)**2 
      enddo
    enddo
    !
    if(mpirank==0)  print *, "Integrate point allocated"
    !
    if(mpirank==0) then
      open(fh,file='pp/SGSintegral.info',form='formatted')
      write(fh,"(2(A9,1x))")'NumL','NumAlpha'
      write(fh,"(2(I9,1x))")num_l,num_alpha
      write(fh,"(2(A9,1x),2(A15,1x))")'i','j','l_lim','l_sqrtalpha'
      do i=1,num_l
        do j=1,num_alphas(i)
        ! Output file of rank information.
          write(fh,"(2(I9,1x),2(E15.7E3,1x))")i,j,l_lim(i),l_sqrtalpha(i,j)
        enddo
      enddo
      !
      close(fh)
      print*,' << SGSintegral.info ... done !'
    endif
    !
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    c_w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
    c_w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
    c_w3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
    c_rhol = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rhol, rhol, [imfftw,jmfftw,kmfftw])
    !
    c_w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1_filted, w1_filted, [imfftw,jmfftw,kmfftw])
    c_w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2_filted, w2_filted, [imfftw,jmfftw,kmfftw])
    c_w3_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3_filted, w3_filted, [imfftw,jmfftw,kmfftw])
    c_rho_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rho_filted, rho_filted, [imfftw,jmfftw,kmfftw])
    !
    c_w1w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w1, w1w1, [imfftw,jmfftw,kmfftw])
    c_w1w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w2, w1w2, [imfftw,jmfftw,kmfftw])
    c_w1w3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w3, w1w3, [imfftw,jmfftw,kmfftw])
    c_w2w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w1, w2w1, [imfftw,jmfftw,kmfftw])
    c_w2w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w2, w2w2, [imfftw,jmfftw,kmfftw])
    c_w2w3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w3, w2w3, [imfftw,jmfftw,kmfftw])
    c_w3w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3w1, w3w1, [imfftw,jmfftw,kmfftw])
    c_w3w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3w2, w3w2, [imfftw,jmfftw,kmfftw])
    c_w3w3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3w3, w3w3, [imfftw,jmfftw,kmfftw])
    !
    c_w1w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w1_filted, w1w1_filted, [imfftw,jmfftw,kmfftw])
    c_w1w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w2_filted, w1w2_filted, [imfftw,jmfftw,kmfftw])
    c_w1w3_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1w3_filted, w1w3_filted, [imfftw,jmfftw,kmfftw])
    c_w2w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w1_filted, w2w1_filted, [imfftw,jmfftw,kmfftw])
    c_w2w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w2_filted, w2w2_filted, [imfftw,jmfftw,kmfftw])
    c_w2w3_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2w3_filted, w2w3_filted, [imfftw,jmfftw,kmfftw])
    c_w3w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3w1_filted, w3w1_filted, [imfftw,jmfftw,kmfftw])
    c_w3w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3w2_filted, w3w2_filted, [imfftw,jmfftw,kmfftw])
    c_w3w3_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3w3_filted, w3w3_filted, [imfftw,jmfftw,kmfftw])
    !
    c_A11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A11, A11, [imfftw,jmfftw,kmfftw])
    c_A12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A12, A12, [imfftw,jmfftw,kmfftw])
    c_A13 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A13, A13, [imfftw,jmfftw,kmfftw])
    c_A21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A21, A21, [imfftw,jmfftw,kmfftw])
    c_A22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A22, A22, [imfftw,jmfftw,kmfftw])
    c_A23 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A23, A23, [imfftw,jmfftw,kmfftw])
    c_A31 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A31, A31, [imfftw,jmfftw,kmfftw])
    c_A32 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A32, A32, [imfftw,jmfftw,kmfftw])
    c_A33 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A33, A33, [imfftw,jmfftw,kmfftw])
    !
    c_tau11_term = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_tau11_term, tau11_term, [imfftw,jmfftw,kmfftw])
    c_tau12_term = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_tau12_term, tau12_term, [imfftw,jmfftw,kmfftw])
    c_tau13_term = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_tau13_term, tau13_term, [imfftw,jmfftw,kmfftw])
    c_tau21_term = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_tau21_term, tau21_term, [imfftw,jmfftw,kmfftw])
    c_tau22_term = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_tau22_term, tau22_term, [imfftw,jmfftw,kmfftw])
    c_tau23_term = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_tau23_term, tau23_term, [imfftw,jmfftw,kmfftw])
    c_tau31_term = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_tau31_term, tau31_term, [imfftw,jmfftw,kmfftw])
    c_tau32_term = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_tau32_term, tau32_term, [imfftw,jmfftw,kmfftw])
    c_tau33_term = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_tau33_term, tau33_term, [imfftw,jmfftw,kmfftw])
    !
    forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    ! wi,wiwj,rhol in physical space, not filted
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
      w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
      w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
      w1w1(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
      w1w2(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
      w1w3(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
      w2w1(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
      w2w2(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
      w2w3(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
      w3w1(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
      w3w2(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
      w3w3(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
      rhol(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T);
      !
    end do
    end do
    end do
    !
    ! wi=rho ui,wiwj =rho ui uj,rhol in spectral space, not filted
    call fftw_mpi_execute_dft(forward_plan,w1,w1)
    call fftw_mpi_execute_dft(forward_plan,w2,w2)
    call fftw_mpi_execute_dft(forward_plan,w3,w3)
    call fftw_mpi_execute_dft(forward_plan,w1w1,w1w1)
    call fftw_mpi_execute_dft(forward_plan,w1w2,w1w2)
    call fftw_mpi_execute_dft(forward_plan,w1w3,w1w3)
    call fftw_mpi_execute_dft(forward_plan,w2w1,w2w1)
    call fftw_mpi_execute_dft(forward_plan,w2w2,w2w2)
    call fftw_mpi_execute_dft(forward_plan,w2w3,w2w3)
    call fftw_mpi_execute_dft(forward_plan,w3w1,w3w1)
    call fftw_mpi_execute_dft(forward_plan,w3w2,w3w2)
    call fftw_mpi_execute_dft(forward_plan,w3w3,w3w3)
    call fftw_mpi_execute_dft(forward_plan,rhol,rhol)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
      w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
      w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
      !
      w1w1(i,j,k)=w1w1(i,j,k)/(1.d0*ia*ja*ka)
      w1w2(i,j,k)=w1w2(i,j,k)/(1.d0*ia*ja*ka)
      w1w3(i,j,k)=w1w3(i,j,k)/(1.d0*ia*ja*ka)
      w2w1(i,j,k)=w2w1(i,j,k)/(1.d0*ia*ja*ka)
      w2w2(i,j,k)=w2w2(i,j,k)/(1.d0*ia*ja*ka)
      w2w3(i,j,k)=w2w3(i,j,k)/(1.d0*ia*ja*ka)
      w3w1(i,j,k)=w3w1(i,j,k)/(1.d0*ia*ja*ka)
      w3w2(i,j,k)=w3w2(i,j,k)/(1.d0*ia*ja*ka)
      w3w3(i,j,k)=w3w3(i,j,k)/(1.d0*ia*ja*ka)
      !
      rhol(i,j,k)=rhol(i,j,k)/(1.d0*ia*ja*ka)
      !
    end do
    end do
    end do
    !
    do m=1,num_l
      !
      if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
      !
      ! Method 1: tauij = rho uiuj - rho ui uj
      ! Filter scale: l
      !
      ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
      ! wi_filted=rho ui,wiwj=rho ui uj,rho_filted in spectral space, filted by l
      do k=1,km
      do j=1,jm
      do i=1,im
        Gl = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_lim(m)**2/2.d0)
        !
        w1_filted(i,j,k)    = w1(i,j,k)   *Gl
        w2_filted(i,j,k)    = w2(i,j,k)   *Gl
        w3_filted(i,j,k)    = w3(i,j,k)   *Gl
        !
        w1w1_filted(i,j,k)  = w1w1(i,j,k) *Gl
        w1w2_filted(i,j,k)  = w1w2(i,j,k) *Gl
        w1w3_filted(i,j,k)  = w1w3(i,j,k) *Gl
        w2w1_filted(i,j,k)  = w2w1(i,j,k) *Gl
        w2w2_filted(i,j,k)  = w2w2(i,j,k) *Gl
        w2w3_filted(i,j,k)  = w2w3(i,j,k) *Gl
        w3w1_filted(i,j,k)  = w3w1(i,j,k) *Gl
        w3w2_filted(i,j,k)  = w3w2(i,j,k) *Gl
        w3w3_filted(i,j,k)  = w3w3(i,j,k) *Gl
        !
        rho_filted(i,j,k)   = rhol(i,j,k) *Gl
      enddo
      enddo
      enddo
      !
      ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
      ! wi_filted=(rho ui)_filted,wiwj=(rho ui uj)_filted,rho_filted in physical space, filted by l 
      call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
      call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
      call fftw_mpi_execute_dft(backward_plan,w1w1_filted,w1w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w1w2_filted,w1w2_filted)
      call fftw_mpi_execute_dft(backward_plan,w1w3_filted,w1w3_filted)
      call fftw_mpi_execute_dft(backward_plan,w2w1_filted,w2w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w2w2_filted,w2w2_filted)
      call fftw_mpi_execute_dft(backward_plan,w2w3_filted,w2w3_filted)
      call fftw_mpi_execute_dft(backward_plan,w3w1_filted,w3w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w3w2_filted,w3w2_filted)
      call fftw_mpi_execute_dft(backward_plan,w3w3_filted,w3w3_filted)
      call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
      !
      ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
      ! wi_filted=(ui)~filted,wiwj=(rho ui uj)_filted,rho_filted in physical space, filted by l 
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
        w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
        w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
        !
        tau(1,1,i,j,k) = dreal(w1w1_filted(i,j,k) - rho_filted(i,j,k) * w1_filted(i,j,k) * w1_filted(i,j,k))
        tau(1,2,i,j,k) = dreal(w1w2_filted(i,j,k) - rho_filted(i,j,k) * w1_filted(i,j,k) * w2_filted(i,j,k))
        tau(1,3,i,j,k) = dreal(w1w3_filted(i,j,k) - rho_filted(i,j,k) * w1_filted(i,j,k) * w3_filted(i,j,k))
        tau(2,1,i,j,k) = dreal(w2w1_filted(i,j,k) - rho_filted(i,j,k) * w2_filted(i,j,k) * w1_filted(i,j,k))
        tau(2,2,i,j,k) = dreal(w2w2_filted(i,j,k) - rho_filted(i,j,k) * w2_filted(i,j,k) * w2_filted(i,j,k))
        tau(2,3,i,j,k) = dreal(w2w3_filted(i,j,k) - rho_filted(i,j,k) * w2_filted(i,j,k) * w3_filted(i,j,k))
        tau(3,1,i,j,k) = dreal(w3w1_filted(i,j,k) - rho_filted(i,j,k) * w3_filted(i,j,k) * w1_filted(i,j,k))
        tau(3,2,i,j,k) = dreal(w3w2_filted(i,j,k) - rho_filted(i,j,k) * w3_filted(i,j,k) * w2_filted(i,j,k))
        tau(3,3,i,j,k) = dreal(w3w3_filted(i,j,k) - rho_filted(i,j,k) * w3_filted(i,j,k) * w3_filted(i,j,k))
        !
        do p=1,3
        do q=1,3
          tau_bis(p,q,i,j,k) = 0.d0
        enddo
        enddo
        !
      enddo
      enddo
      enddo
      !
      ! Method 2: tauij = int_0^l2 （rho_√α  Aik_√α  Ajk_√α_√l2-α 
      ! Filter scale: l
      !
      do n=1,num_alphas(m)
        !
        call date_and_time(values=value) 
        !
        if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alphas(m),',now is ',&
                                value(5), ':', value(6),':',value(7)
        !
        ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
        ! wi_filted=rho ui,rho_filted in spectral space, filted by sqrtalpha
        do k=1,km
        do j=1,jm
        do i=1,im
          Galpha = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_sqrtalpha(m,n)**2/2.d0)
          w1_filted(i,j,k)  = w1(i,j,k)  *Galpha
          w2_filted(i,j,k)  = w2(i,j,k)  *Galpha
          w3_filted(i,j,k)  = w3(i,j,k)  *Galpha
          rho_filted(i,j,k) = rhol(i,j,k)*Galpha
        enddo
        enddo
        enddo
        !
        ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
        ! wi_filted=rho ui,rho_filted in physical space, filted by sqrtalpha
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! wi_filted=(ui)~filted, rho_filted in physical space, filted by sqrtalpha 
        do k=1,km
        do j=1,jm
        do i=1,im
          w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
          w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
          w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
        enddo
        enddo
        enddo
        !
        ! wi_filted=(ui)~filted, Aij = Aij~filted in spectral space, filted by sqrtalpha 
        ! rho_filted in physical space, filted by sqrtalpha 
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j,k)  = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
          w2_filted(i,j,k)  = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
          w3_filted(i,j,k)  = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
          !
          A11(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
          A21(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
          A31(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
          A12(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
          A22(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
          A32(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
          A13(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
          A23(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
          A33(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
          !
        end do
        end do
        end do
        !
        ! wi_filted=(ui)~filted in spectral space, filted by sqrtalpha 
        ! rho_filted, Aij = Aij~filted in physical space, filted by sqrtalpha 
        call fftw_mpi_execute_dft(backward_plan,A11,A11)
        call fftw_mpi_execute_dft(backward_plan,A21,A21)
        call fftw_mpi_execute_dft(backward_plan,A31,A31)
        call fftw_mpi_execute_dft(backward_plan,A12,A12)
        call fftw_mpi_execute_dft(backward_plan,A22,A22)
        call fftw_mpi_execute_dft(backward_plan,A32,A32)
        call fftw_mpi_execute_dft(backward_plan,A13,A13)
        call fftw_mpi_execute_dft(backward_plan,A23,A23)
        call fftw_mpi_execute_dft(backward_plan,A33,A33)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          tau11_term(i,j,k) = rho_filted(i,j,k) * (A11(i,j,k)*A11(i,j,k)+A12(i,j,k)*A12(i,j,k)+A13(i,j,k)*A13(i,j,k))
          tau12_term(i,j,k) = rho_filted(i,j,k) * (A11(i,j,k)*A21(i,j,k)+A12(i,j,k)*A22(i,j,k)+A13(i,j,k)*A23(i,j,k))
          tau13_term(i,j,k) = rho_filted(i,j,k) * (A11(i,j,k)*A31(i,j,k)+A12(i,j,k)*A32(i,j,k)+A13(i,j,k)*A33(i,j,k))
          tau21_term(i,j,k) = rho_filted(i,j,k) * (A21(i,j,k)*A11(i,j,k)+A22(i,j,k)*A12(i,j,k)+A23(i,j,k)*A13(i,j,k))
          tau22_term(i,j,k) = rho_filted(i,j,k) * (A21(i,j,k)*A21(i,j,k)+A22(i,j,k)*A22(i,j,k)+A23(i,j,k)*A23(i,j,k))
          tau23_term(i,j,k) = rho_filted(i,j,k) * (A21(i,j,k)*A31(i,j,k)+A22(i,j,k)*A32(i,j,k)+A23(i,j,k)*A33(i,j,k))
          tau31_term(i,j,k) = rho_filted(i,j,k) * (A31(i,j,k)*A11(i,j,k)+A32(i,j,k)*A12(i,j,k)+A33(i,j,k)*A13(i,j,k))
          tau32_term(i,j,k) = rho_filted(i,j,k) * (A31(i,j,k)*A21(i,j,k)+A32(i,j,k)*A22(i,j,k)+A33(i,j,k)*A23(i,j,k))
          tau33_term(i,j,k) = rho_filted(i,j,k) * (A31(i,j,k)*A31(i,j,k)+A32(i,j,k)*A32(i,j,k)+A33(i,j,k)*A33(i,j,k))
          !
        end do
        end do
        end do
        !
        ! Do filter phi:
        ! F -> product -> F inverse
        call fftw_mpi_execute_dft(forward_plan,tau11_term,tau11_term)
        call fftw_mpi_execute_dft(forward_plan,tau12_term,tau12_term)
        call fftw_mpi_execute_dft(forward_plan,tau13_term,tau13_term)
        call fftw_mpi_execute_dft(forward_plan,tau21_term,tau21_term)
        call fftw_mpi_execute_dft(forward_plan,tau22_term,tau22_term)
        call fftw_mpi_execute_dft(forward_plan,tau23_term,tau23_term)
        call fftw_mpi_execute_dft(forward_plan,tau31_term,tau31_term)
        call fftw_mpi_execute_dft(forward_plan,tau32_term,tau32_term)
        call fftw_mpi_execute_dft(forward_plan,tau33_term,tau33_term)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          Gphi = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_phi(m,n)**2/2.d0) ! Filtre scale :phi
          tau11_term(i,j,k) = tau11_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          tau12_term(i,j,k) = tau12_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          tau13_term(i,j,k) = tau13_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          tau21_term(i,j,k) = tau21_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          tau22_term(i,j,k) = tau22_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          tau23_term(i,j,k) = tau23_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          tau31_term(i,j,k) = tau31_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          tau32_term(i,j,k) = tau32_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          tau33_term(i,j,k) = tau33_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
          !
        enddo
        enddo
        enddo
        !
        !
        call fftw_mpi_execute_dft(backward_plan,tau11_term,tau11_term)
        call fftw_mpi_execute_dft(backward_plan,tau12_term,tau12_term)
        call fftw_mpi_execute_dft(backward_plan,tau13_term,tau13_term)
        call fftw_mpi_execute_dft(backward_plan,tau21_term,tau21_term)
        call fftw_mpi_execute_dft(backward_plan,tau22_term,tau22_term)
        call fftw_mpi_execute_dft(backward_plan,tau23_term,tau23_term)
        call fftw_mpi_execute_dft(backward_plan,tau31_term,tau31_term)
        call fftw_mpi_execute_dft(backward_plan,tau32_term,tau32_term)
        call fftw_mpi_execute_dft(backward_plan,tau33_term,tau33_term)
        !
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          tau_bis(1,1,i,j,k) = tau_bis(1,1,i,j,k) + dreal(tau11_term(i,j,k)) * dl_alpha(m,n)
          tau_bis(1,2,i,j,k) = tau_bis(1,2,i,j,k) + dreal(tau12_term(i,j,k)) * dl_alpha(m,n)
          tau_bis(1,3,i,j,k) = tau_bis(1,3,i,j,k) + dreal(tau13_term(i,j,k)) * dl_alpha(m,n)
          tau_bis(2,1,i,j,k) = tau_bis(2,1,i,j,k) + dreal(tau21_term(i,j,k)) * dl_alpha(m,n)
          tau_bis(2,2,i,j,k) = tau_bis(2,2,i,j,k) + dreal(tau22_term(i,j,k)) * dl_alpha(m,n)
          tau_bis(2,3,i,j,k) = tau_bis(2,3,i,j,k) + dreal(tau23_term(i,j,k)) * dl_alpha(m,n)
          tau_bis(3,1,i,j,k) = tau_bis(3,1,i,j,k) + dreal(tau31_term(i,j,k)) * dl_alpha(m,n)
          tau_bis(3,2,i,j,k) = tau_bis(3,2,i,j,k) + dreal(tau32_term(i,j,k)) * dl_alpha(m,n)
          tau_bis(3,3,i,j,k) = tau_bis(3,3,i,j,k) + dreal(tau33_term(i,j,k)) * dl_alpha(m,n)
          !
        enddo
        enddo
        enddo
        !
      enddo ! loop of integral (alpha)
      !
      do p=1,3
      do q=1,3
        erroravg(p,q)=0.d0
        errormax(p,q)=0.d0
        errorgtr10(p,q)=0.d0
        errorgtr100(p,q)=0.d0
      enddo
      enddo
      errornorm2max=0.d0
      errornorm2avg=0.d0
      errornorm2gtr10=0.d0
      errornorm2gtr100=0.d0
      !
      ! Output comparaison results
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        norm2 = 0.d0
        norm2bis = 0.d0
        do p=1,3
        do q=1,3
          if(abs(tau(p,q,i,j,k))>2.d-4)then
            result = 2*abs(tau_bis(p,q,i,j,k)-tau(p,q,i,j,k))/(abs(tau_bis(p,q,i,j,k))+abs(tau(p,q,i,j,k)))
            norm2 = norm2 + tau(p,q,i,j,k)**2
            norm2bis = norm2bis + tau_bis(p,q,i,j,k)**2
            errormax(p,q) = max(errormax(p,q),result)
            erroravg(p,q) = erroravg(p,q) + result
            if(result > 0.1)then
              errorgtr10(p,q) = errorgtr10(p,q) + 1.d0
            endif
            if(result > 1)then
              errorgtr100(p,q) = errorgtr100(p,q) + 1.d0
              print *, 'p',p,'q',q,tau_bis(p,q,i,j,k),tau(p,q,i,j,k),result, '*'
            endif
          endif
        enddo
        enddo
        result = 2*abs(norm2bis-norm2)/(abs(norm2)+abs(norm2bis))
        errornorm2max = max(errornorm2max,result)
        errornorm2avg = errornorm2avg + result
        if(result > 0.1)then
          errornorm2gtr10 = errornorm2gtr10 + 1.d0
        endif
        if(result > 1)then
          errornorm2gtr100 = errornorm2gtr100 + 1.d0
          print *, norm2, norm2bis, result, '*'
        endif
        !
      enddo
      enddo
      enddo
      !
      !
      do p=1,3
      do q=1,3
        errormax(p,q) = pmax(errormax(p,q))
        erroravg(p,q) = psum(erroravg(p,q))/(1.d0*ia*ja*ka)
        errorgtr10(p,q) = psum(errorgtr10(p,q))/(1.d0*ia*ja*ka)
        errorgtr100(p,q) = psum(errorgtr100(p,q))/(1.d0*ia*ja*ka)
      enddo
      enddo
      errornorm2max = pmax(errornorm2max)
      errornorm2avg = psum(errornorm2avg)/(1.d0*ia*ja*ka)
      errornorm2gtr10 = psum(errornorm2gtr10)/(1.d0*ia*ja*ka)
      errornorm2gtr100 = psum(errornorm2gtr100)/(1.d0*ia*ja*ka)
      !
      if(mpirank==0) then
        write(mname,'(i4.4)')m
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_stress_relative_error_'//stepname//'_'//mname//'.dat'
        else
          outfilename = 'pp/SGS_stress_relative_error_'//mname//'.dat'
        endif
        !
        if(mpirank == 0)then
          open(fh,file=outfilename,form='formatted')
          write(fh,"(A7,1x,2(A20,1x))")'nstep','time','l'
          write(fh,"(I7,1x,2(E20.13E2,1x))")nstep,time,l_lim(m)
          write(fh,"(A8,1x,10(A20,1x))")'type','tau11','tau12','tau13','tau21','tau22','tau23','tau31','tau32','tau33','norm2'
          write(fh,"(A8,1x,10(E20.13E2,1x))")'max',((errormax(p,q),q=1,3),p=1,3),errornorm2max
          write(fh,"(A8,1x,10(E20.13E2,1x))")'avg',((erroravg(p,q),q=1,3),p=1,3),errornorm2avg
          write(fh,"(A8,1x,10(E20.13E2,1x))")'gtr0.1',((errorgtr10(p,q),q=1,3),p=1,3),errornorm2gtr10
          write(fh,"(A8,1x,10(E20.13E2,1x))")'gtr1',((errorgtr100(p,q),q=1,3),p=1,3),errornorm2gtr100
          close(fh)
          print *, '>>>>', outfilename
        endif
        !
        !
        !
      endif
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(loutput)then
        !
        write(mname,'(i4.4)')m
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_stress_'//stepname//'_'//mname//'.h5'
        else
          outfilename = 'pp/SGS_stress_'//mname//'.h5'
        endif
        !
        call h5io_init(trim(outfilename),mode='write')
        !
        do p=1,3
          do q= 1,3
            write (termname, "(A3,I1,I1)") "tau",p,q
            call h5write(var=tau(p,q,1:im,1:jm,1:km),      varname=termname,    mode = modeio) 
            write (termname, "(A3,I1,I1,A3)") "tau",p,q,"bis"
            call h5write(var=tau_bis(p,q,1:im,1:jm,1:km),  varname=termname,    mode = modeio)
          enddo
        enddo
        !
        call h5io_end
        !
      endif
      !
    enddo ! loop of filter point l
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_w1)
    call fftw_free(c_w2)
    call fftw_free(c_w3)
    call fftw_free(c_rhol)
    call fftw_free(c_w1_filted)
    call fftw_free(c_w2_filted)
    call fftw_free(c_w3_filted)
    call fftw_free(c_rho_filted)
    call fftw_free(c_w1w1)
    call fftw_free(c_w1w2)
    call fftw_free(c_w1w3)
    call fftw_free(c_w2w1)
    call fftw_free(c_w2w2)
    call fftw_free(c_w2w3)
    call fftw_free(c_w3w1)
    call fftw_free(c_w3w2)
    call fftw_free(c_w3w3)
    call fftw_free(c_w1w1_filted)
    call fftw_free(c_w1w2_filted)
    call fftw_free(c_w1w3_filted)
    call fftw_free(c_w2w1_filted)
    call fftw_free(c_w2w2_filted)
    call fftw_free(c_w2w3_filted)
    call fftw_free(c_w3w1_filted)
    call fftw_free(c_w3w2_filted)
    call fftw_free(c_w3w3_filted)
    call fftw_free(c_A11)
    call fftw_free(c_A12)
    call fftw_free(c_A13)
    call fftw_free(c_A21)
    call fftw_free(c_A22)
    call fftw_free(c_A23)
    call fftw_free(c_A31)
    call fftw_free(c_A32)
    call fftw_free(c_A33)
    call fftw_free(c_tau11_term)
    call fftw_free(c_tau12_term)
    call fftw_free(c_tau13_term)
    call fftw_free(c_tau21_term)
    call fftw_free(c_tau22_term)
    call fftw_free(c_tau23_term)
    call fftw_free(c_tau31_term)
    call fftw_free(c_tau32_term)
    call fftw_free(c_tau33_term)
    call mpistop
    deallocate(tau,tau_bis,erroravg,errormax,errorgtr10,errorgtr100)
    !
  end subroutine SGSstress3D
  !
  subroutine SGST3D(thefilenumb)
    ! 
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km
    use commarray, only: vel, rho, prs
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
    include 'fftw3-mpi.f03'
    !
    integer,intent(in) :: thefilenumb
    integer :: i,j,k,m,n
    character(len=128) :: infilename,outfilename
    character(len=4) :: stepname
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom,pcom
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
    complex(8) :: imag
    real(8),allocatable,dimension(:) :: sqrtalphas,dalphas
    integer :: num_l,num_alpha,num_alphamin
    integer :: hand_a
    real(8) :: l_min, ratio_max,ratio_min
    real(8) :: Galpha
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted,p_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11_filted,A12_filted,A13_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21_filted,A22_filted,A23_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31_filted,A32_filted,A33_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: p11_filted,p12_filted,p13_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: p21_filted,p22_filted,p23_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: p31_filted,p32_filted,p33_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: p1_filted,p2_filted,p3_filted
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: rho1_filted,rho2_filted,rho3_filted
    complex(8), allocatable, dimension(:,:,:) :: theta_filted,M11_filted,M22_filted,M33_filted
    complex(8), allocatable, dimension(:,:,:) :: M12_filted,M13_filted,M21_filted
    complex(8), allocatable, dimension(:,:,:) :: M23_filted,M31_filted,M32_filted
    real(8) :: Es,Ec,Ts,Tc,Ps1,Pc1,Pc2
    !
    complex(8) :: vxr_D
    !
    type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,c_pcom,forward_plan,backward_plan
    type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted,c_p_filted
    type(C_PTR) :: c_A11_filted,c_A12_filted,c_A13_filted
    type(C_PTR) :: c_A21_filted,c_A22_filted,c_A23_filted
    type(C_PTR) :: c_A31_filted,c_A32_filted,c_A33_filted
    type(C_PTR) :: c_p11_filted,c_p12_filted,c_p13_filted
    type(C_PTR) :: c_p21_filted,c_p22_filted,c_p23_filted
    type(C_PTR) :: c_p31_filted,c_p32_filted,c_p33_filted
    type(C_PTR) :: c_p1_filted,c_p2_filted,c_p3_filted
    type(C_PTR) :: c_rho1_filted,c_rho2_filted,c_rho3_filted
    !
    integer,dimension(8) :: value
    character(len=1) :: modeio
    logical :: loutput
    !
    call readinput
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !!!! Read velocity and density field
    allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km),prs(0:im,0:jm,0:km))
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='p',  var=prs(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    !!!! Prepare initial field in Fourier space
    !! velocity
    c_w1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
    c_w2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
    c_w3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
    c_rhocom = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
    c_pcom = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_pcom, pcom, [imfftw,jmfftw,kmfftw])
    !
    forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T)
      w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T)
      w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T)
      rhocom(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T)
      pcom(i,j,k)=CMPLX(prs(i,j,k),0.d0,C_INTPTR_T)
      !
    end do
    end do
    end do
    !
    !After this bloc, w1 is (rho*u1) in spectral space
    call fftw_mpi_execute_dft(forward_plan,w1,w1)
    call fftw_mpi_execute_dft(forward_plan,w2,w2)
    call fftw_mpi_execute_dft(forward_plan,w3,w3)
    call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
    call fftw_mpi_execute_dft(forward_plan,pcom,pcom)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
      w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
      w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
      rhocom(i,j,k)=rhocom(i,j,k)/(1.d0*ia*ja*ka)
      pcom(i,j,k)=pcom(i,j,k)/(1.d0*ia*ja*ka)
      !
    end do
    end do
    end do

    !
    !
    !! wavenumber
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    do k = 1,km
    do j = 1,jm
    do i = 1,im
      !
      if(im .ne. ia)then
        stop "error! im /= ia"
      endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j,k) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j,k) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if(j <= (ja/2+1)) then
        k2(i,j,k) = real(j-1,8)
      else if(i<=(ia)) then
        k2(i,j,k) = real(j-ja-1,8)
      else
        print *,"Error, no wave number possible, j must smaller than ja-1 !"
      end if
      !
      if((k+k0) <= (ka/2+1)) then
        k3(i,j,k) = real(k+k0-1,8)
      else if((k+k0)<=(ka)) then
        k3(i,j,k) = real(k+k0-ka-1,8)
      else
        print *,"Error, no wave number possible, (k+k0) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    end do
    !
    !! Imaginary number prepare
    imag = CMPLX(0.d0,1.d0,8)
    !
    if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
    !!!! Prepare alpha and others
    call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
    l_min = 2*pi/im
    allocate(sqrtalphas(num_alpha),dalphas(num_alpha))
    !
    do i=1,num_alpha
      sqrtalphas(i) = sqrt( exp(log(ratio_max**2) * (i-1) / (num_alpha-1)) ) * l_min
    enddo
    !
    dalphas(1) = sqrtalphas(1)**2 
    !
    do i=2,num_alpha
      dalphas(i) = sqrtalphas(i)**2 - sqrtalphas(i-1)**2 
    enddo
    !
    if(mpirank==0)  print *, "Integrate point allocated"
    !
    !!!!
    !
    c_w1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw,kmfftw])
    c_w2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw,kmfftw])
    c_w3_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_w3_filted, w3_filted,  [imfftw,jmfftw,kmfftw])
    c_rho_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw,kmfftw])
    c_p_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_p_filted, p_filted,[imfftw,jmfftw,kmfftw])
    !
    c_A11_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A11_filted, A11_filted,[imfftw,jmfftw,kmfftw])
    c_A12_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A12_filted, A12_filted,[imfftw,jmfftw,kmfftw])
    c_A13_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A13_filted, A13_filted,[imfftw,jmfftw,kmfftw])
    c_A21_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A21_filted, A21_filted,[imfftw,jmfftw,kmfftw])
    c_A22_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A22_filted, A22_filted,[imfftw,jmfftw,kmfftw])
    c_A23_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A23_filted, A23_filted,[imfftw,jmfftw,kmfftw])
    c_A31_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A31_filted, A31_filted,[imfftw,jmfftw,kmfftw])
    c_A32_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A32_filted, A32_filted,[imfftw,jmfftw,kmfftw])
    c_A33_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_A33_filted, A33_filted,[imfftw,jmfftw,kmfftw])
    !
    c_p11_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_p11_filted, p11_filted,[imfftw,jmfftw,kmfftw])
    c_p12_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_p12_filted, p12_filted,[imfftw,jmfftw,kmfftw])
    c_p13_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_p13_filted, p13_filted,[imfftw,jmfftw,kmfftw])
    c_p21_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_p21_filted, p21_filted,[imfftw,jmfftw,kmfftw])
    c_p22_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_p22_filted, p22_filted,[imfftw,jmfftw,kmfftw])
    c_p23_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_p23_filted, p23_filted,[imfftw,jmfftw,kmfftw])
    c_p31_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_p31_filted, p31_filted,[imfftw,jmfftw,kmfftw])
    c_p32_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_p32_filted, p32_filted,[imfftw,jmfftw,kmfftw])
    c_p33_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_p33_filted, p33_filted,[imfftw,jmfftw,kmfftw])
    c_p1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_p1_filted, p1_filted,[imfftw,jmfftw,kmfftw])
    c_p2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_p2_filted, p2_filted,[imfftw,jmfftw,kmfftw])
    c_p3_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_p3_filted, p3_filted,[imfftw,jmfftw,kmfftw])
    c_rho1_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rho1_filted, rho1_filted,[imfftw,jmfftw,kmfftw])
    c_rho2_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rho2_filted, rho2_filted,[imfftw,jmfftw,kmfftw])
    c_rho3_filted = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rho3_filted, rho3_filted,[imfftw,jmfftw,kmfftw])
    !
    !
    allocate(theta_filted(1:im,1:jm,1:km),&
    M11_filted(1:im,1:jm,1:km),M22_filted(1:im,1:jm,1:km),M33_filted(1:im,1:jm,1:km),&
    M12_filted(1:im,1:jm,1:km),M13_filted(1:im,1:jm,1:km),M21_filted(1:im,1:jm,1:km),&
    M23_filted(1:im,1:jm,1:km),M31_filted(1:im,1:jm,1:km),M32_filted(1:im,1:jm,1:km))
    !
    !
    !
    Es = 0.d0
    Ec = 0.d0
    Ts = 0.d0
    Tc = 0.d0
    Ps1 = 0.d0
    Pc1 = 0.d0
    Pc2 = 0.d0
    !
    if(mpirank==0)  print *, "Array allocated and initialized"
    !
    do n=1,num_alpha
      !
      call date_and_time(values=value) 
      !
      if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alpha,',now is ',&
                              value(5), ':', value(6),':',value(7)
      !
      !!! Velocity Favre average and density average
      ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
      do i=1,im
      do j=1,jm
      do k=1,km
        Galpha = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*sqrtalphas(n)**2/2.d0) ! Filtre scale :sqrtalpha
        w1_filted(i,j,k)  = w1(i,j,k)    *Galpha
        w2_filted(i,j,k)  = w2(i,j,k)    *Galpha
        w3_filted(i,j,k)  = w3(i,j,k)    *Galpha
        rho_filted(i,j,k) = rhocom(i,j,k)*Galpha
        p_filted(i,j,k)   = pcom(i,j,k)  *Galpha
        !
        rho1_filted(i,j,k) = imag*rho_filted(i,j,k)*k1(i,j,k)
        rho2_filted(i,j,k) = imag*rho_filted(i,j,k)*k2(i,j,k)
        rho3_filted(i,j,k) = imag*rho_filted(i,j,k)*k3(i,j,k)
        !
        p11_filted(i,j,k)  = -p_filted(i,j,k)*k1(i,j,k)*k1(i,j,k)
        p21_filted(i,j,k)  = -p_filted(i,j,k)*k2(i,j,k)*k1(i,j,k)
        p31_filted(i,j,k)  = -p_filted(i,j,k)*k3(i,j,k)*k1(i,j,k)
        p12_filted(i,j,k)  = -p_filted(i,j,k)*k1(i,j,k)*k2(i,j,k)
        p22_filted(i,j,k)  = -p_filted(i,j,k)*k2(i,j,k)*k2(i,j,k)
        p32_filted(i,j,k)  = -p_filted(i,j,k)*k3(i,j,k)*k2(i,j,k)
        p13_filted(i,j,k)  = -p_filted(i,j,k)*k1(i,j,k)*k3(i,j,k)
        p23_filted(i,j,k)  = -p_filted(i,j,k)*k2(i,j,k)*k3(i,j,k)
        p33_filted(i,j,k)  = -p_filted(i,j,k)*k3(i,j,k)*k3(i,j,k)
        p1_filted(i,j,k)   = imag*p_filted(i,j,k)*k1(i,j,k)
        p2_filted(i,j,k)   = imag*p_filted(i,j,k)*k2(i,j,k)
        p3_filted(i,j,k)   = imag*p_filted(i,j,k)*k3(i,j,k)
        !
      enddo
      enddo
      enddo
      !
      ! After this bloc, w1_filted is (rho*u1)_filted in physical space
      call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
      call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
      call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
      call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
      call fftw_mpi_execute_dft(backward_plan,p_filted,p_filted)
      !
      call fftw_mpi_execute_dft(backward_plan,p11_filted,p11_filted)
      call fftw_mpi_execute_dft(backward_plan,p21_filted,p21_filted)
      call fftw_mpi_execute_dft(backward_plan,p31_filted,p31_filted)
      call fftw_mpi_execute_dft(backward_plan,p12_filted,p12_filted)
      call fftw_mpi_execute_dft(backward_plan,p22_filted,p22_filted)
      call fftw_mpi_execute_dft(backward_plan,p32_filted,p32_filted)
      call fftw_mpi_execute_dft(backward_plan,p13_filted,p13_filted)
      call fftw_mpi_execute_dft(backward_plan,p23_filted,p23_filted)
      call fftw_mpi_execute_dft(backward_plan,p33_filted,p33_filted)
      call fftw_mpi_execute_dft(backward_plan,p1_filted,p1_filted)
      call fftw_mpi_execute_dft(backward_plan,p2_filted,p2_filted)
      call fftw_mpi_execute_dft(backward_plan,p3_filted,p3_filted)
      call fftw_mpi_execute_dft(backward_plan,rho1_filted,rho1_filted)
      call fftw_mpi_execute_dft(backward_plan,rho2_filted,rho2_filted)
      call fftw_mpi_execute_dft(backward_plan,rho3_filted,rho3_filted)
      !
      ! After this bloc, w1_filted is u1_filted in physical space
      do i=1,im
      do j=1,jm
      do k=1,km
        !
        w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
        w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
        w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
        !
      enddo
      enddo
      enddo
      !
      ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
      call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
      call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
      call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1_filted(i,j,k)   = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
        w2_filted(i,j,k)   = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
        w3_filted(i,j,k)   = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
        !
        A11_filted(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
        A21_filted(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
        A31_filted(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
        A12_filted(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
        A22_filted(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
        A32_filted(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
        A13_filted(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
        A23_filted(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
        A33_filted(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
        !
      end do
      end do
      end do
      !
      ! After this bloc, A11_filted is A11_filted in physical space
      call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
      call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
      call fftw_mpi_execute_dft(backward_plan,A31_filted,A31_filted)
      call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
      call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
      call fftw_mpi_execute_dft(backward_plan,A32_filted,A32_filted)
      call fftw_mpi_execute_dft(backward_plan,A13_filted,A13_filted)
      call fftw_mpi_execute_dft(backward_plan,A23_filted,A23_filted)
      call fftw_mpi_execute_dft(backward_plan,A33_filted,A33_filted)
      !
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        theta_filted(i,j,k) = A11_filted(i,j,k)+A22_filted(i,j,k)+A33_filted(i,j,k)
        !
        M11_filted(i,j,k) = A11_filted(i,j,k) - 1.d0/3.d0*theta_filted(i,j,k)
        M22_filted(i,j,k) = A22_filted(i,j,k) - 1.d0/3.d0*theta_filted(i,j,k)
        M33_filted(i,j,k) = A33_filted(i,j,k) - 1.d0/3.d0*theta_filted(i,j,k)
        !
        M12_filted(i,j,k) = A12_filted(i,j,k)
        M13_filted(i,j,k) = A13_filted(i,j,k)
        M21_filted(i,j,k) = A21_filted(i,j,k)
        M23_filted(i,j,k) = A23_filted(i,j,k)
        M31_filted(i,j,k) = A31_filted(i,j,k)
        M32_filted(i,j,k) = A32_filted(i,j,k)
        !
      end do
      end do
      end do
      !
      !!!! T terms
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        Es = Es + dreal(0.5d0 * rho_filted(i,j,k)*&
        (M11_filted(i,j,k)*M11_filted(i,j,k) + M12_filted(i,j,k)*M12_filted(i,j,k) + M13_filted(i,j,k)*M13_filted(i,j,k) + &
        M21_filted(i,j,k)*M21_filted(i,j,k) + M22_filted(i,j,k)*M22_filted(i,j,k) + M23_filted(i,j,k)*M23_filted(i,j,k) + &
        M31_filted(i,j,k)*M31_filted(i,j,k) + M32_filted(i,j,k)*M32_filted(i,j,k) + M33_filted(i,j,k)*M33_filted(i,j,k))*dalphas(n))
        !
        Ec = Ec + dreal(1.d0/6.d0 * rho_filted(i,j,k)* theta_filted(i,j,k)*theta_filted(i,j,k) *dalphas(n))
        !
        Ts = Ts + dreal(rho_filted(i,j,k)*&
        (M11_filted(i,j,k)*M11_filted(i,j,k)*M11_filted(i,j,k)+M12_filted(i,j,k)*M21_filted(i,j,k)*M11_filted(i,j,k)+&
        M13_filted(i,j,k)*M31_filted(i,j,k)*M11_filted(i,j,k)+M11_filted(i,j,k)*M12_filted(i,j,k)*M12_filted(i,j,k)+&
        M12_filted(i,j,k)*M22_filted(i,j,k)*M12_filted(i,j,k)+M13_filted(i,j,k)*M32_filted(i,j,k)*M12_filted(i,j,k)+&
        M11_filted(i,j,k)*M13_filted(i,j,k)*M13_filted(i,j,k)+M12_filted(i,j,k)*M23_filted(i,j,k)*M13_filted(i,j,k)+&
        M13_filted(i,j,k)*M33_filted(i,j,k)*M13_filted(i,j,k)+M21_filted(i,j,k)*M11_filted(i,j,k)*M21_filted(i,j,k)+&
        M22_filted(i,j,k)*M21_filted(i,j,k)*M21_filted(i,j,k)+M23_filted(i,j,k)*M31_filted(i,j,k)*M21_filted(i,j,k)+&
        M21_filted(i,j,k)*M12_filted(i,j,k)*M22_filted(i,j,k)+M22_filted(i,j,k)*M22_filted(i,j,k)*M22_filted(i,j,k)+&
        M23_filted(i,j,k)*M32_filted(i,j,k)*M22_filted(i,j,k)+M21_filted(i,j,k)*M13_filted(i,j,k)*M23_filted(i,j,k)+&
        M22_filted(i,j,k)*M23_filted(i,j,k)*M23_filted(i,j,k)+M23_filted(i,j,k)*M33_filted(i,j,k)*M23_filted(i,j,k)+&
        M31_filted(i,j,k)*M11_filted(i,j,k)*M31_filted(i,j,k)+M32_filted(i,j,k)*M21_filted(i,j,k)*M31_filted(i,j,k)+&
        M33_filted(i,j,k)*M31_filted(i,j,k)*M31_filted(i,j,k)+M31_filted(i,j,k)*M12_filted(i,j,k)*M32_filted(i,j,k)+&
        M32_filted(i,j,k)*M22_filted(i,j,k)*M32_filted(i,j,k)+M33_filted(i,j,k)*M32_filted(i,j,k)*M32_filted(i,j,k)+&
        M31_filted(i,j,k)*M13_filted(i,j,k)*M33_filted(i,j,k)+M32_filted(i,j,k)*M23_filted(i,j,k)*M33_filted(i,j,k)+&
        M33_filted(i,j,k)*M33_filted(i,j,k)*M33_filted(i,j,k)+&
        2.d0/3.d0 * &
        (M11_filted(i,j,k)*M11_filted(i,j,k) + M12_filted(i,j,k)*M12_filted(i,j,k) + M13_filted(i,j,k)*M13_filted(i,j,k)+ &
        M21_filted(i,j,k)*M21_filted(i,j,k) + M22_filted(i,j,k)*M22_filted(i,j,k) + M23_filted(i,j,k)*M23_filted(i,j,k) + &
        M31_filted(i,j,k)*M31_filted(i,j,k) + M32_filted(i,j,k)*M32_filted(i,j,k) + M33_filted(i,j,k)*M33_filted(i,j,k)) &
        *theta_filted(i,j,k))*dalphas(n))
        !
        Tc = Tc + dreal(rho_filted(i,j,k)*&
        ((M11_filted(i,j,k)*M11_filted(i,j,k) + M12_filted(i,j,k)*M12_filted(i,j,k) + M13_filted(i,j,k)*M13_filted(i,j,k) + &
        M21_filted(i,j,k)*M21_filted(i,j,k) + M22_filted(i,j,k)*M22_filted(i,j,k) + M23_filted(i,j,k)*M23_filted(i,j,k) + &
        M31_filted(i,j,k)*M31_filted(i,j,k) + M32_filted(i,j,k)*M32_filted(i,j,k) + M33_filted(i,j,k)*M33_filted(i,j,k)) &
        *theta_filted(i,j,k) - 1.d0/3.d0*theta_filted(i,j,k)**3)*dalphas(n))
        !
        Ps1 = Ps1 + dreal((&
        (p11_filted(i,j,k)*A11_filted(i,j,k)+p12_filted(i,j,k)*A12_filted(i,j,k)+p13_filted(i,j,k)*A13_filted(i,j,k) &
        +p21_filted(i,j,k)*A21_filted(i,j,k)+p22_filted(i,j,k)*A22_filted(i,j,k)+p23_filted(i,j,k)*A23_filted(i,j,k) &
        +p31_filted(i,j,k)*A31_filted(i,j,k)+p32_filted(i,j,k)*A32_filted(i,j,k)+p33_filted(i,j,k)*A33_filted(i,j,k))& 
        - 1.d0/3.d0 * (p11_filted(i,j,k)+p22_filted(i,j,k)+p33_filted(i,j,k)) * theta_filted(i,j,k)&
        )*dalphas(n))
        !
        Pc1 = Pc1 + dreal(1.d0/3.d0 * (p11_filted(i,j,k)+p22_filted(i,j,k)+p33_filted(i,j,k)) * &
              theta_filted(i,j,k)*dalphas(n))
        !
        Pc2 = Pc2 + dreal(1.d0/3.d0 / rho_filted(i,j,k) * (rho1_filted(i,j,k)*p1_filted(i,j,k) + &
        rho2_filted(i,j,k)*p2_filted(i,j,k) + rho3_filted(i,j,k)*p3_filted(i,j,k))* theta_filted(i,j,k) *dalphas(n))
        !
      enddo
      enddo
      enddo
      !
    enddo
    !
    Es = psum(Es)/(1.d0*ia*ja*ka)
    Ec = psum(Ec)/(1.d0*ia*ja*ka)
    Ts = psum(Ts)/(1.d0*ia*ja*ka)
    Tc = psum(Tc)/(1.d0*ia*ja*ka)
    Ps1 = psum(Ps1)/(1.d0*ia*ja*ka)
    Pc1 = psum(Pc1)/(1.d0*ia*ja*ka)
    Pc2 = psum(Pc2)/(1.d0*ia*ja*ka)
    !
    if(mpirank==0)  print *, 'Job finish'
    !
    if(mpirank==0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/SGS_T_'//stepname//'.dat'
      else
        outfilename = 'pp/SGS_T.dat'
      endif
      
      call listinit(filename=outfilename,handle=hand_a, &
                    firstline='nstep time Es Ec Ts Tc Ps1 Pc1 Pc2')
      call listwrite(hand_a,Es,Ec,Ts,Tc,Ps1,Pc1,Pc2)
      !
      print *, '>>>>', outfilename
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_w1)
    call fftw_free(c_w2)
    call fftw_free(c_w3)
    call fftw_free(c_rhocom)
    call fftw_free(c_w1_filted)
    call fftw_free(c_w2_filted)
    call fftw_free(c_w3_filted)
    call fftw_free(c_rho_filted)
    call fftw_free(c_A11_filted)
    call fftw_free(c_A12_filted)
    call fftw_free(c_A13_filted)
    call fftw_free(c_A21_filted)
    call fftw_free(c_A22_filted)
    call fftw_free(c_A23_filted)
    call fftw_free(c_A31_filted)
    call fftw_free(c_A32_filted)
    call fftw_free(c_A33_filted)
    call fftw_free(c_p11_filted)
    call fftw_free(c_p12_filted)
    call fftw_free(c_p13_filted)
    call fftw_free(c_p21_filted)
    call fftw_free(c_p22_filted)
    call fftw_free(c_p23_filted)
    call fftw_free(c_p31_filted)
    call fftw_free(c_p32_filted)
    call fftw_free(c_p33_filted)
    call fftw_free(c_p1_filted)
    call fftw_free(c_p2_filted)
    call fftw_free(c_p3_filted)
    call fftw_free(c_rho1_filted)
    call fftw_free(c_rho2_filted)
    call fftw_free(c_rho3_filted)
    call mpistop
    deallocate(theta_filted,M11_filted,M22_filted,M33_filted)
    deallocate(M12_filted,M13_filted,M21_filted)
    deallocate(M23_filted,M31_filted,M32_filted)
    deallocate(k1,k2,k3)
    deallocate(sqrtalphas,dalphas)
    !
  end subroutine SGST3D
  !
  !
  subroutine readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
    !
    use parallel,only: bcast,mpirank
    !
    ! local data
    integer, intent(out) :: num_l,num_alpha,num_alphamin
    real(8), intent(out) :: ratio_max,ratio_min
    logical, intent(out) :: loutput
    character(len=64) :: inputfile
    integer :: fh
    !
    inputfile='datin/SGSinput'
    !
    if(mpirank==0) then
      !
      fh=get_unit()
      !
      open(fh,file=trim(inputfile),action='read')
      read(fh,'(//)')
      read(fh,*)num_l,num_alpha,num_alphamin
      read(fh,'(/)')
      read(fh,*)ratio_max,ratio_min
      read(fh,'(/)')
      read(fh,*)loutput
      close(fh)
      print*,' >> ',trim(inputfile),' ... done'
      print*,' >>> Get: Number of l is',num_l,'Number of alpha is',num_alpha,'Minimum number of alpha is',num_alphamin
      print*,' >>> Ratio max:',ratio_max,'Ratio min',ratio_min
      if(loutput)then
        print *, ' >>> No output stress'
      else
        print *, ' >>> Output stress'
      endif
      !
    endif
    !
    call bcast(num_l)
    call bcast(num_alpha)
    call bcast(num_alphamin)
    call bcast(ratio_max)
    call bcast(ratio_min)
    call bcast(loutput)
    !
  end subroutine readSGSinput
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This function is used to calcuate the spectral energy at any 
  ! wavenumber.
  ! Ref: S. JAMME, et al. Direct Numerical Simulation of the 
  ! Interaction between a Shock Wave and Various Types of Isotropic 
  ! Turbulence, Flow, Turbulence and Combustion, 2002, 68:227每268.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function IniEnergDis(Ac,k0,wnb)
    !
    real(8) :: k0,Ac,var1,wnb,IniEnergDis
    !
    var1=-2.d0*(wnb/k0)**2
    IniEnergDis=Ac*wnb**4*exp(var1)
    !IniEnergDis=Ac*wnb**(-5.d0/3.d0)
    !
    return
    !
  end function IniEnergDis
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function Ek.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  real(8) function wav(i,im)
  ! This function gives the wave number of index i
  ! with the maximum im
    implicit none
    integer, intent(in) :: i, im
    !
    if(i <= (im/2+1)) then
      wav = real(i-1,8)
    else if(i<=im) then
      wav = real(i-im-1,8)
    else
      print *,"Error, no wave number possible, i must smaller than im !"
    end if
  end function wav
  !
  integer function invwav(i,im)
  ! This function gives the wave number of index i
  ! with the maximum im
    implicit none
    integer, intent(in) :: i,im
    !
    if(i < 0) then
      invwav = i + im + 1;
    else
      invwav = i + 1;
    end if
  end function invwav
  !
  integer function kint(k,dk)
  ! This function gives the nearby k
  !!
    implicit none
    real(8), intent(in) :: k,dk
    !
    if(nint(k)>=k-dk) then
      kint = nint(k)
    else if(nint(k)+1<=k+dk) then
      kint = nint(k)+1
    else
      print *,"error in kint!"
      kint = 0
    end if
  end function kint
  !
  real(8) function ProjectP2_2D(i,j,kx,ky)
    !
    !!
    implicit none
    integer, intent(in) :: i,j
    real(8), intent(in) ::kx,ky
    real(8) :: k
    !
    k = dsqrt(kx**2+ky**2+1.d-15)
    !
    if((i .eq. 1) .and. (j .eq. 1)) then
      ProjectP2_2D = 1.d0 - kx*kx/k/k
    else if((i .eq. 1) .and. (j .eq. 2)) then
      ProjectP2_2D = - kx*ky/k/k
    else if((i .eq. 2) .and. (j .eq. 1)) then
      ProjectP2_2D = - kx*ky/k/k
    else if((i .eq. 2) .and. (j .eq. 2)) then
      ProjectP2_2D = 1.d0 - ky*ky/k/k
    else
      stop "ProjectP2_2D error: i,j"
    end if
    !
  end function ProjectP2_2D
  !
  real(8) function ProjectP2_3D(i,j,kx,ky,kz)
  !
  !!
  implicit none
  integer, intent(in) :: i,j
  real(8), intent(in) :: kx,ky,kz
  real(8) :: k
  !
  k = dsqrt(kx**2+ky**2+kz**2+1.d-15)
  !
  if((i .eq. 1) .and. (j .eq. 1)) then
    ProjectP2_3D = 1.d0 - kx*kx/k/k
  else if((i .eq. 1) .and. (j .eq. 2)) then
    ProjectP2_3D = - kx*ky/k/k
  else if((i .eq. 1) .and. (j .eq. 3)) then
    ProjectP2_3D = - kx*kz/k/k
  else if((i .eq. 2) .and. (j .eq. 1)) then
    ProjectP2_3D = - ky*kx/k/k
  else if((i .eq. 2) .and. (j .eq. 2)) then
    ProjectP2_3D = 1.d0 - ky*ky/k/k
  else if((i .eq. 2) .and. (j .eq. 3)) then
    ProjectP2_3D = - ky*kz/k/k
  else if((i .eq. 3) .and. (j .eq. 1)) then
    ProjectP2_3D = - kz*kx/k/k
  else if((i .eq. 3) .and. (j .eq. 2)) then
    ProjectP2_3D = - kz*ky/k/k
  else if((i .eq. 3) .and. (j .eq. 3)) then
    ProjectP2_3D = 1.d0 - kz*kz/k/k
  else
    stop "ProjectP2_3D error: i,j"
  end if
  !
  end function ProjectP2_3D
  !
  real(8) function ProjectP3_2D(i,j,m,kx,ky)
    !
    !!
    implicit none
    integer, intent(in) :: i,j,m
    real(8), intent(in) :: kx,ky
    !
    if((j .eq. 1) .and. (m .eq. 1))then
      ProjectP3_2D = 0.5d0 * kx * ProjectP2(i,j,kx,ky)+ 0.5d0 * kx * ProjectP2(i,m,kx,ky)
    else if((j .eq. 2) .and. (m .eq. 1))then
      ProjectP3_2D = 0.5d0 * kx * ProjectP2(i,j,kx,ky)+ 0.5d0 * ky * ProjectP2(i,m,kx,ky)
    else if((j .eq. 1) .and. (m .eq. 2))then
      ProjectP3_2D = 0.5d0 * ky * ProjectP2(i,j,kx,ky)+ 0.5d0 * kx * ProjectP2(i,m,kx,ky)
    else if((j .eq. 2) .and. (m .eq. 2))then
      ProjectP3_2D = 0.5d0 * ky * ProjectP2(i,j,kx,ky)+ 0.5d0 * ky * ProjectP2(i,m,kx,ky)
    else
      stop "ProjectP3_2D error: j,m"
    endif
  end function ProjectP3_2D
  !
  real(8) function ProjectP3_3D(i,j,m,kx,ky,kz)
  !
  !!
  implicit none
  integer, intent(in) :: i,j,m
  real(8), intent(in) :: kx,ky,kz
  !
  if((j .eq. 1) .and. (m .eq. 1))then
    ProjectP3_3D = 0.5d0 * kx * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * kx * ProjectP2(i,m,kx,ky,kz)
  else if((j .eq. 2) .and. (m .eq. 1))then
    ProjectP3_3D = 0.5d0 * kx * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * ky * ProjectP2(i,m,kx,ky,kz)
  else if((j .eq. 3) .and. (m .eq. 1))then
    ProjectP3_3D = 0.5d0 * kx * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * kz * ProjectP2(i,m,kx,ky,kz)
  else if((j .eq. 1) .and. (m .eq. 2))then
    ProjectP3_3D = 0.5d0 * ky * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * kx * ProjectP2(i,m,kx,ky,kz)
  else if((j .eq. 2) .and. (m .eq. 2))then
    ProjectP3_3D = 0.5d0 * ky * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * ky * ProjectP2(i,m,kx,ky,kz)
  else if((j .eq. 3) .and. (m .eq. 2))then
    ProjectP3_3D = 0.5d0 * ky * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * kz * ProjectP2(i,m,kx,ky,kz)
  else if((j .eq. 1) .and. (m .eq. 3))then
    ProjectP3_3D = 0.5d0 * kz * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * kx * ProjectP2(i,m,kx,ky,kz)
  else if((j .eq. 2) .and. (m .eq. 3))then
    ProjectP3_3D = 0.5d0 * kz * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * ky * ProjectP2(i,m,kx,ky,kz)
  else if((j .eq. 3) .and. (m .eq. 3))then
    ProjectP3_3D = 0.5d0 * kz * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * kz * ProjectP2(i,m,kx,ky,kz)
  else
    stop "ProjectP3_3D error: j,m"
  endif
  end function ProjectP3_3D
  !
  real(8) function ProjectPi2_2D(i,j,kx,ky)
    !
    !!
    implicit none
    integer, intent(in) :: i,j
    real(8), intent(in) :: kx,ky
    real(8) :: k
    !
    k = dsqrt(kx**2+ky**2+1.d-15)
    !
    if((i .eq. 1) .and. (j .eq. 1)) then
      ProjectPi2_2D = kx*kx/k/k
    else if((i .eq. 1) .and. (j .eq. 2)) then
      ProjectPi2_2D = kx*ky/k/k
    else if((i .eq. 2) .and. (j .eq. 1)) then
      ProjectPi2_2D = kx*ky/k/k
    else if((i .eq. 2) .and. (j .eq. 2)) then
      ProjectPi2_2D = ky*ky/k/k
    else
      stop "ProjectPi2_2D error: i,j"
    end if
    !
  end function ProjectPi2_2D
  !
  real(8) function ProjectPi2_3D(i,j,kx,ky,kz)
  !
  !!
  implicit none
  integer, intent(in) :: i,j
  real(8), intent(in) :: kx,ky,kz
  real(8) :: k
  !
  k = dsqrt(kx**2+ky**2+kz**2+1.d-15)
  !
  if((i .eq. 1) .and. (j .eq. 1)) then
    ProjectPi2_3D = kx*kx/k/k
  else if((i .eq. 1) .and. (j .eq. 2)) then
    ProjectPi2_3D = kx*ky/k/k
  else if((i .eq. 1) .and. (j .eq. 3)) then
    ProjectPi2_3D = kx*kz/k/k
  else if((i .eq. 2) .and. (j .eq. 1)) then
    ProjectPi2_3D = ky*kx/k/k
  else if((i .eq. 2) .and. (j .eq. 2)) then
    ProjectPi2_3D = ky*ky/k/k
  else if((i .eq. 2) .and. (j .eq. 3)) then
    ProjectPi2_3D = ky*kz/k/k
  else if((i .eq. 3) .and. (j .eq. 1)) then
    ProjectPi2_3D = kz*kx/k/k
  else if((i .eq. 3) .and. (j .eq. 2)) then
    ProjectPi2_3D = kz*ky/k/k
  else if((i .eq. 3) .and. (j .eq. 3)) then
    ProjectPi2_3D = kz*kz/k/k
  else
    stop "ProjectPi2_3D error: i,j"
  end if
  !
  end function ProjectPi2_3D
  !
  real(8) function ProjectPi3_2D(i,j,m,kx,ky)
    !
    !!
    implicit none
    integer, intent(in) :: i,j,m
    real(8), intent(in) :: kx,ky
    !
    if(m .eq. 1)then
      ProjectPi3_2D = kx * ProjectP2(i,j,kx,ky)
    else if(m .eq. 2)then
      ProjectPi3_2D = ky * ProjectP2(i,j,kx,ky)
    else
      stop "ProjectPi3_2D error: m"
    endif
  end function ProjectPi3_2D
  !
  real(8) function ProjectPi3_3D(i,j,m,kx,ky,kz)
    !
    !!
    implicit none
    integer, intent(in) :: i,j,m
    real(8), intent(in) :: kx,ky,kz
    !
    if(m .eq. 1)then
      ProjectPi3_3D = kx * ProjectP2(i,j,kx,ky,kz)
    else if(m .eq. 2)then
      ProjectPi3_3D = ky * ProjectP2(i,j,kx,ky,kz)
    else if(m .eq. 3)then
      ProjectPi3_3D = kz * ProjectP2(i,j,kx,ky,kz)
    else
      stop "ProjectPi3_3D error: m"
    endif
  end function ProjectPi3_3D
  !
end module pp
!+---------------------------------------------------------------------+
!| The end of the module pp                                            |
!+---------------------------------------------------------------------+