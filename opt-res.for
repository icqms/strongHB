c   170 high -> 1.2 GB
c   215         1.6
c   250         2.5
c   305 high -> 3.7

c       **** opts residues using ONIOM (B97D:PM3MM) one at a time ****
c       **** includes all residues within sphere of opt one 
c       ****   (chains terminated at neigbh CA, single skips forbidden)
c	**** HYD residues put with main, Mg, Ca, Fe kept together 
c	**** chains terminated with H ****
c
c	**** waters done first then chains then cofactors ****

      implicit real*8 (a-h,o-z)
      parameter (mo=10000)
      include 'readhin.cmn'
      character*80 fname

      parameter (mresar=200, mjob=4000, mjq=300,mproc=300)
      character*16 jobnam(mjob)
      integer iord(mresar),iind(mresar),njobs(0:5),njobf(0:5),
     $  jobcof(mmol),imem(mjob)
      logical incl(mres,mmol),lig,hbond
      character*3 ctyp(0:5)
      character*12 resnameuc
      integer imolar(mresar,mjob),iresar(mresar,mjob),nresar(mjob),
     $	      outord(mjob),qlist(mjq,mproc),nlist(mproc),ntot(mjob),nopt(mjob)
      logical incore(mresar)
      real*8 time(mjob),time1(mjob),tlist(mproc)
      common /maxp/ maxptch
      character*2 htype
C      parameter (mwcl=11)
C      integer wclm(mwcl),wclr(mwcl)
C      data wclm /283, 175,265,316,317,320, 211,219,318,  10, 10/
C      data wclr /  1,   1,  1,  1,  1,  1,  1,   1,  1, 119, 42/

      parameter (mwcl=0)
      integer wclm(max(1,mwcl)),wclr(max(1,mwcl))
C      data wclm /283, 175,265,316,317,320, 211,219,318,  10, 10/
C      data wclr /  1,   1,  1,  1,  1,  1,  1,   1,  1, 119, 42/

      common /imsave/ nmol1,natom1,nimage1
      common /chk/ ichkallow,ichkforb
      common /honly/ ifall

C      data rcutin2/9.0/, rcutou2/25.0/, rcutq2/100.0/
C      data rcutin2/9.0/, rcutou2/100.0/, rcutq2/100.0/
      data rcutin2/9.0/, rcutou2/64.0/, rcutq2/64.0/

      data ctyp /'wcl','cof','phy','evn','odd','hoh'/
      data njobs/0,0,0,0,0,0/, njobf/0,0,0,0,0,0/

      ichkallow= 0
      open (87,file='allowed.res',status='old',err=1)
      ichkallow= 1
1     continue

      ichkforb= 0
      open (88,file='forbidden.res',status='old',err=2)
      ichkforb= 1
2     continue

      write (6,*) 'Allowed and forbidden job file flags=',ichkallow,ichkforb

      open (5,file='opt-res.dat',status='old')
      open (4,file='opt-res.out',status='unknown')

c     **** input file name ****

10    write (6,*) 'Enter root name of .hin file'
      read (5,'(a)',err=500,end=500) fname
      if (fname(1:1).eq.'-') goto 500
      nf= 1
      do while (fname(nf:nf).ne.' ')
        nf= nf + 1
      end do
      fname(nf:nf+3)= '.hin'

      if (.not.readhin (fname,0)) goto 500
      rcutou2= min(rcutou2,rinscr2)
      rcutq2= min(rcutq2,rinscr2)
      
c     **** make image atoms manifest as real atoms ****

      call him2real (nmol1,natom1,nimage1)

c     **** test looking at Mg-ligand distances ****

      do imol= 1,nmol1
	i= molind(imol)
	if (nat(i).eq.12) then
	  do j= 1,natom
	    if (amol(j).ne.imol) then
	      r2= hbl2(i,j)
	      if (r2.lt.3.0) then
		write (6,*) 'MG coord error:'
		write (6,*) resname(1,imol),resname(ares(j),amol(j)),atname(j),sqrt(r2)
C		stop 'Mg coord'
	      end if
	    end if
	  end do
	end if
      end do

c     ***********

      write (6,*) 'enter ifall (-1 for local H only, 0 for H only, 1 for X only, 2 for lot)'
      read (5,*) ifall
      write (6,*) ' sel atoms=',ifall
      write (4,*) ' sel atoms=',ifall
      if (ifall.eq.-1) then
	rcutin2= 0.d0
	rcutou2= 0.d0
	rcutq2= 0.d0
      end if

      write (6,*) 'enter 0 to 5 for clusters, cof, phy, evn, odd, hoh'
      write (6,*) '-3, -4 for larger chain clusters, +3 , +4 for smaller'
      read (5,*) njtyp0
      njtyp= abs(njtyp0)
      write (6,*) 'calc type= ',ctyp(njtyp)

      write (6,*) 'enter nber of procs to use'
      read (5,*) nproc
      if (nproc.gt.mproc) stop 'MPROC'

      write (6,'(3(a,f8.1))') ' cutoff**2 internal=',rcutin2,' external=',rcutou2,' pt charge',rcutq2

c     **** set all residues as initially not optimized anywhere ****

      incl= .false.

c     **** clustered residues that need be kept together in joint opts ****

      if (mwcl.gt.0) then
      njob= 1
      njobs(0)= 1

150   continue
	read (5,*,end=160,err=160) nresar(njob),(iresar(i,njob),imolar(i,njob),i=1,nresar(njob))
	write (6,*) 'cluster:', nresar(njob),(iresar(i,njob),imolar(i,njob),i=1,nresar(njob))
	if (nresar(njob).gt.mresar) stop 'MRESAR'
	do i= 1,nresar(njob)
	  ires= iresar(i,njob)
	  imol= imolar(i,njob)
          incl(ires,imol)= .true.
	  if (i.eq.1) then
            call getjname (jobnam(njob),resname(ires,imol),imol)
	  else
	    write (6,810) resname(iresar(1,njob),imolar(1,njob)),resname(ires,imol)
	  end if
	end do
	njob= njob + 1
	goto 150

160   continue

c     **** add to cluster list the asymmetric water cluster at centre of PS1 trimer ****

c     **** first entry is non-rep molecule ****
      imol= wclm(1)
      iresar(1,njob)= 1
      imolar(1,njob)= imol
      incl(1,imol)= .true.
      call getjname (jobnam(njob),resname(1,imol),imol)
      k= 1
c     **** all images of residues in cluster ****
      do i= 2,mwcl
	do iim= 0,2*nmol1,nmol1
	  imol= wclm(i) + iim
	  ires= wclr(i)
	  k= k + 1
	  iresar(k,njob)= ires
	  imolar(k,njob)= imol
          incl(1,imol)= .true.
	  write (6,810) resname(1,wclm(1)),resname(ires,imol)
	end do
      end do
      nresar(njob)= k

      else
c	**** no cluster ****
        njob= 0
      end if

c     **** FS4 and CHL and first part of BCR cofactors ****

      njobf(0)= njob
      njobs(1)= njob+1
      write (6,*) 'looking for FS4, CHL, BCR'

      do imol= 1,nmol1
	resnameuc= resname(1,imol)
	call upcase (resnameuc)
	if (nres(imol).le.3 .and. nares(1,imol).ne.3 .or.
     $      resname(1,imol)(1:2).eq.'CL' .or. resname(1,imol)(1:2).eq.'BC' .or. resname(1,imol)(1:2).eq.'BP' ) then 
C	if (nres(imol).le.3 .and. nares(1,imol).ne.3 ) then
	  njob= njob + 1
	  write (6,*) imol,njob,nres(imol),nares(1,imol),resname(1,imol)
	  if (njob.gt.mjob) stop '1 MJOB'
	  jobcof(imol)= njob
	  nresar(njob)= 1
	  iresar(1,njob)= 1
	  imolar(1,njob)= imol
	  call getjname (jobnam(njob),resname(1,imol),imol)
	  incl(1,imol)= .true.

c	  **** add prev separated fragments around CHL ring ****
	  do ires= 2,nres(imol)
	    resnameuc= resname(ires,imol)
	    call upcase (resnameuc)
	    if (resnameuc(1:2).eq.'CL') then
	      i= nresar(njob) + 1
	      nresar(njob)= i
	      iresar(i,njob)= ires
	      imolar(i,njob)= imol
	      incl(ires,imol)= .true.
	    end if
	  end do

c	  **** search for internal ligands to Ca, Mg, Fe ****
	  i0= iares(1,imol) - 1
	  do i= iares(1,imol),jares(1,imol)
C	    if (nat(i).gt.8) then
	    if (nat(i).gt.12) then

	      do j= 1,natom
	       if (nat(j).gt.1) then
		jres= ares(j)
		jmol= amol(j)
		lig= .false.
		if (.not.incl(jres,jmol)) then
		  r2= hbl2(i,j)
		  lig= r2.lt.rcutin2 
C     $	.or.	    (resname(1,imol)(1:3).eq.'FS4' .and. r2.lt.rcutou2)
		end if

		if (lig) then
		  nresar(njob)= nresar(njob) + 1
		  if (nresar(njob).gt.mresar) stop 'MRESAR'
		  iresar(nresar(njob),njob)= jres
		  imolar(nresar(njob),njob)= jmol
		  incl(jres,jmol)= .true.
		  write (6,810) resname(1,imol),resname(jres,jmol)
810		  format (' to ',a,' is attached ',a)
		end if
	       end if
	      end do

	    end if
	  end do
	end if
      end do	

      njobf(1)= njob
      njobs(2)= njob+1

c     ***** carotene second part, phytyl chains, etc *****
c     ***** LHG 5003 is also incl as it is a ligand to CL1 1801 ****

      do imol= 1,nmol1
	do ires= 1,nres(imol)
	  if (resname(ires,imol)(1:3).eq.'PHY'  .or.
     $	      resname(ires,imol)(1:3).eq.'BCr'  .or.
     $	      resname(ires,imol)(1:3).eq.'PQn'  .or.
     $	      resname(ires,imol)(1:8).eq.'LHG 5003'  .or.
     $	      nres(imol).eq.3 .and. ires.gt.1) then
	    njob= njob + 1
	    if (njob.gt.mjob) stop '2 MJOB'
	    nresar(njob)= 1
	    iresar(1,njob)= ires
	    imolar(1,njob)= imol
	    call getjname (jobnam(njob),resname(ires,imol),imol)
	    incl(ires,imol)= .true.
	  end if
	end do
      end do

      njobf(2)= njob
      njobs(3)= njob+1

c     **** even residues ****

      do imol= 1,nmol1
       do ires= 2,nres(imol),2
	if (nares(ires,imol).gt.3 .and. .not.incl(ires,imol)) then
	  njob= njob + 1
	  if (njob.gt.mjob) stop '3 MJOB'
	  nresar(njob)= 1
	  iresar(1,njob)= ires
	  imolar(1,njob)= imol
	  call getjname (jobnam(njob),resname(ires,imol),imol)
	  incl(ires,imol)= .true.

c         **** look for h-bonded ions, opt all together ****

	  if (abs(qres(ires,imol)).gt.0.99) then
           kresar= 0
           do while (kresar.lt.nresar(njob))
            kresar= kresar + 1
            kmol= imolar(kresar,njob)
            kres= iresar(kresar,njob)
            do jmol= 1,nmol1
             do jres= 1,nres(jmol)
              if (.not.incl(jres,jmol) .and. abs(qres(jres,jmol)).gt.0.99) then
	       do jm= 1,nimage1
                do k= iares(kres,kmol),jares(kres,kmol)
                 do j= iares(jres,jmol),jares(jres,jmol)
                  if (hbond(k,j+(jm-1)*natom1)) then
                    nresar(njob)= nresar(njob) + 1
                    iresar(nresar(njob),njob)= jres
                    imolar(nresar(njob),njob)= jmol+(jm-1)*nmol1
                    incl(jres,jmol)= .true.
		    write (6,*) njob,' ',resname(ires,imol),' optimizing also ',resname(jres,jmol),' in image ',jm
                    goto 166
                  end if
		 end do
                end do
166             continue
               end do
              end if
	     end do
            end do
           end do
          end if

c	  **** search for HYD terminating residues ****

          do kresar= 1,nresar(njob)
            kmol= imolar(kresar,njob)
            kres= iresar(kresar,njob)
	    i0= iares(1,kmol) - 1
	    do i= iares(kres,kmol),jares(kres,kmol)
	      do jj= 1,ncon(i)
	        j= i0 + icon(jj,i)
		jres= ares(j)
		jmol= amol(j)
		if (resname(jres,jmol)(1:3).eq.'HYD') then
		  nresar(njob)= nresar(njob) + 1
		  if (nresar(njob).gt.mresar) stop 'MRESAR'
		  iresar(nresar(njob),njob)= jres
		  imolar(nresar(njob),njob)= jmol
		  incl(jres,jmol)= .true.
		  write (6,810) resname(kres,kmol),resname(jres,jmol)
		end if
	      end do
	    end do
	  end do

	end if
       end do
      end do	

      njobf(3)= njob
      njobs(4)= njob+1

c     **** odd residues ****

      do imol= 1,nmol1
       do ires= 1,nres(imol),2
	if (nares(ires,imol).gt.3 .and. .not.incl(ires,imol)) then
	  njob= njob + 1
	  if (njob.gt.mjob) stop '3 MJOB'
	  nresar(njob)= 1
	  iresar(1,njob)= ires
	  imolar(1,njob)= imol
	  call getjname (jobnam(njob),resname(ires,imol),imol)
	  incl(ires,imol)= .true.

c         **** look for h-bonded ions, opt all together ****

	  if (abs(qres(ires,imol)).gt.0.99) then
           kresar= 0
           do while (kresar.lt.nresar(njob))
            kresar= kresar + 1
            kmol= imolar(kresar,njob)
            kres= iresar(kresar,njob)
            do jmol= 1,nmol1
             do jres= 1,nres(jmol)
              if (.not.incl(jres,jmol) .and. abs(qres(jres,jmol)).gt.0.99) then
               do jm= 1,nimage1
                do k= iares(kres,kmol),jares(kres,kmol)
                 do j= iares(jres,jmol),jares(jres,jmol)
                  if (hbond(k,j+(jm-1)*natom1)) then
                    nresar(njob)= nresar(njob) + 1
                    iresar(nresar(njob),njob)= jres
                    imolar(nresar(njob),njob)= jmol+(jm-1)*nmol1
                    incl(jres,jmol)= .true.
                    write (6,*) njob,' ',resname(ires,imol),' optimizing also ',resname(jres,jmol),' in image ',jm
                    goto 165
                  end if
                 end do
                end do
165             continue
               end do
              end if
	     end do
            end do
           end do
          end if

c	  **** search for HYD terminating residues ****

          do kresar= 1,nresar(njob)
            kmol= imolar(kresar,njob)
            kres= iresar(kresar,njob)
	    i0= iares(1,kmol) - 1
	    do i= iares(kres,kmol),jares(kres,kmol)
	      do jj= 1,ncon(i)
	        j= i0 + icon(jj,i)
		jres= ares(j)
		jmol= amol(j)
		if (resname(jres,jmol)(1:3).eq.'HYD') then
		  nresar(njob)= nresar(njob) + 1
		  if (nresar(njob).gt.mresar) stop 'MRESAR'
		  iresar(nresar(njob),njob)= jres
		  imolar(nresar(njob),njob)= jmol
		  incl(jres,jmol)= .true.
		  write (6,810) resname(kres,kmol),resname(jres,jmol)
		end if
	      end do
	    end do
	  end do

	end if
       end do
      end do	

      njobf(4)= njob
      njobs(5)= njob+1

c     **** waters ****

      do imol= 1,nmol1
	if (nres(imol).eq.1 .and. resname(1,imol)(1:3).eq.'HOH' .and.
     $		.not.incl(1,imol)) then
	  njob= njob + 1
	  if (njob.gt.mjob) stop '5 MJOB'
	  nresar(njob)= 1
	  iresar(1,njob)= 1
	  imolar(1,njob)= imol
	  call getjname (jobnam(njob),resname(1,imol),imol)
	  incl(1,imol)= .true.

c	  **** look for h-bonded waters, opt all together ****
	  kresar= 0
	  do while (kresar.lt.nresar(njob))
	   kresar= kresar + 1
	   kmol= imolar(kresar,njob)
	   do jmol= 1,nmol1
	    if (.not.incl(1,jmol) .and. resname(1,jmol)(1:3).eq.'HOH') then
	      do jm= 1,nimage1
	        do k= molind(kmol),molend(kmol)
		  do j= molind(jmol),molend(jmol)
		    if (hbond(k,j+(jm-1)*natom1)) then
	  	      nresar(njob)= nresar(njob) + 1
	  	      iresar(nresar(njob),njob)= 1
	  	      imolar(nresar(njob),njob)= jmol+(jm-1)*nmol1
	  	      incl(1,jmol)= .true.
		      write (6,*) resname(1,imol),' optimizing also ',resname(1,jmol),' in image ',jm
		      goto 167
		    end if
		  end do  
		end do
167	      continue
	      end do
	    end if
	   end do
	  end do
	end if
      end do	

      njobf(5)= njob

      do i= 0,5
        write (6,*) 'nber of ',ctyp(i),' calcs=',njobf(i)-njobs(i)+1
      end do

c     ***** generate job files ******

      maxsze= 0
      maxhigh= 0
      maxptch= 0

      write (6,*) 'njtyp=',njtyp,njtyp0
      do i= njobs(njtyp),njobf(njtyp)
	if (njtyp.ne.3 .and. njtyp.ne.4 .or. (njtyp0.gt.0 .xor. nresar(i).gt.1)) then
	write (6,*) 'calling job',i,nresar(i),' ',jobnam(i)
        call job (iresar(1,i),imolar(1,i),incore,mresar,nresar(i),
     $	  nnn,nhigh,nopt(i),ifall,jobnam(i),iord,iind,rcutou2,rcutq2,
     $    jobcof,imolar,iresar,nresar)
	imem(i)= nhigh
	ntot(i)= nnn
	time(i)= nhigh**3
	time1(i)= time(i)
	write (6,*) 'job time=',nnn,time(i)/1.e6
	if (nnn.eq.0) write (6,*) 'JOB ABORTED: ',jobnam(i),nnn
	maxsze= max (maxsze,nnn)
	maxhigh= max (maxhigh,nhigh)
	end if
      end do

c     ***** batch scripts *****

      do i= 0,5
        write (6,*) 'nber of ',ctyp(i),' calcs=',njobf(i)-njobs(i)+1
      end do

      write (6,*) 'BATCH JOBS'

      ityp= njtyp
C      do ityp= 0,5
	write (6,*) 'BATCH JOB',ityp,njobs(ityp),njobf(ityp)
	write (fname,'(3a)') 'opt/run-',ctyp(ityp),'.csh'
        open (8,file=fname,status='unknown')
	write (fname,'(3a)') 'opt/reinsert-',ctyp(ityp),'.dat'
        open (3,file=fname,status='unknown')

c	**** order jobs by est time, biggest first ****

	njobout= 0
	do njob= njobs(ityp),njobf(ityp)
	  tmax= -2.d0
	  do kjob= njobs(ityp),njobf(ityp)
	    if (time1(kjob).gt.tmax) then
	      tmax= time1(kjob)
	      kmax= kjob
	    end if
	  end do
	  time1(kmax)= -1.0
	  if (tmax.gt.0.d0) then
	    njobout= njobout + 1
	    outord(njobout)= kmax
C	    write (6,*) 'job order',kmax,tmax
	  end if
	end do

c	**** assign jobs to procs ****

	do iproc= 1,nproc
	  nlist(i)= 0
	  tlist(i)= 0.d0
	end do

	do ijobout= 1,njobout
	  njob= outord(ijobout)

c	  **** find proc with lowest load ****
	  procmin= 1.e30
	  do iproc= 1,nproc
C	    write (6,*) 'tlist=',tlist(iproc),ijobout,njob
	    if (tlist(iproc).lt.procmin) then
	      procmin= tlist(iproc)
	      jproc= iproc
	    end if
	  end do

C	  write (6,*) 'assign ',njob,' to proc',jproc
	  nlist(jproc)= nlist(jproc) + 1
	  tlist(jproc)= tlist(jproc) + time(njob)
	  qlist(nlist(jproc),jproc)= njob
	end do

c	**** output jobs, slowest first ****

	do iproc= 1,nproc
C	  write (6,*) 'proc=',iproc,' jobs=',nlist(iproc)
	  if (nlist(iproc).gt.0) then
	    write (fname,'(a,i3.3)') ctyp(ityp),iproc
C	    write (6,*) 'opening ',fname
	    open (7,file='opt/'//fname,status='unknown')
	    imem1= 0
	    do ijob= 1,nlist(iproc)
	      njob= qlist(ijob,iproc)
	      imem1= max (imem1,imem(njob))
	    end do
	    if (imem1.le.140) then
	      write (7,840) 1000
	    else if (imem1.le.200) then
	      write (7,840) 1600
	    else if (imem1.le.250) then
	      write (7,840) 2800
	    else if (imem1.le.310) then
	      write (7,840) 4000
	    else 
	      write (7,840) 5000
	    end if
	    write (8,'(5hqsub ,a)') fname
	    do ijob= 1,nlist(iproc),1
	      njob= qlist(ijob,iproc)
C	      write (6,'(a,i5,f8.0,i4)') ' ass=',njob,time(njob),iproc
C	      if (nopt(njob).gt.3) write (7,'(3hsh ,a)') jobnam(njob)
	      write (7,'(3hsh ,a)') jobnam(njob)
	      k= 1
	      do while (jobnam(njob)(k:k).ne.'.')
	        k= k + 1
	      end do
C	      if (nopt(njob).gt.3) write (3,'(a)') jobnam(njob)(1:k-1)
	      write (3,'(a)') jobnam(njob)(1:k-1)
	    end do
	    close (7)
	  end if
	end do

840	format (
     $   '#!/bin/bash' /
     $   '#PBS -l walltime=24:00:00' /
     $   '#PBS -l vmem=3000MB' /
     $   '#PBS -q normal' /
     $   '#PBS -l ncpus=4' /
     $   '#PBS -l jobfs=',i4.4,'MB' /
     $   '#PBS -l software=g09' /
     $   '#PBS -wd' //
     $   'module load gaussian' /)

	close (8)
C      end do

      write (6,*) 'nber jobs actually generated=',njobout
      write (4,*) 'nber jobs actually generated=',njobout
      write (6,*) 'max nber atoms in any calc=',maxsze
      write (4,*) 'max nber atoms in any calc=',maxsze
      write (6,*) 'max nber high-level atoms =',maxhigh
      write (4,*) 'max nber high-level atoms =',maxhigh
      write (6,*) 'max nber point charges    =',maxptch
      write (4,*) 'max nber point charges    =',maxptch

500   continue

      end

c     ***************************************************************

      subroutine job (ires0,imol0,incore,mresar,nres0,natsel,nathigh,nopt,ifall,
     $      fname,iord,iind,rcutou2,rcutq2, jobcof,imolar,iresar,nresar)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      character*(*) fname
      character*80 fname1
      character*20 uffat
      character*8 basis(2),cqlink
      character*5 raname
      character*3 flags,res
      character*2 attypeo
      character*1 optfl,optfladj(m)
      integer imol0(mresar),ires0(mresar),iord(mresar),iind(mresar),
     $	      jobcof(mmol)
      integer imolar(mresar,*),iresar(mresar,*),nresar(*)
      logical incore(mresar),hasc

      common /maxp/ maxptch
      common /imsave/ nmol1,natom1,nimage1

      parameter (mo=3000, manion=50)
      common /outp/ qo(m),xo(mo),yo(mo),zo(mo),nato(mo),nnn,freez(mo),
     $  flags(m),hasfs4,anions(manion),aniono(manion),nanions,naniono,
     $  attypeo(mo),nnntoatnumb(m),ibasis(mo)

      parameter (msave=100)
      character*5 atns
      common /save/ xs(msave),ys(msave),zs(msave),
     $		nats(msave),ihvs(msave), ns, atns(msave)

      logical repl(m),afreez(m),freez,incl(mres,mmol),shortdis,
     $		isfs4,hasfs4,ifuff,isptchg,isaug, hbond
      integer anions,aniono
      character*1 opt
      character*12 resnameuc
      character*2 htype
      common /chk/ ichkallow,ichkforb

      character*8 opttyp(2)
      data opttyp /'modredun','z-matrix'/

      data basis /'6-31G*','6-31+G*'/

c     *************
      natsel= 0
      nathigh= 0
      write (6,*) 'Allowed and forbidden job file flags=',ichkallow,ichkforb

      isptchg= rcutou2.lt.rcutq2
      write (6,*) 'Use point charges flag= ',isptchg

c     **** only proceed if its in the allowed list ****
      if (ichkallow.eq.0) goto 1015
      rewind 87
1010  read (87,'(a)',end=1011) fname1
	 l= len(fname)
	 if (fname.eq.fname1(1:l)) goto 1015
	goto 1010
1011  return

1015  continue
      write (6,*) 'allowing fname'

c     **** only proceed if it isnt in the forbidden list ****
      if (ichkforb.eq.0) goto 1020
      rewind 88
1016  read (88,'(a)',end=1020) fname1
	 l= len(fname)
C	 if (fname.eq.fname1(1:l)) write (6,*) 'excluding ',fname
	 if (fname.eq.fname1(1:l)) return
	 goto 1016

1020  continue


      write (6,*) 'Starting output for job ',fname
      hasfs4= .false.
      nanions= 0
      naniono= 0

c	**** indices ****

	do imol= 1,nmol
	  do ires= 1,nres(imol)
	    incl(ires,imol)= .false.
	  end do
	end do

	do i= 1,natom
	  afreez(i)= .false.
	  repl(i)= .false.
	end do

	call hdeselect1 (1,natom1)

c	**** write G98 opt file and selection info to regen .hin afterwards ****

	nf= 1
	do while (fname(nf:nf).ne.' ')
	  nf= nf + 1
	end do

	fname(nf:nf+3)= '.sel'
	open (1,file='opt/'//fname,status='unknown')
	fname(nf:nf+3)= '.com'
	open (2,file='opt/'//fname,status='unknown')

	nnn= 0
	nqtot= 0

c	********* output atoms in selected residues *************

	do kinc= 1,nres0
	  incore(kinc)= .true.
	  ires= ires0(kinc)
	  imol= imol0(kinc)
	  incl(ires,imol)= .true.
	  ist= iares(ires,imol)
	  ifn= jares(ires,imol)
	  call hselres1 (ires,imol,1)
	  nqtot= nqtot + nint (qres(ires,imol))

	  do i= ist,ifn
	    optfl= '-'
	    if (nat(i).gt.0) then
	      afreez(i)= ifall.eq.0 .and. nat(i).gt.1 .and. bfact(i).gt.0.d0 .or.
     $	                 ifall.eq.1 .and. nat(i).le.1 
	      optfl= 'o'
	      if (afreez(i)) optfl= 'F'
	      optfladj(i)= optfl
	    end if
	  end do
	end do

	write (6,*) 'Initial charge on selected residues=',nqtot

c	*********** find residues nearby to these ****************

	nh2onearby= 0
	nrest= nres0

	do kinc= 1,nres0
	  write (6,*) 'Searching around ',ires0(kinc),imol0(kinc)
	  ist= iares(ires0(kinc),imol0(kinc))
	  ifn= jares(ires0(kinc),imol0(kinc))
	
	  do imol= 1,nmol
	   do ires= 1,nres(imol)
	    if (.not.incl(ires,imol)) then
	     i0= molind(imol) - 1
	     do i= iares(ires,imol),jares(ires,imol)
		do j= ist,ifn
C		 if (.not. afreez(j) .or. nat(i).gt.6) then
		 if (.true.) then
		  if ( shortdis(i,j,rcutou2) ) then
		    if (resname(ires,imol)(1:3).eq.'HOH')
     $			nh2onearby= 1 + nh2onearby
		    call hselres1 (ires,imol,1)
			if (resname(ires,imol)(1:2).eq.'PQ') write (6,*) 'has PQN ',fname
		    incl(ires,imol)= .true.
		    nrest= nrest + 1
		    if (nrest.gt.mresar) stop 'MRESAR'
		    ires0(nrest)= ires
		    imol0(nrest)= imol
		    incore(nrest)= .false.

c		    **** if this is HYD then make sure its primary res is included ****
		    if (resname(ires,imol)(1:3).eq.'HYD' ) then
		      jres= ares(i0+icon(1,i))
		      jmol= amol(i0+icon(1,i))
		      if (.not.incl(jres,jmol)) then
		        call hselres1 (jres,jmol,1)
		        incl(jres,jmol)= .true.
		        nrest= nrest + 1
		        if (nrest.gt.mresar) stop 'MRESAR'
		        ires0(nrest)= jres
		        imol0(nrest)= jmol
		        incore(nrest)= .false.
		      end if
		    end if

		    goto 300
		  end if
		 end if
	 	end do
	     end do

300	     continue
	    end if
	   end do
	  end do

	  if (kinc.eq.1) nres1= nrest

	end do

c	***** make sure HYD residues are included to terminate broken chains ***

	nrest0= nrest
	do kinc= 1,nrest0
	  ires= ires0(kinc)
	  imol= imol0(kinc)
	  i0= molind(imol) - 1
	  do j= iares(ires,imol),jares(ires,imol)
	    do kk= 1,ncon(j)
	      k= i0 + icon(kk,j)
	      kres= ares(k)
	      kmol= amol(k)
	      if (.not.incl(kres,kmol) .and.
     $		    resname(kres,kmol)(1:3).eq.'HYD') then
		nrest= nrest + 1
		if (nrest.gt.mresar) stop 'MRESAR'
		ires0(nrest)= kres
		imol0(nrest)= kmol
		incore(nrest)= incore(kinc)
		if (.not.incore(nrest0)) write (6,*) 'HYD not in core ',fname
		call hselres1 (kres,kmol,1)
		incl(kres,kmol)= .true.
	      end if
	    end do
	  end do
	end do

c	******* put all res assoc with hyd bonding, S-S bonding, and all cov bonds in cofactors into core ****

	do kinc= 1+nres0,nres1
	 if (.not.incore(kinc)) then
	  kres= ires0(kinc)
	  kmol= imol0(kinc)
	  k0= molind(kmol) - 1
	  res= resname(kres,kmol)(1:3)
	  incore(kinc)= .true.
C	  if (res.eq.'HOH') goto 346
	  do k= iares(kres,kmol),jares(kres,kmol)
	    do iinc= 1,nres0
	      ires= ires0(iinc)
	      imol= imol0(iinc)
	      do i= iares(ires,imol),jares(ires,imol)
c		**** check for h-bonds ****
	        if (hbond(i,k)) goto 346
c	        **** check for covalent bonds between residues ****
C	  	if (nres(kmol).le.10 .and. kres.ne.10) then
	  	if (nres(kmol).le.10) then
c		  *** add bonded cofactor residues except PHY ****
		  do ll= 1,ncon(k)
		    l= icon(ll,k) + k0
		    if (l.eq.i) goto 346
		  end do
		end if
	      end do
	    end do
	  end do 
	  incore(kinc)= .false.
346	  continue

	  if (incore(kinc)) then
	    write (6,*) 'Putting into core H-bonded res ',
     $		resname(kres,kmol)
	  end if

	 end if
	end do

c	******** if BC:L or BPH, put CLi into core if central ring is also ****

	do ii= 1,nrest
	 if (incore(ii)) then
	  ires= ires0(ii)
	  imol= imol0(ii)
	  if (resname(ires,imol)(1:3).eq.'BP1' .or. resname(ires,imol)(1:2).eq.'BC') then
	    do jres= 1,nres(imol)
	      if (resname(jres,imol)(1:3).eq.'CLi') then
c			**** see if already selected ****
			do kinc= 1,nrest
			  if (ires0(kinc).eq.jres .and. imol0(kinc).eq.imol) then
			    if (.not.incore(kinc)) then
			      write (6,*) 'Putting into core conjugated res ',resname(jres,imol)
			      incore(kinc)= .true.
			    end if
			    goto 314
			  end if
			end do
c			**** select CLi ****
	    		nrest= nrest + 1
	    		if (nrest.gt.mresar) stop 'MRESAR'
	    		ires0(nrest)= jres
	    		imol0(nrest)= imol
	    		incl(ires0(nrest),imol0(nrest))= .true.
	    		call hselres1 (ires0(nrest),imol0(nrest),1)
			write (6,*) 'Selecting and uptting into core sonjugated res ',resname(jres,imol)
	    		incore(nrest)= .true.
			goto 314

	      end if
	    end do
314	    continue
	  end if
	 end if
	end do


c	******** if FS4 is in core, make sure all CyS's are too ' ****

	isfs4= .false.

	do ii= 1,nrest
	 if (incore(ii)) then
	  ires= ires0(ii)
	  imol= imol0(ii)
	  if (resname(ires,imol)(1:3).eq.'FS4') then
            isfs4= .true.
	    write (6,*) 'Searching for CyS bonded to FS4: ', resname(ires,imol)
	    do jmol= 1,nmol
	      do jres= 1,nres(jmol)
		if (resname(jres,jmol)(1:3).eq.'CyS') then
	          do i= iares(ires,imol),jares(ires,imol)
		    do j= iares(jres,jmol),jares(jres,jmol)
		      if (hbl2(i,j).lt.7.0) then
c			**** see if CyS already selected ****
			do kinc= 1,nrest
			  if (ires0(kinc).eq.jres .and. imol0(kinc).eq.jmol) then
			    if (.not.incore(kinc)) then
			      write (6,*) 'Putting into core FS4-bonded res ',resname(jres,jmol)
			      incore(kinc)= .true.
			    end if
			    goto 324
			  end if
			end do
c			**** select CyS ****
	    		nrest= nrest + 1
	    		if (nrest.gt.mresar) stop 'MRESAR'
	    		ires0(nrest)= jres
	    		imol0(nrest)= jmol
	    		incl(ires0(nrest),imol0(nrest))= .true.
	    		call hselres1 (ires0(nrest),imol0(nrest),1)
			write (6,*) 'Putting into core FS4-bonded res ',resname(jres,jmol)
	    		incore(nrest)= .true.
			goto 324
		      end if
                    end do
		  end do
		end if
324		continue
	      end do
	    end do
	  end if
	 end if
	end do

c	***** make sure no single residue is skipped *****

c	**** get residues in order ****
	do kinc= 1,nrest
	  iind(kinc)= ires0(kinc)+imol0(kinc)*mres
	end do
	do i= 1,nrest
	  imin= 1000000000
	  do j= 1,nrest
	    if (imin.gt.iind(j)) then
	      imin= iind(j)
	      jmin= j
	    end if
	  end do
	  iord(i)= jmin
	  iind(jmin)= 1000000001
	  write (6,*) 'order',i,iord(i),' ',resname(ires0(jmin),imol0(jmin))
	end do

	nrest0= nrest

c	**** look for single miss ****
	write (6,*) 'Looking for single residue misses'
	do ii= 2,nrest0
	  i= iord(ii)
	  j= iord(ii-1)
	  resnameuc= resname(ires0(j),imol0(j))
	  call upcase (resnameuc)
	  if (resnameuc(1:2).ne.'CL' .and. resnameuc(1:2).ne.'BC' .and. resnameuc(1:2).ne.'BP' .and. 
     $        imol0(j).eq.imol0(i) .and. ires0(j).eq.ires0(i)-2) then
	    nrest= nrest + 1
	    if (nrest.gt.mresar) stop 'MRESAR'
	    ires0(nrest)= ires0(i) - 1
	    imol0(nrest)= imol0(i)
	    incore(nrest)= .false.
	    incl(ires0(nrest),imol0(nrest))= .true.
	    call hselres1 (ires0(nrest),imol0(nrest),1)
	  end if
	end do
 
c        **** look for single miss in the core ****
        write (6,*) 'Looking for single residue misses in the core'
        do ii= 3,nrest0
           i= iord(ii)
	   ij1= iord(ii-1)
           j= iord(ii-2)
           resnameuc= resname(ires0(j),imol0(j))
           call upcase (resnameuc)
	   if (resnameuc(1:2).ne.'CL' .and. resnameuc(1:2).ne.'BC' .and. resnameuc(1:2).ne.'BP' .and. 
     $         imol0(j).eq.imol0(i) .and. ires0(j).eq.ires0(i)-2.and.incore(i).and.incore(j)
     $          .and.(.not.incore(ij1))) then
               incore(ij1)=.true.     
               call hselres1 (ires0(i)-1,imol0(i),1)
           end if
        end do
     
c       **** look for double misses in the core ****
        write (6,*) 'Looking for double residue misses in the core'
        do ii= 4,nrest0
           i= iord(ii)
	   ij1= iord(ii-1)
	   ij2= iord(ii-2)
           j= iord(ii-3)
           resnameuc= resname(ires0(j),imol0(j))
           call upcase (resnameuc)
	   if (resnameuc(1:2).ne.'CL' .and. resnameuc(1:2).ne.'BC' .and. resnameuc(1:2).ne.'BP' .and. 
     $         imol0(j).eq.imol0(i) .and. ires0(j).eq.ires0(i)-3 .and. incore(i).and. incore(j)
     $         .and. (.not.incore(ij1)) .and. (.not.incore(ij2))) then
              incore(ij1)=.true.
              incore(ij2)=.true.
              call hselres1 (ires0(i)-2,imol0(i),1)
              call hselres1 (ires0(i)-1,imol0(i),1)
           end if
        end do 
  
c       **** look for double misses, eg if 15 and 18 are selected include 16 and 17 as well ****
c        write (6,*) 'Looking for double residue misses'
c        do ii = 2,nrest
c           i = iord(ii)
c           j = iord(ii-1)
c           resnameuc= resname(ires0(j),imol0(j))
c           call upcase (resnameuc)
c	   if (resnameuc(1:2).ne.'CL' .and. resnameuc(1:2).ne.'BC' .and. resnameuc(1:2).ne.'BP' .and. 
c     $		imol0(j).eq.imol0(i) .and. ires0(j).eq.ires0(i)-3) then
c                nrest= nrest + 1
c                if (nrest.gt.mresar) stop 'MRESAR'
c                ires0(nrest)= ires0(i) - 2 
c                imol0(nrest)= imol0(i)
c                incore(nrest)=.false.
c                incl(ires0(nrest),imol0(nrest))=.true.
c                call hselres1 (ires0(nrest),imol0(nrest),1)
               
c                nrest=nrest + 1
c                if (nrest.gt.mresar) stop 'MRESAR'
c                ires0(nrest)= ires0(i) - 1
c                imol0(nrest)= imol0(i)
c                incore(nrest)=.false.
c                incl(ires0(nrest),imol0(nrest))=.true.
c                call hselres1 (ires0(nrest),imol0(nrest),1)            
  
c           end if
c        end do

c	******** put new atoms into output list *********

	nqtoth= nqtot
	do kinc= 1+nres0,nrest
	  kres= ires0(kinc)
	  kmol= imol0(kinc)
	  nqtot= nqtot + nint(qres(kres,kmol))
	  if (incore(kinc)) nqtoth= nqtoth + nint(qres(kres,kmol))
	  do k= iares(kres,kmol),jares(kres,kmol)
	    optfladj(k)= 'f'
	    if (incore(kinc)) optfladj(k)= 'F'
	  end do
	end do
          
c	******** mark link atoms for oniom ********

      write (6,*) 'marking oniom link atoms'
	do ii= 1,nrest
	 if (incore(ii)) then
	  ires= ires0(ii)
	  imol= imol0(ii)
	  i0= molind (imol) - 1
	  do i= iares(ires,imol),jares(ires,imol)
c	    **** i is an atom in the core part ****

	    do kinc= 1+nres0,nrest
	     if (.not.incore(kinc)) then
	      kres= ires0(kinc)
	      kmol= imol0(kinc)
	      do k= iares(kres,kmol),jares(kres,kmol)
c	        **** k is an atom in the non-core part ****

	        do jj= 1,ncon(i)
	          j= i0 + icon(jj,i)
	          if (j.eq.k) then
c		    **** k is connected to the core

c		    **** find link atoms, put other con atoms into high level ****

		   if (atname(k).eq.'  CGD') then
c		    **** patch for CLh residue ****
		    optfladj(k)= 'L'
		    write (6,*) 'set direct link for CLh in ', resname(kres,kmol)

		   else

		    do ll= 1,ncon(k)
		      l= i0 + icon(ll,k)
		      if (l.ne.i) then

			if (atname(l).eq.'   CA') then

c			  **** protein chain, add H and make links ****
			  optfladj(l)= 'F'

c			  **** atoms conn to CA ****
			  do iji= 1,ncon(l)
			    ij= i0 + icon(iji,l)

			    if (nat(ij).eq.1) then
			      optfladj(ij)= 'F'
			    else
			      if ((atname(ij).eq.'    C' .or. atname(ij).eq.'    N')
     $					.and. optfladj(ij).eq.'L') then
c				**** watch out for this res already been made jn from other end ****
				optfladj(ij)= 'F'
			        write (6,844) ij,atname(ij),resname(kres,kmol)
			      else if (optfladj(ij).eq.'f') then
				optfladj(ij)= 'L'
			        write (6,843) ij,atname(ij),resname(kres,kmol)
			      end if
			    end if
			  end do

			else if (nat(l).eq.6) then 
			  if (scon(ll,k).eq.'s') then
			    optfladj(l)= 'L'
			    write (6,843) l,atname(l),resname(kres,kmol)
			  else
c			    **** break at double bond ****
			    optfladj(l)= 'F'

c			    **** atoms conn to end of double bond ****
			    do iji= 1,ncon(l)
			      ij= i0 + icon(iji,l)
			      if (nat(ij).eq.1) then
			        optfladj(ij)= 'F'
			      else
			        optfladj(ij)= 'L'
				write (6,843) ij,atname(ij),resname(kres,kmol)
			      end if
			    end do
			  end if

			else
c			  *** other connected atom (eg, H or O) ****
			  optfladj(l)= 'F'
			end if

843			format (' set  LINK atom ',i6,1x,a,1x,a)
844			format (' DBLE LINK atom ',i6,1x,a,1x,a)

		      end if
		    end do

c		    *** check if atom k needs to be opt ****

     		    optfladj(k)= 'F'
		    if (ifall.ne.1 .and. atname(k).eq.'HXT')
     $			optfladj(k)= 'o'
		    write (6,*) 'next res con atom= ',optfladj(k),k,
     $				atname(k)

		   end if

		  end if

		end do
	      end do
	     end if
	    end do
	  end do
	 end if
	end do

c	**** actually write entire list ****

	do kinc= 1,nrest
	  kres= ires0(kinc)
	  kmol= imol0(kinc)
	  do k= iares(kres,kmol),jares(kres,kmol)
c	    **** k is an atom in the non-core part ****
	    if (nat(k).gt.0) call tooutl (k,optfladj(k))
	  end do
	end do

c	******** search for bonds to not included residues *****

	ns= 0

	do imol= 1,nmol
	 do ires= 1,nres(imol)
	  if (incl(ires,imol)) then
	    i0= molind(imol) - 1
C	    write (6,*) 'checking links to',resname(ires,imol)
	    do i= iares(ires,imol),jares(ires,imol)

	      do jj= 1,ncon(i)
		j= i0 + icon(jj,i)
		jres= ares(j)
		jmol= amol(j)
		j1= iares(jres,jmol) - 1
C		write (6,*) atname(i),' ',atname(j),resname(jres,jmol)
		if (.not.incl(jres,jmol)) then

	          optfl= optfladj(i)
	          if (optfl.eq.'o') optfl= 'F'
	          if (optfl.eq.'L') write (6,*) 'TERMINAL on LINK'
	          if (optfl.eq.'L') optfl= 'f'
		  write (6,*) 'replacing link to ',atname(j),
     $			  ' ',resname(jres,jmol),' flag= ',optfl

	         resnameuc= resname(jres,jmol)
	         call upcase (resnameuc)
		 if (resname(jres,jmol)(1:3).eq.'cYS' .or.
     $		     resnameuc(1:2).eq.'CL'  .or.
     $		     resnameuc(1:2).eq.'BC'  .or.
     $		     resnameuc(1:2).eq.'BP'  .or.
     $		     resnameuc(1:2).eq.'LM' .or.
     $		     resnameuc(1:2).eq.'LH' .or.
     $		     resnameuc(1:2).eq.'PQ' .or.
     $		     resnameuc(1:3).eq.'PHY') then

c		   **** S-S and CHL-PHY jn, just break bond and put in H ****
C		   write (6,*) 'by direct H insertion'
		   call brbond (i,j,repl)
	           call tooutl (j,optfl)

		 else

c		  **** j is chain atype N or C, find CA next to it ****
		  ica= 0
		  do kk= 1,ncon(j)
		    k= i0 + icon(kk,j)
		    if (atname(k).eq.'   CA') ica= k
		    if (k.ne.i) then
c		      **** k is either H, O, CA or CD for PRO ****
		      if (atname(k).eq.'   CD') call brbond (j,k,repl)
		      call tooutl (k,optfl)
		    end if
		  end do
		  if (ica.eq.0) stop 'CA not found'

c		  **** put in atoms connected to CA (incl atom j) ****
		  ihv= 0
		  do kk= 1,ncon(ica)
		    ihv= i0 + icon(kk,ica)
		    if (nat(ihv).gt.1 .and. ihv.ne.j)
     $			call brbond (ica,ihv,repl)
		    call tooutl (ihv,optfl)
		  end do

		 end if
		end if
	      end do
	    end do
	  end if
	 end do
	end do

c	******** add H to any included CyS if FS4 not in core ****

	do ii= 1,nrest
	  ires= ires0(ii)
	  imol= imol0(ii)
	  if (.not.incore(ii) .and. resname(ires,imol)(1:3).eq.'CyS') then
	    i0= molind (imol) - 1
	    is= iares(ires,imol) + 5
c	    **** look for connected Fe ****
	    jres= 1
	    jfe= 0
	    do jmol= 1,nmol
	      if (resname(jres,jmol)(1:3).eq.'FS4') then
		do j= iares(jres,jmol),jares(jres,jmol)
		  if (nat(j).eq.26 .and. hbl2(is,j).lt.6.0) jfe= j
		end do
	      end if
	    end do
	    if (jfe.eq.0) stop 'CyS not near FS4 ?'
	    call brbond (is,jfe,repl)
	    call tooutl (jfe,optfladj(is))
	    write (6,*) 'Adding H to replace FS4 next to non-core CyS'
	    nqtot= nqtot + 1
	    if (optfladj(is).eq.'F') nqtoth= nqtoth + 1
	  end if
	end do

c	******** count nber of frozen atoms and nber electrons *****

	ifkill= 0
	nfreez= 0
	nberel= 0
	do i= 1,nnn
	  if (freez(i)) nfreez= nfreez + 1
	  nberel= nberel + nato(i)
	end do
	nopt= nnn - nfreez
	nberel= nberel - nqtot
	multip= mod(nberel,2) + 1

c	**** look at sample size to see if plane-wave code useable ****
c	**** get partial ONIOM system info ****

	nathigh= 0
	natlink= 0
	xmin= 1.e30
	ymin= 1.e30
	zmin= 1.e30
	xmax= - xmin
	ymax= - xmin
	zmax= - xmin
	do i= 1,nnn
	  if (flags(i).eq.'   ') then
	    nathigh= nathigh+1
	    if (nato(i).eq.6) hasc= .true.
	  end if
	  if (flags(i).eq.'M L') natlink= natlink+1
	  xmin= min (xmin,xo(i))
	  ymin= min (ymin,yo(i))
	  zmin= min (zmin,zo(i))
	  xmax= max (xmax,xo(i))
	  ymax= max (ymax,yo(i))
	  zmax= max (zmax,zo(i))
	end do

	do i= 4,6,2
	write (i,*) 'total  number of atoms=',nnn
	write (i,*) 'frozen number of atoms=',nfreez
	write (i,*) 'opt    number of atoms=',nopt
	write (i,*) 'high   number of atoms=',nathigh
	write (i,*) 'link   number of atoms=',natlink
	write (i,*) 'High-level charge=',nqtoth
	write (i,*) 'total charge=',nqtot
	write (i,*) 'total nberel=',nberel
	write (i,*) 'multiplicity=',multip
	write (i,870) xmax-xmin,ymax-ymin,zmax-zmin,
     $    (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
	end do

870	format (' xyz range=',3f8.3,' vol=',f8.0)


	if (multip.eq.2) stop 'DOUBLET STATE ???'

c	**** nnn=1600 runs in 2400MB RAM, 1800 3000MB ****
	iopttyp= 1
	if (isptchg .or. nnn.gt.1800) iopttyp= 2
	iopttyp= 2

c	**** G03 job header ****

	write (2,900) fname(1:nf-1),fname(1:nf-1)
	write (6,*) 'isfs4,hasfs4=',isfs4,hasfs4

	if (nnn.eq.nathigh) then
c	  **** no mdedium level, all DFT ****
	  write (6,*) 'Full DFT calc'
	  write (2,913)
C	  if (nanions.gt.0 .or. naniono.gt.0) write (2,915)
	  if (isptchg) write (2,916)
     	  write (2,920) opttyp(iopttyp),fname(1:nf-1)
	  write (2,925) nqtot,multip, nqtoth,1, nqtoth,1
	  ifuff= .false.

	else if (isfs4) then
c	  **** FS4 and its 4 CyS's in high level, use UFF not PM3 '****
	  write (2,911)
C	  if (nanions.gt.0 .or. naniono.gt.0) write (2,915)
	  if (isptchg) write (2,916)
     	  write (2,920) opttyp(iopttyp),fname(1:nf-1)
	  write (2,925) nqtot,multip, nqtoth,1, nqtoth,1
	  ifuff= .true.

	else if (hasfs4) then
c	  **** opt of nearby residue to FS4, put FS4 into 3rd layer ****
	  write (2,912)
C	  if (nanions.gt.0 .or. naniono.gt.0) write (2,915)
	  if (isptchg) write (2,916)
     	  write (2,920) opttyp(iopttyp),fname(1:nf-1)
	  write (2,925) nqtot,multip, nqtot-2,1, nqtot-2,1,
     $			nqtoth,1, nqtoth,1, nqtoth,1
	  ifuff= .true.

	else
c	  **** nothing to do with FS4 ****
	  write (2,910)
C	  if (nanions.gt.0 .or. naniono.gt.0) write (2,915)
	  if (isptchg) write (2,916)
     	  write (2,920) opttyp(iopttyp),fname(1:nf-1)
	  write (2,925) nqtot,multip, nqtoth,1, nqtoth,1
	  ifuff= .false.
	end if

900	format (
     $    'g09 << +++ > ',a,'.log' /
     $    '%mem=1900MB' /
     $    '%Nprocshared=4' /
     $    '%Chk=/short/d63/reimers/',a )

C910	format ('#P ONIOM(B971/6-31G*/Auto:PM3) iop(4/28=5)' )
CC911	format ('#P ONIOM(B971/6-31G*/Auto:UFF)=embed' )
CC912	format ('#P ONIOM(B971/6-31G*/Auto:PM3:UFF)=embed iop(4/28=5)' )
C911	format ('#P ONIOM(B971/6-31G*/Auto:UFF)' )
C912	format ('#P ONIOM(B971/6-31G*/Auto:AMBER:UFF) iop(4/28=5)' )
910	format ('#P ONIOM(cam-b3lyp:AMBER=hardfirst)=EmbedCharge iop(4/28=5)' )
C910	format ('#P ONIOM(B971/GenECP:AMBER=hardfirst)=EmbedCharge iop(4/28=5)' )
911	format ('#P ONIOM(B971/6-31G*/Auto:UFF)' )
912	format ('#P ONIOM(B971/6-31G*/Auto:AMBER:UFF) iop(4/28=5)' )
913	format ('#P B971/GenECP')
915	format ('   extrabasis')
916	format ('   charge')
920	format ('   opt(',a,',maxcycle=7,nrscale) nosym scf(conver=6) ' //
     $		'  lysozyme partial opt ',a / )
925	format (6(i3,i2))

c	**** output coordinates ****

	do i= 1,nnn

c	  **** get shortest dist of medium atom to high level hence dielectric constant ****
	  r2short= 1.d30
	  if (flags(i).eq.'M  ') then
	    do j= 1,nnn
	      if (flags(j).eq.'   ') then
		r2= (xo(i)-xo(j))**2 + (yo(i)-yo(j))**2 + (zo(i)-zo(j))**2 
		r2short= min (r2short,r2) 
	      end if
	    end do
	  end if
	  diel= 1.d0
	  rshort= sqrt (r2short)
	  if (r2short.ne.1.d30 .and. rshort.gt.5.) diel= 1.+(rshort-5.)**2

	  iuffat= 0
c	  if (hasfs4) then
c	    if (flags(i)(1:1).eq.'L') iuffat= 1
c	  else if (ifuff .and. flags(i)(1:1).eq.'M') then
c	    iuffat= 1
c	  end if

	  if (iuffat.eq.1) then
	    write (uffat,'(a2,1h-,a2,1h-,f6.3)')
     $		  atsym(nato(i)),atsym(nato(i)),qo(i)/diel
	    call delspace (uffat,nuffat)
	    write (2,961) uffat(1:nuffat),i,i,i,flags(i)
	  else
	    write (uffat,'(a2,1h-,a2,1h-,f6.3)')
     $		  atsym(nato(i)),attypeo(i),qo(i)/diel
	    call delspace (uffat,nuffat)
	    if (flags(i).ne.'M H') then
              if (iopttyp.eq.1) then
C                write (2,930) nato(i),xo(i),yo(i),zo(i),' '//flags(i)
                write (2,931) uffat(1:nuffat),xo(i),yo(i),zo(i),' '//flags(i)
	      else
	        write (2,961) uffat(1:nuffat),i,i,i,flags(i)
	      end if
            else
              iatom= nnntoatnumb(i)
              imol= amol(iatom)
              ires= ares(iatom)
      	      i0= molind(imol) - 1
              do j= 1,ncon(iatom)
                if (optfladj(icon(j,iatom)+i0).eq.'o' .or.
     $                   optfladj(icon(j,iatom)+i0).eq.'F') then
                  if (nat(icon(j,iatom)+i0).eq.7) then
                    htype= ' H'
                  else if (nat(icon(j,iatom)+i0).eq.6) then
                    htype= 'H1'
                  else if (nat(icon(j,iatom)+i0).eq.16) then                       
                    htype= 'HS'
                  else if (nat(icon(j,iatom)+i0).eq.8) then
                    htype= 'HO'
                  end if
                end if   
              end do   

c             **** Charge link atoms ****
              nlink= 0
              qlink= 0.d0

              do il= iares(ires,imol),jares(ires,imol)
                if (optfladj(il).eq.'L') nlink= nlink+1
                if (optfladj(il).eq.'F' .or. optfladj(il).eq.'o') then
                    qlink= qlink+q(il)                  
                end if
              end do
              qlink= -qlink/nlink/diel

              call delspace (htype,nl)
	      if (qlink.ge.0.0) then
	        write (cqlink,'(a,f5.3)') '-',qlink
	      else
	        write (cqlink,'(a,f5.3)') '--',-qlink
	      end if
              if (iopttyp.eq.1) then
C                write (2,930) nato(i),xo(i),yo(i),zo(i),' '//flags(i)
                write (2,931) uffat(1:nuffat),xo(i),yo(i),zo(i),flags(i),'-',htype(1:nl),cqlink
	      else
	        write (2,961) uffat(1:nuffat),i,i,i,flags(i),'-',htype(1:nl),cqlink
	      end if
              htype='  '
	    end if

	  end if
	end do

960	format (i3,' 0 x',i4.4,' y',i4.4,' z',i4.4,1x,a)
961	format (1x,a,' 0 x',i4.4,' y',i4.4,' z',i4.4,1x,4a)

        if (iopttyp.eq.1) then

c	  **** use modreduntant with berne opt ****
        
	  write (2,*)
	  do i= 1,nnn
	    if (freez(i)) write (2,950) i
	  end do

	else

c	  **** use z-matrix opt ****

	  write (2,*) 'Variables:'
	  do i= 1,nnn
	    if (.not.freez(i)) then
	      write (2,970) 'x',i,xo(i),'y',i,yo(i),'z',i,zo(i)
970	      format (a,i4.4,' = ',f10.4)
	    end if
	  end do

	  write (2,*) 'Constants:'
	  do i= 1,nnn
	    if (freez(i)) then
	      write (2,970) 'x',i,xo(i),'y',i,yo(i),'z',i,zo(i)
	    end if
	  end do

	end if

930	format (i2,3f10.4,1x,4a)
931	format (a,3f10.4,1x,4a)
950	format (i4,' F')

c	**** terminator for input coords ****
	write (2,*)
      
c       **** output of MM and basis set functions ****
	if (nnn.ne.nathigh) then
          write (2,988)
988	  format (
     $	  'HrmStr1   * * 331. 1.09' /
     $	  'HrmBnd1   * * * 20. 90.' /
     $	  'AmbTrs    * * * * 0 0 0 0 0.0 0.0 1.15 0.0 3.0' /
     $	  'VDW       MG 1.09 0.25' /
     $	  'VDW       NT 1.8240 0.17' / )
	end if

      if (isptchg) then

c	******** find point charges around outside **********

	qptch= 0.d0
	nptch= 0

	do kinc= 1,nres0
	  write (6,*) 'Charges around ',resname(ires0(kinc),imol0(kinc))
	  ist= iares(ires0(kinc),imol0(kinc))
	  ifn= jares(ires0(kinc),imol0(kinc))
	
	  do imol= 1,nmol
	   do ires= 1,nres(imol)
	    if (.not.incl(ires,imol)) then
	     i0= molind(imol) - 1

	     do i= iares(ires,imol),jares(ires,imol)
	      if (optfladj(i).ne.'F' .and. optfladj(i).ne.'L' .and. .not.repl(i)
     $		  .and. (nat(i).gt.1 .or. nat(icon(1,i)+i0).eq.6) ) then
		r2short= 1.d30
		do j= ist,ifn
		  r2short= min (r2short,hbl2(i,j))
		end do

		if (r2short.lt.rcutq2) then
	          diel= 1.d0
	          rshort= sqrt (r2short)
	          if (r2short.ne.1.d30 .and. rshort.gt.5.) diel= 1.+(rshort-5.)**2
		  diel= diel * 2.d0

		  qch= q(i)
		  if (nat(i).eq.6) then
c		    **** sum CH charges ****
		    do ll= 1,ncon(i)
		      l= icon(ll,i) + i0
		      if (optfladj(l).ne.'F' .and. .not.repl(l) .and. nat(l).eq.1) qch= qch + q(l)
		    end do
		  end if
		  qptch= qptch + qch
		  nptch= nptch + 1
		  write (2,958) x(i),y(i),z(i),qch/diel
C		  write (6,958) x(i),y(i),z(i),qch/diel,atname(i),resname(ires,imol),qptch
958		  format (3f8.3,f7.3,' 0.':1x,a,1x,a,f8.2)
	 	end if

	      end if
	     end do

	    end if
	   end do
	  end do

	end do

	maxptch= max (maxptch,nptch)
	write (6,*) 'total of point charge included=',nptch,qptch
	write (2,*)

      end if

c     **** basis sets ****

	jbasis= ibasis(1)
	i0= 1
	do i= 1,nnn
	  if (jbasis.ne.ibasis(i)) then
	    write (2,989) i0,i-1,basis(jbasis)
	    jbasis= ibasis(i)
	    i0= i
	  end if
	end do
	write (2,989) i0,nnn,basis(jbasis)
989     format (i4,'-',i4.4,' 0' / a / '****' )
	write (2,*)

c	**** output of charges  ****

	isaug=   naniono.gt.0 .or. nanions.gt.0

980	format (20i4)
981	format ('SP 1 1.00' / ' 0.0845 1.0 1.0' / '****' )
982	format ('SP 1 1.00' / ' 0.0405 1.0 1.0' / '****' )
983     format ('C 0' / 'C 3 0' / '  F and up ' / '3' /
     $          ' 2  0.08  -0.001286 ' / ' 2  0.12   0.003475 ' /
     $          ' 2  0.00725 -0.0000019 ' / '   S - F ' /
     $          '1' / ' 2  1.0  0.0 ' / '   P - F ' /
     $          '1' / ' 2  1.0  0.0 ' / '   D - F ' /
     $          '1' / ' 2  1.0  0.0 ' ) 


C        if (hasc) write (2,983) 
C	write (2,*)

	write (2,990)
990	format ('END')

	close (1)
	close (2)

c	********* save coords for viewing, only selected images ********

2000	continue

C      goto 2256
	nimage= nimage1
	natom= natom1
	nmol= nmol1
	do im= 2,nimage
	  do j= 1,natom1
	    hissel(j,im)= hissel(j+(im-1)*nimage1,1)
	  end do
	end do

	fname(nf:nf+3)= '.hin'
C	call writehin ('opt/'//fname,-1)

	nimage= 1
	natom= natom1*nimage1
	nmol= nmol1*nimage1
2256  continue

	fname(nf:nf+3)= '.com'

c	**** reset coords of changed atoms ****

	write (6,*) 'resetting the',ns,' dummy atoms'
	do i= 1,ns
	  j= ihvs(i)
	  x(j)= xs(i)
	  y(j)= ys(i)
	  z(j)= zs(i)
	  nat(j)= nats(i)
	  atname(j)= atns(i)
	  ats(j)= atsym(nats(i))
	end do

c	********* return status *********

        natsel= nnn
C	if (ifkill.eq.1) then
C	  natsel= 0
C	  nathigh= 0
C	end if

C        if (.not.hasfs4) then
C        if (hasfs4) then
C        if (hasfs4 .or. nh2onearby.ne.0) then
c	  **** set nber atoms=0 to skip this job ****
C          natsel= 0
C          nathigh= 0
C	else
C	  natsel= nnn
C	end if

	return

      end

c     ********************************************************

      subroutine getjname (jobnam,resname,imol)

      character*16 jobnam
      character*12 resname

      jobnam= ' '
      i= 0
      do j= 1,12
	if (resname(j:j).ne.' ') then
	  i= i + 1
	  jobnam(i:i)= resname(j:j)
	end if
      end do
C      write (jobnam(i+1:i+3),'(i3.3)') imol

      return
      end


c     ********************************************************

      subroutine delspace (s,n)

      character*(*) s

c     **** removes spaces from a string and returns its length ****

      n1= len (s)
      n= 0
      do j= 1,n1
	if (s(j:j).ne.' ') then
	  n= n + 1
	  s(n:n)= s(j:j)
	end if
      end do

      return
      end

c     ***************************************************************

      logical function shortdis (i,j,rcutou2)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      common /honly/ ifall


C      shortdis= .false.
C      return

      if (nat(i).gt.1 .and. nat(j).gt.1 .or.
     $	attype(i).eq.'HO' .or. attype(i).eq.'HS' .or.
     $  resname(ares(i),amol(i))(1:3).eq.'HOH' .or.
     $	attype(j).eq.'HO' .or. attype(j).eq.'HS' .or.
     $  resname(ares(j),amol(j))(1:3).eq.'HOH') then

c	**** heavy - heavy  or involves HOH or a rotatable H ****
     	shortdis= hbl2(i,j).lt.rcutou2

      else if (ifall.eq.-1) then
c	**** local H opt only **** 
     	shortdis= .FALSE.

      else
     	shortdis= hbl2(i,j).lt.7.84

      end if

      return
      end

c     ***************************************************************

      logical function hbond (i,j)

c     **** returnd true if hyd bond or ligand or pi-stacked interaction ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      real*8 mii(3),miveci(3,3),mivecj(3,3),cm(3)
      logical hbond1,typeok,ether,water,mg

      parameter (mpi=12)
      integer pi1(mpi),pi2(mpi)
      character*3 pires(mpi)
      data pires /'BC1','BP1','CL1','CL2','PQN','HIS','His','HiS','PHE','TYR','TYH','TRP'/
      data pi1   /    1,   1,   1,    1,    1,    6,    6,    6,    7,    6,    6,    6/
      data pi2   /   34,  34,  34,   34,   20,   10,   10,   10,   12,   11,   11,   14/

C      hbond= .false.
C      return

      ipi= 0
      natc= 0
      i0= molind(amol(i)) - 1
      j0= molind(amol(j)) - 1
      i1= iares(ares(i),amol(i)) - 1
      j1= iares(ares(j),amol(j)) - 1
      water= resname(ares(j),amol(j))(1:3).eq.'HOH' .or. resname(ares(i),amol(i))(1:3).eq.'HOH'
      mg= nat(i).eq.12 .or. nat(j).eq.12

C	write (6,*) j,j0
C	write (6,*) ncon(j)
C	write (6,*) icon(2,j)

      if (nat(i).eq.1) then
c	**** i= H, ic= N,O,S con to H, j= acceptor N,O,S ****
	ic= icon(1,i) + i0
	natc= nat(ic)
	natj= nat(j)
	ether= ncon(j).eq.2
	if (ether) ether= nat(j0+icon(2,j)) .ne. 1
	typeok= (natc.eq.7 .or. natc.eq.8 .or. natc.eq.16) .and.
     $          (natj.eq.7 .and. ncon(j).eq.2 .or. natj.eq.8 .and. .not.ether .or. natj.eq.16) 

      else if (nat(j).eq.1) then
c	**** j= H, jc= N,O,S con to H, i= acceptor N,O,S ****
	jc= icon(1,j) + j0
	natc= nat(jc)
	nati= nat(i)
	ether= ncon(i).eq.2
	if (ether) ether= nat(i0+icon(2,i)) .ne. 1
	typeok= (natc.eq.7 .or. natc.eq.8 .or. natc.eq.16) .and.
     $          (nati.eq.7 .and. ncon(i).eq.2 .or. nati.eq.8 .and. .not.ether .or. nati.eq.16) 

      else if (nat(i).eq.12) then
        typeok= nat(j).eq.7 .or. nat(j).eq.8 .or. nat(j).eq.16

      else if (nat(j).eq.12) then
        typeok= nat(i).eq.7 .or. nat(i).eq.8 .or. nat(j).eq.16

      else
c	**** two heavy atoms not Mg, check for pi stacking ****
C	ipi= 1
C	do while (ipi.le.mpi .and. pires(ipi).ne.resname(ares(i),amol(i))(1:3))
C	  ipi= ipi + 1
C	end do
C	if (ipi.le.mpi .and. i-i1.ge.pi1(ipi) .and. i-i1.le.pi2(ipi)) then
C	  jpi= 1
C	  do while (jpi.le.mpi .and. pires(jpi).ne.resname(ares(j),amol(j))(1:3))
C	    jpi= jpi + 1
C	  end do
C	  typeok= jpi.le.mpi .and. j-j1.ge.pi1(jpi) .and. j-j1.le.pi2(jpi)
C	else
C	  typeok= .false.
C	end if

      end if

       hbond1=.false.

      if (typeok) then
	r2= hbl2(i,j)

        if (ipi.gt.0) then
C          hbond1= r2 .lt. 16.0
C	  if (hbond1) then  
Cc	    **** check for alligned pi planes ****
C	    call mominert (i1+pi1(ipi),i1+pi2(ipi),mii,miveci,cm)
C	    call mominert (j1+pi1(jpi),j1+pi2(jpi),mii,mivecj,cm)
C	    cang= miveci(1,3)*mivecj(1,3) + miveci(2,3)*mivecj(2,3) + miveci(3,3)*mivecj(3,3) 
C	    hbond1= abs(cang) .gt. 0.71
C            if (hbond1) write (6,'(a,2i3,2f7.3)') ' PI STACK FOUND:',ipi,jpi,cang,sqrt(r2)
C	  end if
	else if (mg) then
          hbond1= r2 .lt. 6.0
C          hbond1= r2 .lt. 6.0
	else if (water) then
          hbond1= r2 .lt. 5.0
C          hbond1= r2 .lt. 16.0
	else
          hbond1= r2 .lt. 5.0
C          hbond1= r2 .lt. 8.0
	end if

      else
	hbond1= .false.
      end if

C      if (resname(ares(j),amol(j))(1:3).eq.'GLN')
C     $ write (6,888) i,j,nat(i),nat(j),natc,atname(i),atname(j),resname(ares(i),amol(i)),
C     $ resname(ares(j),amol(j)),sqrt(hbl2(i,j)),typeok,hbond1
C888   format (' hbond:',2i7,3i3,4(1x,a),f7.2,2l2)

      hbond= hbond1

      return
      end

c     ***************************************************************

      subroutine brbond (i,j,repl)

c     ***** breaks bond from i to j, replacing j by H ****
c     ***** also saves old coords and addw new to output list ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

      parameter (msave=100)
      character*5 atns
      logical repl(m)
      common /save/ xs(msave),ys(msave),zs(msave),
     $		nats(msave),ihvs(msave), ns, atns(msave)

	    repl(j)= .true.

c	    **** save orig specs ****
	    ns= ns + 1
	    if (ns.gt.msave) stop 'MSAVE'

	    ihvs(ns)= j
	    nats(ns)= nat(j)
	    xs(ns)= x(j)
	    ys(ns)= y(j)
	    zs(ns)= z(j)
	    atns(ns)= atname(j)

c	    **** replace i->j with i->H ****
	    rnow= sqrt (hbl2(i,j))
	    r= 1.1
	    nat(j)= 1
	    ats(j)= ' H'
	    dx= x(j) - x(i)
	    dy= y(j) - y(i)
	    dz= z(j) - z(i)
    	    x(j)= x(i) + dx * r/rnow
    	    y(j)= y(i) + dy * r/rnow
    	    z(j)= z(i) + dz * r/rnow
	    atname(j)= '  dum'

	return
	end

c     ***************************************************************

      subroutine tooutl (j,optfl)

c     ***** adds atom j to output list ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      character*2 attypeo
      character*1 optfl

      parameter (mo=3000, manion=50)
      common /outp/ qo(m),xo(mo),yo(mo),zo(mo),nato(mo),nnn,freez(mo),
     $  flags(m),hasfs4,anions(manion),aniono(manion),nanions,naniono,
     $  attypeo(mo),nnntoatnumb(m),ibasis(mo)
      character*3 flags,s
      logical freez,hasfs4
      integer anions,aniono
      

	    nnn= nnn + 1
      	    if (nnn.gt.mo) stop 'MO'
	    ibasis(nnn)= 1
            nnntoatnumb(nnn)= j
    	    nato(nnn)= nat(j)
	    qo(nnn)= q(j)
	    xo(nnn)= x(j)
	    yo(nnn)= y(j)
	    zo(nnn)= z(j)
	    attypeo(nnn)= attype(j)
	    freez(nnn)= optfl.ne.'o'
	    flags(nnn)= '   '
	    if (optfl.eq.'f') flags(nnn)= 'M  '
	    if (optfl.eq.'L') flags(nnn)= 'M H'
	    kres= ares(j)
	    kmol= amol(j)
	    k1= iares(kres,kmol) - 1
	    if (resname(kres,kmol)(1:3).eq.'FS4' .and. atname(j).ne.'  dum' .and. optfl.eq.'f') then
c	      **** flag hasfs4 indicates UFF layer present containing FS4 fragment ****
	      hasfs4= .true.
	      flags(nnn)= 'L  '
	    end if

	    if (flags(nnn)(1:1).eq.' ') then
	    s= resname(kres,kmol)(1:3)
	    if ( (s.eq.'CyS' .or. s.eq.'FS4') .and. nat(j).eq.16) then
	      nanions= nanions + 1
	      anions(nanions)= nnn
	      ibasis(nnn)= 2
	      write (6,*) 'anion S',nnn
	    end if
	    if (attype(j).eq.'O2' .or. s.eq.'LHG' .and.
     $		(atname(j).eq.'  OO3' .or. atname(j).eq.'  OO6')) then
	      naniono= naniono + 1
	      aniono(naniono)= nnn
	      ibasis(nnn)= 2
	      write (6,*) 'anion O',nnn
	    end if
	    end if

    	    write (1,980) optfl,nnn,j,atname(j),j-k1,
     $		      resname(kres,kmol),kmol,kres
	    call hselatom1 (j,1)

      return

980	    format (1x,a,i5,i6,1x,a,i4,1x,a,i5,i4)
      end

c     ****************************************************************

      subroutine mominert (ia1,ia2,mii,mivec,cm)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      real*8 mii(3),mivec(3,3),cm(3),wk(3)

c     **** get c of mass, all masses equal ****
      cm(1)= 0.d0
      cm(2)= 0.d0
      cm(3)= 0.d0
      do i= ia1,ia2
	cm(1)= cm(1) + x(i)
	cm(2)= cm(2) + y(i)
	cm(3)= cm(3) + z(i)
      end do
      tmass= ia2 - ia1 + 1
      cm(1)= cm(1) / tmass
      cm(2)= cm(2) / tmass
      cm(3)= cm(3) / tmass

c     **** moment of inertia tensor ****
      mivec(1,1)= 0.d0
      mivec(2,1)= 0.d0
      mivec(3,1)= 0.d0
      mivec(2,2)= 0.d0
      mivec(3,2)= 0.d0
      mivec(3,3)= 0.d0
      do i= ia1,ia2
	mivec(1,1)= mivec(1,1) + (y(i)-cm(2))**2 + (z(i)-cm(3))**2
	mivec(2,2)= mivec(2,2) + (x(i)-cm(1))**2 + (z(i)-cm(3))**2
	mivec(3,3)= mivec(3,3) + (x(i)-cm(1))**2 + (y(i)-cm(2))**2
	mivec(2,1)= mivec(2,1) - (x(i)-cm(1)) * (y(i)-cm(2))
	mivec(3,1)= mivec(3,1) - (x(i)-cm(1)) * (z(i)-cm(3))
	mivec(3,2)= mivec(3,2) - (y(i)-cm(2)) * (z(i)-cm(3))
      end do

c     **** diagonalize inertia matrix ****
      call tred2e (3,3,mivec,mii,wk,mivec)
      call tql2e  (3,3,      mii,wk,mivec,ier)

      return
      end

