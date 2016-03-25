c       **** reenters coords from .log file into .hin using the .sel file ****

c	**** USAGE:   reenter <hinname>.hin < reinsert-<OPT_TYPE>.dat

	implicit real*8 (a-h,o-z)
	include 'readhin.cmn'
	real*8 xn(m),yn(m),zn(m),xo(m),yo(m),zo(m)
	integer natn(m)
	logical endopt,notch3
	character*40 fname,hinname
	character*1 opt
	character*5 atnameopt

        parameter (mhist=100)
	integer nhe(0:mhist),nhh(0:mhist),nhx(0:mhist)
	data dele,delh,delx/0.2,0.02,0.005/

C	read (5,*) hinname
	call getarg (1,hinname)
	write (6,*) hinname
	nf= len_trim (hinname)
	hinname(nf+1:nf+4)= '.hin'
	nf= nf + 1

	write (6,*) 're-'//hinname(1:nf)//'out'
	open (4,file='re-'//hinname(1:nf)//'out',status='unknown')
	if (.not.readhin (hinname,0)) stop 'cant open hin file'

	open (10,file='re-'//hinname(1:nf)//'unstarted',status='unknown')
	open (7,file='re-'//hinname(1:nf)//'err',status='unknown')
	open (8,file='re-'//hinname(1:nf)//'qc' ,status='unknown')
	open (9,file='re-'//hinname(1:nf)//'plo',status='unknown')

	hchmax= 0.d0
	xchmax= 0.d0
	echmax= 1.d0
	nfail= 0
	nprem= 0
	nunstarted= 0
	nh= 0
	nx= 0
	ne= 0
	do i= 0,mhist
	  nhe(i)= 0
	  nhh(i)= 0
	  nhx(i)= 0
	end do

10	continue

	write (4,*)
C	write (6,*) 'Enter root name input .log and .sel files'
	read (5,'(a)',err=500,end=500) fname
	write (4,*) 'inserting from ',fname
	if (fname(1:1).eq.'-') goto 500
	nf= 1
	do while (fname(nf:nf).ne.' ')
	  nf= nf + 1
	end do

	fname(nf:nf+3)= '.log'
	if (.not. endopt (fname,natn,xn,yn,zn,xo,yo,zo,nn,m,echange,
     $	   nfail,nprem,nunstarted, echmax,hchmax,xchmax)	) goto 10

	irch= nint (-echange/dele)
	irch= min (mhist, max (0,irch) )
	nhe(irch)= nhe(irch) + 1
	ne= ne + 1
	if (echange.gt.0) stop 'echange +ve'

	fname(nf:nf+3)= '.sel'
C	open (1,file='opt/'//fname,status='unknown')
	open (1,file=fname,status='unknown')

c	*********

	do ii= 1,nn

	read (1,720) opt,iim,atnameopt,natopt
720	    format (1x,a1,5x,i6,1x,a,i4,1x,a,i5,i4)
C720	    format (1x,a1,i5,i6,1x,a,i4,1x,a,2i4)

c	**** watch out for an image being optimized, back transform ****

	im= (iim-1)/natom + 1
	i= iim - (im-1)*natom
C	write (6,*) atnameopt,opt,natopt,i,natn(ii),nat(i)
	if (im.ne.1) then
     	  if (opt.eq.'o') write (6,'(a,3i6,1x,2a)') 'opt of image ',iim,im,i,resname(ares(i),amol(i)),atname(i)
	  xx= xn(ii) - rotim(10,im)
	  yy= yn(ii) - rotim(11,im)
	  zz= zn(ii) - rotim(12,im)
          xn(ii)= rotim(1,im)*xx + rotim(4,im)*yy + rotim(7,im)*zz
          yn(ii)= rotim(2,im)*xx + rotim(5,im)*yy + rotim(8,im)*zz
          zn(ii)= rotim(3,im)*xx + rotim(6,im)*yy + rotim(9,im)*zz
	end if

	if (atnameopt.ne.'  dum' .and. atnameopt.ne.'dummy') then

     	  if (opt.eq.'o' .and. natn(ii).ne.nat(i)) then
	    write (6,*) opt,i,atnameopt,natopt
	    write (6,*) natn(ii)
	    write (6,*) nat(i),atname(i),' ',resname(ares(i),amol(i))
	    write (6,*) 'NAT CHANGED, exiting'
	    goto 490
C	    stop 'NAT changed'
	  end if

	  if (opt.eq.'o') then
	    rch= sqrt ( (x(i)-xn(ii))**2 + (y(i)-yn(ii))**2 +
     $			  (z(i)-zn(ii))**2 )
	    if (nat(i).eq.1) then
	        irch= nint (rch/delh)
		irch= min (mhist, max (0,irch) )
		nhh(irch)= nhh(irch) + 1
		nh= nh + 1
		write (4,'(a,f6.3,a,i5)') ' H change=',rch,atname(i),nh
		hchmax= max (hchmax,rch)
		if (rch.gt.0.5) then
		  if (notch3(i)) then
		    write (4,737) resname(ares(i),amol(i)),rch,atname(i)
		    write (6,737) resname(ares(i),amol(i)),rch,atname(i)
737		    format (1x,a,' LARGE POS change=',f6.3,' atom ',a)
		  end if
		end if
	    else
	        irch= nint (rch/delx)
		irch= min (mhist, max (0,irch) )
		nhx(irch)= nhx(irch) + 1
		nx= nx + 1
		write (4,'(a,f6.3,a,i5)') ' X change=',rch,atname(i),nx
		xchmax= max (xchmax,rch)
		if (rch.gt.0.2) then
		  write (4,737) resname(ares(i),amol(i)),rch,atname(i)
		  write (6,737) resname(ares(i),amol(i)),rch,atname(i)
		end if
	    end if

            x(i)= xn(ii)
            y(i)= yn(ii)
            z(i)= zn(ii)
	  end if
	end if

	end do

490	continue
	close (1)
	goto 10

500   continue

c     ***** write histograms ****

      write (6,*) 'histogram changes of H=',nh
      write (6,*) 'histogram changes of X=',nx
      write (6,*) 'histogram changes of E=',ne

      write (9,900)

      nhmax= 0
      nxmax= 0
      nemax= 0
      do i= 0,mhist
	nhmax= max (nhmax,nhh(i))
	nxmax= max (nxmax,nhx(i))
	nemax= max (nemax,nhe(i))
      end do

      ymax= nhmax/delh/nh
      write (9,910) 1,1,0*delh,mhist*delh,ymax,
     $		'H change / \gA'
      do i= 0,mhist
	yy= nhh(i)/delh/nh
	write (9,920) i*delh,yy
      end do

      if (nx.gt.0) then
      ymax= nxmax/delx/nx
      write (9,910) 1,2,0*delx,mhist*delx,ymax,
     $		'HVY change / \gA'
      do i= 0,mhist
	yy= nhx(i)/delx/nx
	write (9,920) i*delx,yy
      end do
      end if

      ymax= nemax/dele/ne
      write (9,910) 1,3,-mhist*dele,0.,ymax,
     $		'E change / kcal mol\u-1'
      do i= 0,mhist
	yy= nhe(i)/dele/ne
	write (9,920) -i*dele,yy
      end do

900   format ('tr' / 'gr 1 3 1' / 'ch 0.25' )
910   format ('mg',2i2 / 'xr ',2f8.2 / 'yr 0',f8.4 / 'xl ',a /
     $        'yl Probability' )
920   format (f9.3,f8.4)

      write (6,'(a,f6.2,a,2f7.3)') 'Max E change=',echmax,' H,X change=',hchmax,xchmax
      write (6,*) 'ERROR SUMMARY: number of failed/prem/unstarted end calcs=',nfail,nprem,nunstarted

c     ******** final coords *******

      call writehin (hinname,0)

      end

c     *************************************************************

      function endopt (fname,nat,x,y,z,xo,yo,zo,n,m,echange,nfail,nprem,nunstarted,
     $			echmax,hchmax,xchmax)

      implicit real*8 (a-h,o-z)
      real*8 xo(m),yo(m),zo(m),x(m),y(m),z(m)
      integer nat(m)
      logical scan1,scanok,twoc,mp2f,g98,endopt,eof,scanok1,badconv,g09
      character*3 ptgrp
      character*2 atsym(103)
      character*(*) fname
      character*80 line,command,ptline,scfline,mp2line,errline
      character*6 etype(-2:7)

      data etype / 
     $		'PW91' ,'BLYP',
     $		'SCF','MP2','B3LYP','SVWN','AM1' ,'CIS','PM3', 'ONIOM' /
      data au2ang,au2cm,amu2au /0.529177249D0,219477.D0,1822.845D0/
      data atsym	/
     +  'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',
     1	'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',
     2  'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
     3  'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',
     4  'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
     5  'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
     6  'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
     7  'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
     8  'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
     9  'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
     +  'Md', 'No', 'Lr' /

      open (15,file=fname,status='old',err=400)

      twoc= .true.
      mp2f= .true.
      ndata= 0
      emin= 1.d30
      nout= 0

c     **** check for g09 ****
      read (15,'(a)') line
      g09= inline(line,'G09') .gt. 0
      write (6,*) 'g09 flag=',g09

c     **** read job type ****

C1234567890
C #P MP2 opt freq cc-pvdz
      scanok= scan1 (15,1,3,' #P',line)
      write (4,*) line(1:75)
      command= line
      mp2= -1
      if (inline(line,'SCF').gt.0) mp2= 0
      if (inline(line,'RHF').gt.0) mp2= 0
      if (inline(line,'UHF').gt.0) mp2= 0
      if (inline(line,'MP2').gt.0) mp2= 1
      if (inline(line,'B3LYP').gt.0) mp2= 2
      if (inline(line,'SVWN').gt.0) mp2= 3
      if (inline(line,'AM1').gt.0) mp2= 4
      if (inline(line,'CIS').gt.0) mp2= 5
      if (inline(line,'PM3').gt.0) mp2= 6
      if (inline(line,'PW91').gt.0) mp2= -2
      if (inline(line,'BLYP').gt.0) mp2= -1
c	**** ONIOM must come last !! ****
      if (inline(line,'ONIOM').gt.0) mp2= 7
      if (mp2.eq.-1) stop 'cant determine energy type'
C      write (4,800) etype(mp2)
800   format (' energy type= ',a)

50    continue
      scanok= scan1 (15,2,18,'Center     Atomic',line)
C123456789012345678901234567890123456789012345678901234567890
C Center     Atomic              Coordinates (Angstroms)
      if (scanok .and. line(38:43).eq.'Forces') then
C	write (6,*) 'skipping forces'
        scanok= scan1 (15,2,18,'Center     Atomic',line)
      end if
C      if (.not.scanok) write (4,*) 'No more coords'
      if (.not.scanok) goto 200

C Full point group                 CS      NOp   2
C?      scanok= scan1 (15,2,11,'Full point',line)
C?      if (.not.scanok) goto 200
C?      ptline= line
C?      ptgrp= line(35:37)
C      write (4,*) 'point group= ',ptgrp

C      scanok= scan1 (15,2,18,'Center     Atomic',line)
C123456789012345678901234567890123456789012345678901234567890
C Center     Atomic                   Forces (Hartrees/Bohr)
C      if (.not.scanok .or. line(38:43).eq.'Forces') then
C	write (4,*) 'dont understand ... found ptgrp but no coords'
C	goto 200
C	stop 'coords'
C      end if

      g98= line(24:29).eq.'Atomic'
C      write (4,*) 'G98 flag= ',g98
      ndata= ndata + 1

      read (15,*)
      read (15,*)
      n= 0
100   continue
	if (g98 .or. g09) then
	  read (15,*,err=101) i,nati,ii,xx,yy,zz
	else
	  read (15,*,err=101) i,nati,xx,yy,zz
	end if
	if (i.gt.m) stop 'too many atoms'
	if (nati.eq.-1) goto 100
	if (nati.eq.0) goto 101
	n= n + 1
	nat(n)= nati
	xo(n)= xx
	yo(n)= yy
	zo(n)= zz
C	if (n.eq.1) write (4,'(i4,3f10.6)') n,xx,yy,zz
	goto 100
101   continue
C      write (4,*) 'nber atoms=',n
      n3= n*3

c     **** read scf energy ****

      if (mp2.le.3) then
c       **** PW91, BLYP, SCF, B3LYP, SVWN, and MP2 start ****
C1234567890123456789012345678901234567890123456789012345678901234567890
C SCF Done:  E(RB+HF-LYP) =  -209.190586604     A.U. after    1 cycles
C SCF Done:  E(RHF) =  -131.932796539     A.U. after    1 cycles
      scanok= scan1 (15,2,7,'SCF Do',line)
      if (.not.scanok) goto 200
      scfline= line
      ener= getval (line,1)

      else if (mp2.eq.4 .or. mp2.eq.6) then
c	**** AM1 or PM3 ****
C1234567890123456789012345678901234567890123456789012345678901234567890
C Energy=    0.039289784681 NIter=  13.
        scanok= scan1 (15,2,10,'Energy=  ',line)
        if (.not.scanok) goto 200
	read (line,'(8x,f18.12)') ener
        scfline= line

      else if (mp2.eq.7) then
c	**** ONIOM ****
C1234567890123456789012345678901234567890123456789012345678901234567890
C ONIOM: extrapolated energy =     -77.049133827538
        scanok= scan1 (15,2,10,'ONIOM: ex',line)
        if (.not.scanok) goto 200
	ener= getval (line,1)
        scfline= line

      else if (mp2.eq.5) then
c	**** CIS ****
C1234567890123456789012345678901234567890123456789012345678901234567890
C This state for optimization and/or second-order correction.
C Total Energy, E(Cis) =  -262.500471977    
        scanok= scan1 (15,2,10,'This stat',line)
        if (.not.scanok) goto 200
	read (15,'(a)') line
        scfline= line
        ener= getval (line,1)

      end if

      if (mp2.eq.1) then
c	**** read MP2 energy ****
C1234567890123456789012345678901234567890123456789012345678901234567890
C E2=       -0.4297119674D+00 EUMP2=       -0.13236250850661D+03
C E2 =    -0.1074775235D+01 EUMP2 =    -0.36084832207189D+03
        if (mp2f) scanok= scan1 (15,30,35,'EUMP2=',line)
        if (.not.mp2f .or. .not.scanok) then
	  mp2f= .false.
	  if (.not.scanok) rewind 15
          scanok= scan1 (15,28,34,'EUMP2 =',line)
	end if
        if (.not.scanok) goto 200
	mp2line= line
	ener= getval (line,2)
      end if

c     **** only gets here if all data read ****
      nout= nout + 1
      if (nout.eq.1) estart= ener
      if (ener.lt.emin) then
        emin= ener
C        rewind (16)
        do i= 1,n
	  x(i)= xo(i)
	  y(i)= yo(i)
	  z(i)= zo(i)
C	  if (nat(i).ne.0) write(16,'(i3,3f12.6)') nat(i),x(i),y(i),z(i)
        end do
C        write (16,*)
C        write (16,'(a)') command,ptline,scfline
C        if (mp2.eq.1) write (16,'(a)') mp2line
      end if

c     **** read next set of coords and energies ****
      goto 50

200   continue

c     **** gets here after failed read ****

      echange= emin - estart

      write (4,*) 'nber atoms=',n
      write (4,*) 'Nber of coord sets read=',nout
      write (4,'(a,f12.5)') ' Minimum energy=',emin
      write (4,'(a,f12.5)') ' Change  energy=',echange
      echange= echange * 627.51
      echmax= min (echmax,echange)
      if (echange.lt.-10.) then
	write (4,727) fname,echange
	write (6,727) fname,echange
727	format (1x,a,' LARGE E CHANGE=',f7.1,' kcal/mol')
      end if

c     **** check for sensible energy found ****

      endopt= .true.
      scanok= .false.
      eof= .false.
      if (nout.eq.0 .or. emin.gt.1.d29) then
	write (6,*) 'TOTAL FAIL for ',fname
	write (4,*) 'TOTAL FAIL for ',fname
	nfail= nfail + 1
	endopt= .false.
	goto 280
      end if

c     ****** scan for opt completed ****

      rewind 15
      scanok= scan1 (15,1,25,' Optimization completed.',line)
      eof= .not. scanok

280   continue

c     **** skip to end of file, print last 4 lines ****

      if (.not.eof) scanok1= scan1 (15,1,3,'!@#',line)
      backspace 15
      backspace 15
      backspace 15
      backspace 15
      backspace 15
      do i= 1,4
        read (15,'(a)') line
        write (4,'(a)') line
	if (.not.scanok .and. nout.lt.4) write (6,'(a)') line
	if (i.eq.1) badconv= line(1:38).eq.' Convergence failure -- run terminated'
	if (i.eq.1) errline= line
      end do

      close (15)

C      if (emin.gt.1.d3) then
C	write (6,*) 'ERROR: bad energy ',fname,nout
C	stop 'ERROR: bad energy'
      if (scanok) then
	write (4,*) 'opt completed OK ',fname
      else
	write (4,*) 'OPT COMPLETED FLAG NOT FOUND ',fname
	write (6,*) 'OPT COMPLETED FLAG NOT FOUND ',fname
	if (nout.lt.4) then
	  write (4,864) fname,nout,errline
	  write (6,864) fname,nout,errline
864	  format (' ERROR for ',a12,' OPT INCOMPL & CYCLES=',i2,1x,a40)
	  write (7,'(a)') fname
	  if (badconv) write (8,'(a)') fname
	  nprem= nprem + 1
	end if
      end if

      return

400   continue
      write (6,*) 'ERROR- file not found: ',fname
      write (4,*) 'ERROR- file not found: ',fname
      write (10,'(a)') fname
      endopt= .false.
      nunstarted= nunstarted + 1

      return

      end

      function scan1 (nu,i,j,str,line)

c     **** function to scan1 input till a pattern is matched ****

      character*80 line
      character*(*) str
      integer i,j,nu
      logical scan1

      scan1= .true.
100   read (nu,'(a80)',end=200) line
C	write (6,'(2i5,5a)') i,j,' "',line(i:j),'" "',str,'"'
	if (line(i:j).eq.str) return
	goto 100

200   scan1= .false.
      return

      end

c     ********************************************************************

      function getval (line,n)

c     **** reads the n-th value from a line in format ???? = <VALUE> ****

      real*8		getval,e
      integer		j,k,i,n,ll
      character*(*)	line

C      write (4,*) line
      ll= len (line)
      j= 1
      i= 1
      do while (j.le.ll .and. (line(j:j).ne.'=' .or. i.lt.n) )
	if (line(j:j).eq.'=') i= i + 1
	j= j + 1
      end do
      if (j.ge.ll) stop 'GETVAL: value not found'

      k= j + 1
      do while (line(k:k).eq.' ')
	k= k + 1
      end do

      do while (line(k:k).ne.' ')
	k= k + 1
      end do

C      write (4,*) j,k,'"',line(j+1:k-1),'"'
      read (line(j+1:k-1),*,err=200) e
      getval= e
      return

200   continue
      getval= 2.d29
      return

      end

c     ********************************************************************

      function inline (line,s)

c     *** returns the start of substring s in line, s in upper case ****

      character*80 line
      character*(*) s
      character*1 c
      logical found

      n= len (s)
      i= 1
      found= .false.
      do while (.not.found .and. i.le.80-n+1)
	found= .true.
	do j= 1,n
	  c= line(i+j-1:i+j-1)
	  if (c.ge.'a' .and. c.le.'z') c= char ( ichar(c) + ichar('A')
     $		- ichar('a') )
	  found= found .and. s(j:j).eq.c
	end do
	i= i + 1
      end do

      inline= 0
      if (found) inline= i
      return
      end

c     **************************************************************

      logical function notch3 (i)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

c     **** returns true if H atom i is not a CH3 hydrogen ****

      ires= ares(i)
      imol= amol(i)
      i0= molind(imol) - 1
      j= i0 + icon(1,i)
      notch3= .true.
      if (nat(j).ne.6) return
      if (ncon(j).ne.4) return

      nh= 0
      do kk= 1,ncon(j)
	k= icon(kk,j) + i0
	if (nat(k).eq.1) nh= nh + 1
      end do
      if (nh.ne.3) return

      notch3= .false.
      return

      end

