
c     **** read templates ****

      if (windows) then
	write (6,*) 'WINDOWS op system'
        call system ('dir /o residues > residues.list')
        open (40,file='residues.list',status='old')
        do i= 1,7
	  read (40,*)
        end do
        resfname= 'residues\'
      else
	write (6,*) 'UNIX op system'
        call system ('ls -l residues > residues.list')
        open (40,file='residues.list',status='old')
        resfname= 'residues/'
      end if
      
C12345678901234567890123456789012345678901234567890
C22/06/2009  09:24 PM                14 aLA.LUU

      i= 1
29    continue
      if (windows) then
	read (40,'(39x,a)',end=33) resfname(10:16)
      else
	read (40,'(a)',end=33) resfname(10:16)
      end if
      if (resfname(10:10).eq.' ') goto 33
      if (resfname(13:13).eq.' ') goto 29
        open (2,file=resfname,status='old')
	if (i.gt.mt2) stop 'mt2 too small in template.cmn'
	name(i)= resfname(10:12)
	if (name(i)(3:3).eq.'.') name(i)= ' ' // name(i)(1:2)
	write (6,'(i4,a,a30,a,a)') i,' opening template: ',
     $		resfname, 'for ',name(i)
	read (2,*)
	j= 1
	qtot= 0.d0
	nhdef(i)= 0
30	read (2,32,end=31) kk,dname(j,i),ds(j,i),dtype(j,i),dflag(j,i),
     $	  dq(j,i),dx(j,i),dy(j,i),dz(j,i),
     $	  dncon(j,i),(dicon(k,j,i),dscon(k,j,i),k=1,dncon(j,i))
C	write (6,32) kk,dname(j,i),ds(j,i),dtype(j,i),dflag(j,i),
C     $	  dq(j,i),dx(j,i),dy(j,i),dz(j,i),
C     $	  dncon(j,i),(dicon(k,j,i),dscon(k,j,i),k=1,dncon(j,i))
32	  format (i4,1x,a5,1x,a2,1x,a2,1x,a6,f10.6,3f9.3,i4,20(i4,1x,a1))
	  qtot= qtot + dq(j,i)
	  natdef(j,i)= -1000
	  do k= -1,103
	    if (ds(j,i).eq.atsym(k)) natdef(j,i)= k
	  end do
	  if (natdef(j,i).eq.-1000) then
	    write (6,*) name(i),j,' ',ds(j,i)
	    stop 'unknown symbol in def'
	  end if
	  if (natdef(j,i).gt.1) nhdef(i)= nhdef(i) + 1
	  j= j + 1
	  if (j.eq.ma) stop 'too many atoms in residue'
	  goto 30
31	continue
	ndef(i)= j-1
	close (2)
C	if (abs(qtot).gt.1.d-5) write (6,'(1x,a,a,f12.6)')
	if (abs(qtot-nint(qtot)).gt.1.d-5) write (6,'(1x,a,a,f12.6)')
     $		'WARNING: TEMPLATE CHARGE= ',name(i),qtot
	i= i + 1
	goto 29

33      continue
	close (40)
	nt2= i - 1
	write (6,*) 'Nber template files read=',nt2

c     **** test of template start and end residues ****

      do it= mt1,nt2
	c= name(it)(1:1)
	if (c.ge.'A' .and. c.le.'Z') then
	  if (name(it)(2:2).ge.'a' .and. name(it)(3:3).ge.'A' .and. name(it)(3:3).le.'Z') then
	    uufile= name(it)(1:3)
	    call upcase (uufile)
	    jt= mt1
	    do while (jt.le.nt2 .and. name(jt)(1:3).ne.uufile)
	      jt= jt + 1
	    end do
	    if (jt.le.nt2) then
c	      **** match found, check atom number ****
	      if (uufile.ne.'LHG' .and. uufile.ne.'BCR' .and. uufile.ne.'PQN' .and.uufile.ne.'LMG' .and.
     $		 ndef(jt).ne.ndef(it)-2) then
		write (6,*) 'ERROR in TEMPLATE ',name(it),' ',name(jt)
		stop 'template read error in template.for'
	      end if
	    end if

	  else if (name(it)(3:3).ge.'a' .and. name(it)(2:2).ge.'A' .and. name(it)(2:2).le.'Z') then
	    uufile= name(it)(1:3)
	    call upcase (uufile)
	    jt= mt1
	    do while (jt.le.nt2 .and. name(jt)(1:3).ne.uufile)
	      jt= jt + 1
	    end do
	    if (jt.le.nt2) then
c	      **** match found, check atom number ****
	      if (uufile.ne.'LHG' .and. uufile.ne.'BCR' .and. uufile.ne.'PQN' .and.uufile.ne.'LMG' .and.
     $		 ndef(jt).ne.ndef(it)-1) then
		write (6,*) 'ERROR in TEMPLATE ',name(it),' ',name(jt)
		stop 'template read error in template.for'
	      end if
	    end if
	  end if
	end if
      end do

