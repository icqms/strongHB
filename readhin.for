C     **** subroutines to read/write/modify in a .hin data file ****

c  **** 
c  **** atomic parameters are:
c  ****    atnumb atname ats attype flag q x y z ncon (icon,scon,i=1,ncon)
c  ****			also nat
c  **** 
c  **** indices from atom to residue and molecule:
c  ****    ares(n), amol(n)
c  **** 
c  **** indices from mol, res to atom:
c  **** 
c  ****    nres(imol)		... nber residues in this molecule
c  ****    molind(imol)		... first atom in molecule
c  ****    molend(imol)		... last atom in molecule
c  ****    iares(ires,imol)	... first atom of residue
c  ****    jares(ires,imol)	... last atom of residue
c  ****    nares(ires,imol)	... nber atom in residue
c  ****    q(ires,imol)	... residue charge
c  ****    issel(iatom,image)	... if this atom is selected in this image
c  **** 
c  **** entry points are:
c  ****   logical function readhin (filenm,if_fix_charge)
c  ****   subroutine writehin (filenm,if_images)
c  ****   subroutine hswap (imol,ind,nind)
c  ****   subroutine hcopy (i1,i2,n)
c  ****   subroutine hdel (ist,n)
c  ****   subroutine hrotate (conto,vecto,ang)
c  ****   subroutine hselres (ires,imol,n,image)   or hselres1 (ires,imol,n)
c  ****   subroutine hselatom (i1,n,image)	   or hselatom1 (i1,n)
c  ****   subroutine hdeselect (i1,n,image)	   or hdeselect (i1,n)
c  ****   subroutine hcopyim (i1,n)   ... creates image coords 
c  ****   subroutine him2real (nmol1,natom1,nimage1)   ... converts images into real atoms
c  ****   function hbl2i (i,j,image)		   or hbl2 (i,j)
c  ****   logical function hisbond (i,j)
c  **** 
c  **** internal routines:
c  ****   subroutine extract (line,param,m,n,je)
c  ****   subroutine wtrim (s)
c  ****   subroutine wpack (s)
c  **** 

      function readhin (filenm,iffixch)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      logical ff
      character*(*) filenm
      character*(mc) line
      character*6    fixch(0:1)

      data fixch / 'kept  ','fixed'/
      data natom/0/, nmol/0/
      data (rotim(i,1),i=1,12) / 1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0,
     $	0.d0,0.d0,1.d0, 0.d0,0.d0,0.d0 /

      write (4,'(2a)') '  READING: ',filenm
      write (6,'(2a)') '  READING: ',filenm
      open (90,file=filenm,status='old',err=400)

      hissel= .false.
      nsel= 0
      qtot= 0.D0
      ncres= 0
      ntres= 0
      ifhead= 1
      nhead= 0
      nimage= 0
      n= natom
      rinscr2= 0.d0
      alen= 0.d0
      pi= acos (0.d0) * 2.d0
      pre= 180.d0 / pi
      rt32= sqrt (3.d0) / 2.d0

      do i= nmol+1,mmol
	molind(i)= 0
	molend(i)= 0
      end do       

100   continue
	read (90,'(a)',end=200) line
C	write (6,'(a)') line

	if (line(1:4).eq.'mol ') then
	  call extract (line,param,2,nextract,kk)
C	  write (6,*) param(2)
	  read (param(2),'(i15)') imol1
	  nmol= nmol + 1
	  imol= nmol
	  if (imol.gt.mmol) stop 'max nber mols'
	  molind(imol)= n+1
	  nres(imol)= 0
	  ifhead= 0

	else if (line(1:7).eq.'endmol ') then
	  if (ifresn.eq.0) jares(ires,imol)= n
	  molend(imol)= n
C	  write (6,'(1x,a,i4,a,i4)') 'mol ',imol,' res=',nres(imol)
	  k= molend(imol) - molind(imol) + 1
	  do i= molind(imol),molend(imol)
	    do j= 1,ncon(i)
	      if (icon(j,i).le.0 .or. icon(j,i).gt.k) then
		write (6,*) imol,i,j,ncon(i),icon(j,i),k
		stop 'bond range in writehin'
	      end if
	    end do
	  end do

	else if (line(1:4).eq.'res ') then
	  call extract (line,param,2,nextract,kk)
	  nres(imol)= nres(imol) + 1
	  ires= nres(imol)
	  if (ires.gt.mres) stop 'mres too small'
	  read (param(2),'(i15)') ires1
	  resname(ires,imol)= line(kk+2:kk+13)
C	  write(6,*) 'resname= "',resname(ires,imol),'" ',kk,line(kk:kk)
	  qres(ires,imol)= 0.D0
	  nares(ires,imol)= 0
	  iares(ires,imol)= n + 1
	  ntres= ntres + 1
	  ifresn= 1

	else if (line(1:7).eq.'endres ') then
	  jares(ires,imol)= n

	else if (line(1:5).eq.'atom ') then

	  if (nres(imol).eq.0) then
	    write (6,*) 'WARNING: atom before residue number',imol
	    nres(imol)= nres(imol) + 1
	    ires= nres(imol)
	    if (ires.gt.mres) stop 'mres too small'
	    resname(ires,imol)= 'UNK 1 - -'
C	    write(6,*) 'resname= "',resname(ires,imol),'" ',kk,line(kk:kk)
	    qres(ires,imol)= 0.D0
	    nares(ires,imol)= 0
	    iares(ires,imol)= n + 1
	    ntres= ntres + 1
	    ifresn= 0
C	    stop 'missing residue number'
	  end if
	  n= n + 1
	  if (n.gt.m) stop 'nber atoms'
	  nares(ires,imol)= nares(ires,imol) + 1
	  ares(n)= ires
	  amol(n)= imol

	  call extract (line,param,mext,nextract,kk)
C	  write (6,*) 'nextract=',nextract
C	  do i= 1,nextract
C	    write (6,'(i3,1x,1h",A,1H")') I,PARAM(I)
C	  end do

	  ats(n)= param(4)(14:15)
	  if (ats(n)(1:1).eq.' ') ats(n)= param(4)(15:15)//' '
	  k= -1
	  do while (k.le.103 .and. ats(n).ne.atsym(k))
	    k= k + 1
	  end do
	  if (k.gt.103) then
	    write (6,*) n,k,' "',ats(n),'" "',param(4),'"'
	    stop 'unknown atom'
	  end if
	  nat(n)= k

	  read (param(2),'(i15)') atnumb(n)
	  if (n.eq.molind(imol)) then
c	    **** watch out for molecule split into separate ones ****
	    iatsft= 1 - atnumb(n)
	    if (iatsft.ne.0) write (6,*) imol,' atom nber shift by',iatsft
	  end if
	  read (param(3),'(10x,a5)') atname(n)
	  read (param(5),'(13x,a2)') attype(n)
	  read (param(6),'(9x,a6)') flag(n)
	  read (param(7),'(f15.0)') q(n)
	  read (param(8),'(f15.0)') x(n)
	  read (param(9),'(f15.0)') y(n)
	  read (param(10),'(f15.0)') z(n)
	  read (param(11),'(i15)') ncon(n)
	  if (ncon(n).gt.100000) then
	    do k= nextract,12,-1
	      param(k)= param(k-1)
	    end do
	    nextract= nextract + 1
	    param(12)(10:10)= ' '
	    ncon(n)= ncon(n) / 100000
	  end if
	  if (ncon(n).gt.mcon) then
	    write (6,*) imol,ires,atnumb(n),atname(n),ncon(n),nextract
	    write (6,*) resname(ires,imol)
	    write (6,*) line
	    stop 'too many conections for mcon'
	  end if
	  if (11+2*ncon(n) .gt. nextract) then
	    write (6,*) imol,ires,atnumb(n),atname(n),ncon(n),nextract
	    write (6,*) resname(ires,imol)
	    write (6,*) line
	    stop 'too few data entries on line'
	  end if
	  do i= 1,ncon(n)
	    read (param(10+2*i),'(i15)') icon(i,n)
	    icon(i,n)= icon(i,n) + iatsft
	    read (param(11+2*i),'(14x,a1)') scon(i,n)
	  end do
          if (11+2*ncon(n)+1 .le. nextract) then
c           **** read B value from hin file ****
            read (param(11+2*ncon(n)+1),'(f15.0)') bfact(n)
          else
            bfact(n)= 0.e0
          end if
          if (11+2*ncon(n)+7 .le. nextract) then
c           **** read anisotropies ****
            do ii= 2,7
              read (param(11+2*ncon(n)+ii),'(i15)') aniso(ii-1,n)
            end do
          else
            do ii= 1,6
              aniso(ii,n)= 0
            end do
          end if

	  j= 6
	  ff= .true.
	  do while (j.gt.0 .and. ff)
	    ff= flag(i)(j:j).ne.'s'
	    if (ff) j= j - 1
	  end do
	  hissel(n,1)= j.ne.0
	  if (hissel(n,1)) flag(i)(1:j)= ' '//flag(i)(1:j-1)

	  qtot= qtot + q(n)
	  qres(ires,imol)= qres(ires,imol) + q(n)
C	  write (6,'(a,1x,4i4,4f8.3)') line(1:30),n,imol,ires,k,
C     $		  q(n),x(n),y(n),z(n)

	else if (line(1:8).eq.';samemol') then
	  iatsft= n - molind(imol) + 1
	  write (6,*) ';samemol found ... changing atom numbers by ',iatsft

	else if (ifhead.eq.1) then
c	  if (nhead.ne.0 .or. line(1:1).ne.';') then
	    nhead= nhead + 1
	    if (nhead.gt.mh) stop 'header too big'
	    header(nhead)= line
	    if (line(1:5).eq.'; box') then
	      read (line(7:mc),*) alen,blen,clen,alpha,beta,gamma
	    end if
	    if (line(1:8).eq.'; rinscr') then
	      read (line(10:mc),*) rinscr
	      rinscr2= rinscr**2
	    end if
	    if (line(1:7).eq.'; image') then
	      nimage= nimage + 1
	      if (nimage.gt.mimage) stop 'MIMAGE too small'
	      read (line(9:mc),*) (rotim(j,nimage),j=1,12)
	      do j= 1,9
		if (abs(abs(rotim(j,nimage))-rt32).lt.1.d-3)
     $		  rotim(j,image)= sign (rt32,rotim(j,image))
	      end do
	    end if
c	  end if

	else
	  write (6,*) 'last mol,res=',nmol,nres(nmol)
	  write (6,*) line
	  stop 'unknown line after header'

	end if

	goto 100

200   continue

      close (90)

      natom= n
      readhin= .true.
      iqtot= 0

c     **** check residue charge ****

      do imol= 1,nmol
	do ires= 1,nres(imol)
	  iqres= nint (qres(ires,imol))
	  iqtot= iqtot + iqres
	  if (iqres.ne.0) ncres= ncres + 1

	  if (abs(qres(ires,imol)-iqres).gt.1.D-10) then

c	    **** fix or keep ****
	    if (iffixch.eq.1) stop 'CHARGE FIX NOT CODED'
	    ii= iffixch
	    if (ii.eq.1) then
C	      do j= 1,4
C		if (mol(j).eq.imol .and. nabcl.eq.nares(ires,imol)) ii=0
C		if (resmol(j).eq.imol .and. reslig(j).eq.ires .and.
C     $		  nabcl.eq.nares(1,mol(j)) ) ii= 0
C	      end do
	    end if

     	    write (4,650) imol,ires,resname(ires,imol),qres(ires,imol),
     $		fixch(ii)
c     	    write (6,650) imol,ires,resname(ires,imol),qres(ires,imol),
c     $		fixch(ii)
650	    format(' mol',i4,' nres',i4,1x,a,' charge',f10.6,1x,a)

	    if (ii.eq.0) then
c	      **** ignore ****

	    else if (abs(qres(ires,imol)+0.001D0).lt.0.001D0 .and.
     $		    resname(ires,imol)(1:3).eq.'MET') then
	      i= iares(ires,imol) + 6
	      q(i)= q(i) - qres(ires,imol)
	      qres(ires,imol)= 0.D0

	    else
	      delq= (qres(ires,imol) - iqres) / nares(ires,imol)
	      qmax= 0.D0
	      imax= iares(ires,imol)
	      iq= 0
	      do i= iares(ires,imol),jares(ires,imol)
	        q(i)= q(i) - delq
		iq= iq + nint(q(i)*1.D6)
		if (abs(q(i)).gt.qmax) then
		  qmax= abs(q(i))
		  imax= i
		end if
	      end do
	      q(imax)= q(imax) - iq/1.D6
	      qres(ires,imol)= qres(ires,imol) - delq
	    end if

	  end if

	end do
      end do

      if (alen.eq.0.d0) then
	alen= 1.d3
	blen= 1.d3
	clen= 1.d3
	alpha= 90.d0
	beta= 90.d0
	gamma= 90.d0
      end if

      ca= cos (alpha/pre)
      sa= sin (alpha/pre)
      cb= cos (beta/pre)
      sb= sin (beta/pre)
      cg= cos (gamma/pre)
      sg= sin (gamma/pre)
      boxvecs= 0.d0
      boxvecs(1)= alen
      boxvecs(4)= blen * cg
      boxvecs(5)= blen * sg
      boxvecs(7)= clen * cb
      boxvecs(8)= clen * (ca-cb*cg)/sg
      boxvecs(9)= clen * sqrt (1.d0-ca**2-cb**2-cg**2+2.d0*ca*cb*cg)
     $				/ sg
      boxvecsinv= 0.d0
      boxvecsinv(1)= 1.d0 / boxvecs(1)
      boxvecsinv(5)= 1.d0 / boxvecs(5)
      boxvecsinv(9)= 1.d0 / boxvecs(9)
      boxvecsinv(4)= - boxvecs(4) / boxvecs(1) / boxvecs(5)
      boxvecsinv(8)= - boxvecs(8) / boxvecs(9) / boxvecs(5)
      boxvecsinv(7)= (boxvecs(4)*boxvecs(8)-boxvecs(7)*boxvecs(5))
     $		    	/ boxvecs(1) / boxvecs(5) / boxvecs(9)

      if (rinscr2.eq.0.d0) rinscr2= 1.e6
      if (nimage.eq.0) nimage= 1
      call hcopyim (1,natom)

      write (6,880) n,nmol,ntres,ncres,qtot,iqtot,nimage,rinscr
      write (4,880) n,nmol,ntres,ncres,qtot,iqtot,nimage,rinscr
880   format (' number of atoms=',i7 / ' number of mols= ',i5 /
     $    ' number of res=  ',i5 /
     $    ' number of charged res=',i4,
     $    ' tot charge=',f11.6,' without roundoff=',i4 /
     $    ' nber images=',i5,' inscribed radius=',f8.3)

      return

c     ******************************************

400   continue
      readhin= .false.
      return

      end

c     **********************************************************************

      subroutine extract (line,param,m,n,je)

      character*(*) line
      character*15 param(m)

C      write (6,*) line

      je= 0
      do i= 1,m
	js= je + 1
	do while (line(js:js).eq.' ')
	  js= js + 1
	  if (js.gt.len(line)) return
	end do
	je= js
	do while (je.le.len(line) .and. line(je:je).ne.' ')
	  je= je + 1
	end do
	je= je - 1
	param(i)= ' '
	if (je-js.gt.14) stop 'param too long'
	param(i)(15-je+js:15)= line(js:je)
	n= i
C	write (6,*) i,js,je,' "',param(i),'"'
      end do

      return
      end

c     **********************************************************************

      subroutine wtrim (s)

c     **** trims blancks and writes line ****

      character*(*) s

      j= len(s)
      do while (j.gt.1 .and. s(j:j).eq.' ')
	j= j - 1
      end do

      write (90,'(a)') s(1:j)

      return
      end

c     ********************************************************************

      subroutine wpack (s)

      character*256 sout
      character*(*) s

      j= 1
      sout(1:1)= s(1:1)
      do i= 2,len(s)
	if (s(i-1:i).ne.'  ' ) then
	  j= j + 1
	  sout(j:j)= s(i:i)
	end if
      end do

      if (j.gt.1 .and. sout(j:j).eq.' ') j= j - 1
      write (90,'(a)') sout(1:j)

      return
      end

c     **********************************************************************

      subroutine writehin (filenm,if_images)

c     **** if_images= 0 (main coords only) or 1 (all images as well) ****
c     **** or -1 for main coords + selected atoms in images ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      integer nselinmol(mmol,mimage),iselinmol(m,mimage),ninc(mcon)
      character*6 flag1
      character*(*) filenm
      character*(256) line
      character*80 f800

      write (4,'(2a)') '  WRITING: ',filenm
      write (4,*) 'natoms,nimages,nmol=',natom,nimage,nmol
      write (6,'(2a)') '  WRITING: ',filenm
      write (6,*) 'natoms,nimages,nmol=',natom,nimage,nmol
      open (90,file=filenm,status='unknown')

c     **** if if_images < 0 then find indices of nber atoms sel in each mol and image ****
 
      if (if_images.lt.0) then
      nselinmol= 0
      iselinmol= 0
      do im= 1,nimage
	do imol= 1,nmol
	  do i= molind(imol),molend(imol)
	    if (hissel(i,im)) then
	      nselinmol(imol,im)= nselinmol(imol,im) + 1
	      iselinmol(i,im)= nselinmol(imol,im)
C	      write (6,'(a,5i6)') 'sel',im,imol,i,
C     $		nselinmol(imol,im),iselinmol(i,im)
	    end if
	  end do
	end do
      end do
      end if

      do i= 1,nhead
	if (header(i)(1:5).eq.'; box') write (header(i),'(a,6f9.3)')
     $		'; box',alen,blen,clen,alpha,beta,gamma
	call wtrim (header(i))
      end do

      imolout= 0

      do im= 1,max(1,abs(if_images)*nimage)

      do imol= 1,nmol

	i= iares(1,imol)

	if (nres(imol).gt.0) then
	  imolout= imolout + 1
	  write (line,810) imolout
	  call wpack (line)
	end if

	if (if_images.lt.0 .and. nselinmol(imol,im).eq.0) goto 334

	do ires= 1,nres(imol)
	  write (line,811) ires,resname(ires,imol)
	  call wpack (line)

	  do i= iares(ires,imol),jares(ires,imol)

C	    if (if_images.lt.0 .and. .not.hissel(i,im)) goto 300

	    iq= nint(abs(q(i))*1.D6)

	    if (flag(i).eq.'      ') flag(i)= '     -'
	    flag1= flag(i)

	    if (hissel(i,im)) then
	      if (flag1(6:6).eq.'-') then
		flag1= '     s'
	      else
		j= 6
		do while (j.gt.0 .and. flag1(j:j).ne.' ' .and.
     $		  flag1(j:j).ne.'s')
		  j= j - 1
		end do
		if (j.eq.0) then
		  write (6,*) i,resname(ires,imol),atname(i),' ',flag1
		  stop 'flag overflow'
		end if
		flag1(j:j)= 's'
	      end if

	    else
	      do j= 6,1,-1
		if (flag1(j:j).eq.'s') then
		  do k= j,2,-1
		    flag1(k:k)= flag1(k-1:k-1)
		  end do
		  flag1(1:1)= ' '
		end if
	      end do
	      if (flag1(6:6).eq.' ') flag1(6:6)= '-'
	    end if
C	    write (6,*) i,' flag= "',flag(i),'" "',flag1,'"'

	    if (mod(iq,1000).eq.0) then
C	      write (f800,8000)
	      write (f800,800)
	    else
C	      write (f800,8010)
	      write (f800,801)
	    end if
	    atnumb(i)= i - molind(imol) + 1

	    if (if_images.eq.0) then
	      write (line,f800) atnumb(i),atname(i),ats(i),attype(i),
     $		flag1,q(i),x(i),y(i),z(i),
     $		ncon(i),(icon(j,i),scon(j,i),j=1,ncon(i))
	    else
C	    else if (if_images.gt.0) then
	      write (line,f800) atnumb(i),atname(i),ats(i),attype(i),
     $		flag1,q(i),xi(i,im),yi(i,im),zi(i,im),
     $		ncon(i),(icon(j,i),scon(j,i),j=1,ncon(i))
C	    else
Cc	      **** delete links to unselected atoms ****
C	      ncon1= 0
C	      do j= 1,ncon(i)
C		if (iselinmol(icon(j,i)+amol(i)-1,im).ne.0) then
C		  ncon1= ncon1 + 1
C		  ninc(ncon1)= j
C		end if
C	      end do
C	      write (line,f800) iselinmol(i,im),atname(i),
C     $		ats(i),attype(i),
C     $		flag1,q(i),xi(i,im),yi(i,im),zi(i,im),
C     $		ncon1,(iselinmol(icon(ninc(j),i),im),
C     $		scon(ninc(j),i),j=1,ncon1)
	    end if

	    ii= len_trim(line)
	    write (line(ii+1:ii+7),'(f7.2)') bfact(i)
	    write (line(ii+8:ii+50),'(6i7)') (aniso(j,i),j=1,6)
	    call wpack (line)
300	    continue
	  end do

	  write (line,812) ires
	  call wpack (line)
	end do

334     continue

	if (nres(imol).gt.0) then
	  write (line,813) imolout
	  call wpack (line)
	end if

      end do
      end do

8000  format ('(4hatom,i7,2i5,4(1x,a),f8.3,3f11.4,i3,12(i7,1x,a))')
8010  format ('(4hatom,i7,2i5,4(1x,a),f11.6,3f11.4,i3,12(i7,1x,a))')
800   format ('(4hatom,i7,4(1x,a),f8.3,3f11.4,i3,12(i7,1x,a))')
801   format ('(4hatom,i7,4(1x,a),f11.6,3f11.4,i3,12(i7,1x,a))')

810   format ('mol ',i4)
811   format ('res ',i4,1x,a)
812   format ('endres ',i4)
813   format ('endmol ',i4)

      close (90)

      return

      end

c     *********************************************************************

      function hbl2 (i,j)

c     **** bond length squared between two atoms in main image ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

      hbl2= (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2

      return
      end

c     *********************************************************************

      function hbl2i (i,j,image)

c     **** bond length squared between two atoms, j in image ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

      hbl2= (x(i)-xi(j,image))**2 + (y(i)-yi(j,image))**2
     $    + (z(i)-zi(j,image))**2
      hbl2i= hbl2

      return
      end

c     *********************************************************************

      function hisbond (i,j)

c     **** is there a bond between atoms i and j ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

      imol= amol(i)
      jmol= amol(j)
      if (imol.ne.jmol) then
	hisbond= .false.
	return
      else
	i1= molind(imol)-1
	do k= 1,ncon(i)
	  if (icon(k,i).eq.j-i1) then
	    hisbond= .true.
	    return
	  end if
	end do
	hisbond= .false.
	return
      end if

      end

c     *********************************************************************

      subroutine hselres1 (ires,imol,n)
      call hselres (ires,imol,n,1)
      end

      subroutine hselres (ires,imol,n,image)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

c     **** select n residues starting at ires in imol in image ****

      if (nres(imol).lt.ires+n-1) stop 'HSELRES: residue too big'
      do i= ires,ires+n-1
	write (4,800) resname(i,imol),i,imol,qres(i,imol),image
	write (6,800) resname(i,imol),i,imol,qres(i,imol),image
800    format (' selecting res ',a,2i4,' q=',f6.3,' image=',i3)
        call hselatom (iares(i,imol),nares(i,imol),image)
      end do

      return
      end

c     *********************************************************************

      subroutine hselatom1 (i1,n)
      call hselatom (i1,n,1)
      end

      subroutine hselatom (i1,n,image)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

c     **** select n atoms starting at i1 in image ****

C      write (6,*) 'selecting atom',i1,n,image

      do i= i1,i1+n-1
	if (.not.hissel(i,image)) then
	  nsel= nsel + 1
	  hissel(i,image)= .true.
	end if
      end do

      return
      end

c     *********************************************************************

      subroutine hdeselect1 (i1,n)
      call hdeselect (i1,1,1)
      end

      subroutine hdeselect (i1,n,image)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

c     **** deselect atoms n atoms starting at i1 in image ****

      do i= i1,i1+n-1
	if (hissel(i,image)) then
	  nsel= nsel - 1
	  hissel(i,image)= .false.
	end if
      end do

C      nsel= 0
      return
      end

c     *********************************************************************

      subroutine hswap (imol,ind,nind)

c     **** swaps order of entries in .hin file, must be all same residue ****
c     **** ind(1,?) is from list,  ind(2,?) is to list ****
c     **** ? is number from START of MOLECULE ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

      integer ind(2,nind)

      j0= molind(imol) - 1

C      write (6,*) 'enter hswap:',nind

C	  do k= 1,nind
C	    write (6,'(5i5)') k,ind(1,k),ind(2,k),ares(j0+ind(1,k)),
C     $		ares(j0+ind(2,k))
C	  end do

c     **** direct move, change atom nber ****

      do i= 1,nind
	if (ares(j0+ind(1,i)) .ne. ares(j0+ind(2,i))) then
	  write (6,*) 'hswap error: residue different'
	  do k= 1,nind
	    write (6,'(5i5)') k,ind(1,k),ind(2,k),ares(j0+ind(1,k)),
     $		ares(j0+ind(2,k))
	  end do
	  stop 'residue diff'
	end if
	call hcopy (j0+ind(1,i),natom+i,1)
      end do

      do i= 1,nind
	call hcopy (natom+i,j0+ind(2,i),1)
	atnumb(j0+ind(2,i))= ind(2,i)
      end do

c     **** rearrange connection list ****

      do j= molind(imol),molend(imol)
	do k= 1,ncon(j)
	  i= 1
	  do while (i.le.nind .and. icon(k,j).ne.ind(1,i))
	    i= i + 1
	  end do
	  if (i.le.nind) icon(k,j)= ind(2,i)
	end do
      end do

      return

      end

c     *********************************************************************

      subroutine hcopy (i1,i2,n)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

c     **** subroutine to directly copy n atoms from i1 to i2 ****

c     **** detect direction of copy, run to avoid overwriting ****
      if (i1.lt.i2) then
	ist= i1+n-1
	ifn= i1
	idel= -1
      else
	ist= i1
	ifn= i1+n-1
	idel= 1
      end if

C      write (6,'(a,6i7)') ' hcopy:',i1,i2,n,ist,ifn,idel
      if (i2+n-1 .gt. m) stop 'hcopy: insufficient space'

      do i= ist,ifn,idel
	j= i2 - i1 + i
	ares(j)=   ares(i)
	amol(j)=   amol(i)
	atnumb(j)= atnumb(i)
C	if (i.eq.ist) write (6,*) 'first =',amol(i),ares(i),atnumb(i)
	atname(j)= atname(i)
	ats(j)=    ats(i)
	nat(j)=    nat(i)
	attype(j)= attype(i)
	flag(j)=   flag(i)
	q(j)=      q(i)
	x(j)=      x(i)
	y(j)=      y(i)
	z(j)=      z(i)
	bfact(j)=      bfact(i)
	do kk= 1,6
	  aniso(kk,j)= aniso(kk,i)
	end do
	do im= 1,nimage
	  xi(j,im)= xi(i,im)
	  yi(j,im)= yi(i,im)
	  zi(j,im)= zi(i,im)
	  hissel(j,im)= hissel(i,im)
	end do
	ncon(j)=   ncon(i)
	do k= 1,ncon(j)
	  icon(k,j)= icon(k,i)
	  scon(k,j)= scon(k,i)
	end do
      end do

      end

c     *********************************************************************

      subroutine hdel (idel,ndel)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

c     **** subroutine to delete atoms ****

	n= ndel
        natom= natom - n

c	**** redo indices for this mol ****
        ires= ares(idel)
        imol= amol(idel)
	write (6,865) n,idel,ires,imol,resname(ires,imol),atname(idel)
865	format (' del',i7,' atoms, start at',i7,' in res/mol ',
     $		2i5,1x,a,1x,a)
        jares(ires,imol)= jares(ires,imol) - n
        nares(ires,imol)= nares(ires,imol) - n
        molend(imol)= molend(imol) - n
	do i= idel,idel+n-1
          qres(ires,imol)= qres(ires,imol) - q(i)
          qtot= qtot - q(i)
	end do
	do i= ires+1,nres(imol)
          iares(i,imol)= iares(i,imol) - n
          jares(i,imol)= jares(i,imol) - n
	end do

c	**** move coords ****
	call hcopy (idel+n,idel,natom-idel+1)

c	**** renumber molecule ****

	do i= molind(imol),molend(imol)
	  atnumb(i)= i - molind(imol) + 1
	end do

c       **** rearrange connection list of this molecule ****

	idel0= idel - molind(imol) + 1
	ideln= idel - molind(imol) + n
        do i= molind(imol),molend(imol)
	  k= 1
	  do while (k.le.ncon(i))

	    if (icon(k,i).ge.idel0 .and. icon(k,i).le.ideln) then
c	      **** delete entry ****
	      ncon(i)= ncon(i) - 1
	      do k1= k,ncon(i)
C		icon(k1,i)= icon(k1,i+1)
C		scon(k1,i)= scon(k1,i+1)
		icon(k1,i)= icon(k1+1,i)
		scon(k1,i)= scon(k1+1,i)
	      end do

	    else if (icon(k,i).gt.ideln) then
c	      **** modify atom nber ****
	      icon(k,i)= icon(k,i) - n
	      k= k + 1

	    else
c	      **** no action needed ****
	      k= k + 1
	    end if

	  end do
        end do

c       **** redo indices for later molecules ****

	do jmol= imol+1,nmol
	  molind(jmol)= molind(jmol) - n
	  molend(jmol)= molend(jmol) - n
	  do i= 1,nres(jmol)
            iares(i,jmol)= iares(i,jmol) - n
            jares(i,jmol)= jares(i,jmol) - n
	  end do
        end do

      return

      end

c     *********************************************************************

      subroutine hdelmol (imol)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

c     **** deletes molecule imol ****

C      write (6,*) 'hdelmol',imol

      if (imol.gt.nmol) return

      nadel= molend(imol) - molind(imol) + 1
      call hdel (molind(imol),nadel)
      do i= imol+1,nmol
	molind(i-1)= molind(i)
	molend(i-1)= molend(i)
	nres(i-1)= nres(i)
	do k= 1,nres(i)
	  resname(k,i-1)= resname(k,i)
	  qres(k,i-1)= qres(k,i)
	  iares(k,i-1)= iares(k,i)
	  jares(k,i-1)= jares(k,i)
	  nares(k,i-1)= nares(k,i)
	end do
        do j= molind(i),molend(i)
	  amol(j)= amol(j) - 1
        end do
      end do

      nmol= nmol - 1

      return

      end

c     *********************************************************************

      subroutine haddmol (imol)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

c     **** adds a molecule before locn imol, without atoms ****

C      write (6,*) 'haddmol',imol,nmol,mmol

      if (imol.gt.nmol+1) stop 'must add mol connected to others'
      if (nmol.ge.mmol) stop 'cant add another molecule'

      do i= nmol,imol,-1
	molind(i+1)= molind(i)
	molend(i+1)= molend(i)
	nres(i+1)= nres(i)
	do k= 1,nres(i)
	  resname(k,i+1)= resname(k,i)
	  qres(k,i+1)= qres(k,i)
	  iares(k,i+1)= iares(k,i)
	  jares(k,i+1)= jares(k,i)
	  nares(k,i+1)= nares(k,i)
	end do
        do j= molind(i),molend(i)
	  amol(j)= amol(j) + 1
        end do
      end do

c     **** add new molecule, no atoms ****

      if (imol.eq.nmol+1) then
	molind(imol)= natom+1
	amol(natom+1)= imol
	ares(natom+1)= 1
      end if
      molend(imol)= molend(imol-1)
      nres(imol)= 1
      resname(1,imol)= ' '
      qres(1,imol)= 0.D0
      iares(1,imol)= molind(imol)
      jares(1,imol)= molind(imol)-1
      nares(1,imol)= 0

      nmol= nmol + 1

      return

      end

c     *********************************************************************

      subroutine hmermol (ifrom,ito)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

c     **** merges mol ifrom into mol ito at end ****

      write (6,*) 'hmermol',ifrom,ito
      write (6,*) 'nber atoms at start=',natom

      naadd= molend(ifrom) - molind(ifrom) + 1
      inew= molend(ito) + 1
      iresold= nres(ito)
      i0= molend(ito)-molind(ito)+1
      ia= jares(nres(ito),ito)
      na= nares(nres(ito),ito)
      call hadd (molend(ito),naadd,1)
      jares(nres(ito),ito) = ia
      nares(nres(ito),ito) = na
      call hcopy (molind(ifrom),inew,naadd)
      write (6,*) 'inew=',inew,molend(ito),i0,iresold
      do i= inew,molend(ito)
	amol(i)= ito
	ares(i)= iresold + nres(i)
	atnumb(i)= atnumb(i) + i0
	do j= 1,ncon(i)
	  icon(j,i)= icon(j,i) + i0
	end do
C	write (6,'(5i7)') i,atnumb(i),ncon(i),icon(1,i)
      end do
      nres(ito)= iresold + nres(ifrom)
      do k= 1,nres(ifrom)
	resname(iresold+k,ito)= resname(k,ifrom)
	qres(iresold+k,ito)= qres(k,ifrom)
	iares(iresold+k,ito)= jares(iresold+k-1,ito) + 1
	jares(iresold+k,ito)= jares(iresold+k-1,ito) + nares(k,ifrom)
	nares(iresold+k,ito)= nares(k,ifrom)
	write (6,*) resname(iresold+k,ito),
     $		iares(iresold+k,ito),jares(iresold+k,ito)
      end do

c     **** del atoms but keep molecule place ****
      call hdel (molind(ifrom),molend(ifrom)-molind(ifrom)+1)
      nres(ifrom)= 1
      resname(1,ifrom)= 'xxx 001 - -'

      write (6,*) 'nber atoms at end=',natom

      return

      end


c     *********************************************************************

      subroutine hcopymol (ifrom,ito)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

c     **** replaces mol ito with a copy of ifrom ****

C      write (6,*) 'hcopymol',ifrom,ito

      call hdelmol (ito)
      call haddmol (ito)
      naadd= molend(ifrom) - molind(ifrom) + 1
      call hadd (molind(ito),naadd,0)
      call hcopy (molind(ifrom),molind(ito),naadd)
      do i= molind(ito),molend(ito)
	amol(i)= ito
      end do
      nres(ito)= nres(ifrom)
      do k= 1,nres(ifrom)
	resname(k,ito)= resname(k,ifrom)
	qres(k,ito)= qres(k,ifrom)
	iares(k,ito)= iares(k,ifrom) + molind(ito) - molind(ifrom)
	jares(k,ito)= jares(k,ifrom) + molind(ito) - molind(ifrom)
	nares(k,ito)= nares(k,ifrom)
      end do

      return

      end

c     *********************************************************************

      subroutine hmovmol (ifrom,ito)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

c     **** moves a molecule from locn ifrom to locn ito ****

      write (6,*) 'moving molecule ',ifrom,ito

c     **** first, move ifrom to end ****
      call hcopymol (ifrom,nmol+1)

c     **** move molecules in the middle ****
      if (ito.gt.ifrom) then
	do i= ifrom+1,ito
	  call hcopymol (i-1,i)
	end do
      else
	do i= ifrom-1,ito,-1
	  call hcopymol (i,i+1)
	end do
      end if

c     **** replace final ****
      call hcopymol (nmol,ito)
      call hdelmol(nmol)

      call writehin ('jnk4.hin',0)
      stop 'first move not debugged'

      end

c     *********************************************************************

      subroutine hadd (ist0,n,iplace)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

c     **** subroutine to add atoms, not including specific properties ****
c     **** creates space for n atoms 
c     **** iplace=1 after atom ist0 ****
c     **** iplace=0 before atom ist0 ****  (use this after a haddmol call)
	write (6,*) 'HADD:',iplace,amol(ist0),amol(ist0-1)

	if (iplace.eq.0 .and. ist0.gt.natom) then
c	  **** just called "haddmol (nmol+1)" to put new mol at end ****
	  ires= 1
	  imol= nmol
          ist= ist0 - 1
	else if (iplace.eq.0 .and. amol(ist0).eq.amol(ist0-1)+2) then
c	  **** just called "haddmol (imol)" to put new mol before imol ****
	  ires= 1
	  imol= amol(ist0) - 1
          ist= ist0 - 1
	else if (iplace.eq.1) then
          ires= ares(ist0)
          imol= amol(ist0)
          ist= ist0
	else
          ires= ares(ist0)
          imol= amol(ist0)
          ist= ist0 - 1
	end if

        write (6,*) 'hadd: after ',ist,' nber=',n,' res mol=',ires,imol

        natom= natom + n
	if (natom.gt.m) stop 'not enough space for new atoms'

c       **** redo indices for later molecules ****

	do jmol= nmol,imol+1,-1
	  molind(jmol)= molind(jmol) + n
	  molend(jmol)= molend(jmol) + n
	  do i= 1,nres(jmol)
            iares(i,jmol)= iares(i,jmol) + n
            jares(i,jmol)= jares(i,jmol) + n
	  end do
        end do

c	**** redo indices for this mol ****

        jares(ires,imol)= jares(ires,imol) + n
        nares(ires,imol)= nares(ires,imol) + n
        molend(imol)= molend(imol) + n
	do i= ires+1,nres(imol)
          iares(i,imol)= iares(i,imol) + n
          jares(i,imol)= jares(i,imol) + n
	end do

c	**** shift later atoms ****
	call hcopy (ist+1,ist+1+n,natom-ist-n)

c	**** some new atom indices ****
	do i= ist+1,ist+n
	  amol(i)= imol
	  ares(i)= ires
	  ncon(i)= 0
	  atname(i)= '?????'
	  ats(i)= '?'
	  attype(i)= '?'
	  atnumb(i)= i - molind(imol) + 1
	  nat(i)= 0
	  q(i)= 0.
	  x(i)= -1000.
	  y(i)= -1000.
	  z(i)= -1000.
	  bfact(i)= 0.
	  do kk= 1,6
	    aniso(kk,i)= 0
	  end do
	  do im= 1,nimage
	    xi(i,im)= x(i)
	    yi(i,im)= y(i)
	    zi(i,im)= z(i)
	    hissel(i,im)= .false.
	  end do
	end do

c	**** renumber this molecule ****

	i0= molind(imol) - 1
	do i= molind(imol),molend(imol)
	  atnumb(i)= i - i0
	end do

c       **** rearrange connection list of this molecule ****

        do i= molind(imol),molend(imol)
	  k= 0
	  do while (k.lt.ncon(i))
	    k= k + 1
	    if (icon(k,i).gt.ist-i0) then
c	      **** modify atom nber ****
	      icon(k,i)= icon(k,i) + n
	    end if
	  end do
        end do

      return

      end

c     *********************************************************************

      subroutine hbadd (i,j,type)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      character*1 type

C      write (6,*) 'add bond',i,j,' ',type

c     **** adds a bond between atoms i and j ****

      i0= molind(amol(i)) - 1

      do k= 1,ncon(i)
	if (icon(k,i)+i0 .eq. j) return
      end do

      ncon(i)= ncon(i) + 1
      if (ncon(i).gt.mcon) stop 'hbadd: too many bonds'
      icon(ncon(i),i)= j-i0
      scon(ncon(i),i)= type

      ncon(j)= ncon(j) + 1
      if (ncon(j).gt.mcon) stop 'hbadd: too many bonds'
      icon(ncon(j),j)= i-i0
      scon(ncon(j),j)= type

      return
      end

c     *********************************************************************

      subroutine hbdel (i,j)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

c     **** deletes a bond between atoms i and j ****

      i0= molind(amol(i)) - 1
      k= 1
      do while (k.le.ncon(i) .and. icon(k,i).ne.j-i0)
	k= k + 1
      end do
      if (k.le.ncon(i)) then
	  do k1= k+1,ncon(i)
	    icon(k1-1,i)= icon(k1,i)
	    scon(k1-1,i)= scon(k1,i)
	  end do
	  ncon(i)= ncon(i) - 1
      end if

      do k= 1,ncon(j)
	if (icon(k,j).eq.i-i0) then
	  do k1= k+1,ncon(j)
	    icon(k1-1,j)= icon(k1,j)
	    scon(k1-1,j)= scon(k1,j)
	  end do
	  ncon(j)= ncon(j) - 1
	end if
      end do

      return
      end

c     *************************************************************************

      subroutine hadda2 (i1o,io,ico,xo,yo,zo,ncono,icono,atnameo,i,ic)

c     **** adds atom i connected to atom ic based on geom of a reference system ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      real*8 xo(*),yo(*),zo(*),rot(3,3)
      integer ncono(*),icono(mcon,*)
      character*5 atnameo(*)
      logical found


c     **** find two more atoms to get angle and torsion from ****

      found= .false.
      ijox= 1
      i1= molind(amol(ic)) - 1
      write (6,*) 'hadda2 ref=',i1o,io,ico
      write (6,*) 'hadda2 new=',i1,i,ic

      do while (.not.found .and. ijox.le.ncono(ico))
c	**** ijo= atom con to ico, ij same atom con to ic ****
	ijo= icono(ijox,ico) + i1o
	if (ijo.ne.0 .and. ijo.ne.io) then
          write (6,'(a,3i6,1x,a)') 'search for first def atom',
     $		ncono(ico),ijox,ijo,atnameo(ijo)
c	  **** search for same atom type in new coords ****
	  ijj= 0
	  do ijx= 1,ncon(ic)
	    ij= icon(ijx,ic) + i1
	    if (atname(ij).eq.atnameo(ijo)) ijj= ij
	    write (6,*) ijx,icon(ijx,ic),i1,ij,atname(ij),ijj
	  end do

	  if (ijj.ne.0) then
	    ij= ijj
	    write (6,*) 'match for atom ',atname(ij),atnameo(ijo)

c	    **** search for a second atom on ico and ic ****

	    ikox= ijox + 1
            do while (.not.found .and. ikox.le.ncono(ico))
	      iko= icono(ikox,ico) + i1o
	      if (iko.gt.0) then
              write (6,'(a,3i6,1x,a)') 'search for sec default atom',
     $			ncono(ico),ikox,iko,atnameo(iko)
	      if (iko.ne.ico) then
c	        **** search for same atom type in new coords ****
	        ikk= 0
	        do ikx= 1,ncon(ic)
	          ik= icon(ikx,ic) + i1
	          if (atname(ik).eq.atnameo(iko)) ikk= ik
		  write (6,*) ikx,ik,ikk,atname(ik),atnameo(iko)
	        end do
		ik= ikk
		found= ik .ne. 0
	      end if
	      end if
	      if (.not.found) ikox= ikox + 1
	    end do

	    if (.not.found) then

c	      **** search for connected atom to ijo  ****
	      ikox= 1
              do while (.not.found .and. ikox.le.ncono(ijo))
	        iko= icono(ikox,ijo) + i1o
	        if (iko.gt.0) then
                write (6,*) 'sec',ncono(ijo),ikox,iko
	        if (iko.ne.ico) then
c	          **** search for same atom type in new coords ****
	          ikk= 0
	          do ikx= 1,ncon(ij)
	            ik= icon(ikx,ij) + i1
	            if (atname(ik).eq.atnameo(iko)) ikk= ik
	          end do
		  ik= ikk
		  found= ik .ne. 0
	        end if
	        end if
	        if (.not.found) ikox= ikox + 1
	      end do
	    end if

	  end if
	end if
	if (.not.found) ijox= ijox + 1
      end do

      if (.not.found) then
	call writehin ('err.hin',0)
	stop 'hadda2: cant form z matrix'
      end if
      write (6,'(a,4i6,a,4i6)') 'replicating ',io,ico,ijo,iko,
     $		' to ',i,ic,ij,ik

c     **** get bond length, angle and torsion from ref coords ****

      jo= ico
      ko= ijo
      lo= iko

      xji= xo(io) - xo(jo)
      yji= yo(io) - yo(jo)
      zji= zo(io) - zo(jo)
      rji= sqrt ( xji**2 + yji**2 + zji**2 )
      xji= xji / rji
      yji= yji / rji
      zji= zji / rji

      xjk= xo(ko) - xo(jo)
      yjk= yo(ko) - yo(jo)
      zjk= zo(ko) - zo(jo)
      rjk= sqrt ( xjk**2 + yjk**2 + zjk**2 )
      xjk= xjk / rjk
      yjk= yjk / rjk
      zjk= zjk / rjk
      cosphi= xji*xjk + yji*yjk + zji*zjk
      sinphi= sqrt (1.D0 - cosphi**2)
      phi= acos (cosphi)

      xkl= xo(lo) - xo(ko)
      ykl= yo(lo) - yo(ko)
      zkl= zo(lo) - zo(ko)
      rkl= sqrt ( xkl**2 + ykl**2 + zkl**2 )
      xkl= xkl / rkl
      ykl= ykl / rkl
      zkl= zkl / rkl
      cosph3= - ( xjk*xkl + yjk*ykl + zjk*zkl )
      sinph3= sqrt (1.D0 - cosph3**2)
	    
      xijjk= - yji*zjk + zji*yjk
      yijjk= - zji*xjk + xji*zjk
      zijjk= - xji*yjk + yji*xjk
      xjkkl=   yjk*zkl - zjk*ykl
      yjkkl=   zjk*xkl - xjk*zkl
      zjkkl=   xjk*ykl - yjk*xkl
      costor= (xijjk*xjkkl + yijjk*yjkkl + zijjk*zjkkl )
     $			/ sinphi / sinph3
      if (costor.gt. 1.D0 .and. costor.lt. 1.00001) 
     $			  costor= 1.D0
      if (costor.lt.-1.D0 .and. costor.gt.-1.00001)
     $			  costor=-1.D0
      tor= acos (costor)
      torsign= ( yijjk*zjkkl - zijjk*yjkkl ) * xjk
     $		     + ( zijjk*xjkkl - xijjk*zjkkl ) * yjk
     $		     + ( xijjk*yjkkl - yijjk*xjkkl ) * zjk
      if (torsign.lt.0.D0) tor= - tor
      write (4,'(2a,3f8.3)') atname(ic),' tor,costor,torsign=',tor*pre,costor,torsign
      write (6,'(2a,3f8.3)') atname(ic),' tor,costor,torsign=',tor*pre,costor,torsign

c     **** apply these to new coords ****

      call haddzrad (i,ic,ij,ik,rji,phi,tor)
      atname(i)= atnameo(io)

      return
      end

c     *************************************************************************

      subroutine haddz (i,j,k,l,rji,phi,tor)

c     **** adds atom in z-matrix coords ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      call haddzrad (i,j,k,l,rji,phi/pre,tor/pre)
      end

      subroutine haddzrad (i,j,k,l,rji,phi,tor)

c     **** adds atom in z-matrix coords (angles in radians) ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      real*8 rot(3,3)

	write (6,*) 'haddzrad',i,j,k,l,rji,phi,tor

C	**** rotate JK bond to z axis and KL to in xz plane ****
	xkl= x(l) - x(k)
	ykl= y(l) - y(k)
	zkl= z(l) - z(k)
	rkl= sqrt ( xkl**2 + ykl**2 + zkl**2 )
	xkl= xkl / rkl
	ykl= ykl / rkl
	zkl= zkl / rkl
	xjk= x(k) - x(j)
	yjk= y(k) - y(j)
	zjk= z(k) - z(j)
	rjk= sqrt ( xjk**2 + yjk**2 + zjk**2 )
	xjk= xjk / rjk
	yjk= yjk / rjk
	zjk= zjk / rjk

	rot(3,1)= - xjk
	rot(3,2)= - yjk
	rot(3,3)= - zjk

	cosjkl= - (xjk*xkl + yjk*ykl + zjk*zkl)
	rot(1,1)= xkl + cosjkl*xjk
	rot(1,2)= ykl + cosjkl*yjk
	rot(1,3)= zkl + cosjkl*zjk
	xxx= sqrt (rot(1,1)**2 + rot(1,2)**2 + rot(1,3)**2)
	rot(1,1)= rot(1,1) / xxx
	rot(1,2)= rot(1,2) / xxx
	rot(1,3)= rot(1,3) / xxx

c	**** 3rd vector (y) is cross prod of first two *****
	rot(2,1)= rot(3,2)*rot(1,3) - rot(3,3)*rot(1,2)
	rot(2,2)= rot(3,3)*rot(1,1) - rot(3,1)*rot(1,3)
	rot(2,3)= rot(3,1)*rot(1,2) - rot(3,2)*rot(1,1)

c	**** coords in rotated frame ****
	xxx= rji * sin(phi) * cos(tor)
	yyy= rji * sin(phi) * sin(tor)
	zzz= - rji * cos(phi)

c	**** back transform ****
	x(i)= x(j) + rot(1,1)*xxx + rot(2,1)*yyy + rot(3,1)*zzz
	y(i)= y(j) + rot(1,2)*xxx + rot(2,2)*yyy + rot(3,2)*zzz
	z(i)= z(j) + rot(1,3)*xxx + rot(2,3)*yyy + rot(3,3)*zzz

	bfact(i)= 0.d0
	do kk= 1,6
	  aniso(kk,i)= 0
	end do
	call hbadd (i,j,'s')

	write (6,'(a,4i6,3f8.3,a)') 'z-mat added atom at ',
     $		i,j,k,l,x(i),y(i),z(i),atname(i)
	write (6,*) atname(j),atname(k),atname(l)

      end

c     *************************************************************************

      subroutine hadda1 (n,ic,ihybx,r,bang)

c     **** adds coords for one or more atoms in a residue ****
c     **** added after atom ninit, connected to atom ic ****
c     **** n points to last added on exit ****
c     **** use negative value for ihyb to delet H atoms before adding ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      parameter (mclose=100)
      integer iclose(mclose)
      
      write (6,863) n,ic,ihybx,r,bang
863   format (' adding atom after loc',i7,' to atom',i7,' hyb=',i2,
     $		' r=',f6.3,' mult angle=',f8.2)

      ihyb= abs (ihybx)

      if (ihybx.lt.0) then
c	**** delete bond list for any connected H atoms ****
	do jj= ncon(ic),1,-1
	  j= icon(jj,ic) + molind(amol(ic)) - 1
	  if (nat(j).eq.1) then
	    ncon(j)= 0
	    ncon(ic)= ncon(ic) - 1
	    do kk= jj,ncon(ic)
	      icon(kk,ic)= icon(kk+1,ic)
	      scon(kk,ic)= scon(kk+1,ic)
	    end do
	  end if
	end do
      end if

      if (ihyb.gt.3) then
	write (6,*) 'ERROR: SP 1,2,3 only allowed, sp=',ihyb
	stop 'addatm1: valence'
      end if

      ninit= n
      imol= amol(ic)
      ires= ares(ic)
      ic0= molind(imol) - 1

c     **** determine bonds already on center atom ic ****
      nb= ncon(ic)
      nlp= 0
      if (nat(ic).eq.8 .or. nat(ic).eq.16) nlp= 2
      nb= nb + nlp
      write (6,*) '# bonds + (O,S lp) already on atom=',nb
      if (nb.eq.0 .or. nb.gt.ihyb+1) then
	write (6,*) 'ERROR: number of bond on atom ',ic, ' is',nb
	stop 'addatm1: nber bonds'
      end if
      ij= icon(1,ic) + ic0
      ik= icon(2,ic) + ic0
      il= icon(3,ic) + ic0

      if (ihyb+1.eq.nb) then
c	**** full compliment on this atom, do nothing ****
	new= 0

      else if (ihyb.eq.1) then
c	**** SP1 (linear), add one atom ****
	new= 1
	n= n + 1
	r1= r / sqrt ( hbl2(ij,ic) )
	x(n)= x(ic) + (x(ic) - x(ij)) * r1
	y(n)= y(ic) + (y(ic) - y(ij)) * r1
	z(n)= z(ic) + (z(ic) - z(ij)) * r1
	call hcopyim (n,1)

      else if (ihyb.eq.2) then
	if (nb.gt.1) then
c	  **** SP2 (planar triangular), two atoms found, add one atom ****
	  call sp2a (n,r,ic,ij,ik)
	  new= 1
	else
c	  **** SP2 (planar triangular), one atom found, add two atoms ****
c	  **** determine bonds already on center atom ii ****
	  ii= ij
	  ij= icon(1,ii) + ic0
	  ik= icon(2,ii) + ic0
	  call sp2b (n,r,ic,ii,ij,ik,bang)
	  new= 2
	end if

      else if (nb-nlp.eq.1) then
c	**** SP3 (tetrahedral), add 1 (if Oxy) or else 3 atoms ****
	call tetra3 (n,new,r,ic,ij,bang,nlp)

      else if (nb.eq.2) then
c	**** SP3 (tetrahedral), add 2 atoms ****
	call tetra2 (n,new,r,ic,ij,ik,bang)

      else
c	**** SP3 (tetrahedral), add 1 atom at average of 3 poss posns ****
	call tetra1 (n,new,r,ic,ij,ik,il)

      end if

c     **** add bonds ****

      j0= molind(imol) - 1
      do j= ninit+1,n
	ncon(j)= 0
	call hbadd (j,ic,'s')
        if (ihybx.gt.0) then
c	  **** only insert default atomic details if not replacing ****
	  atname(j)= atname(ic)
	  k= 1
	  do while (k.lt.5 .and. atname(j)(k:k).eq.' ')
	    k= k + 1
	  end do
	  if (nat(ic).eq.7) then
	    attype(j)= ' H'
	  else
	    attype(j)= 'H'//atname(j)(k:k)
	  end if
	  atname(j)(k:k)= 'H'
	  if (ninit+1.ne.n) 
     $	    write (atname(j)(k-1:k-1),'(i1)') j-ninit
	  ats(j)= 'H'
	  nat(j)= 1
	  q(j)= 0.0
	  bfact(j)= 0.0
	  do kk= 1,6
	    aniso(kk,j)= 0
	  end do
	  write (6,'(i5,1x,a,3f8.3)') j,atname(j),x(j),y(j),z(j)
	end if
      end do

c     **** if added 3 then rotate to minimize interactions ****
      if (n-ninit.eq.3) then
        write (6,*) 'rotating CH3'
	nclose= 0
	do k= 1,natom
	  rr= hbl2(k,j)
	  if (rr.lt.5.0) then
	    nclose= nclose + 1
	    if (nclose.gt.mclose) stop 'HADDA1: MCLOSE'
	    iclose(nclose)= k
	  end if
	end do
	rrmin= 1.e30
	do irot= 0,36
	  ang= irot * 10.
	  call hrotate (ic,ic0+icon(1,ic),ang)
	  do kclose= 1,nclose
	    k= iclose(kclose)
	    do j= ninit+1,n
	      rr= hbl2(k,j)
	      if (rr.lt.rrmin) then
		rrmin= rr
		angmin= ang
	      end if
	    end do
	  end do
	end do
	call hrotate (ic,ic0+icon(1,ic),angmin)
      end if

c     **** check distances to other atoms ****

      do j= ninit+1,n
	do k= 1,natom
	  rr= hbl2(k,j)
	  if (k.ne.j .and. k.ne.ic .and.
     $		  (k.le.ninit .or. k.gt.n) .and. rr.lt.3.0) then
	    isamec= 0
	    do k1= 1,ncon(ic)
	      if (j0+icon(k1,ic).eq.k) isamec= 1
	    end do
	    rr= sqrt (rr)
	    if (nat(k).le.1 .and. rr.gt.1.4) isamec= 1
	    if (isamec.eq.0 .or. rr.lt.1.3) then 
	      write (4,'(a,i7,1x,a,f6.3,i7,1x,a,1x,a)') ' SHORT INT:',j,
     $	      resname(ires,imol),rr,k,atname(k),resname(ares(k),amol(k))
	      write (6,'(a,i7,1x,a,f6.3,i7,1x,a,1x,a)') ' SHORT INT:',j,
     $	      resname(ires,imol),rr,k,atname(k),resname(ares(k),amol(k))
C	      stop 'SHORT'
	    end if
	  end if
	end do
      end do

      return

      end

c     ******************************************************************

      subroutine sp2a	(n,r,ic,ij,ik)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      real*8		r,a(3),ra

c     **** add atom along the bisector of cj and ck vectors ****

C     write (6,*) 'called SP2A'

      a(1)= (x(ij) + x(ik)) * 0.5 - x(ic)
      a(2)= (y(ij) + y(ik)) * 0.5 - y(ic)
      a(3)= (z(ij) + z(ik)) * 0.5 - z(ic)
      ra= sqrt (a(1)**2 + a(2)**2 + a(3)**2)
      n= n + 1
      x(n)= x(ic) - a(1) / ra * r
      y(n)= y(ic) - a(2) / ra * r
      z(n)= z(ic) - a(3) / ra * r
      call hcopyim (n,1)

      return
      end

c     ******************************************************************

      subroutine sp2b	(n,r,ic,ii,ij,ik,bang)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      real*8		a(3),b(3),r,ra,rb,bang,ang1,cang,sang,dot

C      write (6,*) 'called SP2B'

c     **** calc ci vector ****
      a(1)= x(ii) - x(ic)
      a(2)= y(ii) - y(ic)
      a(3)= z(ii) - z(ic)
      ra= sqrt (a(1)**2 + a(2)**2 + a(3)**2)

c     **** calc jk vector ****
      b(1)= x(ij) - x(ik)
      b(2)= y(ij) - y(ik)
      b(3)= z(ij) - z(ik)

c     **** orthogonalize jk wrt ci ****
      dot= (a(1)*b(1) + a(2)*b(2) + a(3)*b(3)) / ra**2
      b(1)= b(1) - dot*a(1)
      b(2)= b(2) - dot*a(2)
      b(3)= b(3) - dot*a(3)
      rb= sqrt (b(1)**2 + b(2)**2 + b(3)**2)

c     **** add two geminal hydrogens ****
      ang1= bang * 3.14159265 / 360.
      cang= cos (ang1) * r / ra
      sang= sin (ang1) * r / rb
      n= n + 2
      x(n-1)= x(ic) - a(1)*cang + b(1)*sang
      y(n-1)= y(ic) - a(2)*cang + b(2)*sang
      z(n-1)= z(ic) - a(3)*cang + b(3)*sang
      x(n  )= x(ic) - a(1)*cang - b(1)*sang
      y(n  )= y(ic) - a(2)*cang - b(2)*sang
      z(n  )= z(ic) - a(3)*cang - b(3)*sang
      call hcopyim (n-1,2)

      return
      end

c     ******************************************************************

      subroutine tetra1 (n,new,r,ic,ij,ik,il)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      real*8		r,a(3),b(3),ra,rb

C      write (6,*) 'called TETRA1'

c     **** sum the vector from ic to the 3 atoms ****

      a(1)= x(ij) - x(ic)
      a(2)= y(ij) - y(ic)
      a(3)= z(ij) - z(ic)
      ra= sqrt (a(1)**2 + a(2)**2 + a(3)**2)
      b(1)= a(1) / ra
      b(2)= a(2) / ra
      b(3)= a(3) / ra

      a(1)= x(ik) - x(ic)
      a(2)= y(ik) - y(ic)
      a(3)= z(ik) - z(ic)
      ra= sqrt (a(1)**2 + a(2)**2 + a(3)**2)
      b(1)= b(1) + a(1) / ra
      b(2)= b(2) + a(2) / ra
      b(3)= b(3) + a(3) / ra

      a(1)= x(il) - x(ic)
      a(2)= y(il) - y(ic)
      a(3)= z(il) - z(ic)
      ra= sqrt (a(1)**2 + a(2)**2 + a(3)**2)
      b(1)= b(1) + a(1) / ra
      b(2)= b(2) + a(2) / ra
      b(3)= b(3) + a(3) / ra

      rb= sqrt (b(1)**2 + b(2)**2 + b(3)**2)
      x(n+1)= x(ic) - b(1) * r/rb
      y(n+1)= y(ic) - b(2) * r/rb
      z(n+1)= z(ic) - b(3) * r/rb
      call hcopyim (n+1,1)
      n= n + 1
      new= 1
      
      return
      end

c     ******************************************************************

      subroutine tetra2 (n,new,r,ic,ij,ik,bang)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      real*8		r,a(3),b(3),c(3),ra,rb,rc,fa,fc,bang

      write (6,*) 'called TETRA2'

c     **** calc ik-ic-ij bisector and normal to ic-ij-ik plane ***
      a(1)= x(ij) - x(ic)
      a(2)= y(ij) - y(ic)
      a(3)= z(ij) - z(ic)
      ra= sqrt (a(1)**2 + a(2)**2 + a(3)**2)
      b(1)= x(ik) - x(ic)
      b(2)= y(ik) - y(ic)
      b(3)= z(ik) - z(ic)
      rb= sqrt (b(1)**2 + b(2)**2 + b(3)**2)
      a(1)= a(1)/ra + b(1)/rb
      a(2)= a(2)/ra + b(2)/rb
      a(3)= a(3)/ra + b(3)/rb
      ra= sqrt (a(1)**2 + a(2)**2 + a(3)**2)
      c(1)= a(2)*b(3) - a(3)*b(2)
      c(2)= a(3)*b(1) - a(1)*b(3)
      c(3)= a(1)*b(2) - a(2)*b(1)
      rc= sqrt (c(1)**2 + c(2)**2 + c(3)**2)

c     **** add two hydrogens ** new h-c-h angle is bang ****
      fa= -cos (bang/2.*3.14159265/180.) /ra*r
      fc=  sin (bang/2.*3.14159265/180.) /rc*r
      x(n+1)= x(ic) + a(1)*fa + c(1)*fc
      y(n+1)= y(ic) + a(2)*fa + c(2)*fc
      z(n+1)= z(ic) + a(3)*fa + c(3)*fc
      x(n+2)= x(ic) + a(1)*fa - c(1)*fc
      y(n+2)= y(ic) + a(2)*fa - c(2)*fc
      z(n+2)= z(ic) + a(3)*fa - c(3)*fc
      call hcopyim (n+1,2)
      n= n + 2
      new= 2

      return
      end

c     ******************************************************************

      subroutine tetra3 (n,new,r,ic,ij,bang,nlp)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      real*8		rot(3,3),t(3,3),r,ra,sinth,costh,sinph,
     $			cosph,bang

C      write (6,*) 'called TETRA3',nlp

c     **** if nlp=2 then adds just 1 H, else adds 3 H ****

c     **** standard tetrahedronon with variable angle ****
      t(1,1)= sqrt (2.D0*(1.D0-cos(bang*3.14159265/180.))/3.D0) * r
      t(1,2)= 0.D0
      t(1,3)= - sqrt (r**2 - t(1,1)**2)
      t(2,1)= -t(1,1) / 2.D0
      t(2,2)= sqrt(3.D0) * 0.5D0 * t(1,1)
      t(2,3)= t(1,3)
      t(3,1)= t(2,1)
      t(3,2)= -t(2,2)
      t(3,3)= t(2,3)

c     **** calc given bond vector and a rotn matrix from (0,0,1) to this ****

      rot(3,1)= x(ij) - x(ic)
      rot(3,2)= y(ij) - y(ic)
      rot(3,3)= z(ij) - z(ic)
      ra= sqrt (rot(3,1)**2 + rot(3,2)**2 + rot(3,3)**2)
      rot(3,1)= rot(3,1) / ra
      rot(3,2)= rot(3,2) / ra
      rot(3,3)= rot(3,3) / ra
      sinth= sqrt (rot(3,1)**2 + rot(3,2)**2)
      costh= rot(3,3)
      if (sinth.gt.1.e-4) then
	sinph= rot(3,2) / sinth
	cosph= rot(3,1) / sinth
      else
	sinph= 0.D0
	cosph= 1.D0
      end if
      rot(1,1)= costh*cosph
      rot(1,2)= costh*sinph
      rot(1,3)= -sinth
      rot(2,1)= -sinph
      rot(2,2)= cosph
      rot(2,3)= 0.D0
c     write (6,'(1x,3f10.5)') ((rot(i,j),j=1,3),i=1,3)

c     **** rotn of standard tetrahedron ****
      do i= 1,3-nlp
	n= n + 1
	x(n)= x(ic) + rot(1,1)*t(i,1) + rot(2,1)*t(i,2) +rot(3,1)*t(i,3)
	y(n)= y(ic) + rot(1,2)*t(i,1) + rot(2,2)*t(i,2) +rot(3,2)*t(i,3)
	z(n)= z(ic) + rot(1,3)*t(i,1) + rot(2,3)*t(i,2) +rot(3,3)*t(i,3)
        call hcopyim (n,1)
      end do

      new= 3-nlp
      return
      end

c     ******************************************************************

      function hisc2mg (ires,imol)

c     **** determine if input residue is HIS connected to Mg  ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

      hisc2mg= .false.
      if (resname(ires,imol)(1:3).ne.'His'
     $      .and. resname(ires,imol)(1:3).ne.'HiS' ) return

      stop 'HISc2MG not CODED'

      end

c     ****************************************************************

      subroutine hrotate (rcon,rvec,ang)

c     **** rotates all atoms connected to rcon about vec to rvec ****
c     **** through angle ang ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      parameter (mlist=100)
      integer rcon,rvec,list(mlist)
      real*8 r(3,3),r2(3,3),r3(3,3),xx(3)

      ires= ares(rcon)
      imol= amol(rcon)
      i0= iares(1,imol) - 1

c     **** find list of atoms to rotate ****

      nlist= 1
      list(1)= rcon
      ilist= 1
      do while (ilist.le.nlist)
C	write (6,*) 'checking cons on ',list(ilist),nlist
	do j= 1,ncon(list(ilist))
	  jj= icon(j,list(ilist)) + i0
	  isinc= 0
	  do k= 1,nlist
	    if (list(k).eq.jj) isinc= 1
	  end do
	  if (isinc.eq.0 .and. jj.ne.rvec) then
	    nlist= nlist + 1
	    if (nlist.gt.mlist) stop 'HROTATE: MLIST'
	    list(nlist)= jj
	  end if
	end do
	ilist= ilist + 1
      end do

C      write (6,*) 'Rotate list:',nlist
C      write (6,'(10i7)') (list(i),i=1,nlist)

c	**** vector to rotate about ****

      r(3,1)= x(rcon) - x(rvec)
      r(3,2)= y(rcon) - y(rvec)
      r(3,3)= z(rcon) - z(rvec)

c     **** transformation from vector in r(3,i) to the z axis ****
c     **** determine polar angles theta and phi ****

      st= sqrt ( r(3,1)**2 + r(3,2)**2 + r(3,3)**2)
      r(3,1)= r(3,1) / st
      r(3,2)= r(3,2) / st
      r(3,3)= r(3,3) / st
      ct= r(3,3)
      st= sqrt (1.D0-ct**2)
      if (abs(st).lt.1.e-5) then
	cp= 1.D0
	sp= 0.D0
      else
	cp= r(3,1)/st
	sp= r(3,2)/st
      end if

c     **** rest of rotation matrix of vector to Z ****
      r(1,1)= - sp
      r(1,2)=   cp
      r(1,3)= 0.D0
      r(2,1)= - ct*cp
      r(2,2)= - ct*sp
      r(2,3)=   st
C      write (6,820) 'r',((r(i,j),j=1,3),i=1,3)

c     **** rotation matrix about vector ****
      nowpos= ipos
      st= ang/180.d0*3.14159265d0
      ct= cos (st)
      st= sin (st)
      r2(1,1)= ct
      r2(1,2)= -st
      r2(1,3)= 0.D0
      r2(2,1)= st
      r2(2,2)= ct
      r2(2,3)= 0.D0
      r2(3,1)= 0.D0
      r2(3,2)= 0.D0
      r2(3,3)= 1.D0
C      write (6,820) 'r2',((r2(i,j),j=1,3),i=1,3)

c     **** final rotation matrix ****
      call mmult   (r3,r2,r,3,3)
      call mmultta (r2,r,r3,3,3)
C      write (6,820) 'r2',((r2(i,j),j=1,3),i=1,3)
C      write (6,*) 'nratm=',nratm

c     **** rotate connected atoms ****
      do i1= 1,nlist
	i= list(i1)
	xx(1)= x(i) - x(rcon)
	xx(2)= y(i) - y(rcon)
	xx(3)= z(i) - z(rcon)
C	write (6,'(3f10.4)') xx
	call crotn (xx,1,r2)
C	write (6,'(3f10.4)') xx
	x(i)= xx(1) + x(rcon)
	y(i)= xx(2) + y(rcon)
	z(i)= xx(3) + z(rcon)
        call hcopyim (i,1)
      end do

820   format (1x,a2 / (1x,3f10.6) )

      return
      end

c     *********************************************************************

      subroutine hcopyim (i1,n)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

c     **** subroutine to create images for n atoms starting at i1 ****

      do i= i1,i1+n-1
	do j= 1,nimage
	  xi(i,j)= rotim(1,j)*x(i) + rotim(2,j)*y(i) +
     $		   rotim(3,j)*z(i) + rotim(10,j) 
	  yi(i,j)= rotim(4,j)*x(i) + rotim(5,j)*y(i) +
     $		   rotim(6,j)*z(i) + rotim(11,j) 
	  zi(i,j)= rotim(7,j)*x(i) + rotim(8,j)*y(i) +
     $		   rotim(9,j)*z(i) + rotim(12,j) 
	  if (j.gt.1) then
	    rr= (xi(i,j)-x(i))**2 + (yi(i,j)-y(i))**2 + (zi(i,j)-z(i))**2 
	    if (rr.lt.0.1) then
	      imol= amol(i)
	      ires= ares(i)
	      write (6,*) 'atom related by symmetry:',i,j,ires,imol,
     $		resname(ires,imol),' ',atname(i)
	      xi(i,j)= xi(i,j) + max(alen,blen,clen)
	    end if
	  end if
	end do
      end do

      return
      end

c     *********************************************************************

      subroutine him2real (nmol1,natom1,nimage1)

c     **** make image atoms manifest as real atoms ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

      natom1= natom
      nmol1= nmol
      nimage1= nimage

      nimage= 1
      if (nmol1*nimage1.gt.mmol) then
        write (6,*) 'Error: including images gives',
     $		nmol1*nimage1,' mols'
        stop 'mmol too small for expansion'
      end if
      if (natom1*nimage1.gt.m) then
        write (6,*) 'Error: including images gives',
     $		natom1*nimage1,' atoms'
        stop 'm too small for expansion'
      end if
      do im= 2,nimage1
        j0= natom1*(im-1)
        call hcopy (1,j0+1,natom1)
        do j= 1,natom1
          x(j0+j)= xi(j,im)
          y(j0+j)= yi(j,im)
          z(j0+j)= zi(j,im)
          amol(j0+j)= amol(j) + (im-1)*nmol1
        end do
      end do

      nmol= nmol1 * nimage1
      natom= natom1 * nimage1

      do im= 2,nimage1
        do jmol= 1,nmol1
	  imol= jmol + (im-1)*nmol1
	  nres(imol)= nres(jmol)
	  molind(imol)= molind(jmol) + (im-1)*natom1
	  molend(imol)= molend(jmol) + (im-1)*natom1
	  do ires= 1,nres(imol)
	    resname(ires,imol)= resname(ires,jmol)
	    qres(ires,imol)= qres(ires,jmol)
	    nares(ires,imol)= nares(ires,jmol)
	    iares(ires,imol)= iares(ires,jmol) + (im-1)*natom1
	    jares(ires,imol)= jares(ires,jmol) + (im-1)*natom1
	  end do
	end do
      end do

      call hcopyim (1,natom)
      write (6,*) 'Images expanded as separate molecules, new natom=',natom,' nmol=',nmol

      return
      end

c     *************************************************************************

      subroutine	mmult (c,a,b,m,n)

      real*8		a(m,n),b(m,n),c(m,n),xx
      integer		i,j,k,m,n

c     **** matrix multiply **** c= a x b ****

      do i= 1,n
	do j= 1,n
	  xx= 0.D0
	  do k= 1,n
	    xx= xx + a(i,k) * b(k,j)
	  end do
	  c(i,j)= xx
	end do
      end do

      return
      end

c     *************************************************************************

      subroutine	mmultta (c,a,b,m,n)

      real*8		a(m,n),b(m,n),c(m,n),xx
      integer		i,j,k,m,n

c     **** matrix multiply **** c= a(transp) x b ****

      do i= 1,n
	do j= 1,n
	  xx= 0.D0
	  do k= 1,n
	    xx= xx + a(k,i) * b(k,j)
	  end do
	  c(i,j)= xx
	end do
      end do

      return
      end


c     *************************************************************************

      subroutine	mmulttb (c,a,b,m,n)

      real*8		a(m,n),b(m,n),c(m,n),xx
      integer		i,j,k,m,n

c     **** matrix multiply **** c= a x b(transp) ****

      do i= 1,n
	do j= 1,n
	  xx= 0.D0
	  do k= 1,n
	    xx= xx + a(i,k) * b(j,k)
	  end do
	  c(i,j)= xx
	end do
      end do

      return
      end

c     *************************************************************************

      subroutine	crotn (v,n,r)

      real*8		v(3,n),x,y,z,r(3,3)
      integer j,n

c     **** r is rotation matrix for the coordinates ****

      do j= 1,n
	x= r(1,1)*v(1,j) + r(1,2)*v(2,j) + r(1,3)*v(3,j)
	y= r(2,1)*v(1,j) + r(2,2)*v(2,j) + r(2,3)*v(3,j)
	z= r(3,1)*v(1,j) + r(3,2)*v(2,j) + r(3,3)*v(3,j)
	v(1,j)= x
	v(2,j)= y
	v(3,j)= z
      end do

      return
      end

c     ***********************************************************************

      subroutine upcase (s)

      character*(*) s

      n= len (s)
      do i= 1,n
	if (s(i:i).le.'z' .and. s(i:i).ge.'a') s(i:i)= char (
     $		ichar(s(i:i)) + ichar('A') - ichar ('a') )
      end do

      return
      end

c     ***********************************************************************

      subroutine lowcase (s)

      character*(*) s

      n= len (s)
      do i= 1,n
	if (s(i:i).le.'Z' .and. s(i:i).ge.'A') s(i:i)= char (
     $		ichar(s(i:i)) + ichar('a') - ichar ('A') )
      end do

      return
      end

c     ******************************************************************

      block data

      implicit real*8 (a-h,o-z)
      character*2 atsym(-1:103)
      common /atsymb/ atsym
      common /const/ au,au2ev

      data au,au2ev /0.529177D0,27.2107D0/

c     **** Lp = lone-pair electrons modelled by additional charge ****
      data atsym	/ 'Lp', '  ',
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

      end
