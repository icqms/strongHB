c       **** converts pdb file to hin format, seps mols ****

	implicit real*8 (a-h,o-z)
	include 'readhin.cmn'
	parameter (mstr=20, mimage1=mimage*10)

	real*8 wk(3),boxvecs1(12),rotim1(12,mimage1)
	character*80 fname,c,cnew,cani
	character*12 nresname
	character*15 oldatom(m),oldatom1
	character*7 lastres
	character*4 aan
	character*2 atss
        character*1 s
        character*6 strsel(0:mstr)
        character*(mstr) str,strlist(0:mstr)
	integer imol0(100),ires0(100),istr(mres,mmol),iwk(3),imatom(m)
	integer*2 aniso1(6)
	logical isbond,isok,isaniso
        
        logical inmol
        integer nres1(m) 
 
	open (4,file='pdb2hin-aniso.out',status='unknown')
C	open (8,file='8images.pdb',status='unknown')

        if (nargs().lt.2) then
          write (6,*) 'ENTER root name of .pdb file'
          read (5,'(a)') fname
        else
	  call getarg (1,fname)
        end if
          
	nf= len_trim(fname)+1
	fname(nf:nf+3)= '.pdb'
	open (9,file=fname,status='old')


c	**** first read to get structural options *****

	alen= 100.d0
	blen= 100.d0
	clen= 100.d0
	alpha= 90.d0
	beta= 90.d0
	gamma= 90.d0
	nimage1= 0

	imol= -1
	nmol= 0
	natom= 0
	lastres= ' '
	str= ' '
	istr= 0
	nstr= 0
	strlist(nstr)= ' '
	strsel(nstr)= ' '

100   continue
	read (9,'(a)',end=120) c
	len= len_trim(c)

c	**** read header to get list of molecules in unit cell ****

c12345678901234567890123456789012345678901234567890123456789012345678901234567890
cREMARK 290   SMTRY1   1  1.000000  0.000000  0.000000        0.00000 
	if (c(1:18).eq.'REMARK 290   SMTRY') then
	  read (c(19:70),*) ixyz,nimage1,
     $	    (rotim1(ixyz*3-3+j,nimage1),j=1,3),rotim1(9+ixyz,nimage1)
	end if
	
c	**** read header to get xtal translation vectors ****

c12345678901234567890123456789012345678901234567890123456789012345678901234567890
cCRYST1  169.100  169.100  169.100  90.00  90.00  90.00 
cSCALE1      0.003560  0.002050  0.000000        0.00000 
	
	if (c(1:5).eq.'CRYST') then
	  read (c(7:55),*) alen,blen,clen,alpha,beta,gamma
	end if

C	if (c(1:5).eq.'SCALE') then
C	  read (c(6:70),*) ixyz,
C     $		(boxvecs1(ixyz*3-3+j),j=1,3),boxvecs1(9+ixyz)
C	  if (boxvecs1(9+ixyz).ne.0.d0)
C     $		stop 'dont know how to proc box vec shifts'
C	  write (6,'(12f10.6)') (boxvecs1(i),i=1,9)
C	end if

c	**** read molecule and residue names ****

	if (c(1:6).eq.'SEQRES') then

c12345678901234567890123456789012345678901234567890123456789012345678901234567890
CSEQRES   1 A  365  ALA LEU PHE GLY SER ASN ASP VAL THR THR ALA HIS SER 

	  if (c(11:17).ne.lastres) then
	    lastres= c(11:17)
	    nmol= nmol + 1
	    if (nmol.gt.mmol) stop 'mmol too small'
	    ires= 0
	    read (c(13:17),*) nres(nmol)
	    write (6,*) 'SEQRES:  start molecule',nmol,nres(nmol)
	  end if
	  i2= min (nres(nmol),ires+13)
	  read (c,'(18x,13(1x,a3))') (resname(i,nmol)(1:3),i=ires+1,i2)
	  do i= ires+1,i2
	    write (resname(i,nmol)(4:12),'(i5,a,a)') i,' - ',c(12:12)
C	    write (6,*) i,' "',resname(i,nmol),'"'
	  end do
	  ires= i2
	end if


	if (c(1:4).eq.'ATOM'.or.c(1:4).eq.'HETA') then

c12345678901234567890123456789012345678901234567890123456789012345678901234567890
cATOM      1  N  ATHR A   1      16.885  14.078   3.427  0.50  4.48           N 
cATOM     49 HG21ATHR A   2      11.570  11.067   3.373  0.50  8.54           H  

	  s= c(17:17)
	  c(17:17)= ' '
	  if (c(22:22).eq.' ') c(22:22)= '-'

	  if (c(21:26).ne.lastres) then
	    write (nresname,'(a3,a5,3h - ,a1)') c(18:20),c(23:26)
     $          ,c(22:22)

c	    **** update previous residue ****

	    if (str.ne.' ') then
	        do k= 1,nstr
		  if (str.eq.strlist(k)) istr(ires,imol)= k
	        end do
	        if (istr(ires,imol).eq.0) then
		  nstr= nstr + 1
	 	  if (nstr.gt.mstr) stop 'mstr too small'
		  istr(ires,imol)= nstr
		  strlist(nstr)= str
		  write (6,'(2i4,1x,3a,i3)') ires,imol,resname(ires,imol),
     $		     ' STR varients: ',str,nstr
	        end if
		str= ' '
	    end if

c	    **** identity of new residue and mol ****

	    if (c(1:4).eq.'HETA') then
		nmol= nmol + 1
		if (nmol.gt.mmol) stop 'mmol too small for cofactors'
		imol= nmol
		ires= 1
		nres(imol)= 1
		resname(1,imol)= nresname
C		write (6,*) 'COFACTOR: ',imol,' ',resname(ires,imol)
 	    else if (imol.eq.-1) then
	      imol= 1
	      ires= 1
	    else if (c(21:22).ne.lastres(1:2)) then
	      imol= imol + 1
	      if (imol.gt.nmol) stop 'cant find atom residue name in list'
	      ires= 1
	    end if

	    do while (ires.le.nres(imol) .and. nresname.ne.resname(ires,imol))
	      ires= ires + 1
	    end do

	    if (ires.gt.nres(imol)) then
	      write (6,*) c
	      write (6,*) imol,nres(imol),' "',nresname,'"'
	      stop 'residue not found in SEQRES list'
	    end if

	    lastres= c(21:26)
	  end if

	  if (s.ne.' ') then
	    l= len_trim (str)
	    k= 0
	    do i= 1,l
	      if (str(i:i).eq.s) k= i
	    end do
	    if (k.eq.0) str(l+1:l+1)= s
	  end if

	end if

	goto 100
120	continue

c	**** final residue ****

	if (str.ne.' ') then
	  do k= 1,nstr
	    if (str.eq.strlist(k)) istr(ires,imol)= k
	  end do
	  if (istr(ires,imol).eq.0) then
	    nstr= nstr + 1
	    if (nstr.gt.mstr) stop 'mstr too small'
	    istr(ires,imol)= nstr
	    strlist(nstr)= str
	  end if
	  str= ' '
	end if

c	**** box vecs ****

	  write (6,'(a,3f8.3)') 'cell lengths=',alen,blen,clen
	  write (6,'(a,3f8.3)') 'cell angles =',alpha,beta,gamma
	  pre= acos (0.d0) * 2.d0 / 180.d0
	  boxvecs= 0.d0
          ca= cos (alpha*pre)
          sa= sin (alpha*pre)
          cb= cos (beta*pre)
          sb= sin (beta*pre)
          cg= cos (gamma*pre)
          sg= sin (gamma*pre)
          boxvecs= 0.d0
          boxvecs(1)= alen
          boxvecs(4)= blen * cg
          boxvecs(5)= blen * sg
          boxvecs(7)= clen * cb
          boxvecs(8)= clen * (ca-cb*cg)/sg
          boxvecs(9)= clen * sqrt (1.d0-ca**2-cb**2-cg**2+2.d0*ca*cb*cg)
     $				 / sg
	  boxvecsinv= 0.d0
	  boxvecsinv(1)= 1.d0 / boxvecs(1)
	  boxvecsinv(5)= 1.d0 / boxvecs(5)
	  boxvecsinv(9)= 1.d0 / boxvecs(9)
	  boxvecsinv(4)= - boxvecs(4) / boxvecs(1) / boxvecs(5)
	  boxvecsinv(8)= - boxvecs(8) / boxvecs(9) / boxvecs(5)
	  boxvecsinv(7)= (boxvecs(4)*boxvecs(8)-boxvecs(7)*boxvecs(5))
     $				/ boxvecs(1) / boxvecs(5) / boxvecs(9)
	  write (6,'(a,12f10.3)') 'Unit cell vecs=',(boxvecs(i),i=1,9)
	  write (6,'(a,12f10.6)') 'Unit cell inv =',(boxvecsinv(i),i=1,9)

      write (6,*) 'istr(1,3)=',istr(1,3)

c	************** ask which patterns to include **************

	if (nstr.gt.0) then

	  write (6,*) 'Select from alternate structure possibilities:'
	  do i= 1,nstr
	    if (nargs().eq.2+nstr) then
		call getarg (1+i,c)
		strsel(i)= c
	    else
	      write (6,*) ii,' For options: ',strlist(i),' Enter choice or -'
	      read (5,'(a)') strsel(i)
	    end if
	  end do
	  write (6,*)
	  do i= 1,nstr
	    write (6,*) i,' For options: ',strlist(i),
     $		' Selected ',strsel(i)
	  end do
	end if

c     **** replications of data within range of first molecule ****

      if (nimage1.eq.0) then
	nimage1= 1
	rotim1= 0.d0
	rotim1(1,1)= 1.d0
	rotim1(5,1)= 1.d0
	rotim1(9,1)= 1.d0
      end if
      write (6,*) 'Nber of molecules per unit cell=',nimage1

C      call matinv (boxvecs1,3,wk,0,det,3,wk,iwk)
C      alen1= sqrt (boxvecs1(1)**2 + boxvecs1(4)**2 + boxvecs1(7)**2)
C      blen1= sqrt (boxvecs1(2)**2 + boxvecs1(5)**2 + boxvecs1(8)**2)
C      clen1= sqrt (boxvecs1(3)**2 + boxvecs1(6)**2 + boxvecs1(9)**2)
C      write (6,'(a,3f9.3)') 'transf boxvecs in real space, length=',
C     $		alen1,blen1,clen1
C      write (6,'(3f10.4,f15.4)')
C     $	   ((boxvecs1(ixyz*3-3+j),j=1,3),boxvecs1(9+ixyz),ixyz=1,3)

      write (6,*) 'Nber of molecules per unit cell=',nimage1
      rinscr1= min (alen,blen,clen)/2.
      write (6,*) 'Inscribed radius=',rinscr1
      rinscr1= min (rinscr1,15.d0)
      write (6,*) 'Inscribed radius truncated to',rinscr1

c	************** second read to select str option ***************

	rewind 9

	nares= 0
	iares= 0
	jares= -1
	imol= 1
	ires= 1
	natom= 0
	molind(1)= 1
	oldatom= ' '
	imatom= 0
	nmatom= 0

	read (9,'(a)',end=220) cnew

200   continue
	c= cnew
	len= len_trim(c)

	if (c(1:4).ne.'ATOM'.and.c(1:4).ne.'HETA') then
C	  write (8,'(a)') c(1:len)
	  read (9,'(a)',end=220) cnew
	else

	read (9,'(a)',end=220) cnew
	isaniso= cnew(1:6).eq.'ANISOU'
	aniso1= 0
	if (isaniso) then
c12345678901234567890123456789012345678901234567890123456789012345678901234567890
CANISOU    1  N   LYS A   1     1123    684    525    -83     10     56       N 
	  cani= cnew
	  read (cani,'(28x,6i7)') aniso1
	  read (9,'(a)',end=220) cnew
	end if

	oldatom1= c(12:26)
	oldatom1(6:6)= ' '
	j= 1
	do while (j.le.nmatom .and. oldatom(j).ne.oldatom1)
	  j= j + 1
	end do
	if (j.gt.nmatom) then
	  nmatom= j
	  if (nmatom.gt.m) stop 'nber unique atoms m too small'
	  oldatom(j)= oldatom1
	end if
 

C	read (c(57:60),*) wt
C	iwt= nint (wt*8)
C	wt= iwt/8.d0
CC	write (6,*) oldatom1,' match to',j,imatom(j),sngl(wt),iwt
CC	write (8,*) oldatom1,' match to',j,imatom(j),sngl(wt),iwt
C	do i= 1,iwt
C	  ia= ichar('A') + i - 1 + imatom(j)
C	  if (char(ia).eq.'I') stop 'config I'
C	  write (c(56:60),'(f5.3)') wt/iwt
C	  write (8,'(3a)') c(1:16),char(ia),c(18:len)
CC	  write (8,'(3a,f5.3,a)') c(1:16),char(ia),c(18:55),wt,c(61:len)
C	  if (isaniso) write (8,'(4a)') 'ANISOU',c(7:16),char(ia),cani(18:len)
C	end do
	imatom(j)= imatom(j) + iwt

c0        1         2         3         4         5         6         7         8
c12345678901234567890123456789012345678901234567890123456789012345678901234567890
cATOM      1  N  ATHR A   1      16.885  14.078   3.427  0.50  4.48           N 
cATOM     49 HG21ATHR A   2      11.570  11.067   3.373  0.50  8.54           H  

	  s= c(17:17)
	  c(17:17)= ' '
	  if (c(22:22).eq.' ') c(22:22)= '-'

	  write (nresname,'(a3,a5,3h - ,a1)') c(18:20),c(23:26)
     $          ,c(22:22)
	  do while (imol.le.nmol .and. nresname.ne.resname(ires,imol))
	    do while (ires.le.nres(imol) .and. nresname.ne.resname(ires,imol))
	      ires= ires + 1
	    end do
	    if (ires.gt.nres(imol)) then
	      molend(imol)= natom
	      imol= imol + 1
	      ires= 1
	      molind(imol)= natom+1
	    end if
	  end do
	  if (imol.gt.mmol) stop 'second find of resname failed'

c	  **** see if str variant is in the selected list ****

	  isok= s.eq.' '
	  if (.not.isok .and. istr(ires,imol).gt.0) then
	    do j= 1,len_trim(strsel(istr(ires,imol)))
	      if (s.eq.strsel(istr(ires,imol))(j:j)) isok= .true.
	    end do
C	    write (6,'(l2,2i4,3a,i3,3a)') isok,ires,imol,' "',s,'"',
C     $		istr(ires,imol),' "',strsel(istr(ires,imol)),'"'
	  end if
	  if (.not.isok) goto 200

	  natom= natom + 1
	  if (nares(ires,imol).eq.0) iares(ires,imol)= natom

	  do j= 1,6
	    aniso(j,natom)= aniso1(j)
	  end do

c	  **** atomic number ****

	  if (c(13:13).ne.' ') then
	    aan= c(13:16)
	    c(14:17)= aan
	  end if

	  nat(natom)= 0
	  if (c(77:77).ne.' ') then
	    atss= c(77:78)
	    atss(2:2)= char(ichar(atss(2:2))-'A'+'a')
	  else if (c(78:78).ne.' ') then
	    atss= c(78:78)//' '
	  else
	    atss= c(14:14)//' '
	    if (c(14:15).eq.'CU') atss= 'Cu'
	  end if
	  do i= -1,103
	    if (atsym(i).eq.atss) nat(natom)= i
	  end do
	  if (nat(natom).eq.0) then
	    write (6,*) c
	    write (6,*) '"',atss,'"'
	    stop 'unknown atom'
	  end if

c	  **** atomic name ****

	  ll= len_trim(c(1:17))
	  if (c(14:14).eq.'H' .and. c(ll:ll).ge.'0' .and. c(ll:ll).le.'9') then
	    atname(natom)= c(ll:ll)//c(14:ll-1)
	  else
	    atname(natom)= c(14:17)
	  end if
	  if (atname(natom)(2:2).eq.' ' .and. atname(natom)(3:3).ne.' ')
     $		atname(natom)(2:2)= '_'
	  write (6,*) '"',c(13:17),'" "',atname(natom)

	  c(15:15)= ' '

	  jares(ires,imol)= natom
	  nares(ires,imol)= nares(ires,imol) + 1
	  ats(natom)= atsym(nat(natom))
	  attype(natom)= atsym(nat(natom))
	  flag(natom)= '     -'
	  q(natom)= 0.0
	  atnumb(natom)= natom - molind(imol) + 1
	  ares(natom)= ires
	  amol(natom)= imol
	  read (c(30:54),*) x(natom),y(natom),z(natom)
	  read (c(61:66),*) bfact(natom)
    
	end if

	goto 200
220     continue
    
	molend(imol)= natom

c	******** write out and reread **********

	header(1)= 'forcefield amber94'
	header(2)= 'sys 0 0 1'
	header(3)= 'view 40 0.060047 81.4 41.4 1 0 0 0 1 0 0 0 1'//
     $		   ' -11.658 -10.449 -86.03'
	header(4)= 'seed -1111'
        write (header(5),'(a,f8.3)') '; rinscr=',rinscr1
        write (header(6),'(a,6f9.3)') '; box=',alen,blen,clen,
     $		alpha,beta,gamma
	nhead= 6

	ncon= 0
	fname(nf:nf+3)= '.hin'
	nimage= 1
	boxvecs1= boxvecs
	call hcopyim (1,natom)
C	call writehin (fname,0)

C	natom= 0
C	nmol= 0
C	if (.not.readhin (fname,0)) stop 'reread failed'

c	****  periodic images ****

	rinscr= rinscr1
        rinscr2= min(10.d0,rinscr1)**2
	nimage= nimage1
	do i= 1,nimage1
	  do k= 1,12
	    rotim(k,i)= rotim1(k,i)
	  end do
	  write (6,'(i3,12f8.3)') i,(rotim(k,i),k=1,12)
	end do
	write (6,*) 'nimage=',nimage,nimage1
	call hcopyim (1,natom)
	write (6,*) 'copy done'

c	**** find connections ****

	do imol= 1,nmol

C	write (6,*) 'bonds in mol',imol,molind(imol),molend(imol)
C	do ires= 1,nres(imol)
C	  write (6,'(2i4,1x,a,3i6)') ires,imol,resname(ires,imol),
C     $		iares(ires,imol),jares(ires,imol),nares(ires,imol)
C	end do

	do i= molind(imol),molend(imol)
	  do j= molind(imol),i-1
	      r2= hbl2(i,j)
	      isbond= .false.
	      if (r2.lt.9.0) then
		r= sqrt (r2)
		if (nat(i).eq.1 .or. nat(j).eq.1) then
		  isbond= r.lt.1.4
		else if (nat(i).ge.12 .or. nat(j).ge.12) then
		  isbond= r.lt.2.4
		else if (nat(i).eq.8 .or. nat(j).eq.8) then
		  isbond= r.lt.1.9
		else if (nat(i).eq.7 .or. nat(j).eq.7) then
		  isbond= r.lt.1.9
		else
		  isbond= r.lt.1.9
		end if
	      end if
	      if (isbond) call hbadd (i,j,'s')

	  end do
	end do
	end do

c	**** make sure H connected only once ****

	do i= 1,natom
	  if (nat(i).eq.1 .and. ncon(i).gt.1) then
	    write (6,*) 'multiple bonded H: ',resname(ares(i),amol(i))
	    r2min= 1.d30
	    do jj= 1,ncon(i)
	      j= icon(jj,i) + molind(i)-1
	      r2= hbl2 (i,j)
	      write (6,*) atname(i),atname(j),sqrt(r2)
	      if (r2.lt.r2min) then
		r2min= r2
		jmin= jj
	      end if
	    end do
	    do jj= 1,ncon(i)
	      if (jj.ne.jmin) call hbdel (i,icon(jj,i)+molind(i)-1)
	    end do
	  end if
	end do

c       **** H atom names ****

        do i= 1,natom
          if (nat(i).ne.1) then
	  imol= amol(i)
            inst= 1
            do while (inst.lt.6 .and. atname(i)(inst:inst).eq.' ')
              inst= inst + 1
            end do

c           **** nber of Hs attached to it ****
            nh= 0
            do k1= 1,ncon(i)
              k= icon(k1,i)+molind(imol)-1
              if (nat(k).eq.1) nh= nh + 1
            end do

c           **** set label for each H ****
            ih= 0
            do k1= 1,ncon(i)
              k= icon(k1,i)+molind(imol)-1
              if (nat(k).eq.1) then
                ih= ih + 1
                if (nh.eq.1) then
                  atname(k)= atname(i)(inst:5)
                  atname(k)(1:1)= 'H'
                else
                  write (atname(k),'(i1,a,a)')
     $			ih,'H',atname(i)(inst+1:4)
                end if
              end if
            end do

          end if
        end do
       
c     **** check for unattached atoms ****

      do i= 1,natom
	imol= amol(i)
	ires= ares(i)
	if (ncon(i).eq.0 .and. molind(imol).ne.molend(imol)) then
	  write (6,*) 'ERROR: unattached atom',i,
     $		resname(ires,imol),atname(i)
	  j= iares(ires,imol)
	  write (6,*) j,sqrt(hbl2(i,j))
	end if
      end do

c     **** replications of data within range of first molecule ****

      write (6,*) 'replications'
C      goto 2222
      nimage= 1
C      do ia= 0,1
C	do ib= 0,1
C	  do ic= -1,0
C      do ia= 0,0
C	do ib= 0,0
C	  do ic= 0,0
      do ia= -1,1
	do ib= -1,1
	  do ic= -1,1
	    dx= boxvecs1(1)*ia + boxvecs1(4)*ib + boxvecs1(7)*ic
	    dy= boxvecs1(2)*ia + boxvecs1(5)*ib + boxvecs1(8)*ic
	    dz= boxvecs1(3)*ia + boxvecs1(6)*ib + boxvecs1(9)*ic
	    do im= 1,nimage1
	      if (im.ne.1 .or. ia.ne.0 .or. ib.ne.0 .or. ic.ne.0) then
		do i= 1,natom
		  do j= 1,natom
		    r2= (x(i)-xi(j,im)-dx)**2 + (y(i)-yi(j,im)-dy)**2 +
     $			(z(i)-zi(j,im)-dz)**2
C		    if (r2.lt.rinscr2) then
		      nimage= nimage + 1
		      write (6,'(a,5i4)') 'incl image',nimage,ia,ib,ic,im
		      if (nimage.gt.mimage) stop 'mimage too small'
		      do k= 1,9
		        rotim(k,nimage)= rotim1(k,im)
		      end do
		      rotim(10,nimage)= rotim1(10,im) + dx
		      rotim(11,nimage)= rotim1(11,im) + dy
		      rotim(12,nimage)= rotim1(12,im) + dz
		      goto 350
C		    end if
		  end do
		end do
	      end if
350	      continue
	    end do
	  end do
	end do
      end do
2222  continue

	write (6,*) 'nimage=',nimage,nimage1,natom
      call hcopyim (1,natom)

c     **** check on bond lengths to images ****

      rrmin= 1.d30
      do i= 1,natom
	do j= 1,natom
	 if (nat(i).ne.1 .or. nat(j).ne.1) then
	  do jm= 2,nimage
	    rr= hbl2i (i,j,jm)
	    if (rr.lt.8.0) then
	      ires= ares(i)
	      jres= ares(j)
	      imol= amol(i)
	      jmol= amol(j)
	      write (6,'(a,i4,4a,i3,i4,4a,f8.3)')
     $		 'Short dist: ',imol,atname(i),' ',
     $	         resname(ires,imol),' to image ',jm,jmol,atname(j),' ',
     $		 resname(jres,jmol),' r=',sqrt(rr)
	    end if
	    if (rr.lt.rrmin) then
	      rrmin= rr
	      imin= i
	      jmin= j
	    end if
	  end do
	 end if
	end do
      end do
      write (6,*) natom,nimage,' shortest heavy intermol dist=',
     $		sqrt(rrmin),atname(imin),atname(jmin)

      nhead= 6+nimage
      if (nhead.gt.mh) stop 'mh too small for images'
      do i= 1,nimage
        write (header(6+i),'(a,9f8.4,3f9.3)')
     $		'; image=',(rotim(k,i),k=1,12)
      end do


      fname(nf:nf+3)= '.hin'
      call writehin (fname,0)
      fname(nf:nf+10)= '-images.hin'
      call writehin (fname,1)

      end 
