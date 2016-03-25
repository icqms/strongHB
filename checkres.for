c	**** First routine to run after pdb2hin creates .hin from .pdb ****
c	**** and blank residues are added to get residue numbers correct ****
c	**** nb- delete endmol & mol lines and insert ;samemol in .hin file to join two mols,
c	****		say to span missing chain residue atoms
c
c	**** checks all residues match their templates ***
c       **** adds missing H if all others present ****
c
c	**** UUU = normal complete residue 
c	**** lll = error in residue, including no atoms at all
c	**** UlU = chain starting residue
c	**** UUl = chain end residue
c	**** Ull = chain only, nothing past CB
c	**** lUU = first varient of a proper residue,
c		eg cYS for CYS disulfide (CYX),
c		   hIS for HIS deprotonated on E2 (HID),
c		   aSP for protonated ASP  (ASH)
c		   gLU for protonated GLU
c	**** lUl = second varient of a proper residue,
c		eg hIs for HIS deprotonated on D1 (HIE)
c		   cYs for CYS anion (CYH)
c	**** llU = ????

c    bug ... scon isnt reset to default, causes problems say with reading PHY as must have double bond marked for opt-res to run


      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      include 'template.cmn'

      real*8 xfs4(8),yfs4(8),zfs4(8)
      integer nrok(mt1:mt2),nrchain(mt1:mt2),nrbreak(mt1:mt2),nrmiss(mt1:mt2),
     $     ind(2,ma),indlig(17),itmp(mres,mmol),match(ma)
      logical fixchain,specialres
      character*1 c
      character*2 resname12
      character*3 resname13,resnameuc,bclres(10),chlres(10),bphres(10)
      character*3 uufile
      character*5 newnm,dnamew,hxtrealname(m)
      character*7 namext
      character*40 fname

      data chlres /'CL1','CLa','CLb','CLc','CLd','CLe','CLf','CLg','CLh','PHY'/
      data bclres /'BC1','CLa','CLb','CLc','CLi','CLe','CLf','CLg','CLh','PHY'/
      data bphres /'BP1','CLa','CLb','CLc','CLi','CLe','CLf','CLg','CLh','PHY'/

	open (4,file='checkres.out',status='unknown')

c	********* open .hin file *************

	if (nargs().lt.2) then
	  write (6,*) 'Enter root name of .hin file'
	  fname= '3ENItr-1Hcofphyoddevn'
C	  read (5,'(a)') fname
	else
	  call getarg (1,fname)
	end if

	if (fname(1:1).eq.'-') goto 500
	nf= 1
	do while (fname(nf:nf).ne.' ')
	  nf= nf + 1
	end do
	fname(nf:nf+3)= '.hin'

	if (.not.readhin (fname,0)) goto 500

	call hdeselect1 (1,natom)

c     **** read templates ****

      include 'template.for'

      itala= 0
      itgly= 0
      do it= mt1,mt2
	if (name(it).eq.'ALA') itala= it
	if (name(it).eq.'GLY') itgly= it
      end do
      write (6,*) 'ALA and GLY templates=',itala,itgly

      nnewtemp= 0
      naddh= 0
      naddheavy= 0
      novercon= 0
      noldtermh= 0

      do i= mt1,mt2
        nrok(i)= 0
        nrmiss(i)= 0
	nrchain(i)= 0
	nrbreak(i)= 0
      end do

c     ************************ scan all residues **********************


      do imol= 1,nmol
       ires= 0
       do while (ires.lt.nres(imol))
	ires= ires + 1
	write (6,*) 'STARTING ',ires,imol,resname(ires,imol)

	j0= iares(ires,imol) - 1
	j1= molind(imol) - 1

C      if (resname(ires,imol)(1:3).eq.'PHY') then
C	do i= jares(ires,imol),iares(ires,imol),-1
C	  if (nat(i).eq.1) call hdel (i,1)
C	end do
C      end if

c	**** delete lone pairs ****

	do i= jares(ires,imol),iares(ires,imol),-1
	  if (atname(i)(3:4).eq.'LP') then
	    j= icon(1,i) + j1
	    q(j)= q(j) + q(i)
	    write (6,*) 'deleting lone pair ',resname(ares(i),amol(i)),i,j
	    call hdel (i,1)
	  end if
	end do

c	**** convert names of previouly error residues ****
	c= resname(ires,imol)(1:1)
	if (resname(ires,imol)(2:2).ge.'a' .and. resname(ires,imol)(2:2).ge.'a'
     $		.and. resname(ires,imol)(1:3).ne.'Hoh'
     $		  ) call upcase (resname(ires,imol)(1:3))
	resname13= resname(ires,imol)(1:3)


c	**** nber of heavy atoms ****

	nheavy= 0
	do j= iares(ires,imol),jares(ires,imol)
	  if (nat(j).gt.1) nheavy= nheavy + 1
	end do

c	**** check for alternate residue names ****

	ic= 0
	ih= 0
	ioxt= 0
	do i= iares(ires,imol),jares(ires,imol)
	  if (nat(i).eq.6) ic= ic + 1
	  if (nat(i).eq.1) ih= ih + 1
	  if (atname(i).eq.'  OXT') ioxt= ioxt + 1
	end do

c	**** residue name correction for first residue in a chain ****

	resnameuc= resname(ires,imol)(1:3)
	call upcase (resnameuc)
	resname12= resnameuc(1:2)
c	**** specialres are cofactors split into residues ****
	specialres=
     $		resname12.eq.'CL' .or. resname12.eq.'PH' .or.
     $		resname12.eq.'BC' .or. resname12.eq.'BP' .or.
     $		resname12.eq.'LM' .or.	resname12.eq.'LH' .or. resname12.eq.'PQ'

     	if (nares(ires,imol).gt.1 .and. nres(imol).gt.1 .and.
     $	    	ires.eq.1) then
     	  if (resname(ires,imol)(2:2).ge.'A' .and. resname(ires,imol)(2:2).le.'Z'  .and. .not.specialres) then
             
	    call lowcase (resname(ires,imol)(2:2))
        
c	    **** make sure charge is correct (hypechem gets it wrong) ****
	    qold= 0.d0
	    do i= iares(ires,imol),jares(ires,imol)
	      qold= qold + q(i)
	    end do
	    iqold= nint(qold)
	    q(iares(ires,imol))= q(iares(ires,imol)) + iqold-qold
	    write (6,*) 'first residue charge set to',iqold,' ',resname(ires,imol)

Cc	    **** check to see if this has any H ****
C	    in= iares(ires,imol)
C	    if (nares(ires,imol).gt.nheavy .and. nat(in).eq.7) then
Cc	      **** check to see that first N has 3 H connected ****
C	      if (ncon(in).ne.4) then
C		nmissh= 4 - ncon(in)
C		write (6,*) 'adding missing H to terminal N',nmissh
C		iresend= jares(ires,imol)
C		call hadd (iresend,nmissh,1)
C		call hadda1 (iresend,in,3,1.0d0,109.47d0)
C	      end if
C	    end if

	  end if

	end if
	  
c	**** terminating residue name correction ****

     	if (.not.specialres .and. nres(imol).gt.1 .and. ioxt.gt.0) call lowcase (resname(ires,imol)(3:3))

c	*** cysteine ligand or disulfide correction ****

	if (resname(ires,imol)(1:3).eq.'CYS') then 
          do i= iares(ires,imol),jares(ires,imol)
             if (nat(i).eq.16) then
	        if (ncon(i).eq.1) then
                  do j= 1,natom
	            if (j.ne.i .and. nat(j).gt.8) then
		      r2= hbl2(i,j)
		      if (nat(j).gt.16 .and. r2.lt.2.0**2) then
			write (6,*) 'CYS intermol disupfide link to ',resname(ares(j),amol(j)),atname(j),sqrt(r2)
		 	resname(ires,imol)(1:3)= 'cYS'
		      else if (r2.lt.2.6**2) then
			write (6,*) 'CYS ligated to ',resname(ares(j),amol(j)),atname(j),sqrt(r2)
		 	resname(ires,imol)(1:3)= 'cYs'
		      end if
		    end if
		  end do
		else
                  do j= 1,ncon(i)
                    if (nat(icon(j,i)+j1).eq.16) resname(ires,imol)(1:3)= 'cYS'
	          end do
		end if
             end if
          end do
	  write (6,*) 'cyst residue name= ',resname(ires,imol)
        end if

	if (resname(ires,imol)(1:3).eq.'ASP') then 
	  if (nares(ires,imol).eq.13) resname(ires,imol)(1:3)= 'aSP'
	end if

	if (resname(ires,imol)(1:3).eq.'GLU') then 
	  if (nares(ires,imol).eq.16) resname(ires,imol)(1:3)= 'gLU'
	end if

	if (resname(ires,imol)(1:3).eq.'HOH') then 
	  if (nares(ires,imol).eq.1) resname(ires,imol)(1:3)= 'Hoh'
	  if (nares(ires,imol).eq.1 .and. resname(ires,imol)(1:3).eq.'HOH') then
c	    **** add two hydrogens ****
	    j= iares(ires,imol)
	    call hadd (j,2,1)

	    attype(j)= 'OW'
	    q(j)= -0.834
	    flag(j)= '     h'
	    ncon(j)= 2
	    icon(1,j)= 2
	    icon(2,j)= 3
	    scon(1,j)= 's'
	    scon(2,j)= 's'

	    nat(j+1)= 1
	    atname(j+1)= '   H1'
	    ats(j+1)= 'H '
	    attype(j+1)= 'HW'
	    q(j+1)= 0.417
	    flag(j+1)= '     -'
	    ncon(j+1)= 1
	    icon(1,j+1)= 1
	    scon(1,j+1)= 's'
	    x(j+1)= x(j) + 0.97 * cos (52./pre)
	    y(j+1)= y(j) + 0.97 * sin (52./pre)
	    z(j+1)= z(j)

	    nat(j+2)= 1
	    atname(j+2)= '   H2'
	    ats(j+2)= 'H '
	    attype(j+2)= 'HW'
	    q(j+2)= 0.417
	    flag(j+2)= '     -'
	    ncon(j+2)= 1
	    icon(1,j+2)= 1
	    scon(1,j+2)= 's'
	    x(j+2)= x(j) + 0.97 * cos (52./pre)
	    y(j+2)= y(j) - 0.97 * sin (52./pre)
	    z(j+2)= z(j)

	  end if
	end if

	if (resname(ires,imol)(1:3).eq.'HIS') then 
	  ihashd1= 0
	  ihashe2= 0
	  if (ires.eq.110 .and. imol.eq.1) then
	    write (6,*) 'res 110-A'
	  end if

          do i= jares(ires,imol),iares(ires,imol),-1
            if (atname(i).eq.'  ND1') then
c	      **** look for bond to Mg ****
	      do jmol= 1,nmol
		j= molind(jmol)
		if (nat(j).eq.12) then
		  do jm= 1,nimage
		    if (hbl2i(i,j,jm).lt.2.5**2) ihashd1= -1
		  end do
		end if
	      end do
            else if (atname(i).eq.'  NE2') then
c	      **** look for bond to Mg ****
	      do jmol= 1,nmol
		j= molind(jmol)
		if (nat(j).eq.12) then
		  do jm= 1,nimage
		    if (hbl2i(i,j,jm).lt.2.5**2) ihashe2= -1
		  end do
		end if
	      end do
	    end if
          end do

          do i= jares(ires,imol),iares(ires,imol),-1
            if (atname(i).eq.'  HD1') then
	      if (ihashd1.eq.0) then
		ihashd1= 1
	      else
C		call hdel (i,1)
	      end if
            else if (atname(i).eq.'  HE2') then
	      if (ihashe2.eq.0) then
	        ihashe2= 1
	      else
C		call hdel (i,1)
	      end if
	    end if
          end do
	  if (ihashd1.le.0 .and. ihashe2.eq.1) resname(ires,imol)(1:3)= 'hIs'
	  if (ihashd1.eq.1 .and. ihashe2.le.0) resname(ires,imol)(1:3)= 'hIS'
C	if (ires.eq.110 .or.ires.eq.145 .or. ires.eq.289 .or. ires.eq.297) then
C		resname(ires,imol)(1:3)= 'hIS'
C		call hdel (iares(ires,imol)+17,1)
C	end if
C	if (ires.eq.296) then
C		resname(ires,imol)(1:3)= 'hIs'
C		call hdel (iares(ires,imol)+14,1)
C	end if
	  write (6,*) 'histidine residue name= ',resname(ires,imol)
        end if

c	**** fix up errors in original atom names ****
	if (resname(ires,imol)(1:3).eq.'Lmg' .or.
     $	    resname(ires,imol)(1:3).eq.'Lhg') then

	  do j= iares(ires,imol),jares(ires,imol)
	    write (newnm,'(i5)') j-j0
	    jj= 2
	    if (j-j0.lt.10) jj=3
	    if (ats(j)(2:2).ne.' ') then
	      newnm(jj:jj+1)= ats(j)
	    else
	      newnm(jj+1:jj+1)= ats(j)(1:1)
	    end if
	    atname(j)= newnm
	    write (6,*) j-j0,' newname= ',newnm
	  end do
	end if


Cc	**** reorder H to atom 5 for amino acids ****
C
C	if (nres(imol).gt.1 .and. nares(ires,imol).gt.5 .and.
C     $	  atname(j0+1).eq.'    N' .and. atname(j0+4).eq.'    O') then
C	  k= iares(ires,imol)+5
C	  nswap= 0
C	  do while (k.le.jares(ires,imol) .and. atname(k).ne. '    H')
C	    nswap= nswap + 1
C	    ind(1,nswap)= k - iares(1,imol) + 1
C	    ind(2,nswap)= ind(1,nswap) + 1
C	    k= k + 1
C	  end do
C	  if (k.le.jares(ires,imol)) then
C	    nswap= nswap + 1
C	    ind(1,nswap)= k - iares(1,imol) + 1
C	    ind(2,nswap)= ind(2,1) - 1
C	    call hswap (imol,ind,nswap)
C	  end if
C	end if
 
c	*************** match residue name to template names ***************

	it= mt1
	do while (it.le.nt2 .and. resname(ires,imol)(1:3).ne.name(it))
	  it= it + 1
	end do

	if (it.gt.nt2) then

	  write (4,*) 'residue not found: ',resname(ires,imol),
     $		ires,imol
	  write (6,*) 'residue not found: ',resname(ires,imol),
     $		ires,imol

c	  ******* not found, create new template file *******

	  nt2= nt2 + 1
	  if (nt2.gt.mt2) stop 'NO ROOM TO EXPND TEMPLATES'
	  i= nt2
	  name(i)= resname(ires,imol)(1:3)
	  if (name(i)(1:1).eq.' ') name(i)= resname(ires,imol)(2:3)
	  namext= name(i)
	  l= len_trim(namext)
	  namext(l+1:l+1)= '.'
	  do j= 1,l
	    if (namext(j:j).le.'Z' .and. namext(j:j).ge.'A') then
	      namext(j+l+1:j+l+1)= 'U'
	    else
	      namext(j+l+1:j+l+1)= 'L'
	    end if
	  end do

	  if (nares(ires,imol).gt.0) then

	  write (6,*) 'CREATING TEMPLATE: ',name(i)
	  write (4,*) 'CREATING TEMPLATE: ',name(i)
	  nnewtemp= nnewtemp + 1
          open (2,file='residues\'//namext,status='unknown')
	  write (2,'(a)') resname(ires,imol)
	  if (nares(ires,imol).gt.ma) stop 'too many atoms in new residue'

	  i= nt2
	  do j= 1,nares(ires,imol)

	    dname(j,i)= atname(j0+j)
	    ds(j,i)= ats(j0+j)
	    dtype(j,i)= attype(j0+j)
	    dflag(j,i)= flag(j0+j)
c	    **** make sure that first residue flag isnt copied to template ****
	    if (j.eq.1 .and. dname(j,i).eq.'    N') dflag(j,i)= '     i'
	    dq(j,i)= q(j0+j)
	    dx(j,i)= x(j0+j)
	    dy(j,i)= y(j0+j)
	    dz(j,i)= z(j0+j)
	    dncon(j,i)= ncon(j0+j)
	    k= 1
	    do k1= 1,ncon(j0+j)
	      dicon(k,j,i)= icon(k1,j0+j) -iares(ires,imol) +iares(1,imol)
	      dscon(k,j,i)= scon(k1,j0+j)
	      if (dicon(k,j,i).le.0 .or. dicon(k,j,i).gt.nares(ires,imol))
     $		  then
	        dicon(k,j,i)= 0
	      end if
	      k= k + 1
	    end do

	    write (2,3211) j,dname(j,i),ds(j,i),dtype(j,i),dflag(j,i),
     $	       dq(j,i),dx(j,i),dy(j,i),dz(j,i),dncon(j,i),
     $		(dicon(k,j,i),dscon(k,j,i),k=1,dncon(j,i))
3211	  format (i4,1x,a5,1x,a2,1x,a2,1x,a6,f10.6,3f9.3,i4,20(i4,1x,a1))

	    natdef(j,i)= -1000
	    do k= -1,103
	        if (ds(j,i).eq.atsym(k)) natdef(j,i)= k
	      end do
	    if (natdef(j,i).eq.-1000) then
		write (6,*) ds(j,i)
		stop 'unknown symbol in def'
	    end if
	  end do

	  ndef(i)= j-1
	  close (2)
	  it= nt2

	  end if

	end if

c	***** continue checking ***********

	write (6,*)
	write (6,*) resname(ires,imol),nares(ires,imol),' match to',it,name(it),ndef(it)
	nmissh= 0
	itorig= it

	if (nares(ires,imol).eq.0) then
	  if (.not.specialres) call lowcase (resname(ires,imol)(1:3))
	  nrmiss(it)= nrmiss(it) + 1
	  write (4,*) 'MISSING RESIDUE: ',resname(ires,imol)
	  write (6,*) 'MISSING RESIDUE: ',resname(ires,imol)
	  goto 300
	end if

c	*********** fixups if wrong number of atoms **********

	if (nares(ires,imol).eq.ndef(it)) then
	  nrok(it)= nrok(it) + 1

	else
        
	  if (nheavy.eq.nhdef(it).and.nares(ires,imol).lt.ndef(it)) then
c	    **** missing H atoms only ****
	    naddh= naddh + 1
	    nmissh= ndef(it) - nares(ires,imol)
	    nrok(it)= nrok(it) + 1
	    write (6,*) 'ADDING',nmissh,' MISSING H to ',resname(ires,imol)
	    write (4,*) 'ADDING',nmissh,' MISSING H to ',resname(ires,imol)
             	    
	  else if (fixchain (ires,imol)) then
	    nheavy0= nheavy
	    nheavy= 0
	    do j= iares(ires,imol),jares(ires,imol)
	      if (nat(j).gt.1) nheavy= nheavy + 1
	    end do

c	    **** use template for ALA if truncated chain ****
	    if (resname(ires,imol)(1:3).ne.'ALA' .and. resname(ires,imol)(1:3).ne.'GLY') then
	      nrchain(it)= nrchain(it) + 1
     	      call lowcase (resname(ires,imol)(2:3))
	    else
	      nrok(it)= nrok(it) + 1
	    end if

	    if (nheavy.eq.5) then
	      itorig= it
	      it= itala
	    end if
	    if (nheavy.ne.nheavy0) then
	      naddheavy= naddheavy + 1
	      write (6,*) 'COMPLETED chain of ',resname(ires,imol),' to make ',name(it)
	      write (4,*) 'COMPLETED chain of ',resname(ires,imol),' to make ',name(it)
	    end if
	    call hselres1 (ires,imol,1)

	  else
c	    ****** unfixable atom number mismatch ********
	    nheavy0= nheavy
	    nheavy= 0
	    do j= iares(ires,imol),jares(ires,imol)
	      if (nat(j).gt.1) nheavy= nheavy + 1
	    end do
	    if (nheavy.ne.nheavy0) naddheavy= naddheavy + 1
	    write (4,*) 'CHAIN BROKEN IN ',resname(ires,imol)
	    write (6,*) 'CHAIN BROKEN IN ',resname(ires,imol)
	    if (.not.specialres) call lowcase (resname(ires,imol)(1:3))
	    nrbreak(it)= nrbreak(it) + 1
	    call hselres1 (ires,imol,1)
            
	  end if

	end if

   	nswap= 0

	if (name(it).eq.'FS4') then
c	  *** fix erroneous S connectivity on FS4 ****
	  write (6,*) 'Fixing FS4 connectivity'
	  write (4,*) 'Fixing FS4 connectivity'
	  j4= j0+4

	  do i= 1,4
	    xfs4(i)= x(i+j4)
	    yfs4(i)= y(i+j4)
	    zfs4(i)= z(i+j4)
	    k= 0
	    ncon(i+j4)= 3
	    do j= 1,4
	      if (j.ne.5-i) then
		k= k + 1
		icon(k,i+j4)= j
		icon(k+3,i+j0)= j+4
	      end if
	    end do
	  end do

	  do i= 1,4
	    write (6,765) i,(icon(j,i+j0),j=1,6)
	  end do
	  do i= 1,4
	    write (6,765) i+4,(icon(j,i+j4),j=1,3)
765	    format (i4,4x,6i4)
	    j= 1
	    do while (j.le.4 .and. hbl2(j4+i,j0+j) .lt. 2.5**2)
	      j= j + 1
	    end do
	    if (j.gt.4) stop 'FS4 has strange structure'
	    x(i+j4)= xfs4(5-j)
	    y(i+j4)= yfs4(5-j)
	    z(i+j4)= zfs4(5-j)
	  end do
	end if

	if (resname(ires,imol)(1:3).eq.'HYD') noldtermh= noldtermh + 1

c	********** match default atoms to actual atoms ***************

	nmiss= 0

	do j= 1,ndef(it)
4444	  continue
C	  write (6,*) 'searching default atom',dname(j,it)

	  ifound= 0
	  k= 0
	  do while (ifound.eq.0 .and. k.lt.nares(ires,imol))
	    k= k + 1

C	    write (6,*) k,atname(k+j0),dname(j,it)
	    if (atname(k+j0).eq.dname(j,it) .or. it.ne.itorig .and. dname(j,it).eq.'  3HB'.and.atname(k+j0).eq.'  HXT') then
		write (6,*) k,' match',atname(k+j0)

c	      *** check for symbol, attype, flag, charge ****
	      if (ats(k+j0).ne.ds(j,it)) then
	        ats(k+j0)= ds(j,it)
	      end if

	      if (attype(k+j0).ne.dtype(j,it)) then
C	        write (6,835) resname(ires,imol),
C     $		      k,dname(j,it),attype(k+j0),dtype(j,it)
C	        write (4,835) resname(ires,imol),
C     $		      k,dname(j,it),attype(k+j0),dtype(j,it)
	        attype(k+j0)= dtype(j,it)
835	        format (1x,a,' FIXED atom type ',i4,3(1x,a))
	      end if

	      if (flag(k+j0).ne.dflag(j,it)) then
C		write (6,*) 'changed flags: ',flag(k+j0),dflag(j,it)
	        flag(k+j0)= dflag(j,it)
	      end if

	      if (abs(q(k+j0)-dq(j,it)).gt.1.D-7) then
C	        write (4,984) resname(ires,imol),dname(j,it),
C     $		      q(k+j0),dq(j,it)
C	        write (6,984) resname(ires,imol),dname(j,it),
C     $		      q(k+j0),dq(j,it)
984	        format (' FIXED atom q: ',a,1x,a,2f8.4)
	        if (atname(k+j0).ne.'  HTM' .and. atname(k+j0).ne.'  HXT') q(k+j0)= dq(j,it)
	      end if

c	      **** if nber connections wrong, delete any H that is on there ****

	      if (ncon(k+j0).lt.dncon(j,it)) then
		write (6,*) 'nber con wrong',k+j0,ncon(k+j0),dncon(j,it)
		kkk= 1
		do while (kkk.le.ncon(k+j0))
		  kkkk= icon(kkk,k+j0) + j1
		  if (nat(kkkk).eq.1) then
		    write (6,*) 'deleting existing H',kkkk,' on incomplete heavy atom ',kkk,atname(kkk)
		    call hdel (kkkk,1)
		  end if
		  kkk= kkk + 1
		end do
	      end if

c	      **** form list of atoms to swap ****

	      ifound= k
	      match(j)= k+j0
	      write (6,*) 'template match atom ',j,' to actual atom',k+j0

	    end if
	  end do

c	  **** atom not found ****

	  if (ifound.eq.0) then

	    match(j)= 0

	    if (natdef(j,it).ne.1)  then
	      nmiss= nmiss + 1
	      write (6,*) 'atom name not found ',resname(ires,imol),
     $		     ' ' ,dname(j,it)

	    else
c	      **** add missing hydrogens, kdc is default heavy atom, kc real one ****
	      kdc= dicon(1,j,it)
	      kc= match(kdc)
	      if (kc.gt.0) then
	        write (6,*) 'missing ',j,dname(j,it),' on ',name(it),' heavy atom ',dname(kdc,it),kdc,kc
		do kk= 1,ndef(it)
		  write (6,*) kk,dname(kk,it)
		end do
		do kk= iares(ires,imol),jares(ires,imol)
		  write (6,'(i3,a,i3,6i6)') kk,atname(kk),ncon(kk),(icon(kkk,kk),kkk=1,ncon(kk))
		end do
	        iresend= jares(ires,imol)
		if (atname(kc).eq.'    N' .and. ncon(kc).eq.1) then
	          call hadd (iresend,2,1)
		  call hadda1 (iresend,kc,2,1.01d0,120.d0)
		  call hdel (iresend,1)
		  iresend= iresend - 1
		  atname(iresend)= '    H'
		elseif (atname(kc).eq.'    N' .and. ncon(kc).eq.2) then
	          call hadd (iresend,1,1)
		  call hadda1 (iresend,kc,2,1.01d0,120.d0)
		elseif (nat(kc).eq.6) then
c		  **** check only H missing ****
		  nhv= 0
		  do kk= 1,ncon(kc)
		    write (6,*) kk,icon(kk,kc),j1,nat(icon(kk,kc+j1))
		    if (nat(icon(kk,kc)+j1).gt.1) nhv= nhv + 1
		  end do
		  nhvd= 0
		  do kk= 1,dncon(kdc,it)
		    if (dicon(kk,kdc,it).gt.0) then
		      if (natdef(dicon(kk,kdc,it),it).gt.1) nhvd= nhvd + 1
		    else
		      nhvd= nhvd + 1
		    end if
		  end do
		  write (6,*) 'H on C: nheavy=',nhv,nhvd,ncon(kc),dncon(kdc,it)
		  if (nhv.eq.nhvd) then
		    ihyb= dncon(kdc,it)-1
		    aa= 120.
		    if (ihyb.eq.3) aa= 109.5
	            call hadd (iresend,dncon(kdc,it)-ncon(kc),1)
		    call hadda1 (iresend,kc,ihyb,1.1d0,aa)
		  else
	            call hadd (iresend,1,1)
	            iresend= iresend + 1
		    write (6,*) 'calling hadda2',j,kdc,iresend,kc
	            call hadda2 (0,j,kdc,dx(1,it),dy(1,it),dz(1,it),dncon(1,it),
     $		      dicon(1,1,it),dname(1,it), iresend,kc)
		  end if
		else
	          call hadd (iresend,1,1)
	          iresend= iresend + 1
		  write (6,*) 'calling hadda2',j,kdc,iresend,kc
	          call hadda2 (0,j,kdc,dx(1,it),dy(1,it),dz(1,it),dncon(1,it),
     $		    dicon(1,1,it),dname(1,it), iresend,kc)
		end if
		atname(jares(ires,imol))= dname(j,it)
	        goto 4444
	      end if
	    end if
	  end if

c	  **** change name of fictitious H3 of ALA to HXT ***
	  if (itorig.ne.it .and. atname(j).eq.'  3HB') then
	    write (6,*) 'changing 3HB to HXT'
	    atname(j)= '  HXT'
	    hxtrealname(j)= '  3HB'
	  end if

	  k= ifound
	  if (k.ne.j) then
	    nswap= nswap + 1
	    ind(1,nswap)= k+j0-j1
	    ind(2,nswap)= j+j0-j1
C	    write (6,*) 'swap',nswap,k,j,atname(k+j0),dname(j,it)
	  end if

	end do

c	********** swap atoms to correct order in residue **********

	if (nmiss.eq.0 .and. nswap.ne.0) then
	  write (6,*) 'reorder ',resname(ires,imol),imol,ires
	  call hswap (imol,ind,nswap)
	end if

C	if (ires.eq.23) then
C	  call writehin ('jnk.hin',0)
C	  stop 'test 23'
C	end if

300	continue
	itmp(ires,imol)= it

       end do
      end do

      call writehin ('jnk.hin',0)

c     ******************** second pass to check connectivities ***************************

      do i= 4,6,2
	write (i,*)
	write (i,*) 'SECOND PASS TO ADD TERMINATING HYD RESIDUES AND HXT ATOMS'
	write (i,*)
      end do

      do imol= 1,nmol
       do ires= 1,nres(imol)

	it= itmp(ires,imol)
	j0= iares(ires,imol) - 1
	j1= molind(imol) - 1

	do i= iares(ires,imol),jares(ires,imol)
	 if (resname(ires,imol)(1:8).ne.'BCR 4009' .and. resname(ires,imol)(1:3).ne.'LHG' .and.
     $	     resname(ires,imol)(1:3).ne.'HYD' .and. atname(i).ne.'  HXT') then

	 write (6,'(a,a,a,i3,10i6)') 'checking: ',resname(ires,imol),atname(i),ncon(i),(icon(ii,i),ii=1,ncon(i))

c	  **** find atom in template ****
	  j= 1
	  do while (j.le.ndef(it) .and. atname(i).ne.dname(j,it))
	    j= j + 1
	  end do
	  if (j.gt.ndef(it)) then
	    write (6,*) 'ERROR: atom type not found in template ',resname(ires,imol),atname(i)
	    stop 'unknown atom name'
	  end if

	  if (ncon(i).gt.dncon(j,it)) then
	    write (4,*) 'WARNIMNG: more conenctions on atom than in template ',resname(ires,imol),atname(i)
	    write (6,*) 'WARNIMNG: more conenctions on atom than in template ',resname(ires,imol),atname(i)
	    novercon= novercon + 1

	  else if (ncon(i).lt.dncon(j,it)) then
c	    **** nber connections too small ****

	    do jcc= 1,dncon(j,it)
	      jc= dicon(jcc,j,it)
	      if (jc.gt.0) then
	        dnamew= dname(jc,it)
	      else if (dname(j,it).eq.'    N') then
		dnamew= '    C'
	      else if (dname(j,it).eq.'    C') then
		dnamew= '    N'
	      else if (dname(j,it).eq.'    S') then
		dnamew= '    S'
		stop 'broken S-S link not programmed'
	      else
		stop 'dont know type of external connection'
	      end if
	      write (6,*) 'searching for con to ',jc,dnamew
	      ic= 0
	      do icc= 1,ncon(i)
	        if (atname(icon(icc,i)+j1).eq.dnamew) ic= icc
	      end do

	      if (ic.eq.0) then
		write (6,*) 'connection missing',ncon(i),jc
c		**** missing atom identified ****

		if (jc.eq.0) then
c	 	  **** add terminating residue ****
		  ie= molend(imol)
		  call hadd (ie,1,1)
		  iresinto= ares(ie+1)
		  jares(iresinto,imol)= jares(iresinto,imol) - 1
		  nares(iresinto,imol)= nares(iresinto,imol) - 1
		  nres(imol)= nres(imol) + 1
		  if (nres(imol).gt.mres) stop 'mres too small to add HXT'
		  ll= len_trim (resname(nres(imol)-1,imol))
		  write (resname(nres(imol),imol),'(a,i5,a)') 'HYD',nres(imol),resname(nres(imol)-1,imol)(ll-3:ll)
		  write (6,*) 'ADDING HYD terminating residue ',resname(nres(imol),imol),' to ',resname(ires,imol)
	          call hadda1 (ie,i,2,1.05D0,120.d0)
		  iares(nres(imol),imol)= ie
		  jares(nres(imol),imol)= ie
		  nares(nres(imol),imol)= 1
		  ares(ie)= nres(imol)
		  nat(ie)= 1
		  ats(ie)= 'H '
		  attype(ie)= ' H'
		  atname(ie)= '  HTM'
		  flag(ie)= '     -'
	          hxtrealname(ie)= dnamew
	          call hselres1 (ires,imol,1)
		  call hselres1 (nres(imol),imol,1)

		else
c		  **** add terminating HTM inside residue ****
		  ie= jares(ires,imol)
		  call hadd (ie,1,1)
		  ie= ie + 1
		  nat(ie)= 1
		  ats(ie)= 'H '
		  attype(ie)= ' H'
	          call hadda2 (0,jc,j,dx(1,it),dy(1,it),dz(1,it),dncon(1,it),
     $		    dicon(1,1,it),dname(1,it), ie,i)
		  atname(ie)= '  HXT'
		  flag(ie)= '     -'
		  r= sqrt (hbl2(ie,i))
		  rnew= 1.1
		  if (nat(i).eq.7) rnew= 1.01
		  if (nat(i).eq.8) rnew= 0.97
		  x(ie)= x(i) + (x(ie)-x(i))/r*rnew
		  y(ie)= y(i) + (y(ie)-y(i))/r*rnew
		  z(ie)= z(i) + (z(ie)-z(i))/r*rnew
	          hxtrealname(ie)= dnamew
	          call hselres1 (ires,imol,1)
		  write (6,*) 'ADDING HXT terminating atom  inside ',resname(ires,imol)

		end if

	      end if
	    end do

	  end if

	 end if
	end do

       end do
      end do

c     ******************** third pass to split cofactors ***************************

      call writehin ('jnk1.hin',0)

      do imol= 1,nmol

	ires= 1
	resname13= resname(ires,imol)(1:3)
        j1= molind(imol) - 1

	if ((resname13(1:2).eq.'BC' .or. resname13(1:2).eq.'BP' .or. resname13.eq.'CHL') .and. nres(imol).eq.1) then
	  write (6,*) 'SPLITTING ',resname(ires,imol)
	  nares0= nares(1,imol)

c	  **** insert residue names ****
	  do jres= 1,10
	    if (resname13.eq.'BCL') then
	      resname(jres,imol)= bclres(jres) // resname(1,imol)(4:12)
	    else if (resname13.eq.'BPH') then
	      resname(jres,imol)= bphres(jres) // resname(1,imol)(4:12)
	    else if (resname13(1:2).eq.'CHL') then
	      resname(jres,imol)= bclres(jres) // resname(1,imol)(4:12)
	    end if
	    write (6,*) resname(jres,imol)
	  end do

c	  **** find template atoms and match to real ones  ****

	  nswap= 0
	  do jres= 1,10
	    write (6,*) 'searching inside template residue ',resname(jres,imol)
	    it= mt1
	    do while (it.le.nt2 .and. name(it).ne.resname(jres,imol)(1:3))
	      it= it + 1
	    end do
C	if (jres.eq.1) then
C	dname(ndef(it)+1,it)= '   HB'
C	dname(ndef(it)+2,it)= '   HD'
C	dname(ndef(it)+1,it)= '  H2C'
C	dname(ndef(it)+2,it)= '  H3C'
C	ndef(it)= ndef(it) + 2
C	end if
	    if (it.gt.nt2) stop 'template not found for splitting'
	    do j= 1,ndef(it)
	      match(j)= 0
	      k= iares(1,imol)
	      do while (k.le.jares(1,imol) .and. atname(k).ne.dname(j,it) .and. hxtrealname(k).ne.dname(j,it))
		k= k + 1
	      end do
	      if (k.le.jares(1,imol)) then
	        nswap= nswap + 1
	        ind(1,nswap)= k-j1
	        ind(2,nswap)= nswap
		match(j)= k
		nares(1,imol)= nares(1,imol) - 1
		nares(jres,imol)= nares(jres,imol) + 1
	        write (6,'(4a,i6,a,i4,2i6)') 'match template ',name(it),dname(j,it),
     $		   ' to atom nber',k, ' swapping: ',nswap,ind(1,nswap),ind(2,nswap)
	      end if
	    end do
	  end do

c	  *** make sure all real atoms found ****

	  if (nswap.ne.nares0) then
	    write (6,*) 'swapping',nswap,' but nares=',nares0
	    do k= iares(1,imol),iares(1,imol)+nares0-1
	      jj= 0
	      do j= 1,nswap
		if (ind(1,j).eq.k-j1) jj= j
	      end do
	      if (jj.eq.0) then
		write (6,*) 'failed to match atom ',k,atname(k)
	      end if
	    end do
	    stop 'error in splitting'
	  end if

c	  **** move to correct atom order ****

	  call hswap (imol,ind,nswap)

c	  **** create indices ****

	  nres(imol)= 10
	  jares(1,imol)= iares(1,imol) + nares(1,imol) - 1
	  do jres= 2,10
	    iares(jres,imol)= jares(jres-1,imol) + 1
	    jares(jres,imol)= jares(jres-1,imol) + nares(jres,imol)
	    do j= iares(jres,imol),jares(jres,imol)
	      ares(j)= jres
	    end do
	  end do

	end if
      end do

c     *******************************************************************************

c     **** check H distances ****

      nhlenerr= 0
      do i= 1,natom
	  if (nat(i).gt.1) then
	    imol= amol(i)
	    ires= ares(i)
	    i1= molind(imol) - 1
	    do jj= 1,ncon(i)
	      j= icon(jj,i) + i1
	      if (nat(j).eq.1) then
		r= sqrt (hbl2(i,j))
		if (nat(i).eq.16 .and. r.gt.1.55 .or. nat(i).ne.16 .and. r.gt.1.15 .or. r.lt. 0.9) then
		  write (6,*) 'H length error:',r,
     $		    resname(ires,imol),atname(i)
		  call hselres (ires,imol,1,1)
		  nhlenerr= nhlenerr + 1
		  rnew= 1.1
		  if (nat(i).eq.7) rnew= 1.01
		  if (nat(i).eq.8) rnew= 0.97
		  if (nat(i).eq.16) rnew= 0.97
		  x(j)= x(i) + (x(j)-x(i))/r*rnew
		  y(j)= y(i) + (y(j)-y(i))/r*rnew
		  z(j)= z(i) + (z(j)-z(i))/r*rnew
		end if
	      end if
	    end do
	  end if
      end do

c     *******************************************************************************

      do io= 4,6,2
      write (io,*)
      write (io,*) 'nber of complete, chain OK, chain broken, missing residues:'
      nok= 0
      nmiss= 0
      nchain= 0
      nbreak= 0
      do i= mt1,nt2
	write (io,'(1x,a,4i5)') name(i),nrok(i),nrchain(i),nrbreak(i),nrmiss(i)
	nok= nok + nrok(i)
	nmiss= nmiss + nrmiss(i)
	nchain= nchain + nrchain(i)
	nbreak= nbreak + nrbreak(i)
      end do
      write (io,'(a,4i6)') 'total nber OK=                              ',nok
      write (io,'(a,4i6)') 'total nber chain OK but incompete sidechain=',nchain
      write (io,'(a,4i6)') 'total nber residues with inclomplet chain=  ',nbreak
      write (io,'(a,4i6)') 'total nber missing residues=                ',nmiss
      write (io,'(a,4i6)') 'total nber initail HYD term residues=       ',noldtermh
      write (io,'(a,3i6)') 'nber new templates=                         ',nnewtemp
      write (io,'(a,3i6)') 'nber residues with added hydrogens=         ',naddh
      write (io,'(a,3i6)') 'nber residues with added heavy atoms=       ',naddheavy
      write (io,'(a,3i6)') 'nber atoms with more cons than in template= ',novercon
      write (io,'(a,3i6)') 'nber H length errors=                       ',nhlenerr
      end do

c     ****** write output ******

      fname(nf:nf+4)= 'c.hin'
      call writehin (fname,0)

850   format (10(1x,a3,i4))
851   format (/' nber sel atoms=',i6,' total nber of residues in ',a/)
852   format (1x,a,50i4)
853   format (/' total nber of residues in all mutants:'/ (i4,1x,a))
854   format (4x,50i4)

500   continue
      end

c     *********************************************************************

      logical function fixchain (ires,imol)

c     **** completes chain atoms if possible, returnd true if all heavy atoms are N CA C O and CB(if not GLY) ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      real*8 xsav(4),ysav(4),zsav(4)
      integer ka(5,-1:1)
      logical gly,fixchain1,not13
      data kn,kca,kc,ko,kcb/1,2,3,4,5/

      j0= iares(ires,imol) - 1
      j1= molind(imol) - 1
      gly= resname(ires,imol)(1:3) .eq. 'GLY'
      iresend= jares(ires,imol)
      iresend0= iresend

      write (6,*) 'attempting to fix chain ',resname(ires,imol)

c     **** find N CA C O CB atoms in this, previous, and next residues ****

      ka= 0
      do jress= -1,1
	jres= jress + ires
	  write (6,*) 'scanning residue ',jres,' ',resname(jres,imol)
	  do j= iares(jres,imol),jares(jres,imol)
	    write (6,*) j,' "',atname(j),'"'
	    if (atname(j).eq.'    N') ka(kn ,jress)= j
	    if (atname(j).eq.'   CA') ka(kca,jress)= j
	    if (atname(j).eq.'    C') ka(kc ,jress)= j
	    if (atname(j).eq.'    O') ka(ko ,jress)= j
	    if (atname(j).eq.'   CB') ka(kcb,jress)= j
C	    if (jress.eq.0 .and. nat(j).eq.1) stop 'incomplete res has a H ??'
	  end do
	  write (6,'(i3,a,5i6)') jres,' ka=',(ka(i,jress),i=1,5)
      end do

      if (ka(kn,0).eq.0) then
c       **** missing N ****
	if (ka(kca,-1).ne.0 .and. ka(kc,-1).ne.0 .and. ka(ko,-1).ne.0) then
	  call hadd (iresend,1,1)
	  call hadda1 (iresend,ka(kc,-1),2,1.37d0,120.d0)
	  ka(kn,0)= iresend
	else if (ka(kca,0).ne.0 .and. ka(kc,0).ne.0 .and. ka(kcb,0).ne.0) then
	  call hadd (iresend,1,1)
	  iresend= iresend + 1
	  ka(kn,0)= iresend
	  call haddz (ka(kn,0),ka(kca,0),ka(kc,0),ka(kcb,0),1.43d0,109.5d0,-120.d0)
	end if
	if (ka(kn,0).ne.0) then
	  atname(iresend)= '    N'
	  nat(iresend)= 7
	  write (6,*) 'ADDED N'
	end if
      end if

      if (ka(kc,0).eq.0) then
c       **** missing C ****
	if (ka(kca,0).ne.0 .and. ka(kn,0).ne.0 .and. ka(kcb,0).ne.0) then
	  call hadd (iresend,1,1)
	  iresend= iresend + 1
	  ka(kc,0)= iresend
	  call haddz (ka(kc,0),ka(kca,0),ka(kn,0),ka(kcb,0),1.49d0,109.5d0,120.d0)
	  atname(iresend)= '    C'
	  nat(iresend)= 6
	  write (6,*) 'ADDED C'
	end if
      end if

      if (ka(kca,0).eq.0) then
c       **** missing CA ****
	if (ka(kc,0).ne.0 .and. ka(kn,1).ne.0 .and. ka(ko,0).ne.0) then
	  call hadd (iresend,1,1)
	  iresend= iresend + 1
	  ka(kca,0)= iresend
	  call haddz (ka(kca,0),ka(kc,0),ka(kn,1)+1,ka(ko,0),1.49d0,120.d0,180.d0)
	  atname(iresend)= '   CA'
	  nat(iresend)= 6
	  write (6,*) 'ADDED CA'
	end if
      end if

      if (ka(kc,0).eq.0 .and. ka(kca,0).eq.0) then
c       **** missing C and CA ****
	if (ka(kn,0).ne.0 .and. ka(kc,-1).ne.0 .and. ka(kn,1).ne.0 .and. ka(kca,1).ne.0) then
	  call hcopyim (1,natom)
	  nadd= 4
	  if (gly) nadd= 3
	  write (6,*) 'missing C and CA, adding ',nadd
	  call hadd (iresend,nadd,1)
	  ka(kca,0)= iresend+1
	  ka(kc,0)= iresend+2
	  ka(ko,0)= iresend+3
	  if (.not.gly) ka(kcb,0)= iresend+4
	  do ii= 1,5
	    ka(ii,1)= ka(ii,1) + nadd
	  end do

c	  **** grid search on the 2 torsional angles ****
	  errmin= 1.d30
	  do torc= 0.,350.,10.
	    call haddz (ka(kc,0),ka(kn,1),ka(kca,1),1,1.37d0,120.d0,torc)
	    do torca= 0.,350.,10.
	      call haddz (ka(kca,0),ka(kn,0),ka(kc,-1),1,1.43d0,120.d0,torca)
	      call hbadd (ka(kc,0),ka(kca,0),'s')
	      call hbdel (ka(ko,0),ka(kc,0))
	      kk= ka(ko,0)-1
	      call hadda1 (kk,ka(kc,0),2,1.33d0,120.d0)
	      if (.not.gly) call haddz (ka(kcb,0),ka(kca,0),ka(kn,0),ka(kc,0),1.53d0,109.5d0,-120.d0)

c	      **** check C-CA distance ****
	      rloc= sqrt (hbl2 (ka(kc,0),ka(kca,0)))
	      err= (rloc - 1.48)**2

c	      **** check bond angles ****
	      rcca0= sqrt (hbl2 (ka(kc,-1),ka(kca,0)))
	      err= err + (rcca0-2.46)**2 / 8.
	      rca0n= sqrt (hbl2 (ka(kca,0),ka(kn,1)))
	      err= err + (rca0n-2.43)**2 / 8.
	      rc0ca= sqrt (hbl2 (ka(kc,0),ka(kca,1)))
	      err= err + (rcca0-2.46)**2 / 8.
	      rn0c0= sqrt (hbl2 (ka(kn,0),ka(kc,0)))
	      err= err + (rn0c0-2.43)**2 / 8.

	      if (err.lt.0.2) then
c	        **** check that C CA O and CB dont hit anything ****
	        call hcopyim (iresend+1,nadd)
	        extmin= 1.d30
	        do i= iresend+1,iresend+nadd
		  do j= 1,natom
		    if (nat(j).gt.1) then
c		      **** exclude 1-2 and 1-3 ints ****
		      not13= amol(j).ne.imol
		      if (.not.not13) then
			not13= .true.
		        do jj= 1,ncon(j)
			  if (icon(jj,j)+j1.eq.i) not13= .false.
			  do ii= 1,ncon(i)
			    if (icon(ii,i).eq.icon(jj,j)) not13= .false.
			  end do
			end do
		      end if
		      if (not13) then
		        do jm= 1,nimage
		          r2= hbl2i (i,j,jm)
			  if (extmin.gt.r2) then
		            extmin= r2
			    write (6,'(2hE=,2f8.1,i6,a,2i6,f12.6,a,a)')
     $				torc,torca,i,atname(i),j,jm,sqrt(r2),resname(ares(j),amol(j)),atname(j)
			  end if
		        end do
		      end if
		    end if
		  end do
	        end do

		extmin= sqrt (extmin)
		if (extmin.lt.3.0) err= err + (extmin-3.0)**2 / 10.
		err= sqrt (err)
	        if (err.lt.errmin) then
		  errmin= err
		  write (6,'(a,7f8.3)') 'error=',err,rloc,rcca0,rn0c0,rca0n,rc0ca,extmin
		  do i= 1,nadd
		    xsav(i)= x(iresend+i)
		    ysav(i)= y(iresend+i)
		    zsav(i)= z(iresend+i)
		  end do
	        end if
	      end if

	    end do
	  end do

C	i= iares(ires+1,imol)+1
C	call hdel (i,natom-i+1)
C	nmol= 1
C	nres(1)= ires+1

	  if (errmin.lt.0.4) then
c	    **** added OK  ****
	    write (6,*) 'ADDED C and CA',nadd,' len err=',errmin
	    do i= 1,nadd
	      x(iresend+i)= xsav(i)
	      y(iresend+i)= ysav(i)
	      z(iresend+i)= zsav(i)
	    end do
	    atname(iresend+1)= '   CA'
	    atname(iresend+2)= '    C'
	    atname(iresend+3)= '    O'
	    nat(iresend+1)= 6
	    nat(iresend+2)= 6
	    nat(iresend+3)= 8
	    if (.not.gly) then
	      atname(iresend+4)= '   CB'
	      nat(iresend+4)= 6
	    end if
	    iresend= iresend + nadd
	  else
	    write (6,*) 'ADDING C and CA failed',nadd,' err=',errmin
	    call hdel (iresend+1,nadd)
	    ka(kca,0)= 0
	    ka(kc,0)= 0
	  end if

	end if
      end if

c     **** check to see if chain completed, add O and CB if not done ****

      fixchain1= ka(kc,0).gt.0 .and. ka(kn,0).gt.0 .and. ka(kca,0).gt.1
      fixchain= fixchain1
      if (fixchain) then
	if (ka(ko,0).eq.0) then
	  kk= ka(ko,0)-1
	  call hadda1 (kk,ka(kc,0),2,1.33d0,120.d0)
	end if
	if (.not.gly .and. ka(kcb,0).eq.0) then
	  call haddz (ka(kcb,0),ka(kca,0),ka(kn,0),ka(kc,0),1.53d0,109.5d0,-120.d0)
	end if
      end if

      write (4,*) 'ADDED ',iresend-iresend0,' HEAVY ATOMS to ',resname(ires,imol),fixchain1
      call hcopyim (1,natom)
	do ii= iares(ires,imol),jares(ires,imol)
	  write (6,'(i3,a,i3,12i6)') ii,atname(ii),ncon(ii),(icon(iii,ii),iii=1,ncon(ii))
	end do

      return
      end




