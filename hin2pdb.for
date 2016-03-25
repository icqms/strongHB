c       **** writes .hin file out in .pdb format ****

	implicit real*8 (a-h,o-z)
	include 'readhin.cmn'
	parameter (morig=10000, mcopy=8)

	real*8 rotimo(12,mcopy)
	integer resnumbi,resnumbj(morig),resnumb,newnumb(m)
	integer iresorig(mcopy,morig),imolorig(mcopy,morig)
	character*128 fname,line
	character*4 moln,atname1
	character*3 resn
        character*1 chain,str
        character*6 typ,typ0

	common nchainorig

	open (4,file='hin2pdb.out',status='unknown')

c	**** input file name ****

	if (nargs().lt.2) then
	  write (6,*) 'Enter root name of .hin file'
	  read (5,'(a)') fname
	else
	  call getarg (1,fname)
	end if

	if (nargs().eq.3) then
c	  **** July 2011 JRR created firl 2VB1-8spec.hin as data for this ****
	  call getarg (2,line)
	  write (6,*) 'reading image file: ',line
	  open (1,file=line,status='old')
C	  read (1,'(a)') line
	  nrep= 0
20	  read (1,'(a)',end=30) line
	    write (6,*) line(1:70)
	    nrep= nrep + 1
	    if (nrep.gt.mcopy) stop 'too many images specified'
	    read (line(9:128),*) (rotimo(j,nrep),j=1,12)
	    goto 20
30	  close (1)
	  write (6,*) nrep,' images specified'
	else
	  nrep= 1
	  rotimo= 0.d0
	  rotimo(1,1)= 1.d0
	  rotimo(5,1)= 1.d0
	  rotimo(9,1)= 1.d0
	  str= ' '
	end if
	write (6,*) 'Nber chains in orig str=',nchainorig

	nf= 1
	do while (fname(nf:nf).ne.' ')
	  nf= nf + 1
	end do
	fname(nf:nf+3)= '.hin'

	if (.not.readhin (fname,0)) stop 'file not found'

	fname(nf:nf+3)= '.pdb'
	open (2,file=fname,status='unknown')

c	**** find largest chain ****
	nchainorig= 1000
	if (nrep.gt.1) then
	  ic= 0
	  do imol= 1,nmol
	    do ires= 1,nres(imol)
	      k= len_trim (resname(ires,imol))
	      ic= max (ic, ichar(resname(ires,imol)(k:k)))
	    end do
	  end do
	  write (6,*) 'ic=',ic,ichar('A')
	  nchainorig= (ic-ichar('A')) / nrep + 1
	  write (6,*) 'Nber chains in original system=',nchainorig
	end if

c	**** delete lone pairs ****
	write (6,*) 'nmol,natom=',nmol,natom
	do i= natom,1,-1
	  if (nat(i).le.0) call hdel (i,1)
	end do
	write (6,*) 'Nber atoms after lone pair deletion=',natom

	iout= 0
	orig= 0
	norig= 0
	imolorig= 0

c	**** sort out residues into originals + copies ****

	do imol= 1,nmol
	  do ires= 1,nres(imol)
	    resnumbi= resnumb (resname(ires,imol),imnumbi)
	    do jorig= 1,norig
	      if (resnumbi.eq.resnumbj(jorig)) then
	        iresorig(imnumbi,jorig)= ires
	        imolorig(imnumbi,jorig)= imol
		goto 130
	      end if
	    end do
	    norig= norig + 1
	    if (norig.gt.morig) stop 'morig too small'
	    resnumbj(norig)= resnumbi
	    iresorig(imnumbi,norig)= ires
	    imolorig(imnumbi,norig)= imol
	    write (6,*) resname(ires,imol),resnumbi,imnumbi,norig
130	    continue
	  end do
	end do

	write (6,*) 'nber original residues=',norig

	weight= 1.d0 / nrep
	typ= 'ATOM  '

	do jorig= 1,norig
	 do jcopy= 1,nrep

	  imol= imolorig(jcopy,jorig)
	  ires= iresorig(jcopy,jorig)
	  if (imol.eq.0) goto 3333

	  typ0= typ
	  typ= 'ATOM  '  
	  if (nres(imol).le.10) typ= 'HETATM'
	  if (typ.eq.'HETATM' .and. typ0.eq.'ATOM  ') then
		iout= iout + 1
     		write (2,800) 'TER   ',iout,'    ',' ',resn,chain,moln
	  end if

c	  **** parse residue name ****

	  i= 12
	  do while (resname(ires,imol)(i:i).eq.' ')
	    i= i - 1
	  end do
	  chain= resname(ires,imol)(i:i)
	  write (6,*) resname(ires,imol),i,' ',chain
	  if (chain.eq.'-') then
	    chain= ' '
	    str= ' '
	  else
	    ic= ichar (chain) - 'A'
	    chain= char ('A' + mod (ic,nchainorig))
	    if (nrep.gt.1) str= char ('A' + ic/nchainorig)
	  end if
C	  if (resname(ires,imol)(1:3).eq.'HOH') chain= 'S'
	  write (6,*) 'str= ',str,chain,nchainorig

	  resn= resname(ires,imol)(1:3)
	  call upcase (resn)
	  i= 5
	  if (resn(3:3).eq.' ') then
	    resn= ' ' // resname(ires,imol)(1:2)
	    i= 4
	  end if

	  j= i+1
	  do while (resname(ires,imol)(j:j).ne.' ')
	    j= j + 1
	  end do
	  k= j - i
	  moln= '     '
	  moln(4-k+1:4)= resname(ires,imol)(i:j-1)

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
CHETATM17417 CA    CA A1001    -118.957  85.889  14.889  1.00 28.65          Ca
CATOM  25247 MG   CL2 A1011     -84.633  96.228   4.358  1.00 29.49          MG
CTER   17416      ALA X  35                                                      
CATOM     42 HG11 VAL A   2       5.251  11.378  12.259  1.00  7.64           H 

	  do i= iares(ires,imol),jares(ires,imol)
		iout= iout + 1
		newnumb(i)= iout
		k= 1
		do while (atname(i)(k:k).eq.' ')
		  k= k + 1
		end do
		if (atname(i)(k:k).ge.'0' .and. atname(i)(k:k).le.'9') then
		  if (k.ge.3) then
		    atname1= ' ' // atname(i)(k+1:5) // atname(i)(k:k)
		  else
		    atname1= atname(i)(k+1:5) // atname(i)(k:k)
		  end if
		else
		  if (k.ge.3) then
		    atname1= ' ' // atname(i)(k:5)
		  else
		    atname1= atname(i)(k:5)
		  end if
		end if
		xx= x(i) - rotimo(10,jcopy)
		yy= y(i) - rotimo(11,jcopy)
		zz= z(i) - rotimo(12,jcopy)
	  	x(i)= rotimo(1,jcopy)*xx + rotimo(4,jcopy)*yy + rotimo(7,jcopy)*zz
	  	y(i)= rotimo(2,jcopy)*xx + rotimo(5,jcopy)*yy + rotimo(8,jcopy)*zz
	  	z(i)= rotimo(3,jcopy)*xx + rotimo(6,jcopy)*yy + rotimo(9,jcopy)*zz
		write (2,800) typ,iout,atname1,str,resn,chain,moln,
     $		  x(i),y(i),z(i),weight,bfact(i),ats(i)
	  end do

3333	  continue
	 end do
	end do

        if (typ.eq.'ATOM  ') then
		iout= iout + 1
     		write (2,800) 'TER   ',iout,'    ',' ',resn,chain,moln
        end if

c       ******** connectivity list **********

      do i= 1,natom
	i0= molind(amol(i))-1
	if (ncon(i).ne.0) write (2,850) newnumb(i),
     $		(newnumb(icon(j,i)+i0),j=1,ncon(i))
      end do
      write (2,'(a)') 'END'

850   format ('CONECT',50i5)


800   format (a6,i5,1x,a,a,a,1x,2a,4x,3f8.3,f6.3,f6.1,10x,a)
810   format (a6,i5,6x,a,1x,2a)

	end


	integer function resnumb (resname,imnumb)

	character*(*) resname
	character*1 c
	common nchainorig

	i= 4
	do while (resname(i:i).ne.'-')
	  i= i + 1
	end do
	read (resname(4:i-1),*) j
	c= resname(i+2:i+2)
	ic= ichar (c) - ichar ('A')
	k= mod (ic,nchainorig)
	resnumb= j + k*10000
	imnumb= ic/nchainorig + 1
C	write (6,*) 'resnumb: ',c,ic,nchainorig,imnumb

	return
	end
