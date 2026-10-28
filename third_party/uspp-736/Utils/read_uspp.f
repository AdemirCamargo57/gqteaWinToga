c     *********************************************************
      program read_uspp
c     *********************************************************
c
c     Useless program for reading in a binary unformatted
c       pseudopotential file of version 732 or later
c
c     It can form the basis for a conversion program
c       or simply serve as a summary of the file format
c
c     ----------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
      parameter(idim1=1000,idim3=10,idim4=5,idim5=4,idim6=2*idim5-1)
      parameter(idim8=20)
c
      dimension r(idim1),rab(idim1)
      dimension rsatom(idim1),rspsco(idim1)
      dimension nnlzps(idim4),eeps(idim4),wwnlps(idim4),rc(idim5)
      dimension wf(idim1,idim4)
      dimension beta(idim1,idim3),qfunc(idim1,idim3,idim3),
     +  qfcoef(idim8,idim6,idim3,idim3),vloc(idim1),vloc0(idim1),
     +  rinner(idim6)
      dimension lll(idim3),eee(idim3),iptype(idim3),qqq(idim3,idim3),
     +  ddd0(idim3,idim3),ddd(idim3,idim3)
      dimension iver(3),idmy(3)
c
      character*20 title
c
      data iops /14/
c
c     ----------------------------------------------------------
      open (iops,file='sample.uspp',form='unformatted',status='old')
c     ----------------------------------------------------------
c
      read (iops) (iver(i),i=1,3),(idmy(i),i=1,3)
      read (iops) title,zps,zvps,exftps,nvalps,mesh,etot
      read (iops) (nnlzps(i),wwnlps(i),eeps(i),i=1,nvalps)
      read (iops) keyps,ifpcor,rinner(1)
      read (iops) nang,lloc,eloc,ifqopt,nqf,qtryc
      read (iops) (rinner(i),i=1,nang*2-1)
      read (iops) irel
      read (iops) (rc(i),i=1,nang)
      read (iops) nbeta,kkbeta
c
      do j=1,nbeta
        read (iops) lll(j),eee(j),(beta(i,j),i=1,kkbeta)
        do k=j,nbeta
          read (iops) ddd0(j,k),ddd(j,k),qqq(j,k),
     +    (qfunc(i,j,k),i=1,kkbeta),
     +    ((qfcoef(i,lp,j,k),i=1,nqf),lp=1,nang*2-1)
        end do
      end do
c
      read (iops) (iptype(j),j=1,nbeta),npf,ptryc
      read (iops) rcloc,(vloc0(i),i=1,mesh)
c
      if (ifpcor.gt.0) then
         read (iops) rpcor
         read (iops) (rspsco(i),i=1,mesh)
      endif
c
      read (iops) (vloc(i),i=1,mesh)
      read (iops) (rsatom(i),i=1,mesh)
      read (iops) (r(i),i=1,mesh)
      read (iops) (rab(i),i=1,mesh)
      read (iops) nchi
      read (iops) ((wf(i,j), i=1,mesh),j=1,nchi)
c
c     ----------------------------------------------------------
      write (6,'(a)') ' Pseudopotential successfully read'
c     ----------------------------------------------------------
c
      end
