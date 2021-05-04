      program Projmod
      implicit none

cKS:======================================================================= 
cKS:  Thu Sep 22 11:25:39 CEST 2005
cKS:  added writing of a binary eigenvector file on first run
cKS:  this will be called upon subsequent runs
cKS:  all modifications limitted by cKS:===================================
cKS:======================================================================= 
c======================================================================= 

c     PROJection of a difference-vector onto a set of vectors (MODES),
c     and/or:
c    (qdepl=T) displacement of the structure along one of them.

c     Vectors: CERFACS format (as produced by DIAGSTD, DIAGRTB, BLZPACK).
c     Difference-vector: difference of two sets of coordinates 
c    (files in PDB format)

c     Since version 1.41:
c     If vectors were calculated for a set of separated bodies (with no
c     interactions), the difference-vector can be buildt from two
c     files with the coordinates of a single of these bodies.

c     Also, vectors are annotated: "A" for motions of the largest body,
c    "B" for motions of the other ones.

c======================================================================= 
      integer nvecmx, natmax, nresmx, nmotsmax

c     NATMAX: Maximum number of atoms.
c     NRESMX: Maximum number of amino-acid residues.
c     NVECMX: Maximum number of vectors.

      parameter(natmax=100000,nresmx=100000,nvecmx=106,
     .        nmotsmax=132)
 
c======================================================================= 
c     A ajouter:
c     La collectivite de chaque mode; celle du changement de conf.
c     Un fichier de commandes.
c     YHS-Fev-1997: Premiere version (Toulouse).
c....................................................................... 
      integer atequiv(natmax), atind(natmax), fatres(nresmx+1), 
     .        i, ii, imax, iresat(natmax), iresatc(natmax), ivec, 
     .        j, jj, jmax, k, 
     .        lmot, lnomeig, lnompdb, lnompdbc, lnompdbr, mdmax, modnum, 
     .        natdif, natom, natomc, natomeff, natvec, nddl, nddleff, 
     .        nmots, nres, nresc, nrsdif, numvec(natmax), nunit, nvec, 
     .        nvecotr, nzero, nzerovp,
     .        uneig, unmor, unout, unpdb, unpdbc, unpdbq, unpdbr
      double precision avemass, bigzero, dfreq, dq, dr, dr2, entr,
     .       freq(3*natmax), gap, lowfreq, massat(natmax), 
     .       matvec(3*natmax,nvecmx), norme, qdiff, qnorm, qp, qproj, 
     .       qtot(nvecmx), recmax, rmsat, rmsmass, small, tot, 
     .       w(nvecmx), wmain(natmax),
     .       xat(natmax), xref(natmax), xconf(natmax), xcurr(natmax),
     .       yat(natmax), yref(natmax), yconf(natmax), ycurr(natmax),
     .       zat(natmax), zref(natmax), zconf(natmax), zcurr(natmax)
      logical qbig, qblock(natmax), qdepl, qdir, qerror, qexist, qinter, 
     .       qmasse, qok, qsubs
      character atonam(natmax)*4, atonamc(natmax)*4, cformat*10, 
     .       cstatus*10, mots(nmotsmax)*132, namfil*64, nom*7, 
     .       nomb*1, nomeig*64, nompdb*64, nompdbc*64, nompdbr*64, 
     .       program*9, progrer*12, progrwn*12, 
     .       resnam(natmax)*4, resnamc(natmax)*4, 
     .       segid(natmax)*4, ssunam(natmax)*4, ssusel*1, version*40
cBegin:
      version=' Version 1.42, February 2005.' 
 
      nom='Projmod'
      program=' '//nom//'>'

      write(6,'(2A)') program,
     .    ' Projection of a difference vector'//
     .    ' on a set of eigenvectors.'
      write(6,'(2A/)') program,version

cDefaults:
      progrwn='%'//nom//'-Wn>'
      progrer='%'//nom//'-Er>'
      small=1e-4
 
c     Ouverture des fichiers:
c     ----------------------
c     En lecture:
 
      call getnam('Name of the file with the (eigen)vectors ?',
     .     nomeig,lnomeig,qok)
      if (.not.qok) stop
 
      cformat="FORMATTED"
      cstatus="old"
 
      nunit=10
      uneig=nunit
      nunit=nunit+1
      call openam(nomeig,cformat,cstatus,uneig,.true.,
     .            qinter,qexist)
      if (qinter.or..not.qexist) stop
 
c     On recherche l'information/sous-unite:
c     Structure de reference:
 
      call getnam('Pdb file with the reference structure ?',
     .     nompdbr,lnompdbr,qok)
      if (.not.qok) stop
 
      nompdb=nompdbr
      lnompdb=lnompdbr
 
      call string_split(nompdb,lnompdb,':',
     .                  mots,nmotsmax,nmots)
      call stringcl(mots(1),lmot)
 
      if (nmots.gt.1) then
          call stringcl(mots(2),lmot)
          ssusel=mots(nmots)
          write(6,*) 'Pdbsel> Subunit to be selected: ',ssusel
          if (nmots.gt.2) then
              write(6,'(4A)') progrwn,
     .      ' The end of pdb name, ',
     .        nompdb(1:lnompdb),', was not understood.'
          endif
      else
          ssusel=' '
      endif
      nompdb=mots(1)
 
      unpdb=nunit
      nunit=nunit+1
      call openam(nompdb,cformat,cstatus,unpdb,.true.,
     .            qinter,qexist)
      if (qinter.or..not.qexist) stop
      unpdbr=unpdb
 
c     Autre conformere (eventuel):
 
      call getnam('Pdb file with the other conformer ?',
     .     nompdbc,lnompdbc,qok)
      if (.not.qok) stop
 
      nompdb=nompdbc
      lnompdb=lnompdbc
 
      call string_split(nompdb,lnompdb,':',
     .                  mots,nmotsmax,nmots)
      call stringcl(mots(1),lmot)
 
      if (nmots.gt.1) then
          call stringcl(mots(2),lmot)
          ssusel=mots(nmots)
          write(6,*) 'Pdbsel> Subunit to be selected: ',ssusel
          if (nmots.gt.2) then
              write(6,'(4A)') progrwn,
     .      ' The end of pdb name, ',
     .        nompdb(1:lnompdb),', was not understood.'
          endif
      else
          ssusel=' '
      endif
      nompdb=mots(1)
 
      unpdbc=-1
      inquire(file=nompdb,exist=qexist)

      if (qexist) then
      unpdb=nunit
      nunit=nunit+1
      call openam(nompdb,cformat,cstatus,unpdb,.true.,
     .            qinter,qexist)
      if (qinter) stop
      unpdbc=unpdb
      endif
 
      call getrep(
     .   ' Are the masses given in the pdb file ? (y/n)',
     .     qmasse,qok)
      if (.not.qok) qmasse=.false.
 
      if (qmasse) then
          write(6,'(/2A)') program,
     .  ' Masses will be picked in the pdb files.'
      else
          write(6,'(/2A)') program,
     .  ' All masses will all be assumed to be of 1.'
      endif
 
      call getrep(
     .    'Displacement along one mode ? (y/n)',
     .     qdepl,qok)
      if (.not.qok) qdepl=.false.
 
      modnum=-1
      qdir=.true.
      if (qdepl) then
          call getchi('Dq ?',dq,qok)
          if (.not.qok) then
              write(6,'(/2A)') progrer,
     .      ' While reading dq. No displacement.'
              qdepl=.false.
          endif
          call getnum('Mode number ?',modnum,1,-1,qok)
          if (.not.qok.or.modnum.le.0) then
            if (unpdbc.gt.0) then
              write(6,'(/2A)') progrwn,
     .      ' Wrong mode number. Displacement along best mode.'
            else
              write(6,'(/2A)') progrer,
     .      ' Wrong mode number. No displacement.'
              stop 
            endif
            modnum=-1
          endif
          if (modnum.lt.0) then
          call getrep(
     .        'Displacement along difference-vector direction ? (y/n)',
     .         qdir,qok)
          if (.not.qok) qdir=.true.
          endif
      endif

c     En ecriture:
 
      if (unpdbc.gt.0) then
      cformat="FORMATTED"
      cstatus="ove"
 
      namfil='projmod.res'
      unout=nunit
      nunit=nunit+1
      call openam(namfil,cformat,cstatus,unout,.true.,
     .            qinter,qexist)
 
      namfil='dr.res'
      unmor=nunit
      nunit=nunit+1
      cstatus="ove"
      call openam(namfil,cformat,cstatus,unmor,.true.,
     .            qinter,qexist)
      endif

      if (qdepl) then
      namfil='projmod_dq.pdb'
      unpdbq=nunit
      nunit=nunit+1
      cstatus="ove"
      call openam(namfil,cformat,cstatus,unpdbq,.true.,
     .            qinter,qexist)
      endif
 
c     Lecture des fichiers:
c     =====================

      call rdmodfacs(uneig,3*natmax,nvecmx,numvec,freq,
     .               matvec,nddl,nvec)
 
      natvec=nddl/3
      write(6,'(/A,I5,A,I6,A)') program, 
     .      nvec,' vectors, ',natvec,' atoms in vector file.'
 
c     Annotation des vecteurs.
c     ------------------------
c     On suppose que le premier vecteur est representatif.

      qbig=.false.
      nzero=0
      do j=1,nddl,3
         jj=(j+2)/3
         qblock(jj)=.false.
         if (dabs(matvec(j,1)).le.small.and.
     .       dabs(matvec(j+1,1)).le.small.and.
     .       dabs(matvec(j+2,1)).le.small) then
             nzero=nzero+1
             qblock(jj)=.true.
         endif
      enddo

      if (nzero.gt.0) then
          write(6,'(/2A)') program,' Several separated bodies.'

c         Le bloc avec le plus de zeros sera le bloc 'A':
          if (nzero.gt.natvec/2) then
              write(6,'(A,I6,A)') 
     .              program,nzero,' atoms in the largest one.'
              qbig=.true.
          else
              write(6,'(A,I6,A)') 
     .        program,natvec-nzero,' atoms in the largest one.'
          endif
      else
          write(6,'(/2A)') program,' Single body modes.'
      endif

c     Les coordonnees:
c     ----------------
      write(6,'(/2A)') program,' First (reference) structure:'

      call rdatompdb(unpdbr,ssusel,xref,yref,zref,massat,
     .     atonam,iresat,resnam,ssunam,segid,natmax,natom,
     .     fatres,nresmx,nres,qerror)

c     Tests:
c     ------

      qsubs=.false.
      if (natvec.eq.natom.or.
     .   (nzero.gt.0.and.(nzero.eq.natom.or.natvec-nzero.eq.natom)))then

          write(6,'(/2A)') program,
     .  ' Cartesian (eigen)vectors will be studied.'

          if (natvec.eq.natom) then
              do i=1,natom
                 atind(i)=i
              enddo
          else
              write(6,'(2A)') program,
     .      ' Motion of a part of the system will be considered.'
              qsubs=.true.
              if (nzero.ne.natom) then
                  do i=1,natvec
                     if (qblock(i)) then
                         qblock(i)=.false.
                     else
                         qblock(i)=.true.
                     endif
                  enddo
              endif
              ii=0
              do i=1,nddl
                 if (qblock(i)) then
                     ii=ii+1
                     atind(ii)=i
                 endif
              enddo
          endif
      else if (natvec.gt.0) then
          write(6,'(/2A)') progrer,
     .  ' Vector and reference pdb file are not consistent.'
          stop '*Wrong files*'
      else if (nddl.le.0.or.nvec.le.0) then
          write(6,'(/2A)') progrer,
     .  ' Nothing can be done.'
          stop '*Wrong files*'
      endif
 
      if (.not.qmasse) then
      do i=1,natom
          massat(i)=1.d0
      enddo
      endif

      if (unpdbc.le.0) then
          write(6,'(/2A)') program,' One conformer only.'
          if (qdepl) then
              goto 500
          else
              write(6,'(/2A)') progrer,' Nothing done.'
              stop
          endif
      endif

      write(6,'(/2A)') program,' Second structure:'

      call rdatompdb(unpdbc,ssusel,xconf,yconf,zconf,wmain,
     .     atonamc,iresatc,resnamc,ssunam,segid,natmax,natomc,
     .     fatres,nresmx,nresc,qerror)

c     Conformite des deux conformeres consideres:
 
      if (natomc.ne.natom.or.nresc.ne.nres) then
      if (natomc.ne.natom) then
          write(6,'(/2A)') progrwn,
     .  ' Different number of atoms for the other conformer.'
      endif
      if (nresc.ne.nres) then
          write(6,'(/2A)') progrwn,
     .  ' Different number of residues for the other conformer.'
      endif
      do i=1,natom
         do j=1,natomc
         if (iresatc(j).eq.iresat(i)) then
         if (atonamc(j).eq.atonam(i).and.resnamc(j).eq.resnam(i)) then
             atequiv(i)=j
             xcurr(i)=xconf(j) 
             ycurr(i)=yconf(j) 
             zcurr(i)=zconf(j) 
             goto 100
         endif
         endif
         enddo
         write(6,'(2A,I6,2A,I6,1X,2A)') progrwn,
     . ' Atom ',i,' of first conformer: ',
     .   resnam(i),iresat(i),atonam(i),' not found in second one.'
         massat(i)=0.d0
 100  continue
      enddo      
      do i=1,natom
         xconf(i)=xcurr(i) 
         yconf(i)=ycurr(i) 
         zconf(i)=zcurr(i) 
         atonamc(i)=atonam(i)
         resnamc(i)=resnam(i)
         iresatc(i)=iresat(i)
      enddo
      write(6,'(2A)') program,
     .  ' All atoms of first conformer were found in second one.'
      endif
 
c     Vecteur difference:
c     ------------------
      natdif=0
      nrsdif=0
      do i=1,natom
 
c        Un des atomes est inconnu ?
         if ((xref(i).gt.9998.and.yref(i).gt.9998.and.
     .        zref(i).gt.9998).or.
     .       (xconf(i).gt.9998.and.yconf(i).gt.9998.and.
     .        zconf(i).gt.9998)) then
              massat(i)=0.d0
         endif

         xat(i)=xconf(i)-xref(i)
         yat(i)=yconf(i)-yref(i)
         zat(i)=zconf(i)-zref(i)

         if (atonamc(i).ne.atonam(i)) then
             natdif=natdif+1
             if (natdif.lt.5) then
                 write(6,'(2A,I6,5A)') progrwn,' Atom ',i,
     .         ' has name ',atonamc(i),' in a file and ',
     .           atonam(i),' in the other.'
             elseif (natdif.eq.5) then
                 write(6,'(2A)') progrwn,' ... '
             endif
         endif
         if (resnamc(i).ne.resnam(i)) then
             nrsdif=nrsdif+1
             if (nrsdif.lt.5) then
                 write(6,'(2A,I6,5A)') progrwn,' Atom ',i,
     .         ' belongs to residue ',resnamc(i),' in a file and ',
     .           resnam(i),' in the other.'
             elseif (nrsdif.eq.5) then
                 write(6,'(2A)') progrwn,' ... '
             endif
         endif
      enddo      

      if (natdif.gt.0)
     .   write(6,'(/A,I6,A)') progrwn,
     .   natdif,' atoms have different names in the two pdb files.'
      if (nrsdif.gt.0)
     .   write(6,'(/A,I6,A)') progrwn,nrsdif,
     . ' atoms belong to different residues in the two pdb files.'
 
      write(6,'(/2A/)') program,
     . ' File dr.res: displacement=f(atom number).'

      avemass=0.d0
      rmsmass=0.d0
      natomeff=0
      rmsat=0.d0

      do i=1,natom
         dr2=0.d0
         if (massat(i).gt.small) then
             avemass=avemass+massat(i)
             rmsmass=rmsmass+massat(i)**2.d0
             dr2=xat(i)**2.d0+yat(i)**2.d0+zat(i)**2.d0
             rmsat=rmsat+dr2
             natomeff=natomeff+1
         endif
         write(unmor,'(I6,F12.4)') i,dsqrt(dr2)
      enddo
 
      write(6,'(A,I6,A)') program,natomeff,' atoms are considered.'
      nddleff=3*natomeff
      if (natomeff.eq.0) then
          write(6,'(/2A)') progrer,
     .  ' All atoms have unknown coordinate(s) or have zero mass.'
          stop
      endif
 
      avemass=avemass/float(natomeff)
      rmsmass=dsqrt(rmsmass/float(natomeff)-avemass**2.d0)
      rmsat=dsqrt(rmsat)
 
      write(6,'(/2A,F8.2)') program,
     .    ' Atomic r.m.s. displacements=  ',
     .      rmsat/sqrt(float(natomeff))
      write(6,'(A,2(A,F8.2))') program,
     .    ' Atomic average masses      =  ',avemass,' +/- ',rmsmass
 
      if (dabs(avemass-1.).gt.small.or.rmsmass.gt.small) then
      norme=0.d0
      do i=1,natom
         if (massat(i).gt.small) 
     .   norme=norme+
     .  (xat(i)**2.d0+yat(i)**2.d0+zat(i)**2.d0)*massat(i)
      enddo
      if (norme.le.small) then
         write(6,'(2A,I6,A)') progrer,
     . ' Difference-vector has null norm. Projection skipped.'
         goto 500
      endif
      norme=dsqrt(norme)

      write(6,'(2A,F8.2)') program,
     .    ' Atomic mass-weighted rmsd  =  ',
     .      norme/sqrt(float(natomeff))
      else
      norme=rmsat
      endif

      if (norme.le.small) then
          write(6,'(/2A)') progrer,' Null difference vector.'
          stop '*No motion*'
      endif

      write(6,'(/2A/)') program,' File projmod.res:'//
     .    ' dr.vector=f(fqcy), and cumulative square sum.'
      
      nvecotr=0
      recmax=0.d0
      mdmax=-1
      tot=0.d0
 
      do ivec=1,nvec
         qnorm=0.d0
         qproj=0.d0
         do i=1,natom
          if (massat(i).gt.small) then
            ii=3*atind(i)-2
            qp=matvec(ii,ivec)*xat(i)*dsqrt(massat(i))+
     .         matvec(ii+1,ivec)*yat(i)*dsqrt(massat(i))+
     .         matvec(ii+2,ivec)*zat(i)*dsqrt(massat(i))
            qproj=qproj+qp
            qnorm=qnorm+
     .         matvec(ii,ivec)**2.d0+
     .         matvec(ii+1,ivec)**2.d0+
     .         matvec(ii+2,ivec)**2.d0
          endif
         enddo
 
c        Cosinus:
c        -------
         if (qnorm.le.small) then
             if (qsubs) then
                 nvecotr=nvecotr+1
             else
                 write(6,'(2A,I6,A)') progrwn,' Eigenvector ',ivec,
     .         ' has null norm. It is skipped.'
             endif
             goto 400
         endif
         if ((natom.eq.natomeff.and.dabs(qnorm-1.).gt.small).or.
     .       qnorm-1.gt.small) then
             write(6,'(2A,I6,A,F12.4)') progrwn,' Eigenvector ',ivec,
     .           ' Norm= ',qnorm
         endif

         qnorm=dsqrt(qnorm)
         qdiff=qproj/qnorm
         qproj=qdiff/norme

         if (dabs(qproj).gt.dabs(recmax)) then
             recmax=qproj
             mdmax=ivec
         endif

c        Somme cumulee des cosinus-carres:
c        ---------------------------------

         qtot(ivec)=qproj**2.d0
         tot=tot+qtot(ivec)

c        Type de vecteur:
c        ----------------

         nomb='A'
         if (nzero.gt.0) then
         k=0
         do j=1,3*natom,3
            jj=(j+2)/3
            if (matvec(j,ivec).eq.0.d0.and.matvec(j+1,ivec).eq.0.d0.and.
     .          matvec(j+2,ivec).eq.0.d0.and.qblock(jj)) k=k+1
         enddo
         if (k.eq.nzero) then
         if (qbig) then
             nomb='A'
         else
             nomb='B'
         endif
         endif
         endif

         write(6,'(A,I6,1X,2A,F8.2,2(4X,A,F6.3),A,F12.4)') 
     .   ' Vector: ',ivec,nomb,' F= ',freq(ivec),
     .   ' Cos= ',qproj,' Sum= ',tot,
     .   ' q= ',qdiff
         write(unout,'(I6,1X,A,1X,3F12.4)') 
     .     ivec,nomb,freq(ivec),qproj,tot

c        Vecteur suivant:
c        ----------------
 400     continue
      enddo
 
      if (nvecotr.gt.0) 
     .   write(6,'(/A,I6,A)') program,nvecotr,
     . ' vectors skipped (they are not involved in subsystem motion).'

c     Compter les modes a frequence "nulle":
c     Chercher le facteur-gap le plus important.

      gap=0.d0
      do i=2,nvec
         if (freq(i-1).eq.0.d0) then
            if (freq(i).gt.0.d0) then
                dfreq=freq(i)-freq(i-1)
                bigzero=dabs(freq(i-1))
                lowfreq=freq(i)
                nzerovp=i-1
                goto 450
            endif
         else
            dfreq=(freq(i)-freq(i-1))/freq(i-1)
            if (dfreq.gt.gap) then
                bigzero=dabs(freq(i-1))
                lowfreq=freq(i)
                nzerovp=i-1
                gap=dfreq
            endif
         endif
      enddo
 450  continue      

      write(6,'(/A,I6,A,F12.6)') program,nzerovp,
     . ' frequencies less than: ',bigzero
      write(6,'(2A,F13.6)') program,
     .' Lowest non-zero frequency  :',lowfreq
 
      write(6,'(/2A,F6.2,A,I6,A,F8.2,A)') program,
     .  ' Best overlap with diff.vect. = ',recmax,' for mode ',mdmax,
     .'   with F= ',freq(mdmax),' cm-1.'

c     On elimine les contributions des frequences nulles:
c     ---------------------------------------------------

      tot=0.d0
      do i=1,nvec
         if (freq(i).lt.lowfreq) then
             tot=tot+qtot(i)
             qtot(i)=0.d0
         endif
      enddo
      if (tot.gt.0.1) then
          write(6,'(/A,F8.2,A)') progrwn,100*tot,
     .   '% of the motion described with zero-frequency modes.'
          if (nzerovp.eq.6)
     .    write(6,'(2A)') progrwn,
     .  ' The second structure was not fitted ?'
      endif

c     Tri des cosinus-carres par ordre decroissant:
c     ---------------------------------------------

      call trier(qtot,nvec,nvecmx,w,.false.)

c     Sommation des contributions (somme des 3N=1):
      tot=0.d0
      do i=1,nvec
         tot=tot+w(i)
         qtot(i)=tot
      enddo

c     Entropie de la description:

      entr=0.d0
      do i=1,nvec
         w(i)=w(i)/tot
         if (w(i).gt.0.d0) entr=entr-w(i)*log(w(i))
      enddo
      entr=exp(entr)

c     Sortie homogene (en colonnes):

      if (nvec.le.nzerovp+60) then
      if (nvec.gt.nzerovp+9) then
         qtot(12)=qtot(nvec)
      else if (nvec.gt.nzerovp+3) then
         qtot(6)=qtot(nvec)
         qtot(9)=qtot(nvec)
      else
         qtot(3)=qtot(nvec)
      endif
      endif

      write(6,'(/2A,6(F5.3,2X))')  program,
     . ' 1-3-6-9-12-all-best contrb.  =  ',
     .   qtot(1),qtot(3),qtot(6),qtot(9),
     .   qtot(12),qtot(nvec)
      
      write(6,'(/2A,F6.1)')  program,
     . ' Effective nb of modes req.   = ',entr

c     Displacement along one mode: 
c     ---------------------------- 

 500  continue
      if (qdepl) then
      if (modnum.le.0) modnum=mdmax
      if (modnum.gt.nvec) then
          write(6,'(/2A,I6)') progrwn,
     .  ' Displacement can not be performed along mode ',modnum
          modnum=nvec
      endif

c     Direction:

      if ((recmax.gt.0.and..not.qdir).or.(recmax.lt.0.and.qdir)) then
          write(6,'(/2A)') progrwn,' Displacement reversed: '//
     .   'along vector-diff. direction, as required.'
          dq=-dq
      endif

      write(6,'(/2A,F12.6,A,I6)') program,
     .  ' Displacement of dq= ',dq,' to be performed along mode ',modnum

      k=0      
      do i=1,natom
         if (massat(i).gt.small) then
            ii=3*atind(i)-2
            xat(i)=xref(i)+(dq*matvec(ii,modnum))/dsqrt(massat(i))
            yat(i)=yref(i)+(dq*matvec(ii+1,modnum))/dsqrt(massat(i))
            zat(i)=zref(i)+(dq*matvec(ii+2,modnum))/dsqrt(massat(i))
         else
            k=k+1
            if (k.le.5) 
     .      write(6,'(2A,I6,A)') progrer,' Atom ',k,
     .    ' has unknown coordinates or zero mass.'
            if (k.eq.5) 
     .      write(6,'(2A)') progrer,' ... '
         endif
      enddo

      if (k.eq.0) then
      call writpdb(unpdbq,xat,yat,zat,massat,
     .     atonam,iresat,resnam,segid,natom)
      else
      write(6,'(/A,I6,A)') progrer,k,
     .   ' atoms have unknown coordinates or zero mass.'
      stop
      endif
      endif

      write(6,'(/2A)') program,' Normal end.'
 
      stop
      end
 
      subroutine getnum(message,numlu,nummin,nummax,qok)
c
c     NUMLU obtenu en reponse au MESSAGE.
c     NUMLU doit etre inferieur a nummax et superieur a nummin.
c
c     NTRYMX essais en cas de probleme.
c     qok=.false. => Probleme a la lecture.
c
c     YHS-oct-1996: version 1.0
c     YHS-nov-2000: version 3.0
c
      implicit none
cI/O:
      integer numlu, nummax, nummin
      logical qok
      character*(*) message
cLocal:
      integer ntry, ntrymx, iread
cBegin:
      ntrymx=5
c
      qok=.false.
      ntry=0
c
 100  continue
      ntry=ntry+1
      if (ntry.ge.ntrymx) return
c
      write(6,'(A,A)') ' Getnum> ',message
      read(5,*,end=200,err=100) iread
      numlu=iread
c
      write(6,*) 'Getnum> ',numlu

      if (nummin.le.nummax) then
      if (numlu.gt.nummax) then
          write(6,'(A,I6,A)') 
     .  '%Getnum-Err: number larger than ',nummax,
     .  ' This is not allowed. Sorry.'
          numlu=nummax
          return
      else if (numlu.lt.nummin) then
          write(6,'(A,I6,A)') 
     .  '%Getnum-Err: number smaller than ',nummin,
     .  ' This is not allowed. Sorry.'
          numlu=nummin
          return
      endif
      endif
c
      qok=.true.
      return
 200  continue
      return
      end 
c-------------------------------------------
      subroutine getchi(message,numlu,qok)
c
c     NUMLU obtenu en reponse au MESSAGE.
c     NTRYMX essais en cas de probleme.
c     YHS-oct-96
c
      implicit none
cI/O:
      double precision numlu
      logical qok
      character*(*) message
cLocal:
      integer ntry, ntrymx
      double precision iread
cBegin:
      ntrymx=5
c
      qok=.false.
      ntry=0
c
 100  continue
      ntry=ntry+1
      if (ntry.ge.ntrymx) return
c
      write(6,'(A,A)') ' Getchi> ',message
      read(5,*,end=200,err=100) iread
      numlu=iread
c
      write(6,*) 'Getchi> ',numlu
c
      qok=.true.
      return
 200  continue
      return
      end 
c---------------------------------------------------
      subroutine getnam(message,nomlu,lnomlu,qok)
c
c     NOMLU obtenu en reponse au MESSAGE.
c     NTRYMX essais en cas de probleme.
c     YHS-oct-96
c
      implicit none
cI/O:
      integer lnomlu
      logical qok
      character*(*) message, nomlu
cLocal:
      integer ntry, ntrymx
cBegin:
      ntrymx=5
c
      qok=.false.
      ntry=0
c
 100  continue
      ntry=ntry+1
      if (ntry.ge.ntrymx) return
c
      write(6,'(A,A)') ' Getnam> ',message
      read(5,'(A)',end=200,err=100) nomlu
c
      call stringcl(nomlu,lnomlu)
      write(6,'(A,A)') ' Getnam> ',nomlu(1:lnomlu)
c
      qok=.true.
      return
 200  continue
      return
      end 
c---------------------------------------------------
      subroutine getrep(message,qinfo,qok)
c
c     qinfo obtenu en reponse au MESSAGE.
c     NTRYMX essais en cas de probleme.
c     YHS-jan-00
c     YHS-oct-00
c
      implicit none
cI/O:
      logical qok, qinfo
      character*(*) message
cLocal:
      integer ntry, ntrymx
      character*1 cread
cBegin:
      ntrymx=2
c
      qinfo=.false.
      qok=.false.
      ntry=0
c
 100  continue
      ntry=ntry+1
      if (ntry.ge.ntrymx) then
          goto 200
          return
      endif
c
      write(6,'(A,A)') ' Getrep> ',message
      read(5,'(A)',end=200,err=100) cread
c
      if (cread.eq.'T'.or.cread.eq.'t'.or.
     .    cread.eq.'Y'.or.cread.eq.'y'.or.
     .    cread.eq.'O'.or.cread.eq.'o') then
          qinfo=.true.
      else
      if (cread.ne.'F'.and.cread.ne.'f'.and.
     .    cread.ne.'N'.and.cread.ne.'n')
     .    write(6,'(3A)') '%Getrep-W> Unexpected answer:',cread,
     . '. Assumed answer is: NO.'
      endif
c
      write(6,*) 'Getrep> ',qinfo
c
      qok=.true.
      return
 200  continue
      write(6,'(A)') '%Getrep-W> No answer.'//
     .    ' Assumed answer is: NO.'
      return
      end 
c
      SUBROUTINE openam(namfil,cformat,cstatus,unit,qverbos,
     .                  qinterr,qexist)
c
c     Ouverture d'un fichier de nom NAMFIL, sur l'unite UNIT.
c
c     input:
c        namfil: nom du fichier a ouvrir. 
c        "stop", "end", "fin", "quit" : arretent le programme.
c        cstatus: mots-cles fortran... ou "OVE" pour overwrite.
c     output: 
c        qexist: flag / existence du fichier 
c        qinterr: Pas de nom pour le fichier cherche.
c
c     YHS-oct-1993: Premiere version.
c     YHS-jan-2000: Derniere modification.
c I/O:
      logical qinterr, qverbos, qexist
      integer unit
      character*10 cformat, cstatus
      character*64 namfil
c Local
      character*132 ordrunix
c begin:
      if (cstatus.eq.'old') cstatus='OLD'
      if (cstatus.eq.'new') cstatus='NEW'
      if (cstatus.eq.'ove') cstatus='OVE'
      if (cstatus.eq.'unknown') cstatus='UNKNOWN'
c
      qinterr=.false.
      qexist=.false.
c
      if (namfil.eq.' ') then 
          qinterr=.true.
          write(6,'(A)') '%Openam-Err> No filename.'
          return
      endif
c
      if (namfil.eq.'stop'.or.namfil.eq.'end'                         
     &    .or.namfil.eq.'fin'.or.namfil.eq.'quit') then 
         write(*,*) 'Openam> Program is stopping on user request.'
         stop                                                                   
      endif 
c
c     Checks if filename is consistent with the opening:
c
      inquire(file=namfil,exist=qexist)
      if (.not.qexist.and.cstatus.eq.'OLD') then
          qinterr=.true.
          write(6,'(A/A)') '%Openam-Err> File not found. Filename: ',
     .    namfil
          return
      endif
c
      if (qexist.and.cstatus.eq.'NEW') then
         write(*,'(/A)') 
     .      '%Openam-Err> This file exists:',namfil
         stop
      else if (qexist.and.cstatus.eq.'OVE') then
         ordrunix='rm '//namfil
         call system(ordrunix)
      endif
      if (cstatus.eq.'OVE') cstatus='NEW'
c                                                                   
      if (qverbos) then
         write(*,'(/A,I6,A)')
     .           ' Openam> file on opening on unit ',unit,':'
         write(*,*) namfil
      endif
      open(file=namfil,form=cformat,
     .     status=cstatus,unit=unit)                
c        
      return                                                                       
      end
c
      subroutine string_split(chaine,taille,delimiteur,
     .                        souschaine,nbremax,nbre)
c
c     "Chaine" est coupee en "nbre" "souschaine" de part et d'autre du
c     "delimiteur"
c      YHS-Sep-93, Uppsala
c I/O:
      integer taille, nbremax, nbre
      character*(*) chaine, souschaine(*), delimiteur
c Local:
      integer icar, iprev
c
      nbre=1
      iprev=1
      souschaine(1)=chaine
      do icar=1,taille
         if (chaine(icar:icar).eq.delimiteur) then
            if (icar-1.ge.iprev) then
               souschaine(nbre)=chaine(iprev:icar-1)
               nbre=nbre+1
               if (nbre.le.nbremax) then
                  if (icar+1.le.taille.and.
     .               chaine(icar+1:taille).ne.' ') then
                     souschaine(nbre)=chaine(icar+1:taille) 
                  else
                     nbre=nbre-1
                     return
                  endif
               else
                  write(6,'(A,I6,A/A)') 
     .               ' %String_split-Err: more than ',nbremax,
     .               ' substrings in : ',chaine
                  return
               endif
            endif
            iprev=icar+1
         endif
      enddo
c
      return
      end
c
      subroutine stringcl(chaine,nonblancs)
c
c     Les caracteres "blancs" de la CHAINE sont retires (a gauche et au milieu).
c     L'entier NONBLANCS donne la position du dernier caractere.
c
c     YHS-Jan-95, Toulouse.
c     YHS-Oct-00, Bordeaux.
c I/O:
      integer nonblancs
      character*(*) chaine
c Local:
      integer icar, ncar, taille
c Begin:
      nonblancs=0
      taille=len(chaine)
      if (taille.le.0) return
c
      if (index(chaine(1:taille),' ').le.0) then
          nonblancs=taille
          return
      endif
c
c*****Nettoyage des blancs a gauche.
c     Premier non-blanc:
c
      do icar=1,taille
         if (chaine(icar:icar).ne.' ') goto 150
      enddo
      icar=taille
 150  continue
      chaine=chaine(icar:taille)
c
c*****Nettoyage des blancs au milieu.
c
          icar=1
          ncar=1
 170      continue
          icar=icar+1
          ncar=ncar+1
          if (chaine(icar:icar).eq.' ') then
              chaine=chaine(1:icar-1)//chaine(icar+1:taille) 
              icar=icar-1
          endif
          if (ncar.lt.taille-1) goto 170
c
      nonblancs=index(chaine,' ')-1
c
      return
      end
c
      subroutine rdmodfacs(uneig,nddlmax,nvecmx,numvec,freq,
     .           matvec,nddl,nvec)
c
c     Lecture de modes CERFACS.
c     Devra remplacer rdcerfacs.
c     Difference: comptage de l'ordre de la matrice.
c    (et pas des atomes)
c
c     Premieres versions (rdcerfacs):
c     YHS-Nov-1996.
c     Dernieres modifications:
c     YHS-Jan-2001.
cKS:======================================================================= 
cKS:  added binary I/O to a binary file proj_modes_bin.IO
cKS:  if this file does not exist it will be created after first reading
cKS:  if this file exists it will be used instead of the ASCII file
cKS:  NOTE: delete this file when using different eigenvector files (NO CHECK)
cKS:======================================================================= 
cI/O:
      integer numvec(*), nvecmx, nddlmax, nvec, nddl, uneig
      double precision freq(*), matvec(nddlmax,*)
cLocal:
      integer nmotsmax
      parameter(nmotsmax=100)
      integer nerr, ivec, indnm_cfacs, nmots,
     .        i, ii, j, jj, k, kk
      double precision wtofreq
      logical qfound, qold, qfirst
      character*1 carnum
      character*132 lign132, mots(nmotsmax)
cKS:======================================================================= 
      logical qbin
cKS:======================================================================= 
cDefaut:
c     Facteur de conversion (2*pi*f)**2 -> f (cm-1):
      wtofreq=108.586698
c
      nerr=0
      qfirst=.true.
      qold=.false.
      qfound=.false.
 100  continue
      read (uneig,'(A)',end=300,err=110) lign132
      goto 120
 110  continue
      nerr=nerr+1
 120  continue
c
      qfound=qfound.or.
     .      (index(lign132,' value ').gt.0.and.
     .       index(lign132,' vector ').gt.0.and.
     .       index(lign132,' residual ').le.0)
      qold=qold.or.
     .      (index(lign132,' VALUE ').gt.0.and.
     .       index(lign132,' VECTOR ').gt.0)
c
      if (.not.qfound.and..not.qold) goto 100
c________________________________________
c
c     Lecture des frequences des modes :
c________________________________________
c
cKS:======================================================================= 
c     test if file proj_modes_bin.IO exists (FLAG qbin)
      write (6, '(A)') "KS> Testing proj_modes_bin.IO"
      qbin = .false.
      open ( 33, file="proj_modes_bin.IO", status="old", 
     .       form="unformatted", err = 991 )
      qbin = .true.
991   continue
      close ( 33 )

c     if ( qbin ) then
      if ( .not. qbin ) goto 666
c     start reading in binary mode

        write (6, '(A)') "KS> Reading modes from proj_modes_bin.IO"
        open ( 33, file="proj_modes_bin.IO", status="old", 
     .         form="unformatted" )
        read (33) nvec, nddl
        write (6, '(A, I6,A)') ">KS", nvec, " modes in file"
        write (6, '(A, I6,A)') ">KS", nddl, " mots"
        write (6, '(A)') "KS> Reading numvec"
        read (33) (numvec(ii), ii=1,nvec)
        write (6, '(A)') "KS> Reading freq"
        read (33) (freq(ii), ii=1,nvec)
        write (6, '(A)') "KS> Reading matvec"
        read (33) ((matvec(ii,ivec), ii=1,nddl),ivec=1,nvec)
        close ( 33 )

c     end reading in binary mode
c     else
      goto 777
666   continue

        write (6, '(A)') "KS> Reading modes in ASCII mode"
cKS:======================================================================= 
      if (qfirst) then
          if (qold) then
          write(6,'(/A)') 
     .  ' Rdmodfacs> Old Blzpack file format detected.'
          else
          write(6,'(/A)') 
     .  ' Rdmodfacs> Blzpack file format detected.'
          endif
          qfirst=.false.
      endif
c
      ivec=0
      nvec=0
 250  continue
      ivec=ivec+1
      if (ivec.gt.nvecmx) then
          write(6,'(/A,I5,A)') 
     .  '%Rdmodfacs-Err> More than ',nvecmx,' vectors in file.'
          return
      endif
c
      read(lign132,'(7X,I5,12X,G12.4)',end=240,err=240)
     .     numvec(ivec), freq(ivec)
      freq(ivec)=wtofreq*dsqrt(abs(freq(ivec)))
c
      goto 255
 240  continue
      write(6,'(/3A)')
     .    '%Rdmodfacs-W> Pb with ligne: ',lign132(1:36),'...'
 255  continue
c
      nvec=ivec
      write(6,'(/A,I6)')
     .    ' Rdmodfacs> Numero du vecteur CERFACS en lecture:',
     .      numvec(ivec)
      write(6,'(A,1PG12.4)')
     .    ' Rdmodfacs> Frequence du vecteur en lecture:',
     .      freq(ivec)
c
      if (numvec(ivec).le.0)
     .    write(6,'(/A/A)')
     .  '%Rdmodfacs-W> Vector number was expected in:',
     .    lign132
c
      read(uneig,'(A)',end=230,err=230) lign132
 230  continue
      read(lign132,'(1X,A1)',end=232,err=232) carnum
 232  continue
      if ((qfound.and.carnum.ne.'=').or.
     .    (qold.and.carnum.ne.'-')) then
          write(6,'(2A/A)')
     .       ' %Rdmodfacs-Warning> Unexpected character ',
     .       ' in second column of line:',
     .    lign132
      endif
c____________________________________________________
c
c     2) Lecture des coordonnees des modes CERFACS :
c        Format libre.
c____________________________________________________
c
      k=0
 257  continue
      if (k.ge.nddlmax) then
          write(6,'(/A,I6,A,I5)') 
     .  '%Rdmodfacs-Err> More than ',nddlmax,
     .  ' coordinates for vector ',ivec
          return
      endif
c
      read(uneig,'(A)',end=300,err=270) lign132
c
c     Nombre de coordonnees par ligne:
      call string_split(lign132,132,' ',
     .                  mots,nmotsmax,nmots)
c
      if (lign132.eq.' ') then
          read(uneig,'(A)',end=300,err=260) lign132
      else if (.not.qold.or.index(lign132,' VALUE ').le.0) then
          read(lign132,*,end=258)
     .   (matvec(k+ii,ivec),ii=1,nmots)
          k=k+nmots
          goto 257
 258      continue
      endif
      nddl=k
c
 260  continue
      indnm_cfacs=index(lign132,'       VALUE')
      if (indnm_cfacs.le.0)
     .    indnm_cfacs=index(lign132,'       value')
      if (indnm_cfacs.gt.0) then
          goto 250
      else
          write(6,'(A,A/A/A)')
     .  ' Rdmodfacs: Lecture des modes CERFACS terminee.',
     .  ' Item VALUE non trouve dans la ligne:',lign132
          goto 300
      endif
c
 270  continue
      write(6,'(A,I6,A)')
     .   ' %Rdmodfacs-Error: durant la lecture de la coordonnee ',
     .      i,' du mode.'
      stop
c
 220  continue
c*****Ligne suivante de la lecture du fichier des modes en cours :
c
      goto 100
c
c     Fin de la lecture du fichier des modes :
c
 300  continue
cKS:======================================================================= 
c     end reading in ASCII mode
c     start writing in binary mode
c     integer numvec(*), nvecmx, nddlmax, nvec, nddl, uneig
c     double precision freq(*), matvec(nddlmax,*)

        write (6, '(A)') "KS> Writing modes to proj_modes_bin.IO"
        write (6, '(A, I6,A)') ">KS", nvec, " modes to write to file"
        write (6, '(A, I6,A)') ">KS", nddl, " mots"
        open ( 33, file="proj_modes_bin.IO", status="new", 
     .         form="unformatted" )
        write (33) nvec, nddl
        write (6, '(A)') "KS> Writing numvec"
        write (33) (numvec(ii), ii=1,nvec)
        write (6, '(A)') "KS> Writing freq"
        write (33) (freq(ii), ii=1,nvec)
        write (6, '(A)') "KS> Writing matvec"
        write (33) ((matvec(ii,ivec), ii=1,nddl),ivec=1,nvec)
        close ( 33 ) 
        write (6, '(A)') "KS> Writing modes to proj_modes_bin.IO done"

c     end writing in binary mode

c     end if 
777   continue
cKS:======================================================================= 


      return
      end
c----------------------------------------------------------------------
      subroutine trier(y,npoint,nmax,ysort,qcrois)
c
c     Tri par ordre croissant ou decroissant.
c     YHS-Jun-2002: premiere version (Bordeaux).

      implicit none
      integer i, j, nmax, npoint
      logical qcrois
      double precision y(*), ycur, ysort(*)
      character progrer*10

      progrer='%Trier-Er>'

      if (npoint.gt.nmax) then
          write(6,'(A,I9,A,I9,A)') progrer,npoint,
     .  ' points to be sorted, i.e., more than ',nmax,' Sorry.'
          stop
      endif

      do i=1,npoint
         ysort(i)=y(i)
      enddo

      if (qcrois) then
      do i=1,npoint
         do j=1,npoint
            if (ysort(i).lt.ysort(j)) then
                ycur=ysort(i)
                ysort(i)=ysort(j)
                ysort(j)=ycur
            endif
         enddo
      enddo
      else
      do i=1,npoint
         do j=1,npoint
            if (ysort(i).gt.ysort(j)) then
                ycur=ysort(i)
                ysort(i)=ysort(j)
                ysort(j)=ycur
            endif
         enddo
      enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine writpdb(unpdb,xat,yat,zat,binfo,
     .           atonam,ires,resnam,segid,natom)
c
c     Ecriture d'un fichier en format pdb.
c     YHS-Mars-1996.
c
      implicit none
cI/O:
      integer unpdb, ires(*), natom
      double precision xat(*), yat(*), zat(*), binfo(*)
      character*4 atonam(*), resnam(*), segid(*)
cLocal:
      integer i, j, k
      character*5 atom
cBegin:
      write(6,'(/A)') ' Writpdb> Writing pdb file.'
c
      do i=1,natom
c
c     Un petit probleme avec le cadrage des noms des atomes:
c    'iHjj' ou ' NNN'
c
      atom=atonam(i)
      if (index(atonam(i),' ').gt.1) atom=' '//atonam(i)
c
c     write(unpdb,'(A,I7,1X,A4,A4,2X,I4,4X,3F8.3,6X,F6.2,6X,A)') 
      write(unpdb,'(A,I7,1X,A4,1X,A4,1X,I4,4X,3F8.3,6X,F6.2,6X,A)') 
     .      'ATOM', i, atom(1:4), resnam(i), ires(i), 
     .      xat(i), yat(i), zat(i), binfo(i), segid(i)
      enddo
      write(unpdb,'(A)') 'END'
c
      write(6,'(A,I6,A)') ' Writpdb> ',natom,' atoms saved.'
      return
      end
c-----------------------------------------------------------------------
      subroutine rdatompdb(unpdb,ssusel,xat,yat,zat,binfo,
     .           atonam,iresat,resnam,ssunam,segid,natmax,natom,
     .           fatres,nresmx,nres,qerror)
c
c     Lecture ligne a ligne d'un fichier pdb.
c     Uniquement les lignes commencant par 'ATOM'.
c     Uniquement ceux de la sous-unite selectionnee.
c
c     fatres(i): numero du premier atome du residu i.
c
c     YHS-nov-1996: premiere version (Toulouse).
c     YHS-mar-2003: derniere version (Lyon).
c
      implicit none
cI/O:
      integer unpdb, natmax, iresat(*), natom,  lnom,
     .        nresmx, nres, fatres(*)
      double precision xat(*), yat(*), zat(*), binfo(*)
      logical qerror
      character*4 atonam(*), resnam(*), segid(*)
      character*1 ssusel, ssunam(*)
cLocal:
      integer nerr, iatom, irs, irsprev,
     .        i, j, k, ii
      double precision x, y, z, bfact
      character*1  ssu
      character*4  ren, segat
      character*5  atncur
      character*80 lign80
cBegin:
      write(6,'(/A)') ' Rdatompdb> Reading pdb file.'
c
      qerror=.false.
      nerr=0
c
      irsprev=-1
      nres=0
      iatom=0
 105  continue   
      read(unpdb,'(A)',end=200,err=110) lign80 
c 
      goto 120                                
 110  continue
      nerr=nerr+1                            
c
 120  continue                              
      if (lign80(1:4).eq.'ATOM') then
      read(lign80,'(12X,A4,1X,A4,A1,I4,4X,3F8.3,6X,F6.2,6X,A4)') 
     .            atncur, ren, ssu, irs, x, y, z, 
     .            bfact, segat 
      if (iatom.lt.natmax) then
          if (ssu.eq.ssusel.or.ssusel.eq.' ') then
          iatom=iatom+1
          xat(iatom)=x
          yat(iatom)=y
          zat(iatom)=z
          binfo(iatom)=bfact
c
          call stringcl(atncur,lnom)
          atonam(iatom)=atncur
          call stringcl(ren,lnom)
          resnam(iatom)=ren
          iresat(iatom)=irs
          ssunam(iatom)=ssu
          segid(iatom)=segat
c
          if (irs.ne.irsprev) then
              nres=nres+1
              if (nres.gt.nresmx) then
                  write(6,'(A/A,I6)') 
     .          '%Rdatompdb-Er> Too many residues in this file.',
     .          ' Maximum allowed is = ',nresmx
                  stop
              endif
              irsprev=irs
              fatres(nres)=iatom
          endif
          endif
      else
          write(6,'(A/A,I6)') 
     .      '%Rdatompdb-Er> Too many atoms in this file.',
     .      ' Maximum allowed is = ',natmax
          stop
      endif
      endif
c
c     2) Ligne suivante du fichier pdb :
c
      goto 105
c
c     3) Fin de la lecture du fichier pdb :
c
 200  continue 
      write(6,*) 'Rdatompdb> End of file reached.'
      write(6,*) 'Rdatompdb> Number of I/O errors: ',
     .            nerr
c
      natom=iatom
      fatres(nres+1)=natom+1
      irs=0
      if (natom.gt.0) irs=iresat(natom)
c
      write(6,'(/(A,I6))') 
     .' Rdatompdb> Number of residues found = ',nres,
     .'            First residue number     = ',iresat(1),
     .'            Last  residue number     = ',irs,
     .'            Number of atoms found    = ',natom
      write(6,'(A,F8.1)') 
     .'            Mean number per residue  = ',float(natom)/float(nres)
c
      if (natom.eq.0) then
          write(6,'(A)')
     .  '%Rdatompdb-Er> No atom found in file.'
          qerror=.true.
      endif
      if (nres.eq.0) then
          write(6,'(A)')
     .  '%Rdatompdb-Er> No residue found in file.'
          qerror=.true.
      endif
c
      return
      end
