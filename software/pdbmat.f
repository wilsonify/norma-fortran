      program pdbmat

c=======================================================================

c     Computes the mass-weighted second derivatives energy matrix,
c     using Tirion's model, that is, an elastic network model.
c     In such models, all neighbors are linked by springs.

c     Goal: Computing the normal modes of vibration of the system.

c     To do so, the matrix produced by pdbmat has to be diagonalized
c     with another program, like diagstd or, for large systems, diagrtb.

c=======================================================================

c     INPUT: 
c     ******
c     A parameter file, called pdbmat.dat.  

c     Note that each run of pdbmat produces a pdbmat.dat_run file,
c     where all parameter values are given (and commented).
c     pdbmat.dat_run can be modified and used as a pdbmat.dat file, 
c     for further runs. 

c     Among the parameters: 
c     The name of a file with the coordinates of 
c     the system, in FREE or PDB (protein data bank) format.

c     Free format: x, y, z, mass.
c     PDB  format: Masses can be given in the Bfactors column.
c     Note that masses are expected, but not required.

c     OUTPUT:
c     *******
c     A matrix in format: i, j, non-zero-i-j-matrix-element

c     Output matrix filename:
c     Formatted file: pdbmat.sdijf (Free format)
c     Binary    file: pdbmat.sdijb 

c     A coordinate file with: x, y, z, mass, block-number
c     For a pdb file, the block-number is the amino-acid number.
c     It is of no use in the present program.

c     Output coordinate filename:
c     pdbmat.xyzm

c.......................................................................

c     MEMORY LIMITS:
c     **************

c     NATMAX: Maximum number of atoms allowed.
c     NRESMX: Maximum number of residues allowed.

      implicit none
      integer natmax, nresmx

      parameter(NATMAX=50000,NRESMX=50000)

c.......................................................................

c     ABOUT Tirion's model:
c     *********************

c     Principe du modele (Tirion, 1996): 
      
c     Tous les atomes a moins de "cutoff" les uns des autres 
c     sont supposes lies par des ressorts, qui ont tous
c     la meme raideur.
c     Simplification supplementaire par rapport au modele initial: 
c     les atomes sont supposes avoir tous la meme taille
c    (le cutoff est le meme pour toutes les paires d'atomes).
c     On peut de plus poser qu'ils ont tous la meme masse.
c     Sinon, celles-ci sont lues dans la colonne des 
c     facteurs B du fichier pdb.
 
c     Principal resultat:
 
c     Les modes de vibration de basse frequence obtenus
c     a partir d'un tel modele sont tres voisins de ceux
c     obtenus avec un modele beaucoup plus detaille, tels
c     ceux utilises lors des etudes de Dynamique Moleculaire.
 
c     Principaux avantages:
 
c     Pas besoin de prendre en compte tous les atomes.
c     Pas besoin de minimisation d'energie, prealablement
c     au calcul des modes de vibration (E=0 par construction).

c.......................................................................

c     REFERENCES:
c     ***********

c     1) M.M. Tirion (1996):
c    "Large amplitude elastic motions in proteins from
c     a single-parameter, atomic analysis",
c     Phys. Rev. letters vol.77(9), p1905-1908.
 
c     2) F. Tama, Y.H. Sanejouand (2001):
c    "Conformational change of proteins arising 
c     from normal modes calculations"
c     Protein Engineering vol.14, p1-6.

c     In case of problem, feel free to contact: 
c     Yves-Henri.Sanejouand@ens-lyon.fr

c.......................................................................
c     YHS-Nov-1996: Version initiale (Toulouse).
c.......................................................................
      integer nmotsmax, ntopmax
      parameter(ntopmax=10*natmax,nmotsmax=100)
 
      integer fatres(nresmx+1), i, idmax, idres(nresmx), ii, 
     .        imax, imin, ires, iresat(natmax), 
     .        j, jangle(ntopmax), jat, jbond(ntopmax), jj,
     .        k, kcom, kk, klist, l, ll, lmot, lnom, lnompdb, 
     .        namax, namin, nangle(ntopmax), nangles, natom, 
     .        nbig, nbmax, nbmin, nbond(ntopmax), nbonds, ndat, nl,
     .        nmax, nmin, nmots, nntr, nnzero, nres, nunit, nunknown, 
     .        nvoisat(natmax), prtlev, uninp, unout, unpdb
      double precision cutbnd, cutoff, ddf, der2(3,3*natmax), 
     .        dist, dist2, elemnt, elmax, kangle, kbond, kij, knonb, 
     .        levelshft, massat(natmax), nmoy, nrms,
     .        rave, rbig, rdev, rinput, rmax, rmin, rsmall, rx, ry, rz, 
     .        trace, unknown, xat(natmax), yat(natmax), zat(natmax) 
      logical qbinary, qerror, qexist, qinter, qmasse, qok, qpdb
      character atonam(natmax)*4, cformat*32, cstatus*32, lign80*80, 
     .        motinp*80, mots(nmotsmax)*132, nomfich*80, 
     .        nompdb*80, program*8, progrer*11, progrwn*11, 
     .        residus_standards*132, residus_stshort*21, 
     .        resnam(natmax)*4, ssusel*1, typbond*80, typmas*80, 
     .        typout*80, version*32
      parameter(rbig=1e10,rsmall=1e-10,unknown=9999.d9)
c.......................................................................
      version=' Version 3.46, Fevrier 2004.'
 
cBegin:
      idmax=21
      residus_standards='   ILE PHE TRP LEU CYS VAL MET TYR ALA HIS '//
     .                     'GLY THR SER PRO ARG GLN ASN ASP GLU LYS '
      residus_stshort='IFWLCVMYAHGTSPRQNDEKX'
 
      program=' Pdbmat>'
      progrer='%Pdbmat-Er>'
      progrwn='%Pdbmat-Wn>'

      write(6,'(2A)') program,
     .' Computes the Hessian matrix, using an Elastic Network Model.'
      write(6,'(2A)') program,version
c____________________________________________________
c
c     Ouverture et Lecture du fichier d'instructions:
c____________________________________________________
 
c     Valeurs par default:
      nompdb='pdbmat.structure'
      typbond='NONE'
      typmas='CONS'
      typout='  FREE'
      kbond=1000.0d0
      cutbnd=2.d0
      kangle=100.0d0
      knonb=1.0d0
      cutoff=12.d0
      levelshft=1e-9
      qmasse=.false.
      prtlev=0
      qbinary=.false.
 
      nunit=10
 
      uninp=nunit
      nunit=nunit+1
      nomfich='pdbmat.dat'
      cformat="FORMATTED"
      cstatus="old"
      call openam(nomfich,cformat,cstatus,uninp,.false.,
     .     qinter,qexist)
      if (qinter.or..not.qexist) then 
          write(6,'(/2A)') progrwn,
     .  ' No pdbmat.dat file found. Defaults assumed for all options. '
          goto 110
      else
          write(6,'(/2A)') program,
     .  ' Options to be read in pdbmat.dat file.'
      endif
 
 50   continue
      read(uninp,'(A)',end=100) lign80
 
      k=index(lign80,'=') 
      motinp=' '
      if (k.gt.0) then 
          motinp=lign80(1:k)
      else
          kcom=index(lign80,'!')
          if (kcom.ne.1) then
          write(6,'(/2A/A)') progrwn,
     .  ' No separator (=) in command ligne:',
     .    lign80
          write(6,'(2A)') progrwn,' This ligne is skipped.'
          endif
          goto 50
      endif
      call mintomaj(motinp)

      kcom=index(lign80(k+1:80),'!')
      if (kcom.gt.0) lign80(k+kcom:80)=' '
      klist=index(lign80(k+1:80),'?')

      if (index(motinp,' FILENAME').gt.0) then 
          nompdb=lign80(k+1:80)
      else if (index(motinp,' DEFINITION').gt.0) then
          typbond=lign80(k+1:80)
          call mintomaj(typbond)
          call stringcl(typbond,lnom)
          if (typbond(1:3).eq.'ALL') then 
              typbond=' ALL'
          else if (typbond(1:3).eq.'NON') then 
              typbond='NONE'
          else if (typbond(1:3).eq.'CON') then
              typbond='CONS'
          else
              write(6,'(/3A)') progrwn,' Bond definition :',
     .        typbond(1:4)
              if (klist.le.0)
     .        write(6,'(2A)') progrwn,' This is not a known keyword.'
              write(6,'(2A)') progrwn,
     .      ' Valid options are: NONe, ALL, CONsecutive.'
              write(6,'(A)') ' Default assumed.'
              typbond='NONE'
          endif
      else if (index(motinp,'MASS').gt.0) then
          typmas=lign80(k+1:80)
          call mintomaj(typmas)
          call stringcl(typmas,lnom)
          if (typmas(1:3).eq.'PDB'.or.typmas(1:3).eq.'COO') then
              qmasse=.true.
              typmas='COOR'
          else if (typmas(1:3).ne.'CON') then
              write(6,'(/3A)') progrwn,' Origin of mass values :',
     .        typmas(1:3)
              if (klist.le.0)
     .        write(6,'(2A)') progrwn,' This is not a known keyword.'
              write(6,'(2A)') progrwn,
     .      ' Valid options are: CONstant, COOr, PDB.'
              write(6,'(A)') ' Default assumed.'
              qmasse=.false.
              typmas='CONS'
          endif
      else if (index(motinp,'FORMAT').gt.0) then
          typout=lign80(k+1:80)
          call mintomaj(typout)
          call stringcl(typout,lnom)
          if (typout(1:1).eq.'B'.or.typout(1:1).eq.'U') then
              qbinary=.true.
              typout='BINARY'
          else if (typout(1:1).ne.'F') then
              write(6,'(/3A)') progrwn,' Kind of matrix format :',
     .        typout(1:1)
              if (klist.le.0)
     .        write(6,'(2A)') progrwn,' This is not a known keyword.'
              write(6,'(2A)') progrwn,
     .      ' Valid options are: Free, Binary, Formatted, Unformatted.'
              write(6,'(A)') ' Default assumed.'
              qbinary=.false.
              typout='  FREE'
          else
              qbinary=.false.
              typout='  FREE'
          endif
      else 
          qok=.false.
          read(lign80(k+1:80),*) rinput
          if (index(motinp,'SHIFT ').gt.0) then 
               qok=.true.
               levelshft=rinput
          else if (index(motinp,' CUTOF').gt.0.or.
     .             index(motinp,' DIST').gt.0) then
               qok=.true.
               cutoff=rinput
          else if (index(motinp,'INTERAC').gt.0) then
               if (index(motinp,' FORCE ').gt.0.or.
     .             index(motinp,' CONST').gt.0) then
                   qok=.true.
                   knonb=rinput
               endif
               if (index(motinp,' CUTOF').gt.0.or.
     .             index(motinp,' DIST').gt.0) then
                   qok=.true.
                   cutoff=rinput
               endif 
          else if (index(motinp,' BOND').gt.0.and.
     .        (index(motinp,' FORCE ').gt.0.or.
     .         index(motinp,' CONST').gt.0)) then
               qok=.true.
               kbond=rinput
          else if (index(motinp,' LENGTH').gt.0) then
               qok=.true.
               cutbnd=rinput
          else if (index(motinp,'PRINT').gt.0) then
               qok=.true.
               prtlev=int(rinput)
          else if (index(motinp,'ANGLE').gt.0) then
               if (index(motinp,' FORCE ').gt.0.or.
     .             index(motinp,' CONST').gt.0) then
                   qok=.true.
                   kangle=rinput
               endif
          endif
          if (.not.qok) then
               write(6,'(/2A/A)') progrwn,' No known keyword in ligne:',
     .         motinp
               write(6,'(2A/)') progrwn,
     .       ' This command ligne is skipped.'
          endif
      endif
      goto 50
 
 100  continue
      close(uninp)
 110  continue
      call stringcl(nompdb,lnompdb)
 
      write(6,'(/3A)') program,' Coordinate filename     = ',
     .      nompdb(1:lnompdb)
      write(6,'(/(A,F10.2))') 
     .' Pdbmat> Distance cutoff         = ',cutoff,
     .'         Force constant          = ',knonb
      if (typbond.ne.'NONE') then
      write(6,'(A,6X,A)') 
     .'         Kind of bond definition = ',typbond(1:4)
      write(6,'(A,F10.2)') 
     .'         Maximum bond length     = ',cutbnd,
     .'         Bond force constant     = ',kbond,
     .'         Angle force constant    = ',kangle
      endif
      write(6,'(A,6X,A)') 
     .'         Origin of mass values   = ',typmas(1:4)
      if (prtlev.gt.0) then
      write(6,'(A,1PG10.1)') 
     .' Pdbmat> Levelshift              = ',levelshft
      write(6,'(A,3X,I7)') 
     .'         PRINTing level          = ',prtlev
      endif

c     Sauvegarde du fichier de commandes complet:

      uninp=nunit
      nunit=nunit+1
      nomfich='pdbmat.dat_run'
      cformat="FORMATTED"
      cstatus="ove"
      call openam(nomfich,cformat,cstatus,uninp,.false.,
     .     qinter,qexist)

c     Pour etre plus clair:

      if (typbond.eq.'NONE') then
          cutbnd=0.d0
          kangle=0.d0
          kbond=0.d0
      endif

      write(uninp,'(2A)') 
     .'! This file can be modified and used as a command file',
     .' (named pdbmat.dat) for pdbmat.'
      write(uninp,'(2A)') ' Coordinate FILENAME        = ',
     .      nompdb(1:lnompdb)
      write(uninp,'(A,F10.3)') ' INTERACtion DISTance CUTOF = ',cutoff
      write(uninp,'(A,F10.3)') ' INTERACtion FORCE CONStant = ',knonb
      write(uninp,'(A,6X,2A)') ' Origin of MASS values      = ',
     .      typmas(1:4),' ! CONstant, or from COOrdinate file.'
      write(uninp,'(A,8X,I2,A)') ' Output PRINTing level      = ',
     .      prtlev,' ! =1: more detailled. =2: debug level.'
      write(uninp,'(A,6X,2A)') ' Bond DEFINITION            = ',
     .      typbond(1:4),' ! NONe, ALL, or between CONsecutive atoms.'
      write(uninp,'(A,F10.3)') ' Maximum bond LENGTH        = ',cutbnd
      write(uninp,'(A,F10.3)') ' BOND FORCE CONStant        = ',kbond
      write(uninp,'(A,F10.3)') ' ANGLE FORCE CONStant       = ',kangle
      write(uninp,'(A,1PG10.1,A)') ' LevelSHIFT                 = ',
     .      levelshft,
     .  ' ! Non-zero value often required (numerical reasons).'
      write(uninp,'(A,4X,2A)') ' Matrix FORMAT              = ',
     .      typout(1:6),' ! Free, or Binary, matrix saved.'
      close(uninp)

c     Test:

      if (cutoff.lt.0.d0.or.knonb.lt.0.d0.or.
     .   (typbond(1:4).ne.'NONE'.and.(cutbnd.lt.0.d0.or.kbond.lt.0.d0
     .   .or.kangle.lt.0.d0))) then
          write(6,'(/2A)') progrer,
     .  ' Distances and force constants can not be negative values !' 
          stop
      endif

c     On recherche l'information/sous-unite:
 
      call string_split(nompdb,lnompdb,":",
     .                  mots,nmotsmax,nmots)
      call stringcl(mots(1),lnom)
 
      if (nmots.gt.1) then
          call stringcl(mots(2),lnom)
          ssusel=mots(nmots)
          write(6,*) 'Pdbmat> Subunit to be selected: ',ssusel
          if (nmots.gt.2) then
              write(6,'(A)') 
     .      '%Pdbmat-W> The end of filename, ',
     .        nompdb(1:lnompdb),', was not understood.'
          endif
      else
          ssusel=' '
      endif
      nompdb=mots(1)
      call stringcl(nompdb,lnompdb)
c_______________________________________
c                                
c     Lecture du fichier de coordonnees:
c_______________________________________
      if (prtlev.gt.0)
     .  write(6,'(/(3A))') 
     .' Pdbmat> Coordinate file ',nompdb(1:lnompdb),' to be opened.'

      unpdb=nunit
      nunit=nunit+1
      cformat="FORMATTED"
      cstatus="old"
      call openam(nompdb,cformat,cstatus,unpdb,.true.,
     .     qinter,qexist)
      if (qinter) stop '*Input error*'

c     Format pdb ?

      nl=0
      qpdb=.false.
 120  continue
      read(unpdb,'(A)',end=130) lign80
      if (lign80(1:5).eq.'ATOM '.or.lign80(1:6).eq.'HETATM ') then
          qpdb=.true.
          goto 130
      else
          nl=nl+1
      endif
      goto 120
 130  continue
      rewind(unpdb)

      do i=1,natmax
         xat(i)=unknown
         yat(i)=unknown
         zat(i)=unknown
         massat(i)=unknown
         iresat(i)=i
      enddo

      if (qpdb) then
          write(6,'(/2A)') program,
     .  ' Coordinate file in PDB format.'
          call rdatompdb(unpdb,ssusel,xat,yat,zat,massat,
     .         atonam,iresat,resnam,natmax,natom,
     .         fatres,nresmx,nres,qerror,prtlev)
      else
          if (nl.eq.0) then
              write(6,'(/2A)') progrer,' Empty coordinate file.'
              stop
          endif
          write(6,'(/2A)') program,
     .  ' Coordinate file in Free format.'

          call readxyz(unpdb,xat,yat,zat,massat,iresat,natmax,natom,
     .         ndat,qerror,prtlev)

          if (qmasse.and.ndat.lt.4) then
              write(6,'(/2A)') progrer,
     .      ' Masses were not all found, as expected.'
              qmasse=.false.
          endif
      endif

c     Tests:

      if (qerror) stop
      if (natom.le.1) then
          write(6,'(2A)') progrer,
     .  ' Not enough atoms found in file. Nothing done.'
          stop
      endif

      write(6,'(/2A)') program,' Coordinate statistics: '

      call vecstat(xat,natom,rmin,rmax,rave,rdev)
      write(6,'(4(A,F12.6))') 
     .' <x>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

      call vecstat(yat,natom,rmin,rmax,rave,rdev)
      write(6,'(4(A,F12.6))') 
     .' <y>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

      call vecstat(zat,natom,rmin,rmax,rave,rdev)
      write(6,'(4(A,F12.6))') 
     .' <z>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

      if (qmasse) then
          write(6,'(/2A)') program,' Mass statistics: '
          call vecstat(massat,natom,rmin,rmax,rave,rdev)
          write(6,'(4(A,F12.6))') 
     .  ' <m>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

          if (rmin.le.0.d0) then
              write(6,'(2A)') progrer,
     .      ' Negative or null masses found !'
              qmasse=.false.
          endif
      endif

      if (.not.qmasse) then
          write(6,'(2A)') program,' Masses are all set to one.'
          do i=1,natom
             massat(i)=1.d0
          enddo
      endif
 
c     Test/identification des residus.
 
      if (qpdb) then
      nunknown=0
      do i=1,nres
         ires=fatres(i)
         idres(i)=index(residus_standards,resnam(ires))/4
         if (idres(i).le.0) then
             nunknown=nunknown+1
             if (nunknown.lt.10) then
                 write(6,'(3A)') 
     .         "%Pdbmat-W> residue:'",
     .           resnam(ires),"' is not a well known amino-acid."
                 idres(i)=idmax
             else if (nunknown.eq.10) then
                 write(6,'(A)') 
     .         '%Pdbmat-W> ........'
                 idres(i)=idmax
             endif
         endif
      enddo
 
      if (nunknown.gt.0) then
          write(6,'(/A,I6,A)') 
     .  '%Pdbmat-W> ',nunknown,' residue(s) not known.'
      endif
      endif
 
c     Detection des liaisons covalentes:
 
      nbonds=0
      if (typbond.ne.' ALL'.and.typbond.ne.'CONSECUTIF') goto 200
 
      if (typbond.eq.' ALL') then
         nbmax=0
         nbmin=999
         imax=-1
         imin=-1
         k=1
         nbond(1)=1
         do i=1,natom
            do j=1,natom
               if (i.ne.j) then
                   rx=xat(i)-xat(j)
                   ry=yat(i)-yat(j)
                   rz=zat(i)-zat(j)
                   dist=dsqrt(rx*rx + ry*ry + rz*rz)
                   if (dist.le.cutbnd) then
                       jbond(k)=j
                       k=k+1
                   endif
               endif
            enddo
            nbond(i+1)=k
            if (nbond(i+1)-nbond(i).gt.nbmax) then
                nbmax=nbond(i+1)-nbond(i)
                imax=i
            endif
            if (nbond(i+1)-nbond(i).lt.nbmin) then
                nbmin=nbond(i+1)-nbond(i)
                imin=i
            endif
            if (k-1.gt.ntopmax) then
                write(6,'(/A,I12)') 
     .        '%Pdbmat-Err> Too many bonds. Maximum is: ',
     .          ntopmax
                stop
            endif
         enddo
         nbonds=k-1
      else if (typbond.eq.'CONSECUTIF') then
c
c        On fait attention aux distances...
c        Il peut y avoir plusieurs molecules,
c        plusieurs chaines, dans le systeme.
c       (A nettoyer: deux fois le meme calcul)
c
         nbond(1)=1
         jbond(1)=2
         k=2
         do i=2,natom-1
            nbond(i)=k
            j=i-1
            rx=xat(i)-xat(j)
            ry=yat(i)-yat(j)
            rz=zat(i)-zat(j)
            dist=dsqrt(rx*rx + ry*ry + rz*rz)
            if (dist.le.cutbnd) then
                jbond(k)=j
                k=k+1
            endif
            j=i+1
            rx=xat(i)-xat(j)
            ry=yat(i)-yat(j)
            rz=zat(i)-zat(j)
            dist=dsqrt(rx*rx + ry*ry + rz*rz)
            if (dist.le.cutbnd) then
                jbond(k)=j
                k=k+1
            endif
            if (k.gt.ntopmax) then
                write(6,'(/A,I12)') 
     .        '%Pdbmat-Err> Too many bonds. Maximum is: ',
     .          ntopmax
                stop
            endif
         enddo
         nbond(natom)=k
         jbond(k)=natom-1
         nbond(natom+1)=k+1
         imax=2
         imin=1
         nbmin=1
         nbmax=2
         nbonds=k
      endif
 
      if (nbonds.eq.0) then
          write(6,'(/A/)') 
     .  '%Pdbmat-W> No bond found.'
          goto 200
      endif
 
      if (prtlev.gt.0)
     .write(6,'(A,I6,A,F5.2,A)') program,
     .  nbonds/2,' covalent bonds, i.e.,',
     .  float(nbonds)/float(2*natom),' per atom.'

      if (prtlev.gt.1)
     .write(6,'(A,I3,A,I6)') 
     .'         Maximum number found =',nbmax,' for atom ',imax,
     .'         Minimum number found =',nbmin,' for atom ',imin
 
c     Angles de valence:
 
      namax=0
      namin=9999
      imax=-1
      imin=-1
      ii=1
      nangle(1)=1
      do i=1,natom
         if (nbond(i+1).gt.nbond(i)) then
         do jj=nbond(i),nbond(i+1)-1
            j=jbond(jj)
            if (nbond(j+1).gt.nbond(j)) then
            do kk=nbond(j),nbond(j+1)-1
               k=jbond(kk)
               if (k.ne.i) then
                   jangle(ii)=k
                   ii=ii+1
               endif
            enddo
            endif
         enddo
         endif
         nangle(i+1)=ii
         if (nangle(i+1)-nangle(i).gt.namax) then
             namax=nangle(i+1)-nangle(i)
             imax=i
         endif
         if (nangle(i+1)-nangle(i).lt.namin) then
             namin=nangle(i+1)-nangle(i)
             imin=i
         endif
         if (ii.gt.ntopmax) then
             write(6,'(/A,I12)') 
     .      '%Pdbmat-Err> Too many angles. Maximum is: ',
     .       ntopmax
             stop
         endif
      enddo
      nangles=ii-1
 
      if (prtlev.gt.0)
     .write(6,'(A,I6,A,F5.2,A)') program,
     .  nangles/2,' valence angles, i.e.,',
     .  float(nangles)/float(2*natom),' per atom.'

      if (prtlev.gt.1)
     .write(6,'(A,I3,A,I6)') 
     .'         Maximum number found =',namax,' for atom ',imax,
     .'         Minimum number found =',namin,' for atom ',imin

c______________________________________________________________
c
c     Matrice des derivees secondes:
c______________________________________________________________
c
 200  continue

c     Coordonnees et masses utilisees, sauvegardees:

      unout=nunit
      nunit=nunit+1
      nomfich="pdbmat.xyzm"
      cformat="FORMATTED"
      cstatus="ove"
      call openam(nomfich,cformat,cstatus,unout,.true.,
     .     qinter,qexist)

      do i=1,natom
         write(unout,'(4(1PG20.12),I9)')  
     .   xat(i), yat(i), zat(i), massat(i), iresat(i)
      enddo
      close(unout) 
      if (prtlev.gt.0)
     .write(6,'(2A)') program,
     .    ' Coordinates and masses considered are saved.'

      if (qbinary) then
      nomfich="pdbmat.sdijb"
      cformat="UNFORMATTED"
      else
      nomfich="pdbmat.sdijf"
      cformat="FORMATTED"
      endif
      cstatus="ove"
      call openam(nomfich,cformat,cstatus,unout,.true.,
     .     qinter,qexist)
 
c     ========================================
c     Les atomes sont tous lies deux a deux,
c     par un potentiel "universel" (M.Tirion).
c     ========================================
 
      elmax=0.d0
      trace=0.d0
      nnzero=0
      nntr=0
      ll=0
      nbig=0

      do i=1,natom
         ii=3*i-2
         nvoisat(i)=0
 
c        On calcule trois lignes de la matrice a la fois:
c        -----------------------------------------------
         do j=1,3*natom
            der2(1,j)=0.d0
            der2(2,j)=0.d0
            der2(3,j)=0.d0
         enddo
c
         do j=1,natom
            if (i.ne.j) then
            jj=3*j-2
            kij=knonb
 
            rx=xat(i)-xat(j)
            ry=yat(i)-yat(j)
            rz=zat(i)-zat(j)
            dist2=rx*rx + ry*ry + rz*rz
            dist=dsqrt(dist2)
 
            if (dist.lt.rsmall) then
                write(6,'(/2A,1PG10.4,A/2(I6,2A,I6,1X,2A))') 
     .          progrer,' Too small distance = ',dist,
     .        ' between following atoms.',
     .          i,': ',resnam(i),iresat(i),atonam(i),' and ',
     .          j,': ',resnam(j),iresat(j),atonam(j)
                stop '*Wrong coordinates*'
            endif
 
c           Constante de force topologique:
 
            if (nbonds.gt.0) then
            if (nbond(i+1).gt.nbond(i)) then
                do k=nbond(i),nbond(i+1)-1
                   if (jbond(k).eq.j) then
                       kij=kbond
                       goto 300
                   endif
                enddo
            else
            goto 300
            endif
c
            if (nangle(i+1).gt.nangle(i)) then
                do k=nangle(i),nangle(i+1)-1
                   if (jangle(k).eq.j) then
                       kij=kangle
                       goto 300
                   endif
                enddo
            endif
            endif
 300        continue
 
c           Calcul des elements: (potentiel harmonique)
c           -------------------------------------------
            if (dist.le.cutoff) then
                ll=ll+1
                nvoisat(i)=nvoisat(i)+1
 
c               Elements diagonaux des blocs i et j:
c               -----------------------------------
                ddf=kij/dist2
                elemnt=rx*rx*ddf
                der2(1,ii)=der2(1,ii)+elemnt
                der2(1,jj)=der2(1,jj)-elemnt
                elemnt=ry*ry*ddf
                der2(2,ii+1)=der2(2,ii+1)+elemnt
                der2(2,jj+1)=der2(2,jj+1)-elemnt
                elemnt=rz*rz*ddf
                der2(3,ii+2)=der2(3,ii+2)+elemnt
                der2(3,jj+2)=der2(3,jj+2)-elemnt
 
c               Elements extra-diagonaux des deux blocs:
c               ---------------------------------------
                elemnt=rx*ry*ddf
                der2(1,ii+1)=der2(1,ii+1)+elemnt
                der2(2,ii)=der2(2,ii)+elemnt
                der2(1,jj+1)=der2(1,jj+1)-elemnt
                der2(2,jj)=der2(2,jj)-elemnt
                elemnt=rx*rz*ddf
                der2(1,ii+2)=der2(1,ii+2)+elemnt
                der2(3,ii)=der2(3,ii)+elemnt
                der2(1,jj+2)=der2(1,jj+2)-elemnt
                der2(3,jj)=der2(3,jj)-elemnt
                elemnt=ry*rz*ddf
                der2(2,ii+2)=der2(2,ii+2)+elemnt
                der2(3,ii+1)=der2(3,ii+1)+elemnt
                der2(2,jj+2)=der2(2,jj+2)-elemnt
                der2(3,jj+1)=der2(3,jj+1)-elemnt
            endif
            endif
         enddo
 
c        Sortie de la matrice-bande calculee:
c        -----------------------------------
c       (Uniquement la demi-matrice superieure)
 
c        Level-shift, pour eviter les zeros numeriques,
c        lors de la diagonalisation a venir:
c       (la minimisation est parfaite, par definition)
 
         der2(1,ii)=der2(1,ii)+levelshft
         der2(2,ii+1)=der2(2,ii+1)+levelshft
         der2(3,ii+2)=der2(3,ii+2)+levelshft
 
         do j=ii,3*natom
            jat=(j-1)/3+1
            if (der2(1,j).ne.0.d0) then
                nnzero=nnzero+1
                if (qbinary) then
                write(unout) 
     .          ii,j,der2(1,j)/dsqrt(massat(i)*massat(jat))
                else
                write(unout,'(2I6,1PG20.12)') 
     .          ii,j,der2(1,j)/dsqrt(massat(i)*massat(jat))
                endif
                if (dabs(der2(1,j)).gt.rbig) nbig=nbig+1
                if (dabs(der2(1,j)).gt.elmax) elmax=dabs(der2(1,j))
            endif
         enddo
         do j=ii+1,3*natom
            jat=(j-1)/3+1
            if (der2(2,j).ne.0.d0) then
                nnzero=nnzero+1
                if (qbinary) then
                write(unout) 
     .          ii+1,j,der2(2,j)/dsqrt(massat(i)*massat(jat))
                else
                write(unout,'(2I6,1PG20.12)') 
     .          ii+1,j,der2(2,j)/dsqrt(massat(i)*massat(jat))
                endif
                if (dabs(der2(2,j)).gt.rbig) nbig=nbig+1
                if (dabs(der2(2,j)).gt.elmax) elmax=dabs(der2(2,j))
            endif
         enddo
         do j=ii+2,3*natom
            jat=(j-1)/3+1
            if (der2(3,j).ne.0.d0) then
                nnzero=nnzero+1
                if (qbinary) then
                write(unout) 
     .          ii+2,j,der2(3,j)/dsqrt(massat(i)*massat(jat))
                else
                write(unout,'(2I6,1PG20.12)') 
     .          ii+2,j,der2(3,j)/dsqrt(massat(i)*massat(jat))
                endif
                if (dabs(der2(3,j)).gt.rbig) nbig=nbig+1
                if (dabs(der2(3,j)).gt.elmax) elmax=dabs(der2(3,j))
            endif
         enddo
         elemnt=(der2(1,ii)+der2(2,ii+1)+der2(3,ii+2))/massat(i)
         if (elemnt.eq.0.d0) then
             write(6,'(A,I6,A)')
     .     '%Pdbmat-W> Atom ',i,
     .     ' has a null second derivatives...'
         else
             nntr=nntr+1
         endif
         trace=trace+elemnt
      enddo
 
      nmoy=0.d0
      nrms=0.d0
      nmin=natom
      nmax=0
      do i=1,natom
         if (nvoisat(i).gt.nmax) nmax=nvoisat(i)
         if (nvoisat(i).lt.nmin) nmin=nvoisat(i)
         nmoy=nmoy+nvoisat(i)
         nrms=nrms+nvoisat(i)**2.d0
      enddo
      nmoy=nmoy/float(natom)
      nrms=nrms/float(natom)-nmoy**2.d0
      if (nrms.gt.0.d0) nrms=dsqrt(nrms)
 
      if (ll.eq.0) then
          write(6,'(/2A,I12,A)') progrer,
     .  ' No atom-atom interaction found. Too short cutoff ?'
          stop '*Empty matrix*'
      endif

      if (prtlev.gt.0)
     .write(6,'(/A)') ' Pdbmat> Matrix statistics:'

      write(6,'(/2A,F8.4,A)') program,' The matrix is ',
     .  100.d0*dfloat(nnzero)/dfloat(3*natom*(3*natom+1)/2),' % Filled.'
      write(6,'(A,I12,A)') 
     .' Pdbmat> ',nnzero,'  non-zero elements.'

      if (prtlev.gt.0) then
      write(6,'(A,I12,A)') 
     .' Pdbmat> ',ll/2,' atom-atom interactions.'
      write(6,'(/A,F9.2,A,F9.2/(A,I6))') 
     .' Pdbmat> Number per atom= ',nmoy,' +/- ',nrms,
     .'         Maximum number = ',nmax,
     .'         Minimum number = ',nmin
      endif
 
      write(6,'(/A,1PG12.6)') 
     .' Pdbmat> Matrix trace   = ',trace

      if (prtlev.gt.0) then
      write(6,'(A,1PG12.6)') 
     .' Pdbmat> Larger element = ',elmax
      write(6,'(A,I6,A,1PG8.1)') 
     .' Pdbmat> ',nbig,' elements larger than +/- ',rbig
      endif
 
      write(6,'(/A)') ' Pdbmat> Normal end.'
      stop
      end
c-----------------------------------------------------------------------------
      subroutine stringcl(chaine,nonblancs)
c
c     Les caracteres "blancs" de la CHAINE sont retires (a gauche et au milieu).
c     L'entier NONBLANCS donne la position du dernier caractere.
c
c     YHS-Jan-95, Toulouse.
c I/O:
      integer nonblancs
      character*(*) chaine
c Local:
      integer icar, ncar, taille
c Begin:
      taille=len(chaine)
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
c---------------------------------------------------------------------
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
c--------------------------------------------------------------------
      SUBROUTINE openam(namfil,cformat,cstatus,unit,qverbos,
     .                  qinterr,qexist)
c
c     Ouverture d'un fichier de nom NAMFIL, sur l'unite UNIT,
c     a priori suite a une interrogation...
c
c     input:
c        namfil: nom du fichier a ouvrir. 
c        "stop", "end", "fin", "quit" : arretent le programme.
c        cstatus: mots-cles fortran... ou "OVE" pour overwrite.
c     output: 
c        qexist: flag / existence du fichier 
c        qinterr: Pas de nom pour le fichier cherche.
c I/O:
      logical qinterr, qverbos, qexist
      integer unit
      character*(*) namfil, cformat, cstatus
c Local
      integer lnom
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
         write(6,'(2A)') 'Openam> Program is stopping on user request.'
         stop                                                                   
      endif 
 
c     Checks if filename is consistent with the opening:

      call stringcl(namfil,lnom)

      inquire(file=namfil,exist=qexist)
      if (.not.qexist.and.cstatus.eq.'OLD') then
          qinterr=.true.
          if (qverbos) then
              write(6,'(/2A)') '%Openam-Err> File: ',namfil(1:lnom)
              write(6,'(A)') 
     .      ' Expected in the current directory, but not found.'
          endif
          return
      endif
c
      if (qexist.and.cstatus.eq.'NEW') then
         write(6,'(/2A)') 
     .      '%Openam-Err> This file exists:',namfil(1:lnom)
         stop
      else if (qexist.and.cstatus.eq.'OVE') then
         ordrunix='rm '//namfil
         call system(ordrunix)
      endif
      if (cstatus.eq.'OVE') cstatus='NEW'
c                                                                   
      open(file=namfil,form=cformat,
     .     status=cstatus,unit=unit)                
c        
      if (qverbos) then
         write(6,'(/2A)') ' Openam> File opened: ',namfil(1:lnom)
      endif

      return                                                                       
      end
c--------------------------------------------------------------------
      subroutine rdatompdb(unpdb,ssusel,xat,yat,zat,binfo,
     .           atonam,iresat,resnam,natmax,natom,
     .           fatres,nresmx,nres,qerror,prtlev)
c
c     Lecture ligne a ligne d'un fichier pdb.
c     Uniquement les lignes commencant par 'ATOM'.
c     Uniquement ceux de la sous-unite selectionnee.
c
c     fatres(i): numero du premier atome du residu i.
c     YHS-novembre-1996.
c
      implicit none
cI/O:
      integer unpdb, natmax, iresat(*), natom, 
     .        nresmx, nres, fatres(*), prtlev
      double precision xat(*), yat(*), zat(*), binfo(*)
      logical qerror
      character*4 atonam(*), resnam(*)
      character*1 ssusel
cLocal:
      integer nerr, iatom, irs, irsprev, lmot, ntit
      double precision x, y, z, bfact
      character*1  ssu
      character*4  ren
      character*5  atncur
      character*80 lign80
cBegin:
      if (prtlev.gt.0)
     .write(6,'(/A)') ' Rdatompdb> Reading pdb file.'
c
      qerror=.false.
      nerr=0
c
      irsprev=-1
      nres=0
      iatom=0
      ntit=0
 105  continue   
      read(unpdb,'(A)',end=200,err=110) lign80 
c 
      goto 120                                
 110  continue
      nerr=nerr+1                            
c
 120  continue                              
      if (lign80(1:4).eq.'ATOM'.or.lign80(1:6).eq.'HETATM') then
      read(lign80,'(12X,A4,1X,A4,A1,I4,4X,3F8.3,6X,F6.2)') 
     .            atncur, ren, ssu, irs, x, y, z, 
     .            bfact
      if (iatom.lt.natmax) then
          if (ssu.eq.ssusel.or.ssusel.eq.' ') then
          iatom=iatom+1
          xat(iatom)=x
          yat(iatom)=y
          zat(iatom)=z
          binfo(iatom)=bfact
c
          call stringcl(atncur,lmot)
          atonam(iatom)=atncur
          call stringcl(ren,lmot)
          resnam(iatom)=ren
          iresat(iatom)=irs
c
          if (irs.ne.irsprev) then
              nres=nres+1
              if (nres.gt.nresmx) then
                  write(6,'(A/A,I6)') 
     .          '%Rdatompdb-Err> Too many residues in this file.',
     .          ' Maximum allowed is = ',nresmx
                  stop
              endif
              irsprev=irs
              fatres(nres)=iatom
          endif
          endif
      else
          write(6,'(A/A,I6)') 
     .      '%Rdatompdb-Err> Too many atoms in this file.',
     .      ' Maximum allowed is = ',natmax
          stop
      endif
      else if (lign80(1:6).eq.'REMARK'.and.prtlev.gt.0) then
          ntit=ntit+1
          if (ntit.le.10) then
              write(6,'(A)') lign80
          else if (ntit.eq.11) then
              write(6,'(A)') ' .... '
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
      if (prtlev.gt.1) then
      write(6,'(A)') ' Rdatompdb> End of file reached.'
      write(6,'(A,I6)') ' Rdatompdb> Number of I/O errors: ',nerr
      endif
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
      if (prtlev.gt.0)
     .write(6,'(A,F8.1)') 
     .'            Mean number per residue  = ',float(natom)/float(nres)
c
      if (natom.eq.0) then
          write(6,'(A)')
     .  '%Rdatompdb-Err> No atom found in file.'
          qerror=.true.
      endif
      if (nres.eq.0) then
          write(6,'(A)')
     .  '%Rdatompdb-Err> No residue found in file.'
          qerror=.true.
      endif
c
      return
      end
c-----------------------------------------------------------------------------
      subroutine mintomaj(chaine)
 
c     Les caracteres minuscules sont mis en MAJUSCULES.
c     Les autres ne sont pas touches.
 
c     YHS-Oct-98: Premiere version (Toulouse).
c     YHS-Sep-03: Dernieres modifications (Lyon).
 
      character*(*) chaine
c Local:
      integer icar, ilettre, taille
      character*26  carmaj, carmin
 
      carmin='qwertyuiopasdfghjklzxcvbnm'
      carmaj='QWERTYUIOPASDFGHJKLZXCVBNM'
 
      taille=len(chaine)
      if (taille.le.0) return

      do icar=1,taille
         ilettre=index(carmin,chaine(icar:icar))
         if (ilettre.gt.0) then
             chaine(icar:icar)=carmaj(ilettre:ilettre)
         endif
      enddo
 
      return
      end
c-----------------------------------------------------------------------------
      subroutine readxyz(uninp,x,y,z,w,ic,nmax,ncoor,ndat,qerror,prtlev)

c     Reads at most NMAX coordinates in free format. 

c     Either:
c     x, y, z
c     or:
c     x, y, z, w
c     or:
c     x, y, z, w, ic

c     If first word in ligne is not a number, the whole ligne is
c     assumed to be a title or a commentary.

c     YHS-Sep-03: First version (Lyon).

cI/O:
      logical qerror
      integer ic(*), ncoor, ndat, nmax, prtlev, uninp
      double precision w(*), x(*), xc, y(*), yc, z(*), zc
cLocal:
      integer nmotsmax
      parameter(nmotsmax=255)

      integer i, lmot, nchi, nlmax, nl, nlu, nmots, stats(0:nmotsmax)
      double precision rlu
      character chiffres*15, lignlg*(nmotsmax), 
     .        mots(nmotsmax)*(nmotsmax), program*9, progrer*11, 
     .        progrwn*11

cBegin:
      program=' Readxyz>'
      progrer='%Readxyz-Er>'
      progrwn='%Readxyz-Wn>'

      chiffres='1234567890.eE+-'

      do i=1,nmotsmax
         stats(i)=0
      enddo

c     Lecture ligne a ligne:
c     ----------------------

      if (prtlev.gt.1) write(6,'(/2A)') program,
     .  ' Comments, or lignes with less than three numbers: '

      qerror=.false.
      ncoor=0
      nl=0
 100  continue
      read(uninp,'(A)',end=200) lignlg 
      nl=nl+1
      call string_split(lignlg,nmotsmax," ",mots,nmotsmax,nmots)

      nfound=0
      do i=1,nmots
         call stringcl(mots(i),lmot)
         if (lmot.le.0) goto 150

c        Commentaire ?
         if (mots(i)(1:1).eq.'!') goto 150

c        Chiffre ?

         do k=1,lmot
            if (index(chiffres,mots(i)(k:k)).le.0) goto 110
         enddo

         nfound=nfound+1

         if (nfound.le.4) then
             read(mots(i)(1:lmot),*,err=110) rlu
             if (nfound.eq.1) xc=rlu
             if (nfound.eq.2) yc=rlu
             if (nfound.eq.3) zc=rlu
             if (nfound.eq.4) wc=rlu
         else if (nfound.eq.5) then 
             read(mots(i)(1:lmot),*,err=110) nlu
         endif

c        Mot suivant:
 110     continue

c        Le premier mot n'est pas un chiffre => ligne de commentaires

         if (nfound.eq.0) goto 150
      enddo
 150  continue

c     Stockage des coordonnees:
c     -------------------------

      stats(nfound)=stats(nfound)+1

      if (nfound.ge.3) then
          ncoor=ncoor+1
          if (ncoor.le.nmax) then
              x(ncoor)=xc
              y(ncoor)=yc
              z(ncoor)=zc
              if (nfound.eq.4) then
                  w(ncoor)=wc
              else
                  w(ncoor)=wc
                  ic(ncoor)=nlu 
              endif
          else
              write(6,'(/2A,I9,A)') progrer,' More than ',
     .        nmax,' particles in file.'
              write(6,'(2A)') progrer,
     .      ' Please increase program memory limits (Sorry for that).'
              stop
          endif
      else
          if (prtlev.gt.1) then
              write(6,'(2A)') lignlg(1:72),'...'
          endif
      endif
      
c     Ligne suivante:
      goto 100
 200  continue
      if (prtlev.gt.1.and.nl.eq.ncoor) write(6,'(2A)') program,' None.'

      write(6,'(/2A,I7)') program,
     .' Number of particles in file (with x,y,z coordinates): ',ncoor
 
      if (ncoor.eq.0) then
          write(6,'(/2A)') progrer,' No coordinate found in file.'
          qerror=.true.
      endif

      nchi=0
      ndat=0
      nlmax=0
      do i=1,nmotsmax
         if (stats(i).gt.0) nchi=nchi+1
         if (stats(i).gt.nlmax) then
             nlmax=stats(i)
             ndat=i
         endif
      enddo

      do i=0,nmotsmax
         if (stats(i).gt.0.and.(prtlev.gt.1.or.nchi.gt.1)) then
            write(6,'(A,I6,A,I7,A)') program,i,
     .    ' numbers found in ',stats(i),' lignes.'  
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------------
      subroutine vecstat(vect,nmax,rmin,rmax,rave,rdev)

c     Statistics for vector Vect(NMAX):
c     minimum, maximum, average (rave) and standard deviation (rdev).

c     YHS-Sep-03: First version (Lyon).

cI/O:
      integer nmax
      double precision rave, rdev, rmax, rmin, vect(*)
cLocal:
      integer i
      character program*9, progrer*11, progrwn*11

cBegin:
      program=' Vecstat>'
      progrer='%Vecstat-Er>'
      progrwn='%Vecstat-Wn>'

      rave=0.d0
      rdev=0.d0
      rmin=-9999.d0
      rmax=9999.d0

      if (nmax.le.0) then
          write(6,'(2A)') progrer,' Zero-length vector.'
          return
      endif

      do i=1,nmax
         if (vect(i).gt.rmax.or.i.eq.1) rmax=vect(i)
         if (vect(i).lt.rmin.or.i.eq.1) rmin=vect(i)
         rave=rave+vect(i)
         rdev=rdev+vect(i)**2.0
      enddo

      rave=rave/dfloat(nmax) 
      rdev=rdev/dfloat(nmax)-rave*rave
      if (rdev.gt.0.d0) rdev=dsqrt(rdev)

      return
      end
