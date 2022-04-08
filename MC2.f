      program montecarlo2

c     LLUÍS TORRES HUGAS
c     darrera modificació: 10/01/2022

c     SIMULACIÓ MODEL D'ISING 2D UTILITZANT EL MÈTODE MONTE-CARLO
c     Aquest codi permet calcular una estimació estadística de 
c     l'energia, la magnetització, la capacitat calorífica i
c     la susceptibilitat en funció de la temperatura per una
c     xarxa 2d LxL amb spins uniaxials.
      
      implicit none
c-----------------------------------------------------------------------
c     DEFINICIÓ DE TOTES LES VARAIBLES
      integer i, j, k
      integer L, N, MCTOT, MCINI, MCD, NSEED, SEED, SEED0, NDAT 
      integer DE, IMC, IPAS

c     valor del nombre de temperatures 
      parameter (NDAT = 50) 

      real*8 TSTEP, T_vect(0:NDAT),E_vect(0:NDAT),dE_vect(0:NDAT)
      real*8 TEMP, W(-8:8), ENERG, ENE, MAGNE, Cv, X

      real*8 genrand_real2

      logical valid

      real*8 SUM, SUME, SUME2, SUMM, SUMM2, SUMAM, VARE, VARM

c     paràmetres que determinen la simulació
      parameter (L = 64)
      parameter (N = L*L)
      parameter (MCTOT = 30000) 
      parameter (MCINI = 2000) 
      parameter (MCD = 20) 
      parameter (NSEED = 200) 
      parameter (SEED0 = 69420)

      integer*2 S(1:L,1:L)
      integer*4 PBC(0:L+1)

      real*4 TIME1, TIME2

      character*48 NOM

      real*8 maxCv, maxX, TmaxCv, TmaxX
      
      external genrand_real2, dfun

c     s'assigna el nom del fitxer de dades en funció dels paràmetres de simulació
      NOM = "SIM-MC-L-64-MCTOT-30000-NSEED-200-T(200-250)-050"

c     guardem temps inicial de l'ordinador per posteriorment conèixer el temps de càlcul
      call cpu_time(TIME1)

c-----------------------------------------------------------------------
c     CONDICIONS DE CONTORN + VECTOR TEMPERATURES
      PBC(0)=L
      PBC(L+1)=1
      do i=1,L
            PBC(i)=i
      end do

      T_vect(0) = 2.0d0
      TSTEP = (2.5d0 - 2.0d0)/NDAT
      do i = 1,NDAT
            T_vect(i) = T_vect(i-1) + TSTEP
      end do
    
c-----------------------------------------------------------------------
c     SIMULACIÓ MONTE-CARLO

c     preparem el fitxer de dades
      open(unit = 1,file = NOM//".dat")
100   FORMAT(8(A20,8X))
      write(1,100) 'T','<e>','σ_e','<m>','<|m|>','σ_m','Cv','X'
            
      DO k = 0,NDAT

      TEMP = T_vect(k)

c     valors de e^{-βΔH} necessaris per comprovar si s'accepta o no el canvi
      do DE = -8,8
            W(DE) = exp(-DE/TEMP)
      end do

c     inicialització de les variables per fer l'estadística
      SUM=0.0d0
      SUME=0.0d0
      SUME2=0.0d0
      SUMM=0.0d0
      SUMM2=0.0d0
      SUMAM=0.0d0

      do SEED = SEED0,SEED0+NSEED-1,1
            
            call init_genrand(SEED)

           !GENERACIÓ DE LA MATRIU S(i,j)
            do i = 1,L
                  do j = 1,L
                        if (genrand_real2() < 0.5d0) then
                              S(i,j) = 1
                        else
                              S(i,j) = -1
                        end if
                  end do
            end do
            !càlcul de l'energia
            ENE = ENERG(S,L,PBC)

            do IMC = 1, MCTOT
                  do IPAS = 1, N
                        call metropolis(S,L,DE,valid,PBC,W)
                        if (valid) then
                              ENE = ENE + DE
                        end if
                  end do
                  

                  if ((IMC.gt.MCINI).and.(MCD*(IMC/MCD).eq.IMC)) then
                        SUM = SUM + 1.0d0
                        SUME = SUME + ENE
                        SUME2 = SUME2 + ENE**2
                        SUMM = SUMM + MAGNE(S,L)
                        SUMM2 = SUMM2 + MAGNE(S,L)**2
                        SUMAM = SUMAM + abs(MAGNE(S,L))
                  end if

            end do
      end do
c-----------------------------------------------------------------------
c     CALCUL DE LES VARIABLES TERMODINÀMIQUES

      SUME  = SUME/SUM !<E>
      SUME2 = SUME2/SUM !<E^2>
      VARE  = (SUME2 - SUME*SUME)/SUM !Var(E)

      SUMM  = SUMM/SUM !<M>
      SUMAM = SUMAM/SUM !<|M|>
      SUMM2 = SUMM2/SUM !<M^2>=<|M|^2>
      VARM  = (SUMM2 - SUMM*SUMM)/SUM !Var(M)

c     es guarden els valors de l'energia per posteriorment fer la derivada
      E_vect(k) = SUME/dble(N)

c     capacitat calorífica
      Cv = (SUME2 - SUME**2)/(TEMP**2)/dble(N)

c     susceptibilitat
      X = (SUMM2 - SUMAM**2)/TEMP/dble(N)

c     es busca el màxim de cada magnitud termodinàmica i es guarden les corresponents T
      if (Cv > maxCv) then
            maxCv = Cv
            TmaxCv = TEMP
      end if

      if (X > maxX) then
            maxX = X
            TmaxX = TEMP
      end if

200   FORMAT(8(F20.10,8X))

c     anotació dels valors obtinguts 
      write(1,200) TEMP, SUME/N, sqrt(VARE)/N,              
     .                    SUMM/N, SUMAM/N, sqrt(VARM)/N,
     .                    Cv, X
      END DO

      write(1,*)
      write(1,*)
300   FORMAT(8(A20,8X))
      write(1,300) "T","d<e>/dT"
c     càlcul de la derivada per comparar amb C_v
      call dfun(NDAT, T_vect, E_vect, dE_vect)
400   FORMAT(2(F20.10,8X))
      do i = 1, NDAT
            write(1,400) T_vect(i), dE_vect(i)
      end do
c     guardem temps inicial de l'ordinador per coneixer el temps de càlcul
      call cpu_time(TIME2)
      
c     escrivim al document el temps total de càlcul i les Tmax de C_v i X
      write(1,*)
      write(1,*) "CPU TIME [s]: ", TIME2-TIME1
      write(1,*)
      write(1,*)
      write(1,*) "Temperatura Crítica de la Cv: ", TmaxCv
      write(1,*) "Temperatura Crítica de la X: ", TmaxX
      close(1) 


      end program montecarlo2

c-----------------------------------------------------------------------
c     FUNCIONS
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     funció que calcula la magnetització d'una configuració S(i,f)
      real*8 function MAGNE(S,L)
      integer*2 S(1:L,1:L)
      integer*4 i,j,L
      real*8 MAG

      MAG=0.0d0
      
      do i = 1,L
            do j = 1,L
                  MAG = MAG + S(i,j)
            end do
      end do
      
      MAGNE = MAG
      
      return
      end

c-----------------------------------------------------------------------
c     funció que calcula l'energia d'una configuració S(i,j)
      real*8 function ENERG(S,L,PBC)
      integer*2 S(1:L,1:L)
      integer*4 i,j,L
      integer*4 PBC(0:L+1)
      real*8 ENE

      ENE = 0.0d0
      
      do i = 1,L
            do j = 1,L
                  ENE = ENE-S(i,j)*S(PBC(i+1),j)-S(i,j)*S(i,PBC(j+1))
            end do
      end do

      ENERG = ENE
      
      return
      end


c-----------------------------------------------------------------------
c     SUBRUTINES
c-----------------------------------------------------------------------
c     subrutina que càlcula la derivada numèrica a partir de les llistes
c     dels valors de f(x) a certs punts x
      subroutine dfun(ndat,x,f,df)
      implicit none
      real*8 h, x(ndat), f(ndat), df(ndat)
      integer ndat, i
      
      do i = 1,ndat  
            h = abs(x(2) - x(1)) 

            if (i == ndat) then
                  df(i) =  (f(i) - f(i-1))/h
            else if (i == 1) then
                  df(i) = (f(i+1)-f(i))/h
            else
                  df(i) = (f(i+1)-f(i-1))/(2*h)
            end if
      end do
      return
      end



c-----------------------------------------------------------------------
c     subrutina que permet guardar la configuració d'un S(i,j)
      subroutine WRITECONFIG(S,L)
      implicit none
      integer i, j, L
      integer*2 S(1:L,1:L)
      open(2, file = 'P5-configuration.conf')
      do i = 1,L
            do j = 1,L
                  if (S(i,j) == 1) then
                        write(2,*) i,j
                  end if
            end do
      end do
      close(2)
      return
      end

c-----------------------------------------------------------------------
c     funció metropolis, encarregada de determinar si s'accepta o no el
c     canvi d'un S(i,j) determinat, a partir de l'algoridme explicat a 
c     l'informe
      subroutine metropolis(S,L,DE,valid,PBC,W)
      implicit none
      integer i, j, L, DE, PBC(0:L+1)
      integer*2 S(1:L,1:L)
      real*8 genrand_real2, W(-8:8), suma, DELTA
      logical valid

      external genrand_real2

      valid = .false.

      i = int(genrand_real2()*L) + 1
      j = int(genrand_real2()*L) + 1

c     càlcul de ΔH
      suma = 0
      suma = S(i,PBC(j+1))+S(i,PBC(j-1))+S(PBC(i+1),j)+S(PBC(i-1),j)
      DE = int(2*S(i,j)*suma)

c     es mira si s'accepta o no el canvi
      if (DE.lt.0) then
            S(i,j) = -S(i,j)
            valid = .true.

      else if (DE.gt.0) then

            DELTA = genrand_real2()
            
            if (DELTA.lt.W(DE)) then
                S(i,j) = -S(i,j)
                valid = .true.
            endif

      end if

      return
      end