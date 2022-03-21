! 
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2003
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
      program sample
c
c     Shows fdf capabilities..
c

      implicit none
      integer  maxa
      parameter ( maxa = 100 )

      character         fname*20, symbol(maxa)*2
      integer           i, ia, isa(maxa), na, na_default, iblk
      real              wmix
      double precision  factor, xa(3,maxa), cutoff, phonon_energy

      logical doit, debug

      include 'fdfdefs.h'
c
      call fdf_init('sample.fdf','sample.out')

      if (fdf_defined('new-style')) write(6,*) 'New-style stuff'

      na_default = 0
      na = fdf_integer('NumberOfAtoms', na_default )
      write(6,*) 'examples: na =', na

      fname = fdf_string('NameOfFile','calimero')
      write(6,*) fname

      cutoff = fdf_physical('MeshCutoff',8.d0,'Ry')
      write(6,*) cutoff

      phonon_energy = fdf_physical('phonon-energy',0.01d0,'eV')
      write(6,*) phonon_energy

      i = fdf_integer('SomeInt',34)
      write(6,*) i

      wmix = fdf_single('WmixValue',0.55)
      write(6,*) wmix

      factor = fdf_double('FactorValue',1.d-10)
      write(6,*) factor

      debug = fdf_boolean('Debug',.true.)
      write(6,*) debug

      doit = fdf_boolean('DoIt',.false.)
      write(6,*) doit
      
      if (fdf_block('AtomicCoordinatesAndAtomicSpecies',iblk)) then
        do ia = 1,na
          read(iblk,*) (xa(i,ia),i=1,3), isa(ia)
        enddo
      endif
      
      if (fdf_block('AtomicSymbolsAndAtomicCoordinates',iblk)) then
         do ia = 1,na
          read(iblk,*) symbol(ia), (xa(i,ia),i=1,3)
        enddo
      endif

      do ia = 1,na
         write(6,*) (xa(i,ia),i=1,3)
      enddo

      if (fdf_block('AtomicInfo',iblk)) then
        do ia = 1,na
          read(iblk,*) (xa(i,ia),i=1,3)
        enddo
      endif

      do ia = 1,na
         write(6,*) (xa(i,ia),i=1,3)
      enddo

      end






