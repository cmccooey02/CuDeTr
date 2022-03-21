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
Program sample

  Use fdf
#ifdef MPI
  use mpi
#endif

  Implicit None

  Integer  maxa
  Parameter ( maxa = 100 )

  Character         fname*20, symbol(maxa)*2
  Integer           i, ia, isa(maxa), na, na_default, iblk
  Real              wmix

  Double Precision  factor, xa(3,maxa), cutoff, phonon_energy
  character*132 line
  Logical doit, debug
  type(block), pointer :: bp
#ifdef MPI
  integer MPIerror, rc
#endif
  integer Node, Nodes

#ifdef MPI
      call MPI_Init( MPIerror)
      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
#else
      Node = 0
      Nodes = 1
#endif

  Call fdf_init('sample.fdf','sample.out')

  If (fdf_defined('new-style')) Write(6,*) 'New-style stuff'

  na_default = 0
  na = fdf_get('NumberOfAtoms', na_default )
  Write(6,*) 'examples: na =', na, Node

  fname = fdf_get('NameOfFile','calimero')
  Write(6,*) fname, Node

  cutoff = fdf_physical('MeshCutoff',8.d0,'Ry')
  Write(6,*) cutoff, Node

  phonon_energy = fdf_get('phonon-energy',0.01d0,'eV')
  Write(6,*) phonon_energy, Node

  i = fdf_integer('SomeInt',34)
  Write(6,*) i, Node

  wmix = fdf_single('WmixValue',0.55)
  Write(6,*) wmix, Node

  factor = fdf_double('FactorValue',1.d-10)
  Write(6,*) factor, Node

  debug = fdf_boolean('Debug',.True.)
  Write(6,*) debug, Node

  doit = fdf_boolean('DoIt',.False.)
  Write(6,*) doit, Node

  If (.not.fdf_block('DFDFFFDDGF',bp)) Then
     write(6,*) 'Node ', Node, " did not find block"
  Endif
!!  call destroy(bp)

  If (fdf_block('AtomicCoordinatesAndAtomicSpecies',bp)) Then
     ia = 0
     loop: do
        if (.not. fdf_bline(bp,line)) exit loop
        ia = ia + 1
        Read(line,*) (xa(i,ia),i=1,3), isa(ia)
     Enddo loop
     write(6,*) 'Node, final ia:', Node, ia
  Endif
!!  call destroy(bp)

  If (fdf_block('AtomicSymbolsAndAtomicCoordinates',iblk)) Then
     Do ia = 1,na
        Read(iblk,*) symbol(ia), (xa(i,ia),i=1,3)
     Enddo
  Endif

  Do ia = 1,na
     Write(6,*) (xa(i,ia),i=1,3)
  Enddo

#ifdef MPI
  call MPI_Finalize(rc)
#endif

End Program sample


