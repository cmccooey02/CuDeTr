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
c
c Copyright Alberto Garcia, Jose Soler, 1996, 1997, 1998
c---
c     I/O variables for the fdf package. 
c
c     In Fortran 90, all this should go in a module...
c
c
c     ndepth: number of open files (maximum maxdepth)
c     fdf_stack holds their unit numbers.

      integer maxdepth
      parameter (maxdepth=5)
      integer ndepth, fdf_stack(maxdepth)
c
c     Unit numbers for input, output, error notification, and
c     debugging output (the latter active if fdf_debug is true)
c
      integer fdf_in, fdf_out, fdf_err, fdf_log
      common /fdf_io/ fdf_in, fdf_out, fdf_err, fdf_log, 
     $                ndepth, fdf_stack
      logical fdf_debug, fdf_debug2, fdf_started, fdf_donothing
      common /fdf_logicals/ fdf_debug, fdf_debug2, fdf_started,
     $                      fdf_donothing

      save /fdf_io/, /fdf_logicals/
c
c     Line just read and parsing info
c
      character*132 line
      integer maxntokens
      parameter (maxntokens=50)
      integer ntokens
      integer first(maxntokens), last(maxntokens)
      common /fdf_line/ line
      common /fdf_parsing/ ntokens, first, last

      save /fdf_line/, /fdf_parsing/
c---


