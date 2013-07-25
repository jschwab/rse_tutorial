! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def

      implicit none

      ! these routines are called by the standard run_star check_model
      contains

      subroutine extras_controls(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


      end subroutine extras_controls


      integer function extras_startup(s, id, restart, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         ierr = 0
         extras_startup = 0
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if
      end function extras_startup


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         extras_check_model = keep_going
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if
      end function extras_check_model


      integer function how_many_extra_history_columns(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         how_many_extra_history_columns = 2
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(s, id, id_extra, n, names, vals, ierr)

         ! MESA provides routines for taking logarithms which will handle
         ! cases where you pass them zero, etc.  Take advantage of this!
         use num_lib, only : safe_log10

         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr

         real(dp), parameter :: frac = 0.90
         integer :: i
         real(dp) :: edot, edot_partial

         ! calculate the total nuclear energy release rate by integrating
         ! the specific rate (eps_nuc) over the star.  using the dot_product
         ! intrinsic is a common idiom for calculating integrated quantities.
         ! note that one needs to explicitly limit the range of the arrays.
         ! never assume that the array size is the same as the number of zones.
         edot = dot_product(s% dm(1:s% nz), s% eps_nuc(1:s% nz))

         ! the center of the star is at i = s% nz and the surface at i = 1 .
         ! so go from the center outward until 90% of the integrated eps_nuc
         ! is enclosed.  exit and then i will contain the desired cell index.
         edot_partial = 0
         do i = s% nz, 1, -1
            edot_partial = edot_partial + s% dm(i) * s% eps_nuc(i)
            if (edot_partial .ge. (frac * edot)) exit
         end do

         ! note: do NOT add these names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.

         ! column 1
         names(1) = "m90"
         vals(1) = s% q(i) * s% star_mass  ! in solar masses

         ! column 2
         names(2) = "log_R90"
         vals(2) = safe_log10(s% R(i) / rsol) ! in solar radii

         ierr = 0
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         how_many_extra_profile_columns = 1
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(s, id, id_extra, n, nz, names, vals, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: k
         ierr = 0

         ! note: do NOT add these names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         if (n /= 1) stop 'data_for_extra_profile_columns'
         names(1) = 'beta'
         do k = 1, nz
            vals(k,1) = s% Pgas(k)/s% P(k)
         end do

      end subroutine data_for_extra_profile_columns


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         integer :: ierr
         integer :: f
         extras_finish_step = keep_going
         call store_extra_info(s)

         ! MESA provides a number of variables that make it easy to get user input.
         ! these are part of the star_info structure and are named
         ! x_character_ctrl, x_integer_ctrl, x_logical_ctrl, and x_ctrl.
         ! by default there are num_x_ctrls, which defaults to 100 of each.
         ! they can be specified in the controls section of your inlist.

         f = s% x_integer_ctrl(1)

         ! MESA also provides a number variables that are useful for implementing
         ! algorithms which require a state. if you just use these variables
         ! restarts, retries, and backups will work without doing anything special.
         ! they are named xtra1 .. xtra30, ixtra1 .. ixtra30, and lxtra1 .. lxtra30.
         ! they are automatically versioned, that is if you set s% xtra1, then
         ! s% xtra1_old will contains the value of s% xtra1 from the previous step
         ! and s% xtra1_older contains the one from two steps ago.

         s% xtra1 = s% log_center_density

         ! this expression will evaluate to true if f times the log center density
         ! has crossed an integer during the last step.  If f = 5, then we will get
         ! output at log center density = {... 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 ... }
         if ((floor(f * s% xtra1_old) - floor(f * s% xtra1) .ne. 0)) then

            ! save a profile & update the history
            s% need_to_update_history_now = .true.
            s% need_to_save_profiles_now = .true.

            ! by default the priority is 1; you can change that if you'd like
            ! s% save_profiles_model_priority = ?

         endif

      end function extras_finish_step


      subroutine extras_after_evolve(s, id, id_extra, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine extras_after_evolve


      ! routines for saving and restoring extra data so can do restarts

         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3


      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info


      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info


      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info


      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op

         integer :: i, j, num_ints, num_dbls, ierr

         i = 0
         ! call move_int or move_flg
         num_ints = i

         i = 0
         ! call move_dbl

         num_dbls = i

         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return

         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if

         contains

         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl

         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int

         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg

      end subroutine move_extra_info

      end module run_star_extras
