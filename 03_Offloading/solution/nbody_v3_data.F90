!*****************************************************************************!
!* Copyright 2013-2015 Intel Corporation All Rights Reserved.                *!
!*                                                                           *!
!* The source code, information and material ("Material") contained herein   *!
!* is owned by Intel Corporation or its suppliers or licensors, and title to *!
!* such Material remains with Intel Corporation or its suppliers or          *!
!* licensors.                                                                *!
!* The Material contains proprietary information of Intel or its suppliers   *!
!* and licensors. The Material is protected by worldwide copyright laws and  *!
!* treaty provisions. No part of the Material may be used, copied,           *!
!* reproduced, modified, published, uploaded, posted, transmitted,           *!
!* distributed or disclosed in any way without Intel's prior express written *!
!* permission. No license under any patent, copyright or other intellectual  *!
!* property rights in the material is granted to or conferred upon you,      *!
!* either expressly, by implication, inducement, estoppel or otherwise. Any  *!
!* license under such intellectual property rights must be express and       *!
!* approved by Intel in writing.                                             *!
!*****************************************************************************!

! Allow selection of floating-point type at compile-time.
#define FPTYPE real
#define SQRT sqrt

module BodyData
  ! Arrays of body data.
  FPTYPE, allocatable :: Position(:), Acceleration(:)
  FPTYPE, allocatable :: Position_X(:), Position_Y(:), Position_Z(:), &
                         Acceleration_X(:), Acceleration_Y(:), Acceleration_Z(:)
  FPTYPE, allocatable :: Mass(:)

  !$omp declare target(Position_X, Position_Y, Position_Z, Mass, &
  !$omp                Acceleration)

  !dir$ attributes align:64 :: Position_X, Position_Y, Position_Z, Mass, Acceleration_X, Acceleration_Y, Acceleration_Z
  integer :: number_of_bodies
  FPTYPE :: epsilon_sqr
  parameter (epsilon_sqr = 0.01)
end module BodyData


subroutine Initialize()
  use BodyData
  implicit none
  integer :: i

  allocate(Mass(1:number_of_bodies))

  ! Aos format.
  allocate(Position(1:3*number_of_bodies), &
           Acceleration(1:3*number_of_bodies))

  ! SoA format.
  allocate(Position_X(1:number_of_bodies), &
           Position_Y(1:number_of_bodies), &
           Position_Z(1:number_of_bodies))
  allocate(Acceleration_X(1:number_of_bodies), &
           Acceleration_Y(1:number_of_bodies), &
           Acceleration_Z(1:number_of_bodies))

  ! Initialize arrays
  do i = 1, 3 * number_of_bodies
    Position(i) = mod(i-1, 5)
    Acceleration(i) = 0
  end do

  do i = 1, number_of_bodies
    Mass(i) = mod(i-1, 4)
    Position_X(i) = Position(3*(i-1)+1)
    Position_Y(i) = Position(3*(i-1)+2)
    Position_Z(i) = Position(3*(i-1)+3)
    Acceleration_X(i) = Acceleration(3*(i-1)+1)
    Acceleration_Y(i) = Acceleration(3*(i-1)+2)
    Acceleration_Z(i) = Acceleration(3*(i-1)+3)
  end do
end subroutine Initialize


subroutine Perform_NBody()
  use BodyData
  implicit none
  integer :: i, j, jj
  integer :: num_bodies
  FPTYPE :: eps_sqr
  integer, parameter :: CHUNK_SIZE = 8192
  FPTYPE :: pos_x, pos_y, pos_z
  FPTYPE :: acc_x, acc_y, acc_z
  FPTYPE :: delta_x, delta_y, delta_z, gamma, s

  eps_sqr = epsilon_sqr
  num_bodies = number_of_bodies

  !$omp target map(Position_X(1:num_bodies), &
  !$omp        Position_Y(1:num_bodies), Position_Z(1:num_bodies), &
  !$omp        Mass(1:num_bodies), Acceleration(1:3*num_bodies))

  !$omp parallel private(i,j,jj,pos_x,pos_y,pos_z,acc_x,acc_y,acc_z,&
  !$omp                  delta_x,delta_y,delta_z,gamma,s)
  do jj = 1, num_bodies, CHUNK_SIZE
    !$omp do
    do i = 1, num_bodies
      pos_x = Position_X(i)
      pos_y = Position_Y(i)
      pos_z = Position_Z(i)
      acc_x = 0
      acc_y = 0
      acc_z = 0

      !$omp simd private(delta_x,delta_y,delta_z,gamma,s) &
      !$omp      reduction(+:acc_x, acc_y, acc_z) &
      !$omp      aligned(Position_X, Position_Y, Position_Z, Mass :64)
      do j = jj, min(jj + CHUNK_SIZE, num_bodies)
        delta_x = Position_X(j) - pos_x
        delta_y = Position_Y(j) - pos_y
        delta_z = Position_Z(j) - pos_z

        gamma = delta_x*delta_x + delta_y*delta_y + &
                delta_z*delta_z + eps_sqr
        s = Mass(j) / (gamma * SQRT(gamma))
        acc_x = acc_x + s * delta_x
        acc_y = acc_y + s * delta_y
        acc_z = acc_z + s * delta_z
      end do

      Acceleration(3*(i-1)+1) = Acceleration(3*(i-1)+1) + acc_x
      Acceleration(3*(i-1)+2) = Acceleration(3*(i-1)+2) + acc_y
      Acceleration(3*(i-1)+3) = Acceleration(3*(i-1)+3) + acc_z
    end do
  end do
  !$omp end parallel

  !$omp end target
end subroutine Perform_NBody


subroutine Checking()
  use BodyData
  implicit none
  integer :: i

  i = (number_of_bodies / 2)+1
  write(*,"(A7,F,F,F)") "Check =", Acceleration(3*(i-1)+1), &
                                   Acceleration(3*(i-1)+2), &
                                   Acceleration(3*(i-1)+3)
end subroutine Checking


program nbody
  use omp_lib
  use BodyData
  implicit none
  character(len=100) :: arg
  double precision :: t0, t1

  if (iargc() .NE. 1) then
    write(*,*) "Usage: nbody [number of bodies]"
    stop
  endif

  call getarg(1, arg)
  read(arg,*) number_of_bodies

  call Initialize()

  !$omp target data map(to: Position_X(1:number_of_bodies), &
  !$omp                     Position_Y(1:number_of_bodies), & 
  !$omp                     Position_Z(1:number_of_bodies), &
  !$omp                     Mass(1:number_of_bodies)) &
  !$omp             map(tofrom:Acceleration(1:3*number_of_bodies))

  ! Reinitialize the output arrays after the warmup 
  ! run to get the right debug sum in the end
  call Perform_NBody()
  Acceleration(1:3*number_of_bodies) = 0
  
  !$omp target update to(Acceleration(1:3*number_of_bodies))

  t0 = omp_get_wtime()
  call Perform_NBody()
  t1 = omp_get_wtime()

  !$omp end target data

  call Checking()

  write(*,'(A,ES12.5,A)') "Run time = ", t1-t0, " seconds."
end program nbody
