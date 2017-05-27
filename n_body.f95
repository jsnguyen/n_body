MODULE t_vector

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: vector, vinit, OPERATOR(+), OPERATOR(-), ASSIGNMENT(=), vdot, vmag, vscalarmult, v_print

  TYPE vector
    REAL :: x,y,z
  END TYPE vector

  INTERFACE OPERATOR (+)
    MODULE PROCEDURE vadd
  END INTERFACE
  INTERFACE OPERATOR (-)
    MODULE PROCEDURE vsub
  END INTERFACE
  INTERFACE ASSIGNMENT (=)
    MODULE PROCEDURE vassign
  END INTERFACE

CONTAINS

  SUBROUTINE vinit(this)
    TYPE(vector), INTENT(INOUT) :: this
    this%x = 0.0
    this%y = 0.0
    this%z = 0.0
  END SUBROUTINE vinit

  FUNCTION vadd(v1,v2) RESULT(v3)
    TYPE(vector), INTENT(IN) :: v1,v2
    TYPE(vector) :: v3
    v3%x = v1%x + v2%x
    v3%y = v1%y + v2%y
    v3%z = v1%z + v2%z
  END FUNCTION vadd

  FUNCTION vsub(v1,v2) RESULT(v3)
    TYPE(vector), INTENT(IN) :: v1,v2
    TYPE(vector) :: v3
    v3%x = v1%x - v2%x
    v3%y = v1%y - v2%y
    v3%z = v1%z - v2%z
  END FUNCTION vsub

  SUBROUTINE vassign(v1,v2)
    TYPE(vector), INTENT(OUT) :: v1
    TYPE(vector), INTENT(IN) :: v2
    v1%x = v2%x
    v1%y = v2%y
    v1%z = v2%z
  END SUBROUTINE vassign

  FUNCTION vdot(v1,v2) RESULT(d)
    TYPE(vector), INTENT(IN) :: v1,v2
    REAL :: d
    d =  v1%x * v2%x + v1%y * v2%y + v1%z * v2%z
  END FUNCTION vdot

  FUNCTION vmag(this) RESULT(m)
    TYPE(vector), INTENT(IN) :: this
    REAL :: m
    m = sqrt(this%x*this%x +this%y*this%y +this%z*this%z)
  END FUNCTION vmag

  FUNCTION vscalarmult(scalar,this) RESULT(v)
    REAL, INTENT(IN) :: scalar
    TYPE(vector), INTENT(IN) :: this
    TYPE(vector) :: v
    v%x = scalar*this%x
    v%y = scalar*this%y
    v%z = scalar*this%z
  END FUNCTION vscalarmult

  SUBROUTINE v_print(this)
    TYPE(vector), INTENT(IN) :: this
    PRINT *,this%x,this%y,this%z
  END SUBROUTINE v_print

END MODULE t_vector

MODULE class_particle
  USE t_vector
  IMPLICIT NONE
  private

  TYPE,PUBLIC :: particle
    TYPE(vector) :: pos,vel
    REAL :: mass
    CONTAINS
      PROCEDURE :: init => p_init
      PROCEDURE :: print => p_print
      PROCEDURE :: set_pos => p_set_pos
      PROCEDURE :: set_vel => p_set_vel
      PROCEDURE :: set_mass => p_set_mass
  END TYPE particle


CONTAINS

  SUBROUTINE p_init(this)
    CLASS(particle), INTENT(INOUT) :: this
    call this%set_pos(0.0,0.0,0.0)
    call this%set_vel(0.0,0.0,0.0)
  END SUBROUTINE p_init

  SUBROUTINE p_set_pos(this,x,y,z)
    CLASS(particle), INTENT(INOUT) :: this
    REAL, INTENT(IN) :: x,y,z
    this%pos%x = x
    this%pos%y = y
    this%pos%z = z
  END SUBROUTINE p_set_pos

  SUBROUTINE p_set_vel(this,x,y,z)
    CLASS(particle), INTENT(INOUT) :: this
    REAL, INTENT(IN) :: x,y,z
    this%vel%x = x
    this%vel%y = y
    this%vel%z = z
  END SUBROUTINE p_set_vel

  SUBROUTINE p_set_mass(this,m)
    CLASS(particle), INTENT(INOUT) :: this
    REAL, INTENT(IN) :: m
    this%mass = m
  END SUBROUTINE p_set_mass

  SUBROUTINE p_print(this)
    CLASS(particle), INTENT(IN) :: this
    PRINT *,'particle position: (',  this%pos%x, this%pos%y, this%pos%z ,')'
    PRINT *,'particle velocity: (',  this%vel%x, this%vel%y, this%vel%z ,')'
    PRINT *,'particle mass: ',this%mass
  END SUBROUTINE p_print

END MODULE class_particle

MODULE eom !equations of motion
  USE class_particle
  USE t_vector

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calc_fg, calc_vf, calc_d, update_image

  REAL, PARAMETER :: G = 6.67408E-11

CONTAINS

  FUNCTION calc_fg(p1,p2) RESULT(acc)
    CLASS(particle), INTENT(IN) :: p1,p2
    TYPE(vector) :: acc
    REAL :: c, mag

    mag = vmag(p2%pos-p1%pos)
    IF (mag < 10000.0) THEN
      mag = 10000.0
    ENDIF

    c = (G*p1%mass*p2%mass) / (mag**3)
    acc = vscalarmult( c , (p2%pos-p1%pos) )
  END FUNCTION calc_fg

  FUNCTION calc_vf(vi,acc,t) RESULT(vf)
    TYPE(vector), INTENT(IN) :: vi,acc
    TYPE(vector) :: vf
    REAL, INTENT(IN) :: t
    vf = vi + vscalarmult(t,acc)
  END FUNCTION calc_vf

  FUNCTION calc_d(vi,acc,t) RESULT(d)
    TYPE(vector), INTENT(IN) :: vi,acc
    REAL, INTENT(IN) :: t
    TYPE(vector) :: d
    d = vscalarmult(t,vi) + vscalarmult(0.5*t*t,acc)
  END FUNCTION calc_d

  SUBROUTINE update_image(p_real,p_img,n_particles,box_size)
    TYPE(particle), DIMENSION(1:n_particles), INTENT(IN) :: p_real
    TYPE(particle), DIMENSION(1:6,1:n_particles), INTENT(OUT) :: p_img
    INTEGER, INTENT(IN) :: n_particles
    INTEGER :: i
    REAL, INTENT(IN) ::  box_size


    DO i=1,n_particles
      p_img(1,i) = p_real(i)
      CALL p_img(1,i)%set_pos(p_img(1,i)%pos%x+box_size/2.0,p_img(1,i)%pos%y,p_img(1,i)%pos%z)
    END DO
    DO i=1,n_particles
      p_img(2,i) = p_real(i)
      CALL p_img(2,i)%set_pos(p_img(2,i)%pos%x-box_size/2.0,p_img(2,i)%pos%y,p_img(2,i)%pos%z)
    END DO
    DO i=1,n_particles
      p_img(3,i) = p_real(i)
      CALL p_img(3,i)%set_pos(p_img(3,i)%pos%x,p_img(3,i)%pos%y+box_size/2.0,p_img(3,i)%pos%z)
    END DO
    DO i=1,n_particles
      p_img(4,i) = p_real(i)
      CALL p_img(4,i)%set_pos(p_img(4,i)%pos%x,p_img(4,i)%pos%y-box_size/2.0,p_img(4,i)%pos%z)
    END DO
    DO i=1,n_particles
      p_img(5,i) = p_real(i)
      CALL p_img(5,i)%set_pos(p_img(5,i)%pos%x,p_img(5,i)%pos%y,p_img(5,i)%pos%z+box_size/2.0)
    END DO
    DO i=1,n_particles
      p_img(6,i) = p_real(i)
      CALL p_img(6,i)%set_pos(p_img(6,i)%pos%x,p_img(6,i)%pos%y,p_img(6,i)%pos%z-box_size/2.0)
    END DO
  END SUBROUTINE update_image


END MODULE eom

PROGRAM n_body
  USE class_particle
  USE t_vector
  USE eom
  IMPLICIT NONE

  REAL :: sim_time = 10
  REAL, PARAMETER :: t_step = 0.01
  INTEGER :: n_time_step

  TYPE(vector) :: f,a,a_sum

  INTEGER, PARAMETER :: n_particles = 2000
  REAL, PARAMETER :: box_size = 1000000.0
  INTEGER :: i,j,k,time, seed
  TYPE(particle), DIMENSION(1:n_particles) :: p_real
  TYPE(particle), DIMENSION(1:6,1:n_particles) :: p_img

  CALL vinit(a_sum)
  n_time_step = int(sim_time/t_step)
  seed =  int(secnds(0.0_4)*1000)
  CALL srand(seed)

  DO i=1,n_particles
    CALL p_real(i)%init
    CALL p_real(i)%set_pos(((rand()*2.0) - 1.0)*box_size/2.0,((rand()*2.0) - 1.0)*box_size/2.0,((rand()*2.0) - 1.0)*box_size/2.0)


    CALL p_real(i)%set_vel(((rand()*2.0) - 1.0)*10000,((rand()*2.0) - 1.0)*10000,((rand()*2.0) - 1.0)*10000)
    CALL p_real(i)%set_mass(rand()*1E24)
  END DO

  CALL p_real(1)%set_mass(1E25)
  CALL p_real(2)%set_mass(1E25)
  CALL p_real(3)%set_mass(1E25)
  CALL p_real(4)%set_mass(1E25)
  CALL p_real(5)%set_mass(1E25)

  CALL update_image(p_real,p_img,n_particles,box_size)


OPEN(UNIT = 1, file = "particles.dat")

DO time=1,n_time_step
  DO i=1,n_particles
    WRITE(1,*) p_real(i)%pos%x, p_real(i)%pos%y, p_real(i)%pos%z

    DO j=1,6
      DO k =1,n_particles
        IF (i.NE.k) THEN
          IF ( vmag( p_real(i)%pos-p_img(j,k)%pos ) .LE. box_size/2.0 ) THEN
            f = calc_fg(p_real(i),p_img(j,k))
            a = vscalarmult(1.0/p_real(i)%mass,f)
            a_sum = a+a_sum
          ENDIF
        ENDIF
      END DO
    END DO

    DO j=1,n_particles
      IF (i.NE.j) THEN
        IF ( vmag( p_real(i)%pos-p_real(j)%pos ) .LE. box_size/2.0 ) THEN
          f = calc_fg(p_real(i),p_real(j))
          a = vscalarmult(1.0/p_real(i)%mass,f)
          a_sum = a+a_sum
        ENDIF
      ENDIF
    END DO

    p_real(i)%pos = p_real(i)%pos + calc_d(p_real(i)%vel, a_sum, t_step)
    p_real(i)%vel = calc_vf(p_real(i)%vel, a_sum, t_step)

    IF (p_real(i)%pos%x > box_size/2.0 .OR.  p_real(i)%pos%x < -box_size/2.0) THEN
      p_real(i)%pos%x = -p_real(i)%pos%x
    ENDIF

    IF (p_real(i)%pos%y > box_size/2.0 .OR.  p_real(i)%pos%y < -box_size/2.0) THEN
      p_real(i)%pos%y = -p_real(i)%pos%y
    ENDIF

    IF (p_real(i)%pos%z > box_size/2.0 .OR.  p_real(i)%pos%z < -box_size/2.0) THEN
      p_real(i)%pos%z = -p_real(i)%pos%z
    ENDIF

    CALL vinit(a_sum)


  END DO
  CALL update_image(p_real,p_img,n_particles,box_size)
  WRITE(1,*) time

END DO
CLOSE(1)


END PROGRAM n_body
