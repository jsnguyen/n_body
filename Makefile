FORT=gfortran
FLAGS=-Wall -fdefault-real-8

all: n_body

n_body:
	$(FORT) $(FLAGS) n_body.f95 -o n_body.exe

clean:
	rm n_body.exe
