test :
	set -e;
	make clean;
	pytest;

unit :
	set -e;
	make clean;
	pytest -v -m unit;

integration :
	set -e;
	make clean;
	pytest -v -m integration;

clean :
	set -e; 	
	if [ -d "scratch" ]; then rm -Rf scratch; fi
	mkdir scratch;
