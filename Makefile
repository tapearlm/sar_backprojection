all: unit_tests.cu 
	nvcc unit_tests.cu -dc -o unit_tests.exe 
