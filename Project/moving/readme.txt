//readme.txt

1. $gcc -c functionname.c $(-lm if using math.h)
//builds the object file for any new functions and links it with its own header file, //functionname.h

2. $add functionname.o to OBJS in Makefile
//adds the object file to the list of OBJS in the library

3. $make
//builds the libgaussq.a library containing all of the object files

4. $make mainname
//links the main program with the library

5. $./mainname
//runs the executable

//note: to just run one program still do gcc -o functionname functionname.c -lm
//-lm links with the math library
