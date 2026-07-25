
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

-----------------------------------------------------------------------------------------

This is the parallel program to simulate the wormhole issue in 3D condition. Of course, you can also use the program to simulate the 2D condition. 

-----------------------------------------------------------------------------------------

Before you run the program on Shaheen, there are some notices you have to pay attention to:

1. You have to set the module environment on Shaheen correctly. You can use the following commands to download proper modules:

module load bluegene

2. You have to create the directory to store the results of your program at first, since there is no system call from FORTRAN on Shaheen. For example, before you run infile1, you can use "mkdir case1" in the terminal at first.

3. You have to set the number of processors coherently in such places: Makefile, .sh file and the input file.  

-----------------------------------------------------------------------------------------

Before you run the program on Neser, there are some notices you have to pay attention to:

1. You have to set the module environment on Neser correctly. You can use the following commands to download proper modules:

module unload openmpi/1.5.4/gcc
module load intel-compilers/11.1
module load openmpi/1.6.4/intel

2. You have to create the directory to store the results of your program at first, since there is no system call from FORTRAN on Shaheen. For example, before you run infile1, you can use "mkdir case1" in the terminal at first.

3. When your program size is large, for example, 160*160 cells, you have to consider the issue of "insufficient virtual memory" on each Neser node. In such condition, you can allocate more nodes to the program and reduce the number of processes on each node.

4. You have to set the number of processors coherently in such places: Makefile, .sh file and the input file. 

-----------------------------------------------------------------------------------------

If you want to simulate with DBF framework, you must make sure that the three parameters isDarcy, isBrinkman and isForchheimer at the head of the file "DBF_resi.F90" are true. If you want to simulate with only Darcy framework, you must make sure the following things: isDarcy = .true., isBrinkman = .false. and isForchheimer = .false. at the head of the file "DBF_resi.F90".

-----------------------------------------------------------------------------------------

There are many solvers to choose. You can choose the solver in the Makefile. For example, if you want to use Hypre solve, you can uncomment the statement "SOLVER = HYPRE" in the Makefile.

-----------------------------------------------------------------------------------------

Before you run the cases, you have to remove the numbers at the tail of the name of the input files at first. 

-----------------------------------------------------------------------------------------

You can see the results in the document "case*", and "*" is a number which depends on your input file. In "case*", there is a file called "matlabplot.m". You can run the file to generate the matlab figures in the document "matlabplots". The matlab figures only show the physical results at the end of the simulation.  

You can also see the results in Tecplot. A series of .plt files have been generated in "case*". The numbers at the tail of the name of the .plt files stand for the serial number of the time step. So in Tecplot, you can see the simulation history results and not only the results at the end of the simulation.
