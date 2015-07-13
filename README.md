# CS325_Max_Sum_Subarray
This program will calculate the maximum sum of a sub-array for a given integer array. The program will implement four different algorithms.

Once the program is compiled, the code can be executed against an input file that has one or many integer arrays. The program will execute four algorithms (enumeration, better enumeration, divide and conquer and linear) on  each array and determine the maximum sum of the sub array. The result is written to output file 'MSS_Results.txt'.

For testing purposes, different input files are provided such as 'test1.txt', 'MSS_100.txt' - 'MSS_1000.txt'. 'test1.txt' consists of a small set of arrays. 'MSS_100.txt' to 'MSS_1000.txt' consists of 10 arrays in each file where each array has 100 to 1000 random integers.

Example:
1) Compile the program - 'g++ -std=c++0x MSS.cpp -o mss'
2) Execute the program against 'test1.txt' file - 'mss test1.txt'
3) Run-time for each algorithm will be displayed on the screen and the results will be stored under 'MSS_Results.txt'.
