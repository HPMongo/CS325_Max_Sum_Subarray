/*
 *      Author: Mike Sharkazy & Huy Pham
 *   	Date Created: 07/11/15
 *      Last Date Modified: 07/12/15
 *      File Name: MSS.cpp
 *     	Assignment: Project 1 - Project Group 18
 *      Overview: A program to calculate the maximum subarray sum of arrays
 *                 	using different algorithms; The program reads input from
 *                 	a text file and output to a text file; The text file
 *                 	contains lines of integers separated by commas with one
 *                 	set of square brackets around each line; The integers on
 *                 	each line are placed in an array; The maximum sum of a
 *                 	subarray of the arrays on each line are calculated; The
 *                 	original array, subarray and the sum is output to a text
 *                 	file for each array; The results of each algorithm for a
 *                 	single array are output to the file before moving to the
 *                 	next array to be calculated
 *
 *                 	The output file is always named MSS_Results.txt and all
 *                 	results are appended to the end of the file. If the
 *                 	program is run more than once, the results for each
 *                 	successive program run will be appended to
 *                 	MSS_Results.txt if it is not deleted
 *
 *                 	Format for command line arguments
 *                 	mss <file name>
 *
 *                 	The file name is the input file being used
 *
 *                 	The program should be compiled with C++11 flags; Example
 *                 	for use with the flip servers below
 *
 *                 	g++ -std=c++0x MSS.cpp -o mss
 *
 *      Input: Text file with bracketed lines of integers
 *      Output: The integers in the arrays, and the sums are output to a
 *                 	text file
 *
 *
 */
 
#include <iostream>
#include <string>
#include <chrono>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <sstream>

/*
 * The enumMSS() function calculates the maximum sum of a subarray of an array
 * through enumeration
 */
void enumMSS(const std::vector<int> &a, std::vector<int> &sa, int &max);

/*
 * The betterEnumMSS() function calculates the maximum sum of a subarray of an array
 * through better enumeration
 */
void betterEnumMSS(const std::vector<int> &a, std::vector<int> &sa, int &max);

/*
 * The dncMSS() function calculates the maximum sum of a subarray of an array
 * using divide and conquer algorithm
 */
void dncMSS(const std::vector<int> &a, std::vector<int> &sa, int &max);

/*
 * Recursive function for the divide and conquer algorithm
 */
int dncMSS_recursive(const std::vector<int> &a, int left, int right, int &startIndex, int &endIndex);

/*
 * Linear function calculate the maximum sum of a subarray of a given array
 * using linear algorithm
 */
void linearMSS(const std::vector<int> &a, std::vector<int> &sa, int &max);

/*
 * The outputToFile() function writes the vectors and sum to the output file
 */
void outputToFile(std::ofstream &out, const std::vector<int> &a,
              	const std::vector<int> &sa, const int &sum);
             	 
/*
 * The fillVector() function fills the vector and with the ints from the input
 * line
 */
void fillVector(std::string &s, std::vector<int> &a);

int main(int argc, char *argv[])
{   
	// Declare input and ouput stream objects and open input stream
	std::ifstream inStream(argv[1]);
	if (inStream.fail())
	{
    		std::cout << "Input file opening failed.\n";
    		exit(1);
	}
	std::ofstream outStream("MSS_Results.txt");//, std::ios::app);
    
	while (!inStream.eof())
	{
    	// Declare a string for the line, get next line, declare vector for ints
    		std::string line;
    		getline(inStream, line);
    		std::vector<int> intVect;
   	 
    	// Fill the vector with the ints in the line
    		fillVector(line, intVect);

    	// Declare subarray vectors and sum variable
    		std::vector<int> subArray1;
    		std::vector<int> subArray2;
		std::vector<int> subArray3;
		std::vector<int> subArray4;
    		int sum = 0;
   	 	if(line.size() > 0) {
   	
	/// Using enumeration method ///	
        		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();  //set first time point
	// Call function to calculate the maximum sum of a subarray
    			enumMSS(intVect, subArray1, sum);
   	 
   			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();  //set second time point
        		auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();      //calculate duration

        		std::cout << "Algorithm 1 took " << duration1 << " milliseconds to run." << std::endl;	
	// Write vector contents and sum to the output file
    			outputToFile(outStream, intVect, subArray1, sum);
    
	/// Using better enumeration method ///
        		std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();  //set first time point
	// Call the next function to calculate the maximum sum of a subarray
    			betterEnumMSS(intVect, subArray2, sum);
    			std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();  //set second time point
        		auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 ).count();      //calculate duration

        		std::cout << "Algorithm 2 took " << duration2 << " milliseconds to run." << std::endl;	  	
    	// Write vector contents and sum to the output file
    			outputToFile(outStream, intVect, subArray2, sum); 	 
		
	/// Using divide and conquer method ///
			std::chrono::high_resolution_clock::time_point t5 = std::chrono::high_resolution_clock::now();  //set first time point
	// Call the next function to calculate the max sum
			dncMSS(intVect, subArray3, sum);
    			std::chrono::high_resolution_clock::time_point t6 = std::chrono::high_resolution_clock::now();  //set second time point
        		auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>( t6 - t5 ).count();      //calculate duration
         		std::cout << "Algorithm 3 took " << duration3 << " milliseconds to run." << std::endl;	  	
    	// Write vector contents and sum to the output file
    			outputToFile(outStream, intVect, subArray3, sum); 	 

	/// Using linear method ///
			std::chrono::high_resolution_clock::time_point t7 = std::chrono::high_resolution_clock::now();  //set first time point
	// Call the next function to calculate the max sum
			linearMSS(intVect, subArray4, sum);
    			std::chrono::high_resolution_clock::time_point t8 = std::chrono::high_resolution_clock::now();  //set second time point
        		auto duration4 = std::chrono::duration_cast<std::chrono::microseconds>( t8 - t7 ).count();      //calculate duration
         		std::cout << "Algorithm 4 took " << duration4 << " milliseconds to run." << std::endl;	  	
    	// Write vector contents and sum to the output file
    			outputToFile(outStream, intVect, subArray4, sum); 	 




		}
    	}
	// Close input and output streams
	inStream.close();
	outStream.close();
 	std::cout << "Main program executed successfully!\n"; 
}

/*   *   *   *   *   *   *
 *
 * Function: enumMSS()
 *
 *	Entry: A const vector by reference, a vector by reference for the
 *       	subarray, and an int by reference
 *
 * 	Exit: Values of parameters sa and sum will be changed by the function        	 
 *
 *  Purpose: Calculate the maximum sum of a subarray using enumeration
 *
 *   *   *   *   *   *   */
void enumMSS(const std::vector<int> &a, std::vector<int> &sa, int &max)
{
	//std::cout << "In enumMSS \n";
	// Set the sum to 0 and declare a variable for the temp sum
	max = 0;
	int tempSum;
	int startIndex;
	int endIndex;
    
	// Loop over the indices, keeping the best sum found
	for (int i = 0; i < a.size(); i++)
	{
    		for (int j = i; j < a.size(); j++)
    		{
        		tempSum = 0;
        		for (int k = i; k <= j; k++)
        		{
            			tempSum += a[k];
            			if (tempSum > max)
            			{
                			max = tempSum;
                			startIndex = i;
                			endIndex = j;
            			}
        		}
    		}
	}
    	//std::cout << "Best sum is: " << max << std::endl;
	// Create the subarray
	for (int i = startIndex; i <= endIndex; i++)
	{
    		sa.push_back(a[i]);
	}
	//std::cout << "Best sub array is: " << std::endl;
	//for (int i = 0; i < sa.size(); i++) 
	//{
	//	std::cout << sa.at(i) << std::endl; 
	//}
}

/*   *   *   *   *   *   *
 *
 * Function: betterEnumMSS()
 *
 *	Entry: A const vector by reference, a vector by reference for the
 *       	subarray, and an int by reference
 *
 * 	Exit: Values of parameters sa and sum will be changed by the function        	 
 *
 *  Purpose: Calculate the maximum sum of a subarray using better enumeration
 *
 *   *   *   *   *   *   */
void betterEnumMSS(const std::vector<int> &a, std::vector<int> &sa, int &max)
{
	// Set the sum to 0 and declare a variable for the temp sum
	max = 0;
	int tempSum;
	int startIndex;
	int endIndex;
    
	// Loop over the indices, keeping the best sum found
	for (int i = 0; i < a.size(); i++)
	{
    		tempSum = 0;
    		for (int j = i; j < a.size(); j++)
    		{
        		tempSum += a[j];
        		if (tempSum > max)
        		{
            			max = tempSum;
            			startIndex = i;
            			endIndex = j;
        		}
    		}
	}
    
	// Create the subarray
	for (int i = startIndex; i <= endIndex; i++)
	{
    		sa.push_back(a[i]);
	}
}

/*   *   *   *   *   *   *
 *
 * Function: dncMSS()
 *
 *	Entry: A const vector by reference, a vector by reference for the
 *       	subarray, and an int by reference
 *
 * 	Exit: Values of parameters sa and sum will be changed by the function        	 
 *
 *  Purpose: Calculate the maximum sum of a subarray using divide and conquer 
 *
 *   *   *   *   *   *   */
void dncMSS(const std::vector<int> &a, std::vector<int> &sa, int &max)
{
	//std::cout << "In main dnc\n";

	int beg = 0;
	int end = a.size() - 1;
	int startIndex = 0;
	int endIndex = 0;
	
	max = dncMSS_recursive(a, beg, end, startIndex, endIndex);

	// Create the subarray
	for (int i = startIndex; i <= endIndex; i++)
	{
		sa.push_back(a[i]);
	}
}

int dncMSS_recursive(const std::vector<int> &a, int left, int right, int &startIndex, int &endIndex)
{
	//std::cout << "In recursive dnc \n";
	//Base case
	if(left == right)
	{
		startIndex = left;
		endIndex = right;
		return a[left];
	}
	int l_idx_start, l_idx_end, r_idx_start, r_idx_end, m_idx_start, m_idx_end;
	int max = 0;
	int mid = (left + right) / 2;
    	int maxLeft = dncMSS_recursive(a, left, mid, l_idx_start, l_idx_end);
	int maxRight = dncMSS_recursive(a, mid + 1, right, r_idx_start, r_idx_end);
	int maxMidLeft = 0, maxMidRight = 0, leftSum = 0, rightSum = 0;

	// Calculate maximum subarray of mid left side
	for (int i = mid; i >= left; i--)
	{
		leftSum += a[i];
		if(leftSum > maxMidLeft)
		{
			maxMidLeft = leftSum;
			m_idx_start = i;
		}	
	}
    	// Calculate maximum subarray of mid right side
    	for (int i = mid + 1; i <= right; i++)
	{
		rightSum += a[i];
		if(rightSum > maxMidRight)
		{
			maxMidRight = rightSum;
			m_idx_end = i;
		}
	}
	
	if(maxLeft > maxRight)
	{
		if(maxLeft > maxMidLeft + maxMidRight)
		{
			startIndex = l_idx_start;
			endIndex = l_idx_end;
			return maxLeft;
		}
		else
		{
			startIndex = m_idx_start;
			endIndex = m_idx_end;
			return maxMidLeft + maxMidRight;
		}
	}
	else 
	{
		if(maxRight > maxMidLeft + maxMidRight)
		{
			startIndex = r_idx_start;
			endIndex = r_idx_end;
			return maxRight;
		}
		else 
		{
			startIndex = m_idx_start;
			endIndex = m_idx_end;
			return maxMidLeft + maxMidRight;
		}
	}
}

/*   *   *   *   *   *   *
 *
 * Function: linearMSS()
 *
 *	Entry: A const vector by reference, a vector by reference for the
 *       	subarray, and an int by reference
 *
 * 	Exit: Values of parameters sa and sum will be changed by the function        	 
 *
 *  Purpose: Calculate the maximum sum of a subarray using linear algorithm
 *
 *   *   *   *   *   *   */
void linearMSS(const std::vector<int> &a, std::vector<int> &sa, int &max)
{
	//std::cout << "In linear\n";
	int startIndex = 0;
	int endIndex = 0;
	int currentMax = 0;
	max = 0;
	int i_low = 0;
	int i_high = 0;;

	for (int i = 1; i < a.size() ; i++)
	{
		i_high = i;
		if(currentMax > 0)
		{
			currentMax += a.at(i); 
		}
		else
		{
			i_low = i;
			currentMax = a.at(i);
		}

		if(currentMax > max)
		{
			max = currentMax;
			startIndex = i_low;
			endIndex = i_high;
		}
	}
	// Create the subarray
	for (int i = startIndex; i <= endIndex; i++)
	{
		sa.push_back(a[i]);
	}
}
/*   *   *   *   *   *   *
 *
 * Function: outputToFile()
 *
 *	Entry: An ofstream object by reference, a const vector by reference, a  
 *       	const vector by reference for the subarray, and an int by reference
 *
 * 	Exit: Values of parameters will be written to the output file        	 
 *
 *  Purpose: Write vectors and sums to the output file
 *
 *   *   *   *   *   *   */
void outputToFile(std::ofstream &out, const std::vector<int> &a,
              	const std::vector<int> &sa, const int &sum)
{
	//std::cout << "In output - \n";
	// Write vector contents and sum to the output stream
    	out << "[";
    	for (int i = 0; i < a.size(); i++)
    	{
        	if (i < a.size() - 1)
        	{
            		out << a[i] << ", ";
        	}
        	else
        	{
            		out << a[i] << "]\n";
        	}
    	}
    	out << "[";
    	for (int i = 0; i < sa.size(); i++)
    	{
        	if (i < sa.size() - 1)
        	{
            		out << sa[i] << ", ";
        	}
        	else
        	{
            		out << sa[i] << "]\n";
        	}
    	}
    	out << sum << "\n\n";
}


void fillVector(std::string &s, std::vector<int> &a)
{
	// Declare stringstream object for line and int for each number
	std::stringstream lineSS(s);
	int num;

	// Check if the first char in the line is the left bracket
	if (lineSS.peek() == '[')
	{
    		lineSS.ignore();
    		while (lineSS >> num)
    		{
        	// Add the number to the vector
        		a.push_back(num);
       	 
        	// Ignore commas or spaces
        		if (lineSS.peek() == ',' || lineSS.peek() == ' ')
        		{
            			lineSS.ignore();
        		}
    		}
		//std::cout << "Vector size: " << a.size() << std::endl;
		//for(int i = 0; i<a.size(); i++){
		//	std::cout << a.at(i) << std::endl;
		//}
	}
	else if(!lineSS.eof())
	{
    		std::cout << "Error: First character in line not a left bracket.\n";
		//std::cout << "'" << lineSS.peek() << "'\n";
    		exit(1);
	}
}
