#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <string>
#include <vector>
#include <queue>
#include <cstdint>
#include "mpi.h"
using namespace std;

//prototype for data reading function
void reader(unsigned long* data, int size_out);
// prototype for multimerge function
int multimerge(unsigned long * start[], const unsigned long lengths[], const int Number, unsigned long newArray[], const int newArrayLength);
//prototypes for some utility functions
int compare_ints(const void *a, const void *b);
string sortedcheck(unsigned long xxx[], int xxxStart, int xxxLength);
string sortedcheck(unsigned long xxx[], int xxxLength);
void dumpArray(int myid, string arrayName, unsigned long array[], int start, int length);
void dumpArray(int myid, string arrayName, unsigned long array[], int length);

struct mmdata {
	int stindex;
	int index;
	int stvalue;

	mmdata(int st=0, int id=0, int stv = 0):stindex(st),index(id),stvalue(stv){}

};

void reader(unsigned long* data, int size_out){
	unsigned long element;
	int counter = 0;
	//cout << "size=" << sizeof(element) << endl;
	ifstream fo("/share/home/shared_dir/psrs_data", ios::binary);
	if (fo.is_open()) {
		cout << "Open successfully" << endl;
		/*
		fo.seekg(0, fo.end);   //to the end of file stream
		unsigned long length = fo.tellg();  //length of the stream
		fo.seekg(0, fo.beg);  //back to the beginning of file stream
		cout << "length=" << length/64<< endl;
		*/
		while (!fo.eof()) {
			counter++;
			fo.read((char*)&element, 8);
			//cout  << element << endl;
			*(data+counter)=element;
		}
		//cout << "length = " << counter << endl;
		size_out=counter;
		return ;
	}
	else
		cout << "Cannot open this file!" << endl;
		return ;
}

// comparison operator
bool operator<( const mmdata & One, const mmdata & Two){
	return One.stvalue > Two.stvalue;
}

int multimerge(int * starts[], const int lengths[], const int Number, 
			   unsigned long newArray[], const int newArrayLength){
	// Create priority queue.  There will be at most one item in the priority queue
	// for each of the Number lists.
	priority_queue< mmdata> priorities;

	// Examine each of the Number start[] lists, place the first location into 
	// the priority 	queue if the list is not empty
	for(int i=0; i<Number;i++){
		if (lengths[i]>0)
		{
			priorities.push(mmdata(i,0,starts[i][0]));
		}
	}


	// As long as priorities is not empty, pull off the top member (the smallest 
	//value from list i), push it into the newArray, and place the next element from 
	// list i in the priority queue
	int newArrayindex = 0;  // index into the merged array
	while (!priorities.empty() && (newArrayindex<newArrayLength))
	{
		// grab the smallest element, and remove it from the priority queue
		mmdata xxx = priorities.top();
		priorities.pop();

		// insert this smallest element into the merged array
		newArray[newArrayindex++] = starts[xxx.stindex][xxx.index];

		// if start[xxx.stindex] is not empty, place the next member into priority
		if ( lengths[xxx.stindex]>(xxx.index+1))
		{
			priorities.push(mmdata(xxx.stindex, xxx.index+1, 
								starts[xxx.stindex][xxx.index+1]));
		}
}

// return the logical size of the merged array
return newArrayindex;
}


// Utility function to help verify that list is sorted
string sortedcheck(int xxx[], int xxxStart, int xxxLength){
	for(int i=xxxStart; i<xxxStart+xxxLength-1;i++){
		if (xxx[i]>xxx[i+1])
			return "NOT Sorted！";
	}

	return "Sorted！";
}


string sortedcheck(int xxx[], int xxxLength){
	return sortedcheck(xxx,0,xxxLength);
}

// ----------------

// Utility to show array values
void dumpArray(int myid, string arrayName, int array[], int start, int length){
	for(int i=start;i<start+length;i++)
	{
		cout << myid << ": " << arrayName << "[" << i << "] = " << array[i] << endl;
	}
	return;
}

void dumpArray(int myid, string arrayName, int array[], int length){ 
	dumpArray(myid, arrayName, array, 0, length);
	return;
}


// The main program that is executed on each processor

int main(){
	// processor rank, and total number of processors
	int myid, numprocs;
	
	// for timing used by root processor (#0)
	double startwtime = 0.0, endwtime;


	// *******************************************
	//
	// PHASE I:  Initialization
	MPI_Init(NULL, NULL);
	MPI_Comm comm=MPI_COMM_WORLD;
	MPI_Comm_size(comm, &numprocs);
	MPI_Comm_rank(comm, &myid);
	//numprocs=MPI_COMM_WORLD.Get_size();
	//myid = MPI_COMM_WORLD.Get_rank();

	unsigned long myDataSize;
	//store read data from file
	unsigned long* myData;
	//length of data of each process
	int myDataLengths[numprocs];
	//offset of data of each process
	int myDataStarts[numprocs];
	// communication buffer used for determination of pivot values
	unsigned long pivotbuffer[numprocs*numprocs];

	// Process #0 reads data
	// start timer!
	if (myid == 0){
		//call data reading function
		/*
		for(int index=0; index<myDataSize; index++){
			myData[index] = random()% 900;
		}
		*/
		reader(myData, myDataSize);
		startwtime = MPI_Wtime();
	}

	// Compute the number of data to receives in each process
	// Remained are given to the last process
	for(int i=0;i<numprocs;i++){
		myDataLengths[i] = myDataSize/numprocs;
		myDataStarts[i]= i*myDataSize/numprocs;
	}
	myDataLengths[numprocs-1]+=(myDataSize%numprocs);
	
	// *******************************************
	//
	// PHASE II:  Scatter data, local sorts and collect regular sample

	//  The data is scattered to all processors from the root processor (#0)
	//  (root processor) does "in place".
	if (myid==0)
	{
		MPI_Scatterv(myData,myDataLengths,myDataStarts,MPI_LONG,
				MPI_IN_PLACE,myDataLengths[myid],MPI_LONG,0,comm);
	}
	else
	{			
		MPI_Scatterv(myData,myDataLengths,myDataStarts,MPI_LONG,
				myData,myDataLengths[myid],MPI_LONG,0,comm);
	}
	  
	// All processors sort their piece of the data using cstdlib::quicksort
	sort(myData,myData+myDataLengths[myid]);

	// All processors collect regular samples from sorted list
	// Consider an offset to the myData[] index
	for(int index=0;index<numprocs;index++){
		pivotbuffer[index]= myData[index*myDataLengths[myid]/numprocs];
	}


	// *******************************************
	//
	// PHASE III:  Gather and merge samples, and broadcast p-1 pivots
	
	// merged list of samples
	unsigned long tempbuffer[numprocs*numprocs];	
	
	// process#0 gathers all samples got from last phase
	if (myid==0){
		MPI_Gather(MPI_IN_PLACE,numprocs,MPI_LONG,
			tempbuffer,numprocs,MPI_LONG,0,comm);
	}
	else{
		MPI_Gather(pivotbuffer,numprocs,MPI_LONG,
			tempbuffer,numprocs,MPI_LONG,0,comm);
	}

	//  Root processor multimerges the lists together and then selects
	//  final pivot values to broadcast
	if (myid == 0){
		//quiksort merged list of samples
		sort(tempbuffer, tempbuffer+numprocs*numprocs);
		// regularly select numprocs-1 of pivot candidates to broadcast
		// as partition pivot values for myData
		for(int i=0; i<numprocs-1; i++){
			pivotbuffer[i] = tempbuffer[(i+1)*numprocs];
		}				
	}

	//process #0 broadcasts partition pivots
	MPI_Bcast(pivotbuffer,numprocs-1,MPI_LONG,0,comm);


	// *******************************************
	//
	// PHASE IV: Local data partitioned

	// Partition information for myData[]: 
	// 		index of beginning of ith class is classStart[i],
	//		length of ith class is classLength[i], and
	// 		members of ith class, myData[j], have the property
	//   		pivotbuffer[i-1]<= myData[j] < pivotbuffer[i]
	int classStart[numprocs];
	int classLength[numprocs];
	
	// need for each processor to partition its list using the values
	// of pivotbuffer
	int dataindex=0;
	for(int classindex=0; classindex<numprocs-1; classindex++){
		classStart[classindex] = dataindex;
		classLength[classindex]=0;

		// as long as dataindex refers to data in the current class
		while((dataindex< myDataLengths[myid]) 
			&& (myData[dataindex]<=pivotbuffer[classindex])){
			classLength[classindex]++;
			dataindex++;
		}		
	}
	// last class
	classStart[numprocs-1] = dataindex;
	classLength[numprocs-1] = myDataLengths[myid] - dataindex;
	
	
	// *******************************************
	//
	// PHASE V:  All ith classes are gathered by processor i 
	
	// buffer to hold merged i^th class
	int recvbuffer[myDataSize];    
	// length of merged i^th class
	int recvLengths[numprocs];     
	int recvStarts[numprocs];   
	 
	for(int iprocessor=0; iprocessor<numprocs; iprocessor++){
		//before merging
		//Gather the length of each class in each process
		MPI_Gather(&classLength[iprocessor], 1, MPI_LONG, 
			recvLengths,1,MPI_LONG,iprocessor,comm);
			
		if (myid == iprocessor){
			recvStarts[0]=0;
			for(int i=1;i<numprocs; i++){
				recvStarts[i] = recvStarts[i-1]+recvLengths[i-1];
			}
		}
		
		MPI_Gatherv((void*)&myData[classStart[iprocessor]],
			classLength[iprocessor],MPI_LONG,
			(void*)recvbuffer,recvLengths,recvStarts,MPI_LONG,iprocessor,comm);
	}
		
	int *mmStarts[numprocs]; // array of list starts
	for(int i=0;i<numprocs;i++){
		mmStarts[i]=recvbuffer+recvStarts[i];
	}
	multimerge(mmStarts,recvLengths,numprocs,myData,myDataSize);
	
	int mysendLength = recvStarts[numprocs-1] + recvLengths[numprocs-1];
	
	// *******************************************
	//
	// PHASE VI:  Root processor collects all the data


	int sendLengths[numprocs]; // lengths of consolidated classes
	int sendStarts[numprocs];  // starting points of classes
	// Root processor gathers up the lengths of all the data to be gathered
	MPI_Gather(&mysendLength,1,MPI_LONG,
		sendLengths,1,MPI_LONG,0,comm);

	// The root processor compute starts from lengths of classes to gather
	if (myid == 0)
	{
		sendStarts[0]=0;
		for(int i=1; i<numprocs; i++)
		{
			sendStarts[i] = sendStarts[i-1]+sendLengths[i-1];
		}	
	}

	// Now we let processor #0 gather the pieces and glue them together in
	// the right order
	int sortedData[myDataSize];
	MPI_Gatherv(myData,mysendLength,MPI_LONG,
		sortedData,sendLengths,sendStarts,MPI_LONG,0,comm);

	// the root processor prints the elapsed clock time
	if (myid == 0)
	{
		endwtime = MPI_Wtime();
		cout << "wall clock time (seconds) = " 
		     << scientific << setprecision(4) << endwtime-startwtime << endl;

		cout << "Data set " << sortedcheck(sortedData,myDataSize) << " sorted:" 
			<< endl;	     
	}
		
	// shutdown MPI on the processor
	MPI_Finalize();
	return 0;
}