#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <algorithm>
#include "mpi.h"
#include <iostream>

using namespace std;

const double RandomMultiplier = 1000.0;
int pnum = 0;  
int prank = -1;  

enum split_st { FirstHalf, SecondHalf };

void copyD(double *pData, int DataSize, double *pDataCopy) {
	copy(pData, pData + DataSize, pDataCopy);
}
 
bool compareD(double *pData1, double *pData2, int DataSize) {
	return equal(pData1, pData1 + DataSize, pData2);
}
 
void SerialBubble(double *pData, int DataSize) {
	double Tmp;
	for (int i = 1; i < DataSize; i++)
		for (int j = 0; j < DataSize - i; j++)
			if (pData[j] > pData[j + 1]) {
				Tmp = pData[j];
				pData[j] = pData[j + 1];
				pData[j + 1] = Tmp;
			}
}
 
void printD(double *pData, int DataSize) {
	for (int i = 0; i < DataSize; i++)
		printf("%7.4f ", pData[i]);
	printf("\n");
}

void RandomData(double *&pData, int& DataSize) {
	srand((unsigned)time(0));
	for (int i = 0; i < DataSize; i++)
		pData[i] = double(rand()) / RAND_MAX * RandomMultiplier;
}

void initialize(double *&pData, int& DataSize, double
	*&pProcData, int& BlockSize) {
	setvbuf(stdout, 0, _IONBF, 0);
	if (prank == 0) {
		do {
			printf("Enter the size of data: ");
			scanf_s("%d", &DataSize);
			if (DataSize < pnum)
				printf("Data size should be greater than number of processes\n");
		} while (DataSize < pnum);
	}
	 
	MPI_Bcast(&DataSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int RestData = DataSize;
	for (int i = 0; i < prank; i++)
		RestData -= RestData / (pnum - i);
	BlockSize = RestData / (pnum - prank);
	pProcData = new double[BlockSize];
	if (prank == 0) {
		pData = new double[DataSize];
		RandomData(pData, DataSize);
		 
	}
}

void term(double *pData, double *pProcData) {
	if (prank == 0)
		delete[]pData;
	delete[]pProcData;
}

void dataD(double *pData, int DataSize, double *pProcData, int
	BlockSize) {
 
	int *pSendInd = new int[pnum];
	int *pSendNum = new int[pnum];
	int RestData = DataSize;
	int CurrentSize = DataSize / pnum;
	pSendNum[0] = CurrentSize;
	pSendInd[0] = 0;
	for (int i = 1; i < pnum; i++) {
		RestData -= CurrentSize;
		CurrentSize = RestData / (pnum - i);
		pSendNum[i] = CurrentSize;
		pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
	}
	MPI_Scatterv(pData, pSendNum, pSendInd, MPI_DOUBLE, pProcData,
		pSendNum[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
	delete[] pSendNum;
	delete[] pSendInd;
}

void collectD(double *pData, int DataSize, double *pProcData, int
	BlockSize) {
	 
	int *pReceiveNum = new int[pnum];
	int *pReceiveInd = new int[pnum];
	int RestData = DataSize;
	pReceiveInd[0] = 0;
	pReceiveNum[0] = DataSize / pnum;
	for (int i = 1; i < pnum; i++) {
		RestData -= pReceiveNum[i - 1];
		pReceiveNum[i] = RestData / (pnum - i);
		pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
	}
	MPI_Gatherv(pProcData, BlockSize, MPI_DOUBLE, pData,
		pReceiveNum, pReceiveInd, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	 
	delete[]pReceiveNum;
	delete[]pReceiveInd;
}

void exchangeD(double *pProcData, int BlockSize, int DualRank,
	double *pDualData, int DualBlockSize) {
	MPI_Status status;
	MPI_Sendrecv(pProcData, BlockSize, MPI_DOUBLE, DualRank, 0,
		pDualData, DualBlockSize, MPI_DOUBLE, DualRank, 0,
		MPI_COMM_WORLD, &status);
}

void testf(double *pData, int DataSize, double *pProcData, int
	BlockSize) {
	MPI_Barrier(MPI_COMM_WORLD);
	if (prank == 0) {
		printf("Initial data:\n");
		printD(pData, DataSize);
		cout << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i = 0; i < pnum; i++) {
		if (prank == i) {
			printf("ProcessRank = %d\n", prank);
			printf("Block:\n");
			printD(pProcData, BlockSize);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

void ParallelPrintD(double *pProcData, int BlockSize) {
 
	for (int i = 0; i < pnum; i++) {
		if (prank == i) {
			cout << endl;
			printf("ProcessRank = %d\n", prank);
			printf("Proc sorted data: \n");
			printD(pProcData, BlockSize);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

void ParallelBubble(double *pProcData, int BlockSize) {
	 
	SerialBubble(pProcData, BlockSize);
	int Offset;
	split_st splitm;
	for (int i = 0; i < pnum; i++) {
		if ((i % 2) == 1) {
			if ((prank % 2) == 1) {
				Offset = 1;
				splitm = FirstHalf;
			}
			else {
				Offset = -1;
				splitm = SecondHalf;
			}
		}
		else {
			if ((prank % 2) == 1) {
				Offset = -1;
				splitm = SecondHalf;
			}
			else {
				Offset = 1;
				splitm = FirstHalf;
			}
		}
		 
		if ((prank == pnum - 1) && (Offset == 1)) continue;

		if ((prank == 0) && (Offset == -1)) continue;

		MPI_Status status;
		int DualBlockSize;

		MPI_Sendrecv(&BlockSize, 1, MPI_INT, prank + Offset, 0,
			&DualBlockSize, 1, MPI_INT, prank + Offset, 0,
			MPI_COMM_WORLD, &status);
		double *pDualData = new double[DualBlockSize];
		double *pMergedData = new double[BlockSize + DualBlockSize];
	 
		exchangeD(pProcData, BlockSize, prank + Offset, pDualData,
			DualBlockSize);
	 
		merge(pProcData, pProcData + BlockSize, pDualData, pDualData +
			DualBlockSize, pMergedData);
		 
		if (splitm == FirstHalf)
			copy(pMergedData, pMergedData + BlockSize, pProcData);
		else
			copy(pMergedData + BlockSize, pMergedData + BlockSize +
				DualBlockSize, pProcData);
		delete[]pDualData;
		delete[]pMergedData;
	}
}

int main(int argc, char *argv[]) {
	double *pData = 0;
	double *pProcData = 0;
	int DataSize = 0;
	int BlockSize = 0;
	double *pSerialData = 0;
	double start, finish;
	double duration = 0.0;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &pnum);
	MPI_Comm_rank(MPI_COMM_WORLD, &prank);
	if (prank == 0)
		printf("Parallel bubble sort program\n");
	initialize(pData, DataSize, pProcData, BlockSize);
	if (prank == 0) {
		pSerialData = new double[DataSize];
		copyD(pData, DataSize, pSerialData);
	}
	start = MPI_Wtime();
	dataD(pData, DataSize, pProcData, BlockSize);
	testf(pData, DataSize, pProcData, BlockSize);
	ParallelBubble(pProcData, BlockSize);
	ParallelPrintD(pProcData, BlockSize);
	collectD(pData, DataSize, pProcData, BlockSize);
	finish = MPI_Wtime();
	duration = finish - start;
	if (prank == 0)
		printf("Time of execution: %f\n", duration);
	if (prank == 0)
		delete[]pSerialData;
	term(pData, pProcData);
	MPI_Finalize();
	system("pause");
	return 0;
}
